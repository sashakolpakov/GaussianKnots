"""Compatibility layer for pyknotid-based knot identification."""

from __future__ import annotations

import importlib
import inspect
import re
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np


@dataclass(frozen=True)
class PyknotidEnvironment:
    """Import and acceleration status for the local pyknotid installation."""

    available: bool
    version: str | None = None
    fast_backend_available: bool | None = None
    error: str | None = None


@dataclass(frozen=True)
class KnotIdentification:
    """Serializable result from one pyknotid identification attempt."""

    status: str
    is_nontrivial: bool | None
    knot_types: tuple[str, ...] = ()
    determinant: str | None = None
    alexander_roots: dict[str, str] = field(default_factory=dict)
    vassiliev_2: str | None = None
    vassiliev_3: str | None = None
    crossing_count: int | None = None
    simplified_crossing_count: int | None = None
    fast_mode_requested: bool = True
    fast_backend_available: bool | None = None
    messages: tuple[str, ...] = ()


def inspect_pyknotid_environment() -> PyknotidEnvironment:
    """Return whether pyknotid imports and whether Cython helpers are present."""

    try:
        _patch_numpy_legacy_aliases()
        module = importlib.import_module("pyknotid")
    except Exception as exc:  # pragma: no cover - depends on local install
        return PyknotidEnvironment(available=False, error=str(exc))

    version = getattr(module, "__version__", None)
    fast_backend_available: bool | None = None
    try:
        spacecurve = importlib.import_module("pyknotid.spacecurves.spacecurve")
        cython_helpers = getattr(spacecurve, "chelpers", None)
        python_helpers = getattr(spacecurve, "helpers", None)
        fast_backend_available = cython_helpers is not None and cython_helpers is not python_helpers
    except Exception:
        fast_backend_available = None

    return PyknotidEnvironment(
        available=True,
        version=str(version) if version is not None else None,
        fast_backend_available=fast_backend_available,
    )


def identify_polygon(vertices: np.ndarray, use_fast: bool = True) -> KnotIdentification:
    """Identify a closed polygon using the best available pyknotid API.

    The adapter is deliberately defensive: pyknotid has had multiple public
    signatures, its catalogue database may be absent, and its optional Cython
    helpers may or may not have compiled at installation time.
    """

    points = np.asarray(vertices, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("vertices must have shape (N, 3)")

    if len(points) < 6:
        return KnotIdentification(
            status="unknot_by_stick_number",
            is_nontrivial=False,
            fast_mode_requested=use_fast,
            messages=("Closed polygonal knots with fewer than 6 sticks are unknotted.",),
        )

    environment = inspect_pyknotid_environment()
    if not environment.available:
        return KnotIdentification(
            status="pyknotid_unavailable",
            is_nontrivial=None,
            fast_mode_requested=use_fast,
            fast_backend_available=environment.fast_backend_available,
            messages=(environment.error or "pyknotid could not be imported",),
        )

    try:
        from pyknotid.spacecurves import Knot
    except Exception as exc:  # pragma: no cover - depends on local install
        return KnotIdentification(
            status="pyknotid_import_error",
            is_nontrivial=None,
            fast_mode_requested=use_fast,
            fast_backend_available=environment.fast_backend_available,
            messages=(str(exc),),
        )

    knot = Knot(points, verbose=False)
    messages: list[str] = []
    raw_kwargs = {"mode": "use_max_jump", "include_closure": True, "try_cython": use_fast}

    crossing_count = None
    simplified_crossing_count = None
    determinant = None
    alexander_roots: dict[str, str] = {}
    vassiliev_2 = None
    vassiliev_3 = None
    knot_types: tuple[str, ...] = ()

    try:
        crossings = _call_supported(knot.raw_crossings, **raw_kwargs)
        crossing_count = int(len(crossings) // 2)
    except Exception as exc:
        messages.append(f"raw_crossings failed: {exc}")

    try:
        gauss_code = _call_supported(knot.gauss_code, **raw_kwargs)
        if hasattr(gauss_code, "simplify"):
            gauss_code.simplify()
        simplified_crossing_count = _safe_crossing_count(gauss_code)
    except Exception as exc:
        messages.append(f"gauss_code simplification failed: {exc}")

    for root in (2, 3, 4):
        try:
            value = knot.alexander_at_root(root)
            alexander_roots[str(root)] = _value_to_string(value)
            if root == 2:
                determinant = alexander_roots[str(root)]
        except Exception as exc:
            messages.append(f"alexander_at_root({root}) failed: {exc}")

    if hasattr(knot, "vassiliev_degree_2"):
        try:
            vassiliev_2 = _value_to_string(
                _call_supported(knot.vassiliev_degree_2, simplify=True, **raw_kwargs)
            )
        except Exception as exc:
            messages.append(f"vassiliev_degree_2 failed: {exc}")

    if hasattr(knot, "vassiliev_degree_3"):
        try:
            vassiliev_kwargs = dict(raw_kwargs)
            vassiliev_kwargs["try_cython"] = use_fast
            vassiliev_3 = _value_to_string(
                _call_supported(
                    knot.vassiliev_degree_3,
                    simplify=True,
                    **vassiliev_kwargs,
                )
            )
        except Exception as exc:
            messages.append(f"vassiliev_degree_3 failed: {exc}")

    try:
        identified = _identify_with_supported_api(knot, use_fast=use_fast)
        knot_types = tuple(_candidate_name(candidate) for candidate in identified)
    except Exception as exc:
        messages.append(f"catalogue identify failed: {exc}")

    classification = _classify(
        knot_types=knot_types,
        determinant=determinant,
        alexander_roots=alexander_roots,
        vassiliev_2=vassiliev_2,
        vassiliev_3=vassiliev_3,
        crossing_count=crossing_count,
        simplified_crossing_count=simplified_crossing_count,
    )

    if classification is True:
        status = "nontrivial_detected" if not knot_types else "identified_nontrivial"
    elif classification is False:
        status = "identified_unknot" if knot_types else "unknot_detected"
    else:
        status = "unknown"

    return KnotIdentification(
        status=status,
        is_nontrivial=classification,
        knot_types=knot_types,
        determinant=determinant,
        alexander_roots=alexander_roots,
        vassiliev_2=vassiliev_2,
        vassiliev_3=vassiliev_3,
        crossing_count=crossing_count,
        simplified_crossing_count=simplified_crossing_count,
        fast_mode_requested=use_fast,
        fast_backend_available=environment.fast_backend_available,
        messages=tuple(messages),
    )


def _identify_with_supported_api(knot: Any, use_fast: bool) -> list[Any]:
    identify = getattr(knot, "identify")
    signature = inspect.signature(identify)
    params = signature.parameters
    kwargs: dict[str, Any] = {}

    if "determinant" in params:
        kwargs["determinant"] = True
    if "vassiliev_2" in params:
        kwargs["vassiliev_2"] = True
    if "vassiliev_3" in params:
        kwargs["vassiliev_3"] = True
    if "alexander" in params:
        kwargs["alexander"] = False
    if "roots" in params:
        kwargs["roots"] = (2, 3, 4)
    if "min_crossings" in params:
        kwargs["min_crossings"] = True
    if "try_cython" in params:
        kwargs["try_cython"] = use_fast

    result = identify(**kwargs)
    if result is None:
        return []
    if isinstance(result, (list, tuple)):
        return list(result)
    return [result]


def _call_supported(func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    signature = inspect.signature(func)
    params = signature.parameters
    if any(param.kind == inspect.Parameter.VAR_KEYWORD for param in params.values()):
        return func(*args, **kwargs)
    filtered = {key: value for key, value in kwargs.items() if key in params}
    return func(*args, **filtered)


def _safe_crossing_count(gauss_code: Any) -> int | None:
    try:
        length = len(gauss_code)
    except Exception:
        return None
    return int(length // 2)


def _classify(
    knot_types: tuple[str, ...],
    determinant: str | None,
    alexander_roots: dict[str, str],
    vassiliev_2: str | None,
    vassiliev_3: str | None,
    crossing_count: int | None,
    simplified_crossing_count: int | None,
) -> bool | None:
    if simplified_crossing_count == 0 or crossing_count == 0:
        return False

    if _any_invariant_detects_nontrivial(determinant, alexander_roots, vassiliev_2, vassiliev_3):
        return True

    if len(knot_types) == 1:
        return not _is_unknot_name(knot_types[0])

    if knot_types and all(_is_unknot_name(name) for name in knot_types):
        return False

    return None


def _any_invariant_detects_nontrivial(
    determinant: str | None,
    alexander_roots: dict[str, str],
    vassiliev_2: str | None,
    vassiliev_3: str | None,
) -> bool:
    values = list(alexander_roots.values())
    if determinant is not None and determinant not in values:
        values.append(determinant)
    for value in values:
        numeric = _as_float(value)
        if numeric is not None and abs(numeric - 1.0) > 1e-9:
            return True
    for value in (vassiliev_2, vassiliev_3):
        numeric = _as_float(value)
        if numeric is not None and abs(numeric) > 1e-9:
            return True
    return False


def _candidate_name(candidate: Any) -> str:
    for attr_name in ("identifier", "name"):
        if hasattr(candidate, attr_name):
            attr = getattr(candidate, attr_name)
            value = attr() if callable(attr) else attr
            if value:
                return str(value)
    text = str(candidate)
    match = re.search(r"(K\d+[an]\d+|\d+_\d+|unknot)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return text


def _is_unknot_name(name: str) -> bool:
    lowered = name.strip().lower()
    return lowered in {"0_1", "unknot", "the unknot"} or "unknot" in lowered


def _value_to_string(value: Any) -> str:
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def _as_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _patch_numpy_legacy_aliases() -> None:
    """Patch aliases used by older pyknotid releases on NumPy 2.x."""

    aliases = {
        "bool": bool,
        "complex": complex,
        "float": float,
        "int": int,
    }
    for name, target in aliases.items():
        if not hasattr(np, name):
            setattr(np, name, target)
