Formal verification
===================

The pinned Lean 4/mathlib project builds with:

.. code-block:: sh

   cd formal
   lake exe cache get
   lake build GaussianKnots

Lean defines the iid Gaussian product measure, row centering, measurable
proper-crossing events, the complete crossing statistic, vertex replacement,
and the six-vertex determinant signature. It proves the exact slot count,
raw--centered invariance, and

.. math::

   |C_N(x[v:=p])-C_N(x)|\leq 2(N-3).

From the concrete four-point crossing probability it derives the exact mean.
From a typed Doob certificate it derives the printed two-sided and normalized
tail bounds using mathlib Azuma--Hoeffding. For six vertices, Lean derives all
six Calvo records from one 15-sign signature and proves the circuit/intersection
criterion from a nonzero affine determinant.

The four-point/orthant formula, construction of the Doob certificate from iid
exposure, Gaussian polar/Wishart transfer, Haar beta law, PL isotopy, and
Calvo's topological theorem remain external inputs. The probability constants
and determinant deductions downstream of those inputs are kernel checked.

There are no unfinished proofs or project axioms. Root ``#print axioms``
audits report only ``propext``, ``Classical.choice``, and ``Quot.sound``.

The maintained declaration-level report is
`FORMAL_VERIFICATION.md <https://github.com/sashakolpakov/GaussianKnots/blob/main/FORMAL_VERIFICATION.md>`_.
