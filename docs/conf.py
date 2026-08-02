"""Sphinx configuration derived from the authoritative manuscript metadata."""

from pathlib import Path
import re


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = REPOSITORY_ROOT / "manuscript" / "knot_projections.tex"
MANUSCRIPT_TEXT = MANUSCRIPT.read_text(encoding="utf-8")


def latex_values(command: str) -> list[str]:
    """Extract simple, unnested LaTeX metadata commands from the manuscript."""

    return re.findall(rf"\\{command}\{{([^{{}}]+)\}}", MANUSCRIPT_TEXT)


titles = latex_values("title")
authors = latex_values("author")
dates = latex_values("date")
if len(titles) != 1 or not authors:
    raise RuntimeError(f"could not read title/authors from {MANUSCRIPT}")

project = titles[0]
author = " and ".join(authors)
release = dates[0] if dates else ""
copyright = f"2026, {author}"

extensions = ["sphinx.ext.mathjax"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_static_path = ["_static"]
html_theme = "alabaster"
html_title = project
html_theme_options = {
    "github_user": "sashakolpakov",
    "github_repo": "GaussianKnots",
    "github_button": True,
    "github_type": "star",
}
rst_epilog = f"""
.. |paper_title| replace:: {project}
.. |paper_authors| replace:: {author}
"""
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"]],
        "displayMath": [["\\[", "\\]"]],
    }
}
