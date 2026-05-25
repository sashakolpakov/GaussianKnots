"""Sphinx configuration for GaussianKnots."""

project = "GaussianKnots"
author = "Alexander Kolpakov"
copyright = "2026, Alexander Kolpakov"

extensions = [
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "alabaster"
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"]],
        "displayMath": [["\\[", "\\]"]],
    }
}
