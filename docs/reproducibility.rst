Reproducibility
===============

Install the numerical dependencies from the repository root:

.. code-block:: sh

   python3 -m venv .venv
   . .venv/bin/activate
   python -m pip install -r requirements-analysis.txt

Fast tracked-data audit
-----------------------

.. code-block:: sh

   python scripts/repro/calvo_n6_exact.py \
     --expected-trefoils 1856 \
     --validate-standard-direct

Exact-law simulation checks
---------------------------

.. code-block:: sh

   python scripts/repro/validate_exact_projection_laws.py \
     --crossing-samples 10000 \
     --beta-samples 50000 \
     --batch-size 500

Full sign-stream replay
-----------------------

.. code-block:: sh

   python scripts/repro/calvo_n6_exact.py \
     --expected-trefoils 1856 \
     --validate-standard-direct \
     --audit-signs

The replay regenerates all 500,000 bucket assignments and checks the 100
closest determinant walls at 100-digit decimal precision.  It can take a few
minutes.

For manuscript, figure, Monte Carlo, documentation, and Lean build commands,
see the repository's
`REPRODUCIBILITY.md <https://github.com/sashakolpakov/GaussianKnots/blob/main/REPRODUCIBILITY.md>`_.
