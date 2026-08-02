Quick usage
===========

Install
-------

.. code-block:: sh

   python3 -m venv .venv
   . .venv/bin/activate
   python -m pip install -r requirements.txt

Sample projected-simplex knots
------------------------------

.. code-block:: sh

   python scripts/run_knot_experiment.py \
     --projection-model haar \
     --vertices 5,6,7,8 \
     --samples 1000 \
     --seed 20260604 \
     --output-dir results/haar_example

Use ``--projection-model gaussian`` for the raw Gaussian comparison.  Every
output directory records its seed, projection model, classifier environment,
per-sample data, summary counts, and knot-label frequencies.

If the optional pyknotid catalogue is unavailable, geometry generation can be
tested without classification:

.. code-block:: sh

   python scripts/run_knot_experiment.py \
     --vertices 6 \
     --samples 3 \
     --allow-missing-pyknotid

See :doc:`reproducibility` for the paper's exact commands and tracked-data
checks.
