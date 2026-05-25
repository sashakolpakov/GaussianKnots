Usage
=====

Install dependencies from the repository root:

.. code-block:: sh

   python3 -m venv .venv
   . .venv/bin/activate
   pip install -r requirements.txt

The requirements install the ``sashakolpakov/pyknotid`` fork with its
``performance`` extra, which enables the Numba JIT backend when supported by the
local Python environment.

Download the optional pyknotid catalogue database if named knot identification
is needed:

.. code-block:: sh

   python3 scripts/run_knot_experiment.py --download-pyknotid-db --samples 1 --vertices 6

Run a Haar experiment:

.. code-block:: sh

   python3 scripts/run_knot_experiment.py \
     --projection-model haar \
     --vertices 5,6,7,8,9,10,11,12 \
     --samples 250 \
     --seed 20260524 \
     --output-dir results/haar_N5-12_250

Run the Gaussian comparison:

.. code-block:: sh

   python3 scripts/run_knot_experiment.py \
     --projection-model gaussian \
     --vertices 5,6,7,8,9,10,11,12 \
     --samples 250 \
     --seed 20260524 \
     --output-dir results/gaussian_N5-12_250

Outputs
-------

Each run writes:

``summary.csv``
   Per-:math:`N` knot-type rates and metric deformation means.

``type_counts.csv``
   Empirical knot-label counts by :math:`N`.

``samples_N*.csv``
   One record per sampled polygon, including Gauss codes, invariants, edge
   lengths, all-pair distance deformation, and embedding sanity checks.
