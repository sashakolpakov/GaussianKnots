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
     --samples 1000 \
     --seed 20260604 \
     --output-dir results/haar_N5-12_1000

Run the Gaussian comparison:

.. code-block:: sh

   python3 scripts/run_knot_experiment.py \
     --projection-model gaussian \
     --vertices 5,6,7,8,9,10,11,12 \
     --samples 1000 \
     --seed 20260604 \
     --output-dir results/gaussian_N5-12_1000

Summarize the knot-type probability estimates with Wilson intervals:

.. code-block:: sh

   python3 scripts/repro/summarize_knot_probabilities.py results/haar_N5-12_1000
   python3 scripts/repro/summarize_knot_probabilities.py results/gaussian_N5-12_1000

The current volume-first grouped order-type workflow is:

.. code-block:: sh

   python3 scripts/repro/order_type_grouped_volume.py \
     --vertices 6,7,8 \
     --samples 500000 \
     --seed 20260604 \
     --classify-top-groups 200 \
     --checks-per-group 3 \
     --direct-classify-samples 1000 \
     --output-dir results/order_type_grouped_volume_N6-8_500k

For :math:`N=6`, this run checks all observed order-type buckets and gives
``1856/500000 = 0.003712`` for the trefoil bucket, with Wilson interval
``0.003547--0.003884``.  It also records direct 1000-sample classifier checks
for :math:`N=6,7,8`.

The lower-level Grassmannian sampler is:

.. code-block:: sh

   python3 scripts/repro/run_n6_grassmannian_probe.py \
     --samples 2048 \
     --sampler halton \
     --seed 20260604 \
     --output-dir results/n6_grassmannian_halton_2048

Exploratory :math:`N=7` and :math:`N=8` wall/order-type probes are available as:

.. code-block:: sh

   python3 scripts/repro/n7_chamber_probe.py \
     --samples 5000 \
     --classify-reps 200 \
     --classify-samples 200 \
     --skip-full-wall-symmetry \
     --output-dir results/n7_chamber_probe_5000

   python3 scripts/repro/n8_chamber_probe.py \
     --samples 5000 \
     --classify-reps 200 \
     --classify-samples 200 \
     --skip-full-wall-symmetry \
     --output-dir results/n8_chamber_probe_5000

These scripts estimate :math:`D_N` order-type bucket masses and direct-label
consistency for a limited classified subset.

The N=6-only wall/order-type diagnostic computes the nine determinant-wall
symmetries, their orthogonal lift, the :math:`D_6` Hamiltonian-cycle subgroup,
and a 15-sign order-type proxy:

.. code-block:: sh

   python3 scripts/repro/n6_chamber_analysis.py \
     --samples 50000 \
     --seed 20260604 \
     --output results/n6_full_chamber_analysis_diagnostic.csv

This verifies the orthogonal lift of the N=6 wall symmetries and records the
N=6 order-type grouping structure.

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
