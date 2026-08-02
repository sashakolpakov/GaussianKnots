Experiments and reference data
==============================

The numerical tables in the :doc:`manuscript` come from 1000 Haar and 1000
Gaussian samples for every :math:`N=5,\ldots,12`, with master seed
``20260604``.  Unknown or ambiguous classifier outputs are retained rather
than silently assigned.

Exact finite :math:`N=6` audit
------------------------------

The tracked bucket table contains 500,000 generic Haar samples in 155 observed
mirror-merged, :math:`D_6`-canonical order-type buckets.  Applying the exact
Calvo sign decision gives

.. math::

   1856/500000=0.003712

trefoils, with Wilson 95 percent interval ``0.003547--0.003884``.  The
standard direct runs contain ``2/1000`` Haar and ``5/1000`` Gaussian trefoils;
their intervals are much wider and overlap the larger estimate.

The word *exact* applies to classification of each generic stored sign vector,
conditional on Calvo's theorem.  Sampling error remains in the estimated
probability.

The compact inputs and provenance notes are in
`data/reference/n6 <https://github.com/sashakolpakov/GaussianKnots/tree/main/data/reference/n6>`_.
The complete commands are in :doc:`reproducibility`.

Historical pilot
----------------

The repository retains an earlier 250-sample report for historical comparison.
Its tables and the retired metric plots are not manuscript data and are not
used as documentation defaults.
