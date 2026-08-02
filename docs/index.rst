|paper_title|
===============

|paper_authors|

GaussianKnots contains the manuscript, experiments, exact six-stick sign
classifier, and formal-verification project for random Hamiltonian cycles
through projected simplex vertices.

.. important::

   ``manuscript/knot_projections.tex`` is the source of truth for mathematical
   statements. These pages explain the model and reproducibility workflow; they
   do not replace the paper.

The model joins

.. math::

   x_i=P e_i\in\mathbb R^3,
   \qquad i=1,\ldots,N,

in cyclic order. Gaussian and Haar row-orthonormal projections have the same
knot-type law. The Lean development checks the concrete Gaussian crossing
statistic and sensitivity bound, derives the exact crossing mean and tail
constants from typed analytic inputs, and derives all six Calvo records from
one determinant signature.

.. toctree::
   :maxdepth: 2

   manuscript
   model
   experiments
   reproducibility
   formal_verification
   usage
