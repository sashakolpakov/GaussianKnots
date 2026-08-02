Projection model
================

This page is a short guide. Definitions and proofs are in the
:doc:`manuscript`.

Gaussian and Haar forms
-----------------------

The Gaussian model samples a matrix

.. math::

   P\in\mathbb R^{3\times N},\qquad P_{ab}\sim\mathcal N(0,1/3),

and takes :math:`x_i=P e_i`. After centering, the polar row factorization is

.. math::

   P\Pi_H=AQ,\qquad QQ^T=I_3,\qquad Q\mathbf 1=0,

with :math:`A` almost surely invertible. The target-space linear map
:math:`A` preserves knot type but not metric statistics, which is why the
repository reports Gaussian and Haar distance deformation separately.

Exact marginal laws
-------------------

For a Haar row-orthonormal projection and :math:`N>4`, the manuscript proves

.. math::

   \frac{\lVert Q(e_i-e_j)\rVert^2}{2}
   \sim \mathrm{Beta}\!\left(\frac32,\frac{N-4}{2}\right).

For a generic planar shadow and :math:`N\geq4`, it proves

.. math::

   \mathbb E C_N=\frac{N(N-3)}{\pi}\arcsin(1/3)

and, for every :math:`t>0`,

.. math::

   \Pr\!\left(|C_N-\mathbb E C_N|\geq t\right)
   \leq 2\exp\!\left(-\frac{t^2}{2N(N-3)^2}\right).

Six vertices: classification versus volume
-------------------------------------------

For a generic embedded hexagon, the 15 tetrahedral orientation signs determine
the segment--triangle piercings in Calvo's three oriented algebraic
intersection numbers. Calvo's cited theorem then gives an exact
unknot/right-trefoil/left-trefoil classification. This does **not** require a
sign realization set to be connected.

Connectivity and exact integration are different questions. The measured
frequency :math:`1856/500000` is a Monte Carlo estimate of the Haar mass of the
trefoil sector. Neither the finite classifier nor Lean turns that frequency
into an exact chamber volume.

For :math:`N=7,8`, order-type buckets remain exploratory organizational tools;
the repository does not claim that they are connected knot chambers or that a
single representative labels an entire bucket.
