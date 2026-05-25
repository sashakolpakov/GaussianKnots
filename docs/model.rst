Projection Models
=================

Gaussian model
--------------

The Gaussian projection model samples

.. math::

   P \in \mathbb{R}^{3\times N}, \qquad
   P_{ab} \sim \mathcal{N}(0,1/3),

and sets

.. math::

   x_i = P e_i.

Thus the projected simplex vertices are independent
:math:`\mathcal{N}(0,I_3/3)` vectors.  For every original simplex edge
:math:`\Delta_{ij}=e_i-e_j`,

.. math::

   \|P\Delta_{ij}\|^2 \sim \frac{2}{3}\chi^2_3,

so the expected squared projected distance equals the original squared
simplex distance :math:`2`.

Haar model
----------

For knot type, the singular values of :math:`P` are irrelevant.  Almost surely
:math:`P` has polar row factorization

.. math::

   P = A Q,\qquad A=(P P^T)^{1/2}\in GL^+(3,\mathbb{R}),\qquad
   Q Q^T = I_3.

The map :math:`A:\mathbb{R}^3\to\mathbb{R}^3` is an orientation-preserving
linear diffeomorphism, so it preserves ambient isotopy class.  Consequently
the knot-type law can be sampled from the Haar model.

In the code, the Haar model samples a random 3-frame in
:math:`\mathbf{1}^{\perp}\subset\mathbb{R}^N`.  Its vertex matrix satisfies

.. math::

   X^T X = I_3,\qquad \sum_{i=1}^N x_i = 0.

This model is globally scaled relative to the Gaussian model: the expected
squared projected distance is :math:`2\cdot 3/(N-1)`.  The experiment
therefore reports both absolute distance ratios and scale-free RMS-normalized
distance ratios.

Chambers
--------

The projection parameter space is divided by discriminant walls where two
non-adjacent polygon edges intersect.  Away from these walls the polygon is an
embedded stick knot and knot type is locally constant.  At a generic
codimension-one wall, exactly one crossing changes.
