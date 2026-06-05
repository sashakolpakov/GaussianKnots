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

For :math:`N=6`, there are nine non-adjacent edge-pair determinant walls.  The
full projective symmetry group of this wall arrangement has order 72 and lifts
orthogonally to :math:`H=\mathbf{1}^{\perp}`.  The subgroup coming from
Hamiltonian six-cycle relabellings is :math:`D_6` of order 12.  This subgroup
preserves the labelled knot model; the other 60 wall symmetries preserve the
wall arrangement but may change knot type.

The current numerical order-type bucket proxy records the full order type of the six
projected vertices:

.. math::

   \operatorname{sgn}\det(x_b-x_a,\;x_c-x_a,\;x_d-x_a),
   \qquad \{a,b,c,d\}\subset\{0,\ldots,5\}.

There are :math:`\binom{6}{4}=15` such signs.  They refine the nine wall signs,
which are too coarse for knot type in the sampled :math:`N=6` data.  The
resulting order-type estimate should be treated as numerical evidence for
sector volume, not as an exact chamber-volume formula: one still has to prove
connectivity of the relevant sign cells in the Stiefel/Grassmannian parameter
space and control the segment-intersection inequalities.


For :math:`N=7` and :math:`N=8`, ``gaussian_knots.general_chamber_geometry``
builds the same determinant-wall and order-type data for exploratory probes.
The corresponding repro scripts sample :math:`D_N` order-type buckets and check
a limited number of direct classifier labels against representative labels;
these diagnostics do not replace a proof of chamber connectivity.
