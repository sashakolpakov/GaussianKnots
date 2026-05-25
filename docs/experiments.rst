First Experiment
================

The first tracked run used 250 samples for each ``N=5,...,12`` with seed
``20260524``.  Full command lines and knot-label details are in
``reports/haar_vs_gaussian_N5-12_250.md``.

Knot frequencies
----------------

Haar row-space model:

.. list-table::
   :header-rows: 1

   * - N
     - classified
     - nontrivial
     - trivial
     - unknown
     - nontrivial rate
   * - 5
     - 250
     - 0
     - 250
     - 0
     - 0.000
   * - 6
     - 250
     - 1
     - 249
     - 0
     - 0.004
   * - 7
     - 250
     - 1
     - 249
     - 0
     - 0.004
   * - 8
     - 250
     - 8
     - 242
     - 0
     - 0.032
   * - 9
     - 250
     - 9
     - 241
     - 0
     - 0.036
   * - 10
     - 250
     - 27
     - 223
     - 0
     - 0.108
   * - 11
     - 249
     - 29
     - 220
     - 1
     - 0.116
   * - 12
     - 248
     - 36
     - 212
     - 2
     - 0.145

Raw Gaussian model:

.. list-table::
   :header-rows: 1

   * - N
     - classified
     - nontrivial
     - trivial
     - unknown
     - nontrivial rate
   * - 5
     - 250
     - 0
     - 250
     - 0
     - 0.000
   * - 6
     - 250
     - 1
     - 249
     - 0
     - 0.004
   * - 7
     - 250
     - 6
     - 244
     - 0
     - 0.024
   * - 8
     - 250
     - 7
     - 243
     - 0
     - 0.028
   * - 9
     - 250
     - 20
     - 230
     - 0
     - 0.080
   * - 10
     - 250
     - 24
     - 226
     - 0
     - 0.096
   * - 11
     - 250
     - 24
     - 226
     - 0
     - 0.096
   * - 12
     - 250
     - 30
     - 220
     - 0
     - 0.120

Metric deformation
------------------

The original simplex has all pairwise distances equal to \(\sqrt{2}\).  The
tables below report the mean Hamiltonian edge max/min distortion, the mean
all-pair max/min distortion, the all-pair RMS distance divided by \(\sqrt{2}\),
and the scale-free all-pair min/max range after RMS normalization.

Haar row-space model:

.. list-table::
   :header-rows: 1

   * - N
     - edge max/min
     - all-pair max/min
     - RMS abs ratio
     - normalized min
     - normalized max
   * - 5
     - 1.920
     - 2.363
     - 0.866
     - 0.539
     - 1.151
   * - 6
     - 2.454
     - 3.573
     - 0.775
     - 0.414
     - 1.269
   * - 7
     - 2.941
     - 4.317
     - 0.707
     - 0.369
     - 1.358
   * - 8
     - 3.459
     - 5.707
     - 0.655
     - 0.322
     - 1.434
   * - 9
     - 3.818
     - 6.066
     - 0.612
     - 0.290
     - 1.504
   * - 10
     - 3.615
     - 6.873
     - 0.577
     - 0.264
     - 1.560
   * - 11
     - 4.576
     - 8.266
     - 0.548
     - 0.234
     - 1.606
   * - 12
     - 4.401
     - 9.123
     - 0.522
     - 0.222
     - 1.645

Raw Gaussian model:

.. list-table::
   :header-rows: 1

   * - N
     - edge max/min
     - all-pair max/min
     - RMS abs ratio
     - normalized min
     - normalized max
   * - 5
     - 3.078
     - 4.322
     - 0.972
     - 0.415
     - 1.488
   * - 6
     - 3.687
     - 5.671
     - 0.990
     - 0.342
     - 1.572
   * - 7
     - 3.773
     - 5.963
     - 0.984
     - 0.321
     - 1.617
   * - 8
     - 4.609
     - 7.432
     - 0.975
     - 0.274
     - 1.713
   * - 9
     - 4.566
     - 8.249
     - 0.987
     - 0.245
     - 1.752
   * - 10
     - 4.983
     - 8.953
     - 0.984
     - 0.237
     - 1.780
   * - 11
     - 4.995
     - 9.366
     - 0.987
     - 0.232
     - 1.808
   * - 12
     - 5.281
     - 10.775
     - 0.983
     - 0.208
     - 1.829
