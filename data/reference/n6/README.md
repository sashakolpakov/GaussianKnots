# N=6 reference data

These compact files are the tracked inputs for the manuscript's six-stick
audits.

- `bucket_volume.csv`: 155 $D_6$-canonical, mirror-merged order-type
  buckets from 500,000 Haar samples.
- `run_summary.json`: seed, sample count, symmetry counts, and grouped-run
  parameters.
- `direct_sample_summary.csv`: the independent 1000-sample direct check
  generated inside the grouped workflow.
- `haar_samples_N6.csv` and `gaussian_samples_N6.csv`: the N=6 rows
  from the standard manuscript-scale runs, used by
  `--validate-standard-direct`.

The grouped run's direct subsample contains seven trefoils.  It is distinct
from the standard manuscript runs, which contain two Haar and five Gaussian
trefoils.  All three are finite Monte Carlo samples.

`calvo_n6_exact.py` applies the exact finite sign decision to
`bucket_volume.csv` and obtains 1,856 trefoils.  This verifies classification
of the stored generic sign vectors; it does not make their frequency an exact
Haar probability.

## SHA-256

```text
0afb3be96a101ef0e711bcb5f0088744a39bf3590831d2b01e25c8ca01f95759  bucket_volume.csv
ca0c41568f6f23cf233493effbb696fe5303e728851777cc0f206a6974b2499d  direct_sample_summary.csv
33ef11b781bed9b7c29c49c6967861f21d0458aabb5b0636d940ce343fba580d  gaussian_samples_N6.csv
5421ae8f141c87b3c4c6c6fe9c95cfef0daa64e676af15bbe56f628239a4a559  haar_samples_N6.csv
9a58b34d215a02df84f77692692d1ca5e6ea6453b139bd8bd0560a6b92c0afb4  run_summary.json
```
