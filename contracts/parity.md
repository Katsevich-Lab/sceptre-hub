# Parity contract

The R package is the reference implementation. pysceptre and the Nextflow
pipeline agree with SCEPTRE when they reproduce R's outputs on a shared set of
small fixtures, within stated tolerances.

## Where things live

| Piece | Repository | Location |
| --- | --- | --- |
| Fixture inputs and expected outputs (canonical copy) | sceptre (R) | `tests/testthat/fixtures/parity/` |
| Generator script that runs R on the inputs and writes the expected outputs | sceptre (R) | A build-ignored developer script, e.g. `tools/parity/generate.R` |
| Reference stability test: current R output matches the committed fixtures | sceptre (R) | `tests/testthat/test-parity-reference.R` |
| Vendored copy of the fixtures at a pinned version, with checksums | pysceptre | `tests/parity/fixtures/` plus a sync script |
| Tolerances file | pysceptre | `tests/parity/tolerances.*` |
| Comparison tests | pysceptre; pipeline | Each repository's test suite, run on every pull request |
| Full-scale agreement checks on real data | hub | Validation and benchmark scripts with a results table |

Fixtures are small (kilobytes to a few megabytes), so they ship with the tests
and need no download during a build.

## Fixture versions

Each fixture set carries a version and records the R package version and commit
that produced it.

- **If R's output changes unintentionally,** the reference stability test fails.
  That is a bug to fix in R.
- **If R's output changes intentionally** (a bug fix or an accepted method
  change), the same pull request regenerates the fixtures and bumps their
  version. The changelog says which outputs changed and why.
- **pysceptre and the pipeline** then sync to the new version and update their
  code until their comparison tests pass again.

## Tolerance tiers

| Tier | Applies to | Agreement required |
| --- | --- | --- |
| Exact | Deterministic steps: QC, guide assignment, pair construction, model fits and test statistics given the same inputs | Equal, or equal within floating-point tolerance stated per quantity |
| Shared draws | Resampling | Both implementations are fed the same saved random draws from the fixture; resulting p-values agree within the exact tier |
| Statistical | Full analyses with each implementation's own random draws, and full-scale checks | p-values agree in distribution and discoveries largely overlap, with thresholds stated in the tolerances file |

## When the reference is wrong

R being the reference does not mean R is always right. If work on another
implementation uncovers a bug in R, fix R first in a bug-fix pull request,
regenerate the fixtures, and then let the other implementations follow.

A deliberate divergence from R (for example, GPU random draws that differ from
R's) needs an accepted proposal. It gets its own tolerance tier, recorded in the
tolerances file.
