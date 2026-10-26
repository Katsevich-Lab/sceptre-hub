# 0001: Development priorities and branching during Bioconductor review

- **Status:** Accepted
- **Author:** Eugene Katsevich
- **Date:** 2026-10-26
- **Affects:** all

## Motivation

In autumn 2026, development across the R and Python packages spread across
many parallel branches and topics at once, and little of it reached a reviewed,
released state. At the same time, the R package entered Bioconductor review. The
project needs a clear order of priorities and a way to keep developing without
disturbing the review.

## Proposal

### Priorities

1. GPU acceleration of association testing.
2. Native data-model foundations: `SingleCellExperiment` in R, AnnData/MuData in
   Python.

Everything else waits unless an accepted proposal says otherwise. Details are in
the [roadmap](../ROADMAP.md).

### R package branching during review

- `main` is the version under Bioconductor review, plus bug fixes and responses
  to reviewer comments. Only maintainers push to the Bioconductor staging
  repository, deliberately, with a 0.99.x version bump.
- `dev` collects the refactor as a series of reviewed pull requests. Fixes on
  `main` are merged into `dev` regularly.
- After acceptance, `dev` is merged into `main`, which then corresponds to
  Bioconductor's development branch.

### Order of R work

On `main`:

1. Bug fixes identified during the internal refactor, in a few small pull
   requests grouped by area, each with a regression test.

On `dev`:

1. A test suite that pins the current behavior of `main`.
2. The internal restructuring, with results unchanged.
3. A design note for the `SingleCellExperiment`-based object, reviewed before
   any code.
4. An internal accessor layer, so the rest of the code no longer reads object
   fields directly.
5. The `SingleCellExperiment`-based object behind those accessors, with
   migration from objects saved by the released version.
6. Disk-backed counts through DelayedArray/HDF5Array.
7. 10x Genomics import.
8. Speedups, one per pull request, each with unchanged results and a timing
   table.

### Ways of working

Adopt the process in [CONTRIBUTING.md](../CONTRIBUTING.md) and the rules for
AI-assisted work in [AGENTS.md](../AGENTS.md), including the limit of one or two
open pull requests per repository.

## Alternatives considered

- **Merge the existing integration branch as a whole.** Rejected: too large to
  review, and it mixes in work outside the current priorities.
- **Move the R package to `SingleCellExperiment` during the Bioconductor
  review.** Rejected: it would change most of the code under review and likely
  restart the review.

## Validation

Each pull request states what it checks. Behavior-preserving steps must give
unchanged results on the pinned test suite; steps that change results document
the change and who is affected.

## Implementation

Tracked on the project board under the two focus areas, with milestones in each
repository.

## Decision

Accepted by Eugene Katsevich on 2026-10-26. Open for review by the other
maintainers.
