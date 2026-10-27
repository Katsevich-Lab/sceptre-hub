# Rules for AI-assisted contributions

These rules apply to AI coding agents working on any SCEPTRE repository. Each
repository's own agent instructions point here. Where a repository has stricter
rules, those apply too.

## Before starting

- Read [ROADMAP.md](ROADMAP.md) and [CONTRIBUTING.md](CONTRIBUTING.md). Work only
  on items that serve the current priorities, unless a maintainer asks
  otherwise.
- Work from an issue. If there is none, ask a maintainer before starting.
- When the scope is unclear, stop and ask on the issue.

## Scope

- One open pull request at a time per repository. Wait for it to be merged
  before starting the next.
- No new public functions, dependencies, object formats or result columns
  without a maintainer's sign-off. Decisions that span languages need an
  accepted proposal.
- Do not add work that was not asked for, even if it seems useful. Suggest it
  as an issue instead.

## Pull requests

- Keep pull requests small and single-purpose. Behavior-preserving changes are
  separate from changes that alter results.
- State what the pull request checks and how that was verified. Report failures
  and skipped checks plainly.
- Before opening a pull request, run the checks and a fresh-context review, and
  address the findings. Follow the review stages and round limits in
  [CONTRIBUTING.md](CONTRIBUTING.md#review).
- Report progress on the issue and pull request, not in separate status or
  handoff documents.
- Write commit messages, comments and documentation in plain language a
  maintainer can read quickly.
- If a change alters R's output on the parity fixtures, say so. Regenerate the
  fixtures and bump their version in the same pull request only when the
  change is intended; otherwise fix the code. See
  [contracts/parity.md](contracts/parity.md).

## Repository hygiene

- Do not commit evidence files, logs, saved outputs or large data to a package
  repository. Benchmark results go in the benchmark table in this repository;
  raw logs stay outside git.
- Do not keep superseded code paths, snapshots or prototypes in the package
  source. Remove them, or ask.

## Never

- Push to the Bioconductor staging repository or publish a release.
- Merge your own pull request.
- Change defaults or statistical methods without an accepted proposal.
