# Rules for AI-assisted contributions

These rules apply to AI coding agents working on any SCEPTRE repository. Each
repository's own agent instructions point here. Where a repository has stricter
rules, those apply too.

## Before starting

- Read [ROADMAP.md](ROADMAP.md) and [CONTRIBUTING.md](CONTRIBUTING.md). Work only
  on items that serve the current priorities, unless a maintainer asks
  otherwise.
- When the scope is unclear, stop and ask.

## Scope

- No new public functions, dependencies, object formats or result columns
  without a maintainer's sign-off.
- Do not add work that was not asked for, even if it seems useful. Suggest it
  instead.

## Changes

- Keep behavior-preserving changes separate from changes that alter results.
- State what a change checks and how that was verified. Report failures and
  skipped checks plainly.
- Write commit messages, comments and documentation in plain language a
  maintainer can read quickly.
- If a change alters R's output on the parity fixtures, say so. Regenerate the
  fixtures and bump their version in the same change only when the change is
  intended; otherwise fix the code. See
  [contracts/parity.md](contracts/parity.md).
- Never print or commit a token.

## Repository hygiene

- Do not commit evidence files, logs, saved outputs or large data to a package
  repository. Benchmark results go in the benchmark table in this repository;
  raw logs stay outside git.
- Do not keep superseded code paths, snapshots or prototypes in the package
  source. Remove them, or ask.

## Never

- Push to the Bioconductor staging repository or publish a release.
- Change defaults or statistical methods without a maintainer's decision.
