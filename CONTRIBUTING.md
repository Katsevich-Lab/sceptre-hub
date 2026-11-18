# How we work

This guide covers what belongs where across the SCEPTRE repositories. Each
repository's own contribution guide covers local setup, style, testing and how
changes are reviewed there.

## Where things go

| What | Where |
| --- | --- |
| A decision that spans languages: new methods, changed defaults, new result columns, public functions in more than one package | Decided by the maintainers and recorded in this repository |
| A design decision within one package, such as a new object layout | A short design note in that repository, reviewed before code |
| A bug, task or feature in one package | An issue in that repository |
| Ranked priorities | [ROADMAP.md](ROADMAP.md) |
| Contracts the implementations must agree on | [contracts/](contracts/) |

## Who decides

| Repository | Maintainer |
| --- | --- |
| sceptre (R) | Eugene Katsevich ([@ekatsevi](https://github.com/ekatsevi)) |
| pysceptre | Eugenio Mattei ([@emattei](https://github.com/emattei)) |

Methodological changes need scientific evidence. Passing implementation tests
alone does not establish calibration, power or biological validity.

## Branches

Each repository documents its own branching.

The R package is under Bioconductor review. Only maintainers push to the
Bioconductor staging repository, deliberately, with a version bump.

## Releases

Each release is tracked by a release-checklist issue made from the template for
[Bioconductor](templates/release-checklist-bioconductor.md) or
[PyPI](templates/release-checklist-pypi.md). Packages release on their own
schedules. The [contracts](contracts/) record which package versions support
which version of the parity fixtures and the result schema.

## Research boundaries

Do not copy unpublished research, private issue contents, manuscript
discussions or internal notes into public issues or this repository without the
relevant collaborators' agreement. The roadmap can describe a capability under
development while the specific approach stays private until it is published.

## AI-assisted contributions

AI-assisted work is welcome and follows the same rules as any other
contribution, plus those in [AGENTS.md](AGENTS.md).
