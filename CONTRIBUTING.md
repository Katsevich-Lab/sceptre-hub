# How we work

This guide covers how work moves from an idea to a merged change across the
SCEPTRE repositories. Each repository's own contribution guide covers local
setup, style and testing.

## Where things go

| What | Where |
| --- | --- |
| A decision that spans languages: new methods, changed defaults, new result columns, public functions in more than one package | A proposal in [proposals/](proposals/) |
| A design decision within one package, such as a new object layout | A short design note in that repository, reviewed before code |
| A bug, task or feature in one package | An issue in that repository |
| Ranked priorities | [ROADMAP.md](ROADMAP.md) |
| Specs both packages follow | [specs/](specs/) |

## Issues, pull requests and milestones

- **Issue:** one problem or piece of work. Aim for one issue per pull request.
  Larger work becomes a parent issue with sub-issues, each closed by one pull
  request.
- **Pull request:** links its issue in the description ("Closes #12", or
  "Part of #12" for an intermediate step). It states what it checks and how
  that was verified, following the
  [pull-request template](templates/pull_request_template.md).
- **Milestone:** a release of one repository, with a target date. Assign
  milestones to issues rather than pull requests.
- **Project board:** one GitHub Project tracks issues and pull requests from all
  repositories. Its fields are Status, Focus area (the roadmap priority the item
  serves), Repository and Iteration.

Closing keywords such as "Closes #12" act only on merges into a repository's
default branch. Issues fixed by pull requests merged into a development branch
are closed by hand, or when that branch is merged.

## Review

Match the depth of review to the risk of the change:

| Kind of change | Review |
| --- | --- |
| Restructuring that should not change behavior | Can be one larger pull request if tests written against the base branch beforehand show unchanged results. Review the design and the tests, and skim the moved code. |
| Bug fixes and anything that changes results | Small, separate pull requests, each read in full, with a regression test and a changelog entry saying who is affected |
| New public functions, dependencies, object formats or result columns | A reviewed design note or accepted proposal first, then a modest pull request |
| Performance | Results unchanged on the tests, plus one before/after timing table |

Methodological changes need scientific evidence. Passing implementation tests
alone does not establish calibration, power or biological validity.

## Keeping work flowing

- **Limit open pull requests to one or two per repository.** Finish and merge
  before starting the next. This keeps review manageable and keeps work landing.
- **Progress is counted in merged pull requests.**
- **Two-week cycles.** At the end of each cycle: what merged, what is next, what
  is waiting on review.
- **Weekly check-in** on the project board and open reviews.

## Branches

Each repository documents its own branching. In general, `main` holds the
released or submitted version, and a `dev` branch collects reviewed work for the
next release. Short-lived feature branches come off `dev`, one per pull request,
and are deleted after merging.

The R package is under Bioconductor review. Only maintainers push to the
Bioconductor staging repository, deliberately, with a version bump.

## Releases

Each release is tracked by a release-checklist issue made from the template for
[Bioconductor](templates/release-checklist-bioconductor.md) or
[PyPI](templates/release-checklist-pypi.md). Packages release on their own
schedules. The specs record which package versions support which version of the
result schema.

## Research boundaries

Do not copy unpublished research, private issue contents, manuscript
discussions or internal notes into public issues or this repository without the
relevant collaborators' agreement. The roadmap can describe a capability under
development while the specific approach stays private until it is published.

## AI-assisted contributions

AI-assisted work is welcome and follows the same rules as any other
contribution, plus those in [AGENTS.md](AGENTS.md). Pull requests disclose
whether and how AI was used.
