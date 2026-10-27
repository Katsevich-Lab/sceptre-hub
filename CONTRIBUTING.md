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
| Contracts the implementations must agree on | [contracts/](contracts/) |

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

Maintainer attention is the scarce resource, so machines do most of the
checking and a person reads the evidence. A bug that reaches a release cannot be
pulled back quickly, and wrong results rarely raise an error, so the checking
happens before merging.

### Stages

Every pull request goes through these stages in order:

1. **Author verification.** The author, person or agent, runs the change and
   shows the output: each fix's regression test fails on the base branch and
   passes on the pull request, plus the full test suite, the package check and
   the formatter.
2. **Agent review.** A reviewer with fresh context, given the diff and the issue
   but not the author's own summary. Prefer a standard multi-agent reviewer that
   verifies each finding before reporting it (for example `/code-review` in
   Claude Code) over a single hand-prompted agent.
3. **Regression check.**
   - Every pull request: the frozen behavior suite and the parity fixtures, in
     CI.
   - Changes to import or analysis code: the base branch and the pull request
     run side by side on small real datasets, comparing every result, run time
     and peak memory. Run at full scale only when the change targets
     large-scale performance (speed or memory on big data).
   - Before each release: full-scale analyses of public datasets on the
     development version, including negative-control calibration checks.
4. **Maintainer approval.** See [Approval](#approval).

If a later stage leads to a code change, the earlier stages no longer describe
the code. Rerun them in proportion to the change: a wording fix needs nothing,
a change to import logic needs review and the regression check again.

Design questions, such as whether a default should change or input that used to
be accepted should now be rejected, go to the approver before the work starts,
not at review.

### Depth by risk

Match the depth of review to the risk of the change:

| Kind of change | Review |
| --- | --- |
| Restructuring that should not change behavior | Can be one larger pull request if tests written against the base branch beforehand show unchanged results. Review the design and the tests, and skim the moved code. |
| Bug fixes and anything that changes results | Small, separate pull requests, each read in full, with a regression test and a changelog entry saying who is affected |
| New public functions, dependencies, object formats or result columns | A reviewed design note or accepted proposal first, then a modest pull request |
| Performance | Results unchanged on the tests, plus one before/after timing table |

Methodological changes need scientific evidence. Passing implementation tests
alone does not establish calibration, power or biological validity.

### Approval

Each package has one approver:

| Repository | Approver |
| --- | --- |
| sceptre (R) | Eugene Katsevich ([@ekatsevi](https://github.com/ekatsevi)) |
| pysceptre | Eugenio Mattei ([@emattei](https://github.com/emattei)) |

- The approver reads the evidence: the commands run and their output, the
  reviewer's findings, the regression report. If it is unclear what the checks
  would miss, ask for a fuller review rather than approving.
- The approver may be the person who directed the agent that wrote the change.
- An agent's approval never counts as the approver's.
- An approval covers one commit. New commits need approval again.
- Every review finding is answered before merging: fixed, or resolved with a
  short reason.

### Review rounds

- At most three rounds of review and fixes. After that, step back and
  reconsider the approach.
- After the first round, review only the new commits.
- Name the class of a bug, so the author fixes every instance, not one line.
- Do not repeat findings already raised.

### Fixing findings at the source

When reviewers keep raising the same finding, fix it where it comes from: a
formatter or lint check, a test, or a rule in the repository's `REVIEW.md`. Each
package repository keeps a short `REVIEW.md` saying what failure looks like
there, such as results that are wrong without an error, results that change
without the changelog saying who is affected, or memory and run time that blow
up at scale. Bugs found in the past are a test set for the reviewer: when its
setup changes, check that it still finds them.

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
schedules. The [contracts](contracts/) record which package versions support
which version of the parity fixtures and the result schema.

## Research boundaries

Do not copy unpublished research, private issue contents, manuscript
discussions or internal notes into public issues or this repository without the
relevant collaborators' agreement. The roadmap can describe a capability under
development while the specific approach stays private until it is published.

## AI-assisted contributions

AI-assisted work is welcome and follows the same rules as any other
contribution, plus those in [AGENTS.md](AGENTS.md). Pull requests disclose
whether and how AI was used.
