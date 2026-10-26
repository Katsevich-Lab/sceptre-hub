<!--
Release checklist for the sceptre R package (Bioconductor).
Open as an issue titled "Release sceptre x.y.z".
-->

## Prepare

- [ ] All issues in the milestone are closed or moved
- [ ] `git pull`
- [ ] Check the current [Bioconductor build report](https://bioconductor.org/checkResults/)
- [ ] Polish NEWS: user-facing changes only, each accurately scoped
- [ ] `devtools::check()` passes locally
- [ ] `BiocCheck::BiocCheck()`: no errors or warnings; notes reviewed
- [ ] CI passes on Linux, macOS and Windows
- [ ] Vignettes build and examples run
- [ ] Objects saved by the previous release still load
- [ ] Reference stability test passes; if R's output changed on purpose, the parity fixtures were regenerated and their version bumped

## Release

- [ ] Bump the version following [Bioconductor versioning](https://contributions.bioconductor.org/versionnum.html)
- [ ] Push to the Bioconductor repository
- [ ] Confirm the next build report is clean

## Afterwards

- [ ] Tag the release on GitHub
- [ ] Update the pkgdown site
- [ ] Update the compatibility table in the hub contracts
- [ ] Close the milestone
