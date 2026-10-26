<!--
Release checklist for pysceptre (PyPI).
Open as an issue titled "Release pysceptre x.y.z".
-->

## Prepare

- [ ] All issues in the milestone are closed or moved
- [ ] `git pull`
- [ ] Polish CHANGELOG: user-facing changes only, each accurately scoped
- [ ] Tests pass on CPU in CI
- [ ] GPU tests pass on a GPU machine; record which one
- [ ] Documentation builds
- [ ] Parity comparison tests pass against the latest R fixture version (sync first if R has released a new one)

## Release

- [ ] Bump the version
- [ ] Build the sdist and wheels; check them with `twine check`
- [ ] Publish to TestPyPI and install from it in a clean environment
- [ ] Publish to PyPI

## Afterwards

- [ ] Tag the release on GitHub
- [ ] Publish the documentation
- [ ] Update the compatibility table in the hub contracts
- [ ] Close the milestone
