# First-release checklist

The repository remains at a development version until these items are complete.

- [ ] Revoke or rotate the FRED API key that appeared in repository history.
- [ ] Confirm that current files contain no credentials or restricted data.
- [ ] Resolve the open items in `DATA_PERMISSIONS.md`.
- [ ] Verify attribution and licence notices for every bundled dataset.
- [ ] Run `R CMD check` with no errors or warnings and review all notes.
- [ ] Run the companion Bookdown examples against the release candidate.
- [ ] Set the package version for the first stable release.
- [ ] Create a Git tag and GitHub release from the checked source archive.
- [ ] Consider archiving the release with Zenodo and adding its DOI to the
      package citation metadata.

Generated source archives should be attached to the GitHub release rather than
committed to the repository.

Baseline verification completed on 2026-08-18: the development package passed
`R CMD check --no-manual`, and the public HTML Bookdown build completed all 255
code chunks without `private_data`. Repeat both checks on the final release
candidate.
