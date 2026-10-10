# Bouchet campaign packages as branches

Written 10 Oct 2026. The campaign packages o to ab were built in Claude session
scratchpads as separate git trees and shipped as tarballs in
`C:\Users\Gregory\Downloads`. Packages x to ab were fetched into this repository
on 10 Oct 2026 as branches `pkg/<letter>` (`git fetch <tree> <branch>:pkg/<letter>`),
with their full history. Nothing was merged, rebased or checked out; nothing
was pushed.

## Branches

| Branch | HEAD | Built | Tarball | sha256 |
|---|---|---|---|---|
| pkg/x  | c5140dc453a2d4864cffcf5e5ed1b864790a0deb | 7 Oct 2026, rebuilt 9 Oct 09:28 after the EPSG fix | tcsmewn_v09_20261007x.tar.gz  | eae75709af83e1e652161f070568b2a38ac8d225a893802c2191de2dec454a9e |
| pkg/y  | 95ed01860e591c25f7540b88d5b6d851b25eaa9f | 8 Oct 2026 19:01 | tcsmewn_v09_20261008y.tar.gz  | e926a04ac7d708fefed8e7f05d056ff86f83e32642d789043cda4bf04bd5f526 |
| pkg/z  | 3436327f9faa75bccda21e96ccfcec04ccc4e0a6 | 8 Oct 2026 19:59 | tcsmewn_v09_20261008z.tar.gz  | 053da4c62a2c013575267f5af4471ba5f8ac8cfb1ebfdaecd3624c3333c51e0c |
| pkg/aa | 5329f235390d6967113362d5e019dd6c60f1cff8 | 9 Oct 2026 09:42 | tcsmewn_v09_20261008aa.tar.gz | 05f2303402347cd7dce299e95beeaa94aa064d44ddcd0f5c8040abbc005b313e |
| pkg/ab | 2a958c47fa1a7d77a02a666b7438db8135ff33f2 | 9 Oct 2026 10:18 | tcsmewn_v09_20261009ab.tar.gz | b0fccf1f2e88767e1fd071d8afb53b4d50e4483f6adfcd793ec3acc75c4e1208 |

Dates are the HEAD commit times (America/New_York); each tarball was written
within two minutes of its commit.

## Verification (10 Oct 2026)

- Every tarball's sha256 matched its `.sha256` file in Downloads.
- Every tarball was extracted and diffed (`diff -rq`) against `git archive` of
  the HEAD listed above. All five are identical except for the empty folder
  `v2/Experimental observations/Kantola nitrogen`, which the build script adds
  because git does not track empty directories.
- The tarball READMEs and RUNBOOKs do not record a commit sha; the content diff
  above is the verification.

## Sources

| Branch | Fetched from | Source branch |
|---|---|---|
| pkg/x  | scratchpad `f08b9c7f-…`/pkg_x/tcsmewn_v09_20261007x | master |
| pkg/y  | scratchpad `df917fe5-…`/k/ky | package-y |
| pkg/z  | scratchpad `df917fe5-…`/k/ky | erw-kernels-on-y |
| pkg/aa | scratchpad `22926be7-…`/pkg_z/tree | field-campaign |
| pkg/ab | scratchpad `df917fe5-…`/k/ky | field2d-kernel-gate |

Scratchpad root: `C:\Users\Gregory\AppData\Local\Temp\claude\C--Users-Gregory-OneDrive-Yale-TC-SMEW-N-T-C-SMEW-N-v0-9\<session>\scratchpad`.

## History

- The package history (44 commits on pkg/ab) is rooted at the orphan commit
  6a4e54f "Package w as delivered (78f65a22...), base for package x". It shares
  no commit with `main` or `bouchet-campaign-gating`; package w was imported as
  files, not as history.
- pkg/y, pkg/z, pkg/aa and pkg/ab form one line: y < z < aa < ab.
- pkg/x (c5140dc) is not on that line. x and y both branch from 869116c
  "Pins for the package-x build with the Field_2D options"; c5140dc is x's
  EPSG fix in Load_Field_DEM (strfind), while y carries its own fix at 6e12e10.
- Package w (tcsmewn_v09_20261005w.tar.gz, sha256 78f65a22f2820e95bb0802d3fb1d015a29b896070f11abb52d2c86c5445788c1)
  and earlier packages have no git tree; their content is the root commit's tree.

## Package ac (pending)

Package ac (revised uncertainty members) was being built by session
local_72650a3c on 10 Oct 2026 in scratchpad `17de9852-…`/ac/tree, branch
`uncertainty-bounds-ac`, at the time of writing still at HEAD 2a958c4 (ab) with
uncommitted edits (README, RUNBOOK, Uncertainty_Members.m,
test_uncertainty_campaign.m, the two Uncertainty READMEs, a new
UNCERTAINTY_BOUNDS_REVISION_2026-10-09.md and precheck_1d.m). No tarball and no
hand-over note existed. Not fetched. When it is committed and shipped, fetch it
the same way as `pkg/ac` and add its row above.
