# Git / WSL release steps

The intended public repository is:

`https://github.com/FrancisRVillatoro/fourth-order-dispersive-regularization-cuspon`

## 1. Unpack in the WSL filesystem

Prefer the Linux filesystem (`~/...`) rather than `/mnt/c/...` for Git metadata and Python environments.

```bash
cd ~
unzip /mnt/c/Users/FRVillatoro/Downloads/fourth-order-dispersive-regularization-cuspon_v1.0.0.zip
cd fourth-order-dispersive-regularization-cuspon
```

Adjust the Windows download path if necessary.

## 2. Verify the frozen release tree

Do **not** rerun the numerical scripts before the first Git commit: regenerated PDFs/logs can differ byte-for-byte even when the numerical results are identical, which would intentionally make the frozen checksums fail. Verify the prepared release instead:

```bash
sha256sum -c SHA256SUMS.txt
bash preflight_release.sh
```

If you want to perform an independent clean scientific rerun, do it after the `v1.0.0` commit/tag (or in a separate clone) with `bash reproduce_all.sh`.

## 3. Create the empty GitHub repository

On GitHub create a **PUBLIC**, empty repository named

`fourth-order-dispersive-regularization-cuspon`

Do not ask GitHub to add a README, `.gitignore`, or license because all three are already present here.

## 4. Initialize local Git and tag v1.0.0

```bash
bash git_init_v1.0.0_wsl.sh
```

The script verifies hashes and numerical checks, initializes `main`, creates the first commit, creates the annotated tag `v1.0.0`, and configures the SSH remote.

Inspect:

```bash
git status
git log --oneline --decorate --graph -5
git tag -n
```

## 5. Push

```bash
git push -u origin main
git push origin v1.0.0
```

If SSH is not configured, replace the remote with HTTPS:

```bash
git remote set-url origin https://github.com/FrancisRVillatoro/fourth-order-dispersive-regularization-cuspon.git
git push -u origin main
git push origin v1.0.0
```

## 6. Publish the GitHub Release

Create a release for the existing tag `v1.0.0`.

- Title: `v1.0.0 — reproducibility release`
- Body: paste `RELEASE_NOTES_v1.0.0.md`
- This is a full release, not a prerelease.

## 7. Zenodo

Before publishing the GitHub release, enable this repository in the Zenodo GitHub integration. Zenodo will use `CITATION.cff` for metadata. After the GitHub release is published and Zenodo finishes ingestion, record the minted **version DOI** and **concept DOI**.

Follow `ZENODO_RELEASE.md` for the post-DOI updates.
