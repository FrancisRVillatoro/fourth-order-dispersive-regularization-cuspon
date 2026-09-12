#!/usr/bin/env bash
set -euo pipefail

REPO_NAME="fourth-order-dispersive-regularization-cuspon"
OWNER="FrancisRVillatoro"
REMOTE_SSH="git@github.com:${OWNER}/${REPO_NAME}.git"

if [[ -d .git ]]; then
  echo "ERROR: .git already exists. Run this script only in a fresh unpacked release tree." >&2
  exit 2
fi

# Local safety gate: the reproducibility repository must not contain manuscript files.
if find . -type f \( -name '*.tex' -o -name '*.bib' -o -iname 'manuscript*.pdf' -o -iname 'paper*.pdf' -o -iname 'cover_letter*' -o -iname 'response_to_referee*' \) -print -quit | grep -q .; then
  echo "ERROR: manuscript/editorial artifact detected; aborting." >&2
  exit 3
fi

sha256sum -c SHA256SUMS.txt
bash preflight_release.sh

git init -b main
if ! git config user.name >/dev/null; then
  git config user.name "Francisco R. Villatoro"
fi
if ! git config user.email >/dev/null; then
  git config user.email "frvillatoro@uma.es"
fi

git add .
git commit -m "Initial reproducibility release v1.0.0"
git tag -a v1.0.0 -m "v1.0.0: reproducibility package"
git remote add origin "$REMOTE_SSH"

echo
echo "Local Git repository and annotated tag v1.0.0 created."
echo "Create an EMPTY PUBLIC repository on GitHub named: $REPO_NAME"
echo "Then push with:"
echo "  git push -u origin main"
echo "  git push origin v1.0.0"
echo
echo "After pushing, publish a GitHub Release for tag v1.0.0 using RELEASE_NOTES_v1.0.0.md."
