#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "== forbidden manuscript/editorial artifacts =="
BAD=$(find . -type f \( -name '*.tex' -o -name '*.bib' -o -iname 'manuscript*.pdf' -o -iname 'paper*.pdf' -o -iname 'cover_letter*' -o -iname 'response_to_referee*' \) -print)
if [[ -n "$BAD" ]]; then
  echo "$BAD"
  exit 1
fi
echo "none"

echo "== Python syntax =="
python -m py_compile code/*.py

echo "== numerical checks =="
python code/check_results.py

echo "== required files =="
for f in README.md LICENSE DATA_LICENSE.md CITATION.cff requirements.txt reproduce_all.sh RELEASE_NOTES_v1.0.0.md ZENODO_RELEASE.md; do
  test -s "$f" || { echo "missing/empty: $f"; exit 1; }
done

echo "PASS: release preflight"
