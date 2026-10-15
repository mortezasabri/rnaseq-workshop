#!/usr/bin/env bash
# Render the presentation: slides/rnaseq_intro.qmd -> docs/rnaseq_intro.html (one self-contained
# file, figures embedded; open it in any web browser, press F for full screen).
# The slides only show figures and tables from results/, so run the analysis first.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
command -v quarto >/dev/null || { echo "ERROR: quarto not found; activate the environment (see README.md)" >&2; exit 1; }
[[ -d results/downstream/3_dea_vis ]] || { echo "ERROR: no results yet; run ./run_all.sh downstream first" >&2; exit 1; }

quarto render slides/rnaseq_intro.qmd
mkdir -p docs
mv slides/rnaseq_intro.html docs/rnaseq_intro.html
rm -rf slides/rnaseq_intro_files

# The slides get shared: refuse to pass a file that contains a local folder path.
if grep -q "$HOME" docs/rnaseq_intro.html; then
  echo "ERROR: docs/rnaseq_intro.html contains the local path '$HOME'; remove it from the slide source." >&2
  exit 2
fi
echo "Slides: docs/rnaseq_intro.html"
