#!/bin/bash
# Link every panel (symlink_figures.sh), compose every figure into composed/ and copy the composed PDFs / PNGs to EFS.
#   ./build_all.sh                                  panels from symlink_config.yml BASE_DIR
#   BASE_DIR=<panel dir> ./build_all.sh             another panel folder (e.g. figures/paper/_panels_before_rerun/)
#   PYTHON=<python> COMPOSED_EFS=<dir> ./build_all.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-$HOME/.conda/envs/gbm_main/bin/python}"
COMPOSED_EFS="${COMPOSED_EFS:-/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/paper/composed}"

mkdir -p "$SCRIPT_DIR/figures"
"$SCRIPT_DIR/symlink_figures.sh" --clean > "$SCRIPT_DIR/figures/symlink_figures.log"
grep -a -A1 "not found" "$SCRIPT_DIR/figures/symlink_figures.log" || true

for fig in figure1 figure2 figure3 figure4 figure5 figure6 \
           suppfigure1 suppfigure2 suppfigure3 suppfigure4 suppfigure5 suppfigure6; do
    "$PYTHON" "$SCRIPT_DIR/compose/$fig.py" "$@"
done

mkdir -p "$COMPOSED_EFS"
cp "$SCRIPT_DIR"/composed/*.pdf "$SCRIPT_DIR"/composed/*.png "$COMPOSED_EFS"/
echo "Copied composed/*.pdf, *.png to $COMPOSED_EFS"
