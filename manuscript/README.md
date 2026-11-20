# manuscript/

Composed full-figure pages of Watson et al. (MOSAIC GBM), built from the individual panels that `analysis/paper/`
writes to EFS `figures/paper/panels/` (one `<Id>.pdf/.png/.csv` per manuscript panel, see `analysis/paper/README.md`).
Same structure as the `manuscript/` folders of bdsc-tds/mesothelioma and bdsc-tds/he-predictability.

- `symlink_config.yml` + `symlink_figures.sh [figureN]`: link panels into `figures/<figureN>/`, their CSVs into
  `figures/<figureN>/data/` (`--nocsv` per line; CSV-only entries go straight to `data/`). `BASE_DIR` in the config
  (EFS `figures/paper/panels/`); `BASE_DIR=<dir> ./symlink_figures.sh` uses another panel folder, e.g.
  `figures/paper/_panels_before_rerun/`. `tables` and `response` are linked only (CSVs of the tables, response-to-reviewer
  panels), not composed.
- `compose/figureN.py`, `compose/suppfigureN.py`: place the linked panels on A4 -> `composed/<fig>.pdf` and `.png`
  (`<fig>_p<N>.png` per page for a multi-page figure), and write `compose/order/<fig>.order.yml` (letter -> panels, from
  the layout; do not edit by hand).
- `compose/_common.py`: shared by the compose scripts (`linked()`, `size()`, `place()`, `place_crop()`, placeholders,
  letter overlay, `write_order()`, `grid()`, `flow()`, `Published`, `compose()` for box layouts, `check_inside()`,
  `check_no_overlap()`). Panels are scaled down, never up, to fit their box; `--final` fails on placeholders, missing
  panels, crops or text below 4 pt.
- `build_all.sh`: `symlink_figures.sh --clean`, every compose script, then copies `composed/*.pdf, *.png` to EFS
  `figures/paper/composed/` (`BASE_DIR`, `PYTHON`, `COMPOSED_EFS` override the defaults; extra arguments go to every
  compose script).

```bash
cd manuscript && ./build_all.sh                       # all figures
BASE_DIR=<panel dir> ./build_all.sh                    # from another panel folder
./symlink_figures.sh figure4 && python compose/figure4.py   # one figure
```

## Main figures (1-6): the published layout

The published figures are raster images on pages 5, 7, 11, 15, 20 and 23 of the manuscript PDF
(`Watson et al MOSAIC 2026 v4.docx.pdf` at the repo root, not tracked; `MANUSCRIPT_PDF=<pdf>` for another copy). Each
`figureN.py` gives its boxes in inches on that published page (top-left origin, letter at the box's top-left);
`Published(page, content)` scales them to the A4 width, so positions and proportions follow the published figure.
Items of a layout:

| item | placed |
|---|---|
| `Fig4h.pdf` | the panel, vector, scaled to fit, centred, top-aligned |
| `part:Fig4b.pdf#x,y,w,h` | that region (inches) of our panel, vector (e.g. one half of a two-part panel) |
| `crop:<pdf>@<page>#x,y,w,h[!x,y,w,h]` (`Published.crop()`) | that region of the published page, for panels not produced by code; the published letter (clip corner, or the region after `!`) is painted white and redrawn |
| `pending: <text>` | grey placeholder |

A letter `c*` keeps the published letter inside its crop (white letters on the Fig 6 microscopy). `grid()` places a
panel made of several files (Fig 1b) as a grid inside one box.

## Supplementary figures (1-6)

Not in the manuscript PDF. `flow()` lays out rows of panels at the width of the page (letters from the panel ids), a new
page when a row does not fit or at a `None` row: Supp Fig 2 has six pages, Supp Fig 4 two.
