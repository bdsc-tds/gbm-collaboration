#!/usr/bin/env python3
"""
Assemble Figure 4 (cell-type and ECM neighborhoods) on A4 in the layout of the published figure (manuscript PDF
page 15): boxes read off that page and scaled to A4 by Published; letters at the published positions.
  a, b  fig4ab_suppfig4abg_Cell_Neighborhoods; b = the cell-type map (left half) of Fig4b.pdf, its ECM map is c
  c, d  fig4cdg_suppfig4f_ECM_niches_v4; g = the contact-frequency heatmap (left half) of Fig4g.pdf, as published
  e     fig4e_suppfig4e_colocalization_plots
  f     fig4f_pathologist_annotation_analyses
  h     fig4h_suppfig4cd_neighborhood_survival
Run (from manuscript/): ./symlink_figures.sh figure4, then from the repo root:
  python manuscript/compose/figure4.py
Output: manuscript/composed/figure4.{pdf,png}; compose/order/figure4.order.yml.
"""
from _common import Published, cli, compose

P = Published(15, (0.55, 1.0, 7.75, 10.95))

LAYOUT = [
    ("a", P.box(0.58, 1.06, 5.40, 2.00), "Fig4a.pdf"),
    ("b", P.box(6.00, 1.00, 1.72, 2.10), "part:Fig4b.pdf#0,0,10.1,8.78"),
    ("c", P.box(0.62, 3.42, 2.25, 3.00), "Fig4c.pdf"),
    ("d", P.box(2.95, 3.42, 4.65, 1.60), "Fig4d.pdf"),
    ("e", P.box(2.95, 5.15, 4.60, 2.85), "Fig4e.pdf"),
    ("f", P.box(0.62, 6.95, 2.25, 4.00), "Fig4f.pdf"),
    ("g", P.box(3.10, 8.45, 2.15, 2.50), "part:Fig4g.pdf#0,0,6.85,6.03"),
    ("h", P.box(5.40, 8.45, 2.30, 2.50), "Fig4h.pdf"),
]

if __name__ == "__main__":
    compose("figure4", LAYOUT, cli("figure4"))
