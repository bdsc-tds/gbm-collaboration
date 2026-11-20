#!/usr/bin/env python3
"""
Assemble Figure 3 (genomic landscape, WES CNV clusters) on A4 in the layout of the published figure (manuscript PDF
page 11): boxes read off that page and scaled to A4 by Published; letters at the published positions.
  a  fig3a_oncoplot
  b  fig3b_suppfig3ac_infercnv_wes_analyses
  c  fig3c_suppfig3b_wes_profile_os_rmst
  d  fig3d_correlation_cnv_wes_ecm
Run (from manuscript/): ./symlink_figures.sh figure3, then from the repo root:
  python manuscript/compose/figure3.py
Output: manuscript/composed/figure3.{pdf,png}; compose/order/figure3.order.yml.
"""
from _common import Published, cli, compose

P = Published(11, (0.62, 1.08, 7.70, 8.80))

LAYOUT = [
    ("a", P.box(0.68, 1.12, 6.97, 3.33), "Fig3a.pdf"),
    ("b", P.box(0.68, 4.52, 6.90, 1.55), "Fig3b.pdf"),
    ("c", P.box(0.68, 6.30, 3.30, 2.50), "Fig3c.pdf"),
    ("d", P.box(4.00, 6.30, 3.60, 2.50), "Fig3d.pdf"),
]

if __name__ == "__main__":
    compose("figure3", LAYOUT, cli("figure3"))
