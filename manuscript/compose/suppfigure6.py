#!/usr/bin/env python3
"""
Assemble Supplementary Figure 6 (NicheNet circos of the other TIL zones) on A4, one panel per row (not in the manuscript
PDF; letters from the panel ids). Panels: fig6a_suppfig6_NicheNet_circos (symlink_config.yml suppfigure6).
Run (from manuscript/): ./symlink_figures.sh suppfigure6, then from the repo root:
  python manuscript/compose/suppfigure6.py
Output: manuscript/composed/suppfigure6.{pdf,png}; compose/order/suppfigure6.order.yml.
"""
from _common import cli, compose, flow

ROWS = [
    [("a", "SuppFig6a.pdf", 5.45)],
    [("b", "SuppFig6b.pdf", 5.45)],
]

if __name__ == "__main__":
    args = cli("suppfigure6")
    compose("suppfigure6", flow(ROWS, args.panel_dir), args)
