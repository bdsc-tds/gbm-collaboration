#!/usr/bin/env python3
"""
Assemble Supplementary Figure 1 (survival model diagnostics) on A4, one panel per row (not in the manuscript PDF; letters
from the panel ids). Panels: fig1b_suppfig1_GBM_MS_Summary (symlink_config.yml suppfigure1).
Run (from manuscript/): ./symlink_figures.sh suppfigure1, then from the repo root:
  python manuscript/compose/suppfigure1.py
Output: manuscript/composed/suppfigure1.{pdf,png}; compose/order/suppfigure1.order.yml.
"""
from _common import cli, compose, flow

ROWS = [
    [("a", "SuppFig1a.pdf", 5.4)],
    [("b", "SuppFig1b.pdf", 2.6)],
    [("c", "SuppFig1c.pdf", 2.6)],
]

if __name__ == "__main__":
    args = cli("suppfigure1")
    compose("suppfigure1", flow(ROWS, args.panel_dir), args)
