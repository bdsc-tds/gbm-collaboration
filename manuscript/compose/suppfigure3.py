#!/usr/bin/env python3
"""
Assemble Supplementary Figure 3 (WES CNV clusters) on A4 (not in the manuscript PDF; letters from the panel ids);
e and f are the two halves of SuppFig3e-f.pdf.
Panels: fig3b_suppfig3ac_infercnv_wes_analyses, fig3c_suppfig3b_wes_profile_os_rmst, suppfig3d_wes_cluster_covariate_rmst,
suppfig3ef_analyses_wes_sc_infercnv (symlink_config.yml suppfigure3).
Run (from manuscript/): ./symlink_figures.sh suppfigure3, then from the repo root:
  python manuscript/compose/suppfigure3.py
Output: manuscript/composed/suppfigure3.{pdf,png}; compose/order/suppfigure3.order.yml.
"""
from _common import cli, compose, flow

ROWS = [
    [("a", "SuppFig3a.pdf", 2.8)],
    [("b", "SuppFig3b.pdf", 2.6)],
    [("c", "SuppFig3c.pdf"), ("d", "SuppFig3d.pdf")],
    [("e", "part:SuppFig3e-f.pdf#0,0,6.3,5.84"), ("f", "part:SuppFig3e-f.pdf#6.3,0,6.26,5.84")],
]

if __name__ == "__main__":
    args = cli("suppfigure3")
    compose("suppfigure3", flow(ROWS, args.panel_dir), args)
