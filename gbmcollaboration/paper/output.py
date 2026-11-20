"""Output folder, figure saving, plot style and package versions for the paper scripts."""

import importlib.metadata as md
import logging
import re
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from gbmcollaboration.viz import savefig, set_mpl_rcparams

from .config import OUT_ROOT, REPO_ROOT


def setup(script_name: str) -> Path:
    """Common notebook settings; returns the script's output folder OUT_ROOT/<script_name>."""
    warnings.filterwarnings("ignore", category=FutureWarning)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    set_mpl_rcparams()
    out = OUT_ROOT / script_name
    out.mkdir(parents=True, exist_ok=True)
    return out


PANELS_DIR = OUT_ROOT / "panels"


def panel_name(panel: str) -> str:
    """'Fig 4h' -> 'Fig4h', 'Supp Fig 5a' -> 'SuppFig5a', 'Fig 6b hBrain' -> 'Fig6b_hBrain'."""
    m = re.match(r"^\s*((?:Supp\s+)?(?:Fig|Table)\s*\d+[a-z]?(?:-[a-z])?)\s*(.*)$", panel)
    if not m:
        return "_".join(panel.split())
    head, rest = "".join(m.group(1).split()), "_".join(m.group(2).split())
    return f"{head}_{rest}" if rest else head


def save(out_dir: Path, stem: str, fig=None, data: pd.DataFrame | None = None, panel: str | None = None) -> None:
    """Write <stem>.pdf/.png (and <stem>.csv with the plotted data) to out_dir.

    panel: the manuscript panel this figure is (e.g. "Fig 4h", "Supp Fig 5a"). The figure is
    then also written to OUT_ROOT/panels/<Fig4h>.pdf/.png/.csv, the one folder holding every
    panel of the paper under its manuscript name. Use a suffix for several files of one panel
    ("Fig 6b hBrain").
    """
    savefig(df=data, out_file=Path(out_dir) / stem, fig=fig)
    if panel is not None:
        PANELS_DIR.mkdir(parents=True, exist_ok=True)
        savefig(df=data, out_file=PANELS_DIR / panel_name(panel), fig=fig)


def save_table(out_dir: Path, stem: str, df: pd.DataFrame, panel: str | None = None, **kw) -> None:
    """Write <stem>.csv to out_dir; with panel (e.g. "Supp Table 3b") also to OUT_ROOT/panels/."""
    df.to_csv(Path(out_dir) / f"{stem}.csv", **kw)
    if panel is not None:
        PANELS_DIR.mkdir(parents=True, exist_ok=True)
        df.to_csv(PANELS_DIR / f"{panel_name(panel)}.csv", **kw)




@contextmanager
def scanpy_style(dpi: int = 200):
    """scanpy figure defaults (axes grid, 14 pt font) for one block, as the notebooks that
    called sc.set_figure_params at startup drew their panels; restores settings afterwards."""
    import scanpy as sc

    with plt.rc_context():
        sc.set_figure_params(dpi=dpi, color_map="viridis_r")
        yield


def versions(py_packages=(), r_packages=()) -> pd.DataFrame:
    """Python and R package versions this run used (D17)."""
    def v(p):
        try:
            return md.version(p)
        except md.PackageNotFoundError:
            if p == "gbmcollaboration":  # used from the repo, not installed: report the commit
                import subprocess
                rev = subprocess.run(["git", "-C", str(REPO_ROOT), "describe", "--always", "--dirty"],
                                     capture_output=True, text=True).stdout.strip()
                return f"repo {rev}" if rev else "repo"
            return "not installed"

    rows = [("python", sys.version.split()[0])] + [(p, v(p)) for p in ("gbmcollaboration", *py_packages)]
    if r_packages:
        import rpy2.robjects as ro

        rows.append(("R", ro.r("paste(R.version$major, R.version$minor, sep = '.')")[0]))
        rows += [(f"R {p}", ro.r(f'as.character(packageVersion("{p}"))')[0]) for p in r_packages]
    return pd.DataFrame(rows, columns=["package", "version"])
