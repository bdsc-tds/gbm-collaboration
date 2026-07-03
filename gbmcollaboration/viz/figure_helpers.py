from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl

def set_mpl_rcparams(
    params={
        "pdf.fonttype": 42,  # embed TrueType fonts
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "savefig.transparent": True,
    },
):
    mpl.rcParams.update(params)


def savefig(df, out_file, fig=None, dpi=300, bbox_inches="tight",close=False):
    out_file = Path(out_file)
    out_file.parent.mkdir(exist_ok=True, parents=True)

    if df is not None:
        df.to_csv(out_file.with_suffix(".csv"))

    if fig is None:
        fig = plt.gcf()  # get current figure

    fig.savefig(out_file.with_suffix(".pdf"), dpi=dpi, bbox_inches=bbox_inches)
    fig.savefig(out_file.with_suffix(".png"), dpi=dpi, bbox_inches=bbox_inches)

    if close:
        plt.close(fig)
