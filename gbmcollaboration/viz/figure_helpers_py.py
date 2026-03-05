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


def savefig(df, out_file, dpi=300, bbox_inches="tight"):
    df.to_csv(Path(out_file).with_suffix(".csv"))
    plt.savefig(Path(out_file).with_suffix(".pdf"), dpi=dpi, bbox_inches=bbox_inches)
    plt.savefig(Path(out_file).with_suffix(".png"), dpi=dpi, bbox_inches=bbox_inches)
    plt.close()
