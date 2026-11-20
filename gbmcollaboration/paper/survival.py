"""Survival helpers through rpy2 (loads R survival and survRM2): Kaplan-Meier panels, Cox PH tests.

Inputs are in years (os_years, os_event).
"""

import numpy as np
import pandas as pd

from . import config  # noqa: F401  (sets R_HOME for rpy2)
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import pandas2ri  # noqa: E402
from rpy2.robjects.conversion import localconverter  # noqa: E402

ro.r("suppressPackageStartupMessages({library(survival); library(survRM2)})")


def to_r(name: str, df: pd.DataFrame) -> None:
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv[name] = ro.conversion.py2rpy(df.reset_index(drop=True))


def cox_zph(df: pd.DataFrame, covariates, time: str = "os_years", event: str = "os_event") -> pd.DataFrame:
    """Cox model (R survival) and cox.zph test of the scaled Schoenfeld residuals."""
    d = df[[time, event, *covariates]].dropna()
    to_r("cd", d.rename(columns={time: "time", event: "event"}))
    f = " + ".join(f"`{c}`" for c in covariates)
    ro.r(f"fit <- coxph(Surv(time, event) ~ {f}, data = cd); z <- cox.zph(fit)")
    tab = np.array(ro.r("z$table"))
    rows = [str(x) for x in ro.r("rownames(z$table)")]
    out = pd.DataFrame(tab, index=rows, columns=["chisq", "df", "p"])
    out.attrs["n"] = len(d)
    return out


def km(ax, df: pd.DataFrame, arm: str, labels: dict, colors: dict, tau: float | None = 3.0,
       time: str = "os_years", event: str = "os_event", ci_alpha: float = 0.15):
    """Kaplan-Meier curves of arm 1 and arm 0 on ax, with n in the legend."""
    from lifelines import KaplanMeierFitter

    for v in (1, 0):
        g = df[df[arm] == v]
        KaplanMeierFitter().fit(g[time], g[event], label=f"{labels[v]} (n={len(g)})").plot_survival_function(
            ax=ax, color=colors[v], ci_alpha=ci_alpha, lw=1.6)
    if tau is not None:
        ax.axvline(tau, ls="--", c="grey", lw=1)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, None)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
