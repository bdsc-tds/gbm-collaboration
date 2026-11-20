"""Survival helpers: survRM2 RMST through rpy2, Kaplan-Meier panels, Cox PH tests.

RMST differences are returned in months; inputs are in years (os_years, os_event).
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


def rmst(df: pd.DataFrame, arm: str, covariates=(), tau: float = 3.0,
         time: str = "os_years", event: str = "os_event") -> dict:
    """survRM2::rmst2 of arm 1 vs arm 0 at tau, adjusted for covariates if given.

    Rows with a missing value in any used column are dropped (complete cases), so code
    covariates beforehand if missing values should count as a level.
    Returns n, n_high (arm 1), n_low, rmst_high, rmst_low, diff_unadj, p_unadj and the
    (adjusted) diff, lo, hi, p, all in months.
    """
    cols = [time, event, arm, *covariates]
    d = df[cols].dropna()
    d = d[d[time] > 0]
    r = d.copy()
    r.columns = ["time", "event", "arm", *[f"c{i}" for i in range(len(covariates))]]
    to_r("sd", r)
    ro.globalenv["tv"] = float(tau)
    ro.r("a <- rmst2(sd$time, sd$event, sd$arm, tau = tv)")
    out = dict(n=len(d), n_high=int(d[arm].sum()), n_low=int((d[arm] == 0).sum()),
               rmst_high=ro.r("a$RMST.arm1$rmst[['Est.']]")[0] * 12,
               rmst_low=ro.r("a$RMST.arm0$rmst[['Est.']]")[0] * 12)
    u = list(ro.r("a$unadjusted.result[1,]"))
    out.update(diff_unadj=u[0] * 12, p_unadj=u[3])
    if covariates:
        cs = ", ".join(f"sd$c{i}" for i in range(len(covariates)))
        a = list(ro.r(f"rmst2(sd$time, sd$event, sd$arm, tau = tv, covariates = cbind({cs}))$adjusted.result[1,]"))
        out.update(diff=a[0] * 12, lo=a[1] * 12, hi=a[2] * 12, p=a[3])
    else:
        out.update(diff=u[0] * 12, lo=u[1] * 12, hi=u[2] * 12, p=u[3])
    return out


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
