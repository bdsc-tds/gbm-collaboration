"""Clinical survival analyses of GBM_MS_Summary (Spencer Watson) used by more than one run.

ph_diagnostics: the "Check Cox HR assumptions" cell (Supp Fig 1a-1b), moved here verbatim
as a function so that it runs on the published table (every subject of the clinical
tables, N = 287) and on the 260-patient survival cohort without duplicating the cell.
Changes from the cell: clin_out -> clin_table argument, SAVE_DIR -> save_dir argument,
`suffix` added to the file names, figures saved with output.save (panel copies when
panels=True), and the Schoenfeld p-values are returned.
"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

from . import config  # noqa: F401  (sets R_HOME for rpy2)
from . import output
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import pandas2ri  # noqa: E402
from rpy2.robjects.conversion import localconverter  # noqa: E402

ro.r("suppressPackageStartupMessages({library(survival); library(survRM2)})")


def ph_diagnostics(clin_table: pd.DataFrame, save_dir, suffix: str = "", panels: bool = False) -> dict:
    """Cox model of MGMT + age + sex, cox.zph on the scaled Schoenfeld residuals, PH
    diagnostic figure (Supp Fig 1a) and tau-selection figure (Supp Fig 1b) for clin_table.

    One row per "Subject Id" (baseline row first); a missing MGMT or sex counts as
    unmethylated or male, as in the cell. Returns n, events and the cox.zph p-values.
    """
    # =============================================================================
    # CONFIG
    # =============================================================================

    TIME_COL  = "os_years"
    EVENT_COL = "os_event"
    MGMT_COL  = "MGMT promoter methylation"
    AGE_COL   = "age_at_diagnosis"
    SEX_COL   = "Administrative gender"

    COL_METH   = "#ff7b02"
    COL_UNMETH = "#007ebc"

    SAVE_DIR = save_dir

    # =============================================================================
    # PREPARE DATA
    # =============================================================================

    clin_ph = clin_table.copy()
    clin_ph["_priority"] = clin_ph["Sample collection chronology"].map(
        {"Baseline": 0, "At recurrence": 1}).fillna(2)
    clin_ph = clin_ph.sort_values(["Subject Id", "_priority"])
    clin_ph = clin_ph.drop_duplicates(subset=["Subject Id"], keep="first")
    clin_ph = clin_ph.drop(columns=["_priority"])

    clin_ph["time"]  = pd.to_numeric(clin_ph[TIME_COL],  errors="coerce")
    clin_ph["event"] = pd.to_numeric(clin_ph[EVENT_COL], errors="coerce")
    clin_ph["mgmt"]  = (clin_ph[MGMT_COL].astype(str) == "Methylated").astype(int)
    clin_ph["sex"]   = (clin_ph[SEX_COL].astype(str) == "Female").astype(int)
    age_vals         = pd.to_numeric(clin_ph[AGE_COL], errors="coerce")
    clin_ph["age_z"] = (age_vals - age_vals.mean()) / age_vals.std()

    df = clin_ph[["time", "event", "mgmt", "age_z", "sex"]].dropna()
    df = df[df["time"] > 0].copy()

    n_meth   = int(df["mgmt"].sum())
    n_unmeth = int((df["mgmt"] == 0).sum())
    median_os = float(df[df["event"] == 1]["time"].median())

    print(f"N = {len(df)},  events = {int(df['event'].sum())}")
    print(f"MGMT methylated: {n_meth} / {len(df)}")
    print(f"Median OS (events): {median_os:.2f} years")

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv["surv_data"] = ro.conversion.py2rpy(df)

    # =============================================================================
    # COX + SCHOENFELD (already run, re-run cleanly here)
    # =============================================================================

    ro.r("""
    cox_fit  <- coxph(Surv(time, event) ~ mgmt + age_z + sex,
                      data = surv_data, ties = "efron")
    zph_test <- cox.zph(cox_fit)
    zph_id   <- cox.zph(cox_fit, transform = "identity")

    # KM curves
    km_meth   <- survfit(Surv(time, event) ~ 1,
                         data = surv_data[surv_data$mgmt == 1,])
    km_unmeth <- survfit(Surv(time, event) ~ 1,
                         data = surv_data[surv_data$mgmt == 0,])
    km_all    <- survfit(Surv(time, event) ~ 1, data = surv_data)

    t_meth   <- km_meth$time;   s_meth   <- km_meth$surv
    t_unmeth <- km_unmeth$time; s_unmeth <- km_unmeth$surv
    t_all    <- km_all$time
    h_all    <- km_all$n.event / (km_all$n.risk *
                 c(diff(c(0, km_all$time))))

    # Schoenfeld residuals
    scho_t    <- zph_id$time
    scho_mgmt <- zph_id$y[, "mgmt"]
    scho_age  <- zph_id$y[, "age_z"]
    beta_mgmt <- zph_id$y[, "mgmt"]

    lo_mgmt  <- lowess(scho_t, scho_mgmt, f = 0.5)
    lo_age   <- lowess(scho_t, scho_age,  f = 0.5)
    lo_beta  <- lowess(scho_t, beta_mgmt, f = 0.5)

    ref_beta_mgmt <- coef(cox_fit)["mgmt"]
    ref_beta_age  <- coef(cox_fit)["age_z"]

    # O/E ratio — compute manually per group
    # Methylated group
    km_log_meth   <- survfit(Surv(time, event) ~ 1,
                              data = surv_data[surv_data$mgmt == 1,])
    km_log_unmeth <- survfit(Surv(time, event) ~ 1,
                              data = surv_data[surv_data$mgmt == 0,])

    # Cumulative observed events
    cum_obs_meth   <- cumsum(km_log_meth$n.event)
    cum_obs_unmeth <- cumsum(km_log_unmeth$n.event)

    # Expected under PH: scale by overall event rate
    total_events <- sum(surv_data$event)
    total_n      <- nrow(surv_data)
    prop_meth    <- sum(surv_data$mgmt) / total_n
    prop_unmeth  <- 1 - prop_meth

    # Cumulative expected (simple proportional)
    cum_exp_meth   <- (1:length(cum_obs_meth))   * prop_meth   * 0 +
                       cumsum(km_log_meth$n.event   / km_log_meth$n.risk   *
                              km_log_meth$n.risk)
    cum_exp_unmeth <- cumsum(km_log_unmeth$n.event / km_log_unmeth$n.risk *
                              km_log_unmeth$n.risk)

    # Use log-rank expected counts approach via survdiff
    sd_mgmt <- survdiff(Surv(time, event) ~ mgmt, data = surv_data)
    """)

    mgmt_zph_p   = float(ro.r("zph_test$table['mgmt',   'p']")[0])
    age_zph_p    = float(ro.r("zph_test$table['age_z',  'p']")[0])
    sex_zph_p    = float(ro.r("zph_test$table['sex',    'p']")[0])
    global_zph_p = float(ro.r("zph_test$table['GLOBAL', 'p']")[0])

    t_meth   = np.array(ro.r("t_meth"));   s_meth   = np.array(ro.r("s_meth"))
    t_unmeth = np.array(ro.r("t_unmeth")); s_unmeth = np.array(ro.r("s_unmeth"))
    t_all    = np.array(ro.r("t_all"));    h_all    = np.array(ro.r("h_all"))
    scho_t   = np.array(ro.r("scho_t"))
    scho_mgmt = np.array(ro.r("scho_mgmt"))
    scho_age  = np.array(ro.r("scho_age"))
    beta_mgmt = np.array(ro.r("beta_mgmt"))
    lo_sx_m  = np.array(ro.r("lo_mgmt$x")); lo_sy_m = np.array(ro.r("lo_mgmt$y"))
    lo_sx_a  = np.array(ro.r("lo_age$x"));  lo_sy_a = np.array(ro.r("lo_age$y"))
    lo_tx    = np.array(ro.r("lo_beta$x")); lo_ty   = np.array(ro.r("lo_beta$y"))
    ref_hr_mgmt = float(np.exp(ro.r("ref_beta_mgmt")[0]))
    ref_hr_age  = float(np.exp(ro.r("ref_beta_age")[0]))

    # O/E via survdiff — extract observed and expected
    ro.r("""
    oe_obs_vals <- sd_mgmt$obs
    oe_exp_vals <- sd_mgmt$exp
    oe_labels   <- names(sd_mgmt$obs)
    """)
    oe_obs_vals = np.array(ro.r("oe_obs_vals"))
    oe_exp_vals = np.array(ro.r("oe_exp_vals"))

    print(f"\nOverall O/E (MGMT):")
    print(f"  Methylated:   O={oe_obs_vals[1]:.0f}, E={oe_exp_vals[1]:.1f}, O/E={oe_obs_vals[1]/oe_exp_vals[1]:.2f}")
    print(f"  Unmethylated: O={oe_obs_vals[0]:.0f}, E={oe_exp_vals[0]:.1f}, O/E={oe_obs_vals[0]/oe_exp_vals[0]:.2f}")

    # For time-varying O/E, compute from KM directly
    ro.r("""
    # Compute time-varying OE using counting process approach
    # Split follow-up into intervals and compute O/E per interval
    breaks  <- quantile(surv_data$time[surv_data$event==1],
                        probs = seq(0, 1, by = 0.1))
    breaks  <- unique(c(0, breaks))
    n_int   <- length(breaks) - 1

    oe_time  <- numeric(n_int)
    oe_ratio_meth   <- numeric(n_int)
    oe_ratio_unmeth <- numeric(n_int)

    for (i in 1:n_int) {
        t_lo <- breaks[i]
        t_hi <- breaks[i + 1]
        oe_time[i] <- (t_lo + t_hi) / 2

        sub <- surv_data[surv_data$time >= t_lo,]
        sub$time2  <- pmin(sub$time, t_hi)
        sub$event2 <- sub$event * (sub$time <= t_hi)

        if (sum(sub$event2) < 2) next

        sd_int <- tryCatch(
            survdiff(Surv(time2, event2) ~ mgmt, data = sub),
            error = function(e) NULL
        )
        if (!is.null(sd_int) && all(sd_int$exp > 0)) {
            oe_ratio_meth[i]   <- sd_int$obs[2] / sd_int$exp[2]
            oe_ratio_unmeth[i] <- sd_int$obs[1] / sd_int$exp[1]
        } else {
            oe_ratio_meth[i]   <- NA
            oe_ratio_unmeth[i] <- NA
        }
    }
    """)

    oe_time         = np.array(ro.r("oe_time"))
    oe_ratio_meth   = np.array(ro.r("oe_ratio_meth"), dtype=float)
    oe_ratio_unmeth = np.array(ro.r("oe_ratio_unmeth"), dtype=float)

    # Replace 0 with NaN (unfilled intervals)
    oe_ratio_meth[oe_ratio_meth == 0]     = np.nan
    oe_ratio_unmeth[oe_ratio_unmeth == 0] = np.nan

    # =============================================================================
    # FIGURE 1 — PH DIAGNOSTICS
    # =============================================================================

    fig1 = plt.figure(figsize=(14, 11))
    gs   = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.35)
    ax1  = fig1.add_subplot(gs[0, 0])
    ax2  = fig1.add_subplot(gs[0, 1])
    ax3  = fig1.add_subplot(gs[1, 0])
    ax4  = fig1.add_subplot(gs[1, 1])

    # Panel 1: Log-log plot
    mask_m = (s_meth   > 0) & (s_meth   < 1) & (t_meth   > 0)
    mask_u = (s_unmeth > 0) & (s_unmeth < 1) & (t_unmeth > 0)
    ax1.plot(np.log(t_meth[mask_m]),   np.log(-np.log(s_meth[mask_m])),
             color=COL_METH,   linewidth=1.5, label=f"Meth. (n={n_meth})")
    ax1.plot(np.log(t_unmeth[mask_u]), np.log(-np.log(s_unmeth[mask_u])),
             color=COL_UNMETH, linewidth=1.5, label=f"Unmeth. (n={n_unmeth})")
    ax1.set_xlabel("log(time)", fontsize=10)
    ax1.set_ylabel("log(−log(S(t)))", fontsize=10)
    ax1.set_title("Log-log survival plot (MGMT)\nParallel lines → PH satisfied",
                  fontsize=10, fontweight="bold")
    ax1.legend(fontsize=9, frameon=False)
    ax1.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel 2: Scaled Schoenfeld residuals — show both MGMT and Age
    ax2.scatter(scho_t, scho_age, color="gray", alpha=0.35, s=12, zorder=2,
                label="Age residuals")
    ax2.plot(lo_sx_a, lo_sy_a, color="#C0392B", linewidth=2.0, zorder=3,
             label=f"Age smooth (p={age_zph_p:.3f})")
    ax2.scatter(scho_t, scho_mgmt, color="lightblue", alpha=0.35, s=12, zorder=2)
    ax2.plot(lo_sx_m, lo_sy_m, color=COL_UNMETH, linewidth=2.0, zorder=3,
             linestyle="--", label=f"MGMT smooth (p={mgmt_zph_p:.3f})")
    ax2.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax2.set_xlabel("Time (years)", fontsize=10)
    ax2.set_ylabel("Scaled Schoenfeld residual", fontsize=10)
    ax2.set_title("Scaled Schoenfeld residuals\nTime trend → PH violated",
                  fontsize=10, fontweight="bold")
    ax2.legend(fontsize=8, frameon=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Panel 3: Time-varying HR (MGMT)
    ax3.scatter(scho_t, np.exp(beta_mgmt), color="gray", alpha=0.4, s=15, zorder=2)
    ax3.plot(lo_tx, np.exp(lo_ty), color="#C0392B", linewidth=2.0, zorder=3,
             label="Time-varying HR (lowess)")
    ax3.axhline(ref_hr_mgmt, color="#2C3E50", linewidth=1.2, linestyle="--",
                label=f"Cox HR = {ref_hr_mgmt:.2f} (constant)")
    ax3.axhline(1.0, color="black", linewidth=0.6, linestyle=":", alpha=0.5,
                label="HR = 1")
    ax3.set_xlabel("Time (years)", fontsize=10)
    ax3.set_ylabel("Hazard ratio (meth. vs unmeth.)", fontsize=10)
    ax3.set_title("Time-varying hazard ratio (MGMT)\nFlat line → PH satisfied",
                  fontsize=10, fontweight="bold")
    ax3.legend(fontsize=9, frameon=False)
    ax3.set_ylim(bottom=0)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    # Panel 4: Time-varying O/E ratio
    valid_m = ~np.isnan(oe_ratio_meth)
    valid_u = ~np.isnan(oe_ratio_unmeth)
    ax4.plot(oe_time[valid_m], oe_ratio_meth[valid_m],
             color=COL_METH,   linewidth=1.5, marker="o", markersize=4,
             label=f"Meth. (n={n_meth})")
    ax4.plot(oe_time[valid_u], oe_ratio_unmeth[valid_u],
             color=COL_UNMETH, linewidth=1.5, marker="o", markersize=4,
             label=f"Unmeth. (n={n_unmeth})")
    ax4.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.7,
                label="O/E = 1 (PH satisfied)")
    ax4.set_xlabel("Time (years)", fontsize=10)
    ax4.set_ylabel("O/E ratio (per interval)", fontsize=10)
    ax4.set_title("Observed / Expected events over time\nConstant ratio → PH satisfied",
                  fontsize=10, fontweight="bold")
    ax4.legend(fontsize=9, frameon=False)
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)

    fig1.suptitle(
        f"Proportional hazards diagnostics — MOSAIC GBM cohort (N={len(df)})\n"
        f"Schoenfeld p: MGMT={mgmt_zph_p:.3f}  Age={age_zph_p:.3f}  "
        f"Sex={sex_zph_p:.3f}  Global={global_zph_p:.3f}",
        fontsize=10, fontweight="bold", y=1.01)

    output.save(SAVE_DIR, f"ph_diagnostics_cohort{suffix}", fig1, panel="Supp Fig 1a" if panels else None)
    plt.show()
    plt.close(fig1)
    print(f"Saved ph_diagnostics_cohort{suffix}.png/.pdf")

    # =============================================================================
    # FIGURE 2 — TAU SELECTION
    # =============================================================================

    h_smooth = uniform_filter1d(h_all, size=max(1, len(t_all) // 15))

    fig2, (ax_km, ax_hz) = plt.subplots(1, 2, figsize=(12, 5))

    for t_arr, s_arr, color, label in [
        (t_meth,   s_meth,   COL_METH,   f"MGMT meth. (n={n_meth})"),
        (t_unmeth, s_unmeth, COL_UNMETH, f"MGMT unmeth. (n={n_unmeth})"),
    ]:
        t_plot = np.concatenate([[0], t_arr])
        s_plot = np.concatenate([[1], s_arr])
        ax_km.step(t_plot, s_plot, where="post",
                   color=color, linewidth=1.5, label=label)

    ax_km.axvline(1.5, color="gray",  linestyle="--", linewidth=1.0, alpha=0.8)
    ax_km.axvline(3.0, color="black", linestyle="--", linewidth=1.0, alpha=0.8)
    ax_km.text(1.55, 0.97, "τ=1.5y", fontsize=8, color="gray",  va="top")
    ax_km.text(3.05, 0.97, "τ=3y",   fontsize=8, color="black", va="top")
    ax_km.set_xlabel("Time (years)", fontsize=10)
    ax_km.set_ylabel("Survival probability", fontsize=10)
    ax_km.set_title("Kaplan-Meier by MGMT status\nwith tau markers",
                    fontsize=10, fontweight="bold")
    ax_km.set_ylim(0, 1.05)
    ax_km.set_xlim(0)
    ax_km.legend(fontsize=9, frameon=False)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)

    ax_hz.plot(t_all, h_smooth, color="#2C3E50", linewidth=2.0,
               label="Smoothed hazard (cohort)")
    ax_hz.axvline(1.5, color="gray",  linestyle="--", linewidth=1.0, alpha=0.8,
                  label=f"τ=1.5y (≈ median OS {median_os:.1f}y)")
    ax_hz.axvline(3.0, color="black", linestyle="--", linewidth=1.0, alpha=0.8,
                  label="τ=3y (inflection)")
    ax_hz.set_xlabel("Time (years)", fontsize=10)
    ax_hz.set_ylabel("Smoothed hazard rate", fontsize=10)
    ax_hz.set_title("Population hazard rate over time\n"
                    "Inflection and plateau justify tau selection",
                    fontsize=10, fontweight="bold")
    ax_hz.legend(fontsize=9, frameon=False)
    ax_hz.spines["top"].set_visible(False)
    ax_hz.spines["right"].set_visible(False)

    fig2.suptitle(
        f"Tau selection justification — MOSAIC GBM cohort (N={len(df)})\n"
        f"Median OS (events) = {median_os:.2f} years",
        fontsize=11, fontweight="bold")
    fig2.tight_layout()
    output.save(SAVE_DIR, f"tau_selection_justification{suffix}", fig2, panel="Supp Fig 1b" if panels else None)
    plt.show()
    plt.close(fig2)
    print(f"Saved tau_selection_justification{suffix}.png/.pdf")

    # =============================================================================
    # SUMMARY
    # =============================================================================

    print("\n" + "="*60)
    print("PROPORTIONAL HAZARDS ASSESSMENT SUMMARY")
    print("="*60)
    print(f"\nN={len(df)}, events={int(df['event'].sum())}, "
          f"median OS={median_os:.2f}y")
    print(f"\nSchoenfeld residuals test:")
    for var, p in [("MGMT",   mgmt_zph_p),
                   ("Age",    age_zph_p),
                   ("Sex",    sex_zph_p),
                   ("Global", global_zph_p)]:
        flag = "*** PH VIOLATED" if p < 0.05 else "(not rejected)"
        print(f"  {var:<8}: p = {p:.4f}  {flag}")
    print(f"\nConclusion:")
    if age_zph_p < 0.05 or global_zph_p < 0.05:
        print(f"  PH assumption violated for at least one covariate.")
        print(f"  RMST is the more appropriate survival estimand for this cohort.")
    else:
        print(f"  Formal Schoenfeld test does not reject PH, but visual")
        print(f"  diagnostics should be inspected for time-varying patterns.")
    print(f"\nTau selection:")
    print(f"  tau=3y  : inflection in population hazard; plateau in KM curves")
    print(f"  tau=1.5y: approximate median OS ({median_os:.2f}y), "
          f"period of maximal event accrual")

    return dict(n=len(df), events=int(df["event"].sum()), mgmt=mgmt_zph_p, age=age_zph_p,
                sex=sex_zph_p, global_=global_zph_p)
