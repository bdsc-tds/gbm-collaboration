"""
Uses continuous relative Hazard Ratios and RMST to report effects at a given biomarker values.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
from lifelines import CoxPHFitter
from scipy.integrate import trapezoid
from scipy.stats import norm, percentileofscore

# =============================================================================
# 1. MODEL FITTING
# =============================================================================


def fit_spline_cox_model(df, time_col, event_col, cell_col, df_spline=3, cov_names=[], strata=[], treatment_col=None):
    required_cols = [time_col, event_col, cell_col] + list(cov_names) + list(strata)
    if treatment_col:
        required_cols.append(treatment_col)

    df_clean = df[required_cols].dropna().copy()

    # Build Patsy Formula
    if treatment_col:
        spline_formula = f"Q('{treatment_col}') * cr(Q('{cell_col}'), df={df_spline}, constraints='center')"
    else:
        spline_formula = f"cr(Q('{cell_col}'), df={df_spline}, constraints='center')"

    if cov_names:
        cov_str = " + ".join([f"Q('{c}')" for c in cov_names])
        formula = f"{spline_formula} + {cov_str}"
    else:
        formula = spline_formula

    X = patsy.dmatrix(formula + " - 1", data=df_clean, return_type="dataframe")
    df_cox = pd.concat([df_clean[[time_col, event_col] + list(strata)], X], axis=1)

    cph = CoxPHFitter()
    cph.fit(df_cox, duration_col=time_col, event_col=event_col, strata=strata if strata else None)

    # Extract overall, non-linearity, and interaction p-values
    cell_cols = [col for col in X.columns if "cr(" in col and ":" not in col]
    interaction_cols = [col for col in X.columns if "cr(" in col and ":" in col]

    p_overall = cph.summary.loc[cell_cols, "p"].min() if cell_cols else np.nan
    p_non_linear = cph.summary.loc[cell_cols[1:], "p"].min() if len(cell_cols) > 1 else np.nan
    p_interaction = cph.summary.loc[interaction_cols, "p"].min() if interaction_cols else np.nan

    return cph, X.design_info, df_clean, p_overall, p_non_linear, p_interaction


# =============================================================================
# 2. CURVE PREDICTIONS (GLOBAL REFERENCE SUPPORTED)
# =============================================================================


def predict_continuous_curves(
    cph,
    design_info,
    df,
    cell_col,
    tau,
    n_points=200,
    cov_names=[],
    strata=[],
    treatment_col=None,
    bootstrap_CIs=True,
    B=100,
    time_col=None,
    event_col=None,
    global_reference=True,
):
    """
    Predicts continuous RMST and Hazard Ratios.
    If global_reference=True and treatment_col is active, both curves are compared
    to a single baseline (Control arm at the minimum biomarker value).
    """
    x_range = np.linspace(df[cell_col].quantile(0.01), df[cell_col].quantile(0.99), n_points)
    time_grid = np.linspace(0, tau, 200)
    coef_cols = cph.params_.index
    cov_mat = cph.variance_matrix_.values

    results = {"x_range": x_range}

    # Pre-build prediction design matrices
    X_preds = {}

    if treatment_col:
        trt_levels = sorted(list(df[treatment_col].unique()))
        trt_values = [trt_levels[-1], trt_levels[0]]  # [Treated/Arm1, Control/Arm0]
    else:
        trt_values = [None]

    # --- A. Pre-build SINGLE Global Reference Profile if requested ---
    if global_reference:
        df_ref_global = pd.DataFrame({cell_col: [x_range[0]]})
        if treatment_col:
            df_ref_global[treatment_col] = trt_levels[0]  # Arm 0 (Control)
            if isinstance(df[treatment_col].dtype, pd.CategoricalDtype):
                df_ref_global[treatment_col] = df_ref_global[treatment_col].astype(df[treatment_col].dtype)

        for cov in cov_names:
            df_ref_global[cov] = df[cov].median() if pd.api.types.is_numeric_dtype(df[cov]) else df[cov].mode()[0]
        for st in strata:
            df_ref_global[st] = df[st].median() if pd.api.types.is_numeric_dtype(df[st]) else df[st].mode()[0]

        X_ref_global = patsy.build_design_matrices([design_info], df_ref_global, return_type="dataframe")[0]
        for st in strata:
            X_ref_global[st] = df_ref_global[st].values
        X_ref_global_row = X_ref_global[coef_cols].iloc[0].values

    # --- B. Calculate Point Estimates ---
    for idx, trt in enumerate(trt_values):
        df_pred = pd.DataFrame({cell_col: x_range})
        if treatment_col:
            df_pred[treatment_col] = trt
            if isinstance(df[treatment_col].dtype, pd.CategoricalDtype):
                df_pred[treatment_col] = df_pred[treatment_col].astype(df[treatment_col].dtype)

        for cov in cov_names:
            df_pred[cov] = df[cov].median() if pd.api.types.is_numeric_dtype(df[cov]) else df[cov].mode()[0]
        for st in strata:
            df_pred[st] = df[st].median() if pd.api.types.is_numeric_dtype(df[st]) else df[st].mode()[0]

        X_p = patsy.build_design_matrices([design_info], df_pred, return_type="dataframe")[0]
        for st in strata:
            X_p[st] = df_pred[st].values

        # Predict RMST
        surv_funcs = cph.predict_survival_function(X_p, times=time_grid)
        rmst_curve = trapezoid(surv_funcs.values, x=surv_funcs.index, axis=0)

        # Predict Hazard Ratio relative to reference
        if global_reference:
            X_diff = X_p[coef_cols].values - X_ref_global_row
        else:
            ref_val = x_range[0]
            df_ref = df_pred.copy()
            df_ref[cell_col] = ref_val
            X_ref = patsy.build_design_matrices([design_info], df_ref, return_type="dataframe")[0]
            for st in strata:
                X_ref[st] = df_ref[st].values
            X_diff = X_p[coef_cols].values - X_ref[coef_cols].values

        log_hr = X_diff.dot(cph.params_.values)
        var_log_hr = (X_diff.dot(cov_mat) * X_diff).sum(axis=1)
        se_log_hr = np.sqrt(np.maximum(0, var_log_hr))

        hr_curve = np.exp(log_hr)
        hr_lower = np.exp(log_hr - 1.96 * se_log_hr)
        hr_upper = np.exp(log_hr + 1.96 * se_log_hr)

        suffix = f"_trt_{1 - idx}" if treatment_col is not None else ""
        results[f"rmst_curve{suffix}"] = rmst_curve
        results[f"hr_curve{suffix}"] = hr_curve
        results[f"hr_lower{suffix}"] = hr_lower
        results[f"hr_upper{suffix}"] = hr_upper

        X_preds[trt] = X_p

    # --- C. Fast & Robust Bootstrap ---
    if bootstrap_CIs and time_col and event_col:
        boot_curves = {trt: [] for trt in trt_values}

        required_cols = [time_col, event_col, cell_col] + list(cov_names) + list(strata)
        if treatment_col:
            required_cols.append(treatment_col)

        df_clean = df[required_cols].dropna().copy()

        groupby_cols = []
        if strata:
            groupby_cols.extend(strata)
        if treatment_col:
            groupby_cols.append(treatment_col)

        for b in range(B):
            if groupby_cols:
                df_boot = df_clean.groupby(groupby_cols, group_keys=False).sample(frac=1.0, replace=True)
            else:
                df_boot = df_clean.sample(frac=1.0, replace=True)

            try:
                X_boot = patsy.build_design_matrices([design_info], df_boot, return_type="dataframe")[0]
                df_cox_boot = pd.concat([df_boot[[time_col, event_col] + list(strata)], X_boot], axis=1)

                cph_b = CoxPHFitter(penalizer=0.01)
                cph_b.fit(df_cox_boot, duration_col=time_col, event_col=event_col, strata=strata if strata else None)

                for trt in trt_values:
                    surv_b = cph_b.predict_survival_function(X_preds[trt], times=time_grid)
                    rmst_b = trapezoid(surv_b.values, x=surv_b.index, axis=0)
                    boot_curves[trt].append(rmst_b)
            except Exception:
                continue

        for idx, trt in enumerate(trt_values):
            curves_list = boot_curves[trt]
            if len(curves_list) > 5:
                curves_arr = np.array(curves_list)
                lower_ci = np.percentile(curves_arr, 2.5, axis=0)
                upper_ci = np.percentile(curves_arr, 97.5, axis=0)

                suffix = f"_trt_{1 - idx}" if treatment_col is not None else ""
                results[f"rmst_lower{suffix}"] = lower_ci
                results[f"rmst_upper{suffix}"] = upper_ci

    return results


# =============================================================================
# 3. CALCULUS / INFLECTION
# =============================================================================


def detect_inflection_point(x_range, y_curve, raw_data_series, edge_trim=0.05):
    dx = x_range[1] - x_range[0]
    n = len(x_range)
    trim_idx = int(n * edge_trim)
    valid_slice = slice(trim_idx, n - trim_idx)

    x_valid = x_range[valid_slice]
    y_valid = y_curve[valid_slice]

    d_y = np.gradient(y_curve, dx)
    d2_y = np.gradient(d_y, dx)

    idx_peak_valid = np.argmax(y_valid)
    x_peak = x_valid[idx_peak_valid]
    pct_peak = percentileofscore(raw_data_series.dropna(), x_peak)

    idx_inflect_valid = np.argmax(np.abs(d2_y[valid_slice]))
    x_inflection = x_valid[idx_inflect_valid]
    pct_inflection = percentileofscore(raw_data_series.dropna(), x_inflection)

    return {"x_peak": x_peak, "pct_peak": pct_peak, "x_inflection": x_inflection, "pct_inflection": pct_inflection}


# =============================================================================
# 4. CONTRASTS (75th vs 25th)
# =============================================================================


def compute_rmst_contrasts(
    cph, design_info, df, cell_col, tau, cov_names=[], strata=[], treatment_col=None, treatment_val=None
):
    val_q75 = df[cell_col].quantile(0.75)
    val_q25 = df[cell_col].quantile(0.25)

    df_contrast = pd.DataFrame({cell_col: [val_q75, val_q25]})
    if treatment_col is not None and treatment_val is not None:
        df_contrast[treatment_col] = treatment_val
        if isinstance(df[treatment_col].dtype, pd.CategoricalDtype):
            df_contrast[treatment_col] = df_contrast[treatment_col].astype(df[treatment_col].dtype)

    for cov in cov_names:
        df_contrast[cov] = df[cov].median() if pd.api.types.is_numeric_dtype(df[cov]) else df[cov].mode()[0]
    for st in strata:
        df_contrast[st] = df[st].median() if pd.api.types.is_numeric_dtype(df[st]) else df[st].mode()[0]

    X_contrast = patsy.build_design_matrices([design_info], df_contrast, return_type="dataframe")[0]
    for st in strata:
        X_contrast[st] = df_contrast[st].values

    time_grid = np.linspace(0, tau, 200)
    surv_contrast = cph.predict_survival_function(X_contrast, times=time_grid)
    rmst_q75 = trapezoid(surv_contrast.iloc[:, 0].values, x=surv_contrast.index)
    rmst_q25 = trapezoid(surv_contrast.iloc[:, 1].values, x=surv_contrast.index)

    rmst_diff = rmst_q75 - rmst_q25

    coef_cols = cph.params_.index
    X_diff = X_contrast[coef_cols].iloc[0].values - X_contrast[coef_cols].iloc[1].values
    var_log_hr = X_diff.dot(cph.variance_matrix_.values).dot(X_diff)
    se_diff = np.abs(rmst_diff) * np.sqrt(np.maximum(0, var_log_hr))

    p_val_diff = 2 * (1 - norm.cdf(np.abs(rmst_diff / se_diff))) if se_diff > 0 else np.nan
    ci_lower = rmst_diff - 1.96 * se_diff
    ci_upper = rmst_diff + 1.96 * se_diff

    return {
        "val_q75": val_q75,
        "val_q25": val_q25,
        "rmst_q75": rmst_q75,
        "rmst_q25": rmst_q25,
        "rmst_diff": rmst_diff,
        "se_diff": se_diff,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "p_val_diff": p_val_diff,
    }


# =============================================================================
# 5. PLOTTING (GLOBAL REFERENCE LABEL SUPPORTED)
# =============================================================================


def plot_continuous_results(
    x_range, curves, inflection_dict, df, cell_col, tau, plot_type="rmst", treatment_col=None, global_reference=True
):
    val_q75 = df[cell_col].quantile(0.75)
    has_trt = treatment_col is not None

    if plot_type == "both":
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 6))

    # --- PLOT RMST ---
    if plot_type in ["rmst", "both"]:
        ax_rmst = ax1
        if has_trt:
            ax_rmst.plot(
                x_range, curves["rmst_curve_trt_1"], color="darkred", lw=2.5, label=f"Treated Arm 1 (tau={tau}$)"
            )
            if "rmst_lower_trt_1" in curves:
                ax_rmst.fill_between(
                    x_range, curves["rmst_lower_trt_1"], curves["rmst_upper_trt_1"], color="red", alpha=0.10
                )

            ax_rmst.plot(
                x_range, curves["rmst_curve_trt_0"], color="darkblue", lw=2.5, label=f"Control Arm 0 (tau={tau}$)"
            )
            if "rmst_lower_trt_0" in curves:
                ax_rmst.fill_between(
                    x_range, curves["rmst_lower_trt_0"], curves["rmst_upper_trt_0"], color="blue", alpha=0.10
                )

            ax_rmst.axvline(
                x=inflection_dict["trt_1"]["x_peak"],
                color="red",
                linestyle=":",
                alpha=0.8,
                label=f"Peak Trt1 ({inflection_dict['trt_1']['x_peak']:.3f})",
            )
            ax_rmst.axvline(
                x=inflection_dict["trt_0"]["x_peak"],
                color="blue",
                linestyle=":",
                alpha=0.8,
                label=f"Peak Trt0 ({inflection_dict['trt_0']['x_peak']:.3f})",
            )
        else:
            ax_rmst.plot(x_range, curves["rmst_curve"], color="darkred", lw=2.5, label=f"Model RMST (tau={tau}$)")
            if "rmst_lower" in curves:
                ax_rmst.fill_between(
                    x_range, curves["rmst_lower"], curves["rmst_upper"], color="red", alpha=0.12, label="95% CI"
                )

            ax_rmst.axvline(
                x=inflection_dict["x_peak"],
                color="darkgreen",
                linestyle="-",
                lw=2,
                label=f"Peak ({inflection_dict['x_peak']:.3f} | {inflection_dict['pct_peak']:.1f}th %)",
            )
            ax_rmst.axvline(
                x=inflection_dict["x_inflection"],
                color="purple",
                linestyle="-.",
                lw=1.8,
                label=f"Inflection ({inflection_dict['x_inflection']:.3f})",
            )

        ax_rmst.axvline(x=val_q75, color="black", linestyle="--", alpha=0.6, label=f"Q75 ({val_q75:.3f})")
        ax_rmst.set_ylabel(f"RMST (Mean Survival up to {tau})", fontsize=11)
        ax_rmst.set_title(
            f"Continuous RMST Interaction Analysis ({cell_col})"
            if has_trt
            else f"Continuous RMST Analysis ({cell_col})",
            fontsize=12,
            fontweight="bold",
        )
        ax_rmst.legend(loc="best", frameon=True)
        ax_rmst.grid(True, alpha=0.3)

    # --- PLOT HAZARD RATIO ---
    if plot_type in ["hazard", "both"]:
        ax_hr = ax2 if plot_type == "both" else ax1

        # --- FIXED: Label changes dynamically based on global_reference ---
        if global_reference and has_trt:
            ref_label = f"rel to Control at {x_range[0]:.3f}"
        else:
            ref_label = f"rel to {x_range[0]:.3f}"

        if has_trt:
            ax_hr.plot(x_range, curves["hr_curve_trt_1"], color="darkred", lw=2, label=f"HR Treated ({ref_label})")
            ax_hr.fill_between(x_range, curves["hr_lower_trt_1"], curves["hr_upper_trt_1"], color="red", alpha=0.10)

            ax_hr.plot(x_range, curves["hr_curve_trt_0"], color="darkblue", lw=2, label=f"HR Control ({ref_label})")
            ax_hr.fill_between(x_range, curves["hr_lower_trt_0"], curves["hr_upper_trt_0"], color="blue", alpha=0.10)
        else:
            ax_hr.plot(x_range, curves["hr_curve"], color="darkblue", lw=2.5, label=f"Hazard Ratio ({ref_label})")
            ax_hr.fill_between(
                x_range, curves["hr_lower"], curves["hr_upper"], color="blue", alpha=0.15, label="95% CI"
            )

        ax_hr.axhline(y=1.0, color="black", linestyle="-", alpha=0.5)
        ax_hr.set_yscale("log")
        ax_hr.set_ylabel("Hazard Ratio (Log Scale)", fontsize=11)
        ax_hr.set_xlabel(f"{cell_col} (Continuous Fraction)", fontsize=11)
        ax_hr.legend(loc="best", frameon=True)
        ax_hr.grid(True, alpha=0.3)

    if plot_type == "rmst":
        ax1.set_xlabel(f"{cell_col} (Continuous Fraction)", fontsize=11)

    plt.tight_layout()
    plt.show()


# =============================================================================
# 6. MASTER ORCHESTRATOR (GLOBAL REFERENCE PASSED CLEANLY)
# =============================================================================


def master_threshold_analysis(
    df,
    time_col,
    event_col,
    cell_col,
    tau,
    cov_names=[],
    strata=[],
    treatment_col=None,
    plot_type="rmst",
    global_reference=True,
):
    """
    Generalized Master Pipeline.
    Supports continuous spline-only models OR Treatment * Spline interaction models.
    Supports global_reference=True to benchmark both treatment hazard curves to Control at x_min.
    """
    # 1. Fit Model
    cph, design_info, df_clean, p_overall, p_non_linear, p_interaction = fit_spline_cox_model(
        df, time_col, event_col, cell_col, cov_names=cov_names, strata=strata, treatment_col=treatment_col
    )

    # 2. Predict Continuous Curves (With Bootstrap & Global Reference options)
    curves = predict_continuous_curves(
        cph,
        design_info,
        df_clean,
        cell_col,
        tau,
        cov_names=cov_names,
        strata=strata,
        treatment_col=treatment_col,
        bootstrap_CIs=True,
        B=100,
        time_col=time_col,
        event_col=event_col,
        global_reference=global_reference,
    )

    # 3. Detect Inflection/Peak Points & Compute Contrasts
    if treatment_col:
        inflection_dict = {
            "trt_1": detect_inflection_point(curves["x_range"], curves["rmst_curve_trt_1"], df_clean[cell_col]),
            "trt_0": detect_inflection_point(curves["x_range"], curves["rmst_curve_trt_0"], df_clean[cell_col]),
        }
        trt_levels = sorted(list(df_clean[treatment_col].unique()))
        trt_1_val = trt_levels[-1]
        trt_0_val = trt_levels[0]

        contrasts_1 = compute_rmst_contrasts(
            cph,
            design_info,
            df_clean,
            cell_col,
            tau,
            cov_names=cov_names,
            strata=strata,
            treatment_col=treatment_col,
            treatment_val=trt_1_val,
        )
        contrasts_0 = compute_rmst_contrasts(
            cph,
            design_info,
            df_clean,
            cell_col,
            tau,
            cov_names=cov_names,
            strata=strata,
            treatment_col=treatment_col,
            treatment_val=trt_0_val,
        )

        # Build Advanced INTERACTION Summary Table
        summary_table = pd.DataFrame(
            {
                "Metric": [
                    f"Treated ({trt_1_val}) Peak Survival (Max RMST)",
                    f"Treated ({trt_1_val}) Peak Survival Percentile",
                    f"Control ({trt_0_val}) Peak Survival (Max RMST)",
                    f"Control ({trt_0_val}) Peak Survival Percentile",
                    f"Treated ({trt_1_val}) RMST Diff (75th vs 25th) @ {contrasts_1['val_q75']:.3f} vs {contrasts_1['val_q25']:.3f}",
                    f"Control ({trt_0_val}) RMST Diff (75th vs 25th) @ {contrasts_0['val_q75']:.3f} vs {contrasts_0['val_q25']:.3f}",
                    "Overall Continuous Spline p-value",
                    "Spline Non-linearity p-value (Threshold Proof)",
                    "Treatment-Spline Interaction p-value (Interaction Proof)",
                ],
                "Estimate": [
                    inflection_dict["trt_1"]["x_peak"],
                    inflection_dict["trt_1"]["pct_peak"],
                    inflection_dict["trt_0"]["x_peak"],
                    inflection_dict["trt_0"]["pct_peak"],
                    contrasts_1["rmst_diff"],
                    contrasts_0["rmst_diff"],
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "SE": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    contrasts_1["se_diff"],
                    contrasts_0["se_diff"],
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "Lower .95": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    contrasts_1["ci_lower"],
                    contrasts_0["ci_lower"],
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "Upper .95": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    contrasts_1["ci_upper"],
                    contrasts_0["ci_upper"],
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "p-value": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    contrasts_1["p_val_diff"],
                    contrasts_0["p_val_diff"],
                    p_overall,
                    p_non_linear,
                    p_interaction,
                ],
            }
        )
    else:
        inflection_dict = detect_inflection_point(curves["x_range"], curves["rmst_curve"], df_clean[cell_col])
        contrasts = compute_rmst_contrasts(
            cph, design_info, df_clean, cell_col, tau, cov_names=cov_names, strata=strata
        )

        summary_table = pd.DataFrame(
            {
                "Metric": [
                    "Peak Survival Value (Max RMST)",
                    "Peak Survival Empirical Percentile",
                    "Interior Inflection Value (Max Bending)",
                    "Interior Inflection Empirical Percentile",
                    f"RMST at 75th Percentile ({contrasts['val_q75']:.3f})",
                    f"RMST at 25th Percentile ({contrasts['val_q25']:.3f})",
                    "RMST Difference (75th vs 25th)",
                    "Overall Continuous Feature p-value",
                    "Non-linearity p-value (Threshold Proof)",
                ],
                "Estimate": [
                    inflection_dict["x_peak"],
                    inflection_dict["pct_peak"],
                    inflection_dict["x_inflection"],
                    inflection_dict["pct_inflection"],
                    contrasts["rmst_q75"],
                    contrasts["rmst_q25"],
                    contrasts["rmst_diff"],
                    np.nan,
                    np.nan,
                ],
                "SE": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, contrasts["se_diff"], np.nan, np.nan],
                "Lower .95": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, contrasts["ci_lower"], np.nan, np.nan],
                "Upper .95": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, contrasts["ci_upper"], np.nan, np.nan],
                "p-value": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    contrasts["p_val_diff"],
                    p_overall,
                    p_non_linear,
                ],
            }
        )

    # 4. Plot (With global_reference passed cleanly)
    plot_continuous_results(
        curves["x_range"],
        curves,
        inflection_dict,
        df_clean,
        cell_col,
        tau,
        plot_type=plot_type,
        treatment_col=treatment_col,
        global_reference=global_reference,
    )

    return summary_table, cph
