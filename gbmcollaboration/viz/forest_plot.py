import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .figure_helpers import savefig

def forest_plot(
    df: pd.DataFrame,
    
    # --- Column Mappings (FLEXIBLE INPUTS) ---
    col_label: str = "cluster",
    col_n_comp: str = "n_cluster",
    col_n_ref: str = "n_rest",
    col_diff: str = "diff_months_adj",
    col_lo: str = "ci_months_low_adj",
    col_hi: str = "ci_months_high_adj",
    col_p: str = "p_adj",
    
    # --- Scaling (Useful if diff is in months but CIs are in years) ---
    diff_multiplier: float = 1.0,
    ci_multiplier: float = 1.0,
    
    # --- Output & Titling ---
    save_path: str = None,
    title: str = "Adjusted RMST difference by cluster",
    xlabel: str = "RMST difference (months)  |  cluster vs rest",
    tau: int = None,
    
    # --- Figure Dimensions & Resolution ---
    figsize: tuple = (15, 5),
    dpi: int = 300,
    subplots_adj_left: float = 0.38,
    subplots_adj_right: float = 0.72,
    
    # --- Fonts ---
    fs_hdr: float = 8.5,
    fs_row: float = 8.5,
    fs_axis: float = 8.5,
    fs_title: float = 9.5,
    
    # --- Marker & Line Styles ---
    marker: str = "D",
    marker_sz: float = 5,
    cap_sz: float = 3,
    lw: float = 1.3,
    elw: float = 1.1,
    
    # --- Colors ---
    col_sig: str = "#C0392B",
    col_ns: str = "#2C3E50",
    col_n_text: str = "gray",
    col_vline: str = "#7F8C8D",
    grid_alpha: float = 0.4,
    
    # --- Column Positioning (transAxes: 0 is left axis edge, 1 is right axis edge) ---
    x_label: float = -0.25,  # covariate label (right-aligned)
    x_n: float     = -0.02,  # n column        (right-aligned)
    x_est: float   = 1.02,   # estimate + CI   (left-aligned)
    x_p: float     = 1.4    # p-value         (left-aligned)
):
    # Helper to format p-values
    def fmt_p(p):
        if pd.isna(p): return "—"
        if p < 0.0001: return f"{p:.2e}"
        return f"{p:.4f}"

    n_rows = len(df)
    
    # Append tau to titles if provided
    if tau is not None:
        title += f" (τ={int(tau)} years)"
        xlabel = xlabel.replace("months)", f"months, τ={int(tau)} years)")

    # 1. Setup Plot
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(n_rows - 1, -1, -1, dtype=float)

    # 2. Alternating row shading
    for i in range(n_rows):
        if i % 2 == 0:
            ax.axhspan(y_pos[i] - 0.45, y_pos[i] + 0.45, color="whitesmoke", zorder=0)

    # Track min/max for x-axis scaling
    all_lo, all_hi = [], []

    # 3. Draw Error bars
    for i, (_, row) in enumerate(df.iterrows()):
        y = y_pos[i]
        
        # Extract and scale values
        p_val = row[col_p]
        is_sig = False if pd.isna(p_val) else p_val < 0.05
        c_color = col_sig if is_sig else col_ns
        
        val_diff = row[col_diff] * diff_multiplier
        val_lo   = row[col_lo] * ci_multiplier
        val_hi   = row[col_hi] * ci_multiplier
        
        if pd.notna(val_diff) and pd.notna(val_lo) and pd.notna(val_hi):
            ax.errorbar(
                x=val_diff, y=y,
                xerr=[[val_diff - val_lo], [val_hi - val_diff]],
                fmt=marker, color=c_color, markersize=marker_sz,
                capsize=cap_sz, linewidth=lw, elinewidth=elw, zorder=3,
            )
            all_lo.append(val_lo)
            all_hi.append(val_hi)

    # 4. Reference line & Axes styling
    ax.axvline(x=0, color=col_vline, linewidth=0.9, linestyle="--", zorder=2)

    # Dynamic x-axis limits based on SCALED values
    x_min_data = min(all_lo) if all_lo else -5
    x_max_data = max(all_hi) if all_hi else 5
    x_pad = max((x_max_data - x_min_data) * 0.08, 1.0)
    
    ax.set_xlim(x_min_data - x_pad, x_max_data + x_pad)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.set_yticks([]) # Hide y-ticks
    
    ax.set_xlabel(xlabel, fontsize=fs_axis)
    ax.set_title(title + "\n", fontsize=fs_title, pad=8)
    
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", labelsize=fs_axis, width=0.6, length=3)
    ax.grid(axis="x", linewidth=0.4, alpha=grid_alpha, zorder=0)

    # 5. Text Annotations
    trans = ax.transAxes
    y_lo_ax, y_hi_ax = ax.get_ylim()
    y_to_frac = lambda y_data: (y_data - y_lo_ax) / (y_hi_ax - y_lo_ax)

    # Headers
    header_frac = y_to_frac(n_rows - 0.05)
    ax.text(x_label, header_frac, "Cluster", transform=trans, fontsize=fs_hdr, fontweight="bold", ha="right", va="bottom")
    ax.text(x_n,     header_frac, "n (clust/rest)", transform=trans, fontsize=fs_hdr, fontweight="bold", ha="right", va="bottom")
    ax.text(x_est,   header_frac, "Estimate [95% CI] (months)", transform=trans, fontsize=fs_hdr, fontweight="bold", ha="left", va="bottom")
    ax.text(x_p,     header_frac, "adj p", transform=trans, fontsize=fs_hdr, fontweight="bold", ha="left", va="bottom")

    # Rows
    for i, (_, row) in enumerate(df.iterrows()):
        yf = y_to_frac(y_pos[i])
        
        p_val = row[col_p]
        is_sig = False if pd.isna(p_val) else p_val < 0.05
        c_color = col_sig if is_sig else col_ns

        # Covariate label (cast to string in case it's an integer like 0, 1, 2)
        lbl_str = row[col_label]
        ax.text(x_label, yf, lbl_str, transform=trans, fontsize=fs_row, ha="right", va="center", color=col_ns)

        # N
        n_str = f"{int(row[col_n_comp])}/{int(row[col_n_ref])}" if pd.notna(row[col_n_comp]) else "—/—"
        ax.text(x_n, yf, n_str, transform=trans, fontsize=fs_row, ha="right", va="center", color=col_n_text)

        # Estimate and P-value
        val_diff = row[col_diff] * diff_multiplier
        val_lo   = row[col_lo] * ci_multiplier
        val_hi   = row[col_hi] * ci_multiplier

        if pd.notna(val_diff):
            est_str = f"{val_diff:+.1f} [{val_lo:+.1f}, {val_hi:+.1f}]"
            p_str = fmt_p(p_val) + (" *" if is_sig else "")
        else:
            est_str, p_str = "—", "—"

        ax.text(x_est, yf, est_str, transform=trans, fontsize=fs_row, ha="left", va="center", color=c_color)
        ax.text(x_p, yf, p_str, transform=trans, fontsize=fs_row, ha="left", va="center", color=c_color)

    # 6. Final Adjustments & Save
    fig.subplots_adjust(left=subplots_adj_left, right=subplots_adj_right)

    if save_path:
        savefig(df=df, out_file=save_path)
    
    plt.show()