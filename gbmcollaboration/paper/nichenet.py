"""NicheNet circos plots of the selected communication axes (Fig 6a, Supp Fig 6a-b).

The R code is the circos cell of analysis/pipelines/NicheNet_v3.ipynb (cells 80, 116 and
152, identical apart from the zone), kept as written and run through rpy2. Only the head of
the cell changed: the input CSV, the zone label of the title and the output PNG are
arguments of the R function `circos_zone()`, and the function returns the plotted axes.
Needs R circlize, dplyr and readr in the active environment.
"""

from pathlib import Path

import pandas as pd

from . import config  # noqa: F401  (sets R_HOME for rpy2)

CIRCOS_R = r'''
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(circlize)
})

options(bitmapType = "cairo")

circos_zone <- function(AXES_SELECTED_CSV, zone_dir, OUT_PNG) {
stopifnot(file.exists(AXES_SELECTED_CSV))

MAX_EDGES <- Inf

LINK_ALPHA <- 0.5

BORDER_COL <- "#777777"
BORDER_LWD <- 0.25

TITLE_CEX <- 1.00

START_DEGREE <- -8
SMALL_GAP <- 1
BIG_GAP <- 5

WIDTH_MIN <- 0.10
WIDTH_MAX <- 2.50

COLOR_CLIP_LO <- 0.02
COLOR_CLIP_HI <- 0.98

CANVAS_PAD <- 0.40

LABEL_CEX <- 0.46
LABEL_RADIUS <- 1.0
LABEL_NICEFACING <- TRUE
LABEL_MIN_DEG <- 0

PNG_SIZE <- 2600
PNG_RES <- 240

axes <- readr::read_csv(AXES_SELECTED_CSV, show_col_types = FALSE)

required_cols <- c(
  "axis_id","direction",
  "sender_celltype","receiver_celltype",
  "ligand","target",
  "rank_score","ligand_activity"
)
missing <- setdiff(required_cols, colnames(axes))
if (length(missing) > 0) stop(paste0("axes_selected.csv missing columns: ", paste(missing, collapse = ", ")))

axes <- axes %>%
  transmute(
    axis_id = as.character(axis_id),
    direction = as.character(direction),
    sender = as.character(sender_celltype),
    receiver = as.character(receiver_celltype),
    ligand = as.character(ligand),
    target = as.character(target),
    rank_score = as.numeric(rank_score),
    ligand_activity = as.numeric(ligand_activity)
  ) %>%
  filter(
    direction %in% c("TIL_high","TIL_low"),
    is.finite(rank_score),
    is.finite(ligand_activity)
  )

if (nrow(axes) == 0) stop("No rows left after filtering. Check direction/rank_score/ligand_activity.")

if (is.finite(MAX_EDGES) && nrow(axes) > MAX_EDGES) {
  axes <- axes %>% slice_head(n = MAX_EDGES)
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) == 0) return(rep(0.5, length(x)))
  lo <- min(x[ok]); hi <- max(x[ok])
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(rep(0.5, length(x)))
  y <- (x - lo) / (hi - lo)
  y[!is.finite(y)] <- 0.5
  pmin(pmax(y, 0), 1)
}

clip01_quantile <- function(x, qlo = 0.02, qhi = 0.98) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) == 0) return(rep(0.5, length(x)))
  lo <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = qlo, na.rm = TRUE, names = FALSE)))
  hi <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = qhi, na.rm = TRUE, names = FALSE)))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(rescale01(x))
  y <- (x - lo) / (hi - lo)
  y[!is.finite(y)] <- 0.5
  pmin(pmax(y, 0), 1)
}

axes <- axes %>%
  mutate(
    width_norm = rescale01(ligand_activity),
    w = WIDTH_MIN + width_norm * (WIDTH_MAX - WIDTH_MIN)
  )

rank_norm <- clip01_quantile(axes$rank_score, COLOR_CLIP_LO, COLOR_CLIP_HI)

col_high <- circlize::colorRamp2(c(0, 1), c("#FEE6CE", "#D55E00"))
col_low  <- circlize::colorRamp2(c(0, 1), c("#DEEBF7", "#0072B2"))

axes$link_col <- vapply(seq_len(nrow(axes)), function(i) {
  if (axes$direction[i] == "TIL_high") col_high(rank_norm[i]) else col_low(rank_norm[i])
}, character(1))
axes$link_col <- adjustcolor(axes$link_col, alpha.f = LINK_ALPHA)

S_PREFIX <- "S|"
R_PREFIX <- "R|"
clean_label <- function(x) sub("^[SR]\\|", "", x)

senders <- axes %>%
  group_by(sender) %>%
  summarise(total = sum(w), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(sender)

receivers <- axes %>%
  group_by(receiver) %>%
  summarise(total = sum(w), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(receiver)

sender_sectors <- paste0(S_PREFIX, senders)
receiver_sectors <- paste0(R_PREFIX, receivers)
sectors <- c(sender_sectors, receiver_sectors)

ct_names <- sort(unique(c(senders, receivers)))
ct_cols <- setNames(grDevices::hcl.colors(length(ct_names), palette = "Dark 3"), ct_names)

grid_cols <- c(
  setNames(ct_cols[clean_label(sender_sectors)], sender_sectors),
  setNames(ct_cols[clean_label(receiver_sectors)], receiver_sectors)
)

axes <- axes %>%
  mutate(
    s_sector = paste0(S_PREFIX, sender),
    r_sector = paste0(R_PREFIX, receiver),
    row_id = seq_len(n())
  )

sender_totals <- tapply(axes$w, axes$s_sector, sum)
receiver_totals <- tapply(axes$w, axes$r_sector, sum)

totals <- setNames(rep(1, length(sectors)), sectors)
if (length(sender_totals) > 0) totals[names(sender_totals)] <- sender_totals
if (length(receiver_totals) > 0) totals[names(receiver_totals)] <- receiver_totals
totals[!is.finite(totals) | totals <= 0] <- 1

xlim_mat <- cbind(rep(0, length(sectors)), as.numeric(totals[sectors]))
rownames(xlim_mat) <- sectors

# Collapse ligands per sender:
# Allocate one contiguous sender-side block per (sender, ligand), then subdivide within block for links.
lig_groups <- axes %>%
  group_by(s_sector, ligand) %>%
  summarise(
    group_w = sum(w),
    group_order = min(row_id),
    .groups = "drop"
  ) %>%
  arrange(s_sector, group_order)

sender_cursor <- setNames(rep(0, length(sender_sectors)), sender_sectors)
lig_groups$lig_start <- NA_real_
lig_groups$lig_end <- NA_real_

for (k in seq_len(nrow(lig_groups))) {
  ss <- lig_groups$s_sector[k]
  ww <- lig_groups$group_w[k]
  st <- sender_cursor[ss]
  en <- st + ww
  lig_groups$lig_start[k] <- st
  lig_groups$lig_end[k] <- en
  sender_cursor[ss] <- en
}

axes <- axes %>%
  left_join(lig_groups %>% select(s_sector, ligand, lig_start, lig_end), by = c("s_sector", "ligand")) %>%
  group_by(s_sector, ligand) %>%
  arrange(row_id, .by_group = TRUE) %>%
  mutate(
    s_x1 = lig_start + dplyr::lag(cumsum(w), default = 0),
    s_x2 = lig_start + cumsum(w)
  ) %>%
  ungroup()

# Receiver allocation stays per-edge (unchanged)
receiver_cursor <- setNames(rep(0, length(receiver_sectors)), receiver_sectors)
r_x1 <- numeric(nrow(axes)); r_x2 <- numeric(nrow(axes))

for (i in seq_len(nrow(axes))) {
  rr <- axes$r_sector[i]
  ww <- axes$w[i]
  r_x1[i] <- receiver_cursor[rr]
  r_x2[i] <- receiver_cursor[rr] + ww
  receiver_cursor[rr] <- r_x2[i]
}

axes$r_x1 <- r_x1; axes$r_x2 <- r_x2
axes$s_mid <- (axes$s_x1 + axes$s_x2) / 2
axes$r_mid <- (axes$r_x1 + axes$r_x2) / 2

gap_after <- c(
  rep(SMALL_GAP, max(length(sender_sectors) - 1, 0)), BIG_GAP,
  rep(SMALL_GAP, max(length(receiver_sectors) - 1, 0)), BIG_GAP
)

thin_by_min_deg <- function(theta, min_deg) {
  if (!is.finite(min_deg) || min_deg <= 0) return(rep(TRUE, length(theta)))
  ord <- order(theta)
  keep <- rep(FALSE, length(theta))
  last <- NA_real_
  for (j in ord) {
    if (!is.finite(theta[j])) next
    if (is.na(last) || abs(theta[j] - last) >= min_deg) {
      keep[j] <- TRUE
      last <- theta[j]
    }
  }
  keep
}

plot_one <- function() {
  circos.clear()
  circos.par(
    start.degree = START_DEGREE,
    gap.after = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01),
    canvas.xlim = c(-(1 + CANVAS_PAD), 1 + CANVAS_PAD),
    canvas.ylim = c(-(1 + CANVAS_PAD), 1 + CANVAS_PAD)
  )

  circos.initialize(factors = sectors, xlim = xlim_mat)

  circos.trackPlotRegion(
    factors = sectors,
    ylim = c(0, 1),
    track.height = 0.08,
    bg.border = NA,
    panel.fun = function(x, y) {
      si <- get.cell.meta.data("sector.index")
      xl <- get.cell.meta.data("xlim")
      circos.rect(xl[1], 0, xl[2], 1, col = grid_cols[si], border = NA)
    }
  )

  for (i in seq_len(nrow(axes))) {
    circos.link(
      sector.index1 = axes$s_sector[i],
      point1 = c(axes$s_x1[i], axes$s_x2[i]),
      sector.index2 = axes$r_sector[i],
      point2 = c(axes$r_x1[i], axes$r_x2[i]),
      col = axes$link_col[i],
      border = BORDER_COL,
      lwd = BORDER_LWD,
      directional = 1,
      arr.type = "big.arrow",
      arr.length = 0.04
    )
  }

  # Sender-side ligand labels: one per (sender, ligand), positioned at the center of the ligand block
  ligand_labels <- lig_groups %>%
    mutate(lig_mid = (lig_start + lig_end) / 2)

  theta_lig <- rep(NA_real_, nrow(ligand_labels))
  for (k in seq_len(nrow(ligand_labels))) {
    ss <- ligand_labels$s_sector[k]
    m <- circlize(ligand_labels$lig_mid[k], 0.5, sector.index = ss, track.index = 1)
    theta_lig[k] <- m[1, 1]
  }

  keep_lig <- rep(TRUE, nrow(ligand_labels))
  if (LABEL_MIN_DEG > 0) {
    for (ss in unique(ligand_labels$s_sector)) {
      idx <- which(ligand_labels$s_sector == ss)
      keep_lig[idx] <- thin_by_min_deg(theta_lig[idx], LABEL_MIN_DEG)
    }
  }

  # Receiver-side target labels remain per-edge (unchanged)
  theta_r <- rep(NA_real_, nrow(axes))
  for (rr in unique(axes$r_sector)) {
    idx <- which(is.finite(axes$r_mid) & axes$r_sector == rr)
    if (length(idx) == 0) next
    m <- circlize(axes$r_mid[idx], rep(0.5, length(idx)), sector.index = rr, track.index = 1)
    theta_r[idx] <- m[, 1]
  }

  keep_r <- rep(TRUE, nrow(axes))
  if (LABEL_MIN_DEG > 0) {
    for (rr in unique(axes$r_sector)) {
      idx <- which(axes$r_sector == rr)
      keep_r[idx] <- thin_by_min_deg(theta_r[idx], LABEL_MIN_DEG)
    }
  }

  op <- par(xpd = NA)

  for (k in seq_len(nrow(ligand_labels))) {
    if (!keep_lig[k] || !is.finite(theta_lig[k])) next

    ang <- theta_lig[k] * pi / 180
    x <- LABEL_RADIUS * cos(ang)
    y <- LABEL_RADIUS * sin(ang)

    srt <- theta_lig[k]
    adjx <- 0
    if (LABEL_NICEFACING && theta_lig[k] > 90 && theta_lig[k] < 270) {
      srt <- theta_lig[k] + 180
      adjx <- 1
    }

    text(
      x = x, y = y,
      labels = ligand_labels$ligand[k],
      srt = srt,
      adj = c(adjx, 0.5),
      cex = LABEL_CEX
    )
  }

  for (i in seq_len(nrow(axes))) {
    if (keep_r[i] && is.finite(theta_r[i])) {
      ang <- theta_r[i] * pi / 180
      x <- LABEL_RADIUS * cos(ang)
      y <- LABEL_RADIUS * sin(ang)

      srt <- theta_r[i]
      adjx <- 0
      if (LABEL_NICEFACING && theta_r[i] > 90 && theta_r[i] < 270) {
        srt <- theta_r[i] + 180
        adjx <- 1
      }

      text(x = x, y = y, labels = axes$target[i], srt = srt, adj = c(adjx, 0.5), cex = LABEL_CEX)
    }
  }

  par(op)

  legend(
    "topleft",
    legend = ct_names,
    fill = ct_cols[ct_names],
    border = NA,
    cex = 0.75,
    bty = "n",
    title = "Cell types"
  )

  legend(
    "bottomleft",
    legend = c("TIL_high axes", "TIL_low axes"),
    fill = c("#D55E00", "#0072B2"),
    border = NA,
    cex = 0.95,
    bty = "n",
    title = "Axis direction"
  )

  title(
    main = paste0(
      "Selected axes circos | ", zone_dir,
      " | edges=", nrow(axes),
      " | width ~ ligand_activity, color ~ rank_score"
    ),
    cex.main = TITLE_CEX
  )

  circos.clear()
}

# Write PNG once
png(OUT_PNG, width = PNG_SIZE, height = PNG_SIZE, res = PNG_RES, type = "cairo")
plot_one()
dev.off()

OUT_PDF <- sub("\\.png$", ".pdf", OUT_PNG)

pdf(OUT_PDF, width = PNG_SIZE / PNG_RES, height = PNG_SIZE / PNG_RES, useDingbats = FALSE)
plot_one()
dev.off()

cat("Wrote:", OUT_PNG, "\n")
cat("Wrote:", OUT_PDF, "\n")
cat("Edges:", nrow(axes), "\n")
cat("ligand_activity range:", min(axes$ligand_activity), "to", max(axes$ligand_activity), "\n")
cat("rank_score range:", min(axes$rank_score), "to", max(axes$rank_score), "\n")
cat("rank_score clip quantiles:", COLOR_CLIP_LO, "to", COLOR_CLIP_HI, "\n")
cat("Label radius:", LABEL_RADIUS, " Label cex:", LABEL_CEX, " Label min deg:", LABEL_MIN_DEG, "\n")
as.data.frame(axes %>% select(-row_id))
}
'''

_defined = False


def circos_zone(axes_csv, out_png, zone_dir: str | None = None) -> pd.DataFrame:
    """Draw the circos plot of one zone's axes_selected.csv to out_png and the .pdf next to it.

    zone_dir is the zone label of the title (default: the folder name of axes_csv, as in the
    notebook). Returns the plotted axes with cord widths, colours and positions.
    """
    global _defined
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    if not _defined:
        ro.r(CIRCOS_R)
        _defined = True
    axes_csv = Path(axes_csv)
    zone_dir = zone_dir or axes_csv.parent.name
    r_axes = ro.r["circos_zone"](str(axes_csv), zone_dir, str(out_png))
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(r_axes)
