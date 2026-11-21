#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 3 (quantification part) - organ volumes across the post-metamorphic juvenile
# stages (one specimen per stage -> descriptive only, no statistics).
#
# Usage: Rscript make_figure3_quant.R [data_dir] [out_dir]   (defaults: data, figures)
# Input (data_dir): Supplementary_Table_2.csv  (from make_supp_table2.R)
#                   pfish_uCT_segmentation_colors.csv  (same colour key as the adults)
# Output (out_dir):
#   Figure3_quant.pdf / .png / .tiff          panel b of Figure 3: absolute volume (mm3, log10) per
#                                             structure, one panel per organ system
#   FigureS_juvenile_relative.pdf/.png/.tiff  supplementary figure: relative volume (% of total body volume),
#                                             calculated by the script from volume_mm3 / total_body_volume_mm3
#   (panel a of Figure 3 = the specimen renderings, assembled separately)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales); library(ggrepel)
})
args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TEMPORARY   <- FALSE  # TRUE adds a red "temporary figure" caption
SHOW_LABELS <- TRUE   # structure names written next to the last data point of each line

t2   <- read_csv(file.path(data_dir, "2026Parobkova-pfish_uCT-Supplementary_Table_2_20261121"), show_col_types = FALSE)
cols <- read_csv(file.path(data_dir, "pfish_uCT_segmentation_colors.csv"), show_col_types = FALSE)
pal  <- setNames(cols$hex, cols$structure)

# relative volume is calculated here: structure volume / total body volume x 100
d <- t2 %>% filter(segmentation_status == "measured", !is.na(volume_mm3)) %>%
  mutate(relative_volume_pct = 100 * volume_mm3 / total_body_volume_mm3)
if (nrow(d) == 0) stop("No measured volumes in Supplementary_Table_2.csv yet - fill the wide sheet first.")
missing_col <- setdiff(unique(d$structure), names(pal))
if (length(missing_col)) {
  message("NOTE: no colour for: ", paste(missing_col, collapse = ", "), " - drawn in grey.")
  pal <- c(pal, setNames(rep("grey60", length(missing_col)), missing_col))
}
if (any(grepl("PLACEHOLDER", cols$source[cols$structure %in% unique(d$structure)])))
  message("NOTE: placeholder colours are used for structures without a supplied colour.")

sys_levels <- c("Brain", "Nervous system", "Eye", "Cardiovascular", "Digestive", "Other")
d <- d %>% mutate(organ_system = factor(organ_system, levels = sys_levels),
                  structure = factor(structure, levels = names(pal)))

# x axis: one tick per stage, labelled with age and standard length
st <- d %>% distinct(specimen_id, age_wpf, standard_length_mm) %>% arrange(age_wpf) %>%
  mutate(x = row_number(),
         lab = paste0(age_wpf, "\n", format(standard_length_mm, nsmall = 1, trim = TRUE)))
xmax <- max(st$x)
d <- d %>% left_join(select(st, specimen_id, x), by = "specimen_id")

theme_fig <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(colour = "black"), axis.text.x = element_text(size = 7, lineheight = 0.95), axis.title.x = element_text(size = 6.5, colour = "grey30"), strip.background = element_blank(),
        strip.text = element_text(size = 7.5), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0.005, 0.995),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))

stage_panel <- function(df, yvar, ylab, tag, free_y = TRUE) {
  lab <- df %>% group_by(structure) %>% slice_max(x, n = 1, with_ties = FALSE) %>% ungroup()
  ggplot(df, aes(x, .data[[yvar]], colour = structure, fill = structure, group = structure)) +
    geom_line(linewidth = 0.45) +
    geom_point(shape = 21, size = 2, stroke = 0.3, colour = "black") +
    { if (SHOW_LABELS) geom_text_repel(data = lab, aes(label = structure), colour = "grey20", hjust = 0, direction = "y",
                                       nudge_x = 0.25, xlim = c(xmax + 0.2, Inf), size = 2.3,
                                       segment.size = 0.2, segment.colour = "grey70",
                                       min.segment.length = 0.1, box.padding = 0.08, seed = 1,
                                       show.legend = FALSE) } +
    facet_wrap(~ droplevels(organ_system), ncol = ifelse(nlevels(droplevels(df$organ_system)) <= 4, 2, 3),
               scales = if (free_y) "free_y" else "fixed") +
    scale_x_continuous(breaks = st$x, labels = st$lab, limits = c(0.7, xmax + 3), expand = c(0, 0)) +
    scale_y_log10(labels = label_number(big.mark = ",", drop0trailing = TRUE),
                  expand = expansion(mult = 0.12)) +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    labs(x = "Age (wpf) / standard length (mm)", y = ylab, tag = tag) + theme_fig
}

temp_caption <- function(p) if (!TEMPORARY) p else p + plot_annotation(
  caption = "TEMPORARY FIGURE - juvenile data are provisional.",
  theme = theme(plot.caption = element_text(size = 6.5, colour = "#B22222", hjust = 0)))

ragg_ok <- requireNamespace("ragg", quietly = TRUE)
save_all <- function(plot, stem, w_mm, h_mm) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = w_mm, height = h_mm, units = "mm",
         device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = w_mm, height = h_mm, units = "mm", dpi = 300,
         bg = "white", device = if (ragg_ok) ragg::agg_png else "png")
  if (ragg_ok) ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot, width = w_mm, height = h_mm, units = "mm",
                      dpi = 300, bg = "white", device = ragg::agg_tiff, compression = "lzw")
}

# main figure: absolute volumes (panel b of Figure 3)
pb <- stage_panel(d, "volume_mm3", expression("Volume (mm"^3*")"), "b")
save_all(temp_caption(pb), "Figure3_quant", 180, 110)

# supplementary figure: relative volumes (common y scale: no exaggeration of small changes)
dr <- d %>% filter(!is.na(relative_volume_pct))
if (nrow(dr) > 0) {
  ps <- stage_panel(dr, "relative_volume_pct", "Relative volume (%)", "", free_y = FALSE)
  save_all(temp_caption(ps), "FigureS_juvenile_relative", 180, 110)
} else {
  message("NOTE: total_body_volume_mm3 is missing in Supplementary Table 2 - relative-volume figure skipped.")
}

message("Done. Files written to: ", normalizePath(out_dir))
