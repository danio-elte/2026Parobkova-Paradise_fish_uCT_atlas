#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 2 - Quantification of the segmentation in adult paradise fish
#            (3 males, 3 females; whole-body micro-CT; relative volumes)
#
# Usage:   Rscript make_figure2.R [data_dir] [out_dir]        (defaults: "data", "figures")
#
# Input (data_dir):
#   Supplementary_Table_1.csv           long format, one row per specimen x structure; includes the
#                                       body measures and relative_volume_pct (used as supplied)
#   pfish_uCT_segmentation_colors.csv   structure -> hex colour (as in the 3D renderings)
#
# Output (out_dir):
#   Figure2.pdf / .png / .tiff          the figure
#   FigureS_facets_volume.*             one panel per structure, ABSOLUTE volumes (supplementary)
#   Figure2_stats_organs.csv            male vs female tests per structure (relative and absolute volumes)
#   Figure2_stats_bodysize.csv          male vs female tests for the body-size measures
#
# Required packages: readr, dplyr, tidyr, ggplot2, patchwork, scales, stringr (ragg optional)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(scales); library(stringr)
})

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- settings ------------------------------------------------------------------
P_ADJUST    <- "none"  # significance marks: raw Welch p-values ("none") or Benjamini-Hochberg ("BH")
ALPHA       <- 0.05
MAKE_FACETS <- TRUE    # also write the one-panel-per-structure supplementary figure

# structures per row of the figure (order as in Table 5 of the manuscript)
row_nerv <- c("telencephalon", "olfactory bulb", "hypothalamus", "pituitary", "optic tectum",
              "cerebellum", "brain stem",                  # Brain
              "spinal cord", "optic nerve",                # Nervous system
              "retina", "choroid", "sclera", "lens")       # Eye
row_rest <- c("heart", "sinus venosus", "bulbus arteriosus",    # Cardiovascular
              "liver", "stomach", "intestines",                 # Digestive
              "swimming bladder", "muscle", "testis", "ovary")  # Other (gonads included)
sex_limited <- c("testis", "ovary")
struct_all  <- c(row_nerv, row_rest)

# ---- read data -----------------------------------------------------------------
tab  <- read_csv(file.path(data_dir, "2026Parobkova-pfish_uCT-Supplementary_Table_1_20261121"), show_col_types = FALSE)
cols <- read_csv(file.path(data_dir, "pfish_uCT_segmentation_colors.csv"), show_col_types = FALSE)

pal <- setNames(cols$hex, cols$structure)
stopifnot(all(struct_all %in% names(pal)), all(struct_all %in% tab$structure))
missing_colour <- cols$structure[grepl("PLACEHOLDER", cols$source) & cols$structure %in% struct_all]
if (length(missing_colour))
  message("NOTE: placeholder colours are used for: ", paste(missing_colour, collapse = ", "),
          " (edit pfish_uCT_segmentation_colors.csv)")

sex_levels <- c("male", "female")
tab <- tab %>% mutate(sex = factor(sex, levels = sex_levels),
                      organ_system = ifelse(organ_system == "Gonad", "Other", organ_system))  # gonads -> "Other"

vol <- tab %>%
  filter(segmentation_status == "measured", structure %in% struct_all) %>%
  mutate(structure = factor(structure, levels = struct_all))
system_of <- tab %>% distinct(structure, organ_system) %>% { setNames(.$organ_system, .$structure) }

body <- tab %>%
  distinct(specimen_id, sex, standard_length_mm, head_length_mm, body_depth_mm, total_body_volume_mm3) %>%
  pivot_longer(-c(specimen_id, sex), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
           levels = c("standard_length_mm", "head_length_mm", "body_depth_mm", "total_body_volume_mm3"),
           labels = c('atop("Standard length","(mm)")', 'atop("Head length","(mm)")',
                      'atop("Body depth","(mm)")', 'atop("Total volume","(mm"^3*")")')))

# ---- statistics: male vs female (Welch t-test; log10 for organ volumes) ---------------
# n = 3 per sex: exploratory only. Sex-limited structures (testis, ovary) are not compared.
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                      ifelse(p < 0.05, "*", ""))))
welch_p <- function(y, sex, log = FALSE) {
  if (log) y <- log10(y)
  tryCatch(t.test(y ~ sex)$p.value, error = function(e) NA_real_)
}
stat_org <- vol %>% filter(!structure %in% sex_limited) %>%
  group_by(structure) %>%
  summarise(mean_rel_male = mean(relative_volume_pct[sex == "male"]),
            mean_rel_female = mean(relative_volume_pct[sex == "female"]),
            p_relative = welch_p(relative_volume_pct, sex, log = TRUE),
            mean_abs_male = mean(volume_mm3[sex == "male"]),
            mean_abs_female = mean(volume_mm3[sex == "female"]),
            p_absolute = welch_p(volume_mm3, sex, log = TRUE),
            .groups = "drop") %>%
  mutate(structure = as.character(structure),
         p_relative_BH = p.adjust(p_relative, "BH"),
         p_absolute_BH = p.adjust(p_absolute, "BH"),
         p_used = if (P_ADJUST == "BH") p_relative_BH else p_relative,
         mark = ifelse(p_used < ALPHA, stars(p_used), ""))
stat_body <- body %>% group_by(measure) %>%
  summarise(p = welch_p(value, sex), .groups = "drop") %>%
  mutate(mark = stars(p))

write_csv(select(stat_org, -p_used), file.path(out_dir, "Figure2_stats_organs.csv"))
write_csv(transmute(stat_body, measure = c("standard length", "head length", "body depth", "total volume"), p, mark),
          file.path(out_dir, "Figure2_stats_bodysize.csv"))
message("Organ comparisons with p < ", ALPHA, " (", ifelse(P_ADJUST == "BH", "BH-adjusted", "uncorrected"), "): ",
        paste(stat_org$structure[stat_org$mark != ""], collapse = ", "))
message("Organ comparisons significant after BH correction: ", sum(stat_org$p_relative_BH < ALPHA, na.rm = TRUE))

# ---- shared look ---------------------------------------------------------------------
shp <- c(male = 21, female = 24)          # filled circle = male, filled triangle = female

theme_fig <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(), strip.text = element_text(size = 7.5),
        legend.position = "none", legend.key.size = unit(3.5, "mm"),
        plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0.005, 0.995),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))

# males left, females right of the tick (sex-limited structures stay centred)
side <- function(df, off) ifelse(df$structure %in% sex_limited, 0, ifelse(df$sex == "male", -off, off))

# ---- dot plot of relative volumes for a set of structures --------------------------------
structure_dotplot <- function(order, tag, legend = FALSE) {
  df <- vol %>% filter(structure %in% order) %>%
    mutate(structure = factor(as.character(structure), levels = order))
  sys_tbl <- tibble(structure = order, system = unname(system_of[order])) %>%
    mutate(i = row_number()) %>% group_by(system) %>%
    summarise(start = min(i), end = max(i), .groups = "drop") %>%
    arrange(start) %>% mutate(mid = (start + end) / 2)
  set.seed(1)
  d  <- df %>% mutate(xc = as.numeric(structure), xp = xc + side(., 0.2) + runif(n(), -0.07, 0.07))
  mu <- df %>% group_by(structure, sex) %>%
    summarise(m = mean(relative_volume_pct), .groups = "drop") %>%
    mutate(xc = as.numeric(structure), x0 = xc + side(., 0.2) - 0.16, x1 = xc + side(., 0.2) + 0.16)
  sg <- df %>% group_by(structure) %>% summarise(top = max(relative_volume_pct), .groups = "drop") %>%
    mutate(structure = as.character(structure)) %>%
    inner_join(filter(stat_org, mark != ""), by = "structure") %>%
    mutate(xc = match(structure, order), y = top * 1.35)   # bracket above the higher group
  ggplot(d, aes(xp, relative_volume_pct)) +
    geom_vline(xintercept = head(sys_tbl$end, -1) + 0.5, colour = "grey80",
               linewidth = 0.25, linetype = "22") +
    geom_segment(data = mu, aes(x = x0, xend = x1, y = m, yend = m), inherit.aes = FALSE, linewidth = 0.5) +
    geom_point(aes(fill = structure, shape = sex), size = 2.1, stroke = 0.3, colour = "black") +
    geom_segment(data = sg, aes(x = xc - 0.2, xend = xc + 0.2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.35) +
    geom_text(data = sg, aes(x = xc, y = y, label = mark), inherit.aes = FALSE, vjust = 0.15, size = 4.2) +
    annotate("text", x = sys_tbl$mid, y = Inf, label = sys_tbl$system, vjust = 1.4, size = 2.5, colour = "grey30") +
    scale_x_continuous(breaks = seq_along(order), labels = order,
                       limits = c(0.5, length(order) + 0.5), expand = c(0, 0)) +
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                  labels = label_number(big.mark = ",", drop0trailing = TRUE),
                  expand = expansion(mult = c(0.05, 0.30))) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = shp, breaks = c("female", "male"),            # legend: triangle on top
                       guide = guide_legend(title = NULL, ncol = 1,
                                            override.aes = list(fill = "grey55", size = 2.4, stroke = 0.3))) +
    labs(x = NULL, y = "Relative volume (%)", tag = tag) +
    theme_fig +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = if (legend) "right" else "none",
          legend.direction = "vertical", legend.text = element_text(size = 8),
          legend.box.background = element_rect(colour = "grey60", linewidth = 0.3),
          legend.box.margin = margin(4, 8, 4, 8))
}

# ---- panel a: body size -----------------------------------------------------------------------
pa <- {
  set.seed(1)
  d  <- body %>% mutate(xc = as.numeric(sex), xp = xc + runif(n(), -0.08, 0.08))
  mu <- d %>% group_by(measure, sex) %>% summarise(m = mean(value), xc = first(xc), .groups = "drop")
  sg <- d %>% group_by(measure) %>% summarise(top = max(value), .groups = "drop") %>%
    inner_join(filter(stat_body, mark != ""), by = "measure") %>% mutate(y = top * 1.07)
  ggplot(d, aes(xp, value)) +
    geom_segment(data = mu, aes(x = xc - 0.28, xend = xc + 0.28, y = m, yend = m),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = sex), fill = "grey55", size = 2.1, stroke = 0.3, colour = "black") +
    geom_segment(data = sg, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
    geom_text(data = sg, aes(x = 1.5, y = y, label = mark), inherit.aes = FALSE, vjust = 0.15, size = 4.2) +
    facet_wrap(~ measure, nrow = 1, scales = "free_y", labeller = label_parsed) +
    scale_x_continuous(breaks = 1:2, labels = sex_levels, limits = c(0.4, 2.6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = label_comma()) +
    expand_limits(y = 0) +
    scale_shape_manual(values = shp, guide = "none") +
    labs(x = NULL, y = NULL, tag = "a") +
    theme_fig
}

# ---- panels b (brain / nervous system / eye) and c (the rest), legend box next to c --------------
pb <- structure_dotplot(row_nerv, "b")
pc <- structure_dotplot(row_rest, "c")
gt  <- ggplotGrob(structure_dotplot(row_rest, "", legend = TRUE))
leg <- wrap_elements(full = gt$grobs[[which(vapply(gt$grobs, function(x)
  grepl("guide-box", x$name) && !inherits(x, "zeroGrob"), logical(1)))[1]]])

note <- paste0("Mean = horizontal bar; n = 3 per sex (testis and ovary: n = 3 of one sex each). ",
  "Male vs female: Welch's t-test (body size: raw values; organs: log10 relative volume); ",
  "* P < 0.05, ** P < 0.01",
  if (P_ADJUST == "none")
    sprintf(" (uncorrected for multiple testing; %d of %d organ comparisons remain significant after Benjamini-Hochberg correction).",
            sum(stat_org$p_relative_BH < ALPHA, na.rm = TRUE), nrow(stat_org))
  else " (Benjamini-Hochberg adjusted).")

fig2 <- (pa / pb / ((pc | leg) + plot_layout(widths = c(10, 3)))) +
  plot_layout(heights = c(1, 1.15, 1.15)) +
  plot_annotation(caption = str_wrap(note, 125),
                  theme = theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30")))

save_all <- function(plot, stem, w_mm, h_mm) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = w_mm, height = h_mm, units = "mm",
         device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(file.path(out_dir, paste0(stem, ".png")),  plot, width = w_mm, height = h_mm, units = "mm",
           dpi = 300, device = ragg::agg_png, bg = "white")
    ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot, width = w_mm, height = h_mm, units = "mm",
           dpi = 300, device = ragg::agg_tiff, bg = "white", compression = "lzw")
  } else {
    ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = w_mm, height = h_mm, units = "mm",
           dpi = 300, bg = "white")
  }
}
save_all(fig2, "Figure2", 180, 190)

# ---- supplementary: one panel per structure, absolute volumes (free y axis) -------------------
if (MAKE_FACETS) {
  set.seed(1)
  d  <- vol %>% mutate(xc = as.numeric(sex), xp = xc + runif(n(), -0.1, 0.1))
  mu <- d %>% group_by(structure, sex) %>% summarise(m = mean(volume_mm3), xc = first(xc), .groups = "drop")
  pf <- ggplot(d, aes(xp, volume_mm3)) +
    geom_segment(data = mu, aes(x = xc - 0.3, xend = xc + 0.3, y = m, yend = m),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_point(aes(fill = structure, shape = sex), size = 1.9, stroke = 0.3, colour = "black") +
    facet_wrap(~ structure, scales = "free_y", ncol = 6, labeller = label_wrap_gen(14)) +
    scale_x_continuous(breaks = 1:2, labels = c("M", "F"), limits = c(0.4, 2.6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + expand_limits(y = 0) +   # axes start at zero
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = shp, guide = "none") +
    labs(x = NULL, y = expression("Volume (mm"^3*")")) + theme_fig
  save_all(pf, "FigureS_facets_volume", 180, 170)
}

message("Done. Files written to: ", normalizePath(out_dir))
