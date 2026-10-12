library(tidyverse)
library(patchwork)
library(grid)

read_clean <- function(file) {
  x <- read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  x[, !is.na(names(x)) & names(x) != "", drop = FALSE]
}

base_family <- "sans"

trait_levels <- c(
  "Overall CH", "DNMT3A-CH", "TET2-CH",
  "ASXL1-CH", "mLOY", "MPN"
)

trait_cols <- c(
  "Overall CH" = "#0072B2",
  "DNMT3A-CH" = "#D55E00",
  "TET2-CH" = "#4D4D4D",
  "ASXL1-CH" = "#009E73",
  "mLOY" = "#CC79A7",
  "MPN" = "#E69F00"
)

theme_forest <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = base_family, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size + 2),
      axis.title = element_text(size = base_size + 4),
      strip.background = element_rect(
        fill = "white", colour = "black", linewidth = 0.5
      ),
      strip.text = element_text(
        family = base_family, colour = "black",
        size = base_size + 3, face = "plain",
        margin = margin(2, 2, 2, 2)
      ),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.5
      ),
      panel.grid.major.x = element_line(
        colour = "grey88", linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      plot.margin = margin(3, 8, 3, 8)
    )
}

infile <- "ed_fig1.txt"

dat <- read_clean(infile) |>
  mutate(
    trait = factor(trait, levels = trait_levels),
    # Reciprocal OR represents the allele opposite to original C.
    # CI endpoints must be inverted AND swapped.
    original_or = as.numeric(or),
    original_l95cl = as.numeric(l95cl),
    original_u95cl = as.numeric(u95cl),
    or = 1 / original_or,
    l95cl = 1 / original_u95cl,
    u95cl = 1 / original_l95cl,
    y_pos = length(trait_levels) - as.numeric(trait) + 1,
    # Original effect allele C is omitted from the title because the
    # reciprocal ORs refer to the other allele (not supplied in file).
    snp_label = "'rs17758695 - T - ' * italic('BCL2')"
  )

stopifnot(
  all(!is.na(dat$trait)),
  all(is.finite(dat$or)),
  all(is.finite(dat$l95cl)),
  all(is.finite(dat$u95cl)),
  all(dat$l95cl <= dat$or & dat$or <= dat$u95cl),
  all(dat$l95cl > 0)
)

p_forest <- ggplot(dat, aes(x = or, y = y_pos, colour = trait)) +
  geom_vline(
    xintercept = 1, linewidth = 0.3,
    linetype = "dashed", colour = "black"
  ) +
  geom_segment(
    aes(x = l95cl, xend = u95cl, yend = y_pos),
    linewidth = 0.55, lineend = "butt"
  ) +
  geom_segment(
    aes(x = l95cl, xend = l95cl,
        y = y_pos - 0.08, yend = y_pos + 0.08),
    linewidth = 0.4, lineend = "butt"
  ) +
  geom_segment(
    aes(x = u95cl, xend = u95cl,
        y = y_pos - 0.08, yend = y_pos + 0.08),
    linewidth = 0.4, lineend = "butt"
  ) +
  geom_point(size = 2.1) +
  facet_wrap(
    ~ snp_label,
    ncol = 1,
    strip.position = "top",
    labeller = label_parsed
  ) +
  scale_colour_manual(values = trait_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(trait_levels),
    labels = rev(trait_levels),
    limits = c(0.5, length(trait_levels) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_log10(
    breaks = c(0.4, 0.6, 0.8, 1, 1.2),
    labels = c("0.4", "0.6", "0.8", "1", "1.2"),
    limits = c(0.4, 1.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Odds ratio (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_forest() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

ggsave(
  "ed_fig1.pdf", p_forest,
  width = 7.2, height = 3.2,
  units = "in", device = cairo_pdf
)

ggsave(
  "ed_fig1.png", p_forest,
  width = 7.2, height = 3.2,
  units = "in", dpi = 600
)

p_forest
