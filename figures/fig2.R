# ============================================================
# Combined and cleaned plotting script for fig2a, fig2b, fig2
# ============================================================

# ---- Packages ----
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(ggrepel)

# ---- GRCh37 / hg19 chromosome lengths (autosomes only) ----
chr_lengths_grch37 <- tibble(
    chr_num = 1:22,
    chr_len = c(
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
        115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566
    )
) %>%
    mutate(
        offset = lag(cumsum(chr_len), default = 0),
        center = offset + chr_len / 2,
        xmin = offset,
        xmax = offset + chr_len,
        stripe = if_else(chr_num %% 2 == 0, "even", "odd")
    )

# ---- Helper: prepare data ----
prepare_plot_data <- function(file, driver_levels) {
    df <- read_tsv(file, show_col_types = FALSE)

    df2 <- df %>%
        mutate(
            chr = as.character(chr),
            chr_num = case_when(
                chr %in% c("X", "x") ~ 23,
                chr %in% c("Y", "y") ~ 24,
                chr %in% c("MT", "M", "mt", "m") ~ 25,
                TRUE ~ suppressWarnings(as.numeric(chr))
            ),
            label = OpenTargets_V2G,
            ngenes = as.character(ngenes),
            status = as.character(status)
        ) %>%
        filter(
            !is.na(chr_num),
            chr_num >= 1,
            chr_num <= 22,
            !is.na(pos),
            !is.na(p_original),
            !is.na(p_sfFDR),
            !is.na(label),
            label != "",
            !is.na(ngenes),
            ngenes != ""
        ) %>%
        arrange(chr_num, pos)

    chr_sizes <- chr_lengths_grch37

    df2 <- df2 %>%
        left_join(chr_sizes %>% select(chr_num, offset, chr_len), by = "chr_num") %>%
        filter(pos <= chr_len) %>%
        mutate(x = pos + offset)

    genome_span <- sum(chr_sizes$chr_len)
    offset_amt <- genome_span * 0.0028

    plot_df <- df2 %>%
        select(CH_phenotype, rsID, chr, chr_num, pos, x, label, ngenes, status, p_original, p_sfFDR) %>%
        pivot_longer(
            cols = c(p_original, p_sfFDR),
            names_to = "p_type",
            values_to = "p_value"
        ) %>%
        mutate(
            neglog10p = -log10(p_value),
            x_plot = if_else(p_type == "p_original", x - offset_amt, x + offset_amt),
            p_type = factor(p_type, levels = c("p_original", "p_sfFDR")),
            ngenes = factor(ngenes, levels = driver_levels)
        )

    label_df <- plot_df %>%
        filter(!is.na(label), label != "") %>%
        group_by(label) %>%
        slice_max(order_by = neglog10p, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(chr_num, x)

    ymax_signal <- max(plot_df$neglog10p, na.rm = TRUE)

    close_threshold_bp <- 20000000
    base_nudge_x <- genome_span * 0.0085

    label_df <- label_df %>%
        group_by(chr_num) %>%
        arrange(x, .by_group = TRUE) %>%
        mutate(
            cluster_id = cumsum(if_else(row_number() == 1, 1L, if_else((x - lag(x)) > close_threshold_bp, 1L, 0L))),
            cluster_n = dplyr::n(),
            y_anchor = neglog10p,
            y_label_base = ymax_signal + 7,
            label_expr = case_when(
                status == "Replicated" ~ paste0("bolditalic('", label, "')^'\\u2020'"),
                status == "Meta" ~ paste0("bolditalic('", label, "')^'\\u2021'"),
                TRUE ~ paste0("bolditalic('", label, "')")
            )
        ) %>%
        group_by(chr_num, cluster_id, .add = TRUE) %>%
        mutate(
            n_in_cluster = n(),
            idx_in_cluster = row_number(),
            nudge_x = if_else(
                n_in_cluster == 1,
                0,
                (idx_in_cluster - (n_in_cluster + 1) / 2) * base_nudge_x
            )
        ) %>%
        ungroup()

    list(
        plot_df = plot_df,
        label_df = label_df,
        chr_sizes = chr_sizes,
        stripe_df = chr_sizes,
        ymax_signal = ymax_signal,
        genome_span = genome_span
    )
}

# ---- Helper: make plot ----
make_plot <- function(
        prep,
        panel_title,
        driver_levels,
        shape_vals,
        driver_label_expr,
        dodge_x_axis = FALSE
) {
    plot_df <- prep$plot_df
    label_df <- prep$label_df
    chr_sizes <- prep$chr_sizes
    stripe_df <- prep$stripe_df
    ymax_signal <- prep$ymax_signal
    genome_span <- prep$genome_span

    thr_genome <- -log10(5e-8)
    thr_suggestive <- -log10(5e-4)

    y_top <- ymax_signal + 14
    base_breaks <- pretty(c(0, y_top), n = 7)
    y_breaks <- sort(unique(c(base_breaks, thr_suggestive, thr_genome)))

    y_labels <- sapply(
        y_breaks,
        function(b) {
            if (abs(b - thr_genome) < 1e-8) {
                format(round(thr_genome, 2), nsmall = 2, trim = TRUE)
            } else if (abs(b - thr_suggestive) < 1e-8) {
                format(round(thr_suggestive, 2), nsmall = 2, trim = TRUE)
            } else {
                if (abs(b - round(b)) < 1e-8) {
                    format(round(b, 0), trim = TRUE, nsmall = 0)
                } else {
                    format(round(b, 1), trim = TRUE, nsmall = 1)
                }
            }
        }
    )

    legend_x <- min(chr_sizes$xmin, na.rm = TRUE) + genome_span * 0.70
    legend_y1 <- y_top - 2.0
    legend_y2 <- y_top - 3.8

    x_scale <- scale_x_continuous(
        breaks = chr_sizes$center,
        labels = chr_sizes$chr_num,
        expand = expansion(mult = c(0.012, 0.02)),
        guide = if (dodge_x_axis) guide_axis(n.dodge = 2) else waiver()
    )

    ggplot() +
        geom_rect(
            data = stripe_df %>% filter(stripe == "even"),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey95",
            color = NA
        ) +
        geom_segment(
            data = plot_df,
            aes(x = x_plot, xend = x_plot, y = 0, yend = neglog10p, color = p_type),
            linewidth = 0.95,
            alpha = 1
        ) +
        geom_point(
            data = plot_df,
            aes(x = x_plot, y = neglog10p, fill = p_type, shape = ngenes),
            color = "black",
            stroke = 0.55,
            size = 6.0,
            alpha = 0.8
        ) +
        geom_hline(
            yintercept = thr_genome,
            linetype = "dashed",
            linewidth = 0.6,
            color = "firebrick3"
        ) +
        geom_hline(
            yintercept = thr_suggestive,
            linetype = "dashed",
            linewidth = 0.6,
            color = "grey20"
        ) +
        annotate(
            "text",
            x = min(chr_sizes$xmin, na.rm = TRUE) + genome_span * 0.02,
            hjust = 0,
            y = y_top - 1,
            label = panel_title,
            size = 10,
            fontface = "plain",
            family = "sans",
            color = "black",
            alpha = 0.9
        ) +
        annotate(
            "text",
            x = legend_x,
            y = legend_y1,
            label = "\u2020 Replicated",
            hjust = 0,
            vjust = 1,
            size = 8,
            fontface = "plain",
            family = "sans",
            color = "black"
        ) +
        annotate(
            "text",
            x = legend_x,
            y = legend_y2,
            label = "\u2021 GWS in meta",
            hjust = 0,
            vjust = 1,
            size = 8,
            fontface = "plain",
            family = "sans",
            color = "black"
        ) +
        geom_text_repel(
            data = label_df,
            aes(x = x, y = y_anchor, label = label_expr),
            inherit.aes = FALSE,
            parse = TRUE,
            angle = 90,
            nudge_x = label_df$nudge_x,
            family = "sans",
            size = 7.2,
            color = "black",
            direction = "both",
            nudge_y = 7,
            force = 2.5,
            force_pull = 0.8,
            box.padding = 0.5,
            point.padding = 0.25,
            min.segment.length = 0,
            segment.color = "black",
            segment.linetype = "dashed",
            segment.size = 0.85,
            max.overlaps = Inf,
            seed = 123
        ) +
        x_scale +
        scale_y_continuous(
            limits = c(0, y_top),
            breaks = y_breaks,
            labels = y_labels,
            expand = expansion(mult = c(0, 0))
        ) +
        scale_color_manual(
            values = c(
                p_original = "#2F5D8A",
                p_sfFDR = "#D96B27"
            ),
            labels = c("Original", "sfFDR"),
            name = NULL
        ) +
        scale_fill_manual(
            values = c(
                p_original = "#2F5D8A",
                p_sfFDR = "#D96B27"
            ),
            guide = "none"
        ) +
        scale_shape_manual(
            values = shape_vals,
            breaks = driver_levels,
            labels = driver_label_expr,
            name = "CH drivers"
        ) +
        labs(
            x = "Chromosome",
            y = expression(-log[10](italic(P)))
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 20) +
        theme(
            axis.line.x = element_line(linewidth = 1.1, color = "black"),
            axis.line.y = element_line(linewidth = 1.1, color = "black"),
            axis.ticks = element_line(linewidth = 1.0, color = "black"),
            axis.ticks.length = unit(3.2, "mm"),
            axis.text.x = element_text(size = 20, color = "black", face = "plain"),
            axis.text.y = element_text(size = 19, color = "black", face = "plain"),
            axis.title.x = element_text(size = 22, face = "plain", margin = margin(t = 14)),
            axis.title.y = element_text(size = 22, face = "plain", margin = margin(r = 16)),

            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.box.just = "center",
            legend.justification = "center",

            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
            legend.margin = margin(6, 12, 6, 12),
            legend.box.margin = margin(2, 2, 2, 2),

            legend.key.width = unit(1.1, "cm"),
            legend.key.height = unit(0.75, "cm"),
            legend.spacing.x = unit(1.5, "cm"),
            legend.spacing.y = unit(0.15, "cm"),

            legend.title = element_text(size = 18, face = "bold", color = "black"),
            legend.text = element_text(size = 17, face = "bold", color = "black"),

            plot.margin = margin(52, 44, 38, 20)
        ) +
        guides(
            color = guide_legend(
                order = 1,
                nrow = 1,
                byrow = TRUE,
                override.aes = list(
                    linewidth = 2.5,
                    alpha = 1,
                    shape = NA
                )
            ),
            shape = guide_legend(
                order = 2,
                nrow = 1,
                byrow = TRUE,
                override.aes = list(
                    size = 5.6,
                    fill = "grey70",
                    color = "black",
                    stroke = 0.55
                )
            )
        )
}

# ---- Panel A settings ----
driver_levels_a <- c(
    "DNMT3A",
    "DNMT3A+TET2",
    "DNMT3A+ASXL1",
    "DNMT3A+TET2+ASXL1"
)

shape_vals_a <- c(
    "DNMT3A" = 21,
    "DNMT3A+TET2" = 22,
    "DNMT3A+ASXL1" = 24,
    "DNMT3A+TET2+ASXL1" = 23
)

driver_label_expr_a <- c(
    expression(italic(DNMT3A)),
    expression(italic(DNMT3A) * "+" * italic(TET2)),
    expression(italic(DNMT3A) * "+" * italic(ASXL1)),
    expression(italic(DNMT3A) * "+" * italic(TET2) * "+" * italic(ASXL1))
)

# ---- Panel B settings ----
driver_levels_b <- c("DNMT3A", "TET2", "ASXL1")

shape_vals_b <- c(
    "DNMT3A" = 21,
    "TET2" = 22,
    "ASXL1" = 24
)

driver_label_expr_b <- c(
    expression(italic(DNMT3A)),
    expression(italic(TET2)),
    expression(italic(ASXL1))
)

# ---- Build plots ----
prep_a <- prepare_plot_data("fig2a.txt", driver_levels_a)
prep_b <- prepare_plot_data("fig2b.txt", driver_levels_b)

p1 <- make_plot(
    prep = prep_a,
    panel_title = "New overall CH risk loci",
    driver_levels = driver_levels_a,
    shape_vals = shape_vals_a,
    driver_label_expr = driver_label_expr_a,
    dodge_x_axis = TRUE
)

p2 <- make_plot(
    prep = prep_b,
    panel_title = "New driver-specific CH risk loci",
    driver_levels = driver_levels_b,
    shape_vals = shape_vals_b,
    driver_label_expr = driver_label_expr_b,
    dodge_x_axis = TRUE
)

# ---- Save original individual figures as separate outputs ----
ggsave(
    "fig2a.png",
    p1,
    width = 18,
    height = 11.5,
    dpi = 600,
    bg = "white"
)

ggsave(
    "fig2b.png",
    p2,
    width = 18,
    height = 11.5,
    dpi = 600,
    bg = "white"
)

# ---- Add panel labels only for combined multi-panel figure ----
p1_labeled <- p1 +
    labs(tag = "a") +
    theme(
        plot.tag = element_text(size = 42, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = c(0.01, 0.99)
    )

p2_labeled <- p2 +
    labs(tag = "b") +
    theme(
        plot.tag = element_text(size = 42, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = c(0.01, 0.99)
    )

# ---- Save combined multi-panel figure ----
png(
    filename = "fig2.png",
    width = 18,
    height = 23,
    units = "in",
    res = 600,
    bg = "white"
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))

print(p1_labeled, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2_labeled, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

dev.off()

# ---- Print plots to current device ----
print(p1)
print(p2)

# ---- PDF outputs using base R Cairo ----

# fig2a
cairo_pdf("fig2a.pdf", width = 18, height = 11.5, bg = "white")
print(p1)
dev.off()

# fig2b
cairo_pdf("fig2b.pdf", width = 18, height = 11.5, bg = "white")
print(p2)
dev.off()

# combined fig2
cairo_pdf("fig2.pdf", width = 18, height = 23, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
print(p1_labeled, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2_labeled, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
