# =========================
# 1. Packages
# =========================
required_pkgs <- c(
  "dplyr", "readr", "stringr",
  "grid", "gridExtra", "gtable"
)

to_install <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(readr)
library(stringr)
library(grid)
library(gridExtra)
library(gtable)

# =========================
# 2. Global settings
# =========================
base_family <- "sans"   
body_fontsize <- 9
header_fontsize <- 9

# =========================
# 3. Read data
# =========================
dat_b <- read_tsv("ed_table1.txt", show_col_types = FALSE)

# =========================
# 4. Helper
# =========================
clean <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x[x == "NA"] <- ""
  x
}

# =========================
# 5. Prepare display table
# =========================
dat_b_disp <- dat_b %>%
  mutate(
    rsid = clean(rsid),
    cs_size = clean(cs_size),
    highprob = clean(highprob),
    ot = clean(ot),
    cs2g = clean(cs2g),
    crispr = clean(crispr),
    scrnaseq = clean(scrnaseq)
  ) %>%
  transmute(
    `Locus` = rsid,
    `CS size` = cs_size,
    `PICS >50%` = highprob,
    `OpenTargets` = ot,
    `cs2G` = cs2g,
    `CRISPR` = crispr,
    `scRNA-Seq` = scrnaseq
  )

wrap_widths <- c(
  "Locus" = 14,
  "PICS >50%" = 14,
  "OpenTargets" = 14,
  "cs2G" = 14,
  "CRISPR" = 14,
  "scRNA-Seq" = 14
)

for (nm in names(wrap_widths)) {
  dat_b_disp[[nm]] <- str_wrap(dat_b_disp[[nm]], width = wrap_widths[[nm]])
}

# =========================
# 6. Build table
# =========================
tbl_theme <- ttheme_minimal(
  base_size = body_fontsize,
  base_family = base_family,
  core = list(
    fg_params = list(
      hjust = 0,
      x = 0.02,
      fontsize = body_fontsize,
      col = "black"
    ),
    bg_params = list(
      fill = rep(c("white", "#EEECE2"), length.out = nrow(dat_b_disp)),
      col = NA
    ),
    padding = unit(c(1.8, 2.4), "mm")
  ),
  colhead = list(
    fg_params = list(
      fontface = "bold",
      hjust = 0,
      x = 0.02,
      fontsize = header_fontsize,
      col = "black"
    ),
    bg_params = list(
      fill = "#D8D5C5",
      col = NA
    ),
    padding = unit(c(2.0, 2.4), "mm")
  )
)

tbl <- tableGrob(
  dat_b_disp,
  rows = NULL,
  theme = tbl_theme
)

tbl$widths <- unit(
  c(1.20, 0.70, 1.00, 1.15, 1.00, 0.95, 1.00),
  "in"
)

tbl$heights <- tbl$heights * 1.00

# =========================
# 7. Italicise gene columns
# =========================

gene_cols <- c(4, 5, 6, 7)

find_cell <- function(tbl, r, c) {
  which(tbl$layout$t == r &
          tbl$layout$l == c &
          tbl$layout$name == "core-fg")
}

n_rows <- nrow(dat_b_disp)

for (r in seq_len(n_rows)) {
  for (c in gene_cols) {
    idx <- find_cell(tbl, r + 1, c)
    if (length(idx) == 1) {
      txt <- tbl$grobs[[idx]]
      if (!is.null(txt$label) && nzchar(txt$label)) {
        tbl$grobs[[idx]]$gp <- gpar(
          fontface = "italic",
          fontsize = body_fontsize,
          fontfamily = base_family,
          col = "black"
        )
      }
    }
  }
}

# =========================
# 8. Add rules
# =========================

# top rule
tbl <- gtable_add_grob(
  tbl,
  segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1,
               gp = gpar(lwd = 1.0, col = "#3A3A3A")),
  t = 1, b = 1, l = 1, r = ncol(tbl)
)

# row rules
for (i in 1:(n_rows + 1)) {
  tbl <- gtable_add_grob(
    tbl,
    segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                 gp = gpar(lwd = 0.7, col = "#3A3A3A")),
    t = i, b = i, l = 1, r = ncol(tbl)
  )
}

# =========================
# 9. Draw and save
# =========================

table_width  <- convertWidth(sum(tbl$widths), "in", valueOnly = TRUE)
table_height <- convertHeight(sum(tbl$heights), "in", valueOnly = TRUE)

pad <- 0.05

# ---- PDF ----
pdf(
  "ed_table1.pdf",
  width = table_width + pad,
  height = table_height + pad,
  family = base_family
)
grid.newpage()
grid.draw(tbl)
dev.off()

# ---- PNG ----
png(
  "ed_table1.png",
  width = table_width + pad,
  height = table_height + pad,
  units = "in",
  res = 600
)
grid.newpage()
grid.draw(tbl)
dev.off()

# ---- Display in R ----
grid.newpage()
grid.draw(tbl)