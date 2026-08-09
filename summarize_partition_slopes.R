#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(tools) })

# Usage:
# Rscript summarize_partition_slopes.R <pairs_dir> <partition_slopes.tsv> <outdir> [r2_min] [len_min]
# - pairs_dir: directory containing files like pairs_###_NAME.tsv (from aa_saturation_partitioned.R)
# - partition_slopes.tsv: the summary with start/end (to get per-partition lengths)
# - r2_min / len_min: optional filters on partitions before aggregating

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript summarize_partition_slopes.R <pairs_dir> <partition_slopes.tsv> <outdir> [r2_min] [len_min]\n")
  quit(status = 1)
}
pairs_dir <- args[1]
slopes_tsv <- args[2]
outdir     <- args[3]
r2_min     <- if (length(args) >= 4) as.numeric(args[4]) else 0.0
len_min    <- if (length(args) >= 5) as.integer(args[5]) else 1L

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- read partition lengths / filter ---
S <- read.table(slopes_tsv, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")

# Back-compat: accept either start/end or start_site/end_site
pick_col <- function(nms, candidates) {
  lnms <- tolower(nms)
  for (c in candidates) {
    hit <- which(lnms == tolower(c))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

cn <- names(S)
col_start <- pick_col(cn, c("start", "start_site", "site_start"))
col_end   <- pick_col(cn, c("end",   "end_site",   "site_end"))

# Keep original stopifnot logic but normalize to start/end first
if (is.na(col_start) || is.na(col_end)) {
  stop("Could not find start/end columns (tried: start,end OR start_site,end_site).")
}

# Create temporary normalized columns so rest of the script stays unchanged
S$start <- as.integer(S[[col_start]])
S$end   <- as.integer(S[[col_end]])

# Ensure required columns exist
stopifnot(all(c("partition","start","end","slope","r2") %in% names(S)))

S$len  <- with(S, end - start + 1L)
S$keep <- with(S, r2 >= r2_min & len >= len_min)

# map partition index -> weight and keep flag
w_map  <- setNames(S$len,  S$partition)
ok_map <- setNames(S$keep, S$partition)

# --- read all per-partition pair tables ---
fs <- list.files(pairs_dir, pattern = "^pairs_\\d+_.*\\.tsv$", full.names = TRUE)
if (!length(fs)) stop("No pairs_*.tsv files found in: ", pairs_dir)

read_one <- function(f){
  # extract integer index from "pairs_###_....tsv"
  m <- regexec("pairs_(\\d+)_", basename(f)); r <- regmatches(basename(f), m)[[1]]
  if (length(r) < 2) return(NULL)
  i <- as.integer(r[2])
  if (is.na(i) || !(i %in% names(ok_map)) || !ok_map[[as.character(i)]]) return(NULL)
  d <- try(read.table(f, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE), silent = TRUE)
  if (inherits(d, "try-error")) return(NULL)
  # expect: taxon1 taxon2 p md
  if (!all(c("taxon1","taxon2","p","md") %in% names(d))) return(NULL)
  d$partition <- i
  d$w <- w_map[[as.character(i)]]
  d$key <- paste(d$taxon1, d$taxon2, sep="||")
  d
}

DL <- lapply(fs, read_one)
DL <- DL[ lengths(DL) > 0 ]
if (!length(DL)) stop("No usable pair tables after filtering.")

D <- do.call(rbind, DL)

# --- length-weighted collapse to one point per taxon pair ---
# For each pair key, compute weighted mean p and md (weights = partition lengths)
wmean <- function(x, w) sum(w * x) / sum(w)

agg <- do.call(rbind, lapply(split(D, D$key), function(df){
  data.frame(
    taxon1 = df$taxon1[1],
    taxon2 = df$taxon2[1],
    p_w    = wmean(df$p,  df$w),
    md_w   = wmean(df$md, df$w)
  )
}))

# --- regression + plot ---
fit <- lm(md_w ~ p_w, data = agg)
slope <- unname(coef(fit)[2]); r2 <- summary(fit)$r.squared

# PNG
png_out <- file.path(outdir, sprintf("aggregate_length_weighted_r2>=%g_len>=%d.png", r2_min, len_min))
png(png_out, width = 1200, height = 900, res = 150)
plot(agg$p_w, agg$md_w, pch = 16, cex = 0.4,
     xlab = "p (length-weighted across partitions)",
     ylab = "model distance (length-weighted across partitions)",
     main = sprintf("AA saturation (length-weighted aggregate)\nSlope=%.3f  R²=%.3f  (parts kept: %d/%d)",
                    slope, r2, sum(S$keep), nrow(S)))
abline(fit, lwd = 2)
grid()
dev.off()

# SVG
svg_out <- sub("\\.png$", ".svg", png_out)
svg(svg_out, width = 12, height = 9)
plot(agg$p_w, agg$md_w, pch = 16, cex = 0.4,
     xlab = "p (length-weighted across partitions)",
     ylab = "model distance (length-weighted across partitions)",
     main = sprintf("AA saturation (length-weighted aggregate)\nSlope=%.3f  R²=%.3f  (parts kept: %d/%d)",
                    slope, r2, sum(S$keep), nrow(S)))
abline(fit, lwd = 2)
grid()
dev.off()

# also write the numbers
write.table(data.frame(
  partitions_total = nrow(S),
  partitions_kept  = sum(S$keep),
  r2_min = r2_min, len_min = len_min,
  slope = slope, r2 = r2
), file.path(outdir, "aggregate_length_weighted_stats.tsv"),
sep = "\t", row.names = FALSE, quote = FALSE)

cat("[OK] Wrote:\n - ", png_out, "\n - ",
    file.path(outdir, "aggregate_length_weighted_stats.tsv"), "\n", sep = "")
