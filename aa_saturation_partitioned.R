#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ape); library(phangorn) })

# ---------- args ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript aa_saturation_partitioned_v2.R <alignment.fasta> <outdir> <partition.txt> [MODEL|AUTO]\n",
      "  MODEL: one of WAG,JTT,LG,Dayhoff,cpREV,mtmam,mtArt,MtZoa,mtREV24,VT,RtREV,HIVw,HIVb,FLU,Blosum62,Dayhoff_DCMut,JTT_DCMut\n",
      "  AUTO : use the model listed per partition line (before the comma)\n")
  quit(status = 1)
}
fa        <- args[1]
outdir    <- args[2]
partfile  <- args[3]
model_arg <- if (length(args) >= 4) args[4] else "JTT"

ok_models <- c("WAG","JTT","LG","Dayhoff","cpREV","mtmam","mtArt","MtZoa","mtREV24",
               "VT","RtREV","HIVw","HIVb","FLU","Blosum62","Dayhoff_DCMut","JTT_DCMut")
use_auto <- toupper(model_arg) == "AUTO"
if (!use_auto && !(model_arg %in% ok_models)) stop("Unsupported model: ", model_arg)

canon_model <- function(x){
  m <- tolower(trimws(x))
  switch(m,
    "wag"="WAG","jtt"="JTT","lg"="LG","dayhoff"="Dayhoff","cprev"="cpREV",
    "mtmam"="mtmam","mtart"="mtArt","mtzoa"="MtZoa","mtrev24"="mtREV24",
    "vt"="VT","rtrev"="RtREV","hivw"="HIVw","hivb"="HIVb","flu"="FLU",
    "blosum62"="Blosum62","dayhoff_dcmut"="Dayhoff_DCMut","jtt_dcmut"="JTT_DCMut",
    NA_character_
  )
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
flatten_dist <- function(D){
  M <- as.matrix(D); ij <- which(upper.tri(M), arr.ind = TRUE)
  data.frame(taxon1 = rownames(M)[ij[,1]],
             taxon2 = colnames(M)[ij[,2]],
             dist   = M[ij],
             stringsAsFactors = FALSE, check.names = FALSE)
}

reg_and_png <- function(df, x, y, title, outfile_png){
  ok <- is.finite(df[[x]]) & is.finite(df[[y]])
  df <- df[ok, , drop = FALSE]
  fit <- lm(df[[y]] ~ df[[x]])
  slope <- unname(coef(fit)[2]); r2 <- summary(fit)$r.squared

  # PNG output
  png(outfile_png, width = 1200, height = 900, res = 150)
  plot(df[[x]], df[[y]], pch = 16, cex = 0.4, xlab = x, ylab = y,
       main = paste0(title, "\nSlope=", round(slope,3), "  R²=", round(r2,3)))
  abline(fit, lwd = 2); grid()
  dev.off()

  # SVG output (same file path, just .svg extension)
  outfile_svg <- sub("\\.png$", ".svg", outfile_png)
  svg(outfile_svg, width = 12, height = 9)
  plot(df[[x]], df[[y]], pch = 16, cex = 0.4, xlab = x, ylab = y,
       main = paste0(title, "\nSlope=", round(slope,3), "  R²=", round(r2,3)))
  abline(fit, lwd = 2); grid()
  dev.off()

  c(slope = slope, r2 = r2)
}

# Parse "MODEL, name = start-end" or "name = start-end". Clamp in SITE space (1..total_sites).
parse_partitions <- function(path, total_sites, report_path = NULL){
  L  <- readLines(path, warn = FALSE)
  rows <- list(); bad <- list()
  norm <- function(x){ trimws(sub("#.*$", "", gsub("\r", "", x, fixed=TRUE))) }
  re1 <- "^\\s*([A-Za-z0-9_.+\\-]+)\\s*,\\s*([^=]+?)\\s*=\\s*([0-9]+)\\s*-\\s*([0-9]+)\\s*$"
  re2 <- "^\\s*([^=]+?)\\s*=\\s*([0-9]+)\\s*-\\s*([0-9]+)\\s*$"

  for (i in seq_along(L)){
    raw <- L[i]; s <- norm(raw)
    if (s == "") { bad[[length(bad)+1L]] <- data.frame(line=i,reason="blank/comment",text=raw); next }
    mdl <- NA_character_; nm <- NA_character_; st <- NA_integer_; en <- NA_integer_

    if (grepl(re1, s, perl=TRUE)){
      r <- regmatches(s, regexec(re1, s, perl=TRUE))[[1]]
      mdl <- trimws(r[2]); nm <- trimws(r[3]); st <- as.integer(r[4]); en <- as.integer(r[5])
    } else if (grepl(re2, s, perl=TRUE)){
      r <- regmatches(s, regexec(re2, s, perl=TRUE))[[1]]
      nm <- trimws(r[2]); st <- as.integer(r[3]); en <- as.integer(r[4])
    } else {
      bad[[length(bad)+1L]] <- data.frame(line=i,reason="regex_mismatch",text=raw); next
    }
    if (!is.finite(st) || !is.finite(en)){
      bad[[length(bad)+1L]] <- data.frame(line=i,reason="non_integer_bounds",text=raw); next
    }
    # clamp in site space
    st0 <- st; en0 <- en
    st <- max(1L, min(st, total_sites)); en <- max(1L, min(en, total_sites))
    if (en < st){
      bad[[length(bad)+1L]] <- data.frame(
        line=i, reason="empty_after_clamp_site_space",
        text=sprintf("%s  (orig %d-%d; total_sites=%d)", raw, st0, en0, total_sites))
      next
    }
    rows[[length(rows)+1L]] <- data.frame(model=mdl, name=nm, start_site=st, end_site=en,
                                          stringsAsFactors=FALSE)
  }

  df <- if (length(rows)) do.call(rbind, rows) else stop("No valid partitions parsed")
  df <- df[order(df$start_site, df$end_site, df$name), , drop=FALSE]

  if (length(bad) && !is.null(report_path)) {
    write.table(do.call(rbind, bad), report_path, sep="\t", row.names=FALSE, quote=FALSE)
  }
  df
}

# Map site index (1..sum(weights)) → pattern index (1..nr)
make_site_to_pat <- function(weights){
  cumw <- cumsum(weights)
  function(s) {
    # returns the first pattern whose cumulative weight >= s
    idx <- which(cumw >= s)
    if (length(idx)) idx[1] else NA_integer_
  }
}

# ---------- read alignment ----------
aa_pd <- phangorn::read.phyDat(fa, format = "fasta", type = "AA")
pat_total   <- attr(aa_pd, "nr")          # patterns
w           <- attr(aa_pd, "weight")      # pattern weights
total_sites <- sum(w)                     # true columns
taxa        <- attr(aa_pd, "names")

cat(sprintf("[Info] Taxa=%d  Patterns=%d  Sites=%d\n", length(taxa), pat_total, total_sites))

# parse partitions in SITE space
parts_site <- parse_partitions(partfile, total_sites,
                               report_path = file.path(outdir, "partition_parse_report.tsv"))
message("[parse] Parsed ", nrow(parts_site), " partitions (site space).")

# site->pattern mapper
site_to_pat <- make_site_to_pat(w)

# add pattern-mapped coordinates
parts_site$pat_start <- vapply(parts_site$start_site, site_to_pat, integer(1))
parts_site$pat_end   <- vapply(parts_site$end_site,   site_to_pat, integer(1))

# audit any that still ended up NA (shouldn't happen unless sites > total_sites etc.)
mask_ok <- is.finite(parts_site$pat_start) & is.finite(parts_site$pat_end) &
           parts_site$pat_end >= parts_site$pat_start

if (!all(mask_ok)) {
  bad2 <- parts_site[!mask_ok, ]
  if (nrow(bad2)){
    write.table(bad2, file.path(outdir, "partition_map_to_patterns_failed.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE)
    message("[warn] ", nrow(bad2), " partitions could not be mapped to patterns (see partition_map_to_patterns_failed.tsv).")
  }
}
parts <- parts_site[mask_ok, , drop=FALSE]
message("[map] Usable partitions after site→pattern mapping: ", nrow(parts), "/", nrow(parts_site))

# Build a stable pair key order from the first *usable* partition
if (!nrow(parts)) stop("No usable partitions after mapping to pattern indices.")
first_idx <- 1L
tmp_subset <- subset(aa_pd, select = parts$pat_start[first_idx]:parts$pat_end[first_idx])
pair_template <- within(flatten_dist(phangorn::dist.hamming(tmp_subset, ratio=TRUE)),
                        { key <- paste(taxon1, taxon2, sep="||") })[, c("key","taxon1","taxon2")]

all_pairs <- list(); slopes <- list()

# ---------- main loop ----------
for (i in seq_len(nrow(parts))) {
  pr  <- parts[i, ]
  stp <- pr$pat_start; enp <- pr$pat_end
  mdl_raw <- pr$model
  mdl <- if (use_auto && !is.na(mdl_raw)) canon_model(mdl_raw) else model_arg
  if (is.na(mdl) || !(mdl %in% ok_models)) {
    warning(sprintf("Partition %d ('%s'): unknown model '%s'; using global model '%s'.",
                    i, pr$name, as.character(mdl_raw), model_arg))
    mdl <- model_arg
  }

  seg <- subset(aa_pd, select = stp:enp)
  if (!identical(attr(seg, "names"), taxa)) stop("Taxon label mismatch in partition ", i)

  # distances (AA p and model-corrected)
  p_raw <- phangorn::dist.hamming(seg, ratio = TRUE)
  m_dst <- phangorn::dist.ml(seg, model = mdl)

  dfp <- transform(flatten_dist(p_raw), p = dist)[, c("taxon1","taxon2","p")]
  dfm <- transform(flatten_dist(m_dst),  md = dist)[, c("taxon1","taxon2","md")]
  df  <- merge(dfp, dfm, by = c("taxon1","taxon2"), sort = FALSE)
  df$key <- paste(df$taxon1, df$taxon2, sep="||")
  df  <- merge(pair_template, df, by = c("key","taxon1","taxon2"), sort = FALSE)

  safe_name <- gsub("[^A-Za-z0-9._-]+", "_", pr$name)
  tsv  <- file.path(outdir, sprintf("pairs_%03d_%s.tsv", i, safe_name))
  pngf <- file.path(outdir, sprintf("plot_%03d_%s.png",  i, safe_name))
  write.table(df[, c("taxon1","taxon2","p","md")], tsv, sep="\t", row.names=FALSE, quote=FALSE)

  title_txt <- sprintf("AA saturation: %s (%s, sites %d-%d → patterns %d-%d)",
                       pr$name, mdl, pr$start_site, pr$end_site, stp, enp)
  stt <- reg_and_png(df, "p", "md", title_txt, pngf)

  slopes[[length(slopes)+1L]] <- data.frame(
    partition      = i,
    name           = pr$name,
    model          = mdl,
    start_site     = pr$start_site,
    end_site       = pr$end_site,
    pat_start      = stp,
    pat_end        = enp,
    slope          = unname(stt["slope"]),
    r2             = unname(stt["r2"]),
    stringsAsFactors = FALSE
  )
  all_pairs[[length(all_pairs)+1L]] <- df[, c("taxon1","taxon2","p","md")]
  cat(sprintf("[Done] %3d: %-35s model=%-8s site_len=%-5d pat_len=%-5d slope=%.3f  R²=%.3f\n",
              i, pr$name, mdl,
              pr$end_site - pr$start_site + 1L,
              enp - stp + 1L,
              stt["slope"], stt["r2"]))
}

# ---------- summaries ----------
slopes_df <- do.call(rbind, slopes)
write.table(slopes_df, file.path(outdir, "partition_slopes.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# pooled regression across all partitions (simple unweighted)
all_df <- do.call(rbind, all_pairs)
agg_png <- file.path(outdir, "aggregate_plot.png")
agg_st  <- reg_and_png(all_df, "p", "md", "AA saturation (all partitions pooled)", agg_png)
write.table(data.frame(model=if (use_auto) "AUTO" else model_arg,
                       slope=unname(agg_st["slope"]), r2=unname(agg_st["r2"])),
            file.path(outdir, "aggregate_slope.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("[OK] Outputs in: ", outdir, "\n", sep = "")
