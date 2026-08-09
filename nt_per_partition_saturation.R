#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(tools)
})

has_gg <- requireNamespace("ggplot2", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript nt_saturation_by_partition.R <alignment.fasta> <outdir> [partition.txt]\n")
  quit(status = 1)
}

fa <- args[1]
outdir <- args[2]
partfile <- if (length(args) >= 3) args[3] else NA
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==================== HELPERS ====================

flatten_dist <- function(D){
  M <- as.matrix(D)
  ij <- which(upper.tri(M), arr.ind = TRUE)
  data.frame(
    taxon1 = rownames(M)[ij[,1]],
    taxon2 = colnames(M)[ij[,2]],
    dist   = M[ij],
    check.names = FALSE
  )
}

# ---------- universal plot writer (PNG + SVG) ----------
plot_dual <- function(plot_fun, outfile, width_px=1200, height_px=900, res=150){
  # PNG
  png(outfile, width=width_px, height=height_px, res=res)
  plot_fun()
  dev.off()

  # SVG
  svgfile <- sub("\\.png$", ".svg", outfile)
  svg(svgfile, width=width_px/100, height=height_px/100)
  plot_fun()
  dev.off()
}

# ---------- regression plots ----------
reg_plot <- function(df, x, y, title, outfile){
  if (!all(c(x,y) %in% names(df))) return(NULL)
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]), , drop=FALSE]
  if (!nrow(df)) return(NULL)

  fit <- lm(df[[y]] ~ df[[x]])
  slope <- unname(coef(fit)[2])
  r2 <- summary(fit)$r.squared

  plot_fun <- function(){
    plot(df[[x]], df[[y]], pch=16, cex=0.5,
         xlab=x, ylab=y,
         main=paste0(title,
                     "\nSlope=", round(slope,3),
                     "  R²=", round(r2,3)))
    abline(fit, lwd=2)
    grid()
  }

  plot_dual(plot_fun, outfile)

  list(slope=slope, r2=r2)
}

# ---------- weighted regression ----------
reg_plot_weighted <- function(df, x, y, w, title, outfile){
  need <- c(x,y,w)
  if (!all(need %in% names(df))) return(NULL)

  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]) &
           is.finite(df[[w]]) & df[[w]] > 0, , drop=FALSE]
  if (!nrow(df)) return(NULL)

  fit <- lm(df[[y]] ~ df[[x]], weights=df[[w]])
  slope <- unname(coef(fit)[2])
  r2 <- summary(fit)$r.squared

  plot_fun <- function(){
    plot(df[[x]], df[[y]], pch=16, cex=0.4,
         xlab=x, ylab=y,
         main=paste0(title,
                     "\nWeighted slope=", round(slope,3),
                     "  R²=", round(r2,3),
                     "\nN=", format(nrow(df), big.mark=",")))
    abline(fit, lwd=2)
    grid()
  }

  plot_dual(plot_fun, outfile)

  list(slope=slope, r2=r2, n=nrow(df))
}

# ==================== PARTITIONS ====================

read_partitions_df <- function(path, nchar_total){
  if (is.na(path)) stop("Partition file required.")

  L <- readLines(path, warn = FALSE)
  norm <- function(x) trimws(gsub("[;\r]", "", sub("#.*","", x)))

  rows <- list()

  for (ln in L){
    s <- norm(ln)
    if (s == "") next

    nm <- NA; st <- en <- NA

    if (grepl(",", s)){
      r <- regmatches(s, regexec("([^=]+)=\\s*([0-9]+)-([0-9]+)", s))[[1]]
      if (length(r)>=4){
        nm <- trimws(r[2]); st <- as.integer(r[3]); en <- as.integer(r[4])
      }
    } else {
      r <- strsplit(s, "\\s+")[[1]]
      if (length(r)>=2){
        st <- as.integer(r[length(r)-1])
        en <- as.integer(r[length(r)])
        nm <- paste(r[-c(length(r)-1, length(r))], collapse="_")
      }
    }

    if (!is.na(st) && !is.na(en) && en>=st){
      st <- max(1, st); en <- min(nchar_total, en)
      rows[[length(rows)+1]] <- data.frame(
        name = ifelse(is.na(nm) || nm=="", paste0("block_",length(rows)+1), nm),
        start = st, end = en,
        stringsAsFactors=FALSE
      )
    }
  }

  df <- do.call(rbind, rows)
  df$label <- gsub("[^A-Za-z0-9._-]+","_",df$name)
  df[order(df$start),]
}

mk_pos <- function(st,en,pos){
  idx <- seq(st+pos-1,en,3)
  idx[idx>=st & idx<=en]
}

# ==================== LOAD DATA ====================

dnab <- read.dna(fa, format="fasta")
nchar <- ncol(dnab)
message(sprintf("[Info] Taxa=%d Sites=%d", nrow(dnab), nchar))

parts <- read_partitions_df(partfile, nchar)
message(sprintf("[Info] Partitions=%d", nrow(parts)))

# ==================== MAIN ====================

combined <- list(all=NULL,pos1=NULL,pos2=NULL,pos3=NULL)
summary_rows <- list()

run_subset <- function(idx, label, weight){
  dsub <- dnab[, idx, drop=FALSE]

  p_raw <- dist.dna(dsub, model="raw", pairwise.deletion=TRUE)
  d_cor <- dist.dna(dsub, model="TN93", pairwise.deletion=TRUE)

  dfp <- transform(flatten_dist(p_raw), p=dist)[,c("taxon1","taxon2","p")]
  dfm <- transform(flatten_dist(d_cor), corr=dist)[,c("taxon1","taxon2","corr")]

  df <- merge(dfp, dfm)
  if (!nrow(df)) return(NULL)

  df$weight_subset_len <- weight

  outfile <- file.path(outdir, paste0("plot_",label,".png"))
  stats <- reg_plot(df, "p","corr",
                    paste0("NT saturation: ",label), outfile)

  list(df=df, stats=stats)
}

for (i in seq_len(nrow(parts))){
  st <- parts$start[i]; en <- parts$end[i]
  lbl <- parts$label[i]

  idx_all <- st:en
  idx1 <- mk_pos(st,en,1)
  idx2 <- mk_pos(st,en,2)
  idx3 <- mk_pos(st,en,3)

  res_all <- run_subset(idx_all, paste0(lbl,"_all"), length(idx_all))
  if (!is.null(res_all)){
    combined$all <- rbind(combined$all, res_all$df)
  }

  res1 <- run_subset(idx1, paste0(lbl,"_pos1"), length(idx1))
  res2 <- run_subset(idx2, paste0(lbl,"_pos2"), length(idx2))
  res3 <- run_subset(idx3, paste0(lbl,"_pos3"), length(idx3))

  if (!is.null(res1)) combined$pos1 <- rbind(combined$pos1, res1$df)
  if (!is.null(res2)) combined$pos2 <- rbind(combined$pos2, res2$df)
  if (!is.null(res3)) combined$pos3 <- rbind(combined$pos3, res3$df)
}

# ==================== WEIGHTED COMBINED ====================

for (nm in names(combined)){
  df <- combined[[nm]]
  if (is.null(df) || !nrow(df)) next

  outfile <- file.path(outdir, paste0("combined_",nm,".png"))
  reg_plot_weighted(df,"p","corr","weight_subset_len",
                    paste("Weighted:",nm), outfile)
}

message("[OK] Finished. PNG + SVG outputs generated.")