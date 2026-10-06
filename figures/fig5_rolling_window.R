# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

suppressMessages({
library(data.table)
library(ggplot2)
library(pals)
library(pwalign)
library(BSgenome.Hsapiens.UCSC.hg38)
})
genome <- BSgenome.Hsapiens.UCSC.hg38

res_dir <- paste0(PROJECT_DIR, "/results/paper_results/")
fig_dir <- paste0(PROJECT_DIR, "/figs/paper_figs/")

alu_consensus <- fread(paste0(PROJECT_DIR, "/results/20240208/20240129_alu_consensus_seq_data.txt"))

# standard text-size + bordered-strip theme used for this style of figure
# elsewhere in bin/ (bin/plot_alu5100_akita.R, bin/plot_alu_pairs_interaction.R,
# figures/fig4_gene_swaps.R, figures/fig5_random_300bp.R)
theme_textsize_border <- theme(
  axis.title = element_text(size = 7),
  strip.text = element_text(size = 7),
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 5),
  legend.key.size = unit(5, "mm"),
  axis.text = element_text(size = 5),
  axis.line  = element_line(linewidth = 0.3),
  axis.ticks = element_line(linewidth = 0.3),
  strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
)

mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)

# ── load & merge one dataset's rolling-window batches (input + scores) ──────
load_rolling <- function(dir_path, prefix, n_batches = 12) {
  dts <- lapply(0:(n_batches-1), function(i) {
    infile  <- sprintf("%s/input/%s_rolling_%d.txt", dir_path, prefix, i)
    scfile  <- sprintf("%s/output/%s_rolling_%d_scores", dir_path, prefix, i)
    if (!file.exists(infile) || !file.exists(scfile)) return(NULL)
    win <- fread(infile)
    sc  <- fread(scfile)
    stopifnot(nrow(win) == nrow(sc))
    cbind(win, sc[, -"var_index"])
  })
  dt <- rbindlist(dts, use.names = TRUE, fill = TRUE)
  mse_cols <- grep("^mse_", names(dt), value = TRUE)
  dt[, mse_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = mse_cols]
  dt
}

# ── align one element's genomic sequence to its subfamily consensus, ────────
# ── return the element(0-based)->consensus(0-based) position map ───────────
align_to_consensus <- function(chrom, start0, end0, cons_seq) {
  elem_seq <- tryCatch(as.character(getSeq(genome, chrom, start = start0 + 1, end = end0)),
                        error = function(e) NA_character_)
  if (is.na(elem_seq) || nchar(elem_seq) == 0) return(NULL)
  aln <- tryCatch(
    pairwiseAlignment(pattern = DNAString(cons_seq), subject = DNAString(elem_seq),
                       type = "global", substitutionMatrix = mat,
                       gapOpening = 5, gapExtension = 0.5),
    error = function(e) NULL)
  if (is.null(aln)) return(NULL)

  c_chars <- strsplit(as.character(pattern(aln)), "")[[1]]
  e_chars <- strsplit(as.character(subject(aln)), "")[[1]]
  e_is_gap <- e_chars == "-"; c_is_gap <- c_chars == "-"
  e_pos <- cumsum(!e_is_gap) - 1L
  c_pos <- cumsum(!c_is_gap) - 1L
  keep <- !e_is_gap & !c_is_gap
  map_dt <- data.table(e_pos = e_pos[keep], c_pos = c_pos[keep])
  setkey(map_dt, e_pos)
  map_dt
}

# ── full per-dataset pipeline: window scores -> consensus-position coverage ─
run_dataset <- function(window_dt, annot_dt, label) {
  annot_el <- unique(annot_dt[, .(orig_idx, CHROM, POS, END, repName, strand)])
  annot_el <- annot_el[repName != "Alu" & strand == "+" & repName %in% alu_consensus$Alu_name]

  wd <- merge(window_dt, annot_el, by = "orig_idx", suffixes = c("", "_elem"))
  cat(sprintf("[%s] %d window rows across %d elements, %d subfamilies\n",
              label, nrow(wd), uniqueN(wd$orig_idx), uniqueN(wd$repName)))

  elements <- unique(wd[, .(orig_idx, CHROM, POS_elem, END_elem, repName)])

  mapped_list <- vector("list", nrow(elements))
  t0 <- Sys.time()
  for (i in seq_len(nrow(elements))) {
    el <- elements[i]
    cons_seq <- alu_consensus[Alu_name == el$repName]$Alu_seq[1]
    map_dt <- align_to_consensus(el$CHROM, el$POS_elem, el$END_elem, cons_seq)
    if (is.null(map_dt)) next

    el_windows <- wd[orig_idx == el$orig_idx]
    win_start <- el_windows$POS - el$POS_elem
    win_end   <- el_windows$END - el$POS_elem
    c_start <- map_dt[.(win_start), c_pos, on = "e_pos"]
    c_end   <- map_dt[.(win_end - 1L), c_pos, on = "e_pos"]
    ok <- !is.na(c_start) & !is.na(c_end)
    if (!any(ok)) next
    mapped_list[[i]] <- data.table(subfamily = el$repName,
                                    consensus_start = c_start[ok],
                                    consensus_end   = c_end[ok],
                                    score = el_windows$mse_mean[ok])
    if (i %% 500 == 0) cat(sprintf("  [%s] %d / %d elements aligned (%.1fs elapsed)\n",
                                    label, i, nrow(elements), as.numeric(Sys.time()-t0, units="secs")))
  }
  mapped_df <- rbindlist(mapped_list, use.names = TRUE, fill = TRUE)
  cat(sprintf("[%s] mapped %d window-score rows total (%.1fs)\n", label, nrow(mapped_df),
              as.numeric(Sys.time()-t0, units="secs")))

  cov_list <- list()
  for (subfam in unique(mapped_df$subfamily)) {
    L <- nchar(alu_consensus[Alu_name == subfam]$Alu_seq[1])
    sub_rows <- mapped_df[subfamily == subfam]
    sum_v <- numeric(L); cnt_v <- numeric(L)
    for (r in seq_len(nrow(sub_rows))) {
      s <- sub_rows$consensus_start[r] + 1L
      e <- min(sub_rows$consensus_end[r] + 1L, L)
      if (s > e || s < 1) next
      sum_v[s:e] <- sum_v[s:e] + sub_rows$score[r]
      cnt_v[s:e] <- cnt_v[s:e] + 1
    }
    avg <- ifelse(cnt_v > 0, sum_v / cnt_v, NA_real_)
    cov_list[[subfam]] <- data.table(position = seq_len(L) - 1L, score = avg, subfamily = subfam)
  }
  coverage <- rbindlist(cov_list, use.names = TRUE, fill = TRUE)

  # which elements actually contributed >=1 successfully-aligned window
  # (used to count how many plotted elements are globally high-disruption)
  used_idx <- elements[sapply(mapped_list, function(x) !is.null(x))]
  list(coverage = coverage, elements_used = used_idx[, .(orig_idx, repName)])
}

# subfamily block order: AluY, AluS, AluJ, then FAM/FLAM/FRAM immediately
# after AluJ (FAM/FLAM/FRAM are the ancestral monomeric precursors AluJ
# derives from, so grouping them adjacent to AluJ reflects that lineage;
# within each block, plain alphabetical order, matching the convention in
# bin/plot_alu5100_akita.R).
order_repName_levels <- function(repName_levels) {
  repName   <- sort(unique(repName_levels))
  fam_names <- sort(repName[repName %in% c('FAM', 'FLAM_A', 'FLAM_C', 'FRAM')])
  j_names   <- sort(repName[substr(repName, 1, 4) == 'AluJ'])
  s_names   <- sort(repName[substr(repName, 1, 4) == 'AluS'])
  y_names   <- sort(repName[substr(repName, 1, 4) == 'AluY'])
  generic   <- sort(repName[repName == 'Alu'])
  stopifnot(length(generic) + length(j_names) + length(s_names) + length(y_names) + length(fam_names) == length(repName))
  c(generic, y_names, s_names, j_names, fam_names)
}

# per-subfamily color convention matching the established figure style
# (SuppFigure_1a.ai panel d, bin/plot_alu5100_akita.R): AluJ = yellow-orange,
# AluS = greens, AluY = blues, FAM/FLAM/FRAM = pink-magenta, light-to-dark
# within each block (color assignment is independent of facet display order).
build_repName_colors <- function(repName_levels) {
  repName <- sort(repName_levels)
  n_fam     <- sum(repName %in% c('FAM', 'FLAM_A', 'FLAM_C', 'FRAM'))
  n_j       <- sum(substr(repName, 1, 4) == 'AluJ')
  n_s       <- sum(substr(repName, 1, 4) == 'AluS')
  n_y       <- sum(substr(repName, 1, 4) == 'AluY')
  n_generic <- sum(repName == 'Alu')
  stopifnot(n_generic + n_j + n_s + n_y + n_fam == length(repName))

  stepped_colors <- stepped(24)
  cols <- c(
    if (n_generic) stepped_colors[22] else NULL,
    brewer.ylorbr(6)[2:5][seq_len(n_j)],
    brewer.greens(20)[4:20][seq_len(n_s)],
    brewer.pubu(32)[8:32][seq_len(n_y)],
    brewer.rdpu(8)[5:8][seq_len(n_fam)]
  )
  setNames(cols, repName)
}

# facet strip label shorthand: "AluYb8 (18/46)" = 18 of the 46 plotted
# elements for that subfamily are globally high-disruption (mse_mean >=
# genome-wide 99th percentile, same threshold used elsewhere in this project,
# e.g. bin/chromHMM_state10_celltype_cooccurrence.R).
build_disruption_labels <- function(elements_used, annot_dt) {
  p99 <- quantile(annot_dt$mse_mean, 0.99, na.rm = TRUE)
  el <- merge(elements_used, annot_dt[, .(orig_idx, mse_mean)], by = "orig_idx")
  el[, high := mse_mean >= p99]
  el[, .(n_total = .N, n_high = sum(high, na.rm = TRUE)), by = repName][
    , label := sprintf("%s (%d/%d)", repName, n_high, n_total)][]
}

make_plot <- function(plot_df, disruption_labels, out_file, title,
                       ncol_facets = 6, width_mm = 180, max_height_mm = 215) {
  level_order <- order_repName_levels(unique(plot_df$subfamily))
  plot_df$subfamily <- factor(plot_df$subfamily, levels = level_order)

  label_lookup <- setNames(disruption_labels$label, disruption_labels$repName)
  plot_df[, facet_label := label_lookup[as.character(subfamily)]]
  facet_levels <- label_lookup[level_order]
  plot_df$facet_label <- factor(plot_df$facet_label, levels = facet_levels)

  n_subfam    <- uniqueN(plot_df$subfamily)
  nrow_facets <- ceiling(n_subfam / ncol_facets)
  panel_mm    <- width_mm / ncol_facets
  # square-ish panels (30mm x ~27mm here) at 6 cols would run to 240mm --
  # over the 215mm cap -- so clamp height rather than add a 7th column,
  # since panel width (x-axis is consensus position, 0-300bp) matters more
  # for legibility here than panel height
  height_mm <- min(panel_mm * nrow_facets, max_height_mm)

  subfam_colors <- build_repName_colors(level_order)
  facet_colors  <- setNames(subfam_colors[level_order], facet_levels)

  p <- ggplot(plot_df, aes(x = position, y = score, color = facet_label)) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(~facet_label, ncol = ncol_facets, scales = "free_x") +
    scale_color_manual(values = facet_colors, guide = "none") +
    labs(x = "Consensus position", y = "Score (mean MSE)", title = title,
         subtitle = "facet labels: (n high-disruption / n total); high = mse_mean ≥ genome-wide 99th percentile") +
    theme_classic() + theme_textsize_border
  ggsave(out_file, p, width = width_mm / 25.4, height = height_mm / 25.4,
         units = "in", limitsize = FALSE)
  cat("wrote", out_file, sprintf(" (%.0fmm x %.0fmm, %d cols x %d rows)\n",
      width_mm, height_mm, ncol_facets, nrow_facets))
}

# ═══════════════════════════════ top100 dataset ═════════════════════════════
# NOTE (2026-08-16, superseded 2026-08-31): the per-batch loader (load_rolling)
# silently dropped batch 11 (10,000 window rows) because its _scores file had
# been split into 4 unexpectedly-named sub-chunks covering only ~62% of that
# batch, with the remaining ~38% never scored at all. That's now fixed on the
# cluster -- results/alu_top100_repName_window_scores_combined.txt is the
# complete, re-run, single-file combined output (116,432 rows, all 5,094
# elements, 0 NA in mse_mean, strand included), so we read it directly instead
# of stitching per-batch input/output files.
#
# The "5100" directory downloaded alongside the original top100 job is
# byte-identical to top100 (same 5094 orig_idx, exactly 100 elements/
# subfamily) -- confirmed via all.equal(). The real 5100-broader-sample
# rolling-window SLURM job (alu5100/alu5100_40640.txt) was cancelled for
# hitting its time limit and never produced output, so only top100 is real
# data right now. When the true 5100 run completes, rerun this same
# run_dataset()/make_plot() pair on it.
window_top100 <- fread(paste0(PROJECT_DIR, "/results/alu_top100_repName_window_scores_combined.txt"))
annot_master <- fread(paste0(res_dir, '20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'),
                       select = c("orig_idx","CHROM","POS","END","repName","strand","mse_mean"))
result_top100 <- run_dataset(window_top100, annot_master, "top100")
fwrite(result_top100$coverage, paste0(fig_dir, '20260817_rolling_window_top100_coverage.csv'))
fwrite(result_top100$elements_used, paste0(fig_dir, '20260817_rolling_window_top100_elements_used.csv'))
labels_top100 <- build_disruption_labels(result_top100$elements_used, annot_master)
make_plot(result_top100$coverage, labels_top100,
          paste0(fig_dir, '20260817_rolling_window_top100.pdf'), "Rolling window (top100 per subfamily)")

cat("DONE\n")
