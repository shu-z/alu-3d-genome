# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")
DOWNLOADS_DIR <- "/Users/shu/Downloads"

library(data.table)
library(ggplot2)
library(pals)
library(dplyr)
library(cowplot)


stepped_colors<-stepped(24)
stepped2_colors<-stepped2(24)
stepped3_colors<-stepped3(24)

fig_dir<-paste0(PROJECT_DIR, "/figs/paper_figs/")
res_dir<-paste0(PROJECT_DIR, "/results/paper_results/")
swap_res<-paste0(PROJECT_DIR, "/results/paper_results/20260511_gene_swaps/")

test1<-fread(paste0(swap_res, list.files(swap_res)[1]))

#theme things
theme_textsize_border<-theme(
  axis.title = element_text(size = 7),     # axis titles ~7 pt
  strip.text = element_text(size = 7),      # facet headers ~7 pt
  plot.title = element_text(size = 8),      # if you add a title
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 5),
  legend.key.size = unit(5, "mm"),
  axis.text = element_text(size = 5),
  axis.line  = element_line(linewidth = 0.3),  # x & y axis lines
  axis.ticks = element_line(linewidth = 0.3),
  strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
  
)
#for linewidth of all plots 
lwd<-0.3


alu_annot_short<-fread(paste0(res_dir, '20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'))
#jank way to remove duplicates
dup_cols<-!duplicated(names(alu_annot_short))
alu_annot_short <- alu_annot_short[,..dup_cols]
alu_annot_short$SVTYPE <-'DEL'

#get MSE percentiles 
p50 <- quantile(alu_annot_short$mse_mean, 0.5)
p99 <- quantile(alu_annot_short$mse_mean, 0.99)
# Classify alu_annot_short by its own mse_mean (same logic that produced alu_category, presumably)
alu_annot_short[, disruption := fcase(mse_mean >= p99, "high", mse_mean <= p50, "low", default = NA_character_)]



genes_set1 <- c("TAF7", "PCDHA1", "SART1", "PCDHGA11", "ATP1B2", "PCDHGA2",
                "HIST1H2BE", "SPG7", "TRIM38", "CPNE7", "TP53", "WRAP53",
                "MAP1B", "ZNF473", "MRPS27", "ITGAL", "NAA38", "AC104532.1",
                "ZNF771", "CAPN1", "PCDHGA1", "DNAH2", "NLRP8", "CFAP74", "NGEF")

genes_set2 <- c("PITPNC1", "ANKRD11", "ABR", "SPIDR", "DNAH2", "SPATS2", "SPG7",
                "ITGAL", "FRMD5", "GNB1", "YWHAE", "MECP2", "NGEF", "TAF7",
                "CRK", "FAM193A", "TUBA1C", "PTPRS", "SRPK2", "CFAP74", "AHNAK",
                "DIAPH1", "TP53", "RASSF3", "TMEM94")

gene_sets <- rbind(
  #data.table(gene = genes_set1, gene_set = "set1"),
  data.table(gene = genes_set2, gene_set = "set2")
)


metric_cols <- c("insert_DEL_mse", "orig_DEL_mse")


process_gene <- function(filepath, gene_name) {
  dt <- fread(filepath)
  
  # Average RC pairs
  rc_cols <- grep("_RC$", names(dt), value = TRUE)
  base_cols <- sub("_RC$", "", rc_cols)
  for (b in base_cols) {
    rc <- paste0(b, "_RC")
    dt[, (paste0(b, "_mean")) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(b, rc)]
  }
  
  # Melt to long
  keep_metrics <- intersect(metric_cols, names(dt))
  if (length(keep_metrics) == 0) return(NULL)
  
  dt_long<-dt
  dt_long <- melt(dt,measure.vars = keep_metrics, variable.name = "metric", value.name = "mse_mean_metric")
  dt_long[, gene := gene_name]
}

######### iterate over all genes 

gene_files <- list.files(swap_res, pattern = "_geneswap\\.txt$", full.names = TRUE)
file_info <- data.table(path = gene_files, gene = sub("_geneswap\\.txt$", "", basename(gene_files)))
file_info <- file_info[gene %in% unique(gene_sets$gene)]

dt_all <- rbindlist(lapply(seq_len(nrow(file_info)),
         function(i) process_gene(file_info$path[i], file_info$gene[i])), use.names = TRUE, fill = TRUE)


#remove orig or inserts <275 bp long 
dt_all <- merge(dt_all, alu_annot_short[, .(orig_idx, SVLEN)], all.x = TRUE)
setnames(dt_all, "SVLEN", "SVLEN_orig")
dt_all <- merge(dt_all, alu_annot_short[, .(orig_idx, SVLEN)],by.x = "alu_insert_orig_idx", by.y = "orig_idx",all.x = TRUE)
setnames(dt_all, "SVLEN", "SVLEN_insert")

dt_all_noshort<-dt_all[(dt_all$SVLEN_orig>=275) & (dt_all$SVLEN_insert>=275) ]

#########  Attach gene set and plot 

dt_all_geneset <- merge(dt_all_noshort, gene_sets, by = "gene", allow.cartesian = TRUE)
dt_all_geneset$metric<-factor(dt_all_geneset$metric, levels=c('orig_DEL_mse', 'insert_DEL_mse'))


# Build a fill key that encodes which color each row should get
dt_all_geneset[, fill_group := fcase(
  metric == "orig_DEL_mse"   & alu_orig_category   == "high", "#E9C4CA",
  metric == "orig_DEL_mse"   & alu_orig_category   == "low",  "#E1ECF4",
  metric == "insert_DEL_mse" & alu_insert_category == "high", "#A43D4E",
  metric == "insert_DEL_mse" & alu_insert_category == "low",  "#326386",
  default = "grey70"  # catches any other metric/category combos
)]

col_vals<-c('#f6f6f6', '#637939')

#calculate summary stats 
stats_by_panel <- dt_all_geneset[, {
  x <- mse_mean_metric[metric == "orig_DEL_mse"]
  y <- mse_mean_metric[metric == "insert_DEL_mse"]
  
  nx <- length(x); ny <- length(y)
  ks <- suppressWarnings(ks.test(x, y))
  mw <- suppressWarnings(wilcox.test(x, y))
  U  <- unname(mw$statistic)
  
  # Effect sizes derived from U
  # Probability of superiority: P(X > Y) + 0.5 * P(X == Y)
  pso       <- U / (nx * ny)
  # Cliff's delta: rescales PSO to [-1, 1]; 0 = stochastic equality
  cliffs_d  <- 2 * pso - 1
  # Rank-biserial correlation: same as Cliff's delta for two independent samples
  rank_bis  <- cliffs_d
  
  .(n_orig        = nx,
    n_insert      = ny,
    median_orig   = median(x, na.rm = TRUE),
    median_insert = median(y, na.rm = TRUE),
    ks_D          = unname(ks$statistic),
    ks_p          = ks$p.value,
    mw_U          = U,
    mw_p          = mw$p.value,
    pso           = pso,
    rank_biserial = rank_bis)
}, by = .(alu_orig_category, alu_insert_category)]

stats_by_panel[, ks_p_adj := p.adjust(ks_p, method = "BH")]
stats_by_panel[, mw_p_adj := p.adjust(mw_p, method = "BH")]


#pdf(paste0(DOWNLOADS_DIR, "/gene_swaps_test_line.pdf"), width=90/25.4, height=60/25.4)
pdf(paste0(DOWNLOADS_DIR, "/gene_swaps_test_line.pdf"), width=120/25.4, height=30/25.4)

#distribution over all genes -- use top 25 disruptive set 
p<-ggplot(dt_all_geneset, aes(x = mse_mean_metric, fill=metric)) +
#p<-ggplot(dt_all_geneset, aes(x = mse_mean_metric, fill=fill_group)) +
  geom_density(alpha = 0.5, lwd=0.3) + 
  #geom_histogram()+
  #scale_fill_identity()+
  
  scale_fill_manual(values = col_vals) + 
  facet_wrap(alu_orig_category ~ alu_insert_category,
             labeller = labeller(
               alu_orig_category   = c(high = "Native: high", low = "Native: low"),
               alu_insert_category = c(high = "Insert: high", low = "Insert: low")
             ), ncol=4) + 
  labs(title='Top 25 genes with high number of disruptive Alus') + 
  theme_classic() + theme_textsize_border
print(p)

#per gene 
for(gene_i in unique(dt_all_geneset$gene)){
  dt_gene<-dt_all_geneset[dt_all_geneset$gene==gene_i]
  p<-ggplot(dt_gene, aes(x = mse_mean_metric, fill=metric)) +
    geom_density(alpha = 0.5, lwd=0.3) + 
    #geom_histogram()+
    #scale_fill_identity()+
    scale_fill_manual(values = col_vals) + 
    # facet_grid(alu_orig_category ~ alu_insert_category,
    #            labeller = labeller(
    #   alu_orig_category   = c(high = "Native: high", low = "Native: low"),
    #   alu_insert_category = c(high = "Insert: high", low = "Insert: low")
    # )) + 
    facet_wrap(alu_orig_category ~ alu_insert_category,
               labeller = labeller(
                 alu_orig_category   = c(high = "Native: high", low = "Native: low"),
                 alu_insert_category = c(high = "Insert: high", low = "Insert: low")
               ), ncol=4) + 
    labs(title=gene_i) + 
    theme_classic() + theme_textsize_border 
  print(p)
}

dev.off()


###################

genes_set3 <- c('RBFOX1', 'CNTNAP2', 'CALN1', 'PRKN', 'WWOX', 'CACNA1A', 'TBCE', 'KBTBD11-OT1')  # fill in

gene_sets <- rbind(
  data.table(gene = genes_set1, gene_set = "set1"),
  data.table(gene = genes_set2, gene_set = "set2"),
  data.table(gene = genes_set3, gene_set = "set3")
)

# --- Modified per-gene MSE plot, pulling from alu_annot_short ------------

plot_gene_mse <- function(annot, gene_sets, out_file = NULL,
                          point_size = 0.8, point_alpha = 0.6, width = 10, height = 4) {
  
  all_genes <- sort(unique(gene_sets$gene))
  
  # Subset annotation to Alus in any of the gene sets
  d <- annot[single_gene_name %in% all_genes]
  
  # Order panels by gene_set then gene, so plots come out grouped
  panel_order <- gene_sets[, .(gene_set = paste(sort(unique(gene_set)), collapse = ",")),
                           by = gene][order(gene_set, gene)]
  
  plots <- lapply(panel_order$gene, function(g) {
    sets_for_g <- panel_order[gene == g, gene_set]
    ggplot(d[single_gene_name == g], aes(x = POS, y = mse_mean, colour = disruption)) +
      geom_point(alpha = point_alpha, size = point_size) +
      scale_colour_manual(values = c(high = "#C15C6D", low = "#6AA1C8"),
      #scale_colour_manual(values = c(high = "#B33E52", low = "#5254A3"),
                          na.value = "grey90") +
      ylim(0, 0.05)+
      labs(x = "Genomic position", y = "MSE mean",
           title = paste0(g, "  [", sets_for_g, "]")) +
      theme_classic() + theme_textsize_border + 
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())
  })
  names(plots) <- panel_order$gene
  
  if (!is.null(out_file)) {
    pdf(out_file, width = width, height = height)
    for (p in plots) print(p)
    dev.off()
  } else {
    for (p in plots) print(p)
  }
  
  invisible(plots)
}

plot_gene_mse(alu_annot_short, gene_sets,out_file = paste0(DOWNLOADS_DIR, "/20260729_all_genes_mse.pdf"), width=90/25.4, height=25/25.4)
