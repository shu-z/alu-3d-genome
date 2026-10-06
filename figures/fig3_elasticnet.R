# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

library(data.table)
library(ggplot2)
library(corrplot)
library(pals)
library(patchwork)

res_dir<-paste0(PROJECT_DIR, "/results/")
fig_dir<-paste0(PROJECT_DIR, "/figs/paper_figs/")



#USE FOR GENERAL TEXT SIZE THINGS 
theme_textsize<-theme(
  axis.title = element_text(size = 7),     # axis titles ~7 pt
  strip.text = element_text(size = 7),      # facet headers ~7 pt
  plot.title = element_text(size = 8),      # if you add a title
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 5),
  axis.text = element_text(size = 5),
  axis.line  = element_line(linewidth = 0.3),  # x & y axis lines
  axis.ticks = element_line(linewidth = 0.3)
)


###############################################
#feature corrs and overlap 
#alu_scores<-fread(paste0(res_dir, 'paper_results/20260311_merged_aluhg38_all_featureannot_NA_blacklist_filter.txt'))
alu_scores<-fread(paste0(res_dir, 'paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'))

#get 1Mb and mse columns
mb_cols <- grep("1Mb", colnames(alu_scores), value = TRUE)
#mb_cols <- grep("_1kb", colnames(alu_scores), value = TRUE)

cols_feat<-c('mse_mean', mb_cols)


substrings_to_remove <- c('30_', 'phyloP100_alu', 'phyloP100_1kb', 'phyloP100_10kb', 
                          'phyloP100_100kb', 'phyloP100_1Mb',
                          'constitutive', 'broadly_shared', 'lineage_restricted', 
                          'cell_type_specific', 'singleton')

pattern <- paste(substrings_to_remove, collapse = '|')
cols_use <- cols_feat[!grepl(pattern, cols_feat, ignore.case = TRUE)]


#make 1mb df and remove suffix
mb_feat<-alu_scores[, ..cols_use]
colnames(mb_feat) <- gsub("_1Mb", "", colnames(mb_feat))
#colnames(mb_feat) <- gsub("_1kb", "", colnames(mb_feat))


#remove double counted mapp
cor_mat <- cor(mb_feat, use = "pairwise.complete.obs", method='spearman')



#put mse_mean first and cluster the rest 
others <- setdiff(colnames(cor_mat), 'mse_mean')
hc <- hclust(as.dist(1 - cor_mat[others, others]), method = "ward.D2")
#reversing so most positive ones on the left
ordered_others <- rev(others[hc$order])

new_order <- c('mse_mean', ordered_others)
cor_mat_reordered <- cor_mat[new_order, new_order]

# Plot
pdf(paste0(PROJECT_DIR, "/figs/paper_figs/202604023_feature_corr_matrix_updatedfeat_1kb_spearman.pdf"), width = 140/25.4, height = 140/25.4)
my_cols <- colorRampPalette(c("#393B79", "white", "#B33E52"))(200)

corrplot(cor_mat_reordered,
  type='upper',
  col=my_cols,
  tl.col = "black",
  tl.cex = 0.4,
  addCoef.col = "black",
  number.cex = 0.4,
)

dev.off()


####################################################

#performance results of elastic net
enet_res<-fread(paste0(res_dir, '20260212_elasticnet/20260507_CVchromsplit_window_EN_performance_logMSE_H1.txt'))
enet_res$window <- factor(enet_res$window, levels = c("alu", "1kb", "10kb", "100kb", "1Mb"))




pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260508_enet_performance_updatedfeat_H1.pdf"), width = 75/25.4, height = 60/25.4)

col_val<-c('#B8C9A6', '#91AA74')
#cols_deg1 <- colorRampPalette(c("#AFCB90", "#344521"))(5)  
#cols_deg2 <- colorRampPalette(c("#C7CDE6", "#2C3763"))(5)  

p_r2 <- ggplot(enet_res, aes(x = window, y = r2, fill = factor(interaction_degree))) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = col_val) +
  labs(x = "Window", y = "R²", fill = "Interaction degree") +
  theme_classic() + theme_textsize

p_mse <- ggplot(enet_res, aes(x = window, y = mse, fill = factor(interaction_degree))) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = col_val) +
  labs(x = "Window", y = "MSE", fill = "Interaction degree") +
  theme_classic() + theme_textsize

cowplot::plot_grid(p_r2, p_mse, ncol=1)

dev.off()



#read in coef 
enet_coeff<-fread(paste0(res_dir, '20260212_elasticnet/20260507_CVchromsplit_window_EN_coeffs_logMSE_HFF.txt'))
#enet_coeff$window <- factor(enet_coeff$window, levels = c("alu", "1kb", "10kb", "100kb", "1Mb"))
#enet_coeff<-enet_coeff[enet_coeff$window=='1Mb']
#remove '1Mb' from the colnames 
enet_coeff$feature <- gsub("_1Mb", "", enet_coeff$feature)

deg1_coeff<-enet_coeff[enet_coeff$degree=='1']
deg2_coeff<-enet_coeff[enet_coeff$degree=='2']
deg1_coeff_nozero<-deg1_coeff[(deg1_coeff$coef!=0)]
deg2_coeff_nozero<-deg2_coeff[(deg2_coeff$coef!=0)]


#actually grab top 20 because otherwise there's too many for deg2
deg1_top20 <- deg1_coeff[order(-abs(deg1_coeff$coef)), ][1:20, ]
deg2_top20 <- deg2_coeff[order(-abs(deg2_coeff$coef)), ][1:20, ]


#replace spaces with ~
deg2_coeff_nozero$feature <- gsub(" ", " ~ ", deg2_coeff_nozero$feature)


pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260508_enet_coeff_1Mb_HFF.pdf"), width = 160/25.4, height = 80/25.4)

max_coef<-max(abs(enet_coeff$coef))
p_deg1 <- ggplot(deg1_coeff_nozero, aes(x = as.factor(reorder(feature, coef)), y = coef, fill=coef)) +
  geom_col(position = position_dodge(width = 0.8)) + coord_flip() + 
  scale_fill_gradientn(colours = c("#5254A3", '#ffffff', "#B33E52"), 
                       values = scales::rescale(c(-max_coef, 0, max_coef)),  # rescale to 0–1
                       limits = c(-max_coef, max_coef)) +       
  geom_hline(yintercept = 0, linewidth = 0.2) +
  labs(x = "Feature", y = "Coefficient Value") +
  theme_classic() + theme_textsize + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())

#p_deg1

p_deg2 <- ggplot(deg2_top20, aes(x = as.factor(reorder(feature, coef)), y = coef, fill=coef)) +
  geom_col(position = position_dodge(width = 0.8)) + coord_flip() + 
  scale_fill_gradientn(colours = c("#5254A3", '#ffffff', "#B33E52"), 
                       values = scales::rescale(c(-max_coef, 0, max_coef)),  # rescale to 0–1
                       limits = c(-max_coef, max_coef)) +       
  geom_hline(yintercept = 0, linewidth = 0.2) +
  labs(x = "Feature", y = "Coefficient Value") +
  theme_classic() + theme_textsize + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())

#p_deg2

cowplot::plot_grid(p_deg1, p_deg2, ncol=2, rel_widths=c(1.25,1.5))

dev.off()


###################################### repeat for logreg 

#performance results of elastic net
enet_res<-fread(paste0(res_dir, '20260212_elasticnet/20260507_logreg_summary_HFF.tsv'))

enet_res$interaction_degree<-rep(c(1, 2),5)
enet_res$window<-rep(c('alu', 'alu', '1kb', '1kb', '10kb', '10kb', '100kb', '100kb', '1Mb', '1Mb'))
enet_res$window <- factor(enet_res$window, levels = c("alu", "1kb", "10kb", "100kb", "1Mb"))

# Palettes per degree
pal_deg1 <- c('#CCD8BF', '#B8C9A6', '#A4BA8D', '#91AA74', '#7E9A5B')  # green
pal_deg2 <- c('#BFD0DE', '#A6B8C9', '#8DA0B4', '#74899F', '#5B728A')  # blue-grey

# Combined palette keyed by window_degree
combined_pal <- c(
  setNames(pal_deg1, paste(levels(enet_res$window), 1, sep = '_')),
  setNames(pal_deg2, paste(levels(enet_res$window), 2, sep = '_'))
)

enet_res$window_degree <- paste(enet_res$window, enet_res$interaction_degree, sep = '_')

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260320_logreg_performance_updatedfeat_1Mb.pdf"),width = 75/25.4, height = 60/25.4)

p_r2 <- ggplot(enet_res, aes(x = window, y = auroc, fill = window_degree)) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = combined_pal, guide = 'none') +
  labs(x = "Window", y = "AUROC") +
  theme_classic() + theme_textsize

p_mse <- ggplot(enet_res, aes(x = window, y = auprc, fill = window_degree)) +1
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = combined_pal, guide = 'none') +
  labs(x = "Window", y = "auprc") +
  theme_classic() + theme_textsize

cowplot::plot_grid(p_r2, p_mse, ncol = 1)
dev.off()



# Load both cell types and combine
load_ct <- function(ct) {
  roc <- fread(paste0(res_dir, '20260212_elasticnet/', '20260507_logreg_roc_curves_', ct, '.txt'))
  summary <- fread(paste0(res_dir, '20260212_elasticnet/', '20260507_logreg_summary_', ct, '.tsv')) %>%
    mutate(cell_type = ct)
  list(roc = roc, summary = summary)
}

hff <- load_ct('HFF')
h1  <- load_ct('H1')

#roc_all     <- bind_rows(hff$roc, h1$roc)
#summary_all <- bind_rows(hff$summary, h1$summary)
roc_all     <- bind_rows(h1$roc)
summary_all <- bind_rows(h1$summary)

# Build labels with AUC included
plot_df <- roc_all %>%
  left_join(summary_all, by = c('cell_type', 'feat_range', 'degree')) %>%
  mutate(label = sprintf('%s (AUC=%.3f)', feat_range, auroc),
         feat_range = factor(feat_range, levels = c('alu', '1kb', '10kb', '100kb', '1Mb')),
         feat_degree = interaction(feat_range, degree, sep = '_'))

pal_deg1 <- c('#CCD8BF', '#B8C9A6', '#A4BA8D', '#91AA74', '#7E9A5B')  # green
pal_deg2 <- c('#BFD0DE', '#A6B8C9', '#8DA0B4', '#74899F', '#5B728A')  # blue-grey

# Combined palette keyed by feat_range_degree
combined_pal <- c(setNames(pal_deg1, paste(levels(plot_df$feat_range), 1, sep = '_')),
                  setNames(pal_deg2, paste(levels(plot_df$feat_range), 2, sep = '_')))

legend_labels <- sub('_[12]$', '', names(combined_pal))

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260508_logreg_performance_1Mb_H1.pdf"),width = 60/25.4, height = 90/25.4)

ggplot(plot_df, aes(x = fpr, y = tpr, color = feat_degree)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.4) +
  facet_wrap(~ degree, ncol=1) + coord_equal() +
  scale_color_manual(values = combined_pal, labels = legend_labels, breaks = names(combined_pal),name = 'Feature range') +
  labs(x = 'False positive rate', y = 'True positive rate') +
  theme_classic() + theme_textsize + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())


dev.off()



#read in coef 
cell_type<-'HFF'
window<-'1Mb'
enet_coeff_all<-fread(paste0(res_dir, '20260212_elasticnet/20260507_logreg_coefs_', cell_type, '.txt'))
#enet_coeff$window <- factor(enet_coeff$window, levels = c("alu", "1kb", "10kb", "100kb", "1Mb"))

#subset to only 1Mb ones 
enet_coeff<-enet_coeff_all[enet_coeff_all$feat_range==window]
#remove '1Mb' from the colnames 
enet_coeff$feature <- gsub(paste0("_", window), "", enet_coeff$feature)

deg1_coeff<-enet_coeff[enet_coeff$degree=='1']
deg2_coeff<-enet_coeff[enet_coeff$degree=='2']
deg1_coeff_nozero<-deg1_coeff[(deg1_coeff$coef!=0)]
deg2_coeff_nozero<-deg2_coeff[(deg2_coeff$coef!=0)]


#actually grab top 20 because otherwise there's too many for deg2
deg1_top20 <- deg1_coeff[order(-abs(deg1_coeff$coef)), ][1:20, ]
deg2_top20 <- deg2_coeff[order(-abs(deg2_coeff$coef)), ][1:20, ] 


#replace spaces with ~
deg2_coeff_nozero$feature <- gsub(" ", " ~ ", deg2_coeff_nozero$feature)


pdf(paste0(fig_dir, "20260508_logreg_coeff_updatedfeat_", window, '_', cell_type, ".pdf"), width = 160/25.4, height = 80/25.4)

max_coef<-max(abs(enet_coeff$coef))
p_deg1 <- ggplot(deg1_coeff_nozero, aes(x = as.factor(reorder(feature, coef)), y = coef, fill=coef)) +
  geom_col(position = position_dodge(width = 0.8)) + coord_flip() + 
  scale_fill_gradientn(colours = c("#5254A3", '#ffffff', "#B33E52"), 
                       values = scales::rescale(c(-max_coef, 0, max_coef)),  # rescale to 0–1
                       limits = c(-max_coef, max_coef)) +       
  geom_hline(yintercept = 0, linewidth = 0.2) +
  labs(x = "Feature", y = "Coefficient Value") +
  theme_classic() + theme_textsize  +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())
#p_deg1

p_deg2 <- ggplot(deg2_top20, aes(x = as.factor(reorder(feature, coef)), y = coef, fill=coef)) +
  geom_col(position = position_dodge(width = 0.8)) + coord_flip() + 
  scale_fill_gradientn(colours = c("#5254A3", '#ffffff', "#B33E52"), 
                       values = scales::rescale(c(-max_coef, 0, max_coef)),  # rescale to 0–1
                       limits = c(-max_coef, max_coef)) +       
  geom_hline(yintercept = 0, linewidth = 0.2) +
  labs(x = "Feature", y = "Coefficient Value") +
  theme_classic() + theme_textsize + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank())

#p_deg2

cowplot::plot_grid(p_deg1, p_deg2, ncol=2, rel_widths=c(1.25,1.5))

dev.off()






########################################  plot individual feature distributions 


#get 1Mb and mse columns
cols_window <- grep("100kb", colnames(alu_scores), value = TRUE)
cols_plot<-c('mse_mean', cols_window)
mb_feat_all<-alu_scores[, ..cols_plot]

top_cut <- quantile(alu_scores$mse_mean, 0.99, na.rm = TRUE)
bottom_cut <- quantile(alu_scores$mse_mean, 0.50, na.rm = TRUE)
mb_feat_subset <- mb_feat_all[mse_mean >= top_cut | mse_mean <= bottom_cut] [, disruption := ifelse(mse_mean >= top_cut, "high", "low")]

dt_long <- melt(mb_feat_subset, id.vars = "disruption", value.name = "value")

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260423_feat_histogram_indivdiual_100kb.pdf"), width = 60/25.4, height = 50/25.4)


plot_feature <- function(feat){
  medians <- dt_long[, .(med = median(value, na.rm = TRUE)), by = .(variable, disruption)]
  median_feat<-medians[medians$variable==feat]
  
  p<-ggplot(dt_long[variable == feat], aes(x = value, fill = disruption)) +
    #geom_density(alpha = 0.5, adjust=1, lwd=0.4) +
    geom_histogram(aes(y=..density..), color='#292929', bins=50, lwd=0.15, alpha=0.9) + 
    facet_wrap(~disruption, ncol=1) + 
    scale_fill_manual(values=c('#7C985D', '#D7E0CC')) + 
    #scale_color_manual(values=c('#73864d', '#c4cbb4')) + 
    geom_vline(data = median_feat, aes(xintercept = med),color='#292929', linetype = "dashed", lwd=0.3) +
    labs(title = feat, x = feat, y = "Density", fill = "Disruption") + 
    theme_classic() + theme_textsize + 
    theme(strip.background = element_blank(),strip.text = element_blank())
  print(p)
}

plots <- lapply(colnames(mb_feat_subset)[2:27], plot_feature)

dev.off()




pdf(paste0(PROJECT_DIR, "/figs/paper_figs/feat_histogram_faceted.pdf"), width = 90/25.4, height = 60/25.4)

medians <- dt_long[, .(med = median(value, na.rm = TRUE)), by = .(variable, disruption)]

ggplot(dt_long, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians, aes(xintercept = med, color = group), linetype = "dashed") +
  facet_wrap(~variable, scales = "free") +
  theme_classic() + theme_textsize
dev.off()




#################################################################################
#plot specifically high/low phyloP vs disruption metrics 

alu_scores_annot<-fread(paste0(res_dir, 'paper_results/20260316_merged_aluhg38_all_featureannot_withCTCF_NA_blacklist_filter_region_alu.txt'))
#alu_scores_annot<-fread(paste0(res_dir, '/20260316_merged_aluhg38_all_featureannot_withCTCF_NA_blacklist_filter_region_alu.txt'))

#get 1Mb and mse columns
cols_window <- grep("1Mb", colnames(alu_scores_annot), value = TRUE)
#cols_plot<-c('mse_mean', 'orig_idx', 'single_gene_name', 'single_gene_type', cols_window)
cols_plot<-c('mse_mean', 'orig_idx', cols_window)

mb_feat_all<-alu_scores_annot[, ..cols_plot]

top_mse <- quantile(alu_scores_annot$mse_mean, 0.99, na.rm = TRUE)
bottom_mse <- quantile(alu_scores_annot$mse_mean, 0.50, na.rm = TRUE)
top_phyloPhigh <- quantile(alu_scores_annot$phyloP100_highfrac_1Mb, 0.8, na.rm = TRUE)
bottom_phyloPhigh <- quantile(alu_scores_annot$phyloP100_highfrac_1Mb, 0.20, na.rm = TRUE)
top_phyloPlow <- quantile(alu_scores_annot$phyloP100_lowfrac_1Mb, 0.8, na.rm = TRUE)
bottom_phyloPlow <- quantile(alu_scores_annot$phyloP100_lowfrac_1Mb, 0.20, na.rm = TRUE)

mse_subset <- mb_feat_all[mse_mean >= top_mse | mse_mean <= bottom_mse] [, disruption := ifelse(mse_mean >= top_mse, "high_mse", "low_mse")]
phyloPhigh_subset <- mse_subset[phyloP100_highfrac_1Mb >= top_phyloPhigh | phyloP100_highfrac_1Mb <= bottom_phyloPhigh] [, phyloPhigh_group := ifelse(phyloP100_highfrac_1Mb >= top_phyloPhigh, "high", "low")]
phyloPlow_subset <- mse_subset[phyloP100_lowfrac_1Mb >= top_phyloPlow | phyloP100_lowfrac_1Mb <= bottom_phyloPlow] [, phyloPlow_group := ifelse(phyloP100_lowfrac_1Mb >= top_phyloPlow, "high", "low")]


phyloPhigh_long <- melt(phyloPhigh_subset, id.vars = c("disruption", 'phyloPhigh_group'), value.name = "value")
phyloPlow_long <- melt(phyloPlow_subset, id.vars = c("disruption", 'phyloPlow_group'), value.name = "value")

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260426_Phighlow_corr_1Mb.pdf"), width = 120/25.4, height = 90/25.4)

phyloPhigh_high <- mse_subset[phyloP100_highfrac_1Mb >= top_phyloPhigh]
phyloPlow_high <- mse_subset[phyloP100_lowfrac_1Mb >= top_phyloPlow]

phyloPhigh_high[, group := "high"]
phyloPlow_high[, group := "low"]
phyloPhighlow <- rbind(phyloPhigh_high, phyloPlow_high, fill = TRUE)

# mark duplicates (rows appearing in both sets) as "both"
# assumes you have a unique row identifier — replace `id` with your key column(s)
phyloPhighlow[, group := ifelse(.N > 1, "both", group), by = orig_idx]
phyloPhighlow <- unique(phyloPhighlow, by = "orig_idx")

p <- ggplot(phyloPhighlow, aes(x = phyloP100_highfrac_1Mb, y = phyloP100_lowfrac_1Mb)) +
  #geom_point(aes(color = group), alpha = 0.4, size = 0.8) +
  stat_density_2d(aes(color = group), geom = "contour", n = 100, bins = 10) +
  scale_color_manual(values = c("high" = "#C0392B", "low" = "#2B3B4F", "both" = "#88A0BF")) +
  labs(x = "phyloP100_highfrac_1Mb",
       y = "phyloP100_lowfrac_1Mb",
       color = "group") +
  theme_classic() +
  theme_textsize
print(p)
dev.off()

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260426_feat_histogram_phyloPhighfacet_1Mb.pdf"), width = 120/25.4, height = 90/25.4)
  plot_feature <- function(feat){
    medians_high <- phyloPhigh_long[, .(med = stats::median(value, na.rm = TRUE)), by = .(variable, disruption, phyloPhigh_group)]
    median_high_feat<-medians_high[medians_high$variable==feat]
    
    p<-ggplot(phyloPhigh_long[variable == feat], aes(x = value, fill = disruption)) +
      #geom_density(alpha = 0.5, adjust=1, lwd=0.4) +
      geom_histogram(aes(y=..density..), color='#292929', bins=50, lwd=0.15, alpha=0.9) + 
      facet_wrap(phyloPhigh_group~disruption, ncol=2) + 
      scale_fill_manual(values=c('#7C985D', '#D7E0CC')) + 
      #scale_color_manual(values=c('#73864d', '#c4cbb4')) + 
      geom_vline(data = median_high_feat, aes(xintercept = med),color='#292929', linetype = "dashed", lwd=0.3) +
      labs(title = feat, x = feat, y = "Density", fill = "Disruption") + 
      theme_classic() + theme_textsize 
      #theme(strip.background = element_blank(),strip.text = element_blank())
    print(p)
  }
  plots <- lapply(colnames(mse_subset)[3:29], plot_feature)
dev.off()


pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260426_feat_histogram_phyloPlowfacet_1Mb.pdf"), width = 120/25.4, height = 90/25.4)
plot_feature <- function(feat){
  medians_low <- phyloPlow_long[, .(med = stats::median(value, na.rm = TRUE)), by = .(variable, disruption, phyloPlow_group)]
  median_low_feat<-medians_low[medians_low$variable==feat]
  
  p<-ggplot(phyloPlow_long[variable == feat], aes(x = value, fill = disruption)) +
    #geom_density(alpha = 0.5, adjust=1, lwd=0.4) +
    geom_histogram(aes(y=..density..), color='#292929', bins=50, lwd=0.15, alpha=0.9) + 
    facet_wrap(phyloPlow_group~disruption, ncol=2) + 
    scale_fill_manual(values=c('#7C985D', '#D7E0CC')) + 
    #scale_color_manual(values=c('#73864d', '#c4cbb4')) + 
    geom_vline(data = median_low_feat, aes(xintercept = med),color='#292929', linetype = "dashed", lwd=0.3) +
    labs(title = feat, x = feat, y = "Density", fill = "Disruption") + 
    theme_classic() + theme_textsize 
  #theme(strip.background = element_blank(),strip.text = element_blank())
  print(p)
}
plots <- lapply(colnames(mse_subset)[3:29], plot_feature)
dev.off()

graphics.off()





#################### look at high MSE alus only
high_mse <- mse_subset[disruption == "high_mse"]

# Assign phyloP groups on each axis (top 20% vs bottom 20%, mid dropped)
high_mse_phyloPhigh <- high_mse[
  phyloP100_highfrac_1Mb >= top_phyloPhigh | phyloP100_highfrac_1Mb <= bottom_phyloPhigh
][, group := ifelse(phyloP100_highfrac_1Mb >= top_phyloPhigh, "high", "low")
][, phyloP_axis := "phyloPhigh"]

high_mse_phyloPlow <- high_mse[
  phyloP100_lowfrac_1Mb >= top_phyloPlow | phyloP100_lowfrac_1Mb <= bottom_phyloPlow
][, group := ifelse(phyloP100_lowfrac_1Mb >= top_phyloPlow, "high", "low")
][, phyloP_axis := "phyloPlow"]

combined <- rbind(high_mse_phyloPhigh, high_mse_phyloPlow, fill = TRUE)

# Long format, excluding the phyloP/phastCon features themselves
feat_cols <- setdiff(cols_window, c('phyloP30_1Mb', 'phyloP100_1Mb',
                                    'phastCon30_1Mb', 'phastCon100_1Mb',
                                    'phyloP100_lowfrac_1Mb', 'phyloP100_highfrac_1Mb'))
combined_long <- melt(combined,
                      id.vars = c("group", "phyloP_axis"),
                      measure.vars = feat_cols,
                      variable.name = "feature",
                      value.name = "value")

################
# Rank-biserial: high vs low phyloP group, per feature × axis
# Positive r_rb => high-phyloP group has larger feature values

rank_biserial <- combined_long[, {
  hi <- value[group == "high"]; hi <- hi[!is.na(hi)]
  lo <- value[group == "low"];  lo <- lo[!is.na(lo)]
  if (length(hi) == 0 || length(lo) == 0) {
    .(r_rb = NA_real_, p = NA_real_, n_hi = length(hi), n_lo = length(lo))
  } else {
    w <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
    U <- unname(w$statistic)
    r_rb <- 2 * U / (length(hi) * length(lo)) - 1
    .(r_rb = r_rb, p = w$p.value,
      n_hi = length(hi), n_lo = length(lo))
  }
}, by = .(feature, phyloP_axis)]

rank_biserial[, p_adj := p.adjust(p, method = "BH"), by = phyloP_axis]

# Order features by max |r_rb| across the two axes
feat_order <- rank_biserial[, .(m = max(abs(r_rb), na.rm = TRUE)),
                            by = feature][order(m), feature]
rank_biserial[, feature := factor(feature, levels = feat_order)]

# rank-biserial dot plot, paired across phyloP axes

rank_biserial <- combined_long[group == "high", {
  hi <- value[phyloP_axis == "phyloPhigh"]; hi <- hi[!is.na(hi)]
  lo <- value[phyloP_axis == "phyloPlow"];  lo <- lo[!is.na(lo)]
  if (length(hi) == 0 || length(lo) == 0) {
    .(r_rb = NA_real_, p = NA_real_, n_hi = length(hi), n_lo = length(lo))
  } else {
    w <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
    U <- unname(w$statistic)
    r_rb <- 2 * U / (length(hi) * length(lo)) - 1
    .(r_rb = r_rb, p = w$p.value,
      n_hi = length(hi), n_lo = length(lo))
  }
}, by = .(feature)]

rank_biserial[, p_adj := p.adjust(p, method = "BH")]

feat_order <- rank_biserial[, .(m = max(abs(r_rb), na.rm = TRUE)),
                            by = feature][order(m), feature]
rank_biserial[, feature := factor(feature, levels = feat_order)]

#p_effect <- ggplot(rank_biserial,aes(x = r_rb, y = feature, color = phyloP_axis)) +
p_effect <- ggplot(rank_biserial,aes(x = r_rb, y = feature)) +
  geom_point(aes(shape = p_adj < 0.05), size = 2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "n.s."),name = NULL) +
  #scale_color_manual(values = c("phyloPhigh" = "#2C5F8D","phyloPlow"  = "#C97B3C")) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(x = "Rank-biserial r (high vs low phyloP, high-MSE Alus only)",
       y = NULL, color = "phyloP group") +
  theme_classic() + theme_textsize

# Plot 2: ECDFs for top axis-dependent features
#rb_wide <- dcast(rank_biserial, feature ~ phyloP_axis, value.var = "r_rb")
#rb_wide <- dcast(rank_biserial, feature ~ phyloP_axis, value.var = "r_rb")

rb_wide[, axis_diff := abs(phyloPhigh - phyloPlow)]
top_features <- rb_wide[order(-axis_diff)][1:6, as.character(feature)]
library(ggh4x)

p_ecdf_top <- ggplot(combined_long[feature %in% top_features][, feature := droplevels(feature)], aes(x = value, color = group)) +
  stat_ecdf(geom = "step", lwd = 0.5) +
  #facet_grid2(feature ~ phyloP_axis, scales = "free", independent = "all") +
  scale_color_manual(values = c("high" = "#7C985D", "low" = "#3a4a2c")) +
  labs(x='overlap', y = "Cumulative fraction", color = "phyloP overlap") +
  theme_classic() + theme_textsize +
  theme(
    strip.text = element_text(size = 4),
    strip.text.y.left = element_text(angle = 0)
  )

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260429_rankbiserial_highMSE_phyloPaxes_1Mb.pdf"),
    width = 90/25.4, height = 140/25.4)
print(p_effect)
print(p_ecdf_top)

dev.off()

graphics.off()


##################################################

# Rank-biserial: high-MSE vs low-MSE Alus, per feature
# Positive r_rb => high-MSE Alus have larger feature values

# Use the same mse_subset built earlier (1Mb features, high/low disruption groups)
# Long format over the 1Mb feature columns, excluding mse_mean and ids

cols_window <- grep("100kb", colnames(alu_scores_annot), value = TRUE)

CTCF_CTS <- grep("^(CTCF_|phyloP30_1|phyloP100_1|phastCon30)", cols_window, value = TRUE)

feat_cols_disr <- setdiff(cols_window,c(CTCF_CTS))
#cols_plot<-c('mse_mean', 'orig_idx', 'single_gene_name', 'single_gene_type', cols_window)
cols_plot<-c('mse_mean', 'orig_idx', cols_window)

mb_feat_all<-alu_scores_annot[, ..cols_plot]

top_mse <- quantile(alu_scores_annot$mse_mean, 0.99, na.rm = TRUE)
bottom_mse <- quantile(alu_scores_annot$mse_mean, 0.50, na.rm = TRUE)
mse_subset <- mb_feat_all[mse_mean >= top_mse | mse_mean <= bottom_mse] [, disruption := ifelse(mse_mean >= top_mse, "high_mse", "low_mse")]


disr_long <- melt(mse_subset,
                  id.vars = "disruption",
                  measure.vars = feat_cols_disr,
                  variable.name = "feature",
                  value.name = "value")

rb_disr <- disr_long[, {
  hi <- value[disruption == "high_mse"]; hi <- hi[!is.na(hi)]
  lo <- value[disruption == "low_mse"];  lo <- lo[!is.na(lo)]
  if (length(hi) == 0 || length(lo) == 0) {
    .(r_rb = NA_real_, p = NA_real_, n_hi = length(hi), n_lo = length(lo))
  } else {
    w <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
    U <- unname(w$statistic)
    n_hi <- as.numeric(length(hi))     # <- force double
    n_lo <- as.numeric(length(lo))
    r_rb <- 2 * U / (n_hi * n_lo) - 1
    .(r_rb = r_rb, p = w$p.value, n_hi = n_hi, n_lo = n_lo)
  }
}, by = feature]

rb_disr[, p_adj := p.adjust(p, method = "BH")]

# Strip _1Mb suffix for display, matching the coefficient plots
rb_disr[, feature := gsub("_100kb", "", feature)]

# Order features by r_rb so the plot reads top-positive to bottom-negative
rb_disr[, feature := factor(feature, levels = feature[order(r_rb)])]

# Symmetric color limits, same palette as the coefficient plots
max_rb <- max(abs(rb_disr$r_rb), na.rm = TRUE)

p_rb_disr <- ggplot(rb_disr, aes(x = r_rb, y = feature)) +
  # geom_segment(aes(x = 0, xend = r_rb, y = feature, yend = feature,
  #                  color = r_rb), linewidth = 0.4) +
  geom_point(aes(fill = r_rb, shape = p_adj < 0.05),
             size = 2, stroke = 0.3, color = "black") +
  scale_fill_gradientn(
    colours = c("#5254A3", "#ffffff", "#B33E52"),values  = scales::rescale(c(-max_rb, 0, max_rb)),
    limits  = c(-max_rb, max_rb),name    = "Rank-biserial r" ) +
  # scale_color_gradientn(
  #   colours = c("#5254A3", "#ffffff", "#B33E52"),
  #   values  = scales::rescale(c(-max_rb, 0, max_rb)),
  #   limits  = c(-max_rb, max_rb),
  #   guide   = "none"
  # ) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q<0.05", `FALSE` = "n.s."),
                     name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Rank-biserial r, highly disruptive vs. neutral Alus", y = NULL) +
  theme_classic() + theme_textsize +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank())


pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260429_rankbiserial_highMSE_100kb.pdf"),
    width = 90/25.4, height = 90/25.4)
print(p_rb_disr)
dev.off()

graphics.off()


##################################################
# Rank-biserial for all window ranges (faceted)

windows         <- c("alu",   "1kb",   "10kb",   "100kb",   "1Mb")
window_patterns <- c("_alu$", "_1kb$", "_10kb$", "_100kb$", "_1Mb$")

top_mse_all    <- quantile(alu_scores_annot$mse_mean, 0.99, na.rm = TRUE)
bottom_mse_all <- quantile(alu_scores_annot$mse_mean, 0.50, na.rm = TRUE)

rb_all <- rbindlist(lapply(seq_along(windows), function(i) {
  win        <- windows[i]
  pat        <- window_patterns[i]
  cols_win   <- grep(pat, colnames(alu_scores_annot), value = TRUE)
  if (length(cols_win) == 0) return(NULL)

  excl      <- grep("^(CTCF_|phyloP30_|phyloP100_|phastCon30_)", cols_win, value = TRUE)
  feat_cols <- setdiff(cols_win, excl)
  if (length(feat_cols) == 0) return(NULL)

  mb_sub  <- alu_scores_annot[, c("mse_mean", feat_cols), with = FALSE]
  mse_sub <- mb_sub[mse_mean >= top_mse_all | mse_mean <= bottom_mse_all][
    , disruption := ifelse(mse_mean >= top_mse_all, "high_mse", "low_mse")]

  dt_long <- melt(mse_sub, id.vars = "disruption", measure.vars = feat_cols,
                  variable.name = "feature", value.name = "value")

  rb <- dt_long[, {
    hi <- value[disruption == "high_mse"]; hi <- hi[!is.na(hi)]
    lo <- value[disruption == "low_mse"];  lo <- lo[!is.na(lo)]
    if (length(hi) == 0 || length(lo) == 0) {
      .(r_rb = NA_real_, p = NA_real_, n_hi = length(hi), n_lo = length(lo))
    } else {
      w   <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
      U   <- as.numeric(unname(w$statistic))
      r_rb <- 2 * U / (as.numeric(length(hi)) * as.numeric(length(lo))) - 1
      .(r_rb = r_rb, p = w$p.value,
        n_hi = as.numeric(length(hi)), n_lo = as.numeric(length(lo)))
    }
  }, by = feature]

  rb[, p_adj   := p.adjust(p, method = "BH")]
  rb[, feature := gsub(pat, "", as.character(feature))]
  rb[, window  := win]
  rb
}), fill = TRUE)

rb_all[, window := factor(window, levels = windows)]

max_rb <- max(abs(rb_all$r_rb), na.rm = TRUE)

rb_subdir <- paste0(fig_dir, "rankbiserial/")
dir.create(rb_subdir, showWarnings = FALSE, recursive = TRUE)

for (win in windows) {
  rb_win <- rb_all[window == win]
  if (nrow(rb_win) == 0) next

  rb_win[, feature := factor(feature, levels = feature[order(r_rb)])]

  p_win <- ggplot(rb_win, aes(x = r_rb, y = feature)) +
    geom_point(aes(fill = r_rb, shape = p_adj < 0.05),
               size = 1.5, stroke = 0.3, color = "black") +
    scale_fill_gradientn(
      colours = c("#5254A3", "#ffffff", "#B33E52"),
      values  = scales::rescale(c(-max_rb, 0, max_rb)),
      limits  = c(-max_rb, max_rb), name = "Rank-biserial r") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                       labels = c(`TRUE` = "BH q<0.05", `FALSE` = "n.s."), name = NULL) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(title = win,
         x = "Rank-biserial r, highly disruptive vs. neutral Alus", y = NULL) +
    theme_classic() + theme_textsize +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.line = element_blank())

  pdf(paste0(rb_subdir, "20260711_rankbiserial_", win, "_width75.pdf"),
      width = 75/25.4, height = 65/25.4)
  print(p_win)
  dev.off()
}

graphics.off()
