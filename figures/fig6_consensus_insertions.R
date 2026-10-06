# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")
DOWNLOADS_DIR <- "/Users/shu/Downloads"

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr) 
library(pals)
library(dplyr)

stepped_colors<-stepped(24)
stepped2_colors<-stepped2(24)
stepped3_colors<-stepped3(24)

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


res_dir<-paste0(PROJECT_DIR, "/results/paper_results/")
fig_dir<-paste0(PROJECT_DIR, "/figs/paper_figs/")
ancestral_alu <- c("FLAM", "FRAM", "FAM", "FLAM_A", "FLAM_C")


alu_score_all<-fread(paste0(res_dir, '20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'))
#alu_ins<-fread(paste0(PROJECT_DIR, "/results/20251022_topalu_insertAlu.txt"))
#alu_ins<-fread(paste0(PROJECT_DIR, "/results/paper_results/20260414_highlow1000_insertconsensus.txt"))
alu_ins<-fread(paste0(PROJECT_DIR, "/results/paper_results/20260414_highlow1000_insertconsensus_withDEL.txt"))
alu_ins[ ,SVLEN:=END-POS]
#rename old mse columns 
cols_to_rename <-  c("mse_mean", "corr_mean", 'mse_median', 'corr_median', "mse_high", "corr_high")
setnames(alu_ins, old = cols_to_rename, new = paste0("DEL_", cols_to_rename), skip_absent=T)

#fix typo on the DEL Columns
setnames(alu_ins, gsub("__", "_", names(alu_ins)))
setnames(alu_ins, gsub("mse_DEL", "mseDEL", names(alu_ins)))
setnames(alu_ins, gsub("corr_DEL", "corrDEL", names(alu_ins)))

#keep meta columns 
#meta_cols <- c("orig_idx", 'repName', 'DEL_mse_mean', 'DEL_corr_mean', 'DEL_mse_high', 'DEL_corr_high', "SVLEN")
meta_cols <- c("orig_idx", 'repName', 'DEL_mse_mean', 'DEL_corr_mean', 'DEL_mse_high', "SVLEN")


# Identify mse and corr columns
score_cols <- grep("^(mse|corr|mseDEL|corrDEL)_", names(alu_ins), value = TRUE)
rc_cols    <- grep("^(mse|corr|mseDEL|corrDEL)_RC_", score_cols, value = TRUE)
fwd_cols   <- sub("_RC_", "_", rc_cols)
#names for mean columns
mean_cols  <- sub("^(mse|corr|mseDEL|corrDEL)_", "\\1_mean_", fwd_cols)
align_cols<-grep("^(alignscore|percent_identity|gaps|)_", names(alu_ins), value = TRUE)

#compuet means, keep only mean columns and align ones 
alu_ins[, (mean_cols) := Map(function(f, r) rowMeans(.SD[, c(f, r), with = FALSE], na.rm = TRUE), fwd_cols, rc_cols)]
alu_ins_means <- alu_ins[, c(meta_cols, mean_cols, align_cols), with = FALSE]


# Grab the mean column names by metric type
mse_mean_cols     <- grep("^mse_mean_",     names(alu_ins_means), value = TRUE)
mseDEL_mean_cols  <- grep("^mseDEL_mean_",  names(alu_ins_means), value = TRUE)
corr_mean_cols    <- grep("^corr_mean_",    names(alu_ins_means), value = TRUE)
corrDEL_mean_cols <- grep("^corrDEL_mean_", names(alu_ins_means), value = TRUE)

# Sort each by window suffix so they align
get_window <- function(x) sub("^[^_]+_mean_", "", x)
mse_mean_cols     <- mse_mean_cols[    order(get_window(mse_mean_cols))]
mseDEL_mean_cols  <- mseDEL_mean_cols[ order(get_window(mseDEL_mean_cols))]
corr_mean_cols    <- corr_mean_cols[   order(get_window(corr_mean_cols))]
corrDEL_mean_cols <- corrDEL_mean_cols[order(get_window(corrDEL_mean_cols))]

# Melt mse and its DEL counterpart
mse_long <- melt(alu_ins_means, id.vars = meta_cols,
                 measure.vars = list(value = mse_mean_cols, value_DEL = mseDEL_mean_cols),variable.name = "window_idx")
mse_long[, metric_type := "mse"]
mse_long[, family := get_window(mse_mean_cols)[window_idx]]

# Same for corr
corr_long <- melt(alu_ins_means,id.vars = meta_cols, 
  measure.vars = list(value = corr_mean_cols, value_DEL = corrDEL_mean_cols), variable.name = "window_idx" )
corr_long[, metric_type := "corr"]
corr_long[, family := get_window(corr_mean_cols)[window_idx]]

alu_ins_long <- rbindlist(list(mse_long, corr_long), use.names = TRUE)



#pull out alignscore -- this time percent identity 
align_cols_save <- grep("^percent_identity_", names(alu_ins), value = TRUE)
align_dt <- melt(alu_ins, id.vars = meta_cols, measure.vars = align_cols_save, 
  variable.name = "align_metric",value.name = "alignscore")[, family := tstrsplit(align_metric, "_", keep = 3)]

#pull out gaps
gaps_alunew_save <- grep("^gaps_alunew", names(alu_ins), value = TRUE)
gaps_alunew_dt <- melt(alu_ins, id.vars = meta_cols, measure.vars = gaps_alunew_save, 
                 variable.name = "gaps_metric",value.name = "gaps_alunew")[, family := tstrsplit(gaps_metric, "_", keep = 3)]
gaps_alucons_save <- grep("^gaps_alucons", names(alu_ins), value = TRUE)
gaps_alucons_dt <- melt(alu_ins, id.vars = meta_cols, measure.vars = gaps_alucons_save, 
                variable.name = "gaps_metric",value.name = "gaps_alucons")[, family := tstrsplit(gaps_metric, "_", keep = 3)]

#merge mse/corr with alignment information 
alu_ins_long <- merge(alu_ins_long, align_dt[, .(orig_idx, family, alignscore)],by = c("orig_idx", "family"), all.x = TRUE)
alu_ins_long <- merge(alu_ins_long, gaps_alunew_dt[, .(orig_idx, family, gaps_alunew)],by = c("orig_idx", "family"), all.x = TRUE)
alu_ins_long <- merge(alu_ins_long, gaps_alucons_dt[, .(orig_idx, family, gaps_alucons)],by = c("orig_idx", "family"), all.x = TRUE)

#edit corr so its 1-corr
alu_ins_long[ ,value:=(ifelse (metric_type=='corr', (1-value), value))]


#add in information on if the exact alus are the same, or if the families match
alu_ins_long[ ,repName_short:=(substr(repName, 1, 4))]
alu_ins_long[ ,family_short:=(substr(family, 1, 4))]
alu_ins_long[ ,repName_family_match:=ifelse(repName==family, 1, 0)]
alu_ins_long[ ,repName_family_short_match:=ifelse(repName_short==family_short, 1, 0)]

alu_ins_long$repName_family_match<-as.factor(alu_ins_long$repName_family_match)
alu_ins_long$repName_family_short_match<-as.factor(alu_ins_long$repName_family_short_match)


#add in diff of how much the DEL score is higher 
alu_ins_long[, DEL_diff:=value_DEL/value]


#merge GC in 
alu_gc <- fread(paste0(DOWNLOADS_DIR, "/alu_dfam_consensus_with_gc.txt"))
alu_ins_long <- merge(alu_ins_long, alu_gc, by.x = "family", by.y = "AluName", all.x = TRUE)
alu_ins_long <- merge(alu_ins_long, alu_score_all[,c('orig_idx', 'GC_alu')], by='orig_idx', all.x = TRUE)




#define colors for repName and family
unique_repName<-unique(alu_score_all$repName)
repName<-unique_repName[order(unique_repName)]
repName_color<-c(stepped_colors[22], brewer.ylorbr(6)[2:5], brewer.greens(20)[4:20], brewer.pubu(32)[8:32], brewer.rdpu(8)[5:8])
levels1<-c(repName[48:51], repName[1:47])
repName_color_df<-data.table(repName, repName_color)
family_color_df<-data.table('family'=repName, 'family_color'=repName_color)

#merge and factor colors 
alus_plot_0<-merge(alu_ins_long, repName_color_df, by='repName')
alus_plot<-merge(alus_plot_0, family_color_df, by='family')
alus_plot$family<-factor(alus_plot$family, levels=levels1)
alus_plot$repName<-factor(alus_plot$repName, levels=levels1)



#split out mse and corr dts separately 
alus_plot_mse<-alus_plot[alus_plot$metric_type=='mse']
alus_plot_corr<-alus_plot[alus_plot$metric_type=='corr']
alus_plot_mse_noanc <- alus_plot_mse[!repName %in% ancestral_alu & !family %in% ancestral_alu]
alus_plot_corr_noanc <- alus_plot_corr[!repName %in% ancestral_alu & !family %in% ancestral_alu]

fwrite(alus_plot, paste0(PROJECT_DIR, "/results/20251022_topalu_insertAlu_longmelt_ancesremoved.txt"), sep='\t')

################################################## plotting

#remove ancestral alus 
#alus_plot <- alus_plot[!repName %in% ancestral_alu & !family %in% ancestral_alu]


alus_plot<-alus_plot_mse_noanc[alus_plot_mse_noanc$DEL_mse_high==1]

#remove alus where initial alu is <275 bp
alus_plot<-alus_plot[alus_plot$SVLEN>=275]


metric_to_plot <- "mse"  # mse or "corr"
outliersize=0.3
linewidth_box<-0.3
linewidth_mean<-0.4



#plot distribution of all disruption 

#plot boxplot of insertion scores by alu name 
box1<-ggplot(alus_plot[metric_type == metric_to_plot],aes(x = family, y = value, color=family_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() +  
  scale_color_identity() +
  geom_hline(data = alus_plot %>% summarize(median_mse = median(value, na.rm = TRUE)),
             aes(yintercept = median_mse), color = "red", linetype = "dashed") +
  labs(x = "Insert Consensus vs REF", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize

box2<-ggplot(alus_plot[metric_type == metric_to_plot],aes(x = family, y = value_DEL, color=family_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() +  
  scale_color_identity() +
  geom_hline(data = alus_plot %>% summarize(median_mse = median(value_DEL, na.rm = TRUE)),
             aes(yintercept = median_mse), color = "red", linetype = "dashed") +
  labs(x = "Insert Consensus vs Native Alu DEL", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize

box3<-ggplot(alus_plot[metric_type == metric_to_plot],aes(x = family, y = DEL_diff, color=family_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() +  
  scale_color_identity() +
  geom_hline(data = alus_plot %>% summarize(median_mse = median(value_DEL, na.rm = TRUE)),
             aes(yintercept = median_mse), color = "red", linetype = "dashed") +
  labs(x = "Insert Consensus vs Native Alu DEL", y = toupper(metric_to_plot)) + 
  scale_y_log10()+
  theme_bw() + theme_textsize


hist1<-ggplot(alus_plot_mse_noanc[metric_type == metric_to_plot], aes(x = value)) +
  geom_histogram() +  facet_wrap(~DEL_mse_high, ncol=1) +
  labs(x = "Histogram of Inserted Alu disruption scores", y="") + 
  theme_bw() + theme_textsize

hist2<-ggplot(alus_plot_mse_noanc[metric_type == metric_to_plot], aes(x = value_DEL)) +
  geom_histogram() + facet_wrap(~DEL_mse_high, ncol=1) +
  labs(x = "Histogram of Inserted Alu disruption scores", y="") + 
  theme_bw() + theme_textsize

pdf(paste0(fig_dir, '20260506_alu_consensusINS_', metric_to_plot, '_noanc.pdf'),  width=120/25.4, height=80/25.4)
cowplot::plot_grid(box1, box2, box3, ncol=3)
cowplot::plot_grid(hist1, hist2, ncol=2)
dev.off()

pdf(paste0(fig_dir, '20260506_alu_consensusINS_scatter', metric_to_plot, '_noanc.pdf'),  width=240/25.4, height=240/25.4)
#plot scatter of insertion scores whether compared to REF or DEL of original Alu 
ggplot(alus_plot[metric_type == metric_to_plot],aes(x = value, y = value_DEL, color=family_color)) + 
  geom_point(size=0.3, outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() + 
  facet_wrap(~family, ncol=6) +
  scale_color_identity() +
  labs(x = "Insert Alu", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize
dev.off()



pdf(paste0(fig_dir, '20260506_alu_consensusINS_features_', metric_to_plot, '_noances.pdf'),  width=12, height=10)

#look at correlation between alignment and alu insertion score 
ggplot(alus_plot[metric_type == metric_to_plot], aes(x = alignscore, y = value)) +
  geom_point(size=0.5) + facet_wrap(~family)+
  stat_cor(method = "pearson", size=3) + 
  labs(x = "Global Alignment b/w Inserted Consensus Alu and Original Alu", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize

ggplot(alus_plot[metric_type == metric_to_plot], aes(x = gaps_alunew, y = value)) +
  geom_point(size=0.5) + facet_wrap(~family)+
  stat_cor(method = "pearson", size=3) + 
  labs(x = "Gaps in original Alu after Global Alignment b/w Inserted Consensus Alu and Original Alu", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize

ggplot(alus_plot[metric_type == metric_to_plot], aes(x = gaps_alucons, y = value)) +
  geom_point(size=0.5) + facet_wrap(~family)+
  stat_cor(method = "pearson", size=3) + 
  labs(x = "Gaps in consensus Alu after Global Alignment b/w Inserted Consensus Alu and Original Alu", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize

ggplot(alus_plot[metric_type == metric_to_plot], aes(x = SVLEN, y = value)) +
  geom_point(size=0.5) + facet_wrap(~family)+
  stat_cor(method = "pearson", size=3) + 
  labs(x = "SV Length", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize


####### look at distribution of matching/non matching alu types 
#add annotation for counts in each group 
counts_family <- alus_plot[metric_type == metric_to_plot] %>%
  group_by(family, repName_family_match) %>%
  summarise(n = n(), .groups = "drop")
# consistent y-position across all facets
y_top <- max(alus_plot[metric_type == metric_to_plot]$value, na.rm = TRUE) * 1.05

ggplot(alus_plot[metric_type == metric_to_plot], 
       aes(x = repName_family_match, y = value, color=family_color)) +
  geom_boxplot(size = 0.5) +  facet_wrap(~family) + scale_color_identity() + 
  geom_text(data = counts_family, aes(x = repName_family_match, y = y_top, label = paste0("n=", n)), 
            inherit.aes=F, hjust = 1.1, vjust = 1, size = 3) +
  labs(x = "Original/Inserted Alu matched on Alu name, facet on inserted Alu", y = toupper(metric_to_plot)) +
  coord_cartesian(clip = "off") + 
  theme_bw() + theme_textsize


####### look at distribution of matching/non matching alu types 
#add annotation for counts in each group 
counts_repName <- alus_plot[metric_type == metric_to_plot] %>%
  group_by(repName, repName_family_match) %>%
  summarise(n = n(), .groups = "drop")
# consistent y-position across all facets
y_top <- max(alus_plot[metric_type == metric_to_plot]$value, na.rm = TRUE) * 1.05
ggplot(alus_plot[metric_type == metric_to_plot], aes(x = repName_family_match, y = value, color=repName_color)) +
  geom_boxplot(size=0.5) + facet_wrap(~repName) + scale_color_identity() + 
  geom_text(data = counts_repName, aes(x = repName_family_match, y = y_top, label = paste0("n=", n)), 
            inherit.aes=F, hjust = 1.1, vjust = 1, size = 3) +
  labs(x = "Original/Inserted Alu matched on Alu name, facet on original Alu", y = toupper(metric_to_plot)) + 
  coord_cartesian(clip = "off") + 
  theme_bw() + theme_textsize


####### 
#repeat for larger Alu families
ggplot(alus_plot[metric_type == metric_to_plot], aes(x = repName_family_short_match, y = value)) +
  geom_boxplot(size=0.5) + facet_wrap(~family_short)+
  labs(x = "matched on broader family, facet wrap on inserted Alu type", y = toupper(metric_to_plot)) + 
  theme_bw() + theme_textsize
dev.off()







############## look at alus with most amount of variance across different insertions 
# compute SD across all inserted alus in each original alu
rank_var<- 'mse' #mse or corr

alus_query<-alus_plot_mse_noanc[alus_plot_mse_noanc$DEL_mse_high==1]
alus_query<-alus_query[alus_query$SVLEN>=275]


if (rank_var=='mse'){
  alus_query_mse<-alus_query[metric_type == metric_to_plot]
  var_rank <- alus_query_mse[, .(var_value = sd(value, na.rm = TRUE)), by = .(repName, DEL_mse_mean, orig_idx, SVLEN)]
  # pick top 50
  top_sets <- var_rank[order(-var_value)][1:50]
  dt_top <- alus_query_mse[top_sets, on = .(repName, DEL_mse_mean)]
  dt_top[, interaction_id := paste(repName, 'mse',  signif(DEL_mse_mean, 4), sep = "_")]
  
}
if (rank_var=='corr'){
  var_rank <- alus_query[, .(var_value = sd(value, na.rm = TRUE)), by = .(repName, DEL_corr_mean)]
  # pick top 50
  top_sets <- var_rank[order(-var_value)][1:50]
  dt_top <- alus_query[top_sets, on = .(repName, DEL_corr_mean)]
  dt_top[, interaction_id := paste(repName, 'corr',  signif(DEL_corr_mean, 4), sep = "_")]
  
}


# plot outliers 
extreme<-0.5 #default is 1.5 
dt_top[, outlier := FALSE]  # initialize
dt_top[, outlier := value < quantile(value, 0.25) - 1.5*IQR(value) |
         value > quantile(value, 0.75) + extreme*IQR(value), 
       by = interaction_id]


pdf(paste0(fig_dir, '20260506_alu_consensusINS_top50', rank_var, '_SD.pdf'), width=8, height=10)

ggplot(var_rank, aes(x = var_value)) +
  geom_histogram() +
  labs(x = "SD across all Alu Insertions", y="") + theme_bw()

#look at top Alu sites with most variable scores based on what Alu is inserted 
ggplot(dt_top, aes(x = interaction_id, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = family_color), width = 0.4, size = 0.8) +
  scale_color_identity() +  
  geom_text_repel(data = dt_top[outlier == TRUE], aes(label = family),
    size = 2.5, nudge_x = 0.2) +
  geom_point(
    data = dt_top[repName_family_match == 1],
    shape = 8, size = 3, color = "red"
  ) +
  theme_bw() + coord_flip()  + 
  labs(x = paste0("Alu DEL ", rank_var), y=paste0("Alu Insertion ", rank_var)) + theme_bw()
dev.off()

pdf(paste0(fig_dir, '20260506_alu_consensusINS_top50', rank_var, '_SD_alignmentcorrs.pdf'), width=12, height=10)
ggplot(dt_top, aes(x = alignscore, y = value, color=family_color)) +
  geom_point(size=2) + facet_wrap(~interaction_id) +
  #stat_cor(method = "pearson", size=3) +
  stat_cor(method = "pearson", size = 3) + scale_color_identity()+
  labs(x = "Global Alignment b/w Inserted Alu and Original Alu", y = toupper(rank_var)) + theme_bw()


dev.off()


######################################################### 
#finally, annotate for phylogenetic distance -- plot distance vs mse_mean

#read in distance matrix
consensus_dist<-fread(paste0(PROJECT_DIR, "/results/phylo/alu_distances_MLE_bestfit.csv"), check.names=T)
distance_mat <- as.matrix(consensus_dist[, -1], rownames=consensus_dist[[1]])


rank_var<- 'mse' #mse or corr 
alus_query<-alus_plot_mse_noanc[alus_plot_mse_noanc$DEL_mse_high==1]
alus_query<-alus_query[alus_query$SVLEN>=250]


#look at phylo dist across all
alus_query[, 'phylo_dist' := mapply(function(r, c) {
  if (r %in% rownames(distance_mat) & c %in% colnames(distance_mat)) {
    distance_mat[r, c]} else {NA_real_}}, 
  as.character(family), as.character(repName))]

pdf(paste0(fig_dir, 'alu_ins_phylodist.pdf'), width = 85/25.4, height = 55/25.4)

p1<-ggplot(alus_query, aes(phylo_dist, value)) +
  geom_hex(bins = 100) +
  scale_fill_gradientn(colours = c("#2B3B4F", '#88A0BF', "#ffffff"), trans = "log10", name = "log10(count)") +
  labs(fill = "log10(count)", x='log10(MSE)', y='log10(CORR)') +
  #scale_x_log10() + scale_y_log10() + 
  theme_classic() + theme_textsize
p1
dev.off()


########################################################

#### split out for top most variable ones 


top_n<-100

if (rank_var=='mse'){
  alus_query_mse<-alus_query[metric_type == metric_to_plot]
  var_rank <- alus_query_mse[, .(var_value = sd(value, na.rm = TRUE)), by = .(repName, DEL_mse_mean)]
  # pick top 50
  top_sets <- var_rank[order(-var_value)][1:top_n]
  dt_top <- alus_query_mse[top_sets, on = .(repName, DEL_mse_mean)]
  dt_top[, interaction_id := paste(repName, 'mse',  signif(DEL_mse_mean, 4), orig_idx, sep = "_")]
  
}
if (rank_var=='corr'){
  var_rank <- alus_query[, .(var_value = sd(value, na.rm = TRUE)), by = .(repName, DEL_corr_mean)]
  # pick top 50
  top_sets <- var_rank[order(-var_value)][1:top_n]
  dt_top <- alus_query[top_sets, on = .(repName, DEL_corr_mean)]
  dt_top[, interaction_id := paste(repName, 'corr',  signif(DEL_corr_mean, 4), orig_idx, sep = "_")]
  
}


#add in distance metric
#some alus don't have consensus seq -- if it's the original Alu, report as NA 



dt_top[, 'phylo_dist' := mapply(function(r, c) {
  if (r %in% rownames(distance_mat) & c %in% colnames(distance_mat)) {
    distance_mat[r, c]} else {NA_real_}}, 
  as.character(family), as.character(repName))]


# plot outliers 
extreme<-0.3 #default is 1.5 
dt_top[, outlier := FALSE]  # initialize
dt_top[, outlier := value < quantile(value, 0.25) - 1.5*IQR(value) |
         value > quantile(value, 0.75) + extreme*IQR(value), by = interaction_id]

#calculate corr
dt_top_corr<-dt_top[, .(cor_phylo = cor(phylo_dist, value, use = "complete.obs")), by = interaction_id]
dt_top_corr<-dt_top[, .(cor_align = cor(alignscore, value, use = "complete.obs")), by = interaction_id]


pdf(paste0(fig_dir, '20260616_alu_consensusINS_top50', rank_var, '_SD_phylodist_noances_highdisruption_morelabel_phylo.pdf'), width = 80/25.4, height = 55/25.4)

#look at top Alu sites with most variable scores based on what Alu is inserted 
ggplot(dt_top, aes(x = alignscore, y = value)) +
  #geom_point(aes(color = family_color),size = 1) +
  #stat_density_2d(aes(color = after_stat(level)), geom = "contour", n=100, bins=20) +
  #scale_color_gradientn(colours = c("#2B3B4F", '#88A0BF', "#f6f6f6"), name = "log10(count)") +
  geom_hex(bins = 50) +
  scale_fill_gradientn(colours = c("#2B3B4F", '#88A0BF', "#ffffff"), trans = "log10", name = "log10(count)") +
  #scale_x_log10() + scale_y_log10() +
  theme_classic() + theme_textsize + 
  labs(x = paste0("Phylogenetic Distance"), y=paste0("MSE INS vs REF")) 
dev.off()


pdf(paste0(fig_dir, '20260618_alu_consensusINS_top50', rank_var, '_SD_phylodist_noances_highdisruption_morelabel.pdf'), width=28, height=20)
#look at top Alu sites with most variable scores based on what Alu is inserted 
ggplot(dt_top, aes(x = phylo_dist, y = value)) +
  geom_point(aes(color = family_color),size = 1) + facet_wrap(~interaction_id, ncol=10) +
  scale_color_identity() +   scale_fill_identity()+
  geom_text_repel(data = dt_top[outlier == TRUE], aes(label = family),size = 2.5, nudge_x = 0.05) +
  geom_point(data = dt_top[repName_family_match == 1], aes(fill = family_color), color='red', shape = 8, size = 2) +
  labs(x = paste0("Phylogenetic Distance between Native and Insert Consensus Alus", rank_var), y=paste0("Alu Insertion ", rank_var)) +
  theme_bw() + theme_textsize 

ggplot(dt_top, aes(x = alignscore, y = value)) +
  geom_point(aes(color = family_color),size = 1) + facet_wrap(~interaction_id, ncol=10) +
  scale_color_identity() +   scale_fill_identity()+
  geom_text_repel(data = dt_top[outlier == TRUE], aes(label = family),size = 2.5, nudge_x = 0.05) +
  geom_point(data = dt_top[repName_family_match == 1], aes(fill = family_color), color='red', shape = 8, size = 2) +
  labs(x = paste0("Align Score between Native Alu and Consensus Insert ", rank_var), y=paste0("Alu Insertion ", rank_var)) +
  theme_bw() + theme_textsize 
dev.off()



pdf(paste0(fig_dir, '20260618_consensus_ins_individual_phylodist.pdf'), width=60/25.4, height=90/25.4)
dt_top_query<-dt_top[dt_top$orig_idx %in% c(308144, 330025)]

ggplot(dt_top_query, aes(x = phylo_dist, y = value)) +
  geom_point(aes(color = family_color),size = 1) + facet_wrap(~interaction_id, ncol=1) +
  scale_color_identity() +   scale_fill_identity()+
  geom_text_repel(data = dt_top_query[outlier == TRUE], aes(label = family),size = 2.5, nudge_x = 0.05) +
  geom_point(data = dt_top_query[repName_family_match == 1], aes(fill = family_color), color='red', shape = 8, size = 2) +
  labs(x = paste0("Phylogenetic Distance between Native and Insert Consensus Alus", rank_var), y=paste0("Alu Insertion ", rank_var)) +
  theme_bw() + theme_textsize 
dev.off()


#############################make sure everything is closed!
graphics.off()
