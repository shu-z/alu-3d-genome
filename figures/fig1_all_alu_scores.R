# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

library(data.table)
library(ggplot2)
library(pals)
library(dplyr)
library(tidyr)
library(cowplot)


stepped_colors<-stepped(24)
stepped2_colors<-stepped2(24)
stepped3_colors<-stepped3(24)


#res_dir<-paste0(PROJECT_DIR, "/results/")
res_dir<-paste0(PROJECT_DIR, "/results/paper_results/")

#fig_dir<-paste0(PROJECT_DIR, "/figs/bmi_rips_remake/")
fig_dir<-paste0(PROJECT_DIR, "/figs/paper_figs/")



#all_alus<-fread(paste0(res_dir, 'aluhg38_all_annot_GC_mapp.txt'))
all_alus<-fread(paste0(res_dir, '20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'))
#jank way to remove duplicates
dup_cols<-!duplicated(names(all_alus))
all_alus <- all_alus[,..dup_cols]

all_alus$SVTYPE <-'DEL'


#edit corr
all_alus$corr_score_mean<-(1-all_alus$corr_mean)
all_alus$corr_score_median<-(1-all_alus$corr_median)

cor(all_alus$mse_mean, all_alus$corr_score_mean, use='complete.obs', method='pearson')

#alus<-rbindlist(svdf_list)
alus<-all_alus
alus$subFamily<-substr(alus$repName,1,4)
#make FLAM/FRAM one section 
alus[subFamily == 'FLAM', subFamily := 'FLAM/FRAM']
alus[subFamily == 'FRAM', subFamily := 'FLAM/FRAM']

#mean(alu_del$mse_mean, na.rm=T); sd(alu_del$mse_mean, na.rm=T)


#make colors for repName and subFamily
unique_repName<-unique(alus$repName)
repName<-unique_repName[order(unique_repName)]
repName_color<-c(stepped_colors[22], brewer.ylorbr(6)[2:5], brewer.greens(20)[4:20], brewer.pubu(32)[8:32], brewer.rdpu(8)[5:8])
repName_color_df<-data.table(repName, repName_color)


unique_subFamily<-unique(alus$subFamily)
subFamily<-c('FAM', 'FLAM/FRAM', 'Alu', 'AluJ', 'AluS', 'AluY')
subFamily_color<-c(brewer.rdpu(8)[5], brewer.rdpu(8)[7], stepped_colors[22], brewer.ylorbr(6)[4], brewer.greens(20)[14], brewer.pubu(32)[22])
subFamily_color_df<-data.table(subFamily, subFamily_color)
subFamily_level<-subFamily


levels1<-c(repName[48:51], repName[1:47])
alus_r<-merge(alus, repName_color_df, by='repName')
alus_plot<-merge(alus_r, subFamily_color_df, by='subFamily')

alus_plot$repName<-factor(alus_plot$repName, levels=levels1)
alus_plot$subFamily<-factor(alus_plot$subFamily, levels=subFamily_level)


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


pdf(paste0(fig_dir, '20250318_Alus_hg38all_MSE_CORR_distribution_noblacklist_w70.pdf'), width = 75/25.4, height = 37/25.4)
#pdf(paste0(fig_dir, 'Alus_akita_repName.pdf'), width=6, height=7)
outliersize=0.3
linewidth_mean<-0.4

all_alus_nona <- all_alus[!is.na(mse_mean)]
#annotate for top 10% most variable alu regions 
q99<-quantile(all_alus_nona$mse_mean, 0.99, na.rm=T)
q90<-quantile(all_alus_nona$mse_mean, 0.90, na.rm=T)
q50<-quantile(all_alus_nona$mse_mean, 0.50, na.rm=T)

# all Alu families 
p<-ggplot(all_alus_nona, aes(x=mse_mean)) + 
  geom_histogram(bins=60, color='navy', fill='#f2f6f8', linewidth=0.18) + 
  geom_vline(aes(xintercept = q50, color = "q50"),linetype = "dashed", linewidth=linewidth_mean) +
  geom_vline(aes(xintercept = q90, color = "q90"),linetype = "dashed", linewidth=linewidth_mean) +
  geom_vline(aes(xintercept = q99, color = "q99"),linetype = "dashed", linewidth=linewidth_mean) +
  scale_color_manual(name = "Statistics", values = c("q50" = "blue", "q90" = "orange", "q99"='red')) +
  xlab('MSE') + ylab('Frequency') + theme_classic() + theme_textsize
p

#repeat for CORR
all_alus_nona <- all_alus[!is.na(corr_score_mean)]
#annotate for top 10% most variable alu regions 
q99<-quantile(all_alus_nona$corr_score_mean, 0.99, na.rm=T)
q90<-quantile(all_alus_nona$corr_score_mean, 0.90, na.rm=T)
q50<-quantile(all_alus_nona$corr_score_mean, 0.50, na.rm=T)


p<-ggplot(all_alus_nona, aes(x=corr_score_mean)) + 
  geom_histogram(bins=60, color='navy', fill='#f2f6f8', linewidth=0.18) + 
  geom_vline(aes(xintercept = q50, color = "q50"),linetype = "dashed", linewidth=linewidth_mean) +
  geom_vline(aes(xintercept = q90, color = "q90"),linetype = "dashed", linewidth=linewidth_mean) +
  geom_vline(aes(xintercept = q99, color = "q99"),linetype = "dashed", linewidth=linewidth_mean) +
  scale_color_manual(name = "Statistics", values = c("q50" = "blue", "q90" = "orange", "q99"='red')) +
  xlab('CORR') + ylab('Frequency') + theme_classic() + theme_textsize
p


dev.off()



pdf(paste0(fig_dir, '20260318_MSE_CORR_percentiles.pdf'), width = 85/25.4, height = 55/25.4)

all_alus_nona <- all_alus_nona %>% mutate(MSE_decile = ntile(mse_mean, 10))
all_alus_nona <- all_alus_nona %>% mutate(CORR_decile = ntile(corr_score_mean, 10))

p1<-ggplot(all_alus_nona, aes(mse_mean, corr_score_mean)) +
  geom_hex(bins = 50) +
  scale_fill_gradientn(colours = c("#2B3B4F", '#88A0BF', "#ffffff"), trans = "log10", name = "log10(count)") +
  labs(fill = "log10(count)", x='log10(MSE)', y='log10(CORR)') +
  scale_x_log10() + scale_y_log10() + 
  theme_classic() + theme_textsize
p1

p2<-ggplot(all_alus_nona, aes(mse_mean, corr_score_mean)) +
  stat_density_2d(aes(color = after_stat(level)), geom = "contour", n=80, bins=14) +
  scale_color_gradientn(colours = c("#2B3B4F", '#88A0BF', "#f6f6f6"), trans = "log10", name = "log10(count)") +
  scale_x_log10() + scale_y_log10() + 
  labs( x='log10(MSE)', y='log10(CORR)') +
  theme_classic() + theme_textsize
p2

p2<-ggplot(all_alus_nona, aes(mse_mean, corr_score_mean)) +
  stat_density_2d(aes(color = after_stat(level)), geom = "contour", n=100, bins=12) +
  scale_color_gradientn(colours = c("#2B3B4F", '#88A0BF', "#f6f6f6"), trans = "log10", name = "log10(count)") +
  scale_x_log10() + scale_y_log10() + 
  labs( x='log10(MSE)', y='log10(CORR)') +
  theme_classic() + theme_textsize
p2
dev.off()


#divide by 25.4 so its in mm 
pdf(paste0(fig_dir, '20260318_Alus_hg38all_DEL_theme_classic.pdf'), width = 85/25.4, height = 100/25.4)
#pdf(paste0(fig_dir, 'Alus_akita_repName.pdf'), width=6, height=7)

outliersize=0.3
linewidth_box<-0.3
linewidth_mean<-0.4
# all Alu families 
p<-ggplot(alus_plot, aes(x=factor(repName), y=mse_mean, color=repName_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() + 
  #facet_wrap(~SVTYPE, nrow=1, scales='free_x') + 
  scale_color_identity()+
  geom_hline(data = alus_plot %>% group_by(SVTYPE) %>% summarize(median_mse = median(mse_mean, na.rm = TRUE)),
             aes(yintercept = median_mse), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  xlab('Alu Name') + ylab('MSE') + theme_classic() + theme_textsize
p

p1<-ggplot(alus_plot, aes(x=factor(repName), y=corr_score_mean, color=repName_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() + 
  #facet_wrap(~SVTYPE, nrow=1, scales='free_x') + 
  scale_color_identity()+
  geom_hline(data = alus_plot %>% group_by(SVTYPE) %>% summarize(median_corr = median(corr_score_mean, na.rm = TRUE)),
             aes(yintercept = median_corr), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  xlab('Alu Name') + ylab('CORR') + theme_classic() + theme_textsize
p1

dev.off()


#PLOT ON ITS SIDE 
pdf(paste0(fig_dir, '20260803_Alus_hg38all_DEL_theme_classic_SIDEWAYS_width180height90.pdf'), width = 180/25.4, height = 90/25.4)
#pdf(paste0(fig_dir, 'Alus_akita_repName.pdf'), width=6, height=7)

outliersize=0.3
linewidth_box<-0.3
linewidth_mean<-0.4
# all Alu families 
p<-ggplot(alus_plot, aes(x=factor(repName), y=mse_mean, color=repName_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) +
  #facet_wrap(~SVTYPE, nrow=1, scales='free_x') + 
  scale_color_identity()+
  geom_hline(data = alus_plot %>% group_by(SVTYPE) %>% summarize(median_mse = median(mse_mean, na.rm = TRUE)),
             aes(yintercept = median_mse), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  xlab('Alu Name') + ylab('MSE') + theme_classic() + theme_textsize + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

p1<-ggplot(alus_plot, aes(x=factor(repName), y=corr_score_mean, color=repName_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + 
  #facet_wrap(~SVTYPE, nrow=1, scales='free_x') + 
  scale_color_identity()+
  geom_hline(data = alus_plot %>% group_by(SVTYPE) %>% summarize(median_corr = median(corr_score_mean, na.rm = TRUE)),
             aes(yintercept = median_corr), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  xlab('Alu Name') + ylab('CORR') + theme_classic() + theme_textsize + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1

dev.off()


#CORR
p<-ggplot(alus_plot, aes(x=factor(repName), y=corr_score_mean, color=repName_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() + 
  facet_wrap(~SVTYPE, nrow=1, scales='free_x') + 
  scale_color_identity()+
  geom_hline(data = alus_plot %>% group_by(SVTYPE) %>% summarize(median_corr = median(corr_score_mean, na.rm = TRUE)),
             aes(yintercept = median_corr), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  xlab('Alu Name') + ylab('Mean 1-CORR (Augmented Scores)') + 
  theme_classic() + theme_textsize
p


#plot2
alus_del_only<-alus_plot[SVTYPE=='DEL']
p<-ggplot(alus_del_only, aes(x=factor(subFamily), y=corr_score_mean, color=subFamily_color)) + 
  geom_boxplot(outlier.size = outliersize, linewidth=linewidth_box) + coord_flip() +
  facet_wrap(~SVTYPE, nrow=1, scales='free_x') + scale_color_identity() +
  geom_hline(data = alus_del_only %>% group_by(SVTYPE) %>% summarize(median_corr = median(corr_score_mean, na.rm = TRUE)),
             aes(yintercept = median_corr), color = "red", linetype = "dashed", linewidth=linewidth_mean) +
  ylab('Mean  1-CORR (Augmented Scores)') + xlab('Alu Subfamily') + #labs(color='Alu Subfamily') +
  theme_classic() + theme_textsize
p
dev.off()


#GC VS WAHTEVER 


pdf(paste0(fig_dir, '20251205_Alus_hg38all_DEL_GC.pdf'), width=10, height=6)
#look at GC content vs MSE, split by alteration type and Alu Subfamily 
p<-ggplot(alus_del_only, aes(x=GC_alu, y=mse_mean, color=subFamily_color)) + 
  #geom_density(aes(group=subFamily)) + 
  geom_point(alpha=0.5) + #coord_flip() + 
  geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
  geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
  #facet_wrap(~subFamily, nrow=2, scales='free_y') + 
  facet_wrap(~subFamily, nrow=2) +  xlab('GC Alu')+
  ylab('Mean MSE (Augmented Scores)') + labs(color='Alu Subfamily') + 
  scale_color_identity() + theme_classic()
p

p<-ggplot(alus_del_only, aes(x=GC_1kb_up, y=mse_mean, color=subFamily_color)) + 
  #geom_density(aes(group=subFamily)) + 
  geom_point(alpha=0.5) + #coord_flip() + 
  geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
  geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
  #facet_wrap(~subFamily, nrow=2, scales='free_y') + 
  facet_wrap(~subFamily, nrow=2) + xlab('GC 1kb upstream')+
  ylab('Mean MSE (Augmented Scores)') + labs(color='Alu Subfamily') + 
  scale_color_identity() + theme_classic()
p

p<-ggplot(alus_del_only, aes(x=GC_1kb_down, y=mse_mean, color=subFamily_color)) + 
  #geom_density(aes(group=subFamily)) + 
  geom_point(alpha=0.5) + #coord_flip() + 
  geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
  geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
  #facet_wrap(~subFamily, nrow=2, scales='free_y') + 
  facet_wrap(~subFamily, nrow=2) +  xlab('GC 1kb downstream')+
  ylab('Mean MSE (Augmented Scores)') + labs(color='Alu Subfamily') + 
  scale_color_identity() + theme_classic()
p

p<-ggplot(alus_del_only, aes(x=GC_1Mb_window, y=mse_mean, color=subFamily_color)) + 
  #geom_density(aes(group=subFamily)) + 
  geom_point(alpha=0.5) + #coord_flip() + 
  geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
  geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
  #facet_wrap(~subFamily, nrow=2, scales='free_y') + 
  facet_wrap(~subFamily, nrow=2) + xlab('GC 1Mb window')+
  ylab('Mean MSE (Augmented Scores)') + labs(color='Alu Subfamily') + 
  scale_color_identity() + theme_classic()
p
dev.off()



# ##### GC FLANKING OVERLAP PLOTS
# 
# pdf(paste0(fig_dir, 'Alu_Akita_GC_flanking_fixedy.pdf'), width=10, height=2)
# 
# ggplot(alu_sv_df, aes(x=GC.x, y=GC_POSflank_100, color=subFamily_color)) + 
#   geom_point(alpha=0.5) + #coord_flip() + 
#   geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
#   geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
#   #facet_wrap(~subFamily, nrow=1, scales='free_y') + 
#   facet_wrap(~subFamily, nrow=1) + 
#   ylab('Start Flank GC (100bp)') + xlab('Alu GC') + 
#   scale_color_identity() + theme_classic()
# 
# 
# ggplot(alu_sv_df, aes(x=GC.x, y=GC_ENDflank_100, color=subFamily_color)) + 
#   geom_point(alpha=0.5) + #coord_flip() + 
#   geom_vline(xintercept=0.45,lwd=0.6, linetype='dashed', colour="red") + 
#   geom_vline(xintercept=0.55,lwd=0.6,linetype='dashed', colour="red") + 
#   #facet_wrap(~subFamily, nrow=1, scales='free_y') + 
#   facet_wrap(~subFamily, nrow=1) + 
#   ylab('End Flank GC (100bp)') + xlab('Alu GC') + 
#   scale_color_identity() + theme_classic()
# 
# dev.off()



##### look at CORR between mappability and score 
#### look at CORR between mappability and var 


#list all pairs that we want to plto 
map<-'map_1kb_up'
pairs <- list(
  c("mse_mean", map),
  c("corr_score_mean", map),
  c("mse_std", map),
  c("corr_std", map)
)

plots <- lapply(pairs, \(p) {
  ggplot(alus_plot, aes_string(p[1], p[2], color="subFamily_color")) +
    geom_point(size=0.2) +
    facet_wrap(~subFamily, nrow=2) +
    scale_x_log10()+
    scale_color_identity()  + theme_classic() +
    ggtitle(paste(p[1], "vs", p[2]))
})

pdf(paste0(fig_dir, '20251205_Alus_hg38all_DEL_mappalu1kbup_score_log10.pdf'), width=6, height=6)
plots  # prints each plot
dev.off()









# TOP BOTTOM MSE Alus  
mse_col<-alus_del_only$mse_mean
top_10_quantile <- quantile(mse_col, 0.9, na.rm=T)
bottom_10_quantile <- quantile(mse_col, 0.1, na.rm=T)

top_10<-mse_col>top_10_quantile
bottom_10<-mse_col<bottom_10_quantile

df_top<-alus_del_only[top_10]
df_bottom<-alus_del_only[bottom_10]

df_combined<-rbind(cbind(df_top, score='top'), cbind(df_bottom, score='bottom'))
head(df_combined)

test1<-data.table(table(alus_del_only$subFamily))
colnames(test1)<-c('subFamily', 'subFamily_count')

df_merged <- merge(df_combined, test1, by = "subFamily", all.x = TRUE)
df_merged$subFamily<-factor(df_merged$subFamily, levels=subFamily_level)


pdf(paste0(fig_dir, 'top_bottom_mse.pdf'), width=7, height=3)

p<-ggplot(df_combined, aes(x=subFamily, fill=score)) + 
  geom_bar(position="dodge", stat="count") +
  ylab('Count') + xlab('Alu Family') + 
  theme_classic()  #+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p

p<-ggplot(df_merged, aes(x=subFamily, fill=score)) + 
  #geom_bar(position="dodge", stat="count", aes(weight = subFamily_count / sum(subFamily_count))) + 
  geom_bar(position="dodge", stat="count", aes(weight = sum(subFamily_count)/subFamily_count)) + 
  
  ylab('Normalized Counts (by Samples/Family)') + xlab('Alu Family') + 
  theme_classic()  #+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p

dev.off()




####### make plots of proportion of each family in top 1% disruptive mse


#USE FOR GENERAL TEXT SIZE THINGS 
theme_textsize1<-theme(
  axis.title = element_text(size = 7),     # axis titles ~7 pt
  strip.text = element_text(size = 7),      # facet headers ~7 pt
  plot.title = element_text(size = 8),      # if you add a title
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.line  = element_line(linewidth = 0.3),  # x & y axis lines
  axis.ticks = element_line(linewidth = 0.3)
)



alus_plot_noNA <- alus_plot[!is.na(mse_mean)]
mse_1perc <- quantile(alus_plot_noNA$mse_mean, 0.99, na.rm = TRUE)
alus_plot_noNA[, is_top1 := mse_mean >= mse_1perc]

repName_summary <- alus_plot_noNA[, .(total_N = .N, top1_N  = sum(is_top1)),by = repName][, prop_top1 := top1_N / total_N]
repName_summary<-merge(repName_summary, repName_color_df, by='repName')
  
  
p_prop <- ggplot(repName_summary,aes(y = reorder(repName, prop_top1),x = prop_top1, fill=repName_color)) +
  geom_col(width = 0.7) +
  scale_x_continuous(labels = scales::percent) + scale_fill_identity() + 
  geom_vline(xintercept=0.01,lwd=0.6, linetype='dashed', colour="red") +
  labs(x = "Alus in Top 1% / Total Alu",y = NULL, title = "Proportion of Alu type in top 1% MSE") +
  theme_classic() + theme_textsize1

#add numbers 
p_prop <- p_prop + geom_text(aes(label = paste0(top1_N, "/", total_N)), hjust = -0.1, size = 2) +
  coord_cartesian(xlim = c(0, max(repName_summary$prop_top1) * 1.25)) + theme_textsize1

#add in total counts 
p_total <- ggplot(repName_summary, aes(y = reorder(repName, prop_top1), x = total_N, fill=repName_color)) +
  geom_col(width = 0.7) +
  scale_fill_identity() + #scale_x_log10()+
  labs(x = "Count", y = NULL, title = "Total Alu") + theme_classic() + theme_textsize1  


pdf(paste0(fig_dir, 'alu_nona_prop_top1.pdf'), height = 210/25.4, width = 180/25.4)
plot_grid(p_prop, p_total, nrow = 1, align = "h", rel_widths = c(1.2, 0.4))
dev.off()


#look at top 1%
alu_top1<-alus_plot_noNA[alus_plot_noNA$is_top1]
repname1<-alu_top1[alu_top1$repName=='AluYb9']

