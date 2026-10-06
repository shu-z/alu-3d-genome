# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shuzhang/pollard_lab/alu")

library(data.table)
library(ggplot2)



#old results folders are 20240118 and 20240129
res_dir<-paste0(PROJECT_DIR, "/results/20250514_GC_mut/")
list.files(res_dir)

#fig_dir<-paste0(PROJECT_DIR, "/figs/20240201/")
#fig_dir<-paste0(PROJECT_DIR, "/figs/20250514_RERUN/")
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
outliersize=0.3

svdf_list<-list()
#for (i in c('DEL', 'DUP', 'INV')){
for (i in c('DEL')){
  
  #read in original alu df and DEL score 
  df_sv<-fread(paste0(res_dir, '20250513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt'), sep='\t')
  setnames(df_sv, "mse_mean", paste0("0_mse_0"))
  
  # Iterate through the dataframes
  for (i in c(10, 20, 30, 40, 50)) {
    # Read dataframe as data.table
    dt <- fread(paste0(res_dir, '20250513_hg38Alus_100sample_repName_DEL_GCmut', i, '_scores'))
    #print(dt)
    
    # Add prefix to "gc" column
    #dt[, mse_mean := paste0("gc", i, "_", mse_mean)]
    setnames(dt, "mse_0", paste0(i, "_mse_0"))
    #print(dt)
    
    # Select relevant columns
    dt <- dt[, 1:2]
    
    # Bind rows to the combined data.table
    df_sv <- merge(df_sv, dt, by='var_index')
  }


  
  
  

  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
   # geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  
  #############plots

  #ALL MSE HISTOGRAMS
  all_mse_df<-df_sv[, c('disruption','0_mse_0', '10_mse_0', '20_mse_0', '30_mse_0', '40_mse_0', '50_mse_0')]
  melted_dt <- melt(all_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")
  
  # Extract the numeric part from the variable names
  melted_dt[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
  melted_dt$gc_mutate<-as.factor(melted_dt$gc_mutate)
  
  
  
  pdf(paste0(fig_dir,'GC_levels_decrease_10_50_MSE_logscale_boxplot.pdf'),  width=67/25.4, height=50/25.4)
    #histogram of original scores, gc (top/bottom colored, facet wrapped old/new)
  p<-ggplot(melted_dt, aes(y = mse_mean, x = gc_mutate, fill=as.factor(disruption))) +
    geom_boxplot(outlier.size=outliersize,alpha = 0.80) +
    scale_y_continuous(trans='log10') + 
    scale_fill_manual(values = c("bottom" = "#D1D2F0", "top" = "#5C5EAD"),  name = "Disruption") + 
    coord_cartesian(ylim = c(1e-9, 1e-2)) + 
    labs(x = "%GC Decrease", y = "log10(MSE)") +
    theme_classic() + theme_textsize
  print(p)
  
  
  mut_mse_df<-df_sv[, c('disruption', '10_mse_0', '20_mse_0', '30_mse_0', '40_mse_0', '50_mse_0')]
  melted_mut_dt <- melt(mut_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")
  p<-ggplot(melted_mut_dt, aes(y = mse_mean, x = gc_mutate, fill=as.factor(disruption))) +
    geom_boxplot(outlier.size=outliersize,alpha = 0.80) +
    scale_y_continuous(trans='log10') + 
    scale_fill_manual(values = c("bottom" = "#D1D2F0", "top" = "#5C5EAD"),  name = "Disruption") + 
    coord_cartesian(ylim = c(1e-9, 1e-2)) + 
    labs(x = "%GC Decrease", y = "log10(MSE)") +
    theme_classic() + theme_textsize
  print(p)
  dev.off()
# 
#  #GC HISTOGRAMS
#  #now add difference metric
#  melted_dt1 <- melt(all_mse_df, id.vars = c("disruption", "0_mse_0"), variable.name = "gc_mutate", value.name = "mse_mean")
#  melted_dt1[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
#  melted_dt1[, gc_mse_diff:=`0_mse_0`-mse_mean]
#  melted_dt1$gc_mutate<-as.factor(melted_dt1$gc_mutate)
#  
#  
#  #histogram of change in scores, gc (top vs bottom)
#  p<-ggplot(melted_dt1, aes(y = gc_mse_diff, x=gc_mutate, fill=disruption)) +
#    geom_boxplot( alpha = 0.75) +
#    scale_y_continuous(trans='log10') + 
#    labs(title = paste0('(Alu MSE) - (Alu GC Loss MSE)'),
#         x = "%GC Loss", y = "MSE Diff") +
#    theme_classic()
#  print(p)
 






}



#################### CHECK REVCOMP

library(data.table)
library(ggplot2)



#old results folders are 20240118 and 20240129
res_dir<-paste0(PROJECT_DIR, "/results/20250524_GC_checkREVCOMP/")
list.files(res_dir)

#fig_dir<-paste0(PROJECT_DIR, "/figs/20240201/")
fig_dir<-paste0(PROJECT_DIR, "/figs/20250514_GCmut_RERUN/")

#read in original alu df and DEL score 
df_sv<-fread(paste0(res_dir, '20250513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt'), sep='\t')
setnames(df_sv, "mse_mean", paste0("0_mse_0"))
df_sv_revcomp<-df_sv


# Iterate through the dataframes
for (i in c(10, 25, 50)) {
  # Read dataframe as data.table
  df_REVCOMP <- fread(paste0(res_dir, '20250524_DEL__MSE_topbottom_GCmut', i, '_shift0_REVCOMP_scores'))
  setnames(df_REVCOMP, "mse_0_revcomp", paste0(i, "_mse_0"))
  df_REVCOMP <- df_REVCOMP[, 1:2]
  
  # Bind rows to the combined data.table
  df_sv_revcomp <- merge(df_sv_revcomp, df_REVCOMP, by.x='V1', by.y='var_index')
  
  # Read dataframe as data.table
  df_noRC <- fread(paste0(res_dir, '20250524_DEL__MSE_topbottom_GCmut', i, '_shift0_scores'))
  setnames(df_noRC, "mse_0", paste0(i, "_mse_0"))
  df_noRC <- df_noRC[, 1:2]
  
  # Bind rows to the combined data.table
  df_sv <- merge(df_sv, df_noRC,  by.x='V1', by.y='var_index')
}



  
  
  

  #############plots
  
  pdf(paste0(res_dir,'GC_levels_10_50_MSE.pdf'), width=8, height=6)
  
  #ALL MSE HISTOGRAMS
  all_mse_df<-df_sv[, c('disruption', '10_mse_0', '25_mse_0', '50_mse_0', 'strand')]
  melted_dt <- melt(all_mse_df, id.vars = c("disruption", "strand"), variable.name = "gc_mutate", value.name = "mse_mean")
  
  # Extract the numeric part from the variable names
  melted_dt[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
  melted_dt$gc_mutate<-as.factor(melted_dt$gc_mutate)

  #histogram of original scores, gc (top/bottom colored, facet wrapped old/new)
  p<-ggplot(melted_dt, aes(y = mse_mean, x = gc_mutate, fill=disruption)) +
    geom_boxplot( alpha = 0.75) +
    facet_wrap(~strand) + 
    scale_y_continuous(trans='log10') + 
    labs(title = paste0('MSE scores with Alu GC Loss'),
         x = "%GC Loss", y = "MSE") +
    theme_classic()
  print(p)
  
  
  
  #ALL MSE HISTOGRAMS
  all_mse_df_RC<-df_sv_revcomp[, c('disruption', '10_mse_0', '25_mse_0', '50_mse_0', 'strand')]
  melted_dt_RC <- melt(all_mse_df_RC, id.vars = c("disruption", 'strand'), variable.name = "gc_mutate", value.name = "mse_mean")
  
  # Extract the numeric part from the variable names
  melted_dt_RC[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
  melted_dt_RC$gc_mutate<-as.factor(melted_dt_RC$gc_mutate)
  
  #histogram of original scores, gc (top/bottom colored, facet wrapped old/new)
  p<-ggplot(melted_dt_RC, aes(y = mse_mean, x = gc_mutate, fill=disruption)) +
    geom_boxplot( alpha = 0.75) +
    facet_wrap(~strand) + 
    scale_y_continuous(trans='log10') + 
    labs(title = paste0('MSE scores with Alu GC Loss (REVCOMP)'),
         x = "%GC Loss", y = "MSE") +
    theme_classic()
  dev.off()
  
  
  
  
  
  
#################### CHECK AT MUTATE (INCREASE GC)
  
res_dir<-paste0(PROJECT_DIR, "/results/20250623_AT_mutate/")
fig_dir<-paste0(PROJECT_DIR, "/figs/20250623_AT_mutate/")

#read in original alu df and DEL score 
df_sv<-fread(paste0(res_dir, '20250513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt'), sep='\t')
setnames(df_sv, "mse_mean", paste0("0_mse_0"))


# Iterate through the dataframes
for (i in c(10, 20, 30, 40, 50)) {
  # Read dataframe as data.table
  df_mut <- fread(paste0(res_dir, '20250623_hg38Alus_100sample_repName_DEL_ATmut', i, '_scores'))
  setnames(df_mut, "mse_0", paste0(i, "_mse_0"))
  df_mut <- df_mut[, 1:2]
  
  # Bind rows to the combined data.table
  df_sv <- merge(df_sv, df_mut, by.x='V1', by.y='var_index')

}


############# plot


#ALL MSE HISTOGRAMS
all_mse_df<-df_sv[, c('disruption','0_mse_0', '10_mse_0', '20_mse_0', '30_mse_0', '40_mse_0', '50_mse_0')]
melted_dt <- melt(all_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")

# Extract the numeric part from the variable names
melted_dt[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
melted_dt$gc_mutate<-as.factor(melted_dt$gc_mutate)

pdf(paste0(fig_dir,'GC_levels_increase_10_50_MSE.pdf'), width=67/25.4, height=50/25.4)
#histogram of original scores, gc (top/bottom colored, facet wrapped old/new)
p<-ggplot(melted_dt, aes(y = mse_mean, x = gc_mutate, fill=as.factor(disruption))) +
  geom_boxplot(outlier.size=outliersize,alpha = 0.80) +
  scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c("bottom" = "#EFD2D6", "top" = "#BC4E60"),  name = "Disruption") + 
  coord_cartesian(ylim = c(1e-9, 1e-2)) + 
  labs(x = "%GC Increase", y = "log10(MSE)") +
  theme_classic() + theme_textsize
print(p)


mut_mse_df<-df_sv[, c('disruption', '10_mse_0', '20_mse_0', '30_mse_0', '40_mse_0', '50_mse_0')]
melted_mut_dt <- melt(mut_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")
p<-ggplot(melted_mut_dt, aes(y = mse_mean, x = gc_mutate, fill=as.factor(disruption))) +
  geom_boxplot(outlier.size=outliersize,alpha = 0.80) +
  scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c("bottom" = "#EFD2D6", "top" = "#BC4E60"),  name = "Disruption") + 
  coord_cartesian(ylim = c(1e-9, 1e-2)) + 
  labs(x = "%GC Increase", y = "log10(MSE)") +
  theme_classic() + theme_textsize
print(p)

dev.off()


#################################### CpG mutate 

#old results folders are 20240118 and 20240129
res_dir<-paste0(PROJECT_DIR, "/results/20250710_CpG_mutate/")
list.files(res_dir)

#fig_dir<-paste0(PROJECT_DIR, "/figs/20240201/")
fig_dir<-paste0(PROJECT_DIR, "/figs/20250710_CpG_mutate/")

svdf_list<-list()
#for (i in c('DEL', 'DUP', 'INV')){
for (i in c('DEL')){
  
  #read in original alu df and DEL score 
  df_sv<-fread(paste0(res_dir, '20250513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt'), sep='\t')
  setnames(df_sv, "mse_mean", paste0("0_mse_0"))
  
  # Iterate through the dataframes
  for (i in c(50, 100, 200)) {
    # Read dataframe as data.table
    dt <- fread(paste0(res_dir, '20250710_hg38Alus_100sample_repName_DEL_MSE_topbottom_CpGmut_', i, '_scores'))
    #print(dt)
    
    # Add prefix to "gc" column
    #dt[, mse_mean := paste0("gc", i, "_", mse_mean)]
    setnames(dt, "mse_0", paste0(i, "_mse_0"))
    #print(dt)
    
    # Select relevant columns
    dt <- dt[, 1:2]
    
    # Bind rows to the combined data.table
    df_sv <- merge(df_sv, dt, by='var_index')
  }
  
  
  
  
  
  
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  # geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  
  #############plots
  
  #ALL MSE HISTOGRAMS
  all_mse_df<-df_sv[, c('disruption','0_mse_0',  '50_mse_0', '100_mse_0', '200_mse_0')]
  melted_dt <- melt(all_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")
  
  # Extract the numeric part from the variable names
  melted_dt[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
  melted_dt$gc_mutate<-as.factor(melted_dt$gc_mutate)
  
  
  
  pdf(paste0(fig_dir,'CpG increase_50_100_200_MSE_logscale_boxplot.pdf'), width=8, height=6)
  #histogram of original scores, gc (top/bottom colored, facet wrapped old/new)
  p<-ggplot(melted_dt, aes(y = mse_mean, x = gc_mutate, fill=disruption)) +
    geom_boxplot( alpha = 0.75) +
    scale_y_continuous(trans='log10') + 
    coord_cartesian(ylim = c(1e-9, 1e-2)) + 
    labs(title = paste0('MSE scores with Alu CpG Increase'),
         x = "%CpG Increase", y = "log10(MSE)") +
    theme_classic()
  print(p)
  
  
  mut_mse_df<-df_sv[, c('disruption', '50_mse_0', '100_mse_0', '200_mse_0')]
  melted_mut_dt <- melt(mut_mse_df, id.vars = c("disruption"), variable.name = "gc_mutate", value.name = "mse_mean")
  p<-ggplot(melted_mut_dt, aes(y = mse_mean, x = gc_mutate, fill=disruption)) +
    geom_boxplot( alpha = 0.75) +
    scale_y_continuous(trans='log10') + 
    coord_cartesian(ylim = c(1e-9, 1e-2)) + 
    labs(title = paste0('MSE scores with Alu CpG Increase'),
         x = "%CpG Increase", y = "log10(MSE)") +
    theme_classic()
  print(p)
  # 
  #  #GC HISTOGRAMS
  #  #now add difference metric
  #  melted_dt1 <- melt(all_mse_df, id.vars = c("disruption", "0_mse_0"), variable.name = "gc_mutate", value.name = "mse_mean")
  #  melted_dt1[, gc_mutate := as.integer(gsub("_.*", "", gc_mutate))]
  #  melted_dt1[, gc_mse_diff:=`0_mse_0`-mse_mean]
  #  melted_dt1$gc_mutate<-as.factor(melted_dt1$gc_mutate)
  #  
  #  
  #  #histogram of change in scores, gc (top vs bottom)
  #  p<-ggplot(melted_dt1, aes(y = gc_mse_diff, x=gc_mutate, fill=disruption)) +
  #    geom_boxplot( alpha = 0.75) +
  #    scale_y_continuous(trans='log10') + 
  #    labs(title = paste0('(Alu MSE) - (Alu GC Loss MSE)'),
  #         x = "%GC Loss", y = "MSE Diff") +
  #    theme_classic()
  #  print(p)
  
  
  
  
  
  dev.off()
  
}


