# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

library(data.table)
library(ggplot2)
require(RIdeogram)
#library(rsvg)
library(dplyr)
library(GenomicRanges)
library(ggrepel)


res_dir<-paste0(PROJECT_DIR, "/results/paper_results/")
#fig_dir<-paste0(PROJECT_DIR, "/figs/bmi_rips_remake/")
#fig_dir<-paste0(PROJECT_DIR, "/figs/20251205_alu_all/")
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


#data 
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
chrom_lengths_hg38<-fread(paste0(PROJECT_DIR, "/chrom_lengths_hg38"))
colnames(chrom_lengths_hg38)<-c('chrom', 'chrom_len')

all_alus<-fread(paste0(res_dir, '20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'))
all_alus$chrom_noprefix <- sub("^chr", "", all_alus$CHROM)

#remove na scores 
all_alus <- all_alus[!is.na(mse_mean)]


#annotate for top 10% most variable alu regions 
#top_alu_std <- all_alus[ mse_std >= quantile(mse_std, 0.99, na.rm=T)]
top_alu_mse <- all_alus[ mse_mean >= quantile(mse_mean, 0.99, na.rm=T)]



##################################################################

# convert alus ranges to granges 
#alu_highvar_grange <- GRanges(seqnames = top_alu_std$chrom_noprefix,
#                      ranges   = IRanges(top_alu_std$POS, top_alu_std$END))
alu_highmse_grange <- GRanges(seqnames = top_alu_mse$chrom_noprefix,
                              ranges   = IRanges(top_alu_mse$POS, top_alu_mse$END))
alu_all_grange <- GRanges(seqnames = all_alus$chrom_noprefix,
                              ranges   = IRanges(all_alus$POS, all_alus$END))


#make 1Mb window ranges 
window_len<-1e6
bin_list <- lapply(1:nrow(chrom_lengths_hg38), function(i) {
  chr <- paste0(chrom_lengths_hg38$chrom[i])
  len <- chrom_lengths_hg38$chrom_len[i]
  starts <- seq(1, len, by = window_len)
  ends   <- pmin(starts + window_len - 1, len)
  GRanges(seqnames = chr, ranges = IRanges(starts, ends))
})

gr_bins <- Reduce(c, bin_list)

#count overlaps between Alus and genomic bins 
#ov_counts_highvar <- countOverlaps(gr_bins, alu_highvar_grange)
ov_counts_highmse <- countOverlaps(gr_bins, alu_highmse_grange)
ov_counts_allalu <- countOverlaps(gr_bins, alu_all_grange)

#make into dt 

alu_all_windows <- data.table(Chr = as.character(seqnames(gr_bins)), Start = start(gr_bins), 
                              End = end(gr_bins), Value = ov_counts_allalu)
#alu_highvar_windows$Value<-as.numeric(alu_highvar_windows$Value)
alu_all_windows$Value <- alu_all_windows$Value 
alu_all_windows$Color<-'green'

############ high SD alus 
# alu_highvar_windows <- data.table(Chr = as.character(seqnames(gr_bins)), Start = start(gr_bins), 
#                  End = end(gr_bins), Value = ov_counts_highvar)
# #alu_highvar_windows$Value<-as.numeric(alu_highvar_windows$Value)
# alu_highvar_windows$Value <- alu_highvar_windows$Value * 1e7
# 
# 
# #normalize high var windows by number of alus in the window
# alu_highvar_windows<-merge(alu_highvar_windows, alu_all_windows, by=c('Chr', 'Start', 'End'))
# alu_highvar_windows$Value<-((alu_highvar_windows$Value.x)/(alu_highvar_windows$Value.y))*1e8
# alu_highvar_windows$Color<-'blue'
# 
# alu_highvar_windows<-alu_highvar_windows[,c('Chr', 'Start', 'End', 'Value', 'Color')]


############ high MSE alus 
alu_highmse_windows <- data.table(Chr = as.character(seqnames(gr_bins)), Start = start(gr_bins), 
                                  End = end(gr_bins), Value = ov_counts_highmse)
#alu_highvar_windows$Value<-as.numeric(alu_highvar_windows$Value)
alu_highmse_windows$Value <- alu_highmse_windows$Value 


#normalize high mse windows by number of alus in the window
alu_highmse_windows<-merge(alu_highmse_windows, alu_all_windows, by=c('Chr', 'Start', 'End'))
alu_highmse_windows$Value<-log10(((alu_highmse_windows$Value.x)/(alu_highmse_windows$Value.y)))
alu_highmse_windows$Value_unlog<-(((alu_highmse_windows$Value.x)/(alu_highmse_windows$Value.y)))

alu_highmse_windows$Color<-'red'


pdf(paste0(PROJECT_DIR, "/figs/paper_figs/1mb_disruption_topwindows_shortlabel.pdf"), width = 70/25.4, height = 65/25.4)
# Identify outliers - using a simple approach based on distance from the diagonal
# or based on residuals. Adjust the threshold as needed.
#alu_highmse_windows <- alu_highmse_windows %>%mutate(coord_label = paste0(Chr, ":", Start, "-", End),
alu_highmse_windows <- alu_highmse_windows %>%mutate(coord_label = paste0('chr', Chr, ":", Start-1),

    # Flag outliers: e.g., points far from the x=y line, or extreme values
    #is_outlier = abs(-(Value_unlog - Value.y)) > quantile(abs(-(Value_unlog - Value.y)), 0.99, na.rm = TRUE))
    is_outlier = abs(Value_unlog) > quantile(abs(Value_unlog), 0.995, na.rm = TRUE))
    

ggplot(alu_highmse_windows, aes(x = Value.y, y = Value_unlog)) +
  geom_point(aes(color = 'darkgreen'), alpha = 0.6, size=0.3) +
  geom_text_repel(data = filter(alu_highmse_windows, is_outlier), aes(label = coord_label),size = 1,
  max.overlaps = Inf, box.padding = 0.5) +
  scale_color_identity() +  # uses the Color column values directly (e.g., "red")
  #scale_x_log10()+
  labs(x = "Total # Alu per bin",y = "Proportion of Alus that are highly disruptive") +
  theme_classic() + theme_textsize

dev.off()


alu_highmse_windows<-alu_highmse_windows[,c('Chr', 'Start', 'End', 'Value', 'Color')]

#check histogram of values
hist(alu_highmse_windows$Value)
#replace Infs and Nans with slightly lower value than lowest value
val_col <- alu_highmse_windows$Value
val_col[!is.finite(val_col)] <- min(val_col[is.finite(val_col)]) - 4
alu_highmse_windows$Value <- val_col



#look at correlation between alu variability and other features 
# alu_highvar_windows<-alu_highvar_windows[,c('Chr', 'Start', 'End', 'Value', 'Color')]
# alu_highmse_windows<-alu_highmse_windows[,c('Chr', 'Start', 'End', 'Value', 'Color')]
#cor(alu_highvar_windows$mse_mean, all_alus$corr_score_mean, use='complete.obs', method='spearman')


##################################################################
#plot!
#svg_data <- readBin("chromosome.svg", what = "raw", n = file.info("chromosome.svg")$size)

# ideogram(karyotype = human_karyotype, overlaid=gene_density,width=300)
#          #label = test_label, label_type = "line")
# convertSVG(svg_data, device = 'png')
#read_data()


# ideogram(karyotype = human_karyotype, overlaid=gene_density, width=250,
#          label = alu_highvar_windows, label_type = "polygon")
# convertSVG("chromosome.svg", device = "pdf")

ideogram(karyotype = human_karyotype, overlaid=alu_highmse_windows, width=250,
         #label = alu_highmse_windows, label_type = "polygon", 
         colorset1 = c("#F5F5F5", '#FFFFFF', '#F2F8F5', '#E4F1EB', '#D0E7DC', '#BBDDCC', '#85C1A5', '#69B590', '#356E53', "#28533F"))
# #FFFFFF #'#356E53',"#1B372A"
convertSVG("chromosome.svg", device = "pdf")
#svg2pdf("chromosome.svg")

