# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

library(ggplot2)
library(data.table)

res_dir<-paste0(PROJECT_DIR, "/results/")


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


svs<-fread(paste0(res_dir, '20260212_1KG/1KG_T2T_scored_SV.txt'))
svs[, MAF:=ifelse(AF>0.5, 1-AF, AF)]
svs[, MAF_quantile=]
pdf(paste0(PROJECT_DIR, "/figs/paper_figs/1kg_AF.pdf"), width = 90/25.4, height = 75/25.4)

col_val<-c('#AFCB90', '#688B41')
p_r2 <- ggplot(svs, aes(x = window, y = r2, fill = which_sv)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = col_val) +
  labs(x = "Window", y = "R²", fill = "Interaction degree") +
  theme_classic() + theme_textsize


cowplot::plot_grid(p_r2, p_mse, ncol=1)

dev.off()