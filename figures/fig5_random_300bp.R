# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)


#theme things
theme_textsize_border<-theme(
    axis.title = element_text(size = 7),     # axis titles ~7 pt
    strip.text = element_text(size = 7),      # facet headers ~7 pt
    plot.title = element_text(size = 8),      # if you add a title
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.line  = element_line(linewidth = 0.3),  # x & y axis lines
    axis.ticks = element_line(linewidth = 0.3),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
  )


alu_random_100kb <- fread(paste0(PROJECT_DIR, "/results/paper_results/alu_300bp_100kb_random_scores_combined.txt"))
alu_random_10kb  <- fread(paste0(PROJECT_DIR, "/results/paper_results/alu_300bp_10kb_random_scores_combined.txt"))

#add window label and combine 
alu_random_100kb[, window := "100kb"]
alu_random_10kb[,  window := "10kb"]
alu_random <- rbindlist(list(alu_random_100kb, alu_random_10kb), fill = TRUE)

#add labels and effect size 
#alu_random[, disruption_category:=ifelse(mse_mean_aluDEL>0.0286, '1% MSE', 'bottom 50% MSE')]
alu_random[, disruption_category:=disruption]
alu_random[, in_gene:=ifelse(shuffle_gene_overlap, 'Gene Overlap', 'No Gene Overlap')]
alu_random[, effect_size := mse_mean - mse_mean_aluDEL]

df_long <- alu_random %>%
  mutate(shuffled_label = paste("Shuffle", window)) %>%
  pivot_longer(cols = c(mse_mean_aluDEL, mse_mean),
    names_to = "score_type", values_to = "score") %>%
  mutate(score_type = case_when(score_type == "mse_mean_aluDEL" ~ "aluDEL",score_type == "mse_mean" ~ shuffled_label
    ))

df_effect <- alu_random %>% mutate(effect_size = mse_mean_aluDEL-mse_mean) 

df_effect_avg <- df_effect %>% 
  group_by(original_idx, window) %>% 
  summarize( avg_effect_size = mean(effect_size, na.rm = TRUE), 
             in_gene = dplyr::first(in_gene), 
             disruption = dplyr::first(disruption_category), 
             mse_original = dplyr::first(mse_mean_aluDEL), .groups = "drop")



pdf(paste0(PROJECT_DIR, "/figs/paper_figs/alu300bp_window_final_2.pdf"), width=110/25.4, height=60/25.4)
lwd<-0.3
ggplot(df_long, aes(x = score_type, y = score)) +
  geom_violin(trim = FALSE,  position=position_dodge(width=0.8), linewidth=lwd) +
  geom_boxplot(width = 0.2, outlier.size=0.3, position=position_dodge(width=0.8), linewidth=lwd) +
  facet_grid(disruption_category ~ in_gene) +
  scale_y_log10() + 
  labs(x = "Disruption Group", y = "Original aluDEL MSE") + 
  theme_classic() + theme_textsize_border 


col_val_window<-c('#C5E4FB', '#bc4e60')
col_val_highlow <- c(
  'low_10kb'   = '#CEE0ED',  # light blue
  'low_100kb'   = '#6AA1C8',  # darker blue
  'high_10kb'  = '#E3B5BD',  # lighter red/pink
  'high_100kb'  = '#C15C6D'   # original red (now the "dark" one)
)

df_effect_avg <- df_effect_avg %>%
  mutate(fill_group = paste(disruption, ifelse(window == '100kb', '100kb', '10kb'), sep = '_'),
    fill_group = factor(fill_group, levels = c('low_10kb', 'low_100kb', 'high_10kb', 'high_100kb')),
    disruption = factor(disruption, levels = c('low', 'high'))
    
  )
df_summary <- df_effect_avg %>%
  group_by(disruption, in_gene, window,fill_group) %>%
  summarize(frac_positive = mean(avg_effect_size > 0, na.rm = TRUE),
    median_effect = median(avg_effect_size, na.rm = TRUE),
    n = n(), .groups = "drop"
  )

ggplot(df_effect_avg, aes(x = disruption, y = avg_effect_size, fill = fill_group, group = interaction(disruption, window))) +
  geom_violin(trim = TRUE, width = 1, position = position_dodge(width = 0.8), linewidth = lwd, alpha = 0.8) +
  geom_boxplot(width = 0.2, outlier.size = 0.3, position = position_dodge(width = 0.8), linewidth = lwd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  scale_fill_manual(values = col_val_highlow) +
  facet_grid(~ in_gene) +
  labs(x = "Disruption Category", y = "MSE_shuffle - MSE_AluDEL") +
  theme_classic() + theme_textsize_border

#this one has text labels 
ggplot(df_effect_avg, aes(x = disruption, y = avg_effect_size, fill = fill_group, group = interaction(disruption, window))) +
  geom_violin(trim = TRUE, width = 1, position = position_dodge(width = 0.8), linewidth = lwd, alpha = 0.8) +
  geom_boxplot(width = 0.2, outlier.size = 0.3, position = position_dodge(width = 0.8), linewidth = lwd) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  geom_text(data = df_summary,aes(y = Inf, label = sprintf("%.0f%%", 100 * frac_positive)),
            position = position_dodge(width = 0.8),vjust = 1.5, size = 2) +
  scale_fill_manual(values = col_val_highlow) +
  facet_grid(~ in_gene) +
  labs(x = "Disruption Category", y = expression( MSE[Alu_DEL]-MSE[random_DEL]),fill = "Window") +
  theme_classic() + theme_textsize_border


dev.off()


######## ######## get stats 
library(dplyr)

grp_cmp <- function(dat, value = "avg_effect_size", group = "disruption",
                    high, low, alt = "two.sided") {
  dat <- dat[!is.na(dat[[value]]), ]
  x <- dat[[value]][dat[[group]] == high]   # high disruption
  y <- dat[[value]][dat[[group]] == low]    # low disruption
  n1 <- length(x); n2 <- length(y)
  tt <- wilcox.test(x, y, alternative = alt, exact = FALSE, conf.int = TRUE)  # Mann-Whitney
  U  <- unname(tt$statistic)
  A  <- U / (n1 * n2)                         # Vargha-Delaney: P(high > low)
  tibble::tibble(
    n_high       = n1,             n_low      = n2,
    median_high  = median(x),      median_low = median(y),
    HL_shift     = unname(tt$estimate),       # location diff (high - low), MSE units
    ci_low       = tt$conf.int[1], ci_high    = tt$conf.int[2],
    W            = U,
    p_value      = tt$p.value,
    cliffs_delta = 2 * A - 1,                 # = rank-biserial; + => high > low
    VD_A         = A
  )
}

HIGH <- "high"; LOW <- "low"   # check: unique(df_effect_avg$disruption)

es_overall <- grp_cmp(df_effect_avg, high = HIGH, low = LOW)  # alt = "greater" if H1 is high > low
es_overall

######## ######## ######## ######## ######## ######## ######## ######## 

######## look at shuffles too 
alu_shuffle  <- fread(paste0(PROJECT_DIR, "/results/paper_results/alu_highlow_shuff_scores.txt"))
alu_shuffle$disruption<-factor(alu_shuffle$disruption, levels=c('low', 'high'))

alu_shuffle_agg <- alu_shuffle %>%
  group_by(orig_idx, disruption) %>%
  summarize(mse_mean = mean(mse_mean, na.rm = TRUE), .groups = "drop")
# then run grp_cmp on alu_shuffle_agg

meds <- alu_shuffle %>%
  group_by(disruption) %>%
  summarize(med = mean(mse_mean, na.rm = TRUE), .groups = "drop")

pct_more <- 100 * (meds$med[meds$disruption == "high"] /
                     meds$med[meds$disruption == "low"] - 1)
pct_more   # -> "shuffling a disruptive Alu yields ~pct_mo

x <- alu_shuffle$mse_mean[alu_shuffle$disruption == "high"]
y <- alu_shuffle$mse_mean[alu_shuffle$disruption == "low"]
x <- x[!is.na(x)]; y <- y[!is.na(y)]
n1 <- length(x); n2 <- length(y)

w      <- wilcox.test(x, y, exact = FALSE)
U_high <- unname(w$statistic)     # U for 'high': #(high > low) pairs (+ 0.5 per tie)
U_low  <- n1 * n2 - U_high        # U for 'low' (the complement)
U_min  <- min(U_high, U_low)      # the value classic tables/textbooks report

c(U_high = U_high, U_low = U_low, U_min = U_min,
  n1 = n1, n2 = n2, max_U = n1 * n2)

pdf(paste0(PROJECT_DIR, "/figs/paper_figs/alu_highlow_shuffle_2.pdf"), width=70/25.4, height=60/25.4)
lwd<-0.3
col_val_highlow <- c('#6AA1C8', '#C15C6D')
  
ggplot(alu_shuffle, aes(x = disruption, y = mse_mean, fill=as.factor(disruption))) +
  geom_violin(trim = T, width=1, position=position_dodge(width=0.8), linewidth=lwd, alpha=0.8) +
  geom_boxplot(width = 0.2, outlier.size=0.3, position=position_dodge(width=0.8), linewidth=lwd) +
  scale_fill_manual(values = col_val_highlow) + scale_y_log10()+ 
  labs(x = "Disruption Category", y = "MSE shuffle") + 
  theme_classic() + theme_textsize_border 
dev.off()
