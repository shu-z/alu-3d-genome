# Paths - edit for your environment
PROJECT_DIR <- Sys.getenv("ALU_PROJECT_DIR", "/Users/shu/pollard_lab/alu")
DOWNLOADS_DIR <- "/Users/shu/Downloads"
PHYLO_OUT_DIR <- "/Users/shu/results"

library(data.table)
library(ape)
library(msa)
library(phangorn)
library(ggplot2)
library(ggtree)
library(cowplot)


stepped_colors<-stepped(24)
stepped2_colors<-stepped2(24)
stepped3_colors<-stepped3(24)


# info in phylo trees in R: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

alu_ins<-fread(paste0(PROJECT_DIR, "/results/20251022_topalu_insertAlu.txt"))
#rename old mse columns 
cols_to_rename <-  c("mse_mean", "corr_mean", 'mse_median', 'corr_median', "mse_high", "corr_high")
setnames(alu_ins, old = cols_to_rename, new = paste0("DEL_", cols_to_rename))

#keep meta columns 
meta_cols <- c("row_idx", "GC", 'repName', 'DEL_mse_mean', 'DEL_corr_mean', 'DEL_mse_high', 'DEL_corr_high', "SVLEN")

#melt columns(only mse_ and corr_)
melt_cols <- grep("^(mse|corr)_", names(alu_ins), value = TRUE)
alu_ins_long <- melt( alu_ins, id.vars = meta_cols, measure.vars = melt_cols,
                      variable.name = "metric", value.name = "value")

# Extract metric type and family name
alu_ins_long[, c("metric_type", "family") := tstrsplit(metric, "_", keep = c(1,2))]

#pull out alignscore
align_cols <- grep("^alignscore_", names(alu_ins), value = TRUE)
align_dt <- melt(alu_ins, id.vars = meta_cols, measure.vars = align_cols, 
                 variable.name = "align_metric",value.name = "alignscore")[, family := tstrsplit(align_metric, "_", keep = 2)]

#merge mse/corr with alignment information 
alu_ins_long <- merge(alu_ins_long, align_dt[, .(row_idx, family, alignscore)],by = c("row_idx", "family"), all.x = TRUE)

#edit corr so its 1-corr
alu_ins_long[ ,value:=(ifelse (metric_type=='corr', (1-value), value))]


#add in information on if the exact alus are the same, or if the families match
alu_ins_long[ ,repName_short:=(substr(repName, 1, 4))]
alu_ins_long[ ,family_short:=(substr(family, 1, 4))]
alu_ins_long[ ,repName_family_match:=ifelse(repName==family, 1, 0)]
alu_ins_long[ ,repName_family_short_match:=ifelse(repName_short==family_short, 1, 0)]

alu_ins_long$repName_family_match<-as.factor(alu_ins_long$repName_family_match)
alu_ins_long$repName_family_short_match<-as.factor(alu_ins_long$repName_family_short_match)


#define colors for repName and family
unique_repName<-unique(alu_ins_long$repName)
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


consensus_seq<-fread(paste0(DOWNLOADS_DIR, "/20240129_alu_consensus_seq_data.txt"), sep='\t')
# Set names
seqs <- consensus_seq$Alu_seq
names(seqs) <- consensus_seq$Alu_name

# Convert character strings to DNAbin
dna <- as.DNAbin(strsplit(seqs, split = ""))

# Convert to DNAStringSet
dna_strings <- DNAStringSet(seqs)
names(dna_strings) <- consensus_seq$Alu_name

# Align sequences
aln <- msa(dna_strings, method = "Muscle")  # or method = "ClustalW"
aln_mat <- as.DNAbin(aln)

dna2 <- as.phyDat(aln_mat)

#############################################
# Parsimony tree 

#look into which model to use 
#d <- dist.dna(dna, model = "raw")   # or model="JC69", "K80", etc.
tre.ini <- nj(dist.dna(aln_mat,model="raw"))
tre.ini
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)

pdf(paste0(DOWNLOADS_DIR, "/tree_parsimony.pdf"), height=10)
# plot(tre.pars, type="unr", show.tip=FALSE, edge.width=2)
# title("Maximum-parsimony tree")
# tiplabels(tre.pars$tip.label, cex=.5, fg="transparent")


tre.pars <- root(tre.pars,'FAM')
node_to_flip <- getMRCA(tre.pars, c("AluSx1", "AluSc5"))
#tre4 <- ladderize(tre4, right=F)
tre.pars <- ape::rotate(tre.pars, node=node_to_flip)

plot(tre.pars, show.tip=FALSE, edge.width=2)

title("Parsimony tree")
tiplabels(tre.pars$tip.label, cex=.5, fg="transparent")
axisPhylo()
dev.off()

#############################################
# MLE tree no selection 

tre.ini <- nj(dist.dna(aln_mat,model="TN93"))
fit.ini <- pml(tre.ini, dna2, k=4)

#optBf adjusts nucleotide frequencies to fit data
#optQ optimizes substitution model parameters
#optGamma optimizes shape parameter 
fit <- optim.pml(fit.ini, optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)

pdf(paste0(DOWNLOADS_DIR, "/tree_MLE.pdf"))

treMLE <- root(fit$tree,'FAM')
#tre4 <- ladderize(tre4)

dist_matrix_MLE<- cophenetic.phylo(treMLE)
write.csv(dist_matrix_MLE, paste0(PHYLO_OUT_DIR, "/phylo/alu_distances_MLE_TN93.csv"), row.names = TRUE)

plot(treMLE, show.tip=FALSE, edge.width=2)
title("Maximum-likelihood tree")
tiplabels(treMLE$tip.label, cex=.5, fg="transparent")
axisPhylo()
dev.off()

##############################################
# MLE tree with model selection

# Starting NJ tree
tre.ini <- nj(dist.dna(aln_mat, model="TN93"))

# Model selection
mt <- modelTest(dna2, tree=tre.ini, model=c("JC","K80","HKY","TN93","GTR"))

# Inspect model results
print(mt[order(mt$BIC),])
cat("Best model (AIC):", mt$Model[which.min(mt$AIC)], "\n")
cat("Best model (BIC):", mt$Model[which.min(mt$BIC)], "\n")

# Fit best model
fit.ini <- as.pml(mt)  # uses best BIC model by default

# Optimize (k=4 rate categories, invariant sites)
fit <- optim.pml(fit.ini, 
                 optNni=TRUE, 
                 optBf=TRUE, 
                 optQ=TRUE, 
                 optGamma=TRUE,
                 optInv=TRUE)

# Inspect fitted model parameters
print(fit)

# Bootstrap support
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)

# Root on FAM outgroup
treMLE <- root(fit$tree, "FAM")

# Cophenetic distance matrix
dist_matrix_MLE <- cophenetic.phylo(treMLE)
#write.csv(dist_matrix_MLE, paste0(PROJECT_DIR, "/results/phylo/20260320_alu_distances_MLE_bestfit.csv"), row.names=TRUE)



# Extract just the family for each tip if your tip labels are more complex
tip_families <- treMLE$tip.label  # or extract substring if needed

# Match colors from your family_color_df
tip_colors <- family_color_df$family_color[match(tip_families, family_color_df$family)]




pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260427_tree_MLE_bestfit_color.pdf"), width=3, height=7)

#bs controls branch length
plotBS(treMLE, bs,type="phylogram", show.tip=FALSE, 
       edge.width = 0.5, cex = 0.5,
       x.lim = c(0, max(node.depth.edgelength(treMLE)) * 2.2))

# Overlay tip labels with colors
# Add tip labels with colored boxes
tiplabels(treMLE$tip.label,
          tip=1:length(treMLE$tip.label),
          cex=0.6,
          frame="rect",         
          bg=adjustcolor(tip_colors[i], alpha.f = 0.7),           # fill the rectangle with your family colors
          col="black", lwd=0.1)        

#title("MLE tree")

dev.off()
graphics.off()





pdf(paste0(PROJECT_DIR, "/figs/paper_figs/20260428_tree_MLE_bestfit_color_rightalign.pdf"),
    width = 3, height = 7)

treMLE1 <- ladderize(treMLE, right = F) 
treMLE1 <- ape::rotate(treMLE1, node = 72) 
#treMLE1 <- ape::rotate(treMLE1, node = 71) 

#treMLE1 <- ape::rotate(treMLE1, node = 71) 
#treMLE1 <- ape::rotate(treMLE1, node = 79) 
#treMLE1 <- ape::rotate(treMLE1, node = 60) 

#relabel 
tip_families <- treMLE1$tip.label  # or extract substring if needed
tip_colors <- family_color_df$family_color[match(tip_families, family_color_df$family)]


#switch ordering of some nodes 

# Plot tree without tip labels, leaving extra horizontal room on the right
plotBS(treMLE1, bs, type = "phylogram", show.tip = FALSE,
       edge.width = 0.7, cex = 0.5,
       x.lim = c(0, max(node.depth.edgelength(treMLE1)) * 2.2),
       col='grey10',
       p=101 #nothing gets labeled this way 
       )

# Get the plotted coordinates of every node/tip
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
ntip   <- length(treMLE1$tip.label)
xx     <- lastPP$xx[1:ntip]   # actual x of each tip
yy     <- lastPP$yy[1:ntip]
x_max  <- max(xx)             # right edge — where we want labels aligned

# Optional: draw faint dotted lines from each tip to the alignment column
segments(x0 = xx, y0 = yy, x1 = x_max, y1 = yy,
         lty = 3, col = "grey20", lwd = 0.3)

# Draw the labels at x_max instead of at each tip's x
# adj = 0 means the LEFT edge of each box sits at x_max (so all left edges align)
# adj = 1 would right-align the right edges instead
text_width <- max(strwidth(treMLE1$tip.label, cex = 0.6)) * 1.05

for (i in seq_len(ntip)) {
  rect(xleft   = x_max,
       xright  = x_max + text_width,
       ybottom = yy[i] - 0.4,
       ytop    = yy[i] + 0.4,
       col     = adjustcolor(tip_colors[i], alpha.f = 0.8),
       border  = "grey20", lwd = 0.1)
  text(x = x_max, y = yy[i],
       labels = treMLE1$tip.label[i],
       adj = c(0, 0.5), cex = 0.6, col = "black")
}

dev.off()
graphics.off()


pdf(paste0(DOWNLOADS_DIR, "/tree_with_nodes.pdf"), width = 20, height = 4)
plot(treMLE1, cex = 0.5)
nodelabels(cex = 0.5, frame = "none", col = "red")
#tiplabels(cex = 0.5, frame = "none", col = "blue")
dev.off()