# Regenerate plots for the paper.
library(Seurat)
library(Matrix)
library(stm)
library(ggplot2)
library(viridis)
library(ggridges)
library(ggthemes)
# Part 1: Plots Matt makes
# Figure 1 A-C are generated from Lalli's place
# and I have tried to 
# Figure 1D
# 77 target genes 
SinaiCols4 <- c('#0CBFF1','#2C2A6F','#D6118A', '#ECEDED')

NPC.target.sub = readRDS('./Lalli_data_store/NPC_target_sub.RDS')
NPC_target_with_guide = readRDS('Lalli_data_store/NPC_target_with_guide.rds')
# 'SUV420H1' is 'KMT5B'
NPC_target_with_guide$gene_level[NPC_target_with_guide$gene_level == 'SUV420H1'] = 'KMT5B'
levels(NPC_target_with_guide$sgRNA)[grep('SUV420H1', levels(NPC_target_with_guide$sgRNA))] = c('KMT5B_G1', 'KMT5B_G2', 'KMT5B_G3')


beta.mimosca = read.csv('First_paper/scRNA_Mimosca/RESULTS/GeneByGeneLevelBetaMimosca.csv')
rownames(beta.mimosca) = beta.mimosca$X
beta.mimosca = beta.mimosca[, -1]
beta.mimosca = data.matrix(beta.mimosca)


# 77 targets and 1 non-targeting 
ASD.targets = c('NONTARGETING', 'ADNP', 'ANKRD11', 'ARID1B', 'ASH1L', 'ASXL3', 'AUTS2', 'BCL11A', 'CHAMP1', 'CHD2', 
                'CHD8', 'CREBBP', 'CTCF', 'CTNNB1', 'CUL3', 'DDX3X', 'DEAF1', 'DNMT3A', 'DYRK1A', 'EBF3', 'EHMT1', 
                'FMR1', 'FOXP1', 'FOXP2', 'GRIN2B', 'IRF2BPL', 'KDM5B', 'KDM6B', 'KMT2A', 'KMT2C', 'KMT2E', 'KMT5B', 
                'MBD5', 'MECP2', 'MED13', 'MED13L', 'MEF2C', 'MYT1L', 'NACC1', 'NF1', 'NIPBL', 'NRXN1', 'NSD1', 'PAX5', 
                'PHF12', 'PHF21A', 'POGZ', 'PPP2R5D', 'PTEN', 'RAI1', 'RFX3', 'RORB', 'SATB1', 'SATB2', 'SCN2A', 'SETD5', 
                'SHANK2', 'SHANK3', 'SIN3A', 'SKI', 'SLC6A1', 'SMARCA2', 'SMARCC2', 'SYNGAP1', 'TBL1XR1', 'TBR1', 'TCF20', 
                'TCF4', 'TCF7L2', 'TLK2', 'TRAF7', 'TRIO', 'TRIP12', 'TSC1', 'VEZF1', 'WAC', 'ZBTB20', 'ZMYND11')
# Among which, TBR1 is not captured by 
load('Lalli_data_store/ASD_target.RData')
# Select 50 non-targeting genes
non.target.genes = sample(setdiff(rownames(z.beta), ASD.targets), 50)

beta.mimosca.sub = beta.mimosca[c(non.target.genes, intersect(ASD.targets, rownames(z.beta))), ASD.targets]
# Z-normalized the coefficients 
z.beta = apply(beta.mimosca.sub, 2, function(x){ (x-mean(x))/sd(x) })
z.beta.sub = z.beta[intersect(ASD.targets, rownames(z.beta)), ASD.targets]
z.beta.sub = rbind(colMeans(z.beta[non.target.genes, ]), z.beta.sub)
# z.beta.sub[z.beta.sub < -8] = -6
rownames(z.beta.sub)[1] = 'NONTARGETING'
pdf('First_paper/Figures/Figure1D.pdf', width = 6.5, height = 6)
aheatmap(z.beta.sub, Rowv = NA, Colv = NA, border_color = NA, breaks = 0)
dev.off()

#ASD_targets = names(which(table(NPC_target_with_guide$gene_level) > 60))
#ASD.targets = c('NONTARGETING', ASD.targets[ASD.targets != 'NONTARGETING'])


load('First_paper/Data_Store/meta_data_cds_full.RData')
NPC_target_with_guide@reductions$umap@cell.embeddings[, 'UMAP_1'] = meta.data.cds.full$UMAP_1[match(colnames(NPC_target_with_guide), meta.data.cds.full$cell)]
NPC_target_with_guide@reductions$umap@cell.embeddings[, 'UMAP_2'] = meta.data.cds.full$UMAP_2[match(colnames(NPC_target_with_guide), meta.data.cds.full$cell)]
NPC_target_with_guide$broad_cluster = meta.data.cds.full$cluster[match(colnames(NPC_target_with_guide), meta.data.cds.full$cell)]
NPC_target_with_guide$pseudotime = meta.data.cds.full$pseudotime[match(colnames(NPC_target_with_guide), meta.data.cds.full$cell)]

pdf('First_paper/Figures/Fig2_umap_npc_broad_cluster.pdf', width = 12, height = 7)
DimPlot(subset(NPC_target_with_guide, broad_cluster != 4), reduction = "umap", group.by="broad_cluster", pt.size = 0.8) + scale_color_manual(values = SinaiCols4[1:3])
dev.off()

pdf('First_paper/Figures/Fig3_umap_npc_pseudotime.pdf', width = 12, height = 7)
FeaturePlot(subset(NPC_target_with_guide, pseudotime != "Inf"), reduction = "umap", features="pseudotime", pt.size = 0.8) + scale_color_viridis(option="inferno")
dev.off()

#### Figure 1D #########
# Working Knockout
NPC_target_all@reductions$umap@cell.embeddings[, 'UMAP_1'] = cds@int_colData$reducedDims$UMAP[match(colnames(NPC_target_all), cell.names), 1]
NPC_target_all@reductions$umap@cell.embeddings[, 'UMAP_2'] = cds@int_colData$reducedDims$UMAP[match(colnames(NPC_target_all), cell.names), 2]
NPC_target_all$broad_cluster = pData(cds)$cluster[match(colnames(NPC_target_all), cell.names)]

NPC.counts = NPC_target_with_guide@assays$RNA@counts
NPC.data = NPC_target_with_guide@assays$RNA@data


non.target.genes = sample(setdiff(rownames(NPC_target_with_guide), ASD.targets), 50)
NPC.data.sub = NPC.data[intersect(c(ASD.targets[-1], non.target.genes), rownames(NPC.data)), ]
Z.genes.floor = function(data, batch, floor = 0){
  batch.names = names(table(batch))
  z.score = matrix(NA, nrow = nrow(data), ncol = ncol(data))
  for(b in batch.names){
    data.sub = data[, batch == b]
    z.sub = (data.sub - rowMeans(data.sub))/(apply(data.sub, 1, sd) + floor)
    z.score[, batch == b] = data.matrix(z.sub)
  }
  rownames(z.score) = rownames(data);
  colnames(z.score) = colnames(data)
  return(z.score)
}
NPC.z.score = Z.genes.floor(NPC.data.sub, NPC_target_with_guide$lane, floor = 0.1)

NPC.z.score = (NPC.data.sub - rowMeans(NPC.data.sub)) /apply(NPC.data.sub, 1, sd)
a = sapply(ASD.targets, function(target){
  rowMeans(NPC.z.score[, NPC_target_with_guide$gene_level == target])
})


cds.cpm.sub = cds.log.cpm[c(ASD.targets[-1], non.target.genes), ]
avg.exprs.mat = sapply(ASD.targets, function(t){
  Matrix::rowMeans(cds.cpm.sub[, which(NPC_target_with_guide$gene_level == t)])
})

z.score.exprs.mat = t(apply(avg.exprs.mat, 1, function(x){
  (x-mean(x))/sd(x)
}))
z.score.exprs.mat = rbind(NONTARGETING = colMeans(z.score.exprs.mat[77:126, ]), z.score.exprs.mat[1:76, ])
heatmap(z.score.exprs.mat, Rowv = NA, Colv = NA)

pdf('plots/heatmap_Lalli_exprs.pdf', width = 8, height = 8)
aheatmap(z.score.exprs.mat, Rowv = NA, Colv = NA, border_color = 'grey')
dev.off()

######## Figure 2B,C,D UMAP of the marker genes ####
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_1'] = meta.data.cds.full$UMAP_1[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_2'] = meta.data.cds.full$UMAP_2[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]

NPC.target.sub$broad_cluster = meta.data.cds.full$cluster[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$pseudotime = meta.data.cds.full$pseudotime[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]

# marker genes
markers.fig2 = c('NKX2-1', 'NKX6-2', 'DCX')
pdf('First_paper/Figures/Fig2_marker.pdf', width = 6, height = 9)
FeaturePlot(NPC.target.sub, reduction = 'umap', features = markers.fig2, ncol = 1)
dev.off()


markers.npc = list(FNG = c('SOX2', 'NES', 'FOXG1'), DCM = c('PAX6'), VF = c('ASCL1', 'DLX1', 'DLX2'), GBA = c('GAD1', 'GAD2'), 
                   GBA.prog = c('GFAP', 'MKI67', 'TUBB3', 'SYP'))

pdf('First_paper/Figures/Supp_NPC_markers.pdf', width = 8, height = 10)
FeaturePlot(NPC.target.sub, reduction = 'umap', features = unlist(markers.npc), ncol = 3)
dev.off()

###### Figure 2E: Barplot of cell clusters #######

# Chisq Test
# Use Effective ASD risk Genes
ASD.Eff = names(table(NPC.target.sub$gene_level))
ASD.Eff = c('NONTARGETING', ASD.Eff[ASD.Eff != 'NONTARGETING'])

chisq.pval.broadcluster = sapply(2:length(ASD.Eff), function(i){
  M = table(NPC.target.sub$gene_level, NPC.target.sub$broad_cluster)[ASD.Eff[c(1,i)], 1:3]
  M = cbind(rowSums(M[, 1:2]), M[, 3])
  Xsq <- chisq.test(M)
  Xsq$p.value})
names(chisq.pval.broadcluster) = ASD.Eff[-1]
adj.p.broadcluster = p.adjust(chisq.pval.broadcluster, method = 'bonferroni')
sum(adj.p.broadcluster < 0.01)
# 13 cells 
names(which(adj.p.broadcluster < 0.01))
# ANKRD11, CHD8(44 cells), DEAF1, KMT2A, KMT2C, MBD5, MED13L, MYT1L(123), NIPBL(270), RFX3, TCF4, TCF7L2(92), ZBTB20
temp = round(cbind(chisq.pval.broadcluster, adj.p.broadcluster), 3)
write.table(temp, quote = F, file = 'First_paper/Data_Store/barplot_pval.txt', sep = '\t')


library(reshape2)
meta.data.npc = NPC.target.sub@meta.data
meta.data.npc = meta.data.npc[meta.data.npc$gene_level %in% ASD.Eff & meta.data.npc$broad_cluster != 4, ]
fracs.clusters = apply(table(meta.data.npc$broad_cluster, meta.data.npc$gene_level), 2, function(x){x/sum(x)})
m.fracs.clusters = melt(fracs.clusters)
colnames(m.fracs.clusters) = c('Cluster', 'Target', 'Fractions')
m.fracs.clusters$Target = factor(m.fracs.clusters$Target, levels = colnames(fracs.clusters)[order(fracs.clusters[3, ])])
m.fracs.clusters$Cluster = factor(m.fracs.clusters$Cluster, levels = c(3, 2, 1))
a = ifelse(colnames(fracs.clusters)[order(fracs.clusters[3, ])] %in% names(which(adj.p.broadcluster < 0.01)), "red", "black")

pdf('First_paper/Figures/Fig2_barplot.pdf', width = 12, height = 7)
ggplot(m.fracs.clusters, aes(x = Target, y = Fractions, fill = Cluster)) +
  geom_bar(stat = "identity") + theme_minimal(base_size = 10) + 
  theme(axis.text.x = element_text(angle = 70, color = a, vjust = 1, hjust = 1)) + 
  ggtitle('Fractions of Cluster') + scale_fill_manual(values = SinaiCols4[3:1]) 
dev.off()

NPC_target_not_zero <- subset(NPC.target.sub, pseudotime != "Inf")

plotData <- as.data.frame(cbind(as.numeric(NPC_target_not_zero$pseudotime), as.character(NPC_target_not_zero$gene_level)))

colnames(plotData) <- c("Pseudotime","Guide")
plotData <- subset(plotData, !(Guide %in% cluster0_genes))  ## already have a story for these genes. BTW this is how to use NOT IN 

pd2 <- plotData %>% group_by(Guide) %>% summarise(mean = mean(as.numeric(Pseudotime)), sem = sd(as.numeric(Pseudotime))/sqrt(n())) 


pd2$zscore <- scale(pd2$mean)[,1] 
pd2$color = 'darkgray'
pd2$color[pd2$zscore < qnorm(0.025)] = '#51127c'
pd2$color[pd2$zscore > -qnorm(0.025)] = '#fde725'
pd2$color[pd2$Guide == 'NONTARGETING'] = 'white'
a = rep('black', nrow(pd2))
a[1:3] = '#51127c'; a[22] = 'white'; a[59:61] = '#fde725'

pdf('First_paper/Figures/Fig3_dots_pseudotime.pdf', width = 10, height = 7)
g <- ggplot(pd2, aes(reorder(Guide, zscore), zscore, fill = color)) #color=zscore))
g + geom_point(shape=21, size=5) + xlab('Target Gene') + ylab('Average Pseudotime (Z-score)') + 
  theme_minimal(base_size = 16)  + scale_fill_manual(values = c('#51127c','#fde725', 'darkgray',  'white' )) + 
  theme(axis.text.x = element_text(angle = 70, hjust=1, color = a), legend.position = "none") 
dev.off()

pt_alter <- c("NF1","MED13", "MBD5", "NONTARGETING","POGZ", "TRAF7")

plotData$Guide <- relevel(as.factor(plotData$Guide), "TRAF7")
plotData$Guide <- relevel(as.factor(plotData$Guide), "POGZ")
plotData$Guide <- relevel(as.factor(plotData$Guide), "NONTARGETING")
plotData$Guide <- relevel(as.factor(plotData$Guide), "MBD5")
plotData$Guide <- relevel(as.factor(plotData$Guide), "MED13")
plotData$Guide <- relevel(as.factor(plotData$Guide), "NF1")

pdf('First_paper/Figures/Fig3_ridgeplot.pdf', width = 6, height = 8)
ggplot(subset(plotData, Guide%in%pt_alter), aes(x=as.numeric(Pseudotime), y=Guide, fill=Guide)) + geom_density_ridges(scale=1.5) + 
  scale_fill_colorblind() + theme_minimal(base_size = 14)#scale_fill_viridis_d(option="magma") 
dev.off()

NPC_target_all_12 = NPC_target_all[, NPC_target_all$broad_cluster %in% c(1, 2)]
npc.avg.pt = sapply(targets, function(t){ 
  mean(NPC_target_all_12$pseudotime[which(NPC_target_all_12$gene_level == t)])
})
z.score.npc = (npc.avg.pt - mean(npc.avg.pt))/sd(npc.avg.pt)
df.in.z.score = data.frame(z.score = z.score.npc, targets = factor(names(z.score.npc), names(z.score.npc)[order(z.score.npc)]), 
                           label = 'na')
df.in.z.score$label[df.in.z.score$targets == 'NONTARGETING'] = 'NONTARGETING'
df.in.z.score$label[df.in.z.score$z.score < qnorm(0.05)] = 'neg'
df.in.z.score$label[df.in.z.score$z.score > qnorm(0.95)] = 'pos'

#SinaiCols4 = "#0CBFF1" "#2C2A6F" "#D6118A" "#ECEDED"

pdf('First_paper/Figures/NPC_inhibitory_pseudotime_dots.pdf', width = 8, height = 5)
ggplot(df.in.z.score, aes(targets, z.score, fill = label)) + geom_point() + 
  geom_point(shape=21, size = 4) + 
  scale_fill_manual(values = c(SinaiCols4[4], 'black', SinaiCols4[1], SinaiCols4[3])) + 
  theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle = 45), legend.position = 'none') + 
  ylab('Average Pseudotime (Z-score)') + xlab('Target Gene') + ggtitle('Inhibitory Neuron')
dev.off()

#### K-S test results
NPC.target.sub = readRDS('./Lalli_data_store/NPC_target_sub.RDS')
load('First_paper/Data_Store/meta_data_cds_full.RData')
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_1'] = meta.data.cds.full$UMAP_1[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_2'] = meta.data.cds.full$UMAP_2[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$broad_cluster = meta.data.cds.full$cluster[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$pseudotime = meta.data.cds.full$pseudotime[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]

# remove lane 2 and 7
NPC.target.sub = NPC.target.sub[, NPC.target.sub$lane %in% c(1, 3:6, 8)]
# store meta data first
dim(NPC.target.sub@meta.data)


target.names.npc = names(table(NPC.target.sub$gene_level))
target.names.npc = c('NONTARGETING', target.names.npc[target.names.npc!= 'NONTARGETING'])

meta.pt = NPC.target.sub@meta.data
meta.pt = meta.pt[meta.pt$pseudotime < Inf, ]
x0 = meta.pt$pseudotime[meta.pt$gene_level == 'NONTARGETING']
ks.pval = sapply(target.names.npc[-1], function(t){
  x = meta.pt$pseudotime[meta.pt$gene_level == t]
  ks.test(x0, x, simulate.p.value = T)$p.value
})
#which(p.adjust(ks.pval) < 0.05)
#ARID1B   BAZ2B  BCL11A    CHD2     CIC   KMT2C    LEO1   MED13  MED13L    PHF3   SETD5 TBL1XR1    TBR1 
ks.pval = sapply(target.names.npc[-1], function(t){
  x = meta.pt$pseudotime[meta.pt$gene_level == t]
  fasano.franceschini.test(data.frame(pt= x0), data.frame(pt = x))$p.value
})
# We can draw a figure for those genes
# ks.sig.genes = names(which(ks.pval <= 0.1))
ks.sig.genes = rev(c('TRAF7', 'RFX3', 'TCF7L2', 'CHAMP1', 'PPP2R5D', 'CTCF', 
                 'ZMYND11', 'NONTARGETING', 'CHD2', 'PHF21A', 'MED13', 'NF1'))
df.ps = meta.pt[meta.pt$gene_level %in% ks.sig.genes, ]
df.ps = df.ps[, c('gene_level', 'pseudotime')]
avg.pt = sapply(target.names.npc, function(t){
  mean(meta.pt$pseudotime[meta.pt$gene_level == t])
})
df.ps$gene_level = factor(df.ps$gene_level, levels = ks.sig.genes)

library(ggridges)
pdf('First_paper/Figures/pseudotime_density_KS_NPC.pdf', width = 6, height = 12)
ggplot(df.ps, aes(pseudotime, gene_level, group = gene_level, fill = gene_level)) + geom_density_ridges() + 
  theme_minimal(base_size = 15) + theme(legend.position = 'none') + ggtitle('K-S test') + ylab('Targets')
dev.off()


# Part 2: The same figure I should make for Li et al. data
library(Seurat)
library(ggplot2)
# 1. load data: 

modules.asd = readRDS('CHOOSE/CHOOSE_ASD_Modules.rds')
# remove clusters = 32
modules.asd = modules.asd[, modules.asd$clusters != 32]
modules.asd = modules.asd[, modules.asd@reductions$umap@cell.embeddings[, 1] > -8]
#pdf('First_paper/Figures/umap_li_celltypes.pdf', width = 10, height = 8)
DimPlot(modules.asd, group.by = 'celltype_cl_coarse2') + ggtitle('Li Dataset')
#dev.off()

meta.data = modules.asd@meta.data
target.names.choose = names(table(meta.data$gRNA))
target.names.choose = c('Control2', target.names.choose[target.names.choose != 'Control2'])


# heatmap: z-score of gene expression
# Get all the cells and the average for gRNA gene expression
# randomly choose 100 genes that are not non-targeted
non.target.genes = sample(rownames(modules.asd), 50)
avg.exprs.mat = sapply(target.names.choose, function(t){
  rowMeans(modules.asd@assays$RNA@data[c(target.names.choose[-1], non.target.genes), modules.asd$gRNA == t])
  })

z.score.exprs.mat = t(apply(avg.exprs.mat, 1, function(x){
  (x-mean(x))/sd(x)
}))
z.score.exprs.mat = rbind(Control2 = colMeans(z.score.exprs.mat[37:86, ]), z.score.exprs.mat[1:36, ])
heatmap(z.score.exprs.mat, Rowv = NA, Colv = NA)

pdf('plots/heatmap_CHOOSE_exprs.pdf', width = 8, height = 8)
aheatmap(z.score.exprs.mat, Rowv = NA, Colv = NA, main = 'All Cells')
dev.off()

# Separate two lineage for Excitatory Neurons and Inhibitory Neurons
avg.exprs.ex.mat = sapply(target.names.choose, function(t){
  rowMeans(modules.asd@assays$RNA@data[c(target.names.choose[-1], non.target.genes), modules.asd$gRNA == t & modules.asd$lineage == 'excitatory'])
})

avg.exprs.in.mat = sapply(target.names.choose, function(t){
  rowMeans( modules.asd@assays$RNA@data[c(target.names.choose[-1], non.target.genes), modules.asd$gRNA == t & modules.asd$lineage == 'inhibitory'] )
})

z.score.ex.exprs.mat = t(apply(avg.exprs.ex.mat, 1, function(x){
  (x-mean(x))/sd(x)
}))
z.score.ex.exprs.mat = rbind(Control2 = colMeans(z.score.ex.exprs.mat[37:86, ]), z.score.ex.exprs.mat[1:36, ])

z.score.in.exprs.mat = t(apply(avg.exprs.in.mat, 1, function(x){
  (x-mean(x))/sd(x)
}))
z.score.in.exprs.mat = rbind(Control2 = colMeans(z.score.in.exprs.mat[37:86, ], na.rm = T), z.score.in.exprs.mat[1:36, ])

pdf('plots/heatmap_CHOOSE_exprs_two_lineages.pdf', width = 8, height = 8)
aheatmap(z.score.ex.exprs.mat, Rowv = NA, Colv = NA, main = 'Excitatory')
aheatmap(z.score.in.exprs.mat, Rowv = NA, Colv = NA, main = 'Inhibitory')
dev.off()

# The Z-scores does not look good overall, or for separate lineages. 


# Plot out UMAPs, colored by lineage, state and cell types.
pdf('First_paper/Figures/CHOOSE_UMAP_celltype.pdf', width = 10, height = 8)
DimPlot(modules.asd, group.by = 'celltype_cl_coarse2') + ggtitle('Cell Type')
dev.off()

pdf('First_paper/Figures/CHOOSE_UMAP_lineage.pdf', width = 10, height = 8)
DimPlot(modules.asd, group.by = 'lineage') + scale_color_brewer(palette = 'Set2') + ggtitle('Lineage')
dev.off()

pdf('First_paper/Figures/CHOOSE_UMAP_state.pdf', width = 10, height = 8)
DimPlot(modules.asd, group.by = 'state') + scale_color_brewer(palette = 'Set3') + ggtitle('State')
dev.off()

# Plot out UMAPs, colored by marker genes. 
marker.genes = c('VIM', 'PAX6', 'ASPM', 'GFAP', 'ASCL1', 'MKI67', 'TOP2A', 'HOPX', 'EOMES', 'TAC3', 'BCL11B', 
                 'SATB2', 'NEUROD6', 'RORB', 'UNC5D', 'NR2F1', 'TLE4', 'FOXP2', 'OLIG1', 'OLIG2', 'PDGFRA', 'ADARB2', 
                 'CALB2', 'DLX2', 'DLX5', 'NR2F2', 'MEIS2', 'GAD1', 'S100B', 'APOE', 'ALDH1L1', 'AUTS2')

# In total we have 32 marker genes, we can plot 8 in 4 pictures



for(i in 1:4){
  pdf(paste0('First_paper/Figures/CHOOSE_marker_gene_expression_', i,'.pdf'), width = 8, height = 12)
  print(FeaturePlot(modules.asd, features = marker.genes[(1:8) + (i-1)*8], ncol = 2))
  dev.off()
}

# Target gene Fraction plot
# I will first contrast the excitatory neurons and inhibtory neurons
# Then I will contrast the neuron, progenitors and others
library(reshape2)
meta.data.choose = modules.asd@meta.data

fracs.state = apply(table(meta.data.choose$state, meta.data.choose$gRNA), 2, function(x){x/sum(x)})
m.fracs.state = melt(fracs.state)
colnames(m.fracs.state) = c('State', 'Target', 'Fractions')
m.fracs.state$Target = factor(m.fracs.state$Target, levels = colnames(fracs.state)[order(fracs.state[1, ])])
m.fracs.state$State = factor(m.fracs.state$State, levels = c('neuron', 'progenitor', 'other'))

pdf('First_paper/Figures/barplot_choose_fraction_of_states.pdf', width = 9, height = 7)
ggplot(m.fracs.state, aes(x = Target, y = Fractions, fill = State)) +
  geom_bar(stat = "identity") + theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle('Fractions of States') + scale_fill_brewer(palette = 'Set3')
dev.off()

# no others
fracs.state.pure = apply(table(meta.data.choose$state, meta.data.choose$gRNA)[-2, ], 2, function(x){x/sum(x)})
m.fracs.state.pure = melt(fracs.state.pure)
colnames(m.fracs.state.pure) = c('State', 'Target', 'Fractions')
m.fracs.state.pure$Target = factor(m.fracs.state.pure$Target, levels = colnames(fracs.state.pure)[order(fracs.state.pure[1, ])])

pdf('First_paper/Figures/barplot_choose_fraction_of_states_no_others.pdf', width = 9, height = 7)
ggplot(m.fracs.state.pure, aes(x = Target, y = Fractions, fill = State)) +
  geom_bar(stat = "identity") + theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle('Fractions of States') + scale_fill_brewer(palette = 'Set3')
dev.off()


# Chi-squre test for the states. 
tab.state = table(meta.data.choose$state, meta.data.choose$gRNA)


fracs.lineage = apply(table(meta.data.choose$lineage, meta.data.choose$gRNA), 2, function(x){x/sum(x)})
m.fracs.lineage = melt(fracs.lineage)
colnames(m.fracs.lineage) = c('Lineage', 'Target', 'Fractions')
m.fracs.lineage$Target = factor(m.fracs.lineage$Target, levels = colnames(fracs.lineage)[order(fracs.lineage[1, ])])


pdf('First_paper/Figures/barplot_choose_fraction_of_lineage.pdf', width = 9, height = 7)
ggplot(m.fracs.lineage, aes(x = Target, y = Fractions, fill = Lineage)) +
  geom_bar(stat = "identity") + theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle('Fractions of Lineage') + scale_fill_brewer(palette = 'Set2')
dev.off()

fracs.cell.type = apply(table(meta.data.choose$celltype_cl_coarse2, meta.data.choose$gRNA), 2, function(x){x/sum(x)} )
m.fracs.cell.type = melt(fracs.cell.type)
colnames(m.fracs.cell.type) = c('CellType', 'Target', 'Fractions')
m.fracs.cell.type$Target = factor(m.fracs.cell.type$Target, levels = colnames(fracs.lineage)[order(fracs.lineage[1, ])] )
m.fracs.cell.type$CellType = factor(m.fracs.cell.type$CellType, levels = colnames() )

pdf('First_paper/Figures/barplot_choose_fraction_of_celltype.pdf', width = 9, height = 7)
ggplot(m.fracs.cell.type, aes(x = Target, y = Fractions, fill = CellType)) +
  geom_bar(stat = "identity") + theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle('Fractions of Cell Type')
dev.off()



# 12-18-2023 (pick up from the line)
# perform (1) Chi-squre test (2) CHM-test
# Chi-square test 
chisq.pval.state = sapply(2:length(target.names.choose), function(i){
  M = table(modules.asd$gRNA, modules.asd$state)[target.names.choose[c(1,i)], ]
  Xsq <- chisq.test(M)
  Xsq$p.value})
names(chisq.pval.state) = target.names.choose[-1]

adj.p.state = p.adjust(chisq.pval.state)

chisq.pval.lineage = sapply(2:length(target.names.choose), function(i){
  M = table(modules.asd$gRNA, modules.asd$lineage)[target.names.choose[c(1, i)], ]
  Xsq <- chisq.test(M)
  Xsq$p.value
})
names(chisq.pval.lineage) = target.names.choose[-1]
adj.p.lineage = p.adjust(chisq.pval.lineage)

save(chisq.pval.state, chisq.pval.lineage, adj.p.lineage, adj.p.state, file = 'CHOOSE/p_value_q_value_fractions.RData')


# Pseudotime Analysis
# In Li et al. paper, they have removed the Astrocytes, ccRG, ccvRG from analysis

modules.asd.pt = modules.asd[, !is.na(modules.asd$pseudotime)]
pdf(file = 'First_paper/Figures/pseudotime_choose.pdf', width = 10, height = 6)
FeaturePlot(modules.asd.pt, features = 'pseudotime', cols = c('yellow', 'purple')) + ggtitle('All Neurons')
dev.off()

pdf(file = 'First_paper/Figures/pseudotime_choose_excitatory.pdf', width = 10, height = 6)
FeaturePlot(modules.asd.pt[, modules.asd.pt$lineage == 'excitatory' & modules.asd.pt@reductions$umap@cell.embeddings[, 2] - modules.asd.pt@reductions$umap@cell.embeddings[, 1]/2 > 1], 
            features = 'pseudotime', cols = c('yellow', 'purple')) + ggtitle('Excitatory Neuron')
dev.off()


pdf(file = 'First_paper/Figures/pseudotime_choose_inhibitory.pdf', width = 10, height = 6)
FeaturePlot(modules.asd.pt[, modules.asd.pt$lineage == 'inhibitory' & modules.asd.pt@reductions$umap@cell.embeddings[, 2] - modules.asd.pt@reductions$umap@cell.embeddings[, 1]/2 < 1 & 
                             modules.asd.pt@reductions$umap@cell.embeddings[, 2] + modules.asd.pt@reductions$umap@cell.embeddings[, 1] > -8], 
            features = 'pseudotime', cols = c('yellow', 'purple')) + ggtitle('Inhibitory Neuron')
dev.off()

modules.asd.pt.in = modules.asd.pt[, modules.asd.pt$lineage == 'inhibitory' & modules.asd.pt@reductions$umap@cell.embeddings[, 2] - modules.asd.pt@reductions$umap@cell.embeddings[, 1]/2 < 1 & 
                                     modules.asd.pt@reductions$umap@cell.embeddings[, 2] + modules.asd.pt@reductions$umap@cell.embeddings[, 1] > -8]
modules.asd.pt.ex = modules.asd.pt[, modules.asd.pt$lineage == 'excitatory' & modules.asd.pt@reductions$umap@cell.embeddings[, 2] - modules.asd.pt@reductions$umap@cell.embeddings[, 1]/2 > 1]



meta.data.pt.ex = modules.asd.pt.ex@meta.data  # 9743 cells
meta.data.pt.in = modules.asd.pt.in@meta.data  # 12108 cells
target.names.choose = names(table(meta.data.pt$gRNA))
target.names.choose = c('Control2', target.names.choose[target.names.choose != 'Control2'])


# The pseudotime for excitatory neurons for each targets
# calculate the average pseudotime for each target and then generate the Z-scores


# Average Pseudotime, Z-score
# Separate the excitatory and inhibitory neurons

ex.avg.pt = sapply(target.names.choose, function(t){ 
  mean(modules.asd.pt.ex$pseudotime[modules.asd.pt.ex$gRNA == t])
  })
z.score.ex = (ex.avg.pt - mean(ex.avg.pt))/sd(ex.avg.pt)
in.avg.pt = sapply(target.names.choose, function(t){
  mean(modules.asd.pt.in$pseudotime[modules.asd.pt.in$gRNA == t])
})
z.score.in = (in.avg.pt - mean(in.avg.pt))/sd(in.avg.pt)

df.in.z.score = data.frame(z.score = z.score.in, targets = factor(names(z.score.in), names(z.score.in)[order(z.score.in)]), 
                           label = 'na')
df.in.z.score$label[df.in.z.score$targets == 'Control2'] = 'control'
df.in.z.score$label[df.in.z.score$z.score < qnorm(0.05)] = 'neg'
df.in.z.score$label[df.in.z.score$z.score > qnorm(0.95)] = 'pos'

SinaiCols4 = c("#0CBFF1", "#2C2A6F", "#D6118A", "#ECEDED")
pdf('First_paper/Figures/choose_inhibitory_pseudotime_dots.pdf', width = 8, height = 6)
ggplot(df.in.z.score, aes(targets, z.score, fill = label)) + geom_point() + 
  geom_point(shape=21, size = 4) + 
  scale_fill_manual(values = c(SinaiCols4[4], 'black', SinaiCols4[1], SinaiCols4[3])) + 
  theme_minimal(base_size = 15, base_line_size = 0.2) + theme(axis.text.x = element_text(angle = 90), legend.position = 'none') + 
  ylab('Average Pseudotime (Z-score)') + xlab('Target Gene') + ggtitle('Inhibitory Neuron')
dev.off()
# ARID1B, CIC, DDX3X

library(ggridges)

df.ps.in = modules.asd.pt.in@meta.data[modules.asd.pt.in$gRNA %in% c('Control2', 'ARID1B', 'CIC', 'DDX3X'), ]
df.ps.in = df.ps.in[, c('gRNA', 'state', 'lineage', 'celltype_cl_coarse2', 'pseudotime', 'pseudotime_ranks')]

df.ps.in$gRNA = factor(df.ps.in$gRNA, levels = c('Control2', 'ARID1B', 'CIC', 'DDX3X'))

library(ggridges)
pdf('First_paper/Figures/pseudotime_density_in_zscore.pdf', width = 8, height = 6)
ggplot(df.ps.in, aes(pseudotime, gRNA, group = gRNA, fill = gRNA)) + geom_density_ridges() + 
  theme_minimal(base_size = 15) + theme(legend.position = 'none')
dev.off()



df.ex.z.score = data.frame(z.score = z.score.ex, targets = factor(names(z.score.ex), names(z.score.ex)[order(z.score.ex)]), 
                           label = 'na')
df.ex.z.score$label[df.ex.z.score$targets == 'Control2'] = 'control'
df.ex.z.score$label[df.ex.z.score$z.score < qnorm(0.05)] = 'neg'
df.ex.z.score$label[df.ex.z.score$z.score > qnorm(0.95)] = 'pos'

pdf('First_paper/Figures/choose_excitatory_pseudotime_dots.pdf', width = 8, height = 6)
ggplot(df.ex.z.score, aes(targets, z.score, fill = label)) + geom_point(shape=21, size = 4) + 
  scale_fill_manual(values = c(SinaiCols4[4], 'black', SinaiCols4[1], SinaiCols4[3])) + 
  theme_minimal(base_size = 15, base_line_size = 0.2) + theme(axis.text.x = element_text(angle = 90), legend.position = 'none') + 
  ylab('Average Pseudotime (Z-score)') + xlab('Target Gene') + ggtitle('Excitatory Neuron')
dev.off()
# KMT2A, SETD5, ADNP

df.ps.ex = modules.asd.pt.ex@meta.data[modules.asd.pt.ex$gRNA %in% c('Control2', 'KMT2A', 'SETD5', 'ADNP'), ]
df.ps.ex = df.ps.ex[, c('gRNA', 'state', 'lineage', 'celltype_cl_coarse2', 'pseudotime', 'pseudotime_ranks')]

df.ps.ex$gRNA = factor(df.ps.ex$gRNA, levels = c('Control2', 'KMT2A', 'SETD5', 'ADNP'))

pdf('First_paper/Figures/pseudotime_density_ex_zscore.pdf', width = 8, height = 6)
ggplot(df.ps.ex, aes(pseudotime, gRNA, group = gRNA, fill = gRNA)) + geom_density_ridges() + 
  theme_minimal(base_size = 15) + theme(legend.position = 'none')
dev.off()



# I can have an overall pseudotime combined both neurons
pt = modules.asd.pt$pseudotime;
meta.pt = modules.asd.pt@meta.data

avg.pt = sapply(target.names.choose, function(t){
  mean(pt[meta.pt$gRNA == t])
})
z.score.pt = (avg.pt - mean(avg.pt))/sd(avg.pt)

df.z.score = data.frame(z.score = z.score.pt, targets = factor(names(z.score.pt), names(z.score.pt)[order(z.score.pt)]), 
                        label = 'na')
df.z.score$label[df.z.score$targets == 'Control2'] = 'control'
df.z.score$label[df.z.score$z.score < qnorm(0.05)] = 'neg'
df.z.score$label[df.z.score$z.score > qnorm(0.95)] = 'pos'

pdf('First_paper/Figures/choose_all_neurons_pseudotime_dots.pdf', width = 8, height = 6)
ggplot(df.z.score, aes(targets, z.score, fill = label)) + geom_point(shape=21, size = 4) + 
  scale_fill_manual(values = c(SinaiCols4[4], 'black', SinaiCols4[1], SinaiCols4[3])) + 
  theme_minimal(base_size = 15, base_line_size = 0.2) + theme(axis.text.x = element_text(angle = 90), legend.position = 'none')  + 
  ylab('Average Pseudotime (Z-score)') + xlab('Target Gene') + ggtitle('All Neurons')
dev.off()
# ARID1B, TBL1XR1, ASXL3, SRSF11

df.ps.all = modules.asd.pt@meta.data[modules.asd.pt$gRNA %in% c('Control2', 'ARID1B', 'TBL1XR1', 'ASXL3', 'SRSF11'), ]
df.ps.all = df.ps.all[, c('gRNA', 'state', 'lineage', 'celltype_cl_coarse2', 'pseudotime', 'pseudotime_ranks')]

df.ps.all$gRNA = factor(df.ps.all$gRNA, levels = c('Control2', 'ARID1B', 'TBL1XR1', 'ASXL3', 'SRSF11'))

pdf('First_paper/Figures/pseudotime_density_all_zscore.pdf', width = 8, height = 6)
ggplot(df.ps.all, aes(pseudotime, gRNA, group = gRNA, fill = gRNA)) + geom_density_ridges() + 
  theme_minimal(base_size = 15) + theme(legend.position = 'none')
dev.off()


save(ex.avg.pt, in.avg.pt, z.score.ex, z.score.in, z.score.pt, avg.pt, file = 'CHOOSE/pseudotime_ex_in_z_score.RData')

# The 4 genes Lalli selected for his data
# MKI67, NES, TUBB3, DOX
marker.gene.exprs = modules.asd.pt[c(marker.genes, 'NES', 'TUBB3', 'DOX'), ]@assays$RNA@data
library(ggformula)
pdf(file = 'First_paper/Figures/choose_pseudotime_spline_marker_genes.pdf', width = 8, height = 5)
for(g in c(marker.genes, 'NES', 'TUBB3', 'DOX') ){
  df.exprs.pt = data.frame(exprs = marker.gene.exprs[g, ], pseudotime_ranks = meta.data.pt$pseudotime_ranks, pseudotime = meta.data.pt$pseudotime)
  gg.temp = ggplot(df.exprs.pt, aes(pseudotime_ranks, exprs, color = pseudotime_ranks)) + geom_point(size = 1) + 
    theme_minimal(base_size = 15) + geom_spline() + 
    scale_color_continuous(low = 'yellow', high = 'purple') + ggtitle(g) + theme(legend.position = 'none')
  print(gg.temp)
}
dev.off()

x0 = meta.pt$pseudotime[meta.pt$gRNA == 'Control2']

ks.pval = sapply(target.names.choose[-1], function(t){
  x = meta.pt$pseudotime[meta.pt$gRNA == t]
  ks.test(x0, x)$p.value
})
#which(p.adjust(ks.pval) < 0.05)
#ARID1B   BAZ2B  BCL11A    CHD2     CIC   KMT2C    LEO1   MED13  MED13L    PHF3   SETD5 TBL1XR1    TBR1 

# We can draw a figure for those genes
ks.sig.genes = names(which(p.adjust(ks.pval) < 0.05))

df.ps = meta.pt[meta.pt$gRNA %in% c('Control2', ks.sig.genes), ]
df.ps = df.ps[, c('gRNA', 'state', 'lineage', 'celltype_cl_coarse2', 'pseudotime', 'pseudotime_ranks')]

df.ps$gRNA = factor(df.ps$gRNA, levels = c('Control2', ks.sig.genes[order(avg.pt[ks.sig.genes])]))

library(ggridges)
pdf('First_paper/Figures/pseudotime_density_all_cells.pdf', width = 6, height = 10)
ggplot(df.ps, aes(pseudotime, gRNA, group = gRNA, fill = gRNA)) + geom_density_ridges() + 
  theme_minimal(base_size = 15) + theme(legend.position = 'none')
dev.off()


################ Dendrograms results. (Do not run anything again)
# NPC with all effective targets
# save(stm.NPC.all, file = 'stm_NPC_all_topic10.RData')
# Refer to stm_trail_choose.R
# Refer to structural_topic_model_for_oct17.R
# 
library(stm)
NPC.target.sub = readRDS('./Lalli_data_store/NPC_target_sub.RDS')
library(Seurat)
load('First_paper/Data_Store/meta_data_cds_full.RData')
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_1'] = meta.data.cds.full$UMAP_1[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_2'] = meta.data.cds.full$UMAP_2[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$broad_cluster = meta.data.cds.full$cluster[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$pseudotime = meta.data.cds.full$pseudotime[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]


load('stm_NPC_all_topic10.RData')
NPC.target.sub = NPC.target.sub[, NPC.target.sub$lane %in% c(1, 3:6, 8)]
# store meta data first
dim(NPC.target.sub@meta.data)

target.names.npc = names(table(NPC.target.sub$gene_level))
target.names.npc = c('NONTARGETING', target.names.npc[target.names.npc!= 'NONTARGETING'])
gene.names.npc = rownames(NPC.target.sub)

raw.counts = NPC.target.sub@assays$RNA@counts
meta.data.npc = NPC.target.sub@meta.data

avg.theta.npc = t(sapply(target.names.npc, function(t){
  colMeans(stm.NPC.all$theta[meta.data.npc$gene_level == t, ])
}))
colnames(avg.theta.npc) = paste0('Topic', 1:10)

topic.prop = data.frame(stm.NPC.all$theta)
colnames(topic.prop) = paste0('Topic', 1:10)
rownames(topic.prop) = colnames(NPC.target.stm)

NPC.target.stm = NPC.target.sub
NPC.target.stm@meta.data = cbind(NPC.target.stm@meta.data, topic.prop)
pdf('First_paper/Figures/UMAP_Lalli_allcells_10topic_props.pdf', width = 15, height = 8)
FeaturePlot(NPC.target.stm, features = paste0('Topic', 1:10))
dev.off()

SinaiCols4 = c("#0CBFF1", "#2C2A6F", "#D6118A", "#ECEDED")
prop.mtx = t(data.matrix(topic.prop)[order(NPC.target.stm$broad_cluster), ])
annotdf <- data.frame(row.names = colnames(prop.mtx), 
                      broad_cluster = NPC.target.stm$broad_cluster[order(NPC.target.stm$broad_cluster)])
mycolors <- SinaiCols4
names(mycolors) <- unique(annotdf$broad_cluster)
mycolors <- list(broad_cluster = mycolors)

library(pheatmap)
pdf('First_paper/Figures/heatmap_Lalli_allcells_10topic_props.pdf', width = 10, height = 6)
pheatmap(prop.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Lalli, Heatmap of topic proportions')
dev.off()

pdf('First_paper/Figures/heatmap_avg_theta_Lalli_allcells_10topics.pdf', width = 10, height = 10)
heatmap(avg.theta.npc, main = 'Lalli data, all cells')
dev.off()

pdf('First_paper/Figures/heatmap_avg_theta_Lalli_allcells_10topics_2.pdf', width = 8, height = 10)
pheatmap(avg.theta.npc, main = 'Lalli data, all cells', scale = 'none', border_color = NA)
dev.off()

pdf('First_paper/Figures/topic_prop_Lalli_allcells_10topics.pdf', width = 8, height = 6)
plot.STM(stm.NPC.all, type = "summary", xlim = c(0, 0.7), n = 5, main = 'Lalli data, all cells')
dev.off()

pdf('First_paper/Figures/graph_topic_Lalli_allcells_10topics.pdf', width = 6, height = 6)
mod.out.corr <- topicCorr(stm.NPC.all)
plot(mod.out.corr, main = 'Lalli all cells, topic correlation')
dev.off()
a = labelTopics(stm.NPC.all, n = 50)
write.csv(a$topics, file = 'First_paper/Data_Store/NPC_topic_top_50_genes.csv')



### Li et al. paper
modules.asd = readRDS('CHOOSE/CHOOSE_ASD_Modules.rds')
# remove clusters = 32
modules.asd = modules.asd[, modules.asd$clusters != 32]
modules.asd = modules.asd[, modules.asd@reductions$umap@cell.embeddings[, 1] > -8]
DimPlot(modules.asd, group.by = 'celltype_cl_coarse2') 

meta.data = modules.asd@meta.data
marker.genes = c('ERBB4', 'SCGN', 'H4C3', 'DLX6-AS1', 'IL1RAPL2', 
                 'AL589740.1', 'TMEM132D', 'NXPH1', 'CST3', 'SPARCL1', 'PDZRN3', 'DTL', 'AL358075.2', 'PTTG1')
npc.genes = read.csv('Lalli_data_store/NPC_genes.csv')
npc.genes = npc.genes$x
modules.genes = intersect(unique(c(npc.genes, marker.genes)), rownames(modules.asd))  # 1050


target.names.choose = names(table(meta.data$gRNA))
target.names.choose = c('Control2', target.names.choose[target.names.choose != 'Control2'])
modules.asd.sub = modules.asd[modules.genes, ]


# all targets for the CHOOSE dataset, 36 targets
load('CHOOSE/list_MODULES_full_topics_sub_genes_stm.RData')
names(list.MODULES.full) = paste0('Topic', 5*(1:6))
MODULES.full = list.MODULES.full$Topic20

gamma.gRNA = MODULES.full$mu$gamma[1:(length(target.names.choose)), ]
rownames(gamma.gRNA) = c('intercept', target.names.choose[-1])
heatmap(gamma.gRNA[-1, ], main = 'Li et al. data')

meta.data.choose = modules.asd.sub@meta.data


topic.prop = data.frame(MODULES.full$theta)
colnames(topic.prop) = paste0('Topic', 1:20)
rownames(topic.prop) = colnames(modules.asd.stm)

modules.asd.stm = modules.asd
modules.asd.stm@meta.data = cbind(modules.asd.stm@meta.data, topic.prop)

pdf('First_paper/Figures/UMAP_Li_allcells_20topic_props.pdf', width = 15, height = 12)
FeaturePlot(modules.asd.stm, features = paste0('Topic', 1:20))
dev.off()

library(RColorBrewer)

lineage.color = c("#66C2A5", "#FC8D62")
prop.lineage.mtx = t(data.matrix(topic.prop)[order(modules.asd.stm$lineage), ])
annotdf <- data.frame(row.names = colnames(prop.lineage.mtx), 
                      lineage = modules.asd.stm$lineage[order(modules.asd.stm$lineage)])
mycolors <- lineage.color
names(mycolors) <- unique(annotdf$lineage)
mycolors <- list(lineage = mycolors)
pdf('First_paper/Figures/heatmap_Li_allcells_lineage_20topic_props.pdf', width = 10, height = 6)
pheatmap(prop.lineage.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()


state.color = c( "#8DD3C7", "#FFFFB3", "#BEBADA")
prop.state.mtx = t(data.matrix(topic.prop)[order(modules.asd.stm$state), ])
annotdf <- data.frame(row.names = colnames(prop.state.mtx), 
                      state = modules.asd.stm$state[order(modules.asd.stm$state)])
mycolors <- state.color
names(mycolors) <- unique(annotdf$state)
mycolors <- list(state = mycolors)
pdf('First_paper/Figures/heatmap_Li_allcells_state_20topic_props.pdf', width = 10, height = 8)
pheatmap(prop.state.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()

library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(15)

prop.celltype.mtx = t(data.matrix(topic.prop)[order(modules.asd.stm$celltype_cl_coarse2), ])
annotdf <- data.frame(row.names = colnames(prop.celltype.mtx), 
                      celltype = modules.asd.stm$celltype_cl_coarse2[order(modules.asd.stm$celltype_cl_coarse2)])
mycolors <- hex
names(mycolors) <- unique(annotdf$celltype)
mycolors <- list(celltype = mycolors)
pdf('First_paper/Figures/heatmap_Li_allcells_celltype_20topic_props.pdf', width = 10, height = 8)
pheatmap(prop.celltype.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()


# There are the topics we want to use
# The average theta from control and targets
avg.theta.choose = t(sapply(target.names.choose, function(t){
  colMeans(MODULES.full$theta[meta.data.choose$gRNA == t, ])
}))
colnames(avg.theta.choose) = paste0('Topic', 1:20)

pdf('First_paper/Figures/heatmap_avg_theta_Li_allcells_20topics.pdf', width = 10, height = 10)
heatmap(avg.theta.choose, main = 'Li et al. data, all cells')
dev.off()
pdf('First_paper/Figures/heatmap_avg_theta_Li_allcells_20topics_2.pdf', width = 8, height = 10)
pheatmap(avg.theta.choose, main = 'Li data, all cells', scale = 'none', border_color = NA)
dev.off()

pdf('First_paper/Figures/topic_prop_Li_allcells_20topics.pdf', width = 8, height = 6)
plot.STM(MODULES.full, type = "summary", xlim = c(0, 0.7), n = 5, main = 'Li data, all cells')
dev.off()

pdf('First_paper/Figures/graph_topic_Li_allcells_20topics.pdf', width = 8, height = 8)
mod.out.corr <- topicCorr(MODULES.full)
plot(mod.out.corr, main = 'Li all cells, topic correlation')
dev.off()

a = labelTopics(MODULES.full, n = 50)
write.csv(a$topics, file = 'First_paper/Data_Store/Li_topic_top_50_genes.csv')

labels.modules <- labelTopics(MODULES.full, n=20)
labels.NPC <- labelTopics(stm.NPC.all, n = 20)
#only keep FREX weighting
topwords <- data.frame("features" = t(labels$frex))
#assign topic number as column name
colnames(topwords) <- paste("Topics", c(1:15))
#Return the result
topwords[1:5]


##### Common 20 targets: get ready for 

load('Lalli_data_store/list_NPC_common_stm_Oct17.RData')
load('CHOOSE/list_MODULES_common_stm_Oct17.RData')

# Refer common_target_stm.R
load('npc_modules_common_stm.RData')

NPC.target.sub = readRDS('./Lalli_data_store/NPC_target_sub.RDS')
load('./Lalli_data_store/ASD_target.RData')
load('./Lalli_data_store/working_KD.RData')
#log.npc.cpm = apply(NPC.target.sub@assays$RNA@counts, 2, function(x){log2(x/sum(x)*1e6+1)})
load('Lalli_data_store/G_mat.RData')
working_ASD = names(table(NPC.target.sub$gene_level))
G_mat = G_mat[colnames(NPC.target.sub), working_ASD]

load('First_paper/Data_Store/meta_data_cds_full.RData')
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_1'] = meta.data.cds.full$UMAP_1[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub@reductions$umap@cell.embeddings[, 'UMAP_2'] = meta.data.cds.full$UMAP_2[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$broad_cluster = meta.data.cds.full$cluster[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]
NPC.target.sub$pseudotime = meta.data.cds.full$pseudotime[match(colnames(NPC.target.sub), meta.data.cds.full$cell)]

# remove lane 2 and 7
NPC.target.sub = NPC.target.sub[, NPC.target.sub$lane %in% c(1, 3:6, 8)]
# store meta data first
dim(NPC.target.sub@meta.data)

target.names.npc = names(table(NPC.target.sub$gene_level))
target.names.npc = c('NONTARGETING', target.names.npc[target.names.npc!= 'NONTARGETING'])

gene.names.npc = rownames(NPC.target.sub)

# Load data from CHOOSE
modules.asd = readRDS('CHOOSE/CHOOSE_ASD_Modules.rds')
# remove clusters = 32
modules.asd = modules.asd[, modules.asd$clusters != 32]

DimPlot(modules.asd, group.by = 'celltype_cl_coarse2') 
DimPlot(modules.asd, group.by = 'seurat_clusters') 

meta.data = modules.asd@meta.data
marker.genes = c('ERBB4', 'SCGN', 'H4C3', 'DLX6-AS1', 'IL1RAPL2', 
                 'AL589740.1', 'TMEM132D', 'NXPH1', 'CST3', 'SPARCL1', 'PDZRN3', 'DTL', 'AL358075.2', 'PTTG1')

npc.genes = read.csv('Lalli_data_store/NPC_genes.csv')
npc.genes = npc.genes$x

modules.genes = intersect(unique(c(npc.genes, marker.genes)), rownames(modules.asd))  # 1050

target.names.choose = names(table(meta.data$gRNA))
target.names.choose = c('Control2', target.names.choose[target.names.choose != 'Control2'])

# Get common targets from two projects
target.names = intersect(target.names.choose, target.names.npc)

NPC.target.common = NPC.target.sub[, NPC.target.sub$gene_level %in% c("NONTARGETING", target.names)]  # only 5298 cells left
modules.asd.common = modules.asd[modules.genes, modules.asd$gRNA %in% c('Control2', target.names)]  # 16292 cells left

# FeaturePlot

topic.prop.modules = MODULES.common$theta
rownames(topic.prop.modules) = colnames(modules.asd.common)
colnames(topic.prop.modules) = paste0('Topic', 1:20)

modules.asd.common@meta.data = cbind(modules.asd.common@meta.data, topic.prop.modules)

pdf('First_paper/Figures/feature_plot_Li_20topics_common.pdf', width = 12, height = 8)
FeaturePlot(modules.asd.common, features = paste0('Topic', 1:20))
dev.off()

topic.prop.npc = NPC.common$theta
rownames(topic.prop.npc) = colnames(NPC.target.common)
colnames(topic.prop.npc) = paste0('Topic', 1:10)

NPC.target.common@meta.data = cbind(NPC.target.common@meta.data, topic.prop.npc)
pdf('First_paper/Figures/feature_plot_Lalli_10topics_common.pdf', width = 15, height = 8)
FeaturePlot(NPC.target.common, features = paste0('Topic', 1:10))
dev.off()


SinaiCols4 = c("#0CBFF1", "#2C2A6F", "#D6118A", "#ECEDED")
prop.mtx = t(data.matrix(topic.prop.npc)[order(NPC.target.common$broad_cluster), ])
annotdf <- data.frame(row.names = colnames(prop.mtx), 
                      broad_cluster = NPC.target.common$broad_cluster[order(NPC.target.common$broad_cluster)])
mycolors <- SinaiCols4
names(mycolors) <- unique(annotdf$broad_cluster)
mycolors <- list(broad_cluster = mycolors)

library(pheatmap)
pdf('First_paper/Figures/heatmap_Lalli_common_10topic_props.pdf', width = 10, height = 6)
pheatmap(prop.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Lalli, Heatmap of topic proportions')
dev.off()

lineage.color = c("#66C2A5", "#FC8D62")
prop.lineage.mtx = t(data.matrix(topic.prop.modules)[order(modules.asd.common$lineage), ])
annotdf <- data.frame(row.names = colnames(prop.lineage.mtx), 
                      lineage = modules.asd.common$lineage[order(modules.asd.common$lineage)])
mycolors <- lineage.color
names(mycolors) <- unique(annotdf$lineage)
mycolors <- list(lineage = mycolors)
pdf('First_paper/Figures/heatmap_Li_common_lineage_20topic_props.pdf', width = 10, height = 6)
pheatmap(prop.lineage.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()


state.color = c( "#8DD3C7", "#FFFFB3", "#BEBADA")
prop.state.mtx = t(data.matrix(topic.prop.modules)[order(modules.asd.common$state), ])
annotdf <- data.frame(row.names = colnames(prop.state.mtx), 
                      state = modules.asd.common$state[order(modules.asd.common$state)])
mycolors <- state.color
names(mycolors) <- unique(annotdf$state)
mycolors <- list(state = mycolors)
pdf('First_paper/Figures/heatmap_Li_common_state_20topic_props.pdf', width = 10, height = 8)
pheatmap(prop.state.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()

library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(15)

prop.celltype.mtx = t(data.matrix(topic.prop.modules)[order(modules.asd.common$celltype_cl_coarse2), ])
annotdf <- data.frame(row.names = colnames(prop.celltype.mtx), 
                      celltype = modules.asd.common$celltype_cl_coarse2[order(modules.asd.common$celltype_cl_coarse2)])
mycolors <- hex
names(mycolors) <- unique(annotdf$celltype)
mycolors <- list(celltype = mycolors)
pdf('First_paper/Figures/heatmap_Li_common_celltype_20topic_props.pdf', width = 10, height = 8)
pheatmap(prop.celltype.mtx, cluster_rows = TRUE, cluster_cols = FALSE, scale = 'none', border_color = NA, 
         color = colorRampPalette(c("lightyellow", "firebrick3"))(50), labels_col = '', 
         annotation_col = annotdf, annotation_colors = mycolors, main = 'Li, Heatmap of topic proportions')
dev.off()

meta.data.choose = modules.asd.common@meta.data
avg.theta.choose = t(sapply(c('Control2', target.names), function(t){
  colMeans(MODULES.common$theta[meta.data.choose$gRNA == t, ])
}))
colnames(avg.theta.choose) = paste0('Topic', 1:20)

meta.data.npc = NPC.target.common@meta.data
avg.theta.npc = t(sapply(c('NONTARGETING', target.names), function(t){
  colMeans(NPC.common$theta[meta.data.npc$gene_level == t, ])
}))
colnames(avg.theta.npc) = paste0('Topic', 1:10)

diff.avg.theta.choose = avg.theta.choose[-1, ] - avg.theta.choose[rep(1, length(target.names)), ]
diff.avg.theta.npc = avg.theta.npc[-1, ] - avg.theta.npc[rep(1, length(target.names)), ]

par(mfrow = c(1, 2))
plot(hclust(dist(diff.avg.theta.npc)), main = 'Lalli dataset')
plot(hclust(dist(diff.avg.theta.choose)), main = 'Li dataset')

pdf('First_paper/Figures/heatmap_diff_theta_Li_common_20topics.pdf', width = 6, height = 4)
pheatmap(t(diff.avg.theta.choose), main = 'Li data, common targets, diff prop', scale = 'none', border_color = NA)
dev.off()

pdf('First_paper/Figures/heatmap_diff_theta_Lalli_common_10topics.pdf', width = 6, height = 4)
pheatmap(t(diff.avg.theta.npc), main = 'Lalli data, common targets, diff prop', scale = 'none', border_color = NA)
dev.off()




pdf('First_paper/Figures/topic_prop_Li_common_20topics.pdf', width = 8, height = 6)
plot.STM(MODULES.common, type = "summary", xlim = c(0, 0.7), n = 5, main = 'Li data, common')
dev.off()

pdf('First_paper/Figures/topic_prop_Lalli_common_10topics.pdf', width = 8, height = 6)
plot.STM(NPC.common, type = "summary", xlim = c(0, 0.7), n = 5, main = 'Lailli data, common')
dev.off()


pdf('First_paper/Figures/graph_topic_Li_common_20topics.pdf', width = 8, height = 8)
mod.out.corr <- topicCorr(MODULES.common)
plot(mod.out.corr, main = 'Li common_targets, topic correlation')
dev.off()

pdf('First_paper/Figures/graph_topic_Lalli_common_10topics.pdf', width = 8, height = 8)
mod.out.corr <- topicCorr(NPC.common)
plot(mod.out.corr, main = 'Lalli common_targets, topic correlation')
dev.off()


a = labelTopics(MODULES.common, n = 50)
write.csv(a$topics, file = 'First_paper/Data_Store/Li_common_target_topic_top_50_genes.csv')

a = labelTopics(NPC.common, n = 50)
write.csv(a$topics, file = 'First_paper/Data_Store/Lalli_common_target_topic_top_50_genes.csv')
