# SUPPLEMENTARY FIGURES -----
# 
# *Fig. S1A -----
seu_list <- lapply(seu_list, function(x) {x[['MT']] <- PercentageFeatureSet(x, pattern = "^MT-"); x})
lapply(seu_list, function(x) {colnames(x[[]])})

### pre-filter -----
lapply(VlnPlot(seu_list, features = c("nFeature_RNA", "nCount_RNA", "MT"), ncol = 3)) 

### filter params
cl <- 500 # min number of total reads per cell 
cl2 <- 25000 # max
fl <- 300 # min number of gene features per cell
fl2 <- 6000 # max
ml <- 15 # % cut off for mt

seu_list <- lapply(seu_list, function(x) {
  x <- subset(x, nCount_RNA > cl &
                nCount_RNA < cl2 &
                nFeature_RNA > fl & 
                nFeature_RNA < fl2 & 
                MT < ml)
})
### post-filter -----
lapply(VlnPlot(seu_list, features = c("nFeature_RNA", "nCount_RNA", "MT"), ncol = 3)) 

# *Fig. S1C-E - dotplots -----
# 
color <- colorRampPalette((c("steelblue", "white", "firebrick")))(50)
### NK cells -----
NKclus <- subset(integrated_tcr, ident = "NK cells")
DefaultAssay(NKclus) <- "RNA"
Idents(NKclus) <- "Cell_Type"
table(Idents(NKclus))

features_NKclus <- c("IL2RB", "GZMK", "IFNG", "CD160", "CXCR6", "CD69", 
                     "CX3CR1", "GNLY", "GZMB", "KLF2", "FGFBP2", "FCGR3A",
                     "NCAM1", "KLRC1", "CAPG", "TNFRSF18", "MKI67")

### Get data running wilcoxon
NKclus.genes <- wilcoxauc(NKclus, "Cell_Type", seurat_assay = "RNA")
head(NKclus.genes)
table(NKclus.genes$group)

# Filter to genes and order correctly
NKclus.genes <- NKclus.genes %>%
  filter(feature %in% features_NKclus) %>%
  arrange(factor(group, levels = levelsNK)) %>%
  arrange(factor(feature, levels = features_NKclus))
head(NKclus.genes, n = 10)

dp2 <- NKclus.genes %>% 
  mutate(`% Expressing` = pct_in) %>% 
  mutate(`Avg. Exp.` = logFC) %>%
  ggplot(aes(x=feature, y = group, color = logFC, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size = 11) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = color,  oob = scales::squish,  
                        limits = c(-2,2), name = "Avg. Exp.") + 
  scale_y_discrete(position = "left") +
  theme(legend.title = element_text(size = 10))
dp2$data$feature <- factor(dp2$data$feature, levels = features_NKclus)
dp2$data$group <- factor(dp2$data$group, levels = rev(levelsNK))
dp2

### Monocytes -----
MYEclus <- subset(integrated_tcr, ident = "Myeloid cells")
DefaultAssay(MYEclus) <- "RNA"
Idents(MYEclus) <- "Cell_Type"
table(Idents(MYEclus))

features_MYEclus <- c("LYZ", "FCN1", "MNDA", "CSF3R", 
                      "THBS1", "EREG", "CXCL10", "TIMP1", "CXCR6",
                      "CX3CR1", "CDKN1C", "LST1", "HMOX1",
                      "SASH1", "GPX1", "IFI30", "PSME2", "PID1",
                      "APOE", "APOC1", "C1QA", "TREM2", "FCGBP",
                      "CLEC9A", "CADM1", "IDO1", "SNX3", "CXCL9",
                      "LILRA4", "GZMB", "JCHAIN", "IRF7", "S100A8", "CXCL8")

### Get data running wilcoxon
MYEclus.genes <- wilcoxauc(MYEclus, "Cell_Type", seurat_assay = "RNA")
head(MYEclus.genes)
table(MYEclus.genes$group)

# Filter to genes and order correctly
MYEclus.genes <- MYEclus.genes %>%
  filter(feature %in% features_MYEclus) %>%
  arrange(factor(group, levels = levelsMYE)) %>%
  arrange(factor(feature, levels = features_MYEclus))
head(MYEclus.genes, n = 10)

dp2 <- MYEclus.genes %>% 
  mutate(`% Expressing` = pct_in) %>% 
  mutate(`Avg. Exp.` = logFC) %>%
  ggplot(aes(x=feature, y = group, color = logFC, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size = 11) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = color,  oob = scales::squish,  
                        limits = c(-2,2), name = "Avg. Exp.") + 
  scale_y_discrete(position = "left") +
  theme(legend.title = element_text(size = 10))
dp2$data$feature <- factor(dp2$data$feature, levels = features_MYEclus)
dp2$data$group <- factor(dp2$data$group, levels = rev(levelsMYE))
dp2

### B cells -----
Bclus <- subset(integrated_tcr, ident = "B cells")
DefaultAssay(Bclus) <- "RNA"
Idents(Bclus) <- "Cell_Type"
table(Idents(Bclus))

features_Bclus <- c("TCL1A", "IGHD", "IGHM", "FCER2", "IL4R", "ISG15", "ISG20", 
                    "IFIT3", "IFI6", "IRF7", "CXCR5", "CD69", "IL6", "ITGB1", 
                    "PLP2", "COCH", "CD82", "VIM", "CXCR6", "FCRL3", "FCRL5", 
                    "FGR", "HCK", "ITGAX", "JCHAIN", "MZB1", "IGKC", "IGHG1")


### Get data running wilcoxon
Bclus.genes <- wilcoxauc(Bclus, "Cell_Type", seurat_assay = "RNA")
head(Bclus.genes)
table(Bclus.genes$group)

# Filter to genes and order correctly
Bclus.genes <- Bclus.genes %>%
  filter(feature %in% features_Bclus) %>%
  arrange(factor(group, levels = levelsB)) %>%
  arrange(factor(feature, levels = features_Bclus))
head(Bclus.genes, n = 10)

dp2 <- Bclus.genes %>% 
  mutate(`% Expressing` = pct_in) %>% 
  mutate(`Avg. Exp.` = logFC) %>%
  ggplot(aes(x=feature, y = group, color = logFC, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size = 11) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = color,  oob = scales::squish, 
                        limits = c(-2,2), name = "Avg. Exp.") + 
  scale_y_discrete(position = "left") +
  theme(legend.title = element_text(size = 10))
dp2$data$feature <- factor(dp2$data$feature, levels = features_Bclus)
dp2$data$group <- factor(dp2$data$group, levels = rev(levelsB))
dp2

### T cells
Idents(integrated_tcr) <- "Main_Cell_Type"
Tclus <- subset(integrated_tcr, ident = "T cells")

features_Tclus <- c("LTB", "CCR7", "SELL", "IL7R", "TNFAIP3", "FOS", "CXCR6", 
                    "CXCR3", "DUSP1", "S100A11", "CISH", "FOXP3", "CTLA4", "TNFRSF4", 
                    "CD52", "LINC02446", "LDHB", "NOSIP", "CD248", "IFNG", "CD69", 
                    "XCL1", "ITGAE", "GZMK", "CD74", "CCL5", "CMC1", "CX3CR1", "GNLY", 
                    "GZMB", "PRF1", "MALAT1", "EOMES", "RNF213", "DDX3X", "IRF7", 
                    "MKI67", "TRAV1-2", "KLRB1", "RORC", "TRDC", "TRDV1")

### *Get avg. exp. and % data running wilcoxon
tclus.genes <- wilcoxauc(Tclus, "Cell_Type", seurat_assay = "RNA")
head(tclus.genes)
table(tclus.genes$group)

# Filter to genes and order correctly
tclus.genes <- tclus.genes %>%
  filter(feature %in% features_Tclus) %>%
  arrange(factor(group, levels = levelsT)) %>%
  arrange(factor(feature, levels = features_Tclus))

dp1 <- tclus.genes %>% 
  mutate(`% Expressing` = pct_in) %>% 
  mutate(`Avg. Exp.` = logFC) %>%
  ggplot(aes(x=feature, y = group, color = logFC, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size = 12) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = color,  oob = scales::squish,  
                        limits = c(-2,2), name = "Avg. Exp.") + 
  scale_y_discrete(position = "left") +
  theme(legend.title = element_text(size = 11))
dp1$data$feature <- factor(dp1$data$feature, levels = features_Tclus)
dp1$data$group <- factor(dp1$data$group, levels = rev(levelsT))
dp1 

# *Fig. S2- correlation heatmaps -----
colorblind <- colorRampPalette(rev(c("#0D0887FF","#BD3786FF", "#F0F921FF")))(15)

head(integrated_tcr@meta.data)
Idents(integrated_tcr) <- "Cell_Type"
table(Idents(integrated_tcr))

counttable.Sample_ID <- table(integrated_tcr$Sample_ID, Idents(integrated_tcr))

counttable.mat.Sample_ID <- as.data.frame.matrix(counttable.Sample_ID)

counttable <- counttable.mat.Sample_ID ### CHANGE ME!!!

prop <- counttable/rowSums(counttable) * 100 
prop.df <- as.data.frame.matrix(prop)

# Correlations 
cor.df <- cor(t(prop.df))
cor_results <- rcorr(as.matrix(t(prop.df)))
# Extract the correlation matrix
cor_matrix <- cor_results$r
# Extract the p-value matrix
p_matrix <- cor_results$P

# Back to props 
row_sums <- rowSums(prop.df) # = 100??

write.xlsx(prop.df, "Cell_Type_Props.xlsx", sheetName = "By_sample_type",
           col.names = TRUE, row.names = TRUE, append = TRUE)

prop.df$Sample_Type <- rownames(prop.df)
prop.df$Sample_Type <- c("Liver","Liver","Liver", "Pancreas", "PBMC", "Liver","Liver","Liver", 
                         "Pancreas", "PBMC", "Liver","Liver","Liver", "Pancreas", "PBMC")
prop.df$Pt_Sample <- c("PT1_LiverA","PT1_LiverB","PT1_LiverC", "PT1_Pancreas", "PT1_PBMC", 
                       "PT2_LiverA","PT2_LiverB","PT2_LiverC", "PT2_Pancreas", "PT2_PBMC", 
                       "PT3_LiverA","PT3_LiverB","PT3_LiverC", "PT3_Pancreas", "PT3_PBMC")

prop.df$Pt_ID <- c("PT1","PT1","PT1", "PT1", "PT1", 
                   "PT2","PT2","PT2", "PT2", "PT2", 
                   "PT3","PT3","PT3", "PT3", "PT3")

melt.df <- reshape2::melt(prop.df)
melt.df <- melt.df %>%
  rename(variable = "Cell_Type") %>%
  rename(value = "% Cell_Type")

### Filter to just liver and make a heatmap for each patient -----
liv.df <- prop.df %>%
  filter(Sample_Type == "Liver")
rownames(liv.df) <- liv.df$Pt_Sample

pt1.df <- liv.df %>%
  filter(Pt_ID == "PT1")
pt2.df <- liv.df %>%
  filter(Pt_ID == "PT2")
pt3.df <- liv.df %>%
  filter(Pt_ID == "PT3")

breaksList1 = seq(0, 20, by = 5)
ComplexHeatmap::pheatmap(t(pt1.df[1:32]), # Change for each patient
                         display_numbers = T, number_color = "grey67", 
                         row_names_side = "left",
                         fontsize_row = 11, fontsize_col = 11,
                         border_color = NA,
                         fontsize = 15, cluster_cols = F, #main = "Cell Type Proportions", 
                         color = rev(colorblind),
                         cluster_rows = F, breaks = breaksList1, 
                         heatmap_legend_param = list(
                           title = "% Cell \nProportion\n"))

### Jaccard overlap per patient -----
# 
subset1 <- subsetClones(combined, name = "patient_id", variables = c("PT1"))
clonalOverlap(subset1, 
              cloneCall = seq_choice,
              chain = chain_choice,
              method = "jaccard",
              group.by = "sample",
              palette = "Plasma") + 
  theme(axis.text.x = element_text(angle = 90)) 

subset2 <- subsetClones(combined, name = "patient_id", variables = c("PT2"))
clonalOverlap(subset2, 
              cloneCall = seq_choice,
              chain = chain_choice,
              method = "jaccard",
              group.by = "sample",
              palette = "Plasma") + 
  theme(axis.text.x = element_text(angle = 90)) 

subset3 <- subsetClones(combined, name = "patient_id", variables = c("PT3"))
clonalOverlap(subset3, 
              cloneCall = seq_choice,
              chain = chain_choice,
              method = "jaccard",
              group.by = "sample",
              palette = "Plasma") + 
  theme(axis.text.x = element_text(angle = 90)) 


# *Fig. S3- site split volcanos -----
# effmemCD8, cytoCD8, Tregs

### CD8 eff mem -----

CD8.eff.markers.PL <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "Liver",
                                  group.by = "sample.type", 
                                  subset.ident = "Effector memory CD8 T cells", ### CHANGE ME!!!
                                  verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.eff.markers.PL$gene <- rownames(CD8.eff.markers.PL)
# Here we use log2FC * adj pval for ranks
#CD8.eff.markers$rank <- CD8.eff.markers$avg_log2FC * CD8.eff.markers$p_val_adj

CD8.eff.markers.PPB <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "PBMC",
                                   group.by = "sample.type", 
                                   subset.ident = "Effector memory CD8 T cells", ### CHANGE ME!!!
                                   verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.eff.markers.PPB$gene <- rownames(CD8.eff.markers.PPB)


CD8.eff.markers.LPB <- FindMarkers(integrated_tcr, ident.1 = "Liver", ident.2 = "PBMC",
                                   group.by = "sample.type", 
                                   subset.ident = "Effector memory CD8 T cells", ### CHANGE ME!!!
                                   verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.eff.markers.LPB$gene <- rownames(CD8.eff.markers.LPB)


### CD8 cyto -----

CD8.cyto.markers.PL <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "Liver",
                                   group.by = "sample.type", 
                                   subset.ident = "Cytotoxic CD8 T cells", ### CHANGE ME!!!
                                   verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.cyto.markers.PL$gene <- rownames(CD8.cyto.markers.PL)
# Here we use log2FC * adj pval for ranks
#CD8.cyto.markers$rank <- CD8.cyto.markers$avg_log2FC * CD8.cyto.markers$p_val_adj

CD8.cyto.markers.PPB <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "PBMC",
                                    group.by = "sample.type", 
                                    subset.ident = "Cytotoxic CD8 T cells", ### CHANGE ME!!!
                                    verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.cyto.markers.PPB$gene <- rownames(CD8.cyto.markers.PPB)


CD8.cyto.markers.LPB <- FindMarkers(integrated_tcr, ident.1 = "Liver", ident.2 = "PBMC",
                                    group.by = "sample.type", 
                                    subset.ident = "Cytotoxic CD8 T cells", ### CHANGE ME!!!
                                    verbose = F, recorrect_umi = FALSE, assay = "SCT")

CD8.cyto.markers.LPB$gene <- rownames(CD8.cyto.markers.LPB)

### Tregs -----

Idents(integrated_tcr) <- "Cell_Type"
treg.markers.PL <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "Liver",
                               group.by = "sample.type", 
                               subset.ident = "Tregs", ### CHANGE ME!!!
                               verbose = F, recorrect_umi = FALSE, assay = "SCT")

treg.markers.PL$gene <- rownames(treg.markers.PL)
# Here we use log2FC * adj pval for ranks
#treg.markers$rank <- treg.markers$avg_log2FC * treg.markers$p_val_adj

treg.markers.PPB <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "PBMC",
                                group.by = "sample.type", 
                                subset.ident = "Tregs", ### CHANGE ME!!!
                                verbose = F, recorrect_umi = FALSE, assay = "SCT")

treg.markers.PPB$gene <- rownames(treg.markers.PPB)


treg.markers.LPB <- FindMarkers(integrated_tcr, ident.1 = "Liver", ident.2 = "PBMC",
                                group.by = "sample.type", 
                                subset.ident = "Tregs", ### CHANGE ME!!!
                                verbose = F, recorrect_umi = FALSE, assay = "SCT")

treg.markers.LPB$gene <- rownames(treg.markers.LPB)

# Volcano plots 3 comparisons -----

fc <- 0.5
pval <- 0.05

markers1 <- treg.markers.PL %>% ### CHANGE ME!!!
  filter(!(gene %in% BC.genes))

degs1 <-  markers1 %>%
  filter(abs(avg_log2FC) > fc) %>% 
  filter(p_val_adj <= pval) %>%
  filter(!(gene %in% BC.genes))

vp1 <- EnhancedVolcano(markers1,
                       lab = rownames(markers1),
                       selectLab = rownames(degs1),
                       subtitle = bquote(italic("Liver <-> Pancreas")),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = "Tregs", ### CHANGE ME!!!
                       titleLabSize = 17,
                       subtitleLabSize = 13,
                       captionLabSize = 13,
                       pCutoff = pval,
                       FCcutoff = fc,
                       legendLabels = NULL,
                       legendIconSize = 0,
                       pointSize = 2.5,
                       labSize = 3,
                       axisLabSize = 13,
                       legendLabSize = 13, 
                       colAlpha = 0.35)

markers2 <- treg.markers.PPB %>% ### CHANGE ME!!!
  filter(!(gene %in% BC.genes))
degs2 <- markers2 %>%
  filter(abs(avg_log2FC) > fc) %>% 
  filter(p_val_adj <= pval) 

vp2 <- EnhancedVolcano(markers2,
                       lab = rownames(markers2),
                       selectLab = rownames(degs2),
                       subtitle = bquote(italic("PBMC <-> Pancreas")),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = "", 
                       titleLabSize = 17,
                       subtitleLabSize = 13,
                       captionLabSize = 13,
                       pCutoff = pval,
                       FCcutoff = fc,
                       legendLabels = c("NS", paste("Log2FC > ", fc),
                                        paste("Adj. pval < ", pval), 
                                        paste("Log2FC > ", fc, "and Adj. pval < ", pval)),
                       #legendIconSize = 0,
                       pointSize = 2.5,
                       labSize = 3,
                       axisLabSize = 13,
                       legendLabSize = 13, 
                       colAlpha = 0.35)

markers3 <- treg.markers.LPB %>% ### CHANGE ME!!!
  filter(!(gene %in% BC.genes))
degs3 <- markers3 %>%
  filter(abs(avg_log2FC) > fc) %>% 
  filter(p_val_adj <= pval)

vp3 <- EnhancedVolcano(markers3,
                       lab = rownames(markers3),
                       selectLab = rownames(degs3),
                       subtitle = bquote(italic("PBMC <-> Liver")),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = "", 
                       titleLabSize = 18,
                       subtitleLabSize = 13,
                       captionLabSize = 13,
                       pCutoff = pval,
                       FCcutoff = fc,
                       legendLabels = NULL,
                       legendIconSize = 0,
                       pointSize = 2.5,
                       labSize = 3,
                       axisLabSize = 13,
                       legendLabSize = 13, 
                       colAlpha = 0.35)
vp <- vp1+vp2+vp3

# Fig. S4- diversity- scRepertoire -----

### Rarefaction Shannon -----
clonalRarefaction(combined,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2,
                  palette = "Zissou 1")
### Diversity- all metrics
clonalDiversity(combined, 
                cloneCall = "aa", 
                group.by = "sample",
                palette = "Roma")

### Disctibution
clonalSizeDistribution(combined, 
                       cloneCall = "aa", 
                       method= "ward.D2",
                       palette = "Roma")
