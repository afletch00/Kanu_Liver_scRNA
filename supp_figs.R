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

# *Fig. S3- site split volcanos -----
# effmemCD8, cytoCD8, Tregs



# Fig. S4- D50 and rarefaction scRepertoire -----