# integrated_tcr is the main Seurat object used for all analysis

# Functions -----
int_drc <- function(integrated){
  
  DefaultAssay(integrated) <- "integrated" # 
  integrated <- ScaleData(integrated, verbose = T) 
  integrated <- RunPCA(object = integrated, verbose = T)
  integrated <- RunUMAP(object = integrated, dims = 1:30, verbose = T)
  integrated <- FindNeighbors(integrated, dims = 1:30, verbose = T)
  integrated <- FindClusters(integrated, resolution = 0.7, 
                             random.seed = 0,
                             verbose = T, method = "igraph", # for large datasets
                             algorithm = 4) # use when working with large object
  # can set multiple resolutions to compare
  return(integrated)
}

assay_drc <- function(sub, assay, res){
    
  DefaultAssay(sub) <- assay
  if (assay == "RNA") {
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 3000)
    sub <- ScaleData(sub, vars.to.regress = "MT") 
    sub <- NormalizeData(sub)
    sub <- RunPCA(object = sub, verbose = TRUE)
    sub <- RunUMAP(object = sub, dims = 1:30, verbose = TRUE, 
                   reduction.name = "umap_unintegrated")
    sub <- FindNeighbors(sub, dims = 1:30, reduction = "pca", verbose = TRUE)
    sub <- FindClusters(sub, resolution = res, verbose = TRUE, algorithm = 4,
                        cluster.name = "unintigrated_clusters")
  } else {
    sub <- ScaleData(sub)
    sub <- RunPCA(object = sub, verbose = TRUE)
    sub <- RunUMAP(object = sub, dims = 1:30, verbose = TRUE)
    sub <- FindNeighbors(sub, dims = 1:30, verbose = TRUE)
    sub <- FindClusters(sub, resolution = res, algorithm = 4, verbose = TRUE)
  }
  return(sub)
}
  
## *Fig. 1B -----

### Set cell order -----
cell_order <- c(# T cells
  "Naïve/central memory CD4 T cells", "Tissue resident memory CD4 T cells",
  "Effector memory CD4 T cells", "Tregs", "Naïve/central memory CD8 T cells", 
  "Tissue resident memory CD8 T cells", "Effector memory CD8 T cells",
  "Cytotoxic CD8 T cells","Follicular CD8 T cells", "IFNg+ CD8 T cells",
  "Proliferating/Cycling T cells","MAIT","Gamma delta T cells", 
  # NK cells
  "Tissue resident NK cells", "Circulating NK cells",
  "CD56BRIGHT transitional NK cells", "Proliferating/Cycling NK cells",
  # Myeloid cells
  "CD14+ blood monocytes", "CD14+ tissue monocytes", "CD16+ blood monocytes", "CD16+ tissue monocytes",
  "Macrophages","Conventional DCs", "Plasmacytoid DCs", "Neutrophils",
  # B cells
  "Naïve B cells", "ISG+ B cells", "CXCR5+ B cells","Circulating memory B cells", "Tissue memory B cells",  
  "Fcrl5+ B cells", "Plasma cells")

integrated_tcr@meta.data$Cell_Type <- factor(integrated_tcr@meta.data$Cell_Type, levels = cell_order)
levels(integrated_tcr@meta.data$Cell_Type)
table((integrated_tcr@meta.data$Cell_Type))

### Dim plots -----
DimPlot(integrated_tcr, group.by = "Cell_Type",label = FALSE, shuffle = TRUE, pt.size = 0.7,
        split.by = NULL) + 
  scale_color_manual(values = colors_32) +
  theme(plot.title = element_text(size = 13, face = "bold")) + 
  ggtitle("")

DimPlot(integrated_tcr, group.by = "Main_Cell_Type",label = FALSE, shuffle = TRUE, pt.size = 0.5,
        split.by = NULL) + 
  scale_color_manual(values = colors_5) +
  theme(plot.title = element_text(size = 13, face = "bold")) + 
  ggtitle("")

DimPlot(integrated_tcr, group.by = "Patient_ID",label = FALSE, shuffle = TRUE, pt.size = 0.5,
        split.by = NULL) + 
  scale_color_manual(values = colors_5) +
  theme(plot.title = element_text(size = 13, face = "bold")) +
  ggtitle("")

### Make integrated_tcr with only liver 1 -----
integrated_tcr_livA <- subset(integrated_tcr, subset = Sample_ID %in%
                                c("PT1-LiverB", "PT1-LiverC","PT2-LiverB", "PT2-LiverC_Cryo",
                                  "PT3-LiverB", "PT3-LiverC"),
                              invert = TRUE) 
integrated_tcr_livA@meta.data$Cell_Type <- factor(integrated_tcr_livA@meta.data$Cell_Type, 
                                                  levels = cell_order)
levels(integrated_tcr_livA@meta.data$Cell_Type)
# Do not re-cluster
DimPlot(integrated_tcr_livA, group.by = "sample.type",label = FALSE, shuffle = TRUE, pt.size = 0.5,
        split.by = NULL) + 
  scale_color_manual(values = colors_5) +
  theme(plot.title = element_text(size = 13, face = "bold")) + 
  ggtitle("")


## *Fig. 1C- Dotplot-----

### Characterization genes -----
genes_3 <- c(
  "IL7R", "CCR7", "THOC3",
  "CD40LG", "TNFAIP3", "CD69",
  "CISH", "ARHGAP15", "PBXIP1",
  "FOXP3", "CTLA4", "CCR4",
  "LDHB", "NOSIP", "CD248",
  "IFNG", "CCL4", "CCL5",
  "GZMK", "CD27", "KLRG1",
  "CX3CR1", "GZMH", "PRF1",
  "TNF", "DNAJB1", "CST7",
  "RORA", "DDX3X", "NKTR",
  "TRAV1-2", "RORC", "DUSP2",
  "TRDV1", "TRGV9", "TRGC2",
  "MKI67", "H3C2", "VDAC3",
  "HIST1H4F", "TMPO-AS1", "UNG",
  "IL2RB", "CD160", "XCL1",
  "GZMB", "FCGR3A", "GNLY",
  "NCAM1", "KLRC1", "KLRF1",
  "IL4R", "FCER2", "IGHD",
  "ISG20", "IFIT3", "TCL1A",
  "CXCR5", "IL6", "CD70",
  "CD82", "IGHA1", "PLP2",
  "MS4A1", "STMN1", "TPSAB1",
  "FCRL5", "CD22", "FCRL3",
  "JCHAIN", "IGKC", "IGHG1",
  "VCAN", "CD14", "CSF3R",
  "THBS1", "EREG", "OSM",
  "S100A8", "CXCL8", "IL1B",
  "CDKN1C", "HMOX1", "C5AR2",
  "CSF1R", "CD68", "TNFRSF8",
  "APOE", "C1QA", "AXL",
  "CLEC9A", "IDO1", "HLA-DQB2",
  "LILRA4", "IRF7", "CLEC4C"
)

Idents(integrated_tcr) <- "Cell_Type"
levels(Idents(integrated_tcr))

### Set order for dotplots -----
cell_order2 <- c(# T cells
  "Naïve/central memory CD4 T cells", "Tissue resident memory CD4 T cells",
  "Effector memory CD4 T cells", "Tregs", "Naïve/central memory CD8 T cells", 
  "Tissue resident memory CD8 T cells", "Effector memory CD8 T cells",
  "Cytotoxic CD8 T cells","Follicular CD8 T cells", "IFNg+ CD8 T cells",
  "MAIT","Gamma delta T cells", "Proliferating/Cycling T cells",
  # NK cells
  "Proliferating/Cycling NK cells", "Tissue resident NK cells", "Circulating NK cells",
  "CD56BRIGHT transitional NK cells", 
  # B cells
  "Naïve B cells", "ISG+ B cells", "CXCR5+ B cells","Circulating Memory B cells", "Tissue memory B cells",  
  "Fcrl5+ B cells", "Plasma cells",
  # Myeloid cells
  "CD14+ blood monocytes", "CD14+ tissue monocytes", "Neutrophils","CD16+ blood monocytes", 
  "CD16+ tissue monocytes","Macrophages","Conventional DCs", "Plasmacytoid DCs")

integrated_tcr@meta.data$Cell_Type_2 <- integrated_tcr@meta.data$Cell_Type
table(integrated_tcr$Cell_Type_2)
integrated_tcr@meta.data$Cell_Type_2 <- factor(integrated_tcr@meta.data$Cell_Type_2, levels = cell_order2)
levels(integrated_tcr@meta.data$Cell_Type_2)

# Make sure data was normalized. If this does not return integers then yes. 
integrated_tcr[['RNA']]@data@x

# log-normalized versions of corrected counts are in @data
Idents(integrated_tcr) <- "Cell_Type_2"
All.genes3 <- wilcoxauc(integrated_tcr, seurat_assay = "RNA", assay = "data")
head(All.genes3)
table(All.genes3$group)

saveRDS(All.genes3, "wilcoxon_stats_RNA_for_dotplot.rds")

### Filter to genes and order correctly -----
plot.genes3 <- All.genes3 %>%
  filter(feature %in% genes_3) %>%
  arrange(factor(group, levels = cell_order2)) %>%
  arrange(factor(feature, levels = genes_3))
#arrange(dplyr::desc(logFC), dplyr::desc(auc)) %>%
head(plot.genes3, n = 10)
table(plot.genes3$feature)

### Color palette -----
RColorBrewer::brewer.pal.info

display.brewer.all(colorblindFriendly = TRUE)   
color_palette <- brewer.pal(11, "RdYlBu")
color <- colorRampPalette((c("steelblue", "white","firebrick")))(50)
color2 <- colorRampPalette((c("white", "purple4")))(50)
color_palette2 <- brewer.pal(9, "Blues")

color3 <- colorRampPalette((c("grey97","mediumpurple1", "blue3")))(50)

### Manual dotplot -----

dp1 <- plot.genes3 %>% 
  mutate(`% Expressing` = pct_in) %>% 
  mutate(`LogFC Exp.` = logFC) %>%
  #filter(avgExpr > 0, `% Expressing` > 0) %>%
  ggplot(aes(x=feature, y = group, color = `LogFC Exp.`, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size = 12) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  #scale_fill_gradientn(colors = c("red"), name = 'Average Expression')
  scale_color_gradientn(colours = color,  oob = scales::squish,  #viridis::magma(50), 
                        limits = c(-2,2), name = "LogFC Exp.") + 
  scale_y_discrete(position = "right") +
  theme(legend.title = element_text(size = 11))
dp1$data$feature <- factor(dp1$data$feature, levels = genes_3)
dp1$data$group <- factor(dp1$data$group, levels = rev(cell_order2))
dp1

### Expanded TCR stacked barplot -----
# 07122023 stacked barplot of unique clones and cell types
# 
# #*Adding "Unique" as a level to the factor variable
TCRs@meta.data$cloneType <- as.factor(TCRs@meta.data$cloneType)
TCRs@meta.data$cloneType <- factor(TCRs@meta.data$cloneType, levels = c(levels(TCRs@meta.data$cloneType), "Unique"))

for (i in 1:nrow(TCRs@meta.data)) {
  if (TCRs@meta.data$CTnt[i] %in% ct_table$Var1 & 
      TCRs@meta.data$Sample_ID[i] %in% ct_table$Var2) {
    TCRs@meta.data$cloneType[i] <- "Unique"
  }
}

TCRmetaCloneNum <- data.frame(TCRs@meta.data)

sub_meta1 <- TCRmetaCloneNum[TCRmetaCloneNum$cloneType != "Unique", ]

sub_table <- table(sub_meta1$Cell_Type, sub_meta1$Sample_ID)
sub_table <- data.frame(sub_table)
sub_table <- sub_table[order(sub_table$Freq, decreasing = FALSE),]

sub_table$Cell_Type <- factor(sub_table$Var1)
sub_table$Sample_ID <- factor(sub_table$Var2)
sub_table$Var1 <- NULL
sub_table$Var2 <- NULL

sub_table <- sub_table %>%
  filter(Cell_Type != "IFNg+ CD8 T cell")

sub_table2 <- sub_table %>%
  pivot_wider(names_from = Cell_Type, values_from = Freq) 
rownames <- sub_table2$Sample_ID
sub_table2$Sample_ID <- NULL
rownames(sub_table2) <- rownames
prop <- sub_table2/rowSums(sub_table2) * 100
openxlsx::write.xlsx(prop, file = "expanded_clones_props.xlsx", asTable = TRUE,
                     colNames = TRUE, rowNames = TRUE)
#selected_rows <- sub_table[grep("PT3-PBMC", df_submeta1$Sample_ID), ]

ggp1 <- ggplot(sub_table, aes(x = Sample_ID, y = Freq, fill = Cell_Type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors_32) + 
  labs(x = "", y = "Expanded Clones") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = "Cell Type"))
ggp1  

# *Fig. 2A- Dimplot-----

DimPlot(integrated_tcr, group.by = "Cell_Type",label = FALSE, shuffle = TRUE, pt.size = 0.7,
        split.by = "sample.type") + 
  scale_color_manual(values = colors_32) +
  theme(plot.title = element_text(size = 13, face = "bold")) + 
  ggtitle("")


# *Fig. 2C- GSEA -----

m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = c("GO:BP")) 
m_df$gs_name <- gsub("GOBP_", "", m_df$gs_name)
head(m_df)

# Subset to pathways of interest

select <- c("T_CELL_ACTIVATION", "T_CELL_DIFFERENTIATION", "T_CELL_MEDIATED_CYTOTOXICITY",
            "T_CELL_MIGRATION", "T_CELL_PROLIFERATION", "T_CELL_SELECTION")
m_df <- m_df %>%
  filter(gs_name %in% select)
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

## Run fgsea

### Use DGE results log2FC to rank genes -----
markers <-  markers
gs.genes <- markers %>%
  na.omit(markers) %>%
  dplyr::arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(gs.genes)
head(ranks)

fgseaRes <- fgsea(fgsea_sets, stats = ranks, minSize  = 5, maxSize  = 500, eps = 0)


# *Fig. 2D- Volcano -----

# DEG analysis
Idents(integrated_tcr) <- "Cell_Type"
mac.markers <- FindMarkers(integrated_tcr, ident.1 = "Pancreas", ident.2 = "Liver",
                           group.by = "sample.type", 
                           subset.ident = "Macrophages", 
                           verbose = F, recorrect_umi = FALSE, assay = "SCT")

mac.markers$gene <- rownames(mac.markers)

# Define cutoffs for gene labels
fc <- 0.5
pval <- 0.05
markers <- mac.markers 

degs <-  markers %>%
  filter(abs(avg_log2FC) > fc) %>% 
  filter(p_val_adj <= pval) 

vp1 <- EnhancedVolcano(markers,
                       lab = markers$gene,
                       selectLab = degs$gene,
                       subtitle = bquote(italic("Liver <-> Pancreas")),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = "Macrophages", ### CHANGE ME!!!
                       titleLabSize = 18,
                       subtitleLabSize = 14,
                       captionLabSize = 14,
                       pCutoff = pval,
                       FCcutoff = fc,
                       legendLabels = c("NS", expression(Log[2] ~ FC ~ '> 0.5'),
                                        "Adj. p-value < 0.05", 
                                        expression('Adj. p-value' ~ and ~ log[2] ~ FC~ '> 0.5')),
                       pointSize = 2.5,
                       labSize = 3.5,
                       axisLabSize = 13,
                       legendLabSize = 13, 
                       drawConnectors = FALSE,
                       widthConnectors = 0.5,
                       typeConnectors = "closed",
                       endsConnectors = "first",
                       lengthConnectors = unit(0.01, "npc"),
                       colConnectors = "grey10",
                       max.overlaps = 17,
                       maxoverlapsConnectors = NULL,
                       min.segment.length = 0,
                       directionConnectors = "both",
                       arrowheads = FALSE,
                       colAlpha = 0.35)
vp1

# *Fig. 2E- Mac signatures-----
library(readr)
Mac_Signature_Scores <- read_csv("/hpc/group/peterallen/af206/scRNA-seq/FINAL_ANALYSIS/for_paper/Mac_Signature_Scores.csv")

Idents(integrated_tcr) <- "Cell_Type"
Macs <- subset(integrated_tcr, ident = "Macrophages")
table(Idents(Macs))

gs1m <- GeneSet(setName = "Azizi M1", setIdentifier="101", 
                unique(Mac_Signature_Scores$"Azizi_M1"))
geneIds <- geneIds(gs1m) 
is(geneIdType(gs1m), "NullIdentifier") == TRUE
geneIdType(gs1m)
geneIdType(gs1m) <- SymbolIdentifier()

gs2m <- GeneSet(setName = "Azizi M2", setIdentifier="102", 
                unique(Mac_Signature_Scores$Azizi_M2))
geneIds <- geneIds(gs2m) 
is(geneIdType(gs2m), "NullIdentifier") == TRUE
geneIdType(gs2m)
geneIdType(gs2m) <- SymbolIdentifier()

try(GeneSetCollection(gs1m, gs2m))
gscM <- GeneSetCollection(gs1m, gs2m)

# ssGSEA on individual cells

subsetMP <- subset(Macs, subset = sample.type %in% "Pancreas")
table(Idents(subsetMP))
gene_exp1 <- GetAssayData(object = subsetMP, assay = "RNA") 
gene_exp1 <- as.matrix(gene_exp1)

subsetML <- subset(Macs, subset = sample.type %in% "Liver")
table(Idents(subsetML))
gene_exp2 <- GetAssayData(object = subsetML, assay = "RNA") 
gene_exp2 <- as.matrix(gene_exp2)

# GSVA needs a gene expression matrix
# Run all samples with both the upregulated and downregulated gene sets
gbmPar1 <- ssgseaParam(gene_exp1, gscM)
#gsvaParam(gene_exp, gsc, method = "ssgsea")
gbm_es1 <- gsva(gbmPar1)

gsva.es.df1 <- data.frame(gbm_es1)
#write.xlsx(gsva.es, "ssGSVA_results2.xlsx")
gsva.es.df1 <- cbind(rownames(gsva.es.df1), gsva.es.df1)
names(gsva.es.df1)[1] <- "Gene_list"
rownames(gsva.es.df1) <- NULL

# do for each group
m.dfP <- reshape2::melt(gsva.es.df1)
m.dfP$Sample_Type <- "Pancreas"
#
#
gbmPar2 <- ssgseaParam(gene_exp2, gscM)
#gsvaParam(gene_exp, gsc, method = "ssgsea")
gbm_es2 <- gsva(gbmPar2)

gsva.es.df2 <- data.frame(gbm_es2)
#write.xlsx(gsva.es, "ssGSVA_results2.xlsx")
gsva.es.df2 <- cbind(rownames(gsva.es.df2), gsva.es.df2)
names(gsva.es.df2)[1] <- "Gene_list"
rownames(gsva.es.df2) <- NULL
m.dfL <- reshape2::melt(gsva.es.df2)
m.dfL$Sample_Type <- "Liver"

# bind 
m.df.allM <- rbind(m.dfL, m.dfP)

# create the ggplot using the merged data
gp <- ggplot(m.df.allM, aes(x = Sample_Type, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(y = value), width = 0.2, height = 0, alpha = 0.3, 
              color = "blue3", size = 0.3) + 
  facet_wrap(~ Gene_list, scales = "free", ncol = 5) +
  labs(x = "", y = "Gene Set Variation Score") +
  #ggtitle("ssGSEA") +
  stat_compare_means(label = "p.signif") + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))  
gp # 8.3x4.3

# *Fig. 3A -----
# https://www.borch.dev/uploads/screpertoire/
# Palettes - https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
# 
### Remove the TCR genes- use this UMAP -----
curr.count <- integrated_tcr@assays$SCT@counts
rownames(curr.count)[1:10]
TCR.genes <- rownames(curr.count)[grepl("^TRAV|^TRAJ|^TRBV|^TRBJ", rownames(curr.count))]
### subset to T cells -----
Tclus <- subset(integrated_tcr, ident = c("Naïve/central memory CD4 T cells", "Tissue resident memory CD4 T cells",
                                           "Effector memory CD4 T cells", "Tregs", "Naïve/central memory CD8 T cells", 
                                           "Tissue resident memory CD8 T cells", "Effector memory CD8 T cells",
                                           "Cytotoxic CD8 T cells","Follicular CD8 T cells", "IFNg+ CD8 T cells",
                                           "Proliferating/Cycling T cells","MAIT","Gamma delta T cells"))

Tclus@meta.data$Cell_Type <- droplevels(Tclus@meta.data$Cell_Type)
levels(Tclus@meta.data$Cell_Type)
table(Idents(Tclus))

DefaultAssay(Tclus) <- "SCT"
counts <- GetAssayData(Tclus, assay = "SCT")
counts2 <- counts[-(which(rownames(counts) %in% TCR.genes)),]
genes <- rownames(counts2)
rownames(counts2)[grepl("^TRAV", rownames(counts2))]

Tclus.noTCRgenes <- subset(Tclus, features = rownames(counts2))
curr.count <- Tclus.noTCRgenes@assays$SCT@counts
rownames(curr.count)[grepl("^TRB", rownames(curr.count))]

### diminsion reduction and cluster -----
DefaultAssay(Tclus.noTCRgenes) <- "integrated"
Tclus.noTCRgenes <- RunPCA(object = Tclus.noTCRgenes, npcs = 20, verbose = TRUE)
Tclus.noTCRgenes <- RunUMAP(object = Tclus.noTCRgenes, dims = 1:20, verbose = TRUE)
Tclus.noTCRgenes <- FindNeighbors(Tclus.noTCRgenes, dims = 1:20, reduction = "pca", verbose = TRUE)
Tclus.noTCRgenes <- FindClusters(Tclus.noTCRgenes, resolution = 0.7, verbose = TRUE)

### scRepertoire -----

# Get file names -----
# https://ncborcherding.github.io/vignettes/vignette.html
tcr.outs <- "/hpc/group/peterallen/af206/scRNA-seq/contig_files/All"
sc_tcr_dirs <- list.files(path = tcr.outs, full.names = T)
sc_tcr_dirs

sc_tcrs <- sc_tcr_dirs[c(2:16)]

# Make metadata file -----
meta_data <- data.frame(path = sc_tcrs) %>%
  dplyr::mutate(sample_fullname = basename(path),
                sample_id = gsub(".csv","",sample_fullname)) %>%
  separate(
    col = sample_id,
    into = c("patient_id", "sample_type"),
    sep = "_",
    remove = F
  ) %>%
  dplyr::mutate(match_id = paste0(patient_id, ":", sample_type))

meta_data$sample_type <- ifelse(
  grepl("^Liver[ABC]$", meta_data$sample_type),
  sub("[ABC]$", "", meta_data$sample_type),
  meta_data$sample_type
)

filtered_contig_t_list <- lapply(meta_data$path,read_csv)
names(filtered_contig_t_list) <- meta_data$match_id

# Get combined TCR object -----
combined <- combineTCR(filtered_contig_t_list,
                       sample = meta_data$match_id,
                       removeNA = FALSE, 
                       removeMulti = FALSE, 
                       filterMulti = FALSE)
# add some metadata variables
combined <- addVariable(combined, 
                        variable.name = "sample_type", 
                        variables = meta_data$sample_type)

combined <- addVariable(combined, 
                        variable.name = "patient_id", 
                        variables = meta_data$patient_id)


# Chain and sequence setting -----
chain_choice <- "TRB"
seq_choice <- "aa"

# FUnction to ilter to non-expanded clones
filter_nonexpand <- function(contiglist){
  tab <- table(contiglist$cdr3_nt)
  tab.df <- as.data.frame(tab)
  remove <- tab.df %>%
    filter(Freq > 1)
  contiglist <- contiglist %>%
    filter(cdr3_nt %in% remove$Var1)
  return(contiglist)
}

contig_list_filt <- lapply(filtered_contig_t_list, filter_nonexpand)

combined <- combineTCR(contig_list_nonExp,
                       removeNA = T, 
                       removeMulti = F)

### Match barcodes
# Define a regular expression pattern to match and remove prefixes
prefix_pattern <- "^(L1_Pt1_|L2_Pt1_|L3_Pt1_|panc_Pt1_|pbmc_Pt1_|L1_Pt2_|L2_Pt2_|L3_Pt2_|panc_Pt2_|pbmc_Pt2_|L1_Pt3_|L2_Pt3_|L3cryo_Pt3_|panc_Pt3_|pbmc_Pt3_|L1cryo_Pt3_|L1NoDigest_Pt3_)"
# Strip the prefixes from the barcodes
stripped_barcodes <- sub(prefix_pattern, "", colnames(TCRs))
new_barcodes <- paste0(TCRs$Sample_ID, "_", stripped_barcodes) 
new_barcodes <- sub("-", "_", new_barcodes)
# Update the barcodes in the Seurat object
TCRs$barcode <- new_barcodes
# Update the barcode names in the Seurat object
colnames(TCRs) <- TCRs$barcode
head(TCRs@meta.data)

### Clonetype -----
slot(TCRs, "meta.data")$cloneSize <- NULL
TCRs <- combineExpression(combined, TCRs, ### OR combined3
                          cloneCall="nt", 
                          group.by = NULL, # ID Calculates across all samples for each patient
                          #leave out to calculate by sample
                          proportion = TRUE, addLabel = FALSE,
                          cloneSize=c(Rare = 0.009, Small = 0.01, Medium = 0.02, 
                                       Large = 0.05, Hyperexpanded = 1))
Idents(TCRs) <- "cloneSize"
table(Idents(TCRs))
slot(TCRs, "meta.data")$cloneSize <- factor(slot(TCRs, "meta.data")$cloneSize, 
                                            levels = c("Hyperexpanded (0.05 < X <= 1)", 
                                                       "Large (0.02 < X <= 0.05)", 
                                                       "Medium (0.01 < X <= 0.02)", 
                                                       "Small (0.009 < X <= 0.01)", 
                                                       "Rare (0 < X <= 0.009)", NA))


#*Adding "Unique" as a level to the factor variable
TCRs@meta.data$cloneSize <- factor(TCRs@meta.data$cloneSize, levels = c(levels(TCRs@meta.data$cloneSize), "Unique"))
TCRs@meta.data$cloneType <- TCRs@meta.data$cloneSize

TCRs@meta.data$cloneType[TCRs@meta.data$clonalFrequency == 1] <- "Unique"

TCRmeta <- data.frame(TCRs@meta.data)

### Add TCR+ column -----
integrated_tcr@meta.data$TCR <- NA
Idents(integrated_tcr) <- "TCR"
table(Idents(integrated_tcr))
integrated_tcr@meta.data$TCR[is.na(integrated_tcr@meta.data$cloneType)] <- "No"
integrated_tcr@meta.data$TCR[!is.na(integrated_tcr@meta.data$cloneType)] <- "Yes"

### Subset to TCR+ cells -----
Idents(Tclus.noTCRgenes) <- "TCR"
table(Idents(Tclus.noTCRgenes))
TCRs <- subset(Tclus.noTCRgenes, ident = "Yes")

### Dimplot -----
DimPlot(TCRs, group.by = "Cell_Type", label = FALSE, shuffle = TRUE, pt.size = 0.5) +
  scale_color_manual(values = colors_32) +
  theme(plot.title = element_text(size = 13, face = "bold")) + 
  ggtitle("")

### Clone dimplot -----
DimPlot(TCRs, group.by = "cloneType", split.by = "Sample_ID", ncol = 5, pt.size = 0.7,
        alpha = 0.5, shuffle = TRUE) +
  scale_color_manual(values = rev(clone.col)) + 
  theme(plot.title = element_blank())

# *Fig. 3B -----
sub_meta1 <- TCRmeta[TCRmeta$cloneType != "Unique", ]

sub_table <- table(sub_meta1$Cell_Type, sub_meta1$Sample_ID)
sub_table <- data.frame(sub_table)
sub_table <- sub_table[order(sub_table$Freq, decreasing = FALSE),]

sub_table$Cell_Type <- factor(sub_table$Var1)
sub_table$Sample_ID <- factor(sub_table$Var2)
sub_table$Var1 <- NULL
sub_table$Var2 <- NULL

sub_table <- sub_table %>%
  filter(Cell_Type != "IFNg+ CD8 T cell")

ggp1 <- ggplot(sub_table, aes(x = Sample_ID, y = Freq, fill = Cell_Type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors_32) + 
  labs(x = "", y = "Expanded Clones") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = "Cell Type"))
ggp1  

# *Fig. 3C -----
### Just use Liver A from each patient
liverAs <- combined[c(1,4,5,6,9,10,11,14,15)]

clonalQuant(liverAs, cloneCall="nt", scale = T, exportTable = F, chain = "TRA") +# , exportTable = T
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('grey','grey','grey','grey','grey','grey','grey','grey','grey'))

# *Fig. 3D -----
### Top 20 PDAC contigs- TRB chain and aa 
### pt1 -----
Panc.1_counts <- liverAs[[2]]
Panc.1_counts2 <- table(Panc.1_counts$CTnt)
Panc.1_counts_df <- data.frame(Panc.1_counts2)
# Get all clones >1 in the pancreas
Panc.1_counts_df2 <- subset(Panc.1_counts_df, Panc.1_counts_df$Freq > 1)
sorted_counts <- sort(Panc.1_counts2, decreasing = TRUE)
top_20_1 <- names(sorted_counts)[1:20]
top_20_1


# amino acid sequences
Panc.1_aa <- table(Panc.1_counts$CTaa)
Panc.1_aa_df <- data.frame(Panc.1_aa)
sorted_aa1 <- sort(Panc.1_aa, decreasing = TRUE)
Panc.1_aa_df <- data.frame(sorted_aa1)
top_20aa_1 <- names(sorted_aa1)[1:20]
top_20aa_1
top_20aa_1.df <- as.data.frame(top_20aa_1)

### pt2 -----
Panc.2_counts <- liverAs[[5]]
Panc.2_counts2 <- table(Panc.2_counts$CTnt)
Panc.2_counts_df <- data.frame(Panc.2_counts2)
Panc.2_counts_df2 <- subset(Panc.2_counts_df, Panc.2_counts_df$Freq > 1)
sorted_counts2 <- sort(Panc.2_counts2, decreasing = TRUE)
top_20_2 <- names(sorted_counts2)[1:20]
top_20_2

# AA
Panc.2_aa <- table(Panc.2_counts$CTaa)
Panc.2_aa_df <- data.frame(Panc.2_aa)
sorted_aa2 <- sort(Panc.2_aa, decreasing = TRUE)
Panc.2_aa_df <- data.frame(sorted_aa2)
top_20aa_2 <- names(sorted_aa2)[1:20]
top_20aa_2
top_20aa_2.df <- as.data.frame(top_20aa_2)

### pt3 -----
Panc.3_counts <- liverAs[[8]]
Panc.3_counts2 <- table(Panc.3_counts$CTnt)
Panc.3_counts_df <- data.frame(Panc.3_counts2)
Panc.3_counts_df2 <- subset(Panc.3_counts_df, Panc.3_counts_df$Freq > 1)
sorted_counts3 <- sort(Panc.3_counts2, decreasing = TRUE)
top_20_3 <- names(sorted_counts3)[1:20]
top_20_3

# aa sequences for annotation
Panc.3_aa <- table(Panc.3_counts$CTaa) 
Panc.3_aa_df <- data.frame(Panc.3_aa)
sorted_aa3 <- sort(Panc.3_aa, decreasing = TRUE)
Panc.3_aa_df <- data.frame(sorted_aa3)
top_20aa_3 <- names(sorted_aa3)[1:20]
top_20aa_3
top_20aa_3.df <- as.data.frame(top_20aa_3)


## Match top 20 per patient to liver and blood TCRs -----
pt1 <- rbind(liverAs[[1]], liverAs[[2]], liverAs[[3]])
pt2 <- rbind(liverAs[[4]], liverAs[[5]], liverAs[[6]])
pt3 <- rbind(liverAs[[7]], liverAs[[8]], liverAs[[9]])

# Filter to top 20 panc aas -----
pt1.top20 <- pt1 %>%
  filter(CTaa %in% top_20aa_1)
pt1.top20.tab <- table(pt1.top20$sample, pt1.top20$ID, pt1.top20$CTaa)
pt1.top20.tab <- as.data.frame(pt1.top20.tab)

pt2.top20 <- pt2 %>%
  filter(CTaa %in% top_20aa_2)
pt2.top20.tab <- table(pt2.top20$sample, pt2.top20$ID, pt2.top20$CTaa)
pt2.top20.tab <- as.data.frame(pt2.top20.tab)

pt3.top20 <- pt3 %>%
  filter(CTaa %in% top_20aa_3)
pt3.top20.tab <- table(pt3.top20$sample, pt3.top20$ID, pt3.top20$CTaa)
pt3.top20.tab <- as.data.frame(pt3.top20.tab)

list.top.20 <- list("pt1.top.20.panc" = pt1.top20.tab,
                    "pt2.top.20.panc" = pt2.top20.tab,
                    "pt3.top.20.panc" = pt3.top20.tab)

openxlsx::write.xlsx(list.top.20, file = "top_20_clones_pancs_tracked.xlsx", asTable = TRUE,
                     colNames = TRUE, rowNames = TRUE)

### Plot ribbons (Sankey) -----
openxlsx::write.xlsx(liverAs, file = "data_for_sankeys.xlsx", asTable = TRUE,
                     colNames = TRUE, rowNames = TRUE)

p1 <- clonalCompare(liverAs, 
                    samples = c("PT1_LiverA", "PT1_Panc","PT1_PBMC"),
                    cloneCall="aa", 
                    clones = c(top_20aa_1),
                    graph = "alluvial") +
  scale_fill_manual(values = reds) +
  theme_grey() 
p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13),
           axis.text.y = element_text(size = 13)) +
  theme(axis.title.x=element_blank()) +
  labs(y = "Relative Abundance\n") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent)


p2 <- clonalCompare(liverAs, 
                    samples = c("PT2_LiverA", "PT2_Panc","PT2_PBMC"),
                    cloneCall="aa", 
                    clones = c(top_20aa_2),
                    graph = "alluvial") +
  scale_fill_manual(values = blues) +
  theme_grey() 
p2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13),
           axis.text.y = element_text(size = 13)) +
  theme(axis.title.x=element_blank()) +
  labs(y = "Relative Abundance\n") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent)


p3 <- clonalCompare(liverAs, 
                    samples = c("PT3_LiverA", "PT3_Panc","PT3_PBMC"),
                    cloneCall="aa", 
                    clones = c(top_20aa_3),
                    graph = "alluvial") +
  scale_fill_manual(values = purples) +
  theme_grey() 
p3 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13),
           axis.text.y = element_text(size = 13)) +
  theme(axis.title.x=element_blank()) +
  labs(y = "Relative Abundance\n") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent)

# *Fig. 4A -----
### Subset to just pancreas TCR+ cells -----
Idents(TCRs) <- "sample.type"
pancT <- subset(TCRs, ident = "Pancreas")
### Run dimensional reduction and clustering function
pancT <- int_drc(pancT)

DimPlot(pancT, group.by = "Cell_Type", split.by = NULL, ncol = 1, pt.size = 0.7,
        alpha = 0.5, shuffle = TRUE) +
  scale_color_manual(values = colors_32) + 
  theme(plot.title = element_blank())

# *Fig. 4B -----
DimPlot(pancT, group.by = "cloneType", split.by = NULL, ncol = 1, pt.size = 0.7,
        alpha = 0.5, shuffle = TRUE) +
  scale_color_manual(values = rev(clone.col5)) + 
  theme(plot.title = element_blank())

# *Fig. 4C
### Make broad cell types for T cells
table(pancT@meta.data$Cell_Type)
pancT@meta.data$Cell_Type <- droplevels(pancT@meta.data$Cell_Type)

pancT@meta.data$Cell_Type_Broad <- NA
Idents(pancT) <- "Cell_Type_Broad"
table(Idents(pancT))
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Naïve/central memory CD4 T cells"] <- "CD4+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Tissue resident memory CD4 T cells"] <- "CD4+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Effector memory CD4 T cells"] <- "CD4+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Tregs"] <- "Tregs"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Naïve/central memory CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Tissue resident memory CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Effector memory CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Cytotoxic CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Follicular CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "IFNg+ CD8 T cells"] <- "CD8+ T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Proliferating/Cycling T cells"] <- "Proliferating/Cycling T cells"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "MAIT"] <- "MAIT"
pancT@meta.data$Cell_Type_Broad[pancT@meta.data$Cell_Type == "Gamma delta T cells"] <- "Gamma delta T cells"
Idents(pancT) <- "Cell_Type_Broad"
table(Idents(pancT))

                                                                           ### Subset to CD4 and CD8 cells -----
pancCD4T <- subset(pancT, ident = c("CD4+ T cells", "Tregs"))
pancCD8T <- subset(pancT, ident = "CD8+ T cells")

### Pancreas CD4 and CD8 TCR gene sets -----
NeoTCR_genes <- read_csv("/hpc/group/peterallen/af206/scRNA-seq/NeoTCR/NeoTCR_genes.csv")
Offringa_TR_Signature <- read_csv("/hpc/group/peterallen/af206/scRNA-seq/NeoTCR/Offringa_TR_Signature.csv")
Eric_Tran_ExRe_CD4 <- read_csv("/hpc/group/peterallen/af206/scRNA-seq/NeoTCR/Eric_Tran_ExRe_CD4.csv")

NeoTCR4 <- unique(NeoTCR_genes$NeoTCR4)
NeoTCR8 <- unique(NeoTCR_genes$NeoTCR8)
OffTR30 <- unique(Offringa_TR_Signature$`Meng.TR-30`)
ExRe4 <- unique(Eric_Tran_ExRe_CD4$Zheng.ExRe.CD4)

### Get expression matrices -----
exprMatrix4 <- GetAssayData(object = pancCD4T, assay = "RNA", layer = "counts")
head(exprMatrix4)

exprMatrix8 <- GetAssayData(object = pancCD8T, assay = "RNA", layer = "counts")
head(exprMatrix8)

### Assign AUC scores for Lowry gene sets -----
TCR4_AUC <- AUCell_run(exprMatrix4, NeoTCR4)
TCR8_AUC <- AUCell_run(exprMatrix8, NeoTCR8)
AUC4 <- TCR4_AUC@assays@data$AUC
AUC4 <- as.data.frame(t(AUC4))
AUC8 <- TCR8_AUC@assays@data$AUC
AUC8 <- as.data.frame(t(AUC8))
AUC4 <- AUC4 %>%
  arrange(desc(geneSet))
AUC8 <- AUC8 %>%
  arrange(desc(geneSet))

### Lowery CD4 and CD8 top 5% -----
top5.4 <- AUC4 %>%
  dplyr::slice(1:floor(0.05 * n()))

pancT@meta.data$NeoTCR4_Top5per <- "NA"
for (i in 1:nrow(pancT@meta.data)) {
  match_row <- which(rownames(top5.4) == rownames(pancT@meta.data)[i])
  if (length(match_row) > 0) {
    pancT@meta.data$NeoTCR4_Top5per[i] <- round(top5.4$geneSet[match_row[1]], digits = 3)
  }
}

head(pancT@meta.data$NeoTCR4_Top5per)

top5.8 <- AUC8 %>%
  dplyr::slice(1:floor(0.05 * n()))

pancT@meta.data$NeoTCR8_Top5per <- "NA"
for (i in 1:nrow(pancT@meta.data)) {
  match_row <- which(rownames(top5.8) == rownames(pancT@meta.data)[i])
  if (length(match_row) > 0) {
    pancT@meta.data$NeoTCR8_Top5per[i] <- round(top5.8$geneSet[match_row[1]], digits = 3)
  }
}

head(pancT@meta.data$NeoTCR8_Top5per)

### Dimplot
DimPlot(pancT, label=F, group.by="sample.type", split.by = NULL, combine = T,
        #scale_color_manual(values = colors_35) +
        cells.highlight= list(Lowry.CD4.Top_5=rownames(top5.4),
                              Lowry.CD8.Top_5=rownames(top5.8)), sizes.highlight = 0.5, #ncol = 3,
        cols.highlight = c("red","blue3"), cols = "grey") + #+ NoLegend() + , cols= "grey"
  ggtitle("") 

# *Fig. 4D -----
### Meng CD8 top 5% -----
MengTCR8_AUC <- AUCell_run(exprMatrix8, OffTR30)
MengTCR8_AUC <- MengTCR8_AUC@assays@data$AUC
MengTCR8_AUC <- as.data.frame(t(MengTCR8_AUC))
MengTCR8_AUC <- MengTCR8_AUC %>%
  arrange(desc(geneSet))

top5.M <- MengTCR8_AUC %>%
  dplyr::slice(1:floor(0.05 * n()))

pancT@meta.data$MengTCR8_Top5per <- "NA"
for (i in 1:nrow(pancT@meta.data)) {
  match_row <- which(rownames(top5.M) == rownames(pancT@meta.data)[i])
  if (length(match_row) > 0) {
    pancT@meta.data$MengTCR8_Top5per[i] <- round(top5.M$geneSet[match_row[1]], digits = 3)
  }
}

head(pancT@meta.data$MengTCR8_Top5per)

### Zheng CD4 top 5% -----
ZhengTCR4_AUC <- AUCell_run(exprMatrix4, ExRe4)
ZhengTCR4_AUC <- ZhengTCR4_AUC@assays@data$AUC
ZhengTCR4_AUC <- as.data.frame(t(ZhengTCR4_AUC))
ZhengTCR4_AUC <- ZhengTCR4_AUC %>%
  arrange(desc(geneSet))

top5.Z <- ZhengTCR4_AUC %>%
  dplyr::slice(1:floor(0.05 * n()))

pancT@meta.data$ZhengTCR4_Top5per <- "NA"
for (i in 1:nrow(pancT@meta.data)) {
  match_row <- which(rownames(top5.Z) == rownames(pancT@meta.data)[i])
  if (length(match_row) > 0) {
    pancT@meta.data$ZhengTCR4_Top5per[i] <- round(top5.Z$geneSet[match_row[1]], digits = 3)
  }
}

head(pancT@meta.data$ZhengTCR4_Top5per)

### Dimplot -----
DimPlot(pancT, label=F, group.by="Cell_Type", split.by = "sample.type", combine = T,
        cells.highlight= list(Zheng.CD4.Top_5=rownames(top5.Z),
                              Meng.CD8.Top_5=rownames(top5.M)), sizes.highlight = 0.5, 
        cols.highlight = c("blue3","red"), cols = "grey") + 
  ggtitle("")

# *Fig. 4E -----
DefaultAssay(pancT) <- "RNA"
genes1 <-  c("CD4","CD8A", "FOXP3", "ITGAE", "TOX", "TOX")
genes2 <- c("CXCL13", "ENTPD1", "PDCD1", "TIGIT", "HAVCR2")

FeaturePlot(pancT, , order = T, features = genes1,
            cols = c("lightgrey", "red4"),
            ncol = 5, label.size = 3) 
FeaturePlot(pancT, , order = T, features = genes2,
            cols = c("lightgrey", "red4"),
            ncol = 5, label.size = 3) 






