if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2","org.At.tair.db","clusterProfiler",
                       "enrichplot","pheatmap","apeglm"))
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

setwd("D:/academic/PhD/WGCNA paper")

#load the matrix
counts_raw <- read.csv(
  "D:/academic/PhD/WGCNA paper/Combined_RAW_Counts.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

dim(counts_raw)
head(counts_raw[,1:5])
#column rename
colnames(counts_raw) <- c(
  "OE_Ctrl_1","OE_Ctrl_2","OE_Ctrl_3",
  "OE_Man_1","OE_Man_2","OE_Man_3",
  "Col_Ctrl_1","Col_Ctrl_2","Col_Ctrl_3",
  "Col_Man_1","Col_Man_2","Col_Man_3"
)
# ============================================
# 2. Sample metadata — YOUR 4 GROUPS
# ============================================
coldata <- data.frame(
  genotype  = factor(c(rep("OE",6), rep("Col0",6)),
                     levels=c("Col0","OE")),          # Col0 = reference
  treatment = factor(c(rep("Control",3), rep("Mannitol",3),
                       rep("Control",3), rep("Mannitol",3)),
                     levels=c("Control","Mannitol")), # Control = reference
  group     = factor(c(
    rep("OE_Ctrl",3), rep("OE_Man",3),
    rep("Col_Ctrl",3), rep("Col_Man",3)
  )),
  row.names = colnames(counts_raw)
)

print(coldata)
# ============================================
# 3. DESeq2 object
# ============================================
dds <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData   = coldata,
  design    = ~ genotype + treatment + genotype:treatment
)
# Filter: keep genes with ≥10 counts in ≥3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")
# Run DESeq2
dds <- DESeq(dds)
# ============================================
# 4. All 5 key comparisons for your design
# ============================================

## --- Comparison 1: Mannitol effect in Col-0 (WT stress response) ---
res_Col_ManVsCtrl <- results(dds,
                             contrast      = c("treatment","Mannitol","Control"),
                             alpha         = 0.05,
                             independentFiltering = TRUE)
## --- Comparison 2: Mannitol effect in OE-7.7 (transgene stress response) ---
res_OE_ManVsCtrl <- results(dds,
                            list(c("treatment_Mannitol_vs_Control",
                                   "genotypeOE.treatmentMannitol")),
                            alpha = 0.05)
## --- Comparison 3: OE-7.7 vs Col-0 under Control (transgene effect alone) ---
res_OEvsCol_Ctrl <- results(dds,
                            contrast      = c("genotype","OE","Col0"),
                            alpha         = 0.05)
## --- Comparison 4: OE-7.7 vs Col-0 under Mannitol (combined effect) ---
res_OEvsCol_Man <- results(dds,
                           list(c("genotype_OE_vs_Col0",
                                  "genotypeOE.treatmentMannitol")),
                           alpha = 0.05)
## --- Comparison 5: Interaction (unique OE mannitol response) ---
res_interaction <- results(dds,
                           name  = "genotypeOE.treatmentMannitol",
                           alpha = 0.05)
# ============================================
# 5. Filter DEGs
# ============================================
get_degs <- function(res, fc_cutoff=1, fdr=0.05, label="") {
  df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj),
           padj < fdr,
           abs(log2FoldChange) >= fc_cutoff) %>%
    mutate(
      regulation = ifelse(log2FoldChange > 0, "Up", "Down"),
      comparison = label
    ) %>%
    arrange(padj)
  cat(label, "- UP:", sum(df$regulation=="Up"),
      "| DOWN:", sum(df$regulation=="Down"), "\n")
  return(df)
}

degs_Col_stress   <- get_degs(res_Col_ManVsCtrl, label="Col0_Mannitol_vs_Control")
degs_OE_stress    <- get_degs(res_OE_ManVsCtrl,  label="OE_Mannitol_vs_Control")
degs_genotype     <- get_degs(res_OEvsCol_Ctrl,  label="OE_vs_Col0_Control")
degs_combined     <- get_degs(res_OEvsCol_Man,   label="OE_vs_Col0_Mannitol")
degs_interaction  <- get_degs(res_interaction,   label="Interaction_OE_x_Mannitol")
# Save all DEG tables
write.csv(degs_Col_stress,  "results/DEGs/Col0_Mannitol_vs_Control.csv",  row.names=FALSE)
write.csv(degs_OE_stress,   "results/DEGs/OE_Mannitol_vs_Control.csv",    row.names=FALSE)
write.csv(degs_genotype,    "results/DEGs/OE_vs_Col0_Control.csv",        row.names=FALSE)
write.csv(degs_combined,    "results/DEGs/OE_vs_Col0_Mannitol.csv",       row.names=FALSE)
write.csv(degs_interaction, "results/DEGs/Interaction_OE_x_Mannitol.csv", row.names=FALSE)
# ============================================
# 6. Visualization
# ============================================

# --- PCA ---
vsd <- vst(dds, blind=FALSE)

pca_plot <- plotPCA(vsd, intgroup=c("genotype","treatment")) +
  scale_color_manual(values=c(
    "OE:Control"   = "#2196F3",
    "OE:Mannitol"  = "#F44336",
    "Col0:Control" = "#4CAF50",
    "Col0:Mannitol"= "#FF9800"
  )) +
  geom_point(size=4) +
  theme_bw(base_size=14) +
  ggtitle("PCA: OE-7.7 vs Col-0 | Control vs Mannitol")

ggsave("results/PCA_all_samples.pdf", pca_plot, width=8, height=6)
# --- Volcano plot function ---
volcano_plot <- function(res_df, title, fc=1) {
  res_df_plot <- as.data.frame(res_df) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = case_when(
        padj < 0.05 & log2FoldChange >  fc ~ "Up-regulated",
        padj < 0.05 & log2FoldChange < -fc ~ "Down-regulated",
        TRUE ~ "Not significant"
      )
    )
  
  p <- ggplot(res_df_plot,
              aes(x=log2FoldChange, y=-log10(padj),
                  color=sig)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c(
      "Up-regulated"   = "#E53935",
      "Down-regulated" = "#1E88E5",
      "Not significant"= "grey70"
    )) +
    geom_vline(xintercept=c(-fc,fc), linetype="dashed",
               color="black", linewidth=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               color="black", linewidth=0.5) +
    annotate("text",
             x=max(res_df_plot$log2FoldChange,na.rm=TRUE)*0.8,
             y=max(-log10(res_df_plot$padj),na.rm=TRUE)*0.95,
             label=paste("Up:",   sum(res_df_plot$sig=="Up-regulated")),
             color="#E53935", size=4) +
    annotate("text",
             x=min(res_df_plot$log2FoldChange,na.rm=TRUE)*0.8,
             y=max(-log10(res_df_plot$padj),na.rm=TRUE)*0.95,
             label=paste("Down:", sum(res_df_plot$sig=="Down-regulated")),
             color="#1E88E5", size=4) +
    theme_bw(base_size=13) +
    labs(title=title, color="",
         x="log2(Fold Change)", y="-log10(adj. p-value)")
  
  return(p)
}

# Generate all volcano plots
v1 <- volcano_plot(res_Col_ManVsCtrl, "Col-0: Mannitol vs Control")
v2 <- volcano_plot(res_OE_ManVsCtrl,  "OE-7.7: Mannitol vs Control")
v3 <- volcano_plot(res_OEvsCol_Ctrl,  "OE-7.7 vs Col-0: Control")
v4 <- volcano_plot(res_interaction,   "Interaction: OE × Mannitol")

ggsave("results/DEGs/Volcano_Col0_stress.pdf",    v1, width=8, height=6)
ggsave("results/DEGs/Volcano_OE_stress.pdf",      v2, width=8, height=6)
ggsave("results/DEGs/Volcano_Genotype_ctrl.pdf",  v3, width=8, height=6)
ggsave("results/DEGs/Volcano_Interaction.pdf",    v4, width=8, height=6)
# --- Heatmap of top DEGs ---
top_degs <- c(
  head(degs_Col_stress$gene_id,  25),
  head(degs_OE_stress$gene_id,   25),
  head(degs_interaction$gene_id, 25)
) %>% unique()

anno_col <- data.frame(
  Genotype  = coldata$genotype,
  Treatment = coldata$treatment,
  row.names = rownames(coldata)
)

anno_colors <- list(
  Genotype  = c(Col0="#4CAF50", OE="#2196F3"),
  Treatment = c(Control="#90CAF9", Mannitol="#EF9A9A")
)

pheatmap(assay(vsd)[top_degs, ],
         annotation_col  = anno_col,
         annotation_colors = anno_colors,
         scale           = "row",
         show_rownames   = FALSE,
         cluster_cols    = TRUE,
         color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main            = "Top DEGs Heatmap",
         filename        = "results/DEGs/Heatmap_TopDEGs.pdf",
         width=10, height=10
)
# --- Venn diagram: shared/unique DEGs ---
install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(
  x = list(
    Col0_Stress = degs_Col_stress$gene_id,
    OE_Stress   = degs_OE_stress$gene_id,
    Interaction = degs_interaction$gene_id
  ),
  filename   = "results/DEGs/Venn_DEGs.png",
  fill       = c("#4CAF50","#2196F3","#FF9800"),
  alpha      = 0.5,
  main       = "DEG Overlap Across Comparisons"
)
#####KEGG enrichment
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tibble)

# ============================================
# Create output directory
# ============================================
dir.create("results/enrichment",
           recursive = TRUE,
           showWarnings = FALSE)

# ============================================
# Universal enrichment function
# ============================================
enrich_arabidopsis <- function(gene_list, label){
  
  cat("\nRunning enrichment for:", label, "\n")
  
  # Clean TAIR IDs
  genes <- unique(gsub("\\.\\d+$", "", gene_list))
  
  # ----------------------------------------
  # GO enrichment
  # ----------------------------------------
  go_res <- enrichGO(
    gene          = genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  # ----------------------------------------
  # TAIR → Entrez conversion
  # ----------------------------------------
  entrez_df <- bitr(
    genes,
    fromType = "TAIR",
    toType   = "ENTREZID",
    OrgDb    = org.At.tair.db
  )
  
  entrez <- unique(entrez_df$ENTREZID)
  
  cat("Genes:", length(genes), "\n")
  cat("Mapped Entrez IDs:", length(entrez), "\n")
  
  # ----------------------------------------
  # KEGG enrichment
  # ----------------------------------------
  kegg_res <- NULL
  
  if(length(entrez) > 0){
    
    kegg_res <- tryCatch({
      
      enrichKEGG(
        gene         = entrez,
        organism     = "ath",
        keyType      = "ncbi-geneid",
        pvalueCutoff = 0.05
      )
      
    }, error=function(e){
      
      message("KEGG failed for ", label)
      NULL
      
    })
  }
  
  # ----------------------------------------
  # Save GO results
  # ----------------------------------------
  if(!is.null(go_res) &&
     nrow(as.data.frame(go_res)) > 0){
    
    write.csv(
      as.data.frame(go_res),
      paste0("results/enrichment/GO_", label, ".csv"),
      row.names = FALSE
    )
    
    p1 <- dotplot(
      go_res,
      showCategory = 20,
      split = "ONTOLOGY"
    ) +
      facet_grid(ONTOLOGY~., scales="free") +
      ggtitle(paste("GO Enrichment -", label))
    
    ggsave(
      paste0("results/enrichment/GO_", label, ".pdf"),
      p1,
      width = 10,
      height = 12
    )
    
    cat("GO terms found:",
        nrow(as.data.frame(go_res)), "\n")
    
  } else {
    
    cat("No significant GO terms.\n")
    
  }
  
  # ----------------------------------------
  # Save KEGG results
  # ----------------------------------------
  if(!is.null(kegg_res) &&
     nrow(as.data.frame(kegg_res)) > 0){
    
    write.csv(
      as.data.frame(kegg_res),
      paste0("results/enrichment/KEGG_", label, ".csv"),
      row.names = FALSE
    )
    
    p2 <- dotplot(
      kegg_res,
      showCategory = 20
    ) +
      ggtitle(paste("KEGG Enrichment -", label))
    
    ggsave(
      paste0("results/enrichment/KEGG_", label, ".pdf"),
      p2,
      width = 10,
      height = 8
    )
    
    cat("KEGG pathways found:",
        nrow(as.data.frame(kegg_res)), "\n")
    
  } else {
    
    cat("No significant KEGG pathways.\n")
    
  }
  
  return(list(
    GO = go_res,
    KEGG = kegg_res
  ))
}

# ============================================
# Run enrichment analyses
# ============================================

enrich_Col_stress <- enrich_arabidopsis(
  degs_Col_stress$gene_id,
  "Col0_Mannitol_vs_Ctrl"
)

enrich_OE_stress <- enrich_arabidopsis(
  degs_OE_stress$gene_id,
  "OE_Mannitol_vs_Ctrl"
)

enrich_genotype <- enrich_arabidopsis(
  degs_genotype$gene_id,
  "OE_vs_Col0_Control"
)

enrich_interaction <- enrich_arabidopsis(
  degs_interaction$gene_id,
  "Interaction_OExMannitol"
)

# ============================================
# Compare KEGG enrichment among DEG sets
# ============================================

gene_clusters <- list(
  Col0_Stress =
    unique(bitr(
      degs_Col_stress$gene_id,
      "TAIR",
      "ENTREZID",
      org.At.tair.db
    )$ENTREZID),
  
  OE_Stress =
    unique(bitr(
      degs_OE_stress$gene_id,
      "TAIR",
      "ENTREZID",
      org.At.tair.db
    )$ENTREZID),
  
  OE_vs_Col0 =
    unique(bitr(
      degs_genotype$gene_id,
      "TAIR",
      "ENTREZID",
      org.At.tair.db
    )$ENTREZID),
  
  Interaction =
    unique(bitr(
      degs_interaction$gene_id,
      "TAIR",
      "ENTREZID",
      org.At.tair.db
    )$ENTREZID)
)

compare_kegg <- tryCatch({
  
  compareCluster(
    geneClusters = gene_clusters,
    fun = "enrichKEGG",
    organism = "ath",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
  
}, error=function(e){
  
  NULL
  
})

# ============================================
# Plot KEGG comparison
# ============================================

if(!is.null(compare_kegg) &&
   nrow(as.data.frame(compare_kegg)) > 0){
  
  p_compare <- dotplot(
    compare_kegg,
    showCategory = 15
  ) +
    ggtitle("KEGG Pathway Comparison Across Groups") +
    theme(
      axis.text.x =
        element_text(angle = 45,
                     hjust = 1)
    )
  
  ggsave(
    "results/enrichment/KEGG_CompareAll.pdf",
    p_compare,
    width = 14,
    height = 10
  )
  
  write.csv(
    as.data.frame(compare_kegg),
    "results/enrichment/KEGG_CompareAll.csv",
    row.names = FALSE
  )
  
} else {
  
  message(
    "No significant KEGG pathways detected across comparisons."
  )
  
}
# ============================================================
# COMPLETE WGCNA PIPELINE
# Design: OE-7.7 vs Col-0 | Control vs Mannitol
# Samples: 12 (3 replicates × 4 groups)
# ============================================================

# ============================================================
# SECTION 1: INSTALL & LOAD PACKAGES
# ============================================================

if (!require("BiocManager")) install.packages("BiocManager")

# Install Bioconductor dependencies
bioc_pkgs <- c("impute", "preprocessCore", "GO.db",
               "AnnotationDbi", "DESeq2", "org.At.tair.db",
               "clusterProfiler", "enrichplot")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}

# Install CRAN packages
cran_pkgs <- c("WGCNA", "dynamicTreeCut", "fastcluster",
               "tidyverse", "pheatmap", "RColorBrewer",
               "ggplot2", "VennDiagram", "igraph")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Load in correct order (impute MUST come before WGCNA)
library(impute)
library(preprocessCore)
library(GO.db)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# WGCNA global settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)  # adjust to your CPU cores

# Confirm loaded correctly
cat("==============================================\n")
cat("WGCNA version    :", as.character(packageVersion("WGCNA")), "\n")
cat("goodSamplesGenes :", exists("goodSamplesGenes"), "\n")
cat("==============================================\n")



# Create output directories
dirs <- c("results/WGCNA/plots",
          "results/WGCNA/hub_genes",
          "results/WGCNA/cytoscape",
          "results/WGCNA/tables")

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("Output folders created\n")

# ============================================================
# SECTION 3: PREPARE EXPRESSION MATRIX FROM DESeq2
# ============================================================

# VST-normalize counts from your DESeq2 object (dds)
vsd     <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

cat("VST matrix dimensions:", dim(vsd_mat), "\n")
# Expected: ~20000+ genes × 12 samples

# ------------------------------------------------------------
# Combine DEGs from all 4 comparisons as WGCNA gene background
# (mirrors paper approach: strict FC filter before WGCNA)
# ------------------------------------------------------------

# Your 4 DEG tables from DESeq2 section
wgcna_genes <- bind_rows(
  degs_Col_stress,    # Col-0 Mannitol vs Control
  degs_OE_stress,     # OE-7.7 Mannitol vs Control
  degs_genotype,      # OE-7.7 vs Col-0 under Control
  degs_interaction    # Interaction: OE × Mannitol
) %>%
  filter(
    abs(log2FoldChange) >= log2(3),  # FC ≥ 3 (paper used FC≥5;
    # relaxed here for 12 samples)
    padj < 0.05
  ) %>%
  pull(gene_id) %>%
  unique()

cat("================================================\n")
cat("Genes entering WGCNA :", length(wgcna_genes), "\n")
cat("Target range         : 500 – 5000 genes\n")
cat("================================================\n")

# Adjust FC cutoff if needed:
# Too few genes (<500)  → lower to log2(2)
# Too many genes (>5000) → raise to log2(5)

# Build expression matrix (samples as ROWS — required by WGCNA)
expr_wgcna <- t(vsd_mat[wgcna_genes, ])

cat("Matrix: samples × genes =", dim(expr_wgcna), "\n")
# Expected: 12 × 3528

# ============================================================
# SECTION 4: SAMPLE QUALITY CHECK
# ============================================================

# Check for bad samples or genes
gsg <- goodSamplesGenes(expr_wgcna, verbose = 3)

if (!gsg$allOK) {
  cat("\n⚠️  Removing problematic samples/genes...\n")
  if (sum(!gsg$goodGenes) > 0)
    cat("   Genes removed   :", sum(!gsg$goodGenes), "\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("   Samples removed  :", sum(!gsg$goodSamples), "\n")
  expr_wgcna <- expr_wgcna[gsg$goodSamples, gsg$goodGenes]
} else {
  cat("✅ All samples and genes passed quality check\n")
}

cat("Final matrix:", nrow(expr_wgcna), "samples ×",
    ncol(expr_wgcna), "genes\n")

# ------------------------------------------------------------
# Sample clustering to visually detect outliers
# ------------------------------------------------------------

sampleTree <- hclust(dist(expr_wgcna), method = "average")

pdf("results/WGCNA/plots/01_SampleClustering.pdf",
    width = 12, height = 6)
par(mar = c(0, 5, 4, 0))
plot(sampleTree,
     main   = "Sample Clustering — Check for Outliers",
     sub    = "",
     xlab   = "",
     ylab   = "Height",
     cex    = 0.9,
     cex.main = 1.2)
# Add colored bar for groups
group_colors <- c(
  rep("#90CAF9", 3),   # OE Control   — light blue
  rep("#1565C0", 3),   # OE Mannitol  — dark blue
  rep("#A5D6A7", 3),   # Col Control  — light green
  rep("#2E7D32", 3)    # Col Mannitol — dark green
)
abline(h = 200, col = "red", lty = 2)  # adjust cutoff if needed
legend("topright",
       legend = c("OE_Control","OE_Mannitol",
                  "Col0_Control","Col0_Mannitol"),
       fill   = c("#90CAF9","#1565C0","#A5D6A7","#2E7D32"),
       cex    = 0.8, bty = "n")
dev.off()

cat("✅ Sample clustering plot saved\n")

# ------------------------------------------------------------
# Trait heatmap alongside sample dendrogram
# ------------------------------------------------------------

# YOUR EXACT 12 SAMPLES in this order:
# OE_Ctrl_1, OE_Ctrl_2, OE_Ctrl_3,
# OE_Man_1,  OE_Man_2,  OE_Man_3,
# Col_Ctrl_1,Col_Ctrl_2,Col_Ctrl_3,
# Col_Man_1, Col_Man_2, Col_Man_3

traits <- data.frame(
  Genotype    = c(1, 1, 1,   # OE Control   → OE=1
                  1, 1, 1,   # OE Mannitol  → OE=1
                  0, 0, 0,   # Col Control  → Col0=0
                  0, 0, 0),  # Col Mannitol → Col0=0
  Mannitol    = c(0, 0, 0,   # OE Control   → no stress
                  1, 1, 1,   # OE Mannitol  → stress
                  0, 0, 0,   # Col Control  → no stress
                  1, 1, 1),  # Col Mannitol → stress
  Interaction = c(0, 0, 0,   # OE Control   → no interaction
                  1, 1, 1,   # OE Mannitol  → OE+stress = 1
                  0, 0, 0,   # Col Control  → no interaction
                  0, 0, 0),  # Col Mannitol → no interaction
  row.names   = rownames(expr_wgcna)
)

traitColors <- numbers2colors(traits, signed = FALSE)

pdf("results/WGCNA/plots/02_SampleDendro_TraitHeatmap.pdf",
    width = 12, height = 8)
plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels   = c("Genotype\n(OE=1)",
                    "Mannitol\n(Stress=1)",
                    "Interaction\n(OE×Stress)"),
  main          = "Sample Clustering and Trait Assignment",
  cex.colorLabels = 0.8,
  cex.dendroLabels = 0.9
)
dev.off()

cat("✅ Trait heatmap saved\n")

# ============================================================
# SECTION 5: SOFT THRESHOLD SELECTION
# ============================================================

library(WGCNA)

powers <- c(1:20)

sft <- pickSoftThreshold(
  expr_wgcna,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# Examine results
fitIndices <- sft$fitIndices

# Try automatic estimate first
chosen_beta <- sft$powerEstimate

# If NA, choose manually
if(is.na(chosen_beta)){
  
  cat("Automatic power selection failed.\n")
  
  # Find first power with R² > 0.7
  candidate <- fitIndices$Power[
    fitIndices$SFT.R.sq > 0.7
  ]
  
  if(length(candidate) > 0){
    
    chosen_beta <- min(candidate)
    
  } else {
    
    # Default for small RNA-seq datasets
    chosen_beta <- 8
    
  }
  
}

cat("Selected soft-threshold power =", chosen_beta, "\n")

# ============================================================
# SECTION 6: BUILD CO-EXPRESSION NETWORK
# ============================================================

dir.create("results/WGCNA",
           recursive = TRUE,
           showWarnings = FALSE)

cat("Building network... this may take several minutes\n")

net <- blockwiseModules(
  expr_wgcna,
  
  power = chosen_beta,
  
  networkType = "signed",
  TOMType = "signed",
  
  minModuleSize = 50,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  
  reassignThreshold = 0,
  pamRespectsDendro = FALSE,
  
  numericLabels = FALSE,
  
  saveTOMs = FALSE,
  
  verbose = 3
)

# ============================================================
# MODULE SUMMARY
# ============================================================

module_table <- sort(
  table(net$colors),
  decreasing = TRUE
)

print(module_table)

write.csv(
  data.frame(
    Module = names(module_table),
    Genes = as.numeric(module_table)
  ),
  "results/WGCNA/Module_Summary.csv",
  row.names = FALSE
)

# Module summary
module_table <- sort(table(net$colors), decreasing = TRUE)
cat("\n================================================\n")
cat("Module detection complete!\n")
cat("Number of modules:", length(unique(net$colors)) - 1,
    "(excluding grey)\n")
cat("\nModule sizes:\n")
print(module_table)
cat("================================================\n")

# Save module assignments
module_df <- data.frame(
  Gene_ID = names(net$colors),
  Module  = net$colors
)
write.csv(module_df,
          "results/WGCNA/tables/Gene_Module_Assignment.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# Plot gene dendrogram with module colors
# ------------------------------------------------------------

pdf("results/WGCNA/plots/04_GeneDendrogram_Modules.pdf",
    width = 14, height = 7)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene Dendrogram and Module Colors\nOE-7.7 vs Col-0 | Mannitol Stress"
)
dev.off()

cat("✅ Gene dendrogram saved\n")

# ============================================================
# SECTION 7: MODULE EIGENGENES & TRAIT CORRELATION
# ============================================================

# Calculate module eigengenes (first PC of each module)
MEs_raw <- moduleEigengenes(expr_wgcna, net$colors)$eigengenes

# Order module eigengenes
MEs <- orderMEs(MEs_raw)

# Correlate eigengenes with traits
MEs_cor  <- cor(MEs, traits, use = "p")
MEs_pval <- corPvalueStudent(MEs_cor, nrow(expr_wgcna))

cat("\nModule-Trait correlation matrix:\n")
print(round(MEs_cor, 3))
cat("\nModule-Trait p-value matrix:\n")
print(round(MEs_pval, 3))

# Save correlation tables
write.csv(
  data.frame(Module = rownames(MEs_cor), MEs_cor),
  "results/WGCNA/tables/ModuleTrait_Correlation.csv",
  row.names = FALSE
)
write.csv(
  data.frame(Module = rownames(MEs_pval), MEs_pval),
  "results/WGCNA/tables/ModuleTrait_Pvalue.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Module-Trait heatmap (key figure — mirrors paper Figure 4B)
# ------------------------------------------------------------

pdf("results/WGCNA/plots/05_ModuleTrait_Heatmap.pdf",
    width = 8, height = 10)

textMatrix <- paste0(
  signif(MEs_cor,  2), "\n",
  "(p=", signif(MEs_pval, 1), ")"
)
dim(textMatrix) <- dim(MEs_cor)

par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix          = MEs_cor,
  xLabels         = c("Genotype\n(OE=1)",
                      "Mannitol\n(Stress=1)",
                      "Interaction\n(OE×Mannitol)"),
  yLabels         = names(MEs),
  ySymbols        = names(MEs),
  colorLabels     = FALSE,
  colors          = blueWhiteRed(50),
  textMatrix      = textMatrix,
  setStdMargins   = FALSE,
  cex.text        = 0.65,
  zlim            = c(-1, 1),
  main            = "Module–Trait Relationships\nOE-7.7 vs Col-0 | Mannitol Stress"
)

dev.off()

cat("✅ Module-Trait heatmap saved\n")

# ------------------------------------------------------------
# Module eigengene expression barplot across your 4 groups
# ------------------------------------------------------------

# Identify significant modules (p < 0.05 for any trait)
sig_modules <- rownames(MEs_pval)[
  apply(MEs_pval, 1, function(x) any(x < 0.05))
]
sig_modules <- gsub("^ME", "", sig_modules)

cat("\n================================================\n")
cat("Significant modules (p<0.05):",
    paste(sig_modules, collapse=", "), "\n")
cat("================================================\n")

# Plot eigengene expression per group for significant modules
group_labels <- c(rep("OE_Control",  3),
                  rep("OE_Mannitol", 3),
                  rep("Col_Control", 3),
                  rep("Col_Mannitol",3))

ME_long <- MEs %>%
  rownames_to_column("sample") %>%
  mutate(group = group_labels) %>%
  pivot_longer(cols      = starts_with("ME"),
               names_to  = "Module",
               values_to = "Eigengene") %>%
  filter(gsub("^ME","", Module) %in% sig_modules)

group_colors <- c(
  "OE_Control"   = "#90CAF9",
  "OE_Mannitol"  = "#1565C0",
  "Col_Control"  = "#A5D6A7",
  "Col_Mannitol" = "#2E7D32"
)

p_eigen <- ggplot(ME_long,
                  aes(x = group, y = Eigengene, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  facet_wrap(~ Module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  ) +
  labs(
    x     = NULL,
    y     = "Module Eigengene",
    title = "Eigengene Expression by Group — Significant Modules"
  )

ggsave("results/WGCNA/plots/06_Eigengene_byGroup.pdf",
       p_eigen, width = 14, height = 10)

cat("✅ Eigengene expression plots saved\n")

# ============================================================
# SECTION 8: HUB GENE IDENTIFICATION
# ============================================================

# Paper criteria: weight > 0.35 AND degree > 5
# We add: Module Membership > 0.8 for robustness

hub_all    <- list()
mm_gs_all  <- list()

for (mod in sig_modules) {
  
  cat("\n--- Processing module:", mod, "---\n")
  
  mod_genes <- names(net$colors)[net$colors == mod]
  cat("    Total genes in module:", length(mod_genes), "\n")
  
  mod_expr  <- expr_wgcna[, mod_genes, drop = FALSE]
  
  # Adjacency matrix
  adj_mat <- adjacency(mod_expr,
                       power = chosen_beta,
                       type  = "signed")
  
  # Intramodular connectivity (degree)
  kIM <- rowSums(adj_mat) - 1
  
  # Module eigengene for this module
  ME_mod <- MEs[, paste0("ME", mod), drop = FALSE]
  
  # Module Membership (MM) — correlation with eigengene
  MM <- cor(mod_expr, ME_mod, use = "p")
  
  # Gene Significance (GS) — correlation with each trait
  GS_genotype    <- cor(mod_expr, traits$Genotype,    use = "p")
  GS_mannitol    <- cor(mod_expr, traits$Mannitol,    use = "p")
  GS_interaction <- cor(mod_expr, traits$Interaction, use = "p")
  
  # Build full MM/GS table for all genes in module
  mm_gs_df <- data.frame(
    Gene_ID        = mod_genes,
    Module         = mod,
    kIM            = round(kIM, 3),
    MM             = round(MM[, 1], 3),
    GS_Genotype    = round(GS_genotype[, 1], 3),
    GS_Mannitol    = round(GS_mannitol[, 1], 3),
    GS_Interaction = round(GS_interaction[, 1], 3)
  ) %>%
    arrange(desc(abs(MM)))
  
  mm_gs_all[[mod]] <- mm_gs_df
  
  # Hub gene filter: kIM > 5 AND |MM| > 0.8
  hubs <- mm_gs_df %>%
    filter(kIM    > 5,
           abs(MM) > 0.8) %>%
    arrange(desc(kIM))
  
  hub_all[[mod]] <- hubs
  
  cat("    Hub genes identified:", nrow(hubs), "\n")
  
  # Save per-module tables
  write.csv(mm_gs_df,
            paste0("results/WGCNA/hub_genes/MM_GS_", mod, ".csv"),
            row.names = FALSE)
  write.csv(hubs,
            paste0("results/WGCNA/hub_genes/HubGenes_", mod, ".csv"),
            row.names = FALSE)
  
  # MM vs GS scatter plot (quality check)
  p_mmgs <- ggplot(mm_gs_df,
                   aes(x = abs(MM),
                       y = abs(GS_Mannitol),
                       color = kIM)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_gradient(low = "grey80", high = "red") +
    geom_vline(xintercept = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 0.2, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(
      x     = "Module Membership |MM|",
      y     = "Gene Significance |GS| for Mannitol",
      color = "Connectivity",
      title = paste("Module:", mod,
                    "| Genes:", length(mod_genes),
                    "| Hubs:", nrow(hubs))
    )
  
  ggsave(paste0("results/WGCNA/plots/07_MM_GS_", mod, ".pdf"),
         p_mmgs, width = 7, height = 6)
}

# Combined hub gene table
hub_combined <- bind_rows(hub_all)
write.csv(hub_combined,
          "results/WGCNA/hub_genes/HubGenes_ALL_Modules.csv",
          row.names = FALSE)

cat("\n================================================\n")
cat("Total hub genes across all modules:",
    nrow(hub_combined), "\n")
cat("Hub genes per module:\n")
print(table(hub_combined$Module))
cat("================================================\n")

# ============================================================
# SECTION 9: EXPORT TO CYTOSCAPE
# ============================================================

for (mod in sig_modules) {
  
  cat("Exporting module to Cytoscape:", mod, "\n")
  
  mod_genes <- names(net$colors)[net$colors == mod]
  mod_expr  <- expr_wgcna[, mod_genes, drop = FALSE]
  
  # Compute TOM for this module
  TOM_mat   <- TOMsimilarityFromExpr(mod_expr,
                                     power       = chosen_beta,
                                     networkType = "signed")
  rownames(TOM_mat) <- colnames(TOM_mat) <- mod_genes
  
  # Export (paper used weight threshold = 0.35)
  exportNetworkToCytoscape(
    TOM_mat,
    edgeFile  = paste0("results/WGCNA/cytoscape/",
                       mod, "_edges.txt"),
    nodeFile  = paste0("results/WGCNA/cytoscape/",
                       mod, "_nodes.txt"),
    weighted  = TRUE,
    threshold = 0.35   # paper criterion
  )
  
  cat("   ✅ Saved:", mod, "_edges.txt and _nodes.txt\n")
}

# ============================================================
# SECTION 10: SUMMARY REPORT
# ============================================================

cat("\n")
cat("╔══════════════════════════════════════════════╗\n")
cat("║         WGCNA ANALYSIS COMPLETE              ║\n")
cat("╠══════════════════════════════════════════════╣\n")
cat("║ Genes analysed      :", ncol(expr_wgcna),          "\n")
cat("║ Samples             :", nrow(expr_wgcna),          "\n")
cat("║ Soft threshold β    :", chosen_beta,               "\n")
cat("║ Total modules       :",
    length(unique(net$colors)) - 1,                       "\n")
cat("║ Significant modules :", length(sig_modules),       "\n")
cat("║ Total hub genes     :", nrow(hub_combined),        "\n")
cat("╠══════════════════════════════════════════════╣\n")
cat("║ OUTPUT FILES                                  \n")
cat("║  plots/      → 7 PDF figures                  \n")
cat("║  tables/     → module assignments + corr      \n")
cat("║  hub_genes/  → hub gene lists per module      \n")
cat("║  cytoscape/  → edge/node files per module     \n")
cat("╚══════════════════════════════════════════════╝\n")

# ============================================
# TF identification (download from PlantTFDB)
# http://planttfdb.gao-lab.org/download.php
# Select: Arabidopsis thaliana → All TFs
# ============================================
tf_db <- read.table("D:\\academic\\PhD\\WGCNA paper\\results\\Ath_TF_list.txt\\Ath_TF_list.txt",
                    header=TRUE, sep="\t",
                    col.names=c("Gene_ID","Family","Species"))

# All genes in significant modules
all_sig_genes <- names(net$colors)[net$colors %in% sig_modules]
all_sig_genes <- gsub("\\.\\d+$","", all_sig_genes)

# Find TFs
tfs_in_modules <- tf_db %>%
  filter(Gene_ID %in% all_sig_genes)

# Summary by family
tf_summary <- tfs_in_modules %>%
  count(Family, sort=TRUE)

# Plot TF families
ggplot(head(tf_summary,15),
       aes(x=reorder(Family,-n), y=n, fill=Family)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-0.3, size=3.5) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none") +
  labs(x="TF Family", y="Gene Count",
       title="Transcription Factors in Stress-Related Modules\nOE-7.7 Mannitol Response")

ggsave("results/WGCNA/TF_families.pdf", width=10, height=6)

write.csv(tfs_in_modules, "results/WGCNA/TF_list_in_modules.csv",
          row.names=FALSE)
cat("Total TFs identified:", nrow(tfs_in_modules), "\n")

# ============================================
# Select validation genes (like paper's 6 DEGs)
# ============================================
top_hubs <- hub_combined %>%
  group_by(Module) %>%
  slice_max(kIM, n=2) %>%  # top 2 per module
  ungroup() %>%
  slice_head(n=6)           # total 6 like paper

cat("\nTop hub genes for qRT-PCR validation:\n")
print(top_hubs[, c("Gene","Module","kIM","GS_Mannitol","GS_Interaction")])

write.csv(top_hubs, "results/WGCNA/qRTPCR_candidates.csv", row.names=FALSE)

# Expression profiles of validation genes
val_genes <- top_hubs$Gene
val_expr  <- assay(vsd)[val_genes, ] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample", values_to="VST_expression") %>%
  left_join(
    rownames_to_column(as.data.frame(coldata), "sample"),
    by="sample"
  )

ggplot(val_expr,
       aes(x=interaction(genotype,treatment),
           y=VST_expression,
           fill=interaction(genotype,treatment))) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.1, size=2) +
  facet_wrap(~gene, scales="free_y", ncol=3) +
  scale_fill_manual(values=c(
    "Col0.Control" = "#A5D6A7",
    "Col0.Mannitol"= "#EF9A9A",
    "OE.Control"   = "#90CAF9",
    "OE.Mannitol"  = "#CE93D8"
  )) +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none") +
  labs(x="Group", y="VST Expression",
       title="Hub Gene Expression — qRT-PCR Validation Candidates")

ggsave("results/WGCNA/ValidationGenes_expression.pdf",
       width=12, height=8)