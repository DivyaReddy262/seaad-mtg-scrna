.libPaths("/N/scratch/dikonda/R_libs")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

cat("Loading DE results...\n")
de_results <- read.csv("/N/scratch/dikonda/seurat_project/results/de_results.csv")

dir.create("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/plots",
           recursive = TRUE, showWarnings = FALSE)
dir.create("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/tables",
           recursive = TRUE, showWarnings = FALSE)

cell_types <- unique(de_results$cell_type)
cat("Cell types:", paste(cell_types, collapse=", "), "\n")

for (ct in cell_types) {
  cat("\nRunning pathway enrichment for:", ct, "\n")

  de_ct <- de_results %>%
    filter(cell_type == ct, p_val_adj < 0.05)

  cat("  Significant DE genes:", nrow(de_ct), "\n")

  if (nrow(de_ct) < 10) {
    cat("  Skipping - not enough genes\n")
    next
  }

  entrez <- bitr(de_ct$gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

  if (nrow(entrez) < 5) {
    cat("  Skipping - not enough mapped genes\n")
    next
  }

  # GO Enrichment
  tryCatch({
    go_result <- enrichGO(gene          = entrez$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

    if (!is.null(go_result) && nrow(go_result) > 0) {
      write.csv(as.data.frame(go_result),
                paste0("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/tables/GO_",
                       gsub("/| ", "_", ct), ".csv"),
                row.names = FALSE)

      # Dotplot using ggplot2
      go_df <- as.data.frame(go_result) %>%
        head(20) %>%
        mutate(Description = factor(Description,
                                    levels = rev(Description)))

      p <- ggplot(go_df, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
           geom_point() +
           scale_color_gradient(low = "red", high = "blue") +
           theme_bw() +
           ggtitle(paste("GO Biological Process:", ct)) +
           xlab("Gene Count") + ylab("") +
           theme(plot.background  = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 axis.text.y      = element_text(size = 8))

      ggsave(paste0("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/plots/GO_",
                    gsub("/| ", "_", ct), ".png"),
             p, width = 10, height = 8, bg = "white")

      cat("  GO terms found:", nrow(go_result), "\n")
    }
  }, error = function(e) cat("  GO error:", e$message, "\n"))

  # KEGG Enrichment
  tryCatch({
    kegg_result <- enrichKEGG(gene          = entrez$ENTREZID,
                               organism      = "hsa",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05)

    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
      write.csv(as.data.frame(kegg_result),
                paste0("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/tables/KEGG_",
                       gsub("/| ", "_", ct), ".csv"),
                row.names = FALSE)

      kegg_df <- as.data.frame(kegg_result) %>%
        head(20) %>%
        mutate(Description = factor(Description,
                                    levels = rev(Description)))

      p <- ggplot(kegg_df, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
           geom_point() +
           scale_color_gradient(low = "red", high = "blue") +
           theme_bw() +
           ggtitle(paste("KEGG Pathways:", ct)) +
           xlab("Gene Count") + ylab("") +
           theme(plot.background  = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 axis.text.y      = element_text(size = 8))

      ggsave(paste0("/N/scratch/dikonda/seurat_project/results/pathway_enrichment/plots/KEGG_",
                    gsub("/| ", "_", ct), ".png"),
             p, width = 10, height = 8, bg = "white")

      cat("  KEGG pathways found:", nrow(kegg_result), "\n")
    }
  }, error = function(e) cat("  KEGG error:", e$message, "\n"))
}

cat("\nPathway enrichment analysis complete!\n")
