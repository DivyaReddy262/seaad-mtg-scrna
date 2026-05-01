.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(speckle)
library(ggplot2)
library(dplyr)

cat("Loading annotated object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds")

# Cell proportion analysis
cat("Running cell proportion analysis...\n")
props <- getTransformedProps(
  clusters = seurat_obj@meta.data$subclass,
  sample   = seurat_obj@meta.data$donor_id,
  transform = "logit"
)

# Plot cell proportions by cognitive status
meta <- seurat_obj@meta.data %>%
  select(donor_id, cogn_status, subclass) %>%
  distinct(donor_id, .keep_all = TRUE)

# Count cells per subclass per cognitive status
prop_df <- seurat_obj@meta.data %>%
  group_by(cogn_status, subclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cogn_status) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Bar plot of proportions
p1 <- ggplot(prop_df, aes(x = cogn_status, y = proportion,
                           fill = subclass)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_bw() +
      ggtitle("Cell Type Proportions by Cognitive Status") +
      xlab("Cognitive Status") +
      ylab("Proportion") +
      theme(plot.background  = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            legend.text      = element_text(size = 8)) +
      guides(fill = guide_legend(ncol = 1))

ggsave("/N/scratch/dikonda/seurat_project/results/cell_proportion_barplot.png",
       p1, width = 10, height = 8, bg = "white")

# Dot plot showing proportion changes
p2 <- ggplot(prop_df, aes(x = cogn_status, y = proportion,
                           color = subclass, group = subclass)) +
      geom_point(size = 3) +
      geom_line() +
      theme_bw() +
      ggtitle("Cell Type Proportion Changes") +
      xlab("Cognitive Status") +
      ylab("Proportion") +
      theme(plot.background  = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/cell_proportion_dotplot.png",
       p2, width = 10, height = 8, bg = "white")

write.csv(prop_df,
          "/N/scratch/dikonda/seurat_project/results/cell_proportions.csv",
          row.names = FALSE)
cat("Cell proportion analysis done!\n")
