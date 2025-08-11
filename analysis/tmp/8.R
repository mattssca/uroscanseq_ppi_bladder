library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)

# =============================================================================
# CREATE URO SUBTYPE VISUALIZATION WITH LUNDTAX HIGHLIGHTING
# =============================================================================

create_uro_lundtax_figure <- function(significant_hubs, hub_results, 
                                      lundtax_genes = LundTax2023Classifier::gene_list) {
  
  # Get LundTax gene symbols
  lundtax_symbols <- lundtax_genes$hgnc_symbol
  
  # Get Uro subtype data
  uro_significant <- significant_hubs[["uro"]]
  uro_all_hubs <- hub_results[["uro"]]
  
  # Prepare comprehensive Uro data
  uro_comprehensive <- uro_all_hubs %>%
    mutate(
      # Check if gene is in LundTax
      in_lundtax = gene %in% lundtax_symbols,
      # Check if gene is significant
      is_significant = gene %in% uro_significant$gene,
      # Create categories
      gene_category = case_when(
        is_significant & in_lundtax ~ "Significant + LundTax",
        is_significant & !in_lundtax ~ "Significant + Novel",
        !is_significant & in_lundtax ~ "LundTax Only",
        TRUE ~ "Other"
      ),
      # For labeling top genes
      label = case_when(
        is_significant & composite_hub_score <= quantile(composite_hub_score, 0.05, na.rm = TRUE) ~ gene,
        in_lundtax & composite_hub_score <= quantile(composite_hub_score, 0.1, na.rm = TRUE) ~ gene,
        TRUE ~ ""
      ),
      # Log transform for better visualization
      log_degree = log10(degree + 1),
      log_expr = log10(mean_expr + 1)
    ) %>%
    filter(!is.na(gene) & gene != "")
  
  # ==========================================================================
  # FIGURE 1: HUB SCORE vs EXPRESSION SCATTER PLOT
  # ==========================================================================
  
  fig1 <- ggplot(uro_comprehensive, aes(x = mean_expr, y = -composite_hub_score)) +
    geom_point(aes(color = gene_category, size = degree), alpha = 0.7) +
    scale_color_manual(values = c(
      "Significant + LundTax" = "red",
      "Significant + Novel" = "orange", 
      "LundTax Only" = "blue",
      "Other" = "lightgrey"
    )) +
    scale_size_continuous(range = c(0.5, 4), name = "Degree") +
    geom_text_repel(aes(label = label), 
                    size = 3, 
                    max.overlaps = 20,
                    box.padding = 0.3,
                    force = 2) +
    theme_minimal() +
    labs(
      title = "Uro Subtype: Hub Genes and LundTax Classifier Overlap",
      subtitle = paste("Total genes:", nrow(uro_comprehensive), 
                       "| LundTax genes:", sum(uro_comprehensive$in_lundtax),
                       "| Significant hubs:", sum(uro_comprehensive$is_significant)),
      x = "Mean Expression",
      y = "Hub Score (negative, higher = better hub)",
      color = "Gene Category"
    ) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom"
    )
  
  # ==========================================================================
  # FIGURE 2: DEGREE vs EXPRESSION WITH LUNDTAX HIGHLIGHTING
  # ==========================================================================
  
  fig2 <- ggplot(uro_comprehensive, aes(x = log_expr, y = log_degree)) +
    geom_point(aes(color = gene_category, alpha = gene_category), size = 2) +
    scale_color_manual(values = c(
      "Significant + LundTax" = "red",
      "Significant + Novel" = "orange", 
      "LundTax Only" = "blue",
      "Other" = "lightgrey"
    )) +
    scale_alpha_manual(values = c(
      "Significant + LundTax" = 1,
      "Significant + Novel" = 1, 
      "LundTax Only" = 0.8,
      "Other" = 0.3
    )) +
    geom_text_repel(aes(label = label), 
                    size = 3, 
                    max.overlaps = 25,
                    box.padding = 0.3) +
    theme_minimal() +
    labs(
      title = "Uro Subtype: Network Degree vs Expression",
      x = "Log10(Expression + 1)",
      y = "Log10(Degree + 1)",
      color = "Gene Category",
      alpha = "Gene Category"
    ) +
    guides(alpha = "none") +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "bottom"
    )
  
  # ==========================================================================
  # FIGURE 3: BAR PLOT OF TOP GENES WITH LUNDTAX STATUS
  # ==========================================================================
  
  # Get top 30 hub genes
  top_uro_genes <- uro_comprehensive %>%
    arrange(composite_hub_score) %>%
    head(30) %>%
    mutate(
      gene_rank = 1:n(),
      gene_ordered = factor(gene, levels = rev(gene))
    )
  
  fig3 <- ggplot(top_uro_genes, aes(x = gene_ordered, y = degree, fill = gene_category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Significant + LundTax" = "red",
      "Significant + Novel" = "orange", 
      "LundTax Only" = "blue",
      "Other" = "lightgrey"
    )) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Top 30 Uro Hub Genes by Degree",
      subtitle = "Colored by LundTax and Significance Status",
      x = "Gene",
      y = "Network Degree",
      fill = "Gene Category"
    ) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
  
  # ==========================================================================
  # FIGURE 4: VENN DIAGRAM-STYLE SUMMARY
  # ==========================================================================
  
  # Create summary counts
  summary_counts <- uro_comprehensive %>%
    group_by(gene_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      percentage = round(count / sum(count) * 100, 1),
      label = paste0(gene_category, "\n(", count, " genes, ", percentage, "%)")
    )
  
  fig4 <- ggplot(summary_counts, aes(x = "", y = count, fill = gene_category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c(
      "Significant + LundTax" = "red",
      "Significant + Novel" = "orange", 
      "LundTax Only" = "blue",
      "Other" = "lightgrey"
    )) +
    theme_void() +
    labs(
      title = "Uro Subtype Gene Categories",
      fill = "Category"
    ) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "right"
    ) +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold")
  
  # ==========================================================================
  # SAVE ALL FIGURES
  # ==========================================================================
  
  ggsave("figs/uro_lundtax_hub_score_vs_expression.pdf", fig1, width = 12, height = 8)
  ggsave("figs/uro_lundtax_degree_vs_expression.pdf", fig2, width = 10, height = 8)
  ggsave("figs/uro_lundtax_top30_genes.pdf", fig3, width = 10, height = 12)
  ggsave("figs/uro_lundtax_summary_pie.pdf", fig4, width = 8, height = 6)
  
  # ==========================================================================
  # CREATE COMBINED FIGURE
  # ==========================================================================
  
  library(patchwork)
  
  combined_fig <- (fig1 | fig2) / (fig3 | fig4)
  
  ggsave("figs/uro_lundtax_comprehensive_figure.pdf", combined_fig, width = 20, height = 16)
  
  # ==========================================================================
  # PRINT DETAILED STATISTICS
  # ==========================================================================
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("URO SUBTYPE - LUNDTAX COMPARISON DETAILED REPORT\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("Total genes in Uro network:", nrow(uro_comprehensive), "\n")
  cat("Genes in LundTax classifier:", sum(uro_comprehensive$in_lundtax), "\n")
  cat("Significant hub genes:", sum(uro_comprehensive$is_significant), "\n\n")
  
  # Category breakdown
  cat("GENE CATEGORIES:\n")
  for(i in 1:nrow(summary_counts)) {
    cat(sprintf("  %s: %d genes (%.1f%%)\n", 
                summary_counts$gene_category[i], 
                summary_counts$count[i], 
                summary_counts$percentage[i]))
  }
  
  # Top LundTax genes
  lundtax_in_uro <- uro_comprehensive %>%
    filter(in_lundtax) %>%
    arrange(composite_hub_score) %>%
    head(10)
  
  cat("\nTOP 10 URO HUBS ALREADY IN LUNDTAX:\n")
  for(i in 1:nrow(lundtax_in_uro)) {
    cat(sprintf("  %s (degree=%d, expr=%.2f, hub_score=%.2f)\n",
                lundtax_in_uro$gene[i], lundtax_in_uro$degree[i],
                lundtax_in_uro$mean_expr[i], lundtax_in_uro$composite_hub_score[i]))
  }
  
  # Top novel significant genes
  novel_significant <- uro_comprehensive %>%
    filter(is_significant & !in_lundtax) %>%
    arrange(composite_hub_score) %>%
    head(10)
  
  if(nrow(novel_significant) > 0) {
    cat("\nTOP 10 NOVEL SIGNIFICANT URO HUBS (NOT in LundTax):\n")
    for(i in 1:nrow(novel_significant)) {
      cat(sprintf("  %s (degree=%d, expr=%.2f, hub_score=%.2f)\n",
                  novel_significant$gene[i], novel_significant$degree[i],
                  novel_significant$mean_expr[i], novel_significant$composite_hub_score[i]))
    }
  }
  
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  return(list(
    fig1 = fig1,
    fig2 = fig2, 
    fig3 = fig3,
    fig4 = fig4,
    combined_fig = combined_fig,
    uro_data = uro_comprehensive,
    summary_counts = summary_counts
  ))
}

# =============================================================================
# RUN URO VISUALIZATION
# =============================================================================

cat("Creating comprehensive Uro subtype visualization with LundTax highlighting...\n")

uro_lundtax_plots <- create_uro_lundtax_figure(
  significant_hubs = significant_hubs,
  hub_results = hub_results
)

cat("\nUro visualization complete! Files created:\n")
cat("1. uro_lundtax_hub_score_vs_expression.pdf - Hub score vs expression scatter\n")
cat("2. uro_lundtax_degree_vs_expression.pdf - Network degree vs expression\n")
cat("3. uro_lundtax_top30_genes.pdf - Top 30 hub genes bar plot\n")
cat("4. uro_lundtax_summary_pie.pdf - Category summary pie chart\n")
cat("5. uro_lundtax_comprehensive_figure.pdf - All plots combined\n")
