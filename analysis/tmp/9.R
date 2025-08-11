library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# CREATE COMPREHENSIVE LUNDTAX VISUALIZATION FOR ALL SUBTYPES
# =============================================================================

create_subtype_lundtax_figure <- function(subtype_name, significant_hubs, hub_results, 
                                          lundtax_genes = LundTax2023Classifier::gene_list) {
  
  # Get LundTax gene symbols
  lundtax_symbols <- lundtax_genes$hgnc_symbol
  
  # Get subtype data
  subtype_significant <- significant_hubs[[subtype_name]]
  subtype_all_hubs <- hub_results[[subtype_name]]
  
  # Prepare comprehensive subtype data
  subtype_comprehensive <- subtype_all_hubs %>%
    mutate(
      # Check if gene is in LundTax
      in_lundtax = gene %in% lundtax_symbols,
      # Check if gene is significant
      is_significant = gene %in% subtype_significant$gene,
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
  
  fig1 <- ggplot(subtype_comprehensive, aes(x = mean_expr, y = -composite_hub_score)) +
    geom_point(aes(color = gene_category, size = degree), alpha = 0.7) +
    scale_color_manual(values = c(
      "Significant + LundTax" = "red",
      "Significant + Novel" = "orange", 
      "LundTax Only" = "blue",
      "Other" = "lightgrey"
    )) +
    scale_size_continuous(range = c(0.5, 4), name = "Degree") +
    geom_text_repel(aes(label = label), 
                    size = 2.5, 
                    max.overlaps = 15,
                    box.padding = 0.3,
                    force = 2) +
    theme_minimal() +
    labs(
      title = paste(toupper(subtype_name), "Subtype: Hub Genes and LundTax Classifier Overlap"),
      subtitle = paste("Total genes:", nrow(subtype_comprehensive), 
                       "| LundTax genes:", sum(subtype_comprehensive$in_lundtax),
                       "| Significant hubs:", sum(subtype_comprehensive$is_significant)),
      x = "Mean Expression",
      y = "Hub Score (negative, higher = better hub)",
      color = "Gene Category"
    ) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  # ==========================================================================
  # FIGURE 2: DEGREE vs EXPRESSION WITH LUNDTAX HIGHLIGHTING
  # ==========================================================================
  
  fig2 <- ggplot(subtype_comprehensive, aes(x = log_expr, y = log_degree)) +
    geom_point(aes(color = gene_category, alpha = gene_category), size = 1.5) +
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
                    size = 2.5, 
                    max.overlaps = 15,
                    box.padding = 0.3) +
    theme_minimal() +
    labs(
      title = paste(toupper(subtype_name), "Subtype: Network Degree vs Expression"),
      x = "Log10(Expression + 1)",
      y = "Log10(Degree + 1)",
      color = "Gene Category",
      alpha = "Gene Category"
    ) +
    guides(alpha = "none") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  # ==========================================================================
  # FIGURE 3: BAR PLOT OF TOP GENES WITH LUNDTAX STATUS
  # ==========================================================================
  
  # Get top 20 hub genes (reduced for better visualization)
  top_genes <- subtype_comprehensive %>%
    arrange(composite_hub_score) %>%
    head(20) %>%
    mutate(
      gene_rank = 1:n(),
      gene_ordered = factor(gene, levels = rev(gene))
    )
  
  fig3 <- ggplot(top_genes, aes(x = gene_ordered, y = degree, fill = gene_category)) +
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
      title = paste("Top 20", toupper(subtype_name), "Hub Genes"),
      x = "Gene",
      y = "Network Degree",
      fill = "Gene Category"
    ) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  # ==========================================================================
  # FIGURE 4: SUMMARY PIE CHART
  # ==========================================================================
  
  # Create summary counts
  summary_counts <- subtype_comprehensive %>%
    group_by(gene_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      percentage = round(count / sum(count) * 100, 1),
      label = paste0(count, "\n(", percentage, "%)")
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
      title = paste(toupper(subtype_name), "Gene Categories"),
      fill = "Category"
    ) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    ) +
    geom_text(aes(label = label), 
              position = position_stack(vjust = 0.5),
              size = 3, fontface = "bold")
  
  return(list(
    fig1 = fig1,
    fig2 = fig2, 
    fig3 = fig3,
    fig4 = fig4,
    subtype_data = subtype_comprehensive,
    summary_counts = summary_counts
  ))
}

# =============================================================================
# GENERATE DETAILED REPORT FOR EACH SUBTYPE
# =============================================================================

generate_subtype_report <- function(subtype_name, subtype_data, summary_counts) {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat(toupper(subtype_name), "SUBTYPE - LUNDTAX COMPARISON DETAILED REPORT\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  cat("Total genes in", toupper(subtype_name), "network:", nrow(subtype_data), "\n")
  cat("Genes in LundTax classifier:", sum(subtype_data$in_lundtax), 
      "(", round(sum(subtype_data$in_lundtax)/nrow(subtype_data)*100, 1), "%)\n")
  cat("Significant hub genes:", sum(subtype_data$is_significant), 
      "(", round(sum(subtype_data$is_significant)/nrow(subtype_data)*100, 1), "%)\n\n")
  
  # Category breakdown
  cat("GENE CATEGORIES:\n")
  for(i in 1:nrow(summary_counts)) {
    cat(sprintf("  %-25s: %4d genes (%5.1f%%)\n", 
                summary_counts$gene_category[i], 
                summary_counts$count[i], 
                summary_counts$percentage[i]))
  }
  
  # Top LundTax genes
  lundtax_in_subtype <- subtype_data %>%
    filter(in_lundtax) %>%
    arrange(composite_hub_score) %>%
    head(5)
  
  if(nrow(lundtax_in_subtype) > 0) {
    cat("\nTOP 5", toupper(subtype_name), "HUBS ALREADY IN LUNDTAX:\n")
    for(i in 1:nrow(lundtax_in_subtype)) {
      cat(sprintf("  %-12s (degree=%3d, expr=%6.2f, hub_score=%6.2f)\n",
                  lundtax_in_subtype$gene[i], lundtax_in_subtype$degree[i],
                  lundtax_in_subtype$mean_expr[i], lundtax_in_subtype$composite_hub_score[i]))
    }
  }
  
  # Top novel significant genes
  novel_significant <- subtype_data %>%
    filter(is_significant & !in_lundtax) %>%
    arrange(composite_hub_score) %>%
    head(5)
  
  if(nrow(novel_significant) > 0) {
    cat("\nTOP 5 NOVEL SIGNIFICANT", toupper(subtype_name), "HUBS (NOT in LundTax):\n")
    for(i in 1:nrow(novel_significant)) {
      cat(sprintf("  %-12s (degree=%3d, expr=%6.2f, hub_score=%6.2f)\n",
                  novel_significant$gene[i], novel_significant$degree[i],
                  novel_significant$mean_expr[i], novel_significant$composite_hub_score[i]))
    }
  } else {
    cat("\nNo novel significant hubs found for", toupper(subtype_name), "subtype.\n")
  }
  
  # Novel genes with high hub scores
  novel_high_hubs <- subtype_data %>%
    filter(!in_lundtax & composite_hub_score <= quantile(composite_hub_score, 0.1, na.rm = TRUE)) %>%
    arrange(composite_hub_score) %>%
    head(5)
  
  if(nrow(novel_high_hubs) > 0) {
    cat("\nTOP 5 NOVEL HIGH-SCORING", toupper(subtype_name), "HUBS (top 10% by hub score):\n")
    for(i in 1:nrow(novel_high_hubs)) {
      cat(sprintf("  %-12s (degree=%3d, expr=%6.2f, hub_score=%6.2f)\n",
                  novel_high_hubs$gene[i], novel_high_hubs$degree[i],
                  novel_high_hubs$mean_expr[i], novel_high_hubs$composite_hub_score[i]))
    }
  }
  
  cat(paste(rep("=", 80), collapse=""), "\n")
}

# =============================================================================
# CREATE COMPREHENSIVE COMPARISON ACROSS ALL SUBTYPES
# =============================================================================

create_all_subtypes_comparison <- function(significant_hubs, hub_results, 
                                           lundtax_genes = LundTax2023Classifier::gene_list) {
  
  all_results <- list()
  
  # Create individual subtype reports
  for(subtype in names(hub_results)) {
    cat("Processing", toupper(subtype), "subtype...\n")
    
    # Create figures
    subtype_results <- create_subtype_lundtax_figure(subtype, significant_hubs, hub_results, lundtax_genes)
    
    # Save individual plots
    ggsave(paste0("figs/", subtype, "_lundtax_hub_score_vs_expression.pdf"), 
           subtype_results$fig1, width = 10, height = 7)
    ggsave(paste0("figs/", subtype, "_lundtax_degree_vs_expression.pdf"), 
           subtype_results$fig2, width = 8, height = 6)
    ggsave(paste0("figs/", subtype, "_lundtax_top20_genes.pdf"), 
           subtype_results$fig3, width = 8, height = 10)
    ggsave(paste0("figs/", subtype, "_lundtax_summary_pie.pdf"), 
           subtype_results$fig4, width = 6, height = 5)
    
    # Create combined figure for this subtype
    combined_fig <- (subtype_results$fig1 | subtype_results$fig2) / 
      (subtype_results$fig3 | subtype_results$fig4)
    ggsave(paste0("figs/", subtype, "_lundtax_comprehensive_figure.pdf"), 
           combined_fig, width = 16, height = 12)
    
    # Generate detailed report
    generate_subtype_report(subtype, subtype_results$subtype_data, subtype_results$summary_counts)
    
    all_results[[subtype]] <- subtype_results
  }
  
  # ==========================================================================
  # CREATE CROSS-SUBTYPE SUMMARY
  # ==========================================================================
  
  cat("\n", paste(rep("=", 100), collapse=""), "\n")
  cat("CROSS-SUBTYPE LUNDTAX COMPARISON SUMMARY\n")
  cat(paste(rep("=", 100), collapse=""), "\n")
  
  # Collect summary statistics
  cross_summary <- data.frame()
  
  for(subtype in names(all_results)) {
    subtype_data <- all_results[[subtype]]$subtype_data
    
    cross_summary <- rbind(cross_summary, data.frame(
      subtype = toupper(subtype),
      total_genes = nrow(subtype_data),
      lundtax_genes = sum(subtype_data$in_lundtax),
      significant_hubs = sum(subtype_data$is_significant),
      sig_and_lundtax = sum(subtype_data$is_significant & subtype_data$in_lundtax),
      sig_and_novel = sum(subtype_data$is_significant & !subtype_data$in_lundtax),
      lundtax_only = sum(!subtype_data$is_significant & subtype_data$in_lundtax),
      pct_lundtax = round(sum(subtype_data$in_lundtax)/nrow(subtype_data)*100, 1),
      pct_significant = round(sum(subtype_data$is_significant)/nrow(subtype_data)*100, 1),
      pct_novel_of_sig = round(sum(subtype_data$is_significant & !subtype_data$in_lundtax)/
                                 sum(subtype_data$is_significant)*100, 1)
    ))
  }
  
  # Print summary table
  cat("SUMMARY TABLE:\n")
  cat(sprintf("%-8s %5s %8s %8s %8s %8s %8s %7s %7s %7s\n",
              "Subtype", "Total", "LundTax", "Sig.Hubs", "Sig+LT", "Sig+Nov", "LT.Only", "%LT", "%Sig", "%Nov"))
  cat(paste(rep("-", 95), collapse=""), "\n")
  
  for(i in 1:nrow(cross_summary)) {
    cat(sprintf("%-8s %5d %8d %8d %8d %8d %8d %6.1f%% %6.1f%% %6.1f%%\n",
                cross_summary$subtype[i], cross_summary$total_genes[i],
                cross_summary$lundtax_genes[i], cross_summary$significant_hubs[i],
                cross_summary$sig_and_lundtax[i], cross_summary$sig_and_novel[i],
                cross_summary$lundtax_only[i], cross_summary$pct_lundtax[i],
                cross_summary$pct_significant[i], cross_summary$pct_novel_of_sig[i]))
  }
  
  # ==========================================================================
  # CREATE CROSS-SUBTYPE COMPARISON PLOTS
  # ==========================================================================
  
  # Plot 1: LundTax overlap percentages
  overlap_plot <- ggplot(cross_summary, aes(x = subtype, y = pct_lundtax, fill = subtype)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = "Percentage of Genes in LundTax Classifier by Subtype",
      x = "Subtype",
      y = "Percentage in LundTax (%)"
    ) +
    theme(legend.position = "none")
  
  # Plot 2: Novel discovery rates
  novel_plot <- ggplot(cross_summary, aes(x = subtype, y = pct_novel_of_sig, fill = subtype)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = "Percentage of Significant Hubs that are Novel (Not in LundTax)",
      x = "Subtype", 
      y = "Percentage Novel (%)"
    ) +
    theme(legend.position = "none")
  
  # Plot 3: Absolute numbers comparison
  comparison_long <- cross_summary %>%
    select(subtype, sig_and_lundtax, sig_and_novel, lundtax_only) %>%
    pivot_longer(-subtype, names_to = "category", values_to = "count") %>%
    mutate(category = factor(category, 
                             levels = c("sig_and_lundtax", "sig_and_novel", "lundtax_only"),
                             labels = c("Significant + LundTax", "Significant + Novel", "LundTax Only")))
  
  comparison_plot <- ggplot(comparison_long, aes(x = subtype, y = count, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("red", "orange", "blue")) +
    theme_minimal() +
    labs(
      title = "Gene Categories Across Subtypes",
      x = "Subtype",
      y = "Number of Genes",
      fill = "Category"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save cross-subtype plots
  ggsave("figs/cross_subtype_lundtax_overlap.pdf", overlap_plot, width = 10, height = 6)
  ggsave("figs/cross_subtype_novel_discovery.pdf", novel_plot, width = 10, height = 6)
  ggsave("figs/cross_subtype_gene_categories.pdf", comparison_plot, width = 12, height = 6)
  
  # Combined cross-subtype figure
  cross_combined <- overlap_plot / novel_plot / comparison_plot
  ggsave("figs/cross_subtype_comprehensive_comparison.pdf", cross_combined, width = 12, height = 15)
  
  return(list(
    all_results = all_results,
    cross_summary = cross_summary,
    overlap_plot = overlap_plot,
    novel_plot = novel_plot,
    comparison_plot = comparison_plot
  ))
}

# =============================================================================
# RUN COMPREHENSIVE ANALYSIS FOR ALL SUBTYPES
# =============================================================================

cat("Creating comprehensive LundTax comparison for all subtypes...\n")

all_subtype_results <- create_all_subtypes_comparison(
  significant_hubs = significant_hubs,
  hub_results = hub_results
)

# Save all results
save(all_subtype_results, file = "data/all_subtypes_lundtax_comparison.Rdata")

cat("\n", paste(rep("=", 100), collapse=""), "\n")
cat("ALL SUBTYPE ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 100), collapse=""), "\n")
cat("Individual subtype files created:\n")
for(subtype in names(hub_results)) {
  cat(sprintf("  %s_lundtax_*.pdf - Individual plots for %s subtype\n", subtype, toupper(subtype)))
}
cat("\nCross-subtype comparison files:\n")
cat("  cross_subtype_lundtax_overlap.pdf\n")
cat("  cross_subtype_novel_discovery.pdf\n") 
cat("  cross_subtype_gene_categories.pdf\n")
cat("  cross_subtype_comprehensive_comparison.pdf\n")
cat("\nDetailed console reports printed above for each subtype.\n")
cat(paste(rep("=", 100), collapse=""), "\n")
