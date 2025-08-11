library(dplyr)

# =============================================================================
# COMPARE SUBTYPE-SPECIFIC GENES WITH LUNDTAX GENE LIST
# =============================================================================

compare_with_lundtax <- function(significant_hubs, specificity_results = NULL, 
                                 lundtax_genes = LundTax2023Classifier::gene_list) {
  
  # Get LundTax gene symbols
  lundtax_symbols <- lundtax_genes$hgnc_symbol
  
  cat("LundTax2023 classifier contains", length(lundtax_symbols), "genes\n\n")
  
  # ==========================================================================
  # 1. CHECK SUBTYPE-SPECIFIC GENES (if specificity analysis was run)
  # ==========================================================================
  
  if(!is.null(specificity_results)) {
    cat("=== SUBTYPE-SPECIFIC GENES vs LUNDTAX ===\n")
    
    specificity_analysis <- specificity_results$specificity_analysis
    
    for(subtype in names(hub_results)) {
      # Get genes specific to this subtype
      specific_genes <- names(specificity_analysis$specific_details)[
        sapply(specificity_analysis$specific_details, function(x) x$subtype == subtype)
      ]
      
      if(length(specific_genes) > 0) {
        # Check overlap with LundTax
        in_lundtax <- specific_genes %in% lundtax_symbols
        novel_specific <- specific_genes[!in_lundtax]
        known_specific <- specific_genes[in_lundtax]
        
        cat("\n", toupper(subtype), "subtype-specific genes:\n")
        cat("  Total specific genes:", length(specific_genes), "\n")
        cat("  Already in LundTax:", length(known_specific), "(", 
            round(length(known_specific)/length(specific_genes)*100, 1), "%)\n")
        cat("  Novel specific genes:", length(novel_specific), "(", 
            round(length(novel_specific)/length(specific_genes)*100, 1), "%)\n")
        
        if(length(known_specific) > 0) {
          cat("  Known genes:", paste(head(known_specific, 5), collapse = ", "))
          if(length(known_specific) > 5) cat(", ...")
          cat("\n")
        }
        
        if(length(novel_specific) > 0) {
          cat("  Novel genes:", paste(head(novel_specific, 5), collapse = ", "))
          if(length(novel_specific) > 5) cat(", ...")
          cat("\n")
        }
      } else {
        cat("\n", toupper(subtype), "subtype: No specific genes found\n")
      }
    }
  }
  
  # ==========================================================================
  # 2. CHECK ALL SIGNIFICANT HUBS vs LUNDTAX
  # ==========================================================================
  
  cat("\n\n=== ALL SIGNIFICANT HUBS vs LUNDTAX ===\n")
  
  all_novel_genes <- list()
  all_known_genes <- list()
  
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      # Get all significant genes for this subtype
      sig_genes <- results$gene[!is.na(results$gene)]
      
      # Check overlap with LundTax
      in_lundtax <- sig_genes %in% lundtax_symbols
      novel_genes <- sig_genes[!in_lundtax]
      known_genes <- sig_genes[in_lundtax]
      
      all_novel_genes[[subtype]] <- novel_genes
      all_known_genes[[subtype]] <- known_genes
      
      cat("\n", toupper(subtype), "significant hubs:\n")
      cat("  Total significant genes:", length(sig_genes), "\n")
      cat("  Already in LundTax:", length(known_genes), "(", 
          round(length(known_genes)/length(sig_genes)*100, 1), "%)\n")
      cat("  Novel significant genes:", length(novel_genes), "(", 
          round(length(novel_genes)/length(sig_genes)*100, 1), "%)\n")
    } else {
      cat("\n", toupper(subtype), "subtype: No significant genes found\n")
      all_novel_genes[[subtype]] <- character(0)
      all_known_genes[[subtype]] <- character(0)
    }
  }
  
  # ==========================================================================
  # 3. GET TOP NOVEL GENES WITH DETAILS
  # ==========================================================================
  
  cat("\n\n=== TOP NOVEL SIGNIFICANT HUBS (NOT in LundTax) ===\n")
  
  novel_with_details <- list()
  
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    novel_genes <- all_novel_genes[[subtype]]
    
    if(length(novel_genes) > 0) {
      # Get details for novel genes
      novel_details <- results[results$gene %in% novel_genes, ] %>%
        arrange(desc(n_methods_significant), expr_fdr)
      
      novel_with_details[[subtype]] <- novel_details
      
      cat("\n", toupper(subtype), "top novel genes:\n")
      
      for(i in 1:min(10, nrow(novel_details))) {
        gene_info <- novel_details[i, ]
        cat(sprintf("  %s: methods=%d, FC=%.2f, expr=%.2f\n",
                    gene_info$gene, gene_info$n_methods_significant,
                    gene_info$expr_fold_change, gene_info$target_expr))
      }
      
      if(nrow(novel_details) > 10) {
        cat(sprintf("  ... and %d more novel genes\n", nrow(novel_details) - 10))
      }
    }
  }
  
  return(list(
    lundtax_genes = lundtax_symbols,
    all_novel_genes = all_novel_genes,
    all_known_genes = all_known_genes,
    novel_with_details = novel_with_details
  ))
}

# =============================================================================
# EXPORT NOVEL GENES TO FILES
# =============================================================================

export_novel_genes <- function(comparison_results, output_dir = "novel_hub_genes") {
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Exporting novel hub genes to", output_dir, "directory...\n")
  
  for(subtype in names(comparison_results$all_novel_genes)) {
    novel_genes <- comparison_results$all_novel_genes[[subtype]]
    
    if(length(novel_genes) > 0) {
      # Export simple gene list
      filename_simple <- file.path(output_dir, paste0(subtype, "_novel_hubs.txt"))
      writeLines(novel_genes, filename_simple)
      
      # Export detailed information
      novel_details <- comparison_results$novel_with_details[[subtype]]
      if(!is.null(novel_details) && nrow(novel_details) > 0) {
        filename_detailed <- file.path(output_dir, paste0(subtype, "_novel_hubs_detailed.csv"))
        write.csv(novel_details, filename_detailed, row.names = FALSE)
      }
      
      cat("  Exported", length(novel_genes), "novel genes for", toupper(subtype), "\n")
    } else {
      cat("  No novel genes for", toupper(subtype), "\n")
    }
  }
}

# =============================================================================
# CREATE SUMMARY VISUALIZATION
# =============================================================================

visualize_novel_vs_known <- function(comparison_results) {
  
  # Create summary data
  summary_data <- data.frame()
  
  for(subtype in names(comparison_results$all_novel_genes)) {
    n_novel <- length(comparison_results$all_novel_genes[[subtype]])
    n_known <- length(comparison_results$all_known_genes[[subtype]])
    n_total <- n_novel + n_known
    
    summary_data <- rbind(summary_data, data.frame(
      subtype = subtype,
      novel = n_novel,
      known = n_known,
      total = n_total,
      pct_novel = ifelse(n_total > 0, n_novel/n_total*100, 0)
    ))
  }
  
  # Reshape for plotting
  summary_long <- summary_data %>%
    select(subtype, novel, known) %>%
    pivot_longer(-subtype, names_to = "gene_type", values_to = "count")
  
  # Create plot
  novel_plot <- ggplot(summary_long, aes(x = subtype, y = count, fill = gene_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("known" = "lightblue", "novel" = "red"),
                      labels = c("Known (in LundTax)", "Novel (not in LundTax)")) +
    theme_minimal() +
    labs(
      title = "Novel vs Known Significant Hub Genes",
      subtitle = "Comparison with LundTax2023 Classifier Gene List",
      x = "Subtype",
      y = "Number of Significant Hub Genes",
      fill = "Gene Status"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("figs/novel_vs_known_hubs.pdf", novel_plot, width = 10, height = 6)
  
  # Create percentage plot
  pct_plot <- ggplot(summary_data, aes(x = subtype, y = pct_novel, fill = subtype)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = "Percentage of Novel Hub Genes by Subtype",
      x = "Subtype",
      y = "Percentage Novel (%)"
    ) +
    theme(legend.position = "none")
  
  ggsave("figs/novel_percentage_by_subtype.pdf", pct_plot, width = 8, height = 6)
  
  return(list(novel_plot = novel_plot, pct_plot = pct_plot))
}

# =============================================================================
# RUN COMPARISON ANALYSIS
# =============================================================================

cat("Comparing hub genes with LundTax2023 classifier...\n\n")

# Run the comparison
comparison_results <- compare_with_lundtax(
  significant_hubs = significant_hubs,
  specificity_results = if(exists("specificity_results")) specificity_results else NULL
)

# Export novel genes
export_novel_genes(comparison_results)

# Create visualizations
comparison_plots <- visualize_novel_vs_known(comparison_results)

# Save results
save(comparison_results, file = "data/lundtax_comparison_results.Rdata")

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("LUNDTAX COMPARISON SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Analysis complete! Check:\n")
cat("1. novel_hub_genes/ directory for novel gene lists\n")
cat("2. figs/ directory for comparison plots\n")
cat("3. Console output above for detailed statistics\n")
cat(paste(rep("=", 60), collapse=""), "\n")
