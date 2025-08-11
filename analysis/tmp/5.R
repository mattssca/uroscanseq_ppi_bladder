library(dplyr)
library(ggplot2)
library(pheatmap)

# =============================================================================
# IDENTIFY SUBTYPE-SPECIFIC HUBS
# =============================================================================

find_subtype_specific_hubs <- function(significant_hubs, min_evidence = 1) {
  
  # Get all significant genes across all subtypes
  all_sig_genes <- unique(unlist(lapply(significant_hubs, function(x) x$gene)))
  all_sig_genes <- all_sig_genes[!is.na(all_sig_genes)]
  
  # Create presence matrix
  presence_matrix <- matrix(0, nrow = length(all_sig_genes), ncol = length(significant_hubs))
  rownames(presence_matrix) <- all_sig_genes
  colnames(presence_matrix) <- names(significant_hubs)
  
  # Fill presence matrix with evidence level
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      for(i in 1:nrow(results)) {
        gene <- results$gene[i]
        if(!is.na(gene) && gene != "" && results$n_methods_significant[i] >= min_evidence) {
          presence_matrix[gene, subtype] <- results$n_methods_significant[i]
        }
      }
    }
  }
  
  # Count how many subtypes each gene is significant in
  gene_subtype_counts <- rowSums(presence_matrix > 0)
  
  # Identify different categories
  subtype_specific <- names(gene_subtype_counts[gene_subtype_counts == 1])
  shared_two <- names(gene_subtype_counts[gene_subtype_counts == 2])
  shared_three <- names(gene_subtype_counts[gene_subtype_counts == 3])
  shared_four <- names(gene_subtype_counts[gene_subtype_counts == 4])
  universal <- names(gene_subtype_counts[gene_subtype_counts == 5])
  
  # Get detailed information for subtype-specific hubs
  specific_details <- list()
  
  for(gene in subtype_specific) {
    # Find which subtype this gene is specific to
    subtype_idx <- which(presence_matrix[gene, ] > 0)
    if(length(subtype_idx) == 1) {
      subtype_name <- colnames(presence_matrix)[subtype_idx]
      
      # Get detailed info from significant_hubs
      gene_info <- significant_hubs[[subtype_name]][significant_hubs[[subtype_name]]$gene == gene, ]
      
      if(nrow(gene_info) > 0) {
        specific_details[[gene]] <- list(
          subtype = subtype_name,
          n_methods = gene_info$n_methods_significant[1],
          expr_fold_change = gene_info$expr_fold_change[1],
          target_expr = gene_info$target_expr[1],
          expr_fdr = gene_info$expr_fdr[1],
          target_rank = gene_info$target_rank[1]
        )
      }
    }
  }
  
  return(list(
    presence_matrix = presence_matrix,
    gene_counts = gene_subtype_counts,
    subtype_specific = subtype_specific,
    shared_two = shared_two,
    shared_three = shared_three,
    shared_four = shared_four,
    universal = universal,
    specific_details = specific_details
  ))
}

# =============================================================================
# ANALYZE SUBTYPE-SPECIFIC PATTERNS
# =============================================================================

analyze_specificity_patterns <- function(hub_results, significant_hubs) {
  
  # Find subtype-specific hubs
  specificity_analysis <- find_subtype_specific_hubs(significant_hubs, min_evidence = 1)
  
  # Print summary
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("SUBTYPE-SPECIFIC HUB ANALYSIS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("Total unique significant genes:", length(specificity_analysis$gene_counts), "\n")
  cat("Subtype-specific hubs (1 subtype):", length(specificity_analysis$subtype_specific), "\n")
  cat("Shared by 2 subtypes:", length(specificity_analysis$shared_two), "\n")
  cat("Shared by 3 subtypes:", length(specificity_analysis$shared_three), "\n")
  cat("Shared by 4 subtypes:", length(specificity_analysis$shared_four), "\n")
  cat("Universal hubs (all 5 subtypes):", length(specificity_analysis$universal), "\n")
  
  # Analyze each subtype's specific hubs
  cat("\n", paste(rep("-", 60), collapse=""), "\n")
  cat("SUBTYPE-SPECIFIC HUBS BY SUBTYPE:\n")
  cat(paste(rep("-", 60), collapse=""), "\n")
  
  subtype_specific_count <- list()
  
  for(subtype in names(significant_hubs)) {
    # Count specific hubs for this subtype
    specific_for_subtype <- names(specificity_analysis$specific_details)[
      sapply(specificity_analysis$specific_details, function(x) x$subtype == subtype)
    ]
    
    subtype_specific_count[[subtype]] <- length(specific_for_subtype)
    
    cat("\n", toupper(subtype), "SUBTYPE (", length(specific_for_subtype), "specific hubs):\n")
    
    if(length(specific_for_subtype) > 0) {
      # Show top 10 specific hubs
      specific_info <- specificity_analysis$specific_details[specific_for_subtype]
      
      # Sort by evidence level and fold change
      specific_sorted <- specific_info[order(
        -sapply(specific_info, function(x) x$n_methods),
        -sapply(specific_info, function(x) abs(log2(x$expr_fold_change)))
      )]
      
      for(i in 1:min(10, length(specific_sorted))) {
        gene <- names(specific_sorted)[i]
        info <- specific_sorted[[i]]
        cat(sprintf("  %s: FC=%.2f, methods=%d, expr=%.2f\n",
                    gene, info$expr_fold_change, info$n_methods, info$target_expr))
      }
      
      if(length(specific_sorted) > 10) {
        cat(sprintf("  ... and %d more\n", length(specific_sorted) - 10))
      }
    } else {
      cat("  No subtype-specific hubs found.\n")
    }
  }
  
  return(list(
    specificity_analysis = specificity_analysis,
    subtype_counts = subtype_specific_count
  ))
}

# =============================================================================
# VISUALIZE SUBTYPE SPECIFICITY
# =============================================================================

visualize_subtype_specificity <- function(specificity_results) {
  
  specificity_analysis <- specificity_results$specificity_analysis
  subtype_counts <- specificity_results$subtype_counts
  
  # 1. Presence/absence heatmap
  if(nrow(specificity_analysis$presence_matrix) > 0) {
    
    # Convert to binary for cleaner visualization
    binary_matrix <- ifelse(specificity_analysis$presence_matrix > 0, 1, 0)
    
    # Only show genes that are specific or shared by few subtypes
    interesting_genes <- names(specificity_analysis$gene_counts[specificity_analysis$gene_counts <= 3])
    
    if(length(interesting_genes) > 0) {
      interesting_matrix <- binary_matrix[interesting_genes, , drop = FALSE]
      
      pheatmap(interesting_matrix,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               color = c("white", "darkblue"),
               breaks = c(-0.5, 0.5, 1.5),
               legend_breaks = c(0, 1),
               legend_labels = c("Not Significant", "Significant Hub"),
               main = "Subtype-Specific and Rarely Shared Hubs",
               filename = "figs/subtype_specific_hubs.pdf",
               width = 8, height = max(6, length(interesting_genes) * 0.2))
    }
  }
  
  # 2. Bar plot of specificity counts
  sharing_summary <- data.frame(
    category = c("Subtype-specific", "Shared by 2", "Shared by 3", "Shared by 4", "Universal"),
    count = c(
      length(specificity_analysis$subtype_specific),
      length(specificity_analysis$shared_two),
      length(specificity_analysis$shared_three),
      length(specificity_analysis$shared_four),
      length(specificity_analysis$universal)
    )
  )
  
  sharing_plot <- ggplot(sharing_summary, aes(x = factor(category, levels = category), y = count)) +
    geom_bar(stat = "identity", fill = c("red", "orange", "yellow", "lightblue", "darkblue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Hub Sharing Patterns Across Subtypes",
      x = "Sharing Category",
      y = "Number of Genes"
    )
  
  ggsave("figs/hub_sharing_summary.pdf", sharing_plot, width = 10, height = 6)
  
  # 3. Subtype-specific counts
  subtype_df <- data.frame(
    subtype = names(subtype_counts),
    specific_hubs = unlist(subtype_counts)
  )
  
  subtype_plot <- ggplot(subtype_df, aes(x = subtype, y = specific_hubs, fill = subtype)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = "Number of Subtype-Specific Hubs",
      x = "Subtype",
      y = "Number of Specific Hubs"
    ) +
    theme(legend.position = "none")
  
  ggsave("figs/subtype_specific_counts.pdf", subtype_plot, width = 8, height = 6)
  
  return(list(sharing_plot = sharing_plot, subtype_plot = subtype_plot))
}

# =============================================================================
# DETAILED ANALYSIS OF SPECIFIC HUBS
# =============================================================================

detailed_specific_hub_analysis <- function(specificity_results, hub_results) {
  
  specificity_analysis <- specificity_results$specificity_analysis
  
  # Get expression data for all specific hubs across all subtypes
  all_hubs_data <- bind_rows(hub_results, .id = "subtype")
  
  specific_expression_analysis <- list()
  
  for(subtype in names(hub_results)) {
    # Get genes specific to this subtype
    specific_genes <- names(specificity_analysis$specific_details)[
      sapply(specificity_analysis$specific_details, function(x) x$subtype == subtype)
    ]
    
    if(length(specific_genes) > 0) {
      # Get expression data for these genes across all subtypes
      specific_expr_data <- data.frame()
      
      for(gene in specific_genes) {
        gene_data <- all_hubs_data[all_hubs_data$gene == gene, ]
        
        if(nrow(gene_data) > 0) {
          specific_expr_data <- rbind(specific_expr_data, data.frame(
            gene = gene,
            target_subtype = subtype,
            subtype = gene_data$subtype,
            expression = gene_data$mean_expr,
            degree = gene_data$degree,
            hub_score = gene_data$composite_hub_score
          ))
        }
      }
      
      specific_expression_analysis[[subtype]] <- specific_expr_data
    }
  }
  
  return(specific_expression_analysis)
}

# =============================================================================
# RUN COMPLETE SUBTYPE SPECIFICITY ANALYSIS
# =============================================================================

cat("Running subtype-specific hub analysis...\n")

# Main analysis
specificity_results <- analyze_specificity_patterns(hub_results, significant_hubs)

# Create visualizations
specificity_plots <- visualize_subtype_specificity(specificity_results)

# Detailed expression analysis
detailed_analysis <- detailed_specific_hub_analysis(specificity_results, hub_results)

# Save results
save(specificity_results, detailed_analysis, file = "data/subtype_specific_analysis.Rdata")

# =============================================================================
# EXPORT SUBTYPE-SPECIFIC HUB LISTS
# =============================================================================

# Create summary tables for each subtype
for(subtype in names(significant_hubs)) {
  specific_genes <- names(specificity_results$specificity_analysis$specific_details)[
    sapply(specificity_results$specificity_analysis$specific_details, function(x) x$subtype == subtype)
  ]
  
  if(length(specific_genes) > 0) {
    # Create detailed table
    specific_table <- data.frame()
    
    for(gene in specific_genes) {
      info <- specificity_results$specificity_analysis$specific_details[[gene]]
      specific_table <- rbind(specific_table, data.frame(
        gene = gene,
        subtype = info$subtype,
        n_methods_significant = info$n_methods,
        expr_fold_change = info$expr_fold_change,
        target_expression = info$target_expr,
        expr_fdr = info$expr_fdr,
        target_rank = info$target_rank
      ))
    }
    
    # Sort by evidence and effect size
    specific_table <- specific_table %>%
      arrange(desc(n_methods_significant), desc(abs(log2(expr_fold_change))))
    
    # Save to CSV
    write.csv(specific_table, file = paste0("data/", subtype, "_specific_hubs.csv"), row.names = FALSE)
    
    cat("Saved", nrow(specific_table), "specific hubs for", toupper(subtype), "subtype\n")
  }
}

cat("\nSubtype-specific analysis complete!\n")
cat("Check data/ directory for individual subtype hub lists (.csv files)\n")
cat("Check figs/ directory for visualization plots\n")
