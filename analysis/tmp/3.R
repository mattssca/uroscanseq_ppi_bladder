library(igraph)
library(dplyr)
library(broom)

# =============================================================================
# FULLY FIXED OUTLIER-BASED SIGNIFICANCE TEST
# =============================================================================

outlier_based_hub_test_fixed <- function(hub_results, outlier_threshold = 2) {
  
  # Combine all hub data
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  outlier_significance <- list()
  
  for(test_subtype in names(hub_results)) {
    
    # Get genes present in multiple subtypes
    multi_subtype_genes <- all_hubs %>%
      group_by(gene) %>%
      filter(n() >= 4) %>%  # Present in at least 4 subtypes (need 3+ for comparison)
      pull(gene) %>%
      unique()
    
    outlier_tests <- data.frame()
    
    cat("  Processing", length(multi_subtype_genes), "genes for", test_subtype, "subtype...\n")
    
    for(gene in multi_subtype_genes) {
      gene_data <- all_hubs[all_hubs$gene == gene, ]
      
      target_expr <- gene_data$mean_expr[gene_data$subtype == test_subtype]
      other_expr <- gene_data$mean_expr[gene_data$subtype != test_subtype]
      
      # Check we have valid data
      if(length(target_expr) == 1 && length(other_expr) >= 3 && 
         !is.na(target_expr) && all(!is.na(other_expr))) {
        
        # Calculate statistics for other subtypes
        other_mean <- mean(other_expr, na.rm = TRUE)
        other_sd <- sd(other_expr, na.rm = TRUE)
        
        # Check for valid standard deviation
        if(!is.na(other_sd) && other_sd > 0 && !is.na(other_mean)) {
          
          # Calculate Z-score relative to other subtypes
          z_score <- (target_expr - other_mean) / other_sd
          
          # Check for valid z-score
          if(!is.na(z_score) && is.finite(z_score)) {
            
            # Is this an outlier?
            is_outlier <- abs(z_score) > outlier_threshold
            
            # Calculate approximate p-value
            p_value <- 2 * (1 - pnorm(abs(z_score)))
            
            # Calculate fold change
            expr_fold_change <- if(other_mean > 0) target_expr / other_mean else NA
            
            outlier_tests <- rbind(outlier_tests, data.frame(
              gene = gene,
              subtype = test_subtype,
              target_expr = target_expr,
              other_expr_mean = other_mean,
              other_expr_sd = other_sd,
              z_score = z_score,
              is_outlier = is_outlier,
              outlier_pvalue = p_value,
              expr_fold_change = expr_fold_change,
              degree = gene_data$degree[gene_data$subtype == test_subtype],
              hub_score = gene_data$composite_hub_score[gene_data$subtype == test_subtype],
              n_other_subtypes = length(other_expr),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    # Adjust for multiple testing
    if(nrow(outlier_tests) > 0) {
      outlier_tests$outlier_fdr <- p.adjust(outlier_tests$outlier_pvalue, method = "fdr")
    }
    
    outlier_significance[[test_subtype]] <- outlier_tests
  }
  
  return(outlier_significance)
}

# =============================================================================
# ALTERNATIVE APPROACH: RANK-BASED SIGNIFICANCE TEST
# =============================================================================

rank_based_hub_test <- function(hub_results, top_percentile = 0.1) {
  
  # Combine all hub data
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  rank_significance <- list()
  
  for(test_subtype in names(hub_results)) {
    
    cat("  Running rank-based test for", test_subtype, "...\n")
    
    # Get all genes in this subtype
    subtype_data <- all_hubs[all_hubs$subtype == test_subtype, ]
    
    # For each gene, compare its rank in this subtype vs others
    rank_tests <- data.frame()
    
    for(gene in subtype_data$gene) {
      if(is.na(gene) || gene == "") next
      
      gene_all_data <- all_hubs[all_hubs$gene == gene, ]
      
      # Need gene in at least 3 subtypes
      if(nrow(gene_all_data) >= 3) {
        
        # Get expression ranks within each subtype
        gene_ranks <- data.frame()
        
        for(st in unique(gene_all_data$subtype)) {
          st_all_genes <- all_hubs[all_hubs$subtype == st, ]
          gene_expr <- gene_all_data$mean_expr[gene_all_data$subtype == st]
          
          if(length(gene_expr) == 1 && !is.na(gene_expr)) {
            # Calculate percentile rank
            percentile_rank <- mean(st_all_genes$mean_expr <= gene_expr, na.rm = TRUE)
            
            gene_ranks <- rbind(gene_ranks, data.frame(
              subtype = st,
              expression = gene_expr,
              percentile_rank = percentile_rank
            ))
          }
        }
        
        if(nrow(gene_ranks) >= 3) {
          target_rank <- gene_ranks$percentile_rank[gene_ranks$subtype == test_subtype]
          other_ranks <- gene_ranks$percentile_rank[gene_ranks$subtype != test_subtype]
          
          if(length(target_rank) == 1 && length(other_ranks) >= 2) {
            
            # Is this gene in top percentile for this subtype?
            is_top_gene <- target_rank >= (1 - top_percentile)
            
            # How does this rank compare to other subtypes?
            rank_difference <- target_rank - mean(other_ranks, na.rm = TRUE)
            
            # Simple significance: is this gene much higher ranked than in other subtypes?
            is_significant <- is_top_gene && rank_difference > 0.3  # 30% higher rank
            
            rank_tests <- rbind(rank_tests, data.frame(
              gene = gene,
              subtype = test_subtype,
              target_rank = target_rank,
              other_ranks_mean = mean(other_ranks, na.rm = TRUE),
              rank_difference = rank_difference,
              is_top_gene = is_top_gene,
              is_rank_significant = is_significant,
              target_expr = gene_ranks$expression[gene_ranks$subtype == test_subtype],
              degree = gene_all_data$degree[gene_all_data$subtype == test_subtype],
              hub_score = gene_all_data$composite_hub_score[gene_all_data$subtype == test_subtype],
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    rank_significance[[test_subtype]] <- rank_tests
  }
  
  return(rank_significance)
}

# =============================================================================
# ROBUST COMPREHENSIVE SIGNIFICANCE TEST
# =============================================================================

robust_comprehensive_significance_test <- function(ppi_networks_final, hub_results, 
                                                   fdr_threshold = 0.05, 
                                                   effect_size_threshold = 1.5,
                                                   outlier_threshold = 2,
                                                   top_percentile = 0.1) {
  
  cat("Running ROBUST comprehensive significance analysis...\n")
  
  # Run expression-based tests with error handling
  cat("1. Expression-based tests (with error handling)...\n")
  expr_sig <- expression_based_hub_test_fixed(hub_results, top_n = 50)
  
  # Run outlier-based tests with fixed error handling
  cat("2. Outlier-based tests (fixed)...\n")
  outlier_sig <- tryCatch({
    outlier_based_hub_test_fixed(hub_results, outlier_threshold = outlier_threshold)
  }, error = function(e) {
    cat("  Outlier-based test failed, skipping...\n")
    list()
  })
  
  # Run rank-based tests as backup
  cat("3. Rank-based tests (backup method)...\n")
  rank_sig <- rank_based_hub_test(hub_results, top_percentile = top_percentile)
  
  # Combine results
  significant_hubs <- list()
  
  for(subtype in names(hub_results)) {
    
    # Get results from different tests
    expr_results <- if(subtype %in% names(expr_sig)) expr_sig[[subtype]] else data.frame()
    outlier_results <- if(subtype %in% names(outlier_sig)) outlier_sig[[subtype]] else data.frame()
    rank_results <- if(subtype %in% names(rank_sig)) rank_sig[[subtype]] else data.frame()
    
    # Find significant genes from each method
    expr_significant_genes <- c()
    if(nrow(expr_results) > 0) {
      if("z_fdr" %in% colnames(expr_results)) {
        expr_significant_genes <- c(expr_significant_genes, expr_results$gene[
          !is.na(expr_results$z_fdr) & 
            expr_results$z_fdr < fdr_threshold & 
            !is.na(expr_results$expr_fold_change) &
            abs(log2(expr_results$expr_fold_change)) > log2(effect_size_threshold)
        ])
      }
      if("t_fdr" %in% colnames(expr_results)) {
        t_significant <- expr_results$gene[
          !is.na(expr_results$t_fdr) & 
            expr_results$t_fdr < fdr_threshold & 
            !is.na(expr_results$expr_fold_change) &
            abs(log2(expr_results$expr_fold_change)) > log2(effect_size_threshold)
        ]
        expr_significant_genes <- c(expr_significant_genes, t_significant)
      }
    }
    
    outlier_significant_genes <- c()
    if(nrow(outlier_results) > 0) {
      outlier_significant_genes <- outlier_results$gene[
        !is.na(outlier_results$outlier_fdr) & 
          outlier_results$outlier_fdr < fdr_threshold & 
          !is.na(outlier_results$z_score) &
          abs(outlier_results$z_score) > outlier_threshold
      ]
    }
    
    rank_significant_genes <- c()
    if(nrow(rank_results) > 0) {
      rank_significant_genes <- rank_results$gene[
        !is.na(rank_results$is_rank_significant) & 
          rank_results$is_rank_significant == TRUE
      ]
    }
    
    # Combine all significant genes
    all_significant_genes <- unique(c(expr_significant_genes, outlier_significant_genes, rank_significant_genes))
    all_significant_genes <- all_significant_genes[!is.na(all_significant_genes)]
    
    # Create combined results
    combined_results <- data.frame()
    
    for(gene in all_significant_genes) {
      expr_row <- if(nrow(expr_results) > 0) expr_results[expr_results$gene == gene, ] else data.frame()
      outlier_row <- if(nrow(outlier_results) > 0) outlier_results[outlier_results$gene == gene, ] else data.frame()
      rank_row <- if(nrow(rank_results) > 0) rank_results[rank_results$gene == gene, ] else data.frame()
      
      combined_row <- data.frame(
        gene = gene,
        subtype = subtype,
        expr_significant = gene %in% expr_significant_genes,
        outlier_significant = gene %in% outlier_significant_genes,
        rank_significant = gene %in% rank_significant_genes,
        n_methods_significant = sum(c(gene %in% expr_significant_genes, 
                                      gene %in% outlier_significant_genes,
                                      gene %in% rank_significant_genes)),
        expr_fold_change = ifelse(nrow(expr_row) > 0, expr_row$expr_fold_change[1], 
                                  ifelse(nrow(outlier_row) > 0, outlier_row$expr_fold_change[1], NA)),
        target_expr = ifelse(nrow(expr_row) > 0, expr_row$target_expr[1], 
                             ifelse(nrow(outlier_row) > 0, outlier_row$target_expr[1],
                                    ifelse(nrow(rank_row) > 0, rank_row$target_expr[1], NA))),
        expr_fdr = ifelse(nrow(expr_row) > 0, 
                          ifelse("z_fdr" %in% colnames(expr_row), expr_row$z_fdr[1], 
                                 ifelse("t_fdr" %in% colnames(expr_row), expr_row$t_fdr[1], NA)), NA),
        outlier_fdr = ifelse(nrow(outlier_row) > 0, outlier_row$outlier_fdr[1], NA),
        target_rank = ifelse(nrow(rank_row) > 0, rank_row$target_rank[1], NA),
        stringsAsFactors = FALSE
      )
      
      combined_results <- rbind(combined_results, combined_row)
    }
    
    # Sort by number of significant methods
    if(nrow(combined_results) > 0) {
      combined_results <- combined_results %>%
        arrange(desc(n_methods_significant), expr_fdr, outlier_fdr)
    }
    
    significant_hubs[[subtype]] <- combined_results
  }
  
  return(significant_hubs)
}

# =============================================================================
# RUN ROBUST SIGNIFICANCE ANALYSIS
# =============================================================================

# Run the robust comprehensive analysis
cat("Starting ROBUST statistical significance analysis...\n")
significant_hubs_robust <- robust_comprehensive_significance_test(ppi_networks_final, hub_results)

# Display results
for(subtype in names(significant_hubs_robust)) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat(toupper(subtype), "SUBTYPE - ROBUST SIGNIFICANT HUBS\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  results <- significant_hubs_robust[[subtype]]
  
  if(nrow(results) > 0) {
    # Multiple methods significant
    multi_sig <- results[results$n_methods_significant >= 2, ]
    if(nrow(multi_sig) > 0) {
      cat("\nGenes significant by MULTIPLE methods:\n")
      for(i in 1:min(10, nrow(multi_sig))) {
        methods <- paste(c(
          if(multi_sig$expr_significant[i]) "Expression",
          if(multi_sig$outlier_significant[i]) "Outlier", 
          if(multi_sig$rank_significant[i]) "Rank"
        ), collapse = ", ")
        cat(sprintf("  %s: FC=%.2f, methods: %s\n",
                    multi_sig$gene[i], multi_sig$expr_fold_change[i], methods))
      }
    }
    
    # Single method significant
    single_sig <- results[results$n_methods_significant == 1, ]
    if(nrow(single_sig) > 0) {
      cat(sprintf("\nAdditional %d genes significant by single method (top 5):\n", nrow(single_sig)))
      for(i in 1:min(5, nrow(single_sig))) {
        method <- if(single_sig$expr_significant[i]) "Expression" else 
          if(single_sig$outlier_significant[i]) "Outlier" else "Rank"
        cat(sprintf("  %s: FC=%.2f, method: %s\n",
                    single_sig$gene[i], single_sig$expr_fold_change[i], method))
      }
    }
  } else {
    cat("No significant hubs found for this subtype.\n")
  }
}

# Save robust results
save(significant_hubs_robust, file = "data/statistically_significant_hubs_robust.Rdata")

cat("\n\nROBUST statistical analysis complete! Results saved to data/statistically_significant_hubs_robust.Rdata\n")

# Update the significant_hubs variable for visualization
significant_hubs <- significant_hubs_robust
