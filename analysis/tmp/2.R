library(igraph)
library(dplyr)
library(broom)

outlier_based_hub_test <- function(hub_results, outlier_threshold = 2) {
  
  # Combine all hub data
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  outlier_significance <- list()
  
  for(test_subtype in names(hub_results)) {
    
    # Get genes present in multiple subtypes
    multi_subtype_genes <- all_hubs %>%
      group_by(gene) %>%
      filter(n() >= 3) %>%  # Present in at least 3 subtypes
      pull(gene) %>%
      unique()
    
    outlier_tests <- data.frame()
    
    for(gene in multi_subtype_genes) {
      gene_data <- all_hubs[all_hubs$gene == gene, ]
      
      target_expr <- gene_data$mean_expr[gene_data$subtype == test_subtype]
      other_expr <- gene_data$mean_expr[gene_data$subtype != test_subtype]
      
      if(length(target_expr) > 0 && length(other_expr) >= 2) {
        
        # Calculate Z-score relative to other subtypes
        other_mean <- mean(other_expr)
        other_sd <- sd(other_expr)
        
        if(other_sd > 0) {
          z_score <- (target_expr - other_mean) / other_sd
          
          # Is this an outlier?
          is_outlier <- abs(z_score) > outlier_threshold
          
          # Calculate approximate p-value
          p_value <- 2 * (1 - pnorm(abs(z_score)))
          
          outlier_tests <- rbind(outlier_tests, data.frame(
            gene = gene,
            subtype = test_subtype,
            target_expr = target_expr,
            other_expr_mean = other_mean,
            other_expr_sd = other_sd,
            z_score = z_score,
            is_outlier = is_outlier,
            outlier_pvalue = p_value,
            expr_fold_change = target_expr / other_mean,
            degree = gene_data$degree[gene_data$subtype == test_subtype],
            hub_score = gene_data$composite_hub_score[gene_data$subtype == test_subtype],
            stringsAsFactors = FALSE
          ))
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
# METHOD 1: PERMUTATION TEST FOR HUB SIGNIFICANCE
# =============================================================================

# Test if a gene's hub status in a subtype is significantly different from random
permutation_hub_test <- function(ppi_networks_final, hub_results, n_permutations = 1000) {
  
  hub_significance <- list()
  
  for(subtype in names(ppi_networks_final)) {
    graph <- ppi_networks_final[[subtype]]
    observed_degrees <- degree(graph)
    
    # Get gene symbols for mapping
    gene_symbols <- V(graph)$gene_symbol
    names(observed_degrees) <- gene_symbols
    
    # Generate null distribution by permuting expression values
    null_distributions <- matrix(0, nrow = n_permutations, ncol = length(observed_degrees))
    colnames(null_distributions) <- gene_symbols
    
    cat("Running permutations for", subtype, "subtype...\n")
    
    for(i in 1:n_permutations) {
      # Permute expression values
      shuffled_expr <- sample(V(graph)$mean_expr)
      V(graph)$mean_expr <- shuffled_expr
      
      # Calculate degrees (topology stays same, but we track this for consistency)
      null_distributions[i, ] <- degree(graph)
    }
    
    # Calculate p-values for each gene
    gene_pvalues <- sapply(gene_symbols, function(gene) {
      if(is.na(gene) || gene == "") return(NA)
      
      observed_degree <- observed_degrees[gene]
      null_degrees <- null_distributions[, gene]
      
      # Two-tailed test: is this gene's degree significantly different?
      p_value <- (sum(null_degrees >= observed_degree) + 1) / (n_permutations + 1)
      return(min(p_value, 1 - p_value) * 2)  # Two-tailed
    })
    
    # Combine with hub results
    hub_data <- hub_results[[subtype]]
    hub_data$permutation_pvalue <- gene_pvalues[hub_data$gene]
    hub_data$permutation_fdr <- p.adjust(hub_data$permutation_pvalue, method = "fdr")
    
    hub_significance[[subtype]] <- hub_data
  }
  
  return(hub_significance)
}

# =============================================================================
# METHOD 2: EXPRESSION-BASED STATISTICAL TEST
# =============================================================================

expression_based_hub_test_fixed <- function(hub_results, top_n = 50) {
  
  # Combine all hub data
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  expression_significance <- list()
  
  for(test_subtype in names(hub_results)) {
    
    # Get top hubs for this subtype
    subtype_top_hubs <- hub_results[[test_subtype]] %>%
      arrange(composite_hub_score) %>%
      head(top_n) %>%
      pull(gene)
    
    # For each top hub, test if its expression is significantly different in this subtype
    hub_expr_tests <- data.frame()
    
    for(gene in subtype_top_hubs) {
      gene_data <- all_hubs[all_hubs$gene == gene, ]
      
      # Need at least 4 subtypes for meaningful statistical test
      if(nrow(gene_data) >= 4) {
        target_expr <- gene_data$mean_expr[gene_data$subtype == test_subtype]
        other_expr <- gene_data$mean_expr[gene_data$subtype != test_subtype]
        
        # Check we have enough observations
        if(length(target_expr) >= 1 && length(other_expr) >= 3) {
          
          # Calculate basic statistics
          target_mean <- mean(target_expr)
          other_mean <- mean(other_expr)
          other_sd <- sd(other_expr)
          
          # Calculate fold change
          expr_fold_change <- target_mean / other_mean
          
          # Method 1: Z-score test (when we have only 1 target observation)
          if(length(target_expr) == 1 && length(other_expr) >= 3) {
            z_score <- (target_mean - other_mean) / (other_sd / sqrt(length(other_expr)))
            z_pvalue <- 2 * (1 - pnorm(abs(z_score)))  # Two-tailed
            
            hub_expr_tests <- rbind(hub_expr_tests, data.frame(
              gene = gene,
              subtype = test_subtype,
              target_expr = target_mean,
              other_expr_mean = other_mean,
              other_expr_sd = other_sd,
              expr_fold_change = expr_fold_change,
              z_score = z_score,
              z_pvalue = z_pvalue,
              test_type = "z_test",
              degree = gene_data$degree[gene_data$subtype == test_subtype],
              hub_score = gene_data$composite_hub_score[gene_data$subtype == test_subtype],
              n_other_subtypes = length(other_expr),
              stringsAsFactors = FALSE
            ))
          }
          
          # Method 2: T-test (when we have multiple observations)
          else if(length(target_expr) > 1 && length(other_expr) > 1) {
            tryCatch({
              t_test <- t.test(target_expr, other_expr)
              wilcox_test <- wilcox.test(target_expr, other_expr)
              
              hub_expr_tests <- rbind(hub_expr_tests, data.frame(
                gene = gene,
                subtype = test_subtype,
                target_expr = target_mean,
                other_expr_mean = other_mean,
                other_expr_sd = other_sd,
                expr_fold_change = expr_fold_change,
                t_pvalue = t_test$p.value,
                wilcox_pvalue = wilcox_test$p.value,
                test_type = "t_test",
                degree = gene_data$degree[gene_data$subtype == test_subtype][1],
                hub_score = gene_data$composite_hub_score[gene_data$subtype == test_subtype][1],
                n_other_subtypes = length(other_expr),
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              cat("Warning: Failed t-test for gene", gene, "in", test_subtype, "\n")
            })
          }
        }
      }
    }
    
    # Adjust for multiple testing
    if(nrow(hub_expr_tests) > 0) {
      # Use appropriate p-value column based on test type
      if("z_pvalue" %in% colnames(hub_expr_tests)) {
        hub_expr_tests$z_fdr <- p.adjust(hub_expr_tests$z_pvalue, method = "fdr")
      }
      if("t_pvalue" %in% colnames(hub_expr_tests)) {
        hub_expr_tests$t_fdr <- p.adjust(hub_expr_tests$t_pvalue, method = "fdr")
      }
      if("wilcox_pvalue" %in% colnames(hub_expr_tests)) {
        hub_expr_tests$wilcox_fdr <- p.adjust(hub_expr_tests$wilcox_pvalue, method = "fdr")
      }
    }
    
    expression_significance[[test_subtype]] <- hub_expr_tests
  }
  
  return(expression_significance)
}

# =============================================================================
# METHOD 3: DEGREE-BASED SIGNIFICANCE TEST
# =============================================================================

# Test if hub genes have significantly higher degree in a specific subtype
degree_significance_test <- function(hub_results, min_degree_threshold = 10) {
  
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  degree_significance <- list()
  
  for(test_subtype in names(hub_results)) {
    
    # Get genes present in multiple subtypes
    multi_subtype_genes <- all_hubs %>%
      group_by(gene) %>%
      filter(n() >= 3) %>%  # Present in at least 3 subtypes
      filter(max(degree) >= min_degree_threshold) %>%  # Has reasonable connectivity
      pull(gene) %>%
      unique()
    
    gene_degree_tests <- data.frame()
    
    for(gene in multi_subtype_genes) {
      gene_data <- all_hubs[all_hubs$gene == gene, ]
      
      target_degree <- gene_data$degree[gene_data$subtype == test_subtype]
      other_degrees <- gene_data$degree[gene_data$subtype != test_subtype]
      
      if(length(target_degree) > 0 && length(other_degrees) > 0) {
        # Test if degree in target subtype is significantly different
        degree_fold_change <- target_degree / mean(other_degrees)
        
        # Z-score approach
        z_score <- (target_degree - mean(other_degrees)) / sd(other_degrees)
        z_pvalue <- 2 * (1 - pnorm(abs(z_score)))  # Two-tailed
        
        gene_degree_tests <- rbind(gene_degree_tests, data.frame(
          gene = gene,
          subtype = test_subtype,
          target_degree = target_degree,
          other_degree_mean = mean(other_degrees),
          degree_fold_change = degree_fold_change,
          z_score = z_score,
          z_pvalue = z_pvalue,
          target_expr = gene_data$mean_expr[gene_data$subtype == test_subtype],
          hub_score = gene_data$composite_hub_score[gene_data$subtype == test_subtype]
        ))
      }
    }
    
    # Adjust for multiple testing
    if(nrow(gene_degree_tests) > 0) {
      gene_degree_tests$z_fdr <- p.adjust(gene_degree_tests$z_pvalue, method = "fdr")
    }
    
    degree_significance[[test_subtype]] <- gene_degree_tests
  }
  
  return(degree_significance)
}

# =============================================================================
# METHOD 4: COMPREHENSIVE SIGNIFICANCE ANALYSIS
# =============================================================================

comprehensive_significance_test_fixed <- function(ppi_networks_final, hub_results, 
                                                  fdr_threshold = 0.05, 
                                                  effect_size_threshold = 1.5,
                                                  outlier_threshold = 2) {
  
  cat("Running FIXED comprehensive significance analysis...\n")
  
  # Run expression-based tests with error handling
  cat("1. Expression-based tests (with error handling)...\n")
  expr_sig <- expression_based_hub_test_fixed(hub_results, top_n = 50)
  
  # Run outlier-based tests
  cat("2. Outlier-based tests...\n")
  outlier_sig <- outlier_based_hub_test(hub_results, outlier_threshold = outlier_threshold)
  
  # Combine results
  significant_hubs <- list()
  
  for(subtype in names(hub_results)) {
    
    # Get expression significance results
    expr_results <- expr_sig[[subtype]]
    outlier_results <- outlier_sig[[subtype]]
    
    # Find significant genes from expression tests
    expr_significant_genes <- c()
    if(nrow(expr_results) > 0) {
      # Use appropriate FDR column based on test type
      if("z_fdr" %in% colnames(expr_results)) {
        expr_significant_genes <- expr_results$gene[
          expr_results$z_fdr < fdr_threshold & 
            abs(log2(expr_results$expr_fold_change)) > log2(effect_size_threshold)
        ]
      }
      if("t_fdr" %in% colnames(expr_results)) {
        t_significant <- expr_results$gene[
          expr_results$t_fdr < fdr_threshold & 
            abs(log2(expr_results$expr_fold_change)) > log2(effect_size_threshold)
        ]
        expr_significant_genes <- c(expr_significant_genes, t_significant)
      }
    }
    
    # Find significant genes from outlier tests
    outlier_significant_genes <- c()
    if(nrow(outlier_results) > 0) {
      outlier_significant_genes <- outlier_results$gene[
        outlier_results$outlier_fdr < fdr_threshold & 
          abs(outlier_results$z_score) > outlier_threshold
      ]
    }
    
    # Combine all significant genes
    all_significant_genes <- unique(c(expr_significant_genes, outlier_significant_genes))
    
    # Create combined results
    combined_results <- data.frame()
    
    for(gene in all_significant_genes) {
      expr_row <- if(nrow(expr_results) > 0) expr_results[expr_results$gene == gene, ] else data.frame()
      outlier_row <- if(nrow(outlier_results) > 0) outlier_results[outlier_results$gene == gene, ] else data.frame()
      
      combined_row <- data.frame(
        gene = gene,
        subtype = subtype,
        expr_significant = gene %in% expr_significant_genes,
        outlier_significant = gene %in% outlier_significant_genes,
        both_significant = gene %in% expr_significant_genes & gene %in% outlier_significant_genes,
        expr_fold_change = ifelse(nrow(expr_row) > 0, expr_row$expr_fold_change[1], 
                                  ifelse(nrow(outlier_row) > 0, outlier_row$expr_fold_change[1], NA)),
        target_expr = ifelse(nrow(expr_row) > 0, expr_row$target_expr[1], 
                             ifelse(nrow(outlier_row) > 0, outlier_row$target_expr[1], NA)),
        z_score = ifelse(nrow(outlier_row) > 0, outlier_row$z_score[1], 
                         ifelse(nrow(expr_row) > 0 && "z_score" %in% colnames(expr_row), expr_row$z_score[1], NA)),
        outlier_fdr = ifelse(nrow(outlier_row) > 0, outlier_row$outlier_fdr[1], NA),
        expr_fdr = ifelse(nrow(expr_row) > 0, 
                          ifelse("z_fdr" %in% colnames(expr_row), expr_row$z_fdr[1], 
                                 ifelse("t_fdr" %in% colnames(expr_row), expr_row$t_fdr[1], NA)), NA),
        target_degree = ifelse(nrow(expr_row) > 0, expr_row$degree[1], 
                               ifelse(nrow(outlier_row) > 0, outlier_row$degree[1], NA)),
        stringsAsFactors = FALSE
      )
      
      combined_results <- rbind(combined_results, combined_row)
    }
    
    # Sort by significance
    if(nrow(combined_results) > 0) {
      combined_results <- combined_results %>%
        arrange(desc(both_significant), expr_fdr, outlier_fdr)
    }
    
    significant_hubs[[subtype]] <- combined_results
  }
  
  return(significant_hubs)
}
# =============================================================================
# RUN SIGNIFICANCE ANALYSIS
# =============================================================================

# Run the fixed comprehensive analysis
cat("Starting FIXED statistical significance analysis...\n")
significant_hubs_fixed <- comprehensive_significance_test_fixed(ppi_networks_final, hub_results)

# Display results
for(subtype in names(significant_hubs_fixed)) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat(toupper(subtype), "SUBTYPE - SIGNIFICANT HUBS (FIXED)\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  results <- significant_hubs_fixed[[subtype]]
  
  if(nrow(results) > 0) {
    # Significant in both tests
    both_sig <- results[results$both_significant == TRUE, ]
    if(nrow(both_sig) > 0) {
      cat("\nGenes significant in BOTH expression and outlier tests:\n")
      for(i in 1:min(10, nrow(both_sig))) {
        cat(sprintf("  %s: expr FC=%.2f, Z-score=%.2f, FDRs: expr=%.3f, outlier=%.3f\n",
                    both_sig$gene[i], both_sig$expr_fold_change[i], both_sig$z_score[i],
                    both_sig$expr_fdr[i], both_sig$outlier_fdr[i]))
      }
    }
    
    # Expression significant only
    expr_only <- results[results$expr_significant == TRUE & results$outlier_significant == FALSE, ]
    if(nrow(expr_only) > 0) {
      cat("\nGenes significant in expression only:\n")
      for(i in 1:min(5, nrow(expr_only))) {
        cat(sprintf("  %s: expr FC=%.2f (FDR=%.3f)\n",
                    expr_only$gene[i], expr_only$expr_fold_change[i], expr_only$expr_fdr[i]))
      }
    }
    
    # Outlier significant only
    outlier_only <- results[results$outlier_significant == TRUE & results$expr_significant == FALSE, ]
    if(nrow(outlier_only) > 0) {
      cat("\nGenes significant as outliers only:\n")
      for(i in 1:min(5, nrow(outlier_only))) {
        cat(sprintf("  %s: Z-score=%.2f (FDR=%.3f)\n",
                    outlier_only$gene[i], outlier_only$z_score[i], outlier_only$outlier_fdr[i]))
      }
    }
  } else {
    cat("No significant hubs found for this subtype.\n")
  }
}

# Save fixed results
save(significant_hubs_fixed, file = "data/statistically_significant_hubs_fixed.Rdata")

cat("\n\nFIXED statistical analysis complete! Results saved to data/statistically_significant_hubs_fixed.Rdata\n")

# Update the significant_hubs variable for visualization
significant_hubs <- significant_hubs_fixed