# =============================================================================
# EXPORT TOP 20 HUB GENES FOR EACH SUBTYPE
# =============================================================================

export_top_hubs_to_txt <- function(hub_results, top_n = 20, output_dir = "hub_gene_lists") {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Exporting top", top_n, "hub genes for each subtype...\n")
  
  for(subtype in names(hub_results)) {
    
    # Get top hub genes based on composite hub score (lower is better)
    top_hubs <- hub_results[[subtype]] %>%
      arrange(composite_hub_score) %>%
      head(top_n) %>%
      pull(gene)
    
    # Remove any NA or empty gene names
    top_hubs <- top_hubs[!is.na(top_hubs) & top_hubs != ""]
    
    # Create filename
    filename <- file.path(output_dir, paste0(subtype, "_top_", top_n, "_hubs.txt"))
    
    # Write genes to file (one per line)
    writeLines(top_hubs, filename)
    
    cat("Exported", length(top_hubs), "genes for", toupper(subtype), "subtype to", filename, "\n")
  }
  
  cat("\nAll hub gene lists exported to", output_dir, "directory\n")
}

# =============================================================================
# ALTERNATIVE: EXPORT BASED ON DIFFERENT CRITERIA
# =============================================================================

export_hubs_multiple_criteria <- function(hub_results, top_n = 20, output_dir = "hub_gene_lists") {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for(subtype in names(hub_results)) {
    
    subtype_data <- hub_results[[subtype]]
    
    # 1. Top by composite hub score
    top_composite <- subtype_data %>%
      arrange(composite_hub_score) %>%
      head(top_n) %>%
      pull(gene)
    
    # 2. Top by degree
    top_degree <- subtype_data %>%
      arrange(desc(degree)) %>%
      head(top_n) %>%
      pull(gene)
    
    # 3. Top by expression-weighted score
    top_expr_weighted <- subtype_data %>%
      arrange(expr_weighted_score) %>%
      head(top_n) %>%
      pull(gene)
    
    # 4. Top by expression level
    top_expression <- subtype_data %>%
      arrange(desc(mean_expr)) %>%
      head(top_n) %>%
      pull(gene)
    
    # Remove NAs and empty strings
    top_composite <- top_composite[!is.na(top_composite) & top_composite != ""]
    top_degree <- top_degree[!is.na(top_degree) & top_degree != ""]
    top_expr_weighted <- top_expr_weighted[!is.na(top_expr_weighted) & top_expr_weighted != ""]
    top_expression <- top_expression[!is.na(top_expression) & top_expression != ""]
    
    # Write different lists
    writeLines(top_composite, file.path(output_dir, paste0(subtype, "_top_", top_n, "_by_composite_score.txt")))
    writeLines(top_degree, file.path(output_dir, paste0(subtype, "_top_", top_n, "_by_degree.txt")))
    writeLines(top_expr_weighted, file.path(output_dir, paste0(subtype, "_top_", top_n, "_by_expr_weighted.txt")))
    writeLines(top_expression, file.path(output_dir, paste0(subtype, "_top_", top_n, "_by_expression.txt")))
    
    cat("Exported 4 different hub lists for", toupper(subtype), "subtype\n")
  }
  
  cat("\nAll hub gene lists exported to", output_dir, "directory\n")
}

# =============================================================================
# EXPORT SIGNIFICANT HUBS (IF AVAILABLE)
# =============================================================================

export_significant_hubs_to_txt <- function(significant_hubs, top_n = 20, output_dir = "significant_hub_lists") {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Exporting significant hub genes for each subtype...\n")
  
  for(subtype in names(significant_hubs)) {
    
    subtype_data <- significant_hubs[[subtype]]
    
    if(nrow(subtype_data) > 0) {
      # Get top significant genes (sorted by number of methods and fold change)
      top_significant <- subtype_data %>%
        arrange(desc(n_methods_significant), desc(abs(log2(expr_fold_change)))) %>%
        head(top_n) %>%
        pull(gene)
      
      # Remove any NA or empty gene names
      top_significant <- top_significant[!is.na(top_significant) & top_significant != ""]
      
      # Create filename
      filename <- file.path(output_dir, paste0(subtype, "_top_", length(top_significant), "_significant_hubs.txt"))
      
      # Write genes to file (one per line)
      writeLines(top_significant, filename)
      
      cat("Exported", length(top_significant), "significant genes for", toupper(subtype), "subtype\n")
    } else {
      cat("No significant genes found for", toupper(subtype), "subtype\n")
    }
  }
  
  cat("\nAll significant hub gene lists exported to", output_dir, "directory\n")
}

# =============================================================================
# EXPORT SUBTYPE-SPECIFIC HUBS (IF AVAILABLE)
# =============================================================================

export_subtype_specific_hubs <- function(specificity_results, output_dir = "subtype_specific_lists") {
  
  if(exists("specificity_results") && !is.null(specificity_results)) {
    
    # Create output directory if it doesn't exist
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    specificity_analysis <- specificity_results$specificity_analysis
    
    for(subtype in names(hub_results)) {
      # Get genes specific to this subtype
      specific_genes <- names(specificity_analysis$specific_details)[
        sapply(specificity_analysis$specific_details, function(x) x$subtype == subtype)
      ]
      
      if(length(specific_genes) > 0) {
        # Sort by evidence level and fold change
        specific_info <- specificity_analysis$specific_details[specific_genes]
        specific_sorted <- specific_info[order(
          -sapply(specific_info, function(x) x$n_methods),
          -sapply(specific_info, function(x) abs(log2(x$expr_fold_change)))
        )]
        
        specific_genes_sorted <- names(specific_sorted)
        
        # Create filename
        filename <- file.path(output_dir, paste0(subtype, "_specific_hubs.txt"))
        
        # Write genes to file
        writeLines(specific_genes_sorted, filename)
        
        cat("Exported", length(specific_genes_sorted), "subtype-specific genes for", toupper(subtype), "\n")
      } else {
        cat("No subtype-specific genes found for", toupper(subtype), "\n")
      }
    }
  }
}

# =============================================================================
# RUN EXPORTS
# =============================================================================

# 1. Export top 20 hubs by composite score (main method)
export_top_hubs_to_txt(hub_results, top_n = 20, output_dir = "hub_gene_lists")

# 2. Export top hubs by multiple criteria (optional)
export_hubs_multiple_criteria(hub_results, top_n = 20, output_dir = "hub_gene_lists_detailed")

# 3. Export significant hubs (if significant_hubs exists)
if(exists("significant_hubs")) {
  export_significant_hubs_to_txt(significant_hubs, top_n = 20, output_dir = "significant_hub_lists")
}

# 4. Export subtype-specific hubs (if specificity analysis was run)
if(exists("specificity_results")) {
  export_subtype_specific_hubs(specificity_results, output_dir = "subtype_specific_lists")
}

# =============================================================================
# PRINT SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("EXPORT SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Files created:\n")
cat("1. hub_gene_lists/ - Top 20 hubs by composite score\n")
cat("2. hub_gene_lists_detailed/ - Top 20 hubs by various criteria\n")
if(exists("significant_hubs")) {
  cat("3. significant_hub_lists/ - Statistically significant hubs\n")
}
if(exists("specificity_results")) {
  cat("4. subtype_specific_lists/ - Subtype-specific hubs only\n")
}
cat("\nEach file contains one gene symbol per line.\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# =============================================================================
# PREVIEW FILES
# =============================================================================

cat("\nPREVIEW OF EXPORTED FILES:\n")
cat(paste(rep("-", 40), collapse=""), "\n")

for(subtype in names(hub_results)) {
  filename <- file.path("hub_gene_lists", paste0(subtype, "_top_20_hubs.txt"))
  
  if(file.exists(filename)) {
    genes <- readLines(filename)
    cat("\n", toupper(subtype), "subtype (", length(genes), "genes):\n")
    cat("  ", paste(head(genes, 5), collapse = ", "))
    if(length(genes) > 5) cat(", ...")
    cat("\n")
  }
}
