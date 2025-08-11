library(igraph)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(scales)

# =============================================================================
# PART 1: NETWORK OVERVIEW AND VALIDATION
# =============================================================================

# Basic network properties
network_overview <- function(ppi_networks_final) {
  
  overview_stats <- data.frame()
  
  for(subtype in names(ppi_networks_final)) {
    graph <- ppi_networks_final[[subtype]]
    
    stats <- data.frame(
      subtype = subtype,
      nodes = vcount(graph),
      edges = ecount(graph),
      density = edge_density(graph),
      avg_degree = mean(degree(graph)),
      max_degree = max(degree(graph)),
      clustering_coeff = transitivity(graph),
      avg_expression = mean(V(graph)$mean_expr, na.rm = TRUE),
      sd_expression = sd(V(graph)$mean_expr, na.rm = TRUE),
      min_expression = min(V(graph)$mean_expr, na.rm = TRUE),
      max_expression = max(V(graph)$mean_expr, na.rm = TRUE)
    )
    
    overview_stats <- rbind(overview_stats, stats)
  }
  
  return(overview_stats)
}

# Run overview
network_stats <- network_overview(ppi_networks_final)
print("Network Overview:")
print(network_stats)

# =============================================================================
# PART 2: COMPREHENSIVE HUB ANALYSIS
# =============================================================================

# Enhanced hub analysis with multiple centrality measures
comprehensive_hub_analysis <- function(ppi_networks_final) {
  
  all_hub_results <- list()
  
  for(subtype in names(ppi_networks_final)) {
    graph <- ppi_networks_final[[subtype]]
    
    # Calculate centrality measures
    degrees <- degree(graph)
    betweenness <- betweenness(graph, normalized = TRUE)
    closeness <- closeness(graph, normalized = TRUE)
    eigenvector <- eigen_centrality(graph)$vector
    pagerank <- page_rank(graph)$vector
    
    # Get node attributes
    gene_symbols <- V(graph)$gene_symbol
    mean_expr <- V(graph)$mean_expr
    string_ids <- V(graph)$name
    
    # Create comprehensive dataframe
    hub_df <- data.frame(
      string_id = string_ids,
      gene = gene_symbols,
      degree = degrees,
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector,
      pagerank = pagerank,
      mean_expr = mean_expr,
      subtype = subtype,
      stringsAsFactors = FALSE
    )
    
    # Remove rows with missing gene symbols
    hub_df <- hub_df[!is.na(hub_df$gene) & hub_df$gene != "", ]
    
    # Calculate ranks for each centrality measure
    hub_df$degree_rank <- rank(-hub_df$degree, ties.method = "min")
    hub_df$betweenness_rank <- rank(-hub_df$betweenness, ties.method = "min")
    hub_df$eigenvector_rank <- rank(-hub_df$eigenvector, ties.method = "min")
    hub_df$pagerank_rank <- rank(-hub_df$pagerank, ties.method = "min")
    
    # Composite hub score (lower is better - top hub)
    hub_df$composite_hub_score <- (hub_df$degree_rank + 
                                     hub_df$betweenness_rank + 
                                     hub_df$eigenvector_rank + 
                                     hub_df$pagerank_rank) / 4
    
    # Expression-weighted hub score
    expr_percentile <- rank(hub_df$mean_expr) / nrow(hub_df)
    hub_df$expr_weighted_score <- hub_df$composite_hub_score * (2 - expr_percentile)
    
    all_hub_results[[subtype]] <- hub_df
  }
  
  return(all_hub_results)
}

# Run comprehensive hub analysis
hub_results <- comprehensive_hub_analysis(ppi_networks_final)

# Save results
save(hub_results, file = "data/comprehensive_hub_results.Rdata")

# =============================================================================
# PART 3: SUBTYPE-SPECIFIC HUB IDENTIFICATION
# =============================================================================

# Find top hubs for each subtype
get_top_hubs <- function(hub_results, top_n = 20, score_type = "composite_hub_score") {
  
  top_hubs_list <- list()
  
  for(subtype in names(hub_results)) {
    df <- hub_results[[subtype]]
    
    top_hubs <- df %>%
      arrange(get(score_type)) %>%
      head(top_n) %>%
      select(gene, degree, betweenness, eigenvector, mean_expr, 
             composite_hub_score, expr_weighted_score)
    
    top_hubs_list[[subtype]] <- top_hubs
  }
  
  return(top_hubs_list)
}

# Get top hubs by different criteria
top_hubs_composite <- get_top_hubs(hub_results, top_n = 25, "composite_hub_score")
top_hubs_expression <- get_top_hubs(hub_results, top_n = 25, "expr_weighted_score")

# Print top hubs for each subtype
cat("=== TOP 10 HUBS BY COMPOSITE SCORE ===\n")
for(subtype in names(top_hubs_composite)) {
  cat(paste("\n", toupper(subtype), "Subtype Top Hubs:\n"))
  print(head(top_hubs_composite[[subtype]], 10))
}

# =============================================================================
# PART 4: CROSS-SUBTYPE HUB COMPARISON
# =============================================================================

# Compare hub rankings across subtypes
compare_hub_rankings <- function(hub_results, top_n = 50) {
  
  # Get top genes from each subtype
  all_top_genes <- unique(unlist(lapply(hub_results, function(x) {
    x %>% arrange(composite_hub_score) %>% head(top_n) %>% pull(gene)
  })))
  
  # Create ranking matrix
  ranking_matrix <- matrix(NA, nrow = length(all_top_genes), ncol = length(hub_results))
  rownames(ranking_matrix) <- all_top_genes
  colnames(ranking_matrix) <- names(hub_results)
  
  # Fill in rankings
  for(subtype in names(hub_results)) {
    df <- hub_results[[subtype]]
    for(gene in all_top_genes) {
      gene_data <- df[df$gene == gene, ]
      if(nrow(gene_data) > 0) {
        ranking_matrix[gene, subtype] <- gene_data$composite_hub_score[1]
      }
    }
  }
  
  return(ranking_matrix)
}

# Create expression matrix for top hubs
create_expression_matrix <- function(hub_results, top_n = 50) {
  
  # Get top genes from each subtype
  all_top_genes <- unique(unlist(lapply(hub_results, function(x) {
    x %>% arrange(composite_hub_score) %>% head(top_n) %>% pull(gene)
  })))
  
  # Create expression matrix
  expr_matrix <- matrix(NA, nrow = length(all_top_genes), ncol = length(hub_results))
  rownames(expr_matrix) <- all_top_genes
  colnames(expr_matrix) <- names(hub_results)
  
  # Fill in expression values
  for(subtype in names(hub_results)) {
    df <- hub_results[[subtype]]
    for(gene in all_top_genes) {
      gene_data <- df[df$gene == gene, ]
      if(nrow(gene_data) > 0) {
        expr_matrix[gene, subtype] <- gene_data$mean_expr[1]
      }
    }
  }
  
  return(expr_matrix)
}

# Create matrices
hub_ranking_matrix <- compare_hub_rankings(hub_results, top_n = 40)
hub_expression_matrix <- create_expression_matrix(hub_results, top_n = 40)

# =============================================================================
# PART 5: SUBTYPE-SPECIFIC VS SHARED HUBS
# =============================================================================

# Identify subtype-specific vs shared hub patterns
analyze_hub_specificity <- function(hub_results, top_n = 30) {
  
  # Get top hubs for each subtype
  subtype_top_hubs <- lapply(hub_results, function(x) {
    x %>% arrange(composite_hub_score) %>% head(top_n) %>% pull(gene)
  })
  
  # Create presence/absence matrix
  all_genes <- unique(unlist(subtype_top_hubs))
  presence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(subtype_top_hubs))
  rownames(presence_matrix) <- all_genes
  colnames(presence_matrix) <- names(subtype_top_hubs)
  
  for(subtype in names(subtype_top_hubs)) {
    presence_matrix[subtype_top_hubs[[subtype]], subtype] <- 1
  }
  
  # Analyze patterns
  hub_counts <- rowSums(presence_matrix)
  
  shared_hubs <- names(hub_counts[hub_counts >= 4])  # Present in 4+ subtypes
  mostly_shared <- names(hub_counts[hub_counts == 3])  # Present in 3 subtypes
  subtype_specific <- names(hub_counts[hub_counts == 1])  # Present in 1 subtype
  
  # Get subtype-specific details
  specific_details <- list()
  for(gene in subtype_specific) {
    subtype_with_gene <- colnames(presence_matrix)[which(presence_matrix[gene, ] == 1)]
    specific_details[[gene]] <- subtype_with_gene
  }
  
  return(list(
    presence_matrix = presence_matrix,
    hub_counts = hub_counts,
    shared_hubs = shared_hubs,
    mostly_shared = mostly_shared,
    subtype_specific = subtype_specific,
    specific_details = specific_details
  ))
}

# Run specificity analysis
hub_specificity <- analyze_hub_specificity(hub_results, top_n = 30)

# Print results
cat("\n=== HUB SPECIFICITY ANALYSIS ===\n")
cat("Shared hubs (4+ subtypes):", length(hub_specificity$shared_hubs), "\n")
cat("Mostly shared hubs (3 subtypes):", length(hub_specificity$mostly_shared), "\n")
cat("Subtype-specific hubs (1 subtype):", length(hub_specificity$subtype_specific), "\n")

cat("\nSubtype-specific hub details:\n")
for(gene in head(hub_specificity$subtype_specific, 10)) {
  cat(paste(gene, "->", hub_specificity$specific_details[[gene]], "\n"))
}

# =============================================================================
# PART 6: EXPRESSION-BASED HUB ANALYSIS
# =============================================================================

# Find hubs with highest expression variability across subtypes
find_variable_expression_hubs <- function(hub_results, min_subtypes = 3) {
  
  # Combine all hub data
  all_hubs <- bind_rows(hub_results, .id = "subtype")
  
  # Calculate expression variability for each gene
  expr_variability <- all_hubs %>%
    group_by(gene) %>%
    filter(n() >= min_subtypes) %>%
    summarise(
      mean_expression = mean(mean_expr, na.rm = TRUE),
      sd_expression = sd(mean_expr, na.rm = TRUE),
      cv_expression = sd_expression / mean_expression,
      max_expression = max(mean_expr, na.rm = TRUE),
      min_expression = min(mean_expr, na.rm = TRUE),
      fold_change = max_expression / (min_expression + 0.001),
      subtypes_present = n(),
      avg_degree = mean(degree, na.rm = TRUE),
      avg_hub_score = mean(composite_hub_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(cv_expression))
  
  return(expr_variability)
}

# Identify expression-variable hubs
variable_hubs <- find_variable_expression_hubs(hub_results, min_subtypes = 4)

cat("\n=== TOP 15 EXPRESSION-VARIABLE HUBS ===\n")
print(head(variable_hubs, 15))

# =============================================================================
# PART 7: VISUALIZATIONS
# =============================================================================

# Create comprehensive visualizations
create_hub_visualizations <- function(hub_expression_matrix, hub_specificity, network_stats) {
  
  # 1. Expression heatmap of top hubs
  # Remove rows with all NA values
  complete_rows <- complete.cases(hub_expression_matrix)
  clean_matrix <- hub_expression_matrix[complete_rows, ]
  
  if(nrow(clean_matrix) > 0) {
    pheatmap(clean_matrix, 
             cluster_rows = TRUE, 
             cluster_cols = FALSE,
             scale = "row",
             main = "Hub Gene Expression Across Subtypes",
             filename = "figs/hub_expression_heatmap.pdf",
             width = 8, height = 12)
  }
  
  # 2. Hub sharing pattern
  pheatmap(hub_specificity$presence_matrix, 
           cluster_rows = TRUE, 
           cluster_cols = FALSE,
           color = c("white", "darkblue"),
           breaks = c(0, 0.5, 1),
           legend_breaks = c(0, 1),
           legend_labels = c("Not Top Hub", "Top Hub"),
           main = "Hub Status Across Subtypes",
           filename = "figs/hub_sharing_pattern.pdf",
           width = 8, height = 12)
  
  # 3. Network properties comparison
  stats_plot <- network_stats %>%
    select(subtype, avg_expression, sd_expression, avg_degree, clustering_coeff) %>%
    pivot_longer(-subtype, names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = subtype, y = value, fill = subtype)) +
    geom_bar(stat = "identity") +
    facet_wrap(~metric, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Network Properties Across Subtypes")
  
  ggsave("figs/network_properties_comparison.pdf", stats_plot, width = 12, height = 8)
  
  # 4. Expression distribution by subtype
  expr_dist_plot <- bind_rows(hub_results, .id = "subtype") %>%
    ggplot(aes(x = mean_expr, fill = subtype)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~subtype, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expression Distribution by Subtype",
         x = "Mean Expression", y = "Density")
  
  ggsave("figs/expression_distribution_by_subtype.pdf", expr_dist_plot, width = 12, height = 8)
}

# Generate all visualizations
create_hub_visualizations(hub_expression_matrix, hub_specificity, network_stats)

# =============================================================================
# PART 8: SUMMARY REPORT
# =============================================================================

# Fixed summary report function
generate_summary_report <- function(network_stats, hub_specificity, variable_hubs, hub_results) {
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
  cat("COMPREHENSIVE SUBTYPE-SPECIFIC HUB ANALYSIS REPORT\n")
  cat(paste(rep("=", 60), collapse=""), "\n", sep="")
  
  cat("\n1. NETWORK OVERVIEW:\n")
  cat("   - All networks have", network_stats$nodes[1], "nodes and", network_stats$edges[1], "edges\n")
  cat("   - Networks share topology but differ in expression patterns\n")
  
  cat("\n2. EXPRESSION PATTERNS:\n")
  for(i in 1:nrow(network_stats)) {
    cat(sprintf("   - %s: avg expr = %.2f (SD = %.2f)\n", 
                toupper(network_stats$subtype[i]), 
                network_stats$avg_expression[i], 
                network_stats$sd_expression[i]))
  }
  
  cat("\n3. HUB SPECIFICITY:\n")
  cat("   - Shared hubs (4+ subtypes):", length(hub_specificity$shared_hubs), "\n")
  cat("   - Subtype-specific hubs:", length(hub_specificity$subtype_specific), "\n")
  
  cat("\n4. TOP SHARED HUBS:\n")
  if(length(hub_specificity$shared_hubs) > 0) {
    for(gene in head(hub_specificity$shared_hubs, 5)) {
      cat("   -", gene, "\n")
    }
  } else {
    cat("   - No shared hubs found\n")
  }
  
  cat("\n5. TOP EXPRESSION-VARIABLE HUBS:\n")
  for(i in 1:min(5, nrow(variable_hubs))) {
    cat(sprintf("   - %s (CV = %.2f, FC = %.2f)\n", 
                variable_hubs$gene[i], 
                variable_hubs$cv_expression[i], 
                variable_hubs$fold_change[i]))
  }
  
  cat("\n6. SUBTYPE-SPECIFIC HUB EXAMPLES:\n")
  if(length(hub_specificity$subtype_specific) > 0) {
    spec_examples <- head(hub_specificity$subtype_specific, 5)
    for(gene in spec_examples) {
      subtype <- hub_specificity$specific_details[[gene]]
      cat(sprintf("   - %s (specific to %s)\n", gene, subtype))
    }
  } else {
    cat("   - No subtype-specific hubs found\n")
  }
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
}

# Run the corrected function
generate_summary_report(network_stats, hub_specificity, variable_hubs, hub_results)

# =============================================================================
# SAVE ALL RESULTS
# =============================================================================

# Save comprehensive results
analysis_results <- list(
  network_stats = network_stats,
  hub_results = hub_results,
  hub_specificity = hub_specificity,
  variable_hubs = variable_hubs,
  top_hubs_composite = top_hubs_composite,
  hub_expression_matrix = hub_expression_matrix
)

save(analysis_results, file = "data/complete_subtype_analysis_results.Rdata")

cat("\nAnalysis complete! Results saved to data/complete_subtype_analysis_results.Rdata\n")
cat("Visualizations saved to figs/ directory\n")
