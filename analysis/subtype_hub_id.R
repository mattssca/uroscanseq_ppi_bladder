# Find top hubs per subtype
find_top_hubs <- function(hub_analysis, top_n = 20) {
  
  top_hubs <- list()
  
  for(subtype in names(hub_analysis)) {
    df <- hub_analysis[[subtype]]
    
    # Get top hubs by composite score
    top_hubs[[subtype]] <- df %>%
      arrange(hub_score) %>%
      head(top_n) %>%
      select(gene, degree, betweenness, eigenvector, mean_expr, hub_score)
  }
  
  return(top_hubs)
}

# Identify subtype-specific vs shared hubs
compare_hubs_across_subtypes <- function(hub_analysis, top_n = 50) {
  
  # Get top hubs for each subtype
  all_top_hubs <- list()
  for(subtype in names(hub_analysis)) {
    all_top_hubs[[subtype]] <- hub_analysis[[subtype]] %>%
      arrange(hub_score) %>%
      head(top_n) %>%
      pull(gene)
  }
  
  # Find shared vs specific hubs
  all_genes <- unique(unlist(all_top_hubs))
  hub_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_top_hubs))
  rownames(hub_matrix) <- all_genes
  colnames(hub_matrix) <- names(all_top_hubs)
  
  for(subtype in names(all_top_hubs)) {
    hub_matrix[all_top_hubs[[subtype]], subtype] <- 1
  }
  
  # Categorize hubs
  hub_counts <- rowSums(hub_matrix)
  shared_hubs <- names(hub_counts[hub_counts >= 4])  # Present in 4+ subtypes
  subtype_specific <- names(hub_counts[hub_counts == 1])  # Present in only 1 subtype
  
  return(list(
    hub_matrix = hub_matrix,
    shared_hubs = shared_hubs,
    subtype_specific = subtype_specific,
    hub_counts = hub_counts
  ))
}

# Run analyses
top_hubs <- find_top_hubs(hub_analysis, top_n = 20)
hub_comparison <- compare_hubs_across_subtypes(hub_analysis, top_n = 50)

# Visualize hub sharing
pheatmap(hub_comparison$hub_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         main = "Hub Sharing Across Subtypes",
         color = c("white", "red"),  # Simple binary color scheme
         breaks = c(0, 0.5, 1),      # Explicit breaks for binary data
         legend_breaks = c(0, 1),    # Legend breaks
         legend_labels = c("Not Hub", "Hub"),
         filename = "figs/hub_sharing_heatmap.pdf")

library(ggplot2)
library(reshape2)

# Method 1: Use ggplot2 for more control
plot_hub_sharing <- function(hub_matrix) {
  
  # Convert to long format
  hub_df <- melt(hub_matrix, varnames = c("gene", "subtype"), value.name = "is_hub")
  
  # Create the plot
  p <- ggplot(hub_df, aes(x = subtype, y = gene, fill = factor(is_hub))) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_manual(values = c("0" = "white", "1" = "steelblue"),
                      labels = c("Not Hub", "Hub"),
                      name = "Hub Status") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Hub Sharing Across Subtypes",
         x = "Subtype", y = "Gene")
  
  return(p)
}

# Method 2: Create a summary heatmap with hub counts
create_hub_summary <- function(hub_analysis, top_n = 30) {
  
  # Get degree values for top hubs across all subtypes
  all_genes <- unique(unlist(lapply(hub_analysis, function(x) {
    x %>% arrange(hub_score) %>% head(top_n) %>% pull(gene)
  })))
  
  # Create matrix with actual degree values
  degree_matrix <- matrix(0, nrow = length(all_genes), ncol = length(hub_analysis))
  rownames(degree_matrix) <- all_genes
  colnames(degree_matrix) <- names(hub_analysis)
  
  for(subtype in names(hub_analysis)) {
    df <- hub_analysis[[subtype]]
    for(gene in all_genes) {
      gene_row <- df[df$gene == gene, ]
      if(nrow(gene_row) > 0) {
        degree_matrix[gene, subtype] <- gene_row$degree[1]
      }
    }
  }
  
  return(degree_matrix)
}

# Method 3: Focus on most variable hubs
plot_variable_hubs <- function(hub_analysis, top_n = 20) {
  
  # Calculate hub variability across subtypes
  all_hubs <- bind_rows(hub_analysis, .id = "subtype")
  
  hub_variability <- all_hubs %>%
    group_by(gene) %>%
    filter(n() >= 3) %>%  # Present in at least 3 subtypes
    summarise(
      degree_cv = sd(degree, na.rm = TRUE) / mean(degree, na.rm = TRUE),
      max_degree = max(degree, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(degree_cv)) %>%
    head(top_n)
  
  # Create matrix for these variable hubs
  variable_genes <- hub_variability$gene
  variable_matrix <- matrix(0, nrow = length(variable_genes), ncol = length(hub_analysis))
  rownames(variable_matrix) <- variable_genes
  colnames(variable_matrix) <- names(hub_analysis)
  
  for(subtype in names(hub_analysis)) {
    df <- hub_analysis[[subtype]]
    for(gene in variable_genes) {
      gene_row <- df[df$gene == gene, ]
      if(nrow(gene_row) > 0) {
        variable_matrix[gene, subtype] <- gene_row$degree[1]
      }
    }
  }
  
  return(variable_matrix)
}

# Run the improved visualizations
hub_plot <- plot_hub_sharing(hub_comparison$hub_matrix)
ggsave("figs/hub_sharing_ggplot.pdf", hub_plot, width = 8, height = 12)

# Create degree-based heatmap
degree_matrix <- create_hub_summary(hub_analysis, top_n = 30)
pheatmap(degree_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         scale = "row",  # Scale by row to show relative differences
         main = "Hub Degree Across Subtypes",
         filename = "figs/hub_degree_heatmap.pdf")

# Create variability-focused heatmap
variable_matrix <- plot_variable_hubs(hub_analysis, top_n = 20)
pheatmap(variable_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         scale = "row",
         main = "Most Variable Hubs Across Subtypes",
         filename = "figs/variable_hubs_heatmap.pdf")


# Fixed hub analysis function
analyze_subtype_hubs_fixed <- function(ppi_networks, top_n = 20) {
  
  hub_results <- list()
  
  for(subtype in names(ppi_networks)) {
    graph <- ppi_networks[[subtype]]
    
    # Calculate multiple centrality measures
    degrees <- degree(graph)
    betweenness <- betweenness(graph, normalized = TRUE)
    closeness <- closeness(graph, normalized = TRUE)
    eigenvector <- eigen_centrality(graph)$vector
    pagerank <- page_rank(graph)$vector
    
    # Get mean expression and gene symbols
    mean_expr <- V(graph)$mean_expr
    gene_symbols <- V(graph)$gene_symbol
    
    # Combine into dataframe using gene symbols as identifiers
    centrality_df <- data.frame(
      string_id = V(graph)$name,
      gene = gene_symbols,  # Use gene symbols, not STRING IDs
      degree = degrees,
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector,
      pagerank = pagerank,
      mean_expr = mean_expr,
      subtype = subtype
    )
    
    # Remove rows with missing gene symbols
    centrality_df <- centrality_df[!is.na(centrality_df$gene), ]
    
    # Rank genes by different measures
    centrality_df$degree_rank <- rank(-centrality_df$degree)
    centrality_df$betweenness_rank <- rank(-centrality_df$betweenness)
    centrality_df$eigenvector_rank <- rank(-centrality_df$eigenvector)
    
    # Calculate composite hub score
    centrality_df$hub_score <- (centrality_df$degree_rank + 
                                  centrality_df$betweenness_rank + 
                                  centrality_df$eigenvector_rank) / 3
    
    hub_results[[subtype]] <- centrality_df
  }
  
  return(hub_results)
}

# Re-run hub analysis with fixed networks
ppi_networks_fixed <- list(
  uro = ppi_uro_fixed,
  gu = ppi_gu_fixed, 
  basq = ppi_basq_fixed,
  mes = ppi_mes_fixed,
  scne = ppi_scne_fixed
)

hub_analysis_fixed <- analyze_subtype_hubs_fixed(ppi_networks_fixed)

# Check if the fix worked
cat("Hub analysis check - top 5 hubs per subtype:\n")
for(subtype in names(hub_analysis_fixed)) {
  cat(paste("\n", toupper(subtype), "subtype:\n"))
  top_5 <- hub_analysis_fixed[[subtype]] %>%
    arrange(hub_score) %>%
    head(5) %>%
    select(gene, degree, mean_expr, hub_score)
  print(top_5)
}
