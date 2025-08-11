library(igraph)
library(dplyr)
library(ggplot2)
library(corrplot)
library(pheatmap)

# Enhanced hub analysis function
analyze_subtype_hubs <- function(ppi_networks, top_n = 20) {
  
  hub_results <- list()
  
  for(subtype in names(ppi_networks)) {
    graph <- ppi_networks[[subtype]]
    
    # Calculate multiple centrality measures
    degrees <- degree(graph)
    betweenness <- betweenness(graph, normalized = TRUE)
    closeness <- closeness(graph, normalized = TRUE)
    eigenvector <- eigen_centrality(graph)$vector
    pagerank <- page_rank(graph)$vector
    
    # Get mean expression
    mean_expr <- V(graph)$mean_expr
    
    # Combine into dataframe
    centrality_df <- data.frame(
      gene = V(graph)$name,
      degree = degrees,
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector,
      pagerank = pagerank,
      mean_expr = mean_expr,
      subtype = subtype
    )
    
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

# Load your PPI networks
ppi_networks <- list(
  uro = ppi_uro,
  gu = ppi_gu, 
  basq = ppi_basq,
  mes = ppi_mes,
  scne = ppi_scne
)

# Run enhanced hub analysis
hub_analysis <- analyze_subtype_hubs(ppi_networks_final)

# Save results
save(hub_analysis, file = "data/hub_analysis_results.Rdata")
