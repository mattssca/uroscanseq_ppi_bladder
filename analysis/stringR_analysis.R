library(STRINGdb)
library(igraph)
library(scales)

run_stringdb <- function(expr_data = NULL, return_cytoscape = FALSE) {

  # Initialize STRINGdb (for human, species = 9606)
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
  
  # Map gene symbols to STRING IDs
  genes <- rownames(expr_data)
  mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  
  # Get interactions for mapped STRING IDs
  interactions <- string_db$get_interactions(mapped$STRING_id)
  
  # Create lookup: STRING ID -> gene symbol
  string_to_symbol <- setNames(mapped$gene, mapped$STRING_id)
  
  if(return_cytoscape){
    # Add gene symbols to the interactions data frame
    interactions$from_gene <- string_to_symbol[interactions$from]
    interactions$to_gene   <- string_to_symbol[interactions$to]
    return(interactions)
  }
  
  # Create igraph object
  graph_obj <- graph_from_data_frame(interactions, directed = FALSE)
  
  # Calculate mean expression for each gene
  mean_expr <- rowMeans(expr_data, na.rm = TRUE)
  
  # Assign mean expression to nodes by mapping STRING IDs to gene symbols
  V(graph_obj)$mean_expr <- mean_expr[string_to_symbol[V(graph_obj)$name]]
  
  return(graph_obj)
}

#run function for each subtype subset
ppi_uro = run_stringdb(expr_data = expr_subtypes$uro)
ppi_gu = run_stringdb(expr_data = expr_subtypes$gu)
ppi_basq = run_stringdb(expr_data = expr_subtypes$basq)
ppi_mes = run_stringdb(expr_data = expr_subtypes$mes)
ppi_scne = run_stringdb(expr_data = expr_subtypes$scne)

#return for cytoscape
cyto_uro = run_stringdb(expr_data = expr_subtypes$uro, return_cytoscape = TRUE)
write.table(cyto_uro, file = "data/cytoscape/uro_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cyto_gu = run_stringdb(expr_data = expr_subtypes$gu, return_cytoscape = TRUE)
write.table(cyto_gu, file = "data/cytoscape/gu_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cyto_basq = run_stringdb(expr_data = expr_subtypes$basq, return_cytoscape = TRUE)
write.table(cyto_basq, file = "data/cytoscape/basq_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cyto_mes = run_stringdb(expr_data = expr_subtypes$mes, return_cytoscape = TRUE)
write.table(cyto_mes, file = "data/cytoscape/mes_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cyto_scne = run_stringdb(expr_data = expr_subtypes$scne, return_cytoscape = TRUE)
write.table(cyto_scne, file = "data/cytoscape/scne_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#write node data with mean expresion values
expr_uro <- expr_subtypes$uro
mean_expr <- rowMeans(expr_uro, na.rm = TRUE)
node_table <- data.frame(gene_symbol = rownames(expr_uro), mean_expr = mean_expr)
write.csv(node_table, "data/cytoscape/uro_node_table.csv", row.names = FALSE)

expr_gu <- expr_subtypes$gu
mean_expr <- rowMeans(expr_gu, na.rm = TRUE)
node_table <- data.frame(gene_symbol = rownames(expr_gu), mean_expr = mean_expr)
write.csv(node_table, "data/cytoscape/gu_node_table.csv", row.names = FALSE)

expr_basq <- expr_subtypes$basq
mean_expr <- rowMeans(expr_basq, na.rm = TRUE)
node_table <- data.frame(gene_symbol = rownames(expr_basq), mean_expr = mean_expr)
write.csv(node_table, "data/cytoscape/basq_node_table.csv", row.names = FALSE)

expr_mes <- expr_subtypes$mes
mean_expr <- rowMeans(expr_mes, na.rm = TRUE)
node_table <- data.frame(gene_symbol = rownames(expr_mes), mean_expr = mean_expr)
write.csv(node_table, "data/cytoscape/mes_node_table.csv", row.names = FALSE)

expr_scne <- expr_subtypes$scne
mean_expr <- rowMeans(expr_scne, na.rm = TRUE)
node_table <- data.frame(gene_symbol = rownames(expr_scne), mean_expr = mean_expr)
write.csv(node_table, "data/cytoscape/scne_node_table.csv", row.names = FALSE)

#igraph objects
ppi_gu = run_stringdb(expr_data = expr_subtypes$gu)
ppi_basq = run_stringdb(expr_data = expr_subtypes$basq)
ppi_mes = run_stringdb(expr_data = expr_subtypes$mes)
ppi_scne = run_stringdb(expr_data = expr_subtypes$scne)

#save PPIs
save(ppi_uro, file = "data/stringR_outs/ppi_uro.Rdata")
save(ppi_gu, file = "data/stringR_outs/ppi_gu.Rdata")
save(ppi_basq, file = "data/stringR_outs/ppi_basq.Rdata")
save(ppi_mes, file = "data/stringR_outs/ppi_mes.Rdata")
save(ppi_scne, file = "data/stringR_outs/ppi_scne.Rdata")

#plot PPI
plot(ppi_uro, vertex.label = NA, vertex.size = 5)

#overlay with gene expression
#uro
pdf("figs/ppi_uro.pdf", width = 10, height = 10)
node_colors <- col_numeric("Blues", domain = range(V(ppi_uro)$mean_expr, na.rm = TRUE))(V(ppi_uro)$mean_expr)
node_colors[is.na(node_colors)] <- "grey"
plot(ppi_uro, vertex.color = node_colors, vertex.size = 5, vertex.label = NA)
dev.off()

#gu
pdf("figs/ppi_gu.pdf", width = 10, height = 10)
node_colors <- col_numeric("Blues", domain = range(V(ppi_gu)$mean_expr, na.rm = TRUE))(V(ppi_gu)$mean_expr)
node_colors[is.na(node_colors)] <- "grey"
plot(ppi_gu, vertex.color = node_colors, vertex.size = 5, vertex.label = NA)
dev.off()

#basq
pdf("figs/ppi_basq.pdf", width = 10, height = 10)
node_colors <- col_numeric("Blues", domain = range(V(ppi_basq)$mean_expr, na.rm = TRUE))(V(ppi_basq)$mean_expr)
node_colors[is.na(node_colors)] <- "grey"
plot(ppi_basq, vertex.color = node_colors, vertex.size = 5, vertex.label = NA)
dev.off()

#mes
pdf("figs/ppi_mes.pdf", width = 10, height = 10)
node_colors <- col_numeric("Blues", domain = range(V(ppi_mes)$mean_expr, na.rm = TRUE))(V(ppi_mes)$mean_expr)
node_colors[is.na(node_colors)] <- "grey"
plot(ppi_mes, vertex.color = node_colors, vertex.size = 5, vertex.label = NA)
dev.off()

#scne
pdf("figs/ppi_scne.pdf", width = 10, height = 10)
node_colors <- col_numeric("Blues", domain = range(V(ppi_scne)$mean_expr, na.rm = TRUE))(V(ppi_scne)$mean_expr)
node_colors[is.na(node_colors)] <- "grey"
plot(ppi_scne, vertex.color = node_colors, vertex.size = 5, vertex.label = NA)
dev.off()
