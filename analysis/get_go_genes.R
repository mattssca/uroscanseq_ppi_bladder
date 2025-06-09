#load packages
library(msigdbr)
library(dplyr)

#retrieve all human GO Biological Process gene sets
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

#filter for your two GO terms by their gs_id (GO ID)
cell_cycle_genes <- msig_go %>% filter(gs_exact_source == "GO:0007049")
cell_cycle_process_genes    <- msig_go %>% filter(gs_exact_source == "GO:0022402")

#combine and deduplicate gene symbols
cell_cycle_go <- unique(c(cell_cycle_genes$gene_symbol, cell_cycle_process_genes$gene_symbol))

#save object
save(cell_cycle_genes,file = "../uroscanseq_ppi_bladder/data/GO/cell_cycle_genes.Rdata")
save(cell_cycle_process_genes,file = "../uroscanseq_ppi_bladder/data/GO/cell_cycle_process_genes.Rdata")
save(cell_cycle_go,file = "../uroscanseq_ppi_bladder/data/GO/cell_cycle_go.Rdata")
