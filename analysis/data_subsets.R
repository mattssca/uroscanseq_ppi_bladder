#load data
load("C:/Users/matts/Desktop/uroscanseq_ppi_bladder/data/GO/cell_cycle_go.Rdata")
load("C:/Users/matts/Desktop/uroscanseq_ppi_bladder/data/expr/expr_mat_uroscanseq.Rdata")

#print head
head(cell_cycle_go)
head(expr_mat_uroscanseq)[1:10]

#subset gene expressiond ata to genes of interest
expr_cell_cycle = expr_mat_uroscanseq %>% 
  rownames_to_column("genes") %>% 
  filter(genes %in% cell_cycle_go) %>% 
  column_to_rownames("genes")

#save expr subset
save(expr_cell_cycle, file = "data/expr/expr_cell_cycle.Rdata")

#run predictor on full data
predicted = LundTax2023Classifier::lundtax_predict_sub(this_data = expr_mat_uroscanseq, 
                                                       log_transform = FALSE, 
                                                       adjust = TRUE, 
                                                       impute = TRUE)

#export prediction set
save(predicted, file = "data/predicted/predicted_full.Rdata")

#get subtype information
these_subtypes = as.data.frame(predicted$predictions_5classes) %>% 
  rownames_to_column("sample_id") %>% 
  rename(subtype = `predicted$predictions_5classes`)

#export subtype information as its own object
save(these_subtypes, file = "data/predicted/these_subtypes.Rdata")

#get sample sets
uro_samples = these_subtypes %>% 
  filter(subtype == "Uro") %>% 
  pull(sample_id)

gu_samples = these_subtypes %>% 
  filter(subtype == "GU") %>% 
  pull(sample_id)

basq_samples = these_subtypes %>% 
  filter(subtype == "BaSq") %>% 
  pull(sample_id)

mes_samples = these_subtypes %>% 
  filter(subtype == "Mes") %>% 
  pull(sample_id)

scne_samples = these_subtypes %>% 
  filter(subtype == "ScNE") %>% 
  pull(sample_id)

subtypes_samples = list(
  uro = uro_samples,
  gu = gu_samples,
  basq = basq_samples,
  mes = mes_samples,
  scne = scne_samples
)

#save subtype samples object
save(subtypes_samples, file = "data/subtypes_samples.Rdata")

#subset to sample groups
epxr_uro = expr_cell_cycle %>% 
  dplyr::select(all_of(uro_samples))

expr_gu = expr_cell_cycle %>% 
  dplyr::select(all_of(gu_samples))

epxr_basq = expr_cell_cycle %>% 
  dplyr::select(all_of(basq_samples))

epxr_mes = expr_cell_cycle %>% 
  dplyr::select(all_of(mes_samples))

epxr_scne = expr_cell_cycle %>% 
  dplyr::select(all_of(scne_samples))

expr_subtypes = list(
  uro = epxr_uro,
  gu = expr_gu,
  basq = epxr_basq,
  mes = epxr_mes,
  scne = epxr_scne
)

#export expr_subtypes
save(expr_subtypes, file = "data/expr/expr_subtypes.Rdata")
  