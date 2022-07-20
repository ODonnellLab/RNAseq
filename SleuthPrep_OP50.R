library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto/sample_table_OP.txt')
metadata <- dplyr::select(metadata, c(genotype, food, additive, replicate, sample)) %>%
  dplyr::filter(sample != 'NA') %>%
  mutate(genotype = fct_relevel(genotype, "N2"),
         food = fct_relevel(food, "OP50"),
         additive = fct_relevel(additive, "none"))
#add file path for abundance data
metadata <- dplyr::mutate(metadata, path = file.path(getwd(),"kallisto",sample),
                          group = interaction(genotype,additive))
ttg <- read_delim('kallisto/transcripts_to_genes.txt',
                  col_names = c('target_id', 'ens_gene', 'ext_gene'))

# soOP <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
soOP <- sleuth_prep(metadata,
                  target_mapping = ttg,
                  aggregation_column = 'ext_gene',
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm  = TRUE,
                  num_cores = 4,
                  gene_mode = TRUE)
##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, soOP let's account for that factor
soOP <- sleuth_fit(soOP, ~ group, 'full')
soOP <- sleuth_fit(soOP, ~ 1, 'reduced') #LRT only useful for nested models
soOP <- sleuth_fit(soOP, ~genotype, 'gt')
soOP <- sleuth_fit(soOP, ~genotype + additive, 'max')
soOP <- sleuth_fit(soOP, ~additive, 'additive')
soOP <- sleuth_lrt(soOP, 'reduced', 'full')
soOP <- sleuth_lrt(soOP, 'reduced', 'gt')
soOP <- sleuth_lrt(soOP, 'gt', 'max')
#####
models(soOP)

### look at LRT for all transcripts:
sleuth_results_table <- sleuth_results(soOP, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(sleuth_results_table, qval <= 0.05)
head(sleuth_significant, 20)

#let's make some wald tests for the results of genotypes and additives:
soOP <- sleuth_wt(soOP, 'additiveOA', which_model = 'max')
soOP <- sleuth_wt(soOP, 'additive5HT', which_model = 'max')

soOP <- sleuth_wt(soOP, 'genotypecest-1.2', which_model = 'gt')
soOP <- sleuth_wt(soOP, 'genotypecest-2.1', which_model = 'gt')
soOP <- sleuth_wt(soOP, 'genotypecest-4', which_model = 'gt')

soOP <- sleuth_wt(soOP, 'additiveOA', which_model = 'additive')
soOP <- sleuth_wt(soOP, 'additive5HT', which_model = 'additive')


soOP <- sleuth_wt(soOP, 'groupcest-1.2.none', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupcest-2.1.none', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupcest-4.none', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupN2.5HT', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupcest-4.5HT', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupN2.OA', which_model = 'full')
soOP <- sleuth_wt(soOP, 'groupcest-2.1.OA', which_model = 'full')


sleuth_cest1.2sig <- dplyr::filter(sleuth_cest1.2table, qval <= 0.05)

sleuth_cest2.1table <- sleuth_results(soOP, 'genotypecest-2.1', 'wt', show_all = FALSE)
sleuth_cest2.1sig <- dplyr::filter(sleuth_cest2.1table, qval <= 0.05)

sleuth_cest4table <- sleuth_results(soOP, 'genotypecest-4', 'wt', show_all = FALSE)
sleuth_cest4sig <- dplyr::filter(sleuth_cest4table, qval <= 0.05)

sleuth_live(soOP)
