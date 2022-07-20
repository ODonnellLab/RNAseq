library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto//sample_table.txt')
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

# so <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
so <- sleuth_prep(metadata,
                  target_mapping = ttg,
                  aggregation_column = 'ext_gene',
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm  = TRUE,
                  num_cores = 4,
                  gene_mode = TRUE)
##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, so let's account for that factor
so <- sleuth_fit(so, ~food + group, 'full')
so <- sleuth_fit(so, ~food, 'reduced') #LRT only useful for nested models
so <- sleuth_fit(so, ~food + genotype, 'gt')
so <- sleuth_fit(so, ~food + genotype + additive, 'max')
so <- sleuth_fit(so, ~food + additive, 'additive')
so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_lrt(so, 'reduced', 'gt')
so <- sleuth_lrt(so, 'gt', 'max')
#####
models(so)


### look at LRT for all transcripts:
sleuth_results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(sleuth_results_table, qval <= 0.05)
head(sleuth_significant, 20)


#let's make some wald tests for the results of genotypes and additives:
so <- sleuth_wt(so, 'additiveTA', which_model = 'max')
so <- sleuth_wt(so, 'additiveBW', which_model = 'max')
so <- sleuth_wt(so, 'additiveOA', which_model = 'max')
so <- sleuth_wt(so, 'additive5HT', which_model = 'max')

so <- sleuth_wt(so, 'genotypecest-1.2', which_model = 'gt')
so <- sleuth_wt(so, 'genotypecest-2.1', which_model = 'gt')
so <- sleuth_wt(so, 'genotypecest-4', which_model = 'gt')

so <- sleuth_wt(so, 'additiveTA', which_model = 'additive')
so <- sleuth_wt(so, 'additiveBW', which_model = 'additive')
so <- sleuth_wt(so, 'additiveOA', which_model = 'additive')
so <- sleuth_wt(so, 'additive5HT', which_model = 'additive')

so <- sleuth_wt(so, 'groupN2.TA:groupcest-1.2.TA', which_model = 'full')

sleuth_cest1.2sig <- dplyr::filter(sleuth_cest1.2table, qval <= 0.05)

sleuth_cest2.1table <- sleuth_results(so, 'genotypecest-2.1', 'wt', show_all = FALSE)
sleuth_cest2.1sig <- dplyr::filter(sleuth_cest2.1table, qval <= 0.05)

sleuth_cest4table <- sleuth_results(so, 'genotypecest-4', 'wt', show_all = FALSE)
sleuth_cest4sig <- dplyr::filter(sleuth_cest4table, qval <= 0.05)

sleuth_live(so)
