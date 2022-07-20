library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto/sample_table_old.txt')
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

# soOld <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
soOld <- sleuth_prep(metadata,
                  target_mapping = ttg,
                  aggregation_column = 'ext_gene',
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm  = TRUE,
                  num_cores = 4,
                  gene_mode = TRUE)


##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, soOld let's account for that factor
soOld <- sleuth_fit(soOld, ~food + group, 'full')
soOld <- sleuth_fit(soOld, ~food, 'reduced') #LRT only useful for nested models
soOld <- sleuth_fit(soOld, ~food + genotype, 'gt')
soOld <- sleuth_fit(soOld, ~food + genotype + additive, 'max')
soOld <- sleuth_fit(soOld, ~food + additive, 'additive')
soOld <- sleuth_lrt(soOld, 'reduced', 'full')
soOld <- sleuth_lrt(soOld, 'reduced', 'gt')
soOld <- sleuth_lrt(soOld, 'gt', 'max')
#####
models(soOld)


### look at LRT for all transcripts:
sleuth_results_table <- sleuth_results(soOld, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(sleuth_results_table, qval <= 0.05)
head(sleuth_significant, 20)


#let's make some wald tests for the results of genotypes and additives:
soOld <- sleuth_wt(soOld, 'additiveTA', which_model = 'max')
soOld <- sleuth_wt(soOld, 'additiveBW', which_model = 'max')
soOld <- sleuth_wt(soOld, 'additiveOA', which_model = 'max')
soOld <- sleuth_wt(soOld, 'additive5HT', which_model = 'max')

soOld <- sleuth_wt(soOld, 'genotypecest-1.2', which_model = 'gt')
soOld <- sleuth_wt(soOld, 'genotypecest-2.1', which_model = 'gt')
soOld <- sleuth_wt(soOld, 'genotypecest-4', which_model = 'gt')

soOld <- sleuth_wt(soOld, 'additiveTA', which_model = 'additive')
soOld <- sleuth_wt(soOld, 'additiveBW', which_model = 'additive')
soOld <- sleuth_wt(soOld, 'additiveOA', which_model = 'additive')
soOld <- sleuth_wt(soOld, 'additive5HT', which_model = 'additive')

soOld <- sleuth_wt(soOld, 'groupN2.TA:groupcest-1.2.TA', which_model = 'full')

sleuth_cest1.2sig <- dplyr::filter(sleuth_cest1.2table, qval <= 0.05)

sleuth_cest2.1table <- sleuth_results(soOld, 'genotypecest-2.1', 'wt', show_all = FALSE)
sleuth_cest2.1sig <- dplyr::filter(sleuth_cest2.1table, qval <= 0.05)

sleuth_cest4table <- sleuth_results(soOld, 'genotypecest-4', 'wt', show_all = FALSE)
sleuth_cest4sig <- dplyr::filter(sleuth_cest4table, qval <= 0.05)

sleuth_live(soOld)
