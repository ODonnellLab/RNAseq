library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto/sample_table_BW.txt')
metadata <- dplyr::select(metadata, c(genotype, food, additive, replicate, sample)) %>%
  dplyr::filter(sample != 'NA') %>%
  mutate(genotype = fct_relevel(genotype, "N2"),
         additive = fct_relevel(additive, "none"))
#add file path for abundance data
metadata <- dplyr::mutate(metadata, path = file.path(getwd(),"kallisto",sample),
                          group = interaction(genotype,additive))
ttg <- read_delim('kallisto/transcripts_to_genes.txt',
                  col_names = c('target_id', 'ens_gene', 'ext_gene'))

# sotnaA <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
sotnaA <- sleuth_prep(metadata,
                    target_mapping = ttg,
                    aggregation_column = 'ext_gene',
                    extra_bootstrap_summary = TRUE,
                    read_bootstrap_tpm  = TRUE,
                    num_cores = 4,
                    gene_mode = TRUE)
##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, sotnaA let's account for that factor
sotnaA <- sleuth_fit(sotnaA, ~ group, 'full')
sotnaA <- sleuth_fit(sotnaA, ~ 1, 'reduced') #LRT only useful for nested models
sotnaA <- sleuth_fit(sotnaA, ~genotype, 'gt')
sotnaA <- sleuth_fit(sotnaA, ~genotype + additive, 'max')
sotnaA <- sleuth_fit(sotnaA, ~additive, 'additive')
sotnaA <- sleuth_lrt(sotnaA, 'reduced', 'full')
sotnaA <- sleuth_lrt(sotnaA, 'reduced', 'gt')
sotnaA <- sleuth_lrt(sotnaA, 'gt', 'max')
#####
models(sotnaA)

### look at LRT for all transcripts:
sleuth_results_table <- sleuth_results(sotnaA, 'gt:max', test_type = 'lrt')
sleuth_significant <- dplyr::filter(sleuth_results_table, qval <= 0.05)
head(sleuth_significant, 20)

#let's make some wald tests for the results of genotypes and additives:
sotnaA <- sleuth_wt(sotnaA, 'additiveTA', which_model = 'max')
sotnaA <- sleuth_wt(sotnaA, 'additiveBW', which_model = 'max')

sotnaA <- sleuth_wt(sotnaA, 'genotypecest-1.2', which_model = 'gt')

sotnaA <- sleuth_wt(sotnaA, 'additiveBW', which_model = 'additive')
sotnaA <- sleuth_wt(sotnaA, 'additiveTA', which_model = 'additive')


sotnaA <- sleuth_wt(sotnaA, 'groupcest-1.2.none', which_model = 'full')
sotnaA <- sleuth_wt(sotnaA, 'groupcest-1.2.BW', which_model = 'full')
sotnaA <- sleuth_wt(sotnaA, 'groupN2.BW', which_model = 'full')
sotnaA <- sleuth_wt(sotnaA, 'groupN2.TA', which_model = 'full')
sotnaA <- sleuth_wt(sotnaA, 'groupcest-1.2.TA', which_model = 'full')


sleuth_cest1.2sig <- dplyr::filter(sleuth_cest1.2table, qval <= 0.05)

sleuth_cest2.1table <- sleuth_results(sotnaA, 'genotypecest-2.1', 'wt', show_all = FALSE)
sleuth_cest2.1sig <- dplyr::filter(sleuth_cest2.1table, qval <= 0.05)

sleuth_cest4table <- sleuth_results(sotnaA, 'genotypecest-4', 'wt', show_all = FALSE)
sleuth_cest4sig <- dplyr::filter(sleuth_cest4table, qval <= 0.05)

sleuth_live(sotnaA)