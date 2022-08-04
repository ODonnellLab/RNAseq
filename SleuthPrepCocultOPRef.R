library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto/sample_table_cocult.txt')
metadata <- dplyr::select(metadata, c(genotype, food, additive, replicate, sample)) %>%
  dplyr::filter(sample != 'NA') %>%
  mutate(genotype = fct_relevel(genotype, "N2"),
         food = fct_relevel(food, "OP50"),
         additive = fct_relevel(additive, "none"))
#add file path for abundance data
metadata <- dplyr::mutate(metadata, path = file.path(getwd(),"kallisto",sample))
ttg <- read_delim('kallisto/transcripts_to_genes.txt',
                  col_names = c('target_id', 'ens_gene', 'ext_gene'))

# soCo <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
soCo <- sleuth_prep(metadata,
                    target_mapping = ttg,
                    aggregation_column = 'ext_gene',
                    extra_bootstrap_summary = TRUE,
                    read_bootstrap_tpm  = TRUE,
                    num_cores = 4,
                    gene_mode = TRUE)


##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, soCo let's account for that factor
soCo <- sleuth_fit(soCo, ~food, 'full')
soCo <- sleuth_fit(soCo, ~1, 'reduced') #LRT only useful for nested models

soCo <- sleuth_lrt(soCo, 'reduced', 'full')

#####
models(soCo)


### look at LRT for all transcripts:
sleuth_results_table <- sleuth_results(soCo, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(sleuth_results_table, qval <= 0.05)
head(sleuth_significant, 20)


#let's make some wald tests for the results of genotypes and additives:
soCo <- sleuth_wt(soCo, 'foodJUb39', which_model = 'full')
soCo <- sleuth_wt(soCo, 'foodMOYb033', which_model = 'full')
soCo <- sleuth_wt(soCo, 'foodCoculture', which_model = 'full')

meta_noQ1 <- filter(metadata, sample != "Q1")

soCo_noD3 <- sleuth_prep(meta_noQ1,
                         target_mapping = ttg,
                         aggregation_column = 'ext_gene',
                         extra_bootstrap_summary = TRUE,
                         read_bootstrap_tpm  = TRUE,
                         num_cores = 4,
                         gene_mode = TRUE)

soCo_noD3 <- sleuth_fit(soCo_noD3, ~food, 'full')
soCo_noD3 <- sleuth_fit(soCo_noD3, ~1, 'reduced') #LRT only useful for nested models

soCo_noD3 <- sleuth_lrt(soCo_noD3, 'reduced', 'full')

soCo_noD3 <- sleuth_wt(soCo_noD3, 'foodJUb39', which_model = 'full')
soCo_noD3 <- sleuth_wt(soCo_noD3, 'foodMOYb033', which_model = 'full')
soCo_noD3 <- sleuth_wt(soCo_noD3, 'foodCoculture', which_model = 'full')

sleuth_live(soCo_noD3)
fr