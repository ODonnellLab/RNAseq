library(tidyverse)
#setwd(file.path('E:/RNAseq/kallisto/')) # replace with kallisto wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read_delim('kallisto/sample_table_GT_retcon.txt')
metadata <- dplyr::select(metadata, c(genotype, food,replicate, sample)) %>%
  dplyr::filter(sample != 'NA') %>%
  mutate(genotype = fct_relevel(genotype, "N2"),
         food = fct_relevel(food, "OP50"))
#add file path for abundance data
metadata <- dplyr::mutate(metadata, path = file.path(getwd(),"kallisto",sample),
                          group = interaction(food,genotype))
ttg <- read_delim('kallisto/transcripts_to_genes.txt',
                  col_names = c('target_id', 'ens_gene', 'ext_gene'))

# soBWGT <- sleuth_prep(metadata, target_mapping = ttg,
#                   aggregation_column = 'ens_gene',
#                   extra_bootstrap_summary = TRUE)


########## full model ############
library(sleuth)
soBWGT <- sleuth_prep(filter(metadata, food == "tnaA"),
                     target_mapping = ttg,
                     aggregation_column = 'ext_gene',
                     extra_bootstrap_summary = TRUE,
                     read_bootstrap_tpm  = TRUE,
                     num_cores = 4,
                     gene_mode = TRUE)

soBWGT_tx <- sleuth_prep(filter(metadata, food == "tnaA"),
                      target_mapping = ttg,
                      aggregation_column = 'ext_gene',
                      extra_bootstrap_summary = TRUE,
                      read_bootstrap_tpm  = TRUE,
                      num_cores = 4,
                      gene_mode = FALSE)

soOPGT <- sleuth_prep(filter(metadata, food == "OP50"),
                      target_mapping = ttg,
                      aggregation_column = 'ext_gene',
                      extra_bootstrap_summary = TRUE,
                      read_bootstrap_tpm  = TRUE,
                      num_cores = 4,
                      gene_mode = TRUE)


##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, soBWGT let's account for that factor
soBWGT <- sleuth_fit(soBWGT, ~ 1, 'reduced') #LRT only useful for nested models
soBWGT <- sleuth_fit(soBWGT, ~ genotype, 'gt')
soBWGT <- sleuth_lrt(soBWGT, 'reduced', 'gt')
soBWGT_tx <- sleuth_fit(soBWGT, ~ genotype, 'gt')
#####
models(soBWGT)

soBWGT <- sleuth_wt(soBWGT, 'genotypecest-1.2', which_model = 'gt')
soBWGT <- sleuth_wt(soBWGT, 'genotypecest-4', which_model = 'gt')

##### including all interaction terms to make pairwise comparisons possible
### food has the largest effect indpendent of indole, soBWGT let's account for that factor
soOPGT <- sleuth_fit(soOPGT, ~ 1, 'reduced') #LRT only useful for nested models
soOPGT <- sleuth_fit(soOPGT, ~ genotype, 'gt')
soOPGT <- sleuth_lrt(soOPGT, 'reduced', 'gt')
#####
models(soOPGT)

soOPGT <- sleuth_wt(soOPGT, 'genotypecest-1.2', which_model = 'gt')
soOPGT <- sleuth_wt(soOPGT, 'genotypecest-2.1', which_model = 'gt')
soOPGT <- sleuth_wt(soOPGT, 'genotypecest-4', which_model = 'gt')


sleuth_live(soBWGT)
sleuth_live(soOPGT)

cest12BW <- sleuth_results(soBWGT, 'genotypecest-1.2', which_model = 'gt')
theme_set(theme_classic())
ggplot(cest12BW) +
  geom_point(aes(x = b, y = -log(pval), group = target_id),colour = "red") +
  gghighlight(-log10(qval) > 5) + #, label_key = target_id
  coord_cartesian(ylim = c(0,100), xlim = c(-2,2)) +
  geom_vline(xintercept = 0, linetype = 'dashed')

cest12OP <- sleuth_results(soOPGT, 'genotypecest-1.2', which_model = 'gt')
ggplot(cest12OP) +
  geom_point(aes(x = b, y = -log(pval), group = target_id),colour = "red") +
  gghighlight(-log10(qval) > 1.3) + #, label_key = target_id
  coord_cartesian(ylim = c(0,100), xlim = c(-2,2)) +
  geom_vline(xintercept = 0, linetype = 'dashed')

p <- sleuth::plot_volcano(soBWGT, 'genotypecest-1.2', test_type = "wt", which_model = "gt")
library(gghighlight)
p + gghighlight(-log10(qval) > 10, label_key = target_id)



p <- sleuth::plot_volcano(soOPGT, 'genotypecest-2.1', test_type = "wt", which_model = "gt")
library(gghighlight)


cest21OP <- sleuth_results(soOPGT, 'genotypecest-2.1', which_model = 'gt')
ggplot(cest21OP) +
  geom_point(aes(x = b, y = -log(pval), group = target_id),colour = "red") +
  gghighlight(-log10(qval) > 1.3) + #, label_key = target_id
  coord_cartesian(ylim = c(0,75), xlim = c(-2,2)) +
  geom_vline(xintercept = 0, linetype = 'dashed')

cest4OP <- sleuth_results(soOPGT, 'genotypecest-4', which_model = 'gt')
ggplot(cest4OP) +
  geom_point(aes(x = b, y = -log(pval), group = target_id),colour = "red") +
  gghighlight(-log10(qval) > 5) + #, label_key = target_id
  coord_cartesian(ylim = c(0,75), xlim = c(-2,2)) +
  geom_vline(xintercept = 0, linetype = 'dashed')

cest4BW <- sleuth_results(soBWGT, 'genotypecest-4', which_model = 'gt')
ggplot(cest4BW) +
  geom_point(aes(x = b, y = -log(pval), group = target_id),colour = "red") +
  gghighlight(-log10(qval) > 10) + #, label_key = target_id
  coord_cartesian(ylim = c(0,75), xlim = c(-3,3)) +
  geom_vline(xintercept = 0, linetype = 'dashed')
  
