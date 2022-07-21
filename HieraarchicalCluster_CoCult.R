library(tidyverse)
# co-culture
# these values are relative to co-culture, so we're looking for unique changes in coculture
JUbdown <- read_csv('~/Downloads/volcano_table (1).csv')
MOYbdown <- read_csv('~/Downloads/volcano_table (2).csv')
JUbup <- read_csv('~/Downloads/volcano_table (3).csv')
MOYbup <- read_csv('~/Downloads/volcano_table (4).csv')

SigDownCoculture <- inner_join(JUbup[,3], MOYbup[,3])
SigUpCoculture <- inner_join(JUbdown[,3], MOYbdown[,3])

#eliminating sample Q1 (outlier)
JUbdown2 <- read_csv('~/Downloads/volcano_table (5).csv')
MOYbdown2 <- read_csv('~/Downloads/volcano_table (6).csv')
JUbup2 <- read_csv('~/Downloads/volcano_table (7).csv')
MOYbup2 <- read_csv('~/Downloads/volcano_table (8).csv')

SigDownCoculture2 <- inner_join(JUbup2[,3], MOYbup2[,3])
SigUpCoculture2 <- inner_join(JUbdown2[,3], MOYbdown2[,3])


c("cpr-5", "nspc-17", "pck-1", "K09C6.9", "far-3","C45B2.1", "grd-5", "nlp-24")


cluster1 <- c("Y39B6A.1", "cpr-4", "mthf-1", "aagr-2", "icl-1", "hach-1", "ech-6", "cpr-5", "mtrr-1", "C05C10.3", "gln-1", "VF13D12L.3", "smd-1", "F54D5.12", "asp-13", "pmp-5", "folt-2", "hphd-1", "clec-50", "clec-41", "ZK6.11", "T01D3.6", "lec-10", "lys-1", "pud-2.2", "pud-2.1", "pud-1.2", "F28B4.3", "clec-63", "nep-17", "asp-5", "pyc-1", "metr-1", "msra-1", "acdh-1")
cluster1WBgene <-ttg %>% filter(ext_gene %in% cluster1) %>% select(ens_gene) %>% unique()

gprofiler2::gost(cluster1WBgene[[1]], organism = "celegans")

gprofiler2::gost(cluster1WBgene[33:35,][[1]], organism = "celegans")
gprofiler2::gost(cluster1WBgene[17:32,][[1]], organism = "celegans")
gprofiler2::gost(cluster1WBgene[1:17,][[1]], organism = "celegans")

cluster1.1 <- c("spp-8", "aagr-1", "cysl-2", "C01B10.6", "abt-4", "F54E2.1", "C05D12.3", "H34I24.2", "R193.2", "gst-13", "gst-10", "clec-83", "clec-186", "asah-2", "K08D8.6", "hpo-15", "ugt-44", "F42A10.7", "asah-1", "dod-24", "F54D5.4", "clec-65", "lys-8", "dod-19", "zip-2", "C12D12.1", "lys-2", "asp-14", "pho-11", "endu-2", "T19D12.4", "C17H12.8")
cluster1.1WBgene <-ttg %>% filter(ext_gene %in% cluster1.1) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster1.1WBgene[[1]], organism = "celegans")


cluster2 <- c("clec-47", "acdh-2", "F21C10.9", "acs-2", "F18E3.11", "F15E6.4", "C09B8.4", "lys-6", "W06A7.4", "cyp-35A3", "ugt-53", "F25D1.5", "nspc-17", "F44A6.5", "F18E3.13", "F23F12.12", "zig-3", "B0205.13", "arf-1.1", "Y38H6C.21", "Y34F4.2", "F18E3.12", "lys-5", "C07G1.7")
cluster2WBgene <-ttg %>% filter(ext_gene %in% cluster2) %>% select(ens_gene) %>% unique()

gprofiler2::gost(cluster2WBgene[[1]], organism = "celegans")

cluster3 <- c("T27C5.8", "comt-2", "F35E8.10", "Y53G8AM.5", "F49H6.13", "K10D11.2", "T24C4.8", "B0348.2", "clec-42", "C49G7.7", "C01G10.5", "tsp-1", "F20G2.5", "Y49G5A.1", "tsp-2", "kreg-1", "lys-3", "C17H12.6", "ech-9", "clec-45")
cluster3WBgene <- ttg %>% filter(ext_gene %in% cluster3) %>% select(ens_gene) %>% unique()

gprofiler2::gost(cluster3WBgene[[1]], organism = "celegans")

gprofiler2::gost(filter(ttg, ext_gene %in% c("dod-21", "C32H11.9", "clec-45"))[[1]], organism = "celegans")

cluster4 <- c("swt-3", "mdl-1", "asp-12", "Y41C4A.32", "dgat-2", "D1054.8", "clec-10", "W02H5.8", "fbxa-72", "ugt-63", "ugt-26", "drd-1", "T07E3.4", "pgp-1", "fat-5", "pho-13", "nhr-68", "K08D8.3", "ugt-17", "W01B11.6", "B0410.3", "cyp-33", "C8ZK673.1", "cth-1ZK593.3", "ugt-19", "mxl-3", "C53A3.2", "DH11.2", "ddo-2")
cluster4WBgene <- ttg %>% filter(ext_gene %in% cluster4) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster4WBgene[[1]], organism = "celegans")

cluster5 <- c("F01D5.1", "asp-17", "F53A9.6", "trx-3", "pho-9", "valv-1", "gst-16", "F53A9.1", "F46A8.7", "fbxa-14", "cln-3.3", "K11H12.5", "hch-1", "oac-6", "K10D11.6", "F39E9.1", "C18H7.11", "gst-22", "F23F12.3", "ugt-39", "oac-14", "F10D2.10")
cluster5WBgene <- ttg %>% filter(ext_gene %in% cluster5) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster5WBgene[[1]], organism = "celegans")

cluster6 <- c("C50F7.5", "Y47H10A.5", "K08D8.4", "clec-4", "K08D8.5", "fbxa-24", "Y17G7B.8", "C08F11.13", "K10D11.5", "mth-1", "dod-17", "cpt-4", "K06A9.1", "T24C4.4", "ech-1.1", "dct-17", "cld-9", "F01D5.1")
cluster6WBgene <- ttg %>% filter(ext_gene %in% cluster6) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster6WBgene[[1]], organism = "celegans")

cluster7 <- c("Y71G12B.2", "W09G12.7", "cyp-13B1", "zip-10", "tag-244", "F33H12.7", "scrm-4", "gst-12", "F27E5.9", "C51E3.10", "K10D11.3", "math-45", "gst-30", "ifa-3", "F01D5.2", "K03H6.2", "cyp-32", "B1C08B6.4", "Y47H9C.1", "F53A9.7", "W02A2.9", "gst-38", "T26H5.9", "F13A7.11", "Y65B4BR.1", "cyp-37B1", "lipl-2", "drd-50", "sdz-24", "lipl-1")
cluster7WBgene <- ttg %>% filter(ext_gene %in% cluster7) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster7WBgene[[1]], organism = "celegans")

cluster8 <- c("irg-5", "B0024.4", "irg-4", "dod-22", "asm-3", "F53A9.8", "C32H11.4", "pud-4", "T24B8.5", "clec-67", "pud-3")
cluster8WBgene <- ttg %>% filter(ext_gene %in% cluster8) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster8WBgene[[1]], organism = "celegans")

cluster9 <- c("F45D11.14", "cpt-5", "hpo-6", "clec-85", "gst-20", "K11H12.4", "cpr-3", "tth-1", "clec-62", "D1086.3", "clec-265", "cbs-1", "C16B8.3", "W04B5.3", "pcp-4", "Y51H4A.25", "F55G11.4", "F45D3.4", "F53C11.1", "C49C3.9")
cluster9WBgene <- ttg %>% filter(ext_gene %in% cluster9) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster9WBgene[[1]], organism = "celegans")

cluster10 <- c("lec-11", "clec-66", "M02H5.8", "dpy-14", "spp-18", "gstk-1", "hbl-1", "Y47G6A.5", "gba-1", "F49C12.7", "clc-1", "clec-84", "scl-2", "F09C8.1", "Y22D7AL.15", "clec-78", "Y69A2AL.2", "asm-2", "gst-28", "C50B6.7", "Y46D2A.2")
cluster10WBgene <- ttg %>% filter(ext_gene %in% cluster10) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster10WBgene[[1]], organism = "celegans")

cluster11 <- c("sodh-1", "T13F3.6", "mpc-1", "dhs-25", "C23H5.8", "lys-7", "mct-3", "asp-10", "mmaa-1", "ugt-62", "Y37A1B.5", "papl-1", "T15B7.1", "T05E12.6", "nhr-114", "C30G12.2", "atf-5", "R08E5.3")
cluster11WBgene <- ttg %>% filter(ext_gene %in% cluster11) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster11WBgene[[1]], organism = "celegans")

cluster12 <- c("Y119C1B.5", "F15B9.10", "F17H10.1", "cyp-35A2", "kel-8", "F32D8.12", "clec-5", "bca-1", "asns-2", "F17A9.4", "ZK822.5", "clec-48", "thn-2", "C35A5.3", "sodh-1")
cluster12WBgene <- ttg %>% filter(ext_gene %in% cluster12) %>% select(ens_gene) %>% unique()
gprofiler2::gost(cluster12WBgene[[1]], organism = "celegans")


GO_all <- gprofiler2::gost(sleuth_significant$ens_gene, organism = 'celegans')


