reticulate::py_install('pandas', pip=TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")
install.packages('msigdbr')

library(msigdbr)
library(fgsea)
library(reticulate)
library(tidyverse)


# gene sets
h_gene_sets = msigdbr(species = 'Homo sapiens', category = 'H')
o_gene_sets = msigdbr(species = 'Homo sapiens', category = 'C6')
i_gene_sets = msigdbr(species = 'Homo sapiens', category = 'C7')

h_list <- split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
o_list <- split(x = o_gene_sets$gene_symbol, f = o_gene_sets$gs_name)
i_list <- split(x = i_gene_sets$gene_symbol, f = i_gene_sets$gs_name)

all_list <- c(h_list, o_list, i_list)

# import pandas
pd <- import('pandas')
ht_df <- pd$read_csv('/Users/nabilazulkapeli/Documents/Honours Thesis 2025/nabs_data/DEG_figures/ht_DEG.csv')
p_df <- pd$read_csv('/Users/nabilazulkapeli/Documents/Honours Thesis 2025/nabs_data/DEG_figures/peri_DEG.csv')

# make stat lists
ht_f <- ht_df %>%
dplyr::select(variable, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(variable) %>%
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))
ht_f
p_f <- p_df %>%
  dplyr::select(variable, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(variable) %>%
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))
p_f

ht_ranks <- deframe(ht_f)
p_ranks <- deframe(p_f)
head(ht_ranks, 20)
# fgsea
set.seed(42)
ht_res <- fgsea(pathways = all_list,
                  stats = ht_ranks,
                  minSize = 1,
                  maxSize = 1000
                  )
ht_res_c <- ht_res %>%
  as_tibble() %>%
  arrange(padj)
ht_res_c
ht_csv <- apply(ht_res_c,2,as.character)
write.csv(ht_csv, '/Users/nabilazulkapeli/Documents/Honours Thesis 2025/nabs_data/DEG_figures/ht_fgsea.csv', row.names=FALSE)

p_res <- fgsea(pathways = h_list,
               stats = p_ranks,
               minSize=1,
               maxSize=1000)
p_res_c <- p_res %>%
  as_tibble() %>%
  arrange(padj)
p_res_c
p_csv <- apply(p_res_c,2,as.character)
write.csv(p_csv, '/Users/nabilazulkapeli/Documents/Honours Thesis 2025/nabs_data/DEG_figures/p_fgsea.csv', row.names=FALSE)