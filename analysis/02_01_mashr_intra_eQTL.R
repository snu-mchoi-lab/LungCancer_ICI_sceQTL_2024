# Lung Cancer eQTL - mashr pipeline
# 20240311 Hyungtai Sim

vlib = c("tidyverse", "mashr", "ashr", "ggpubr", "data.table", 
         "future.apply", "arrow")
lapply(vlib, require, character.only = TRUE, quietly = TRUE) |> suppressMessages()

base.dir = ""

setwd(paste0(base.dir, "202404-sceQTLv7/"))


## create data for one gene

fname_meta = paste0(base.dir, 'analysis/assets/file_meta.txt')
df_meta = read_delim(fname_meta, col_names = c("levels", "time", "cluster_name")) %>% 
  mutate(prefix_file = paste0(levels, "_", time, "_", cluster_name))

df_meta_target = df_meta %>% filter(levels == "anno_l1")

df_meta_significant = df_meta_target %>%
  mutate(fname = paste0("01_calling_eQTL/01_01_tensorqtl_out/", prefix_file,".map_cis.txt.gz"))

# process significant

list_significant = lapply(df_meta_significant$fname, read_delim)
names(list_significant) = df_meta_significant$prefix_file 

# create df_significant, using names
# we don't filter in this scheme.
df_significant = list_significant %>% bind_rows(.id = "cluster") %>% 
  select(phenotype_id, variant_id, cluster, slope, slope_se) %>%
  mutate(name = paste0(phenotype_id, "_", variant_id)) %>%
  distinct(name)


df_meta_target = expand.grid(
              prefix_file = df_meta_target$prefix_file, 
              chr = paste0("chr", c(1:22))) %>%
  left_join(df_meta_target, .) %>%
  mutate(fname = paste0("01_calling_eQTL/01_02_tensorqtl_nominal/", prefix_file, ".cis_qtl_pairs.",chr,".parquet"))

list_meta_target = df_meta_target %>% group_by(prefix_file) %>% group_split()


library(doFuture)
registerDoFuture()
plan(multisession, workers = 16)

list_parquet = foreach(i = c(1:length(list_meta_target))) %dopar% {
  res = lapply(list_meta_target[[i]]$fname, read_parquet) %>%
    bind_rows() %>%
    select(phenotype_id, variant_id, slope, slope_se, pval_nominal) %>%
    mutate(name = paste0(phenotype_id, "_", variant_id),
           condition =  unique(df_meta_target$prefix_file)[i]) %>%
    select(-variant_id)
}

registerDoSEQ()
plan(sequential)
df_parquet = list_parquet %>% bind_rows()



## strong.set 

df_parquet_strong = df_parquet %>% filter(name %in% df_significant$name)

set.seed(42)

df_random = df_parquet %>% distinct(name) %>% sample_n(1e5)

df_parquet_random = df_parquet %>% filter(name %in% df_random$name)
df_parquet_random %>% write_parquet("02_mashr/random_signals.intra_eqtl_top.parquet")

rm(df_parquet, list_parquet); gc()

data_mashr_strong = list()

c_list_strong =
  df_parquet_strong %>% 
  # in this scheme, we really don't filter at all
  #group_by(phenotype_id) %>% 
  # filter by strongest signals
  #mutate(min_pval_nominal = min(pval_nominal)) %>% 
  #filter(pval_nominal == min_pval_nominal) %>% 
  #select(-min_pval_nominal)  %>% 
  # remove signals
  #distinct(phenotype_id, .keep_all= TRUE) %>%
  pull(name)


df_parquet_strong = df_parquet_strong %>% filter(name %in% c_list_strong)


df_parquet_strong %>% write_parquet("02_mashr/strong_signals.intra_eqtl_allposssible.parquet")


df_parquet_strong = read_parquet("02_mashr/strong_signals.intra_eqtl_allposssible.parquet") 
df_parquet_random = read_parquet("02_mashr/random_signals.intra_eqtl_top.parquet")

# create mashr data

mashr_data_strong = list()

mashr_data_strong$Bhat = df_parquet_strong %>% select(name, condition, slope) %>%
  pivot_wider(names_from = condition, values_from = slope, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()

mashr_data_strong$Shat = df_parquet_strong %>% select(name, condition, slope_se) %>%
  pivot_wider(names_from = condition, values_from = slope_se, values_fill = 100) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()

mashr_data_random = list()

mashr_data_random$Bhat = df_parquet_random %>% select(name, condition, slope) %>%
  pivot_wider(names_from = condition, values_from = slope, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()

mashr_data_random$Shat = df_parquet_random %>% select(name, condition, slope_se) %>%
  pivot_wider(names_from = condition, values_from = slope_se, values_fill = 100) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()



# estimate Vhat from random
mashr_data_temp = mash_set_data(mashr_data_random$Bhat,mashr_data_random$Shat)
Vhat = estimate_null_correlation_simple(mashr_data_temp)


#
data.random = mash_set_data(mashr_data_random$Bhat,mashr_data_random$Shat,V=Vhat)
data.strong = mash_set_data(mashr_data_strong$Bhat,mashr_data_strong$Shat, V=Vhat)

U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
c(U.ed,U.c) %>% saveRDS("02_mashr/ulist.intra_eqtl_allpossible.RDS")

# update
Vhat_em = mash_estimate_corr_em(mashr_data_temp, U.c, details = TRUE)

data.random = mash_set_data(mashr_data_random$Bhat,mashr_data_random$Shat,V=Vhat_em$V)
data.strong = mash_set_data(mashr_data_strong$Bhat,mashr_data_strong$Shat, V=Vhat_em$V)


#estimated mixture component proportions, 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1, verbose = TRUE)

m %>% saveRDS("02_mashr/m_mixture_component.em.intra_eqtl_allpossible.RDS")
m = readRDS("02_mashr/m_mixture_component.em.intra_eqtl_allpossible.RDS")


m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE, verbose = TRUE)


m2 %>% saveRDS("02_mashr/m2_intra_eqtl_results.em.intra_eqtl_allpossible.RDS")
m2 = readRDS("02_mashr/m2_intra_eqtl_results.em.intra_eqtl_allpossible.RDS")


mashr_significant = m2 %>% get_lfsr %>% as.data.frame() %>% rownames_to_column() %>%
  pivot_longer(cols = c(2:ncol(.))) %>% filter(value < 0.05) %>%
  rename(name = rowname, condition= name)

df_parquet_strong %>% left_join(mashr_significant) %>% filter(value < 0.05) %>%
  group_by( condition) %>% summarise(n =n())

m2 %>% get_pairwise_sharing(lfsr_thresh = 0.05) %>% pheatmap()

mash_compute_posterior_matrices(g = get_fitted_g(m), m2)


## no_em

data.random = mash_set_data(mashr_data_random$Bhat, mashr_data_random$Shat, V=Vhat)
data.strong = mash_set_data(mashr_data_strong$Bhat, mashr_data_strong$Shat, V=Vhat)

U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
c(U.ed,U.c) %>% saveRDS("02_mashr/ulist.intra_eqtl_allpossible.RDS")

#estimated mixture component proportions, 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1, verbose = TRUE)

m %>% saveRDS("02_mashr/m_mixture_component.intra_eqtl_allpossible.RDS")
m = readRDS("02_mashr/m_mixture_component.intra_eqtl_allpossible.RDS")


m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE, verbose = TRUE)

m2 %>% saveRDS("02_mashr/m2_intra_eqtl_results.intra_eqtl_allpossible.RDS")
m2 = readRDS("02_mashr/m2_intra_eqtl_results.intra_eqtl_allpossible.RDS")


## create map_cis

m2 = readRDS("02_mashr/m2_intra_eqtl_results.intra_eqtl_allpossible.RDS")
mashr_significant = m2 %>% 
  get_lfsr %>% as.data.frame() %>% rownames_to_column() %>%
  pivot_longer(cols = c(2:ncol(.))) %>% 
  rename(lfsr = value, name = rowname, condition= name)

mashr_significant = m2$result$PosteriorMean %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  pivot_longer(cols = c(2:ncol(.))) %>%
  dplyr::rename(PosteriorMean = value, condition = name, name = rowname) %>%
  left_join(mashr_significant, .)

mashr_significant = m2$result$PosteriorSD %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  pivot_longer(cols = c(2:ncol(.))) %>%
  dplyr::rename(PosteriorSD = value, condition = name, name = rowname) %>%
  left_join(mashr_significant, .)

mashr_significant = mashr_significant %>% 
  separate(name, into = c('phenotype_id', 'variant_id'), sep = "_")

process_parquet = function(prefix_file, df_mashr_statistics){
  require(arrow)
  df_mashr_statistics = df_mashr_statistics %>% filter(condition == prefix_file)
  c_parquets = paste0('01_calling_eQTL/01_02_tensorqtl_nominal/', prefix_file, '.cis_qtl_pairs.chr', c(1:22),'.parquet')
  df_parquet = lapply(c_parquets, read_parquet) %>% bind_rows() %>%
    filter(phenotype_id %in% df_mashr_statistics$phenotype_id,
           variant_id %in% df_mashr_statistics$variant_id)
  df_annotated = left_join(df_mashr_statistics, df_parquet) %>% 
    mutate(prefix_file = prefix_file)
  return(df_annotated)
}

library(doFuture)
registerDoFuture()
plan(multisession, workers = 16)

list_all = foreach(i = c(1:length(df_meta_target$prefix_file))) %dopar% {
  res = process_parquet(df_meta_target$prefix_file[i], mashr_significant)
}

registerDoSEQ()
plan(sequential)
df_all = list_all %>% bind_rows() %>% distinct()

df_map_cis_variants = df_all %>% na.omit() %>% 
  group_by(phenotype_id) %>%
  mutate(min_pval_nominal = min(pval_nominal)) %>%
  filter(pval_nominal == min_pval_nominal) %>%
  select(-min_pval_nominal) %>%
  top_n(n = 1, -lfsr) %>%
  add_count() %>%
  mutate(name = paste0(phenotype_id, "_", variant_id)) %>%
  ungroup() %>%
  distinct()

df_map_cis = df_all %>% na.omit() %>%  mutate(name = paste0(phenotype_id, "_", variant_id)) %>%
  filter(name %in% df_map_cis_variants$name)  %>%
  distinct()

df_all %>% write_delim("02_mashr/joined_all_results_unfiltered.txt.gz", delim = "\t")
df_map_cis %>% write_delim("02_mashr/joined_all_results.txt.gz", delim = "\t")