vlib = c("tidyverse", "mashr", "ashr", "ggpubr", "data.table", 
         "future.apply", "arrow", "optparse")
lapply(vlib, require, character.only = TRUE, quietly = TRUE) |> suppressMessages()

# setwd module. 2024.

base.dir = ""

setwd(paste0(base.dir, "202404-sceQTLv7/"))

option_list = list(
  make_option(c("-c", "--cluster"), action="store", default="Mono", type='character',
              help="select a cluster for analysis.")
  
)

opt = parse_args(OptionParser(option_list=option_list))


fname_meta = paste0(base.dir,'codebase-sceQTLv7/file_meta.txt')
df_meta = read_delim(fname_meta, col_names = c("levels", "time", "cluster_name")) %>% 
  mutate(prefix_file = paste0(levels, "_", time, "_", cluster_name))
df_meta$cluster_name %>% unique()

chr_cluster_name = opt$cluster

df_meta_target = df_meta %>% filter(levels == "anno_l1", cluster_name == chr_cluster_name)

df_meta_significant = df_meta_target %>%
  mutate(fname = paste0("01_calling_eQTL/01_01_tensorqtl_out/", prefix_file,".map_cis.txt.gz"))

df_meta_target = expand.grid(
  prefix_file = df_meta_target$prefix_file, 
  chr = paste0("chr", c(1:22))) %>%
  left_join(df_meta_target, .) %>%
  mutate(fname = paste0("01_calling_eQTL/01_02_tensorqtl_nominal/", prefix_file, ".cis_qtl_pairs.",chr,".parquet"))

list_meta_target = df_meta_target %>% group_by(prefix_file) %>% group_split()

library(doFuture)
registerDoFuture()
plan(multisession, workers = 8)

list_parquet = foreach(i = c(1:length(list_meta_target))) %dopar% {
  res = lapply(list_meta_target[[i]]$fname, read_parquet) %>%
    bind_rows() %>%
    select(phenotype_id, variant_id, slope, slope_se) %>%
    mutate(name = paste0(variant_id,"_",phenotype_id) ,
           condition =  unique(df_meta_target$prefix_file)[i]) %>%
    rename(SNPName = variant_id, gene_symbol = phenotype_id, 
           beta = slope, beta_se = slope_se) %>%
    select(name, SNPName, gene_symbol, beta, beta_se)
}

names(list_parquet) = df_meta_significant$prefix_file
registerDoSEQ()
plan(sequential)
df_parquet = list_parquet %>% bind_rows(.id  = "cluster")

# create df_significant, using names

df_all = read_delim("02_mashr/joined_all_results_unfiltered.txt.gz", delim = "\t")

# needs all possible, significant calls
df_significant_from_IO = df_all %>% 
  filter(lfsr < 0.05) %>%
  left_join(df_meta, by = c( "condition" = "prefix_file")) %>%
  select(phenotype_id, variant_id, cluster_name, slope, slope_se) %>%
  filter(cluster_name == chr_cluster_name) %>%
  mutate(name = paste0(variant_id,"_",phenotype_id)) %>%
  distinct(name)

# merging mashr_stats

list_mashr = list.files(
  "02_mashr/sceQTLgen-1MbloodNL",
  "mashr_stats.txt",
  full.names = T)

list_mashr = lapply(list_mashr, read_delim)
names(list_mashr) = list.files(
  "02_mashr/sceQTLgen-1MbloodNL",
  "mashr_stats.txt",
  full.names = F) %>% str_sub(., 1, -20)
list_mashr = c(list_mashr, list_parquet)

# list_strong
list_mashr_strong = list.files(
  "02_mashr/sceQTLgen-1MbloodNL",
  "mashr_stats.strong.txt",
  full.names = T)

list_mashr_strong = lapply(list_mashr_strong, read_delim)


df_mashr_strong = list_mashr_strong %>% bind_rows() %>% distinct(name) %>%
  rbind(df_significant_from_IO) %>% distinct(name) 


BHAT= vector(mode = "list",length =  length(list_mashr))
names(BHAT) = names(list_mashr)
SHAT= vector(mode = "list",length =  length(list_mashr))
names(SHAT) = names(list_mashr)
for (each in c(1:length(list_mashr))){
  BHAT[[each]] = list_mashr[[each]] %>% 
    select(name, beta) 
  SHAT[[each]] = list_mashr[[each]] %>% 
    select(name, beta_se) 
}

set.seed(42)
df_random = read_delim("02_mashr/sceQTLgen-1MbloodNL/240222_random_mashr_sets.names.txt.gz", delim = "\t")

mashr_data_random = list()


mashr_data_random$BHAT = BHAT %>% bind_rows(.id = 'cluster') %>% 
  filter(name %in% df_random$name) %>%
  pivot_wider(names_from = cluster, values_from= beta, values_fill = 0) %>%
  as.data.frame()%>%
  column_to_rownames(var = "name") %>%
  as.matrix()

mashr_data_random$SHAT = SHAT %>% bind_rows(.id = 'cluster') %>% 
  filter(name %in% df_random$name) %>%
  pivot_wider(names_from = cluster, values_from= beta_se, values_fill = 100) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()


mashr_data_strong = list()

mashr_data_strong$BHAT = BHAT %>% bind_rows(.id = 'cluster') %>% 
  filter(name %in% df_mashr_strong$name) %>%
  pivot_wider(names_from = cluster, values_from= beta, values_fill = 0) %>%
  as.data.frame()%>%
  column_to_rownames(var = "name") %>%
  as.matrix()


mashr_data_strong$SHAT = SHAT %>% bind_rows(.id = 'cluster') %>% 
  filter(name %in% df_mashr_strong$name) %>%
  pivot_wider(names_from = cluster, values_from= beta_se, values_fill =100) %>%
  as.data.frame() %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
 
mashr_data_temp = mash_set_data(mashr_data_random$BHAT,mashr_data_random$SHAT)
Vhat = estimate_null_correlation_simple(mashr_data_temp)
data.random = mash_set_data(mashr_data_random$BHAT,mashr_data_random$SHAT,V=Vhat)
data.strong = mash_set_data(mashr_data_strong$BHAT,mashr_data_strong$SHAT, V=Vhat)


U.pca = cov_pca(data.strong, npc = 8)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
list_Ulist = c(U.ed,U.c) 
list_Ulist %>% saveRDS(paste0("02_mashr/sceQTLgen-1MbloodNL/ulist.merged.sceQTLgen.",chr_cluster_name,".RDS"))

list_Ulist = readRDS(paste0("02_mashr/sceQTLgen-1MbloodNL/ulist.merged.sceQTLgen.",chr_cluster_name,".RDS"))
m = mash(data.random, Ulist = list_Ulist, outputlevel = 1, verbose = TRUE)

m %>% saveRDS(paste0("02_mashr/sceQTLgen-1MbloodNL/m_mixture_component.merged.sceQTLgen.",chr_cluster_name,".RDS"))

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE, verbose = TRUE)

m2 %>% saveRDS(paste0("02_mashr/sceQTLgen-1MbloodNL/m2_results.merged.sceQTLgen.",chr_cluster_name,".RDS"))
