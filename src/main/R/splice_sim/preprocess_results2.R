#!/usr/bin/env Rscript
require(data.table)
require(tidyr)
require(dplyr)
require(ggplot2)
require(scales)
require(rjson)
require(stringr)
require(VGAM)
require(cowplot)
require(arrow)
require(tictoc)
require(ggpubr)
require(minpack.lm)
require(readr)
require(testthat)
require(tidylog)
require(RColorBrewer)

# to ensure stripping '\0' (nul) from character vector
options(arrow.skip_nul = TRUE)

# load data from TSV
load_table = function(dataF, append="", header=T, nrows=Inf) {
  dataF = as.character(paste0(dataF, append))
  print(paste("Loading", dataF))
  if ( endsWith(dataF, ".gz") ) {
    return(fread(cmd=paste('gunzip -c', dataF), header=header, sep="\t", na.strings=c("na","NA",".", "None"), nrows=nrows))
  } else {
    return(fread(dataF, header=header, sep="\t", na.strings=c("na","NA",".", "None"), nrows=nrows))
  }
}

# args
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("usage: preprocess_results.R <splice_sim_config>", call.=FALSE)
} 
splice_sim_config=args[1]
home_dir=paste0(dirname(splice_sim_config),'/')
conf=fromJSON(paste(readLines(splice_sim_config), collapse=""))
ref_base=conf$condition$ref
readlen=conf$readlen

# create results dir?
if (!dir.exists(paste0(home_dir,'/results/')){
  dir.create(paste0(home_dir,'/results/')
} 
# results files
data_file=paste0(home_dir,'/results/data.rds')
meta_file=paste0(home_dir,'/results/data.rds')
  
tic("load evaluation data") # ca. 5-10min.

all_data=tibble()
for ( m in names(conf$mappers) ) {
  for (cr in c(0, conf$condition$conversion_rates)) {
    all_data = all_data %>% bind_rows(
      load_table(paste0(home_dir, 'eva/counts/',conf$dataset_name,'.cr',cr,'.',m,'.counts.tsv.gz'))
    )
  }
}

m=list()
m[['tx']]=load_table(paste0(home_dir, 'eva/meta/big3_slamseq_nf.tx.metadata.tsv.gz'))
m[['fx']]=load_table(paste0(home_dir, 'eva/meta/big3_slamseq_nf.fx.metadata.tsv.gz'))
m[['sj']]=load_table(paste0(home_dir, 'eva/meta/big3_slamseq_nf.sj.metadata.tsv.gz'))
m[['ga']]=load_table(paste0(home_dir, 'sim/reference_model/gene_anno.tsv.gz')) 
toc()

tic("data cleaning")  # data cleaning 

# ===================================================================
# tx metadata
# ===================================================================
m[['tx']] = m[['tx']] %>% 
  mutate(chromosome=factor(chromosome),
         len=end-start+1,
         num_exons=case_when(rnk<=5 ~ as.character(rnk), 
                             rnk>5  ~ '>5',
                             TRUE   ~ NA_character_),
         mappability=case_when(mean_map<0.2 ~ 'low',
                               mean_map>0.9 ~ 'high',
                               TRUE ~ 'medium'),
         GC=(G+C)/len,
         frac_convertible := !!rlang::sym(ref_base)/(A+C+T+G), # note: these are RNA bases so already corrected for strand
         convertibility=case_when(frac_convertible>0.3 ~ 'high',
                                  frac_convertible<0.2 ~ 'low',
                                  TRUE ~ 'medium'),
         
  ) %>% mutate(
    num_exons=factor(num_exons, levels=c('1','2','3','4','5','>5')),
    mappability=factor(mappability, levels=c('high', 'medium', 'low')),
    convertibility=factor(convertibility, levels=c('high', 'medium', 'low'))
  ) %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup()

# ===================================================================
# fx metadata
# ===================================================================

m[['fx']] = m[['fx']] %>% 
  mutate(chromosome=factor(chromosome),
         len=end-start+1,
         num_exons=case_when(rnk<=5 ~ as.character(rnk), 
                             rnk>5  ~ '>5',
                             TRUE   ~ NA_character_),
         mappability=case_when(mean_map<0.2 ~ 'low',
                               mean_map>0.9 ~ 'high',
                               TRUE ~ 'medium'),
         GC=(G+C)/len,
         frac_convertible := !!rlang::sym(ref_base)/(A+C+T+G), # note: these are RNA bases so already corrected for strand
         convertibility=case_when(frac_convertible>0.3 ~ 'high',
                                  frac_convertible<0.2 ~ 'low',
                                  TRUE ~ 'medium')
         
  ) %>% mutate(
    num_exons=factor(num_exons, levels=c('1','2','3','4','5','>5')),
    mappability=factor(mappability, levels=c('high', 'medium', 'low')),
    convertibility=factor(convertibility, levels=c('high', 'medium', 'low'))
  ) %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup()

# ===================================================================
# sj metadata
# ===================================================================

m[['sj']] = m[['sj']] %>% 
  mutate(chromosome=factor(chromosome),
         len=end-start+1,
         num_exons=case_when(rnk<=5 ~ as.character(rnk), 
                             rnk>5  ~ '>5',
                             TRUE   ~ NA_character_),
         mappability_acc = case_when(acc_win_map<0.2 ~ 'low',
                                     acc_win_map>0.9 ~ 'high',
                                     TRUE ~ 'medium'),
         mappability_don = case_when(don_win_map<0.2 ~ 'low',
                                     don_win_map>0.9 ~ 'high',
                                     TRUE ~ 'medium'),
         
         don_ex_fc := !!rlang::sym(paste0('don_ex_', ref_base))/(don_ex_A+don_ex_C+don_ex_T+don_ex_G), # note: these are RNA bases so already corrected for strand
         don_ex_fc=case_when(don_ex_fc>0.3 ~ 'high',
                             don_ex_fc<0.2 ~ 'low',
                              TRUE ~ 'medium'),
         don_in_fc := !!rlang::sym(paste0('don_in_', ref_base))/(don_in_A+don_in_C+don_in_T+don_in_G), # note: these are RNA bases so already corrected for strand
         don_in_fc=case_when(don_in_fc>0.3 ~ 'high',
                             don_in_fc<0.2 ~ 'low',
                             TRUE ~ 'medium'),
         acc_ex_fc := !!rlang::sym(paste0('acc_ex_', ref_base))/(acc_ex_A+acc_ex_C+acc_ex_T+acc_ex_G), # note: these are RNA bases so already corrected for strand
         acc_ex_fc=case_when(acc_ex_fc>0.3 ~ 'high',
                             acc_ex_fc<0.2 ~ 'low',
                             TRUE ~ 'medium'),
         acc_in_fc := !!rlang::sym(paste0('acc_in_', ref_base))/(acc_in_A+acc_in_C+acc_in_T+acc_in_G), # note: these are RNA bases so already corrected for strand
         acc_in_fc=case_when(acc_in_fc>0.3 ~ 'high',
                             acc_in_fc<0.2 ~ 'low',
                             TRUE ~ 'medium')
  ) %>% mutate(
    num_exons=factor(num_exons, levels=c('1','2','3','4','5','>5')),
    mappability_don=factor(mappability_don, levels=c('high', 'medium', 'low')),
    mappability_acc=factor(mappability_acc, levels=c('high', 'medium', 'low')),
    don_ex_fc=factor(don_ex_fc, levels=c('high', 'medium', 'low')),
    don_in_fc=factor(don_in_fc, levels=c('high', 'medium', 'low')),
    acc_ex_fc=factor(acc_ex_fc, levels=c('high', 'medium', 'low')),
    acc_in_fc=factor(acc_in_fc, levels=c('high', 'medium', 'low'))
  ) %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% 
    ungroup()
  
# ===================================================================
# ga metadata, see https://www.gencodegenes.org/pages/data_format.html
# ===================================================================

m[['ga']] = m[['ga']] %>% 
  select(tid=transcript_id, gene_type, gene_name, level, tag) %>% 
  mutate(gencode_level=case_when(
    level==1 ~ 'verified',
    level==2 ~ 'manually_annotated',
    level==3 ~ 'automatically_annotated',
    TRUE ~ NA_character_
  )) %>% 
  mutate(gencode_level=factor(gencode_level),
         gene_type=factor(gene_type)) %>% 
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) 

# ===================================================================
# counts
# ===================================================================

all_data = all_data %>% 
  select(-mq_fil) %>% 
  mutate(mapper=factor(mapper),
         conversion_rate=factor(conversion_rate),
         classification=factor(classification),
         class_type=factor(class_type)
  )

d=list()
d[['tx']] = all_data %>% 
  filter(class_type=='tx') %>% 
  select(-class_type) %>% 
  left_join(m[['tx']], by=c('fid'='tid')) %>% # NB too-short tx are lost
  left_join(m[['ga']], by=c('fid'='tid')) 


d[['fx']] = all_data %>% 
  filter(class_type=='fx') %>% 
  select(-class_type) %>% 
  left_join(m[['fx']], by='fid') %>% # NB too-short tx are lost
  left_join(m[['ga']], by='tid')


d[['sj']] = all_data %>% 
  filter(class_type %in% c('spl', 'don', 'acc')) %>% 
  left_join(m[['sj']], by='fid') %>% # NB too-short tx are lost
  left_join(m[['ga']], by='tid')

toc()

# ===================================================================
# Save data
# ===================================================================
tic("save data") # about 2min
saveRDS(d, data_file, compress = FALSE)
saveRDS(m, meta_file, compress = FALSE)
toc()

print(paste0("Done. Load data via d=readRDS('",data_file,"')"))
