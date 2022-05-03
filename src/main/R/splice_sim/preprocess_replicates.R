#!/usr/bin/env Rscript
require(data.table)
require(tidyr)
require(dplyr)
require(rjson)
require(tictoc)
require(readr)

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
if (length(args)!=4) {
  stop("usage: preprocess_replicates.R <r1> <r2> <r3> <outdir", call.=FALSE)
} 
splice_sim_config=list()
splice_sim_config[['rep1']]=args[1]
splice_sim_config[['rep2']]=args[2]
splice_sim_config[['rep3']]=args[3]
outdir=args[4]

home_dir=list()
for ( r in names(splice_sim_config)) {
  home_dir[[r]]=paste0(dirname(splice_sim_config[[r]]),'/')
}

conf=fromJSON(paste(readLines(splice_sim_config[['rep1']]), collapse=""))
ref_base=conf$condition$ref
readlen=conf$readlen

# create results dir?
if (!dir.exists(outdir)) {
  dir.create(outdir)
} 

tic("preprocess metadata tables") 

m=list()
m[['tx']]=load_table(paste0(home_dir[['rep1']], 'eva/meta/', conf$dataset_name, '.tx.metadata.tsv.gz'))
m[['fx']]=load_table(paste0(home_dir[['rep1']], 'eva/meta/', conf$dataset_name, '.fx.metadata.tsv.gz'))
m[['sj']]=load_table(paste0(home_dir[['rep1']], 'eva/meta/', conf$dataset_name, '.sj.metadata.tsv.gz'))
m[['ga']]=load_table(paste0(home_dir[['rep1']], 'sim/reference_model/gene_anno.tsv.gz')) 
if (file.exists(paste0(home_dir_r1, 'tids.tsv'))) {
  m[['all_tids']]=load_table(paste0(home_dir[['rep1']], 'tids.tsv')) %>% pull(transcript_id)  
} else {
  m[['all_tids']]=NA
}
toc()

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
  left_join(m[['tx']] %>% select(tid, tx_rnk=rnk, num_exons, tx_mappability=mappability), by='tid') %>% # get cols from tx data
  mutate(chromosome=factor(chromosome),
         len=end-start+1,
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
  left_join(m[['tx']] %>% select(tid, tx_rnk=rnk, num_exons, tx_mappability=mappability), by='tid') %>% # get cols from tx data
  mutate(chromosome=factor(chromosome),
         len=end-start+1,
         acc_mappability = case_when(acc_win_map<0.2 ~ 'low',
                                     acc_win_map>0.9 ~ 'high',
                                     TRUE ~ 'medium'),
         don_mappability = case_when(don_win_map<0.2 ~ 'low',
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
    don_mappability=factor(don_mappability, levels=c('high', 'medium', 'low')),
    acc_mappability=factor(acc_mappability, levels=c('high', 'medium', 'low')),
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

meta_file=paste0(outdir,'/meta.rds')
saveRDS(m, meta_file, compress = FALSE)

toc()

# ===================================================================
# counts
# ===================================================================
tic("preprocess count tables") 
for (mq in c('','.mq20')) {
  all_data=tibble()
  for ( r in names(splice_sim_config)) {
    for ( m in names(conf$mappers) ) {
      for (cr in c(0, conf$condition$conversion_rates)) {
        tab = load_table(paste0(home_dir[[r]], 'eva/counts/',conf$dataset_name,'.cr',cr,'.',m,'.counts',mq,'.tsv.gz')) %>% 
          mutate(rep=r) %>% 
          group_by(rep, mapper, conversion_rate, class_type, fid, true_isoform, classification) %>% 
          summarise(all=sum(count), 
                    tc1=sum(count[cv1==1]), 
                    tc2=sum(count[cv2==1])) %>% 
          pivot_wider(names_from=classification, values_from=c(all, tc1, tc2), names_sort=T) %>% 
          mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>% 
          mutate(
            all_true =  all_TP+all_FN,
            all_found = all_TP+all_FP,
            tc1_true =  tc1_TP+tc1_FN,
            tc1_found = tc1_TP+tc1_FP,
            tc2_true =  tc2_TP+tc2_FN,
            tc2_found = tc2_TP+tc2_FP,
            true_fcr1 = ifelse(all_true>0,tc1_true/all_true,0),
            found_fcr1 = ifelse(all_found>0,tc1_found/all_found,0),
            true_fcr2 = ifelse(all_true>0,tc2_true/all_true,0),
            found_fcr2 = ifelse(all_found>0,tc2_found/all_found,0),
            # performance
            all_precision=ifelse(all_TP+all_FP>0, all_TP/(all_TP+all_FP), NA),
            all_recall=ifelse(all_TP+all_FN>0,all_TP/(all_TP+all_FN), NA),
            all_F1=ifelse((2*all_TP+all_FP+all_FN)>0,2*all_TP/(2*all_TP+all_FP+all_FN),NA),
            tc1_precision=ifelse(tc1_TP+tc1_FP>0, tc1_TP/(tc1_TP+tc1_FP), NA),
            tc1_recall=ifelse(tc1_TP+tc1_FN>0,tc1_TP/(tc1_TP+tc1_FN), NA),
            tc1_F1=ifelse((2*tc1_TP+tc1_FP+tc1_FN)>0,2*tc1_TP/(2*tc1_TP+tc1_FP+tc1_FN),NA),
            tc2_precision=ifelse(tc2_TP+tc2_FP>0, tc2_TP/(tc2_TP+tc2_FP), NA),
            tc2_recall=ifelse(tc2_TP+tc2_FN>0,tc2_TP/(tc2_TP+tc2_FN), NA),
            tc2_F1=ifelse((2*tc2_TP+tc2_FP+tc2_FN)>0,2*tc2_TP/(2*tc2_TP+tc2_FP+tc2_FN),NA)
          all_data = all_data %>% bind_rows(tab_fx)
      }
    } \
    data_file=paste0(outdir,'/data',mq,'.rds')
    saveRDS(d, data_file, compress = FALSE)
  }
}
toc()

