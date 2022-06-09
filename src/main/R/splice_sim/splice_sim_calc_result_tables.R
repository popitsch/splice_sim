#!/usr/bin/env Rscript
# splice_sim_calc_result_tables.R
# author: niko.popitsch@univie.ac.at
require(data.table)
require(tidyr)
require(dplyr)
require(rjson)
require(tictoc)
require(readr)
require(writexl)

# ARGS-------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3 | length(args)>4) {
  stop("usage: splice_sim_calc_result_tables.R <splice_sim_config> <data.rds> <meta.rds> [<out_dir>]", call.=FALSE)
} 
splice_sim_config=args[1]
conf=fromJSON(paste(readLines(splice_sim_config), collapse=""))
data_file=args[2]
meta_file=args[3]
home_dir=paste0(dirname(splice_sim_config),'/')
if (length(args)<4) {
  outdir=paste0(home_dir,'/results/')
} else {
  outdir=args[4]
}
out_file_rds=paste0(outdir, conf$dataset_name, '.mappability.rds')
out_file_xlsx=paste0(outdir, conf$dataset_name, '.mappability.xlsx')
# __________==========================================================
# DATA ---------------------------------------------------------------
d=readRDS(data_file)
m=readRDS(meta_file)

tx = d[['tx']] %>% left_join(m[['tx']], by=c('fid'='tid')) 
fx = d[['fx']] %>% left_join(m[['fx']], by='fid') 
sj = d[['sj']] %>% 
  group_by(mapper, conversion_rate, class_type, fid, true_isoform, classification) %>% 
  summarise(count=sum(count)) %>% 
  ungroup()

results=list()

# __________===============================================================
# MAP ---------------------------------------------------------------

calc_map_table = function(tab) {
  tab = tab %>% 
    select(fid, ftype, mapper, conversion_rate, classification, count) %>% 
    group_by(fid, ftype, mapper, conversion_rate, classification) %>% 
    summarise(count=sum(count), .groups='drop') %>%
    pivot_wider(names_from=classification, values_from=count) %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>% 
    mutate(F1=ifelse(((2*TP+FP+FN)>0) & (TP+FN>0),2*TP/(2*TP+FP+FN),NA))
  
  # calc F1 from conversion_rate==0 only
  noconv=tab %>% filter(conversion_rate==0) %>% select(fid, mapper, ftype, F1_0=F1)
  # create final table
  tab = tab %>% group_by(fid,ftype,mapper) %>% 
    summarise(TP=sum(TP), FP=sum(FP), FN=sum(FN), .groups='drop') %>% 
    mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>% 
    mutate(F1=ifelse(((2*TP+FP+FN)>0) & (TP+FN>0),2*TP/(2*TP+FP+FN),NA)) %>% # set F1 to NA if no reads were simulated for tx (e.g., because too short)
    left_join(noconv, by=c('fid', 'mapper', 'ftype')) %>% # add F1_0
    pivot_wider(names_from=mapper, values_from=c(F1_0, F1), id_cols=c('fid', 'ftype'))
  
  return(tab)
}


# result tables
results[['map_tx']] = tx %>% calc_map_table() %>% rename(tid='fid') %>% 
  mutate(best_mapper_no_conv = ifelse(F1_0_HISAT3N-F1_0_STAR>0.05, 'HISAT3N',
                                      ifelse(F1_0_HISAT3N-F1_0_STAR< -0.05,'STAR','Both')),
         best_mapper_conv = ifelse(F1_HISAT3N-F1_STAR>0.05, 'HISAT3N',
                                   ifelse(F1_HISAT3N-F1_STAR< -0.05,'STAR','Both'))) %>% 
  left_join(m[['tx']], by=c('tid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid')) %>% 
  select(tid, gene_name, everything())

results[['map_ex']] = fx %>% filter(ftype=='exon') %>% calc_map_table() %>% 
  mutate(best_mapper_no_conv = ifelse(F1_0_HISAT3N-F1_0_STAR>0.05, 'HISAT3N',
                                      ifelse(F1_0_HISAT3N-F1_0_STAR< -0.05,'STAR','Both')),
         best_mapper_conv = ifelse(F1_HISAT3N-F1_STAR>0.05, 'HISAT3N',
                                   ifelse(F1_HISAT3N-F1_STAR< -0.05,'STAR','Both'))) %>% 
  left_join(m[['fx']], by=c('fid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid')) %>% 
  select(tid, gene_name, everything())

results[['map_in']] = fx %>% filter(ftype=='intron') %>% calc_map_table() %>% 
  mutate(best_mapper_no_conv = ifelse(F1_0_HISAT3N-F1_0_STAR>0.05, 'HISAT3N',
                                      ifelse(F1_0_HISAT3N-F1_0_STAR< -0.05,'STAR','Both')),
         best_mapper_conv = ifelse(F1_HISAT3N-F1_STAR>0.05, 'HISAT3N',
                                   ifelse(F1_HISAT3N-F1_STAR< -0.05,'STAR','Both'))) %>% 
  left_join(m[['fx']], by=c('fid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid')) %>% 
  select(tid, gene_name, everything())


# __________===============================================================
# FCR ---------------------------------------------------------------
fcr = 
  tx %>% 
    bind_rows(fx) %>% 
    group_by(fid, mapper, conversion_rate, classification, true_isoform, mappability, ftype) %>% 
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
    ) %>% 
    mutate(ftype=factor(ftype, levels=c('tx', 'exon', 'intron')))

fcr_perf = fcr %>%
    select(fid, conversion_rate,true_isoform,ftype,mappability,true_fcr1, found_fcr1) %>% 
    ungroup() %>% 
    pivot_wider(names_from=mapper, values_from=c(true_fcr1, found_fcr1)) %>% 
    mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>%   # replace NA with 0; NA are resulting from mat not having introns. So if there are no FP values then there are no entries!
    rename(simulated=true_fcr1_HISAT3N, HISAT3N=found_fcr1_HISAT3N, STAR=found_fcr1_STAR) %>%
    select(-true_fcr1_STAR) %>% 
    mutate(diff_HISAT3N=HISAT3N-simulated, 
           diff_STAR=STAR-simulated,
           best_diff=pmin(abs(diff_HISAT3N), abs(diff_STAR)),
           best_mapper=case_when(
             best_diff>0.1 ~ 'None',
             abs(diff_HISAT3N)<0.05 & abs(diff_STAR)<0.05 ~ 'Both',
             abs(diff_HISAT3N)<abs(diff_STAR) ~ 'HISAT3N',
             T ~ 'STAR')
    )

# result tables
fcr_tab = fcr_perf %>% 
  filter(as.numeric(as.character(conversion_rate))>0) %>% 
  group_by(fid, ftype) %>% 
  summarise(m_diff_HISAT3N = mean(diff_HISAT3N,na.rm=T), 
            m_diff_STAR = mean(diff_STAR,na.rm=T)) %>% 
  mutate(best_m_diff=pmin(abs(m_diff_HISAT3N), abs(m_diff_STAR)),
         best_mapper=case_when(
           best_m_diff>0.1 ~ 'None',
           abs(m_diff_HISAT3N)<0.05 & abs(m_diff_STAR)<0.05 ~ 'Both',
           abs(m_diff_HISAT3N)<abs(m_diff_STAR) ~ 'HISAT3N',
           T ~ 'STAR')
  )

results[['fcr_tx']] = fcr_tab %>% filter(ftype=='tx') %>%
  rename(tid=fid) %>% 
  left_join(m[['tx']], by=c('tid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid'))  %>% 
  select(tid, gene_name, everything())

results[['fcr_ex']] = fcr_tab %>% filter(ftype=='exon') %>% 
  left_join(m[['fx']], by=c('fid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid'))  %>% 
  select(tid, gene_name, everything())

results[['fcr_in']] = fcr_tab %>% filter(ftype=='intron') %>% 
  left_join(m[['fx']], by=c('fid', 'ftype')) %>% 
  left_join(m[['ga']], by=c('tid'))  %>% 
  select(tid, gene_name, everything())

# __________===============================================================
# FMAT ---------------------------------------------------------------

sj_per = 
  sj %>% 
    group_by(mapper, conversion_rate, fid, class_type, classification) %>% 
    summarise(count=sum(count)) %>% 
    pivot_wider(names_from=c(class_type, classification), values_from=count) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>% 
    left_join(m[['sj']], by='fid') %>% 
    mutate(
      spl_pre=ifelse(spl_TP+spl_FP>0, spl_TP/(spl_TP+spl_FP), NA),
      spl_rec=ifelse(spl_TP+spl_FN>0, spl_TP/(spl_TP+spl_FN), NA),
      spl_F1=ifelse((2*spl_TP+spl_FP+spl_FN)>0,2*spl_TP/(2*spl_TP+spl_FP+spl_FN),NA),
      acc_pre=ifelse(acc_TP+acc_FP>0, acc_TP/(acc_TP+acc_FP), NA),
      acc_rec=ifelse(acc_TP+acc_FN>0, acc_TP/(acc_TP+acc_FN), NA),
      acc_F1=ifelse((2*acc_TP+acc_FP+acc_FN)>0,2*acc_TP/(2*acc_TP+acc_FP+acc_FN),NA),
      don_pre=ifelse(don_TP+don_FP>0, don_TP/(don_TP+don_FP), NA),
      don_rec=ifelse(don_TP+don_FN>0, don_TP/(don_TP+don_FN), NA),
      don_F1=ifelse((2*don_TP+don_FP+don_FN)>0,2*don_TP/(2*don_TP+don_FP+don_FN),NA),
      true_mat=spl_FN+spl_TP,
      found_mat=spl_FP+spl_TP,
      true_pre=don_FN+don_TP+acc_FN+acc_TP,
      found_pre=don_FP+don_TP+acc_FP+acc_TP,
      true_spliced=ifelse(true_pre>0,true_mat/(true_pre+true_mat), NA),
      found_spliced=ifelse(found_pre>0,found_mat/(found_pre+found_mat), NA),
      spliced_diff=abs(true_spliced-found_spliced)
    )  %>% mutate(spliced_diff_rnk = dense_rank(spliced_diff)) # rank of difference to true spliced value. Bad introns have higher ranks

# calc intron mappability
min_spliced_diff = 0.1 # filter SJ only if diff to true value is greater than this value
min_remaining_readcount = 100 # filter SJs util this remaining_readcount
intron_mappability = 
  sj_per %>% 
    group_by(tid, mapper, conversion_rate) %>% 
    arrange(spliced_diff_rnk) %>% 
    mutate(read_sum=cumsum(found_mat+found_pre)) %>% 
    select(tid, mapper, conversion_rate, fid, spliced_diff_rnk, spliced_diff, 
           true_pre, true_mat, found_pre, found_mat, read_sum) %>% 
    arrange(desc(spliced_diff_rnk)) %>% # sort per tx by decreasing intron rnk, i.e.,bad introns first
    mutate(remaining_reads=lead(read_sum, default=0)) %>% # remaining reads after removing the current intron
    mutate(low_map = (is.na(spliced_diff) | 
                        (spliced_diff>!!min_spliced_diff & remaining_reads>!!min_remaining_readcount))) %>% 
    group_by(mapper, tid, fid) %>% summarise(
      include_converted=ifelse(sum(low_map)>sum(!low_map),1,0),   # filter because in most conditions it is filtered
      unconverted_only=low_map[conversion_rate==0]  # filter because in cond0 it is filtered
    ) %>% 
    ungroup() %>% 
    pivot_longer(c(include_converted, unconverted_only), names_to = 'filter_type', values_to='filtered')

  
sfrac_conv =
  sj_per %>% 
    left_join(intron_mappability, by=c('mapper', 'fid', 'tid' )) %>% 
    filter(filter_type=='include_converted') %>%  # TC part
    group_by(tid, mapper, conversion_rate, filter_type, tx_mappability) %>% 
    summarise(true_mat=sum(spl_FN)+sum(spl_TP),
              found_mat=sum(spl_FP)+sum(spl_TP),
              true_pre=sum(don_FN)+sum(don_TP)+sum(acc_FN)+sum(acc_TP),
              found_pre=sum(don_FP)+sum(don_TP)+sum(acc_FP)+sum(acc_TP),
              # --- filtered -
              true_mat_fil=sum(spl_FN[!filtered])+sum(spl_TP[!filtered]),
              found_mat_fil=sum(spl_FP[!filtered])+sum(spl_TP[!filtered]),
              true_pre_fil=sum(don_FN[!filtered])+sum(don_TP[!filtered])+sum(acc_FN[!filtered])+sum(acc_TP[!filtered]),
              found_pre_fil=sum(don_FP[!filtered])+sum(don_TP[!filtered])+sum(acc_FP[!filtered])+sum(acc_TP[!filtered]),
              n_fil=sum(filtered, na.rm = T),
              .groups = 'drop'
    )

sfrac_unconv = 
  sj_per %>% 
    left_join(intron_mappability, by=c('mapper', 'fid', 'tid' )) %>% 
    filter(filter_type=='unconverted_only', conversion_rate==0) %>%  # no-TC part
    group_by(tid, mapper, conversion_rate, filter_type, tx_mappability) %>% 
    summarise(true_mat=sum(spl_FN)+sum(spl_TP),
              found_mat=sum(spl_FP)+sum(spl_TP),
              true_pre=sum(don_FN)+sum(don_TP)+sum(acc_FN)+sum(acc_TP),
              found_pre=sum(don_FP)+sum(don_TP)+sum(acc_FP)+sum(acc_TP),
              # --- filtered -
              true_mat_fil=sum(spl_FN[!filtered])+sum(spl_TP[!filtered]),
              found_mat_fil=sum(spl_FP[!filtered])+sum(spl_TP[!filtered]),
              true_pre_fil=sum(don_FN[!filtered])+sum(don_TP[!filtered])+sum(acc_FN[!filtered])+sum(acc_TP[!filtered]),
              found_pre_fil=sum(don_FP[!filtered])+sum(don_TP[!filtered])+sum(acc_FP[!filtered])+sum(acc_TP[!filtered]),
              n_fil=sum(filtered, na.rm = T),
              .groups = 'drop'
    )

sfrac = 
  bind_rows(sfrac_conv, sfrac_unconv) %>% 
    filter( found_mat_fil+found_pre_fil>!!min_remaining_readcount ) %>% # drop if due to filtering the overall readcount became too low
    mutate(
      true_spliced=ifelse(true_pre>0,true_mat/(true_pre+true_mat), NA),
      found_spliced=ifelse(found_pre>0,found_mat/(found_pre+found_mat), NA),
      found_spliced_fil=ifelse(found_pre_fil>0,found_mat_fil/(found_pre_fil+found_mat_fil), NA),
    ) %>% 
    mutate(diff_spliced=true_spliced-found_spliced,
           diff_spliced_fil=true_spliced-found_spliced_fil) %>% 
    mutate(improvement=abs(diff_spliced) - abs(diff_spliced_fil)) %>% 
    left_join(m[['tx']] %>% select(tid, n_overlapping, rnk), by='tid') %>% 
    left_join(m[['ga']] %>% select(tid, gene_type), by='tid')

best_mapper_fmat = 
  sfrac %>% 
    select(mapper, tid, conversion_rate,filter_type, n_fil, n_overlapping, 
           true_spliced, found_spliced, found_spliced_fil, tx_mappability) %>% 
    group_by(mapper, tid, n_fil, n_overlapping, true_spliced, tx_mappability, filter_type) %>% 
    summarise( m_found_spliced=median(found_spliced, na.rm=T), 
               m_found_spliced_fil=median(found_spliced_fil)) %>% 
    ungroup() %>% 
    pivot_wider(names_from=mapper, 
                values_from=c(m_found_spliced, m_found_spliced_fil, n_fil), 
                id_cols=c('tid', 'n_overlapping', 'true_spliced', 'tx_mappability', 'filter_type')) %>% 
    mutate(diff_STAR=abs(true_spliced-m_found_spliced_STAR),
           diff_STAR_fil=abs(true_spliced-m_found_spliced_fil_STAR),
           diff_HISAT3N=abs(true_spliced-m_found_spliced_HISAT3N),
           diff_HISAT3N_fil=abs(true_spliced-m_found_spliced_fil_HISAT3N)
    ) %>%   
    mutate(best_HISAT3N=pmin(diff_HISAT3N, diff_HISAT3N_fil),
           best_STAR=pmin(diff_STAR, diff_STAR_fil )) %>% 
    mutate(best_diff=pmin(best_HISAT3N, best_STAR)) %>% 
    mutate(best_mapper_fmat=ifelse(best_HISAT3N<best_STAR, 'HISAT3N', ifelse(best_STAR<best_HISAT3N,'STAR','Both')),
           after_fil=(diff_STAR_fil==best_diff | diff_HISAT3N_fil==best_diff) & (diff_HISAT3N!=best_diff) & (diff_STAR!=best_diff),
           diff_to_other=abs(best_HISAT3N-best_STAR)
    ) 

results[['spl_tx_include_converted']] = best_mapper_fmat %>% 
  filter(filter_type=='include_converted') %>% 
  left_join(m[['tx']] %>% select(-n_overlapping), by='tid') %>% 
  left_join(m[['ga']], by='tid') %>% 
  select(tid, gene_name, everything())

results[['spl_tx_unconverted_only']] = best_mapper_fmat %>% 
  filter(filter_type=='unconverted_only') %>% 
  left_join(m[['tx']] %>% select(-n_overlapping), by='tid') %>% 
  left_join(m[['ga']], by='tid') %>% 
  select(tid, gene_name, everything())

# __________===============================================================
# WRITE -------------------------------------------------------------------
tab_summary =
  results[['map_tx']] %>% select(gene_name, tid, mappability,
                                     map_best_no_conv=best_mapper_no_conv, 
                                     map_best_conv=best_mapper_conv) %>% 
  full_join(
    results[['fcr_tx']] %>% select(gene_name, tid, 
                                       fcr_best=best_mapper), 
    by=c('gene_name', 'tid' )) %>% 
  full_join(
    results[['spl_tx_include_converted']] %>% select(gene_name, tid, 
                                       spl_best_include_converted=best_mapper_fmat), 
    by=c('gene_name', 'tid' )) %>% 
  full_join(
    results[['spl_tx_unconverted_only']] %>% select(gene_name, tid, 
                                                     spl_best_unconverted_only=best_mapper_fmat), 
    by=c('gene_name', 'tid' )) 

results=c(list('summary'=tab_summary), results)

tab_config = tribble(
  ~property, ~value,
  "dataset_name", conf$dataset_name,
  "configured_tx", as.character(length(m[['all_tids']] )),
  "config_file",  splice_sim_config,
  "out_file_rds", out_file_rds,
  "out_file_xlsx", out_file_xlsx
)
results=c(list('configuration'=tab_config), results)

# write rds
saveRDS(d, out_file_rds, compress = FALSE)

# write excel result file
write_xlsx(results, path=out_file_xlsx)

print(paste0("Done. Load data via mappability=readRDS(",out_file_rds,".rds)"))
