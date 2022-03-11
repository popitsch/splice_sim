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
if (length(args)==0) {
  stop("usage: preprocess_results.R <splice_sim_config>", call.=FALSE)
} 
splice_sim_config=args[1]

home_dir=paste0(dirname(splice_sim_config),'/')
conf=fromJSON(paste(readLines(splice_sim_config), collapse=""))
data_file=paste0(home_dir,'/results/data.rds')

tic("load evaluation data") # ca. 5-10min.

d=list()
d[['all_tids']]=load_table(paste0(home_dir, 'tids.tsv')) %>% pull(transcript_id)
d[['tid_performance']]=open_dataset(paste0(home_dir,'/results/overall_performance/'), format="parquet") %>% collect()
d[['splice_site_performance']]=open_dataset(paste0(home_dir,'/results/splice_site_performance/')) %>% collect()
d[['splice_site_mappability']]=open_dataset(paste0(home_dir,'/results/splice_site_mappability/')) %>% collect()
d[['gene_anno']]=open_dataset(paste0(home_dir,'/results/gene_anno/')) %>% collect()
fc_schema_mismapped=schema(
  condition_id = string(), # need to set explicit types as otherwise arrow thinks this is an int
  tid = string(),
  fid = string(),
  start = int64(),
  end = int64(),
  FP = int64(),
  FN = int64(),
  mapper = string(),
  chromosome = string(),
  ftype = string()
)
d[['feature_counts_mismapped']]=open_dataset(paste0(home_dir,'/results/feature_counts_mismapped/'), schema = fc_schema_mismapped) # too big, do not collect
d[['feature_counts']]=open_dataset(paste0(home_dir,'/results/feature_counts/')) # too big, do not collect
d[['tx_metadata']]=open_dataset(paste0(home_dir,'/results/tx_metadata/')) %>% collect()
d[['sj_metadata']]=open_dataset(paste0(home_dir,'/results/sj_metadata/')) %>% collect()
# optional: read info
if (file.exists(paste0(home_dir,'/results/reads/')) ) {
  d[['reads']]=open_dataset(paste0(home_dir,'/results/reads/')) # too big, do not collect
} else {
  d[['reads']]=NA
}
d[['extra_transcript_table']]=load_table(params$extra_transcript_table) %>% select(transcript_id, rnk, GC_frac, mappability_median, mappability_mean, T_Pos)
toc()

tic("data cleaning") 
# data cleaning ~1min
# add mappability, coords, etc.
d[['gene_anno']] = d[['gene_anno']] %>% 
  mutate(transcript_len=end-start+1) %>% 
  rename(tid=transcript_id) %>% 
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) 
d[['extra_transcript_table']] = d[['extra_transcript_table']] %>% 
  rename(tid=transcript_id) %>% 
  select(-rnk)
d[['tid_performance']]=d[['tid_performance']] %>% 
  merge(d[['gene_anno']], by='tid', all.x=T) %>% 
  merge(d[['extra_transcript_table']], by='tid', all.x=T) %>% 
  mutate(mapper=factor(mapper),
         condition_id=factor(condition_id),
         precision=TP/(TP+FP), 
         recall=TP/(TP+FN), 
         F1=2*TP/(2*TP+FP+FN)
  ) %>% 
  mutate(num_exons=case_when(rnk<=5 ~ as.character(rnk), 
                             rnk>5 ~ '>5',
                             T ~ 'NA'),
         mappability=case_when(mappability_mean<0.2 ~ 'low',
                               mappability_mean>0.9 ~ 'high',
                               T ~ 'medium')
  ) %>% mutate(
    num_exons=factor(num_exons, levels=c('1','2','3','4','5','>5')),
    mappability=factor(mappability, levels=c('high', 'medium', 'low'))
  ) %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) 

d[['sj_metadata']] = d[['sj_metadata']] %>% 
  mutate(across(where(is.character), ~na_if(., "None"))) %>% 
  mutate(acc_win_min_map=as.numeric(acc_win_min_map),
         acc_win_max_map=as.numeric(acc_win_max_map),
         don_win_min_map=as.numeric(don_win_min_map),
         don_win_max_map=as.numeric(don_win_max_map))

d[['splice_site_mappability']] = d[['splice_site_mappability']] %>% 
  left_join(d[['sj_metadata']], by=c( "tid","intron_id","chromosome","start","end","strand" )) %>% 
  mutate(
    acc_win_min_map=replace(acc_win_min_map, acc_win_min_map == "None", 0), # set mappability to zero if 'None' reported
    acc_win_max_map=replace(acc_win_max_map, acc_win_max_map == "None", 0),
    don_win_min_map=replace(don_win_min_map, don_win_min_map == "None", 0),
    don_win_max_map=replace(don_win_max_map, don_win_max_map == "None", 0),
  )  %>% 
  separate(intron_id, c(NA, "intron_rnk"), sep = '_') %>%
  mutate(intron_rnk=ifelse(strand=='-',as.numeric(intron_rnk),as.numeric(intron_rnk)-1)) %>% # NB there is a 'bug' when assigning intron ids: + strand introns start counting at 2...
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) # replace NaN with NA

# intron_coords=d[['feature_counts']] %>% 
#   filter(mapper=='HISAT3N', condition_id=='0.01', ftype=='intron') %>% 
#   select(fid, chromosome,start,end) %>% 
#   distinct() %>% 
#   collect() %>% 
#   mutate(intron_len=end-start+1) %>% 
#   rename(intron_id=fid)

d[['splice_site_performance']] = d[['splice_site_performance']] %>% 
  merge(d[['gene_anno']] %>% select(-c(chromosome,start,end,strand)), by='tid', all.x=T) %>% 
  merge(d[['extra_transcript_table']], by='tid', all.x=T) %>% 
  mutate(mapper=factor(mapper),
         condition_id=factor(condition),
         intron_len=end-start+1,
         spl_pre=spl_TP/(spl_TP+spl_FP), 
         spl_rec=spl_TP/(spl_TP+spl_FN), 
         spl_F1=2*spl_TP/(2*spl_TP+spl_FP+spl_FN),
         don_pre=don_TP/(don_TP+don_FP), 
         don_rec=don_TP/(don_TP+don_FN), 
         don_F1=2*don_TP/(2*don_TP+don_FP+don_FN),
         acc_pre=acc_TP/(acc_TP+acc_FP), 
         acc_rec=acc_TP/(acc_TP+acc_FN), 
         acc_F1=2*acc_TP/(2*acc_TP+acc_FP+acc_FN),
         all_F1=2*(spl_TP+don_TP+acc_TP)/(2*(spl_TP+don_TP+acc_TP)+(spl_FP+don_FP+acc_FP)+(spl_FN+don_FN+acc_FN)),
         true_mat=spl_FN+spl_TP,
         found_mat=spl_FP+spl_TP,
         true_pre=don_FN+don_TP+acc_FN+acc_TP,
         found_pre=don_FP+don_TP+acc_FP+acc_TP
  ) %>% 
  mutate(true_spliced=ifelse(true_pre>0,true_mat/(true_pre+true_mat), NA),
         found_spliced=ifelse(found_pre>0,found_mat/(found_pre+found_mat), NA),
         spliced_diff=abs(true_spliced-found_spliced)) %>% 
  mutate(spliced_diff_rnk = dense_rank(spliced_diff)) %>% # rank of difference to true spliced value. Bad introns have higher ranks
  separate(intron_id, c(NA, "intron_rnk"), sep = '_', remove=F) %>%
  mutate(intron_rnk=ifelse(strand=='-',as.numeric(intron_rnk),as.numeric(intron_rnk)-1)) %>% # NB there is a 'bug' when assigning intron ids: + strand introns start counting at 2...
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) # replace NaN with NA


# rename columns and match our id schema for exons/introns
d[['feature_counts']] = d[['feature_counts']] %>% 
  collect() %>% 
  rename(chromosome=Chr,
         start=Start, 
         end=End,
         length=Length,
         strand=Strand) %>% 
  mutate(fid=case_when(
    ftype=='exon' ~ str_replace(str_sub(Geneid, start=nchar('exon:')+1), ':', '_ex'),
    ftype=='intron' ~ str_replace(str_sub(Geneid, start=nchar('intron:')+1), ':', '_'),
    T ~ NA_character_
  )) %>% 
  separate(fid, c('tid', NA), '_', remove=F) %>% 
  ungroup()
# exon: ENSMUST00000193812.1_ex1
# intron: ENSMUST00000161581.1_2

d[['feature_counts_mismapped']] =
  d[['feature_counts_mismapped']] %>%
  collect()

# FIXME: some introns in d[['feature_counts_mismapped']] are misnamed (wrong intron #), so merge w/o fid here
d[['counts_per_feature']] = d[['feature_counts_mismapped']] %>% 
  left_join(d[['feature_counts']], by=c('tid', 'chromosome', 'start', 'end', 'ftype', 'mapper', 'condition_id' ) ) %>% 
  left_join(d[['gene_anno']] %>% select(tid, gene_name, gene_type, strand,rnk,transcript_len), 
            by=c('tid', 'strand')) %>% 
  left_join(d[['extra_transcript_table']], by='tid') %>% 
  distinct() %>% # remove duplicates due to joins
  mutate(mapper=factor(mapper), condition_id=factor(condition_id)) %>% 
  mutate(cov=read_count*conf$readlen/length) %>% 
  mutate(TP=read_count+FN-FP) %>% 
  mutate(F1=2*TP/(2*TP+FP+FN)) %>% 
  mutate(num_exons=case_when(rnk<=5 ~ as.character(rnk), 
                             rnk>5 ~ '>5',
                             T ~ 'NA'),
         mappability=case_when(mappability_mean<0.2 ~ 'low',
                               mappability_mean>0.9 ~ 'high',
                               T ~ 'medium')
  ) %>% mutate(
    num_exons=factor(num_exons, levels=c('1','2','3','4','5','>5')),
    mappability=factor(mappability, levels=c('high', 'medium', 'low'))
  ) %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup()

d[['counts_per_tx']] = d[['counts_per_feature']] %>% 
  group_by(mapper, condition_id, tid, ftype, strand, mappability, GC_frac, gene_name, gene_type) %>% 
  summarise(FP=sum(FP), FN=sum(FN), TP=sum(TP), read_count=sum(read_count), length=sum(length)) %>% 
  mutate(cov=read_count*conf$readlen/length) %>% 
  mutate(F1=2*TP/(2*TP+FP+FN)) %>% 
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% 
  ungroup()

toc()

tic("save data") # about 2min
saveRDS(d, data_file, compress = FALSE)
toc()

print("Done. Load data via d=readRDS('",data_file,"')")
