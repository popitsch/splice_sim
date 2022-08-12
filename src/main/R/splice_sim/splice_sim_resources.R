
# __________===============================================================
# Functions ---------------------------------------------------------------


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

# load RDS file
load_rds = function(data_file) {
  if (file.exists(data_file)) {
    tic(paste0("load ", data_file))
    d=readRDS(data_file)
    toc()
    return(d)
  } else {
    stop("Could not find data_file")
  }
}

# multiple plots with single title
my_plot_grid = function(plots, main='', ncol=NULL, nrow=NULL, labels=NULL, rel_widths = 1, rel_heights = 1) {

  if (!is.null(labels)){
    if (labels==T) { labels=LETTERS[1:length(plots)] }
  }
  plot_row=plot_grid(plotlist=plots, ncol=ncol, nrow=nrow, labels=labels, rel_heights=rel_heights, rel_widths=rel_widths)
  title <- ggdraw() + draw_label( main, fontface = 'bold', x = 0, hjust = 0 ) + theme(plot.margin = margin(0, 0, 0, 7))
  return (plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1)))
}

# calculate confidence interval (ci),  a measure of precision
# call with mtcars %>% calc_ci(mpg) or mtcars %>% group_by(cyl) %>% calc_ci(mpg)
# plot with ...  %>% ggplot(aes(x=cyl, y=col.mean)) + geom_line() + geom_ribbon(aes(ymin=col.lower, ymax=col.upper), alpha=0.1 ...)
calc_ci = function(d, col, min_lower_ci=NA, max_upper_ci=NA) {
  ret = d %>% summarise(
    col.mean := mean({{col}}, na.rm = TRUE),
    #col.median := median({{col}}, na.rm = TRUE),
    col.sd = sd({{col}}, na.rm = TRUE),
    col.n = n(),
    .groups = 'drop') %>%
    mutate(stderr = col.sd / sqrt(col.n),
           col.lower = col.mean - qt(1 - (0.05 / 2), col.n - 1) * stderr,
           col.upper = col.mean + qt(1 - (0.05 / 2), col.n - 1) * stderr) %>%
    ungroup() %>%
    mutate(col.lower = pmax(min_lower_ci, col.lower, na.rm=T),
           col.upper = pmin(max_upper_ci, col.upper, na.rm=T))
  return (ret)
}

# calculate interquartile range (iqr), a measure of dispersion
# call with mtcars %>% calc_iqr(mpg) or mtcars %>% group_by(cyl) %>% calc_iqr(mpg)
# plot with  ... %>% ggplot(aes(x=cyl, y=col.median)) + geom_line() + geom_ribbon(aes(ymin=col.lower, ymax=col.upper), alpha=0.1 ...)
calc_iqr = function(d, col, min_lower_ci=NA, max_upper_ci=NA) {
  d %>% summarise(
    col.median = median({{col}}, na.rm = TRUE),
    #col.mean = mean({{col}}, na.rm = TRUE),
    col.upper = quantile({{col}}, .75, na.rm = TRUE),
    col.lower = quantile({{col}}, .25, na.rm = TRUE),
    col.n = n(),
    .groups = 'drop') %>%
    mutate(col.lower = pmax(min_lower_ci, col.lower, na.rm=T),
           col.upper = pmin(max_upper_ci, col.upper, na.rm=T))
}

# calculate outlier cutoffs based on IQR
# Lower Bound: (Q1 - 1.5 * IQR)
# Upper Bound: (Q3 + 1.5 * IQR)
# usage: mtcars %>% group_by(cyl) %>% calc_outlier(mpg)
calc_outlier = function(d, col) {
  d %>% summarise(
    col.median = median({{col}}, na.rm = TRUE),
    col.iqr = quantile({{col}}, .75, na.rm = TRUE)-quantile({{col}}, .25, na.rm = TRUE),
    col.upper = quantile({{col}}, .75, na.rm = TRUE) + 1.5 * col.iqr,
    col.lower = quantile({{col}}, .25, na.rm = TRUE) - 1.5 * col.iqr,
    col.n = n(),
    .groups = 'drop')
}

#
# write result tibble to BED file
#
write_bed = function(dat, bed_file, title, header=F) {
  sink(bed_file)
  cat(paste0("track name=",title," description=\"",title,"\" useScore=1 itemRgb=\"On\"\n"))
  sink()
  dat %>% write_tsv( bed_file, col_names = F, append = T )
}

# calculates performance and coverage on grouped data.
# required columns: count, classification, len
calc_performance=function(tab, readlen) {
  tab %>%
    summarise(count=sum(count)) %>%
    pivot_wider(names_from=classification, values_from=count, names_sort=T) %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>%
    mutate(
      read_count=TP+FN,
      cov=read_count*!!readlen/len,
      precision=ifelse(TP+FP>0, TP/(TP+FP), NA),
      recall=ifelse(TP+FN>0,TP/(TP+FN), NA),
      F1=ifelse((2*TP+FP+FN)>0,2*TP/(2*TP+FP+FN),NA)
    )
}

# calculate coverage for a grouped table
# calc_coverage = function(grp_tab) {
#   grp_tab %>%
#     summarise(count=sum(count)) %>%
#     pivot_wider(names_from=c(mapper,classification), values_from=count, names_sort=T) %>%
#     mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), 0, .))) %>%
#     filter(len>conf$readlen) %>% # filter too-short introns as this will lead to wrong coverage calc. Example: ENSMUST00000116560.2_in5
#     mutate(
#       simulated=(HISAT3N_TP+HISAT3N_FN)*!!conf$readlen/len, # same as simulated_star!
#       HISAT3N=(HISAT3N_TP+HISAT3N_FP)*!!conf$readlen/len,
#       STAR=(STAR_TP+STAR_FP)*!!conf$readlen/len,
#     ) %>% pivot_longer(c(simulated, HISAT3N, STAR), names_to='mapper') %>%
#     mutate(mapper=factor(mapper, levels=c('simulated', 'HISAT3N', 'STAR')),
#            is_converted=ifelse(conversion_rate==0,'no conversions', 'converted reads')) %>%
#     ungroup()
# }

#see https://github.com/tidyverse/dplyr/issues/4223
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}

# simple exponential decay model
decay_model= function(t, k) {
  return( exp(t * -k) )
}
# fit decay model and calculate halflife
fit_halflife = function(TP, dat) {
  if (any(is.na(dat))) {
    return (list(mod=NA,k=NA,hl=NA,pseudoR2=NA, bic=NA))
  }
  mod = tryCatch({
    nlsLM(dat~decay_model(TP, k),
          start=list(
            k=0),
          lower = c(0),
          upper = c(Inf),
          control = nls.lm.control(maxiter = 1000),
          na.action = na.omit)
  }, error=function(e){
    print(e)
  })
  if ( inherits(mod, "simpleError")) {return (list(mod=NA,k=NA,hl=NA,pseudoR2=NA, bic=NA)) }
  k=coef(mod)[length(coef(mod))]
  hl=ifelse(k>0, log(2)/k, NA)
  rss = sum(residuals(mod)^2)
  tss = sum((dat - mean(dat ,na.rm = TRUE))^2 ,na.rm = TRUE)  # Total sum of squares
  pseudoR2 = 1 - (rss/tss)  # R-squared measure
  bic = BIC(mod)
  return (list(mod=mod,k=k,hl=hl,pseudoR2=pseudoR2, bic=bic))
}

calc_cor = function(d,a,b) {
  return(
    tibble(
      r_pearson=round(cor(d[[a]], d[[b]], use = "complete.obs"), 4),
      r_spearman = round(cor(d[[a]], d[[b]], use = "complete.obs", method="spearman"), 4),
      n =nrow(na.omit(d %>% select(all_of(a),all_of(b))))
    )
  )
}

# simple correlation plot
plot_corr = function(d, a, b, col_, shape_=NULL, main_title=NA, xlog=F, use_hex=F, draw_diag=T, 
                     max_x=NA, show_legend=T, alpha_=0.2, use_scattermore=FALSE, calc_corr=T) {
  thecor=''
  if ( calc_corr) {
    thecor = paste("r_pearson = ", round(cor(d[[a]], d[[b]], use = "complete.obs"), 4), 
                 "\nr_spearman = ", round(cor(d[[a]], d[[b]], use = "complete.obs", method="spearman"), 4),
                 "\nn =",nrow(na.omit(d %>% select(all_of(a),all_of(b))))  )
  }
  if (is.na(main_title)) {
    main_title=paste0("Correlation between ",a," and ",b)
  }

  if ( use_hex) {
    p = ggplot( d, aes_string(x=a, y=b) ) +
      geom_hex(bins=100)
  } else if (use_scattermore) {
    p = ggplot( d, aes_string(x=a, y=b, col=col_, shape=shape_) ) +
      geom_scattermore(aes(alpha=alpha_), pointsize=2) 
  } 
  else {
    p = ggplot( d, aes_string(x=a, y=b, col=col_, shape=shape_) ) +
      geom_point(aes(alpha=alpha_))
  }
  if (main_title!='') {
    p=p+ggtitle(main_title, thecor)
  }
  if ( !is.na(max_x) ) {
    p=p+xlim(0,max_x)+ylim(0, max_x)
  }
  if (draw_diag) {
    p=p+geom_abline(intercept = 0, slope = 1, col="black",linetype="dotted")
  }
  if ( xlog ) {
    p=p+scale_x_log10()+scale_y_log10()
  }
  if ( ! show_legend ) {
    p=p+theme(legend.position="none")
  }
  return(p)
}


gn2tid = function(gn) {
  return(unique( m[['ga']] %>% filter(gene_name == gn) %>% pull(tid) ))
}
tid2coord = function(tid) {
  return(unique( m[['tx']] %>% filter(tid == !!tid) %>% select(chromosome, start, end) ))
}

# __________===============================================================
# Cache ---------------------------------------------------------------

# cache results
# usage: d = cache({x %>% head()}, 'x_head')
cache = function(my_expr, name, rerun=F) {
  xfun::cache_rds(my_expr, rerun=rerun, dir=paste0(params$out_dir, '/cache/'), file=name)
}

clean_cache = function() {
  cached_files = list.files(paste0(params$out_dir, '/cache/'), '_[0-9a-f]{32}[.]rds$', full.names = TRUE)
  unlink(cached_files)
}

gene_type2cat = function(gts) {
  # vectorized
  # see https://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
  ncRNA_cat = c('Mt_rRNA','Mt_tRNA','miRNA','misc_RNA','rRNA','scRNA','snRNA','snoRNA','ribozyme','sRNA','scaRNA')
  lncRNA_cat = c('non_coding','3prime_overlapping_ncRNA','antisense', 'antisense_RNA', 'lincRNA', 'retained_intron','sense_intronic','sense_overlapping','macro_lncRNA','bidirectional_promoter_lncRNA')
  ret=c()
  for(gt in gts) {
    ret=c(ret, case_when(
      gt %in% c('protein_coding','lncRNA') ~ gt,
      grepl('pseudogene', gt) ~ 'pseudogene',
      gt %in% ncRNA_cat ~ 'ncRNA',
      gt %in% lncRNA_cat ~ 'lncRNA',
      gt == 'protein_coding' ~ 'protein_coding',
      is.na(gt) ~ NA_character_,
      TRUE ~ 'other'
    ))
  }
  return(ret)
}


# __________===============================================================
# Colors ------------------------------------------------------------------


# usage p + my_scales()
# all_factors=c("None", "Both", "STAR","HISAT3N","simulated","truth","fast","moderate","slow" )
# pal=brewer.pal(name="Dark2",n=8)
# my_colors=setNames(colorRampPalette(pal)(length(all_factors)), all_factors)
my_colors=c(
  "None"="grey21",
  "Both"="goldenrod4",
  "STAR"="green4",
  "HISAT3N"="darkorange2",
  "simulated"="#551A8B",
  "truth"="black",
  "fast"="#551A8B",
  "moderate"="sienna4",
  "slow"="violetred4",
  "high"="red",
  "medium"="blue4",
  "low"="darkolivegreen",
  "outlier"="red",
  "No Outlier"="grey21"
  "Outlier in both"="goldenrod4",
  "Outlier in STAR"="darkolivegreen",
  "Outlier in HISAT3N"="darkorange4",
  "MERANGS"="violetred4",
  "MOSAIC"="blue4"
)
# show_col(my_colors)
my_scales=function() {
  return(list(
    scale_color_manual(values = my_colors, limits = force),
    scale_fill_manual(values = my_colors, limits = force)
  )
  )
}
