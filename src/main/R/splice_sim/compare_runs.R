# mcompare


d1=readRDS('/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big4_slamseq_nf/eva_can_chrom_only/results/data.rds')
d2=readRDS('/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big4_slamseq_nf/results/data.rds')
m=readRDS('/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big4_slamseq_nf/results/meta.rds')
diff = d1[['tx']] %>% 
  left_join(d2[['tx']], 
            by=c('mapper', 'conversion_rate', 'fid', 'true_isoform', 
                 'cv1', 'cv2','se1','se2', 'classification'))

diff %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>% 
  count(mapper, conversion_rate, true_isoform, classification)

diff %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>% ggplot(aes(count.x, count.y, col=true_isoform)) +
  geom_hex(bins=50) +
  facet_wrap(classification~mapper, scales='free')


diff_fx = d1[['fx']] %>% 
  left_join(d2[['fx']], 
            by=c('mapper', 'conversion_rate', 'fid', 'true_isoform', 
                 'cv1', 'cv2','se1','se2', 'classification'))
diff_fx %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>% count(mapper, conversion_rate, true_isoform, classification)

diff_fx %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>% 
  ggplot(aes(count.x, count.y, col=true_isoform)) +
  geom_hex(bins=50) +
  facet_wrap(classification~mapper) +
  xlab("can_chrom_only") + ylab("all chrom") +
  ggtitle('fx, diff only')

diff_sj = d1[['sj']] %>% 
  left_join(d2[['sj']], 
            by=c('mapper', 'conversion_rate', 'fid', 'true_isoform', 
                 'cv1', 'cv2','se1','se2', 'classification', 'class_type')) %>% ungroup()

diff_sj %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>% 
  count(mapper, conversion_rate, true_isoform, classification)

diff_sj %>% mutate(d=count.x- count.y) %>% filter(d!=0) %>%  
  ggplot(aes(count.x, count.y, col=true_isoform)) +
  geom_hex(bins=50) +
  facet_wrap(classification~mapper) +
  xlab("can_chrom_only") + ylab("all chrom") +
  ggtitle('sj, diff only')
