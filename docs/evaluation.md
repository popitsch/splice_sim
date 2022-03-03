# Extra transcript table

Metadata table for transcripts that contains additional columns such as gene_type, len, ilen, GC_frac, T-Pos and F_TC_5min

| transcrtip_id        | rnk  | GC_frac | mappability_median | mappability_mean | T_Pos |
|----------------------|------|---------|--------------------|------------------|-------|
| ENSMUST00000193812.1 | NA   | 0.34206 | 1.000              | 0.93119792       | 998   |

# gene_anno

Metadata table for genes

| gene_id              | transcript_id        | gene_type      | gene_name | transcript_type | transcript_name | level | transcript_support_level | tag                           | havana_transcript    | chromosome |     start |       end | strand | rnk |
|----------------------|----------------------|----------------|-----------|-----------------|-----------------|-------|--------------------------|-------------------------------|----------------------|------------|-----------|-----------|
| ENSMUSG00000102693.1 | ENSMUST00000000001.4 | protein_coding | Gnai3     | protein_coding  | Gnai3-201       | 2     | 1                        | basic,appris_principal_1,CCDS | OTTMUST00000016610.1 | 2          | 108107280 | 108146146 |     - |   9 |


# transcript_info (meta/*.transcript.metadata.tsv)

Self-calculated table mostly for overlapping transcript info that lists for each transcript the following columns

| tid                  | exon_len  | intron_len | n_overlapping | overapping_tids                           |
|----------------------|-----------|------------|---------------|-------------------------------------------|
| ENSMUST00000161581.1 | 250       | 46717      | 2             | ENSMUST00000192973.1,ENSMUST00000194099.1 |

# sj_info (meta/*.SJ.metadata.tsv)

Self-calculated table characterising exon-intron SJ from both sides and their base composition and genomic mappability.

| tid                  | intron_id              | chromosome | start   | end     | strand | don_ex_A | don_ex_C | don_ex_T | don_ex_G | don_in_A | don_in_C | don_in_T | don_in_G | don_win_min_map | don_win_max_map | acc_ex_A | acc_ex_C | acc_ex_T | acc_ex_G | acc_in_A | acc_in_C | acc_in_T | acc_in_G | acc_win_min_map | acc_win_max_map |
|----------------------|------------------------|------------|---------|---------|--------|----------|----------|----------|----------|----------|----------|----------|-----------------|-----------------|-----------------|----------|----------|----------|----------|----------|----------|----------|----------|-----------------|-----------------|
| ENSMUST00000161581.1 | ENSMUST00000161581.1_2 | 1          | 3466688 | 3513404 | +      | 28       | 14       | 28       | 30       | 26       | 14       | 25       | 35       | 1.000           | 1.000           | 32       | 28       | 24       | 16       | 24       | 11       | 39       | 26       | 1.000           | 1.000           |   

# tid_performance (overall_performance/*.tid_performance.tsv)

Mapping performance across the entire transcript with the following columns:

| is_converted_bam | mapper  | condition_id | iso | tid                  | TP  | FP | FN |
|------------------|---------|--------------|-----|----------------------|-----|----|----|
| 1                | HISAT3N | 0.01         | pre | ENSMUST00000210418.1 | 350 | 0  | 0  |

# splice_site_performance (splice_site_performance/*.splice_site_performance.tsv)

Mapping performance across SJs the following columns:

| is_converted_bam | mapper  | condition_id | tid                  | intron_id              | spl_TP | spl_FP | spl_FN | don_TP | don_FP | don_FN | acc_TP | acc_FP | acc_FN | chromosome | start   | end     | strand |
|------------------|---------|--------------|----------------------|------------------------|--------|--------|--------|--------|--------|--------|--------|--------|------------|------------|---------|---------|--------|
| 1                | HISAT3N | 0.01         | ENSMUST00000161581.1 | ENSMUST00000161581.1_2 | 56     | 1      | 13     | 174    | 5      | 0      | 174    | 8      | 1      | 1          | 3466688 | 3513404 | +      |

# splice_site_mappability (splice_site_mappability/*.mappability.tsv)

Mappability tracks of the mapper on the SJ calculated as follows: 1 - abs(acceptorTruthFractionTP - acceptorMappedFractionTP)

| is_converted_bam | mapper  | condition | tid                  | intron_id              | chromosome | start   | end     | strand | don_sj_mappability_TP | don_sj_mappability_FP | don_TP_reads | don_FP_reads | acc_sj_mappability_TP | acc_sj_mappability_FP | acc_TP_reads | acc_FP_reads |
|------------------|---------|-----------|----------------------|------------------------|------------|---------|---------|--------|-----------------------|------------------------------------------|--------------|--------------|-----------------------|-----------------------|--------------|--------------|
| 1                | HISAT3N | 0.01      | ENSMUST00000161581.1 | ENSMUST00000161581.1_2 | 1          | 3466688 | 3513404 | +      | 0.9478632478632478    | 0.4444444444444444    | 117          | 9            | 0.9521655701754386    | 0.7777777777777778    | 114          | 36           |
