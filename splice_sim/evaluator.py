"""
    @author: @author: niko.popitsch@univie.ac.at

    Evaluation methods for splice_sim data
"""
import concurrent.futures
import statistics
from collections import Counter

import pyBigWig

from splice_sim.genomic_iterators import BlockLocationIterator, AnnotationOverlapIterator
from splice_sim.genomic_iterators import Location, get_chrom_dicts
from splice_sim.genomic_iterators import ReadIterator
from splice_sim.model import *

# RGB codes for read colouring
read_colors = {
    'fn': '200,200,200',
    'fn+fp': '200,200,100',
    'fp+fn': '200,100,0',
    'fp': '200,0,0'
}


def read_aligns_to_loc(loc, read):
    """ Tests whether a read aligns to the passed location """
    # chr is not checked
    for read_start, read_end in read.get_blocks():
        if loc.start <= read_end and read_start < loc.end:
            return True
    return False


def calc_overlap(a, b):
    """ overlap between two intervals """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def read_envelops_loc(iloc, read):
    """ Tests whether a spliced read covers the passed (intronic) location """
    return (read.reference_start < iloc.start) and (read.reference_end >= iloc.end) and ('N' in read.cigarstring)


def read_splices_intron(iloc, read, max_diff=5):
    """ Tests whether the passed read splices the passed intronic) location.
    The maximum difference in coordinates between the intron and the splicing block (N section) can be configured.
    A value>0 allows for some fuzziness """
    os, oe = None, None
    for [start, end] in read.get_blocks():
        if os is not None:
            splice_block = (oe + 1, start - 1)
            diff = abs(splice_block[0] - iloc.start) + abs(splice_block[1] - iloc.end)
            if diff <= max_diff:
                return True
        os, oe = start, end
    return False


def get_class(read, features):
    """ Calculate the classification of the passed read wrt different classification groups """
    classification = {
        'tx': set(),
        'fx': set(),
        'don': set(),
        'acc': set(),
        'spl': set()}
    for feat_loc, (tid, ftype, fid) in features:  # NB: chrom is not checked, must be matched beforehand!
        strand_match = feat_loc.strand == '-' if read.is_reverse else feat_loc.strand == '+'
        if not strand_match:
            continue
        # does read align into feature?
        if read_aligns_to_loc(feat_loc, read):  # NB no need to check chr here
            classification['tx'].add(tid)
            classification['fx'].add(fid)
            if ftype == 'intron':
                if read.reference_start < feat_loc.start <= read.reference_end:
                    classification['don' if feat_loc.strand == '+' else 'acc'].add(fid)
                elif read.reference_start < feat_loc.end <= read.reference_end:
                    classification['acc' if feat_loc.strand == '+' else 'don'].add(fid)
        # does read span over feature?
        if read_envelops_loc(feat_loc, read) and ftype == 'intron' and read_splices_intron(feat_loc, read):
            classification['spl'].add(fid)
    return classification


def set_compare(truth, truth_ext, found):
    """ Compares two sets and returns TP, FP, FNs.
    truth contains the classification of the true read wrt. features of the true tx,
    found contains the classification of the found read wrt. overlapping features and
    truth_ext contains the classification of the true read wrt. the union of features overlapping the found read and true tx features. """
    TP = truth & found
    FP = found.difference(truth_ext)
    FN = truth.difference(found)
    return TP, FP, FN


classtype2cat = {
    'tx': 'tx',
    'fx': 'fx',
    'don': 'sj',
    'acc': 'sj',
    'spl': 'sj'}
empty_class = {'tx': set(), 'fx': set(), 'don': set(), 'acc': set(), 'spl': set()}


def classify_read(found_read, found_overlapping_raw, tid2feat, performance, dict_samout, dict_chr2idx, min_mapq=20):
    """ Classifies a read as TP/FP/FN """
    # recreate true read
    true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = found_read.query_name.split(
        '_')
    true_read = pysam.AlignedSegment()
    true_read.query_name = found_read.query_name
    true_read.cigarstring = true_cigar
    true_read.reference_start = int(true_start)
    true_read.is_reverse = true_strand == '-'
    cv1 = 1 if int(n_converted) > 0 else 0
    cv2 = 1 if int(n_converted) > 1 else 0
    se1 = 1 if int(n_seqerr) > 0 else 0
    se2 = 1 if int(n_seqerr) > 1 else 0
    # which features does the true read overlap with and whats the true iclass
    true_features = set(tid2feat[true_tid])
    true_class = get_class(true_read, true_features)
    true_features_ext = set(tid2feat[true_tid] + found_overlapping_raw)
    true_class_ext = true_class if true_features_ext == true_features else get_class(true_read, true_features_ext)
    found_class = {}
    found_class[0] = empty_class if found_read.is_unmapped else get_class(found_read, found_overlapping_raw)
    found_class[1] = empty_class if ((found_read.is_unmapped) or (found_read.mapping_quality < min_mapq)) else \
    found_class[0]
    # calculate TP,FP,FN (+/- mq filtering)
    for mq_fil in [0, 1]:
        for class_type in true_class.keys():
            TP, FP, FN = set_compare(true_class[class_type], true_class_ext[class_type],
                                     found_class[mq_fil][class_type])
            for fid in TP:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'TP'] += 1.0
            for fid in FN:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FN'] += 1.0
            for fid in FP:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FP'] += 1.0 / len(FP)
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FP_raw'] += 1.0
            # write reads (no mq filtering!)  if mismapped
            if (mq_fil == 0) and (dict_samout is not None) and (len(FN) + len(FP) > 0):
                samout = dict_samout[classtype2cat[class_type]]
                # write FN
                if len(FN) > 0:
                    is_correct_strand = (found_read.is_reverse and true_strand == '-') or (
                                (not found_read.is_reverse) and true_strand == '+')
                    true_read.query_sequence = found_read.query_sequence if is_correct_strand else reverse_complement(
                        found_read.query_sequence)
                    true_read.query_qualities = found_read.query_qualities if is_correct_strand else list(
                        reversed(found_read.query_qualities))
                    true_read.reference_id = dict_chr2idx[true_chr]
                    true_read.tags = (("NM", n_seqerr + n_converted),)
                    true_read.flag = 0 if true_strand == '+' else 16
                    true_read.set_tag(tag='YC', value=read_colors['fn+fp'] if len(FP) > 0 else read_colors['fn'],
                                      value_type="Z")
                    samout.write(true_read)
                # write FP
                if len(FP) > 0:
                    found_read.set_tag(tag='YC', value=read_colors['fp+fn'] if len(FN) > 0 else read_colors['fp'],
                                       value_type="Z")
                    samout.write(found_read)
    return performance


def classify_region(target_region, config, bam_file, chrom2feat, tid2feat, filter_regions, out_file_prefix):
    """ Classify all reads in the passed target region as TP/FP/FN """
    print("Classify region", target_region)
    if target_region == 'unmapped':
        filter_tree = None
    else:
        filter_tree = filter_regions[target_region.chr_idx] if ((filter_regions is not None) and (target_region.chr_idx in filter_regions)) else None

    min_mapq = config['min_mapq'] if 'min_mapq' in config else 20
    write_bam = config['write_bam'] if 'write_bam' in config else True
    samin = pysam.AlignmentFile(bam_file, "rb")
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    dict_samout = {
        'tx': pysam.AlignmentFile(out_file_prefix + '.tx.eval.%s.bam' % slugify(str(target_region)), "wb",
                                  template=samin),
        'fx': pysam.AlignmentFile(out_file_prefix + '.fx.eval.%s.bam' % slugify(str(target_region)), "wb",
                                  template=samin),
        'sj': pysam.AlignmentFile(out_file_prefix + '.sj.eval.%s.bam' % slugify(str(target_region)), "wb",
                                  template=samin)
    } if write_bam else None
    performance = Counter()
    wr = 0
    fr = 0
    if target_region == 'unmapped':
        for unmapped_read in samin.fetch(contig=None, until_eof=True):
            if not unmapped_read.is_unmapped:
                continue
            performance = classify_read(unmapped_read, [], tid2feat, performance, dict_samout, dict_chr2idx, min_mapq)
    else:
        chrom = dict_idx2chr[target_region.chr_idx]
        if chrom in samin.references:  # there are reads on this chrom
            aits = [BlockLocationIterator(
                iter(sorted(chrom2feat[chrom] if chrom in chrom2feat else [], key=lambda x: x[0])))]
            rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=target_region.start, end=target_region.end,
                               max_span=None, flag_filter=0)  # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits, check_read_alignment=True, report_splicing=True)
            for loc, (read, annos) in it:
                if (filter_tree is not None) and (filter_tree.overlap(loc.start, loc.end)):
                    fr += 1
                    continue
                wr += 1
                overlapping_annos = annos[0]  # select 1st entry as there is only one ait
                performance = classify_read(read, overlapping_annos, tid2feat, performance, dict_samout, dict_chr2idx,
                                            min_mapq)
    # close SAM files
    samin.close()
    if dict_samout is not None:
        for cat, samout in dict_samout.items():
            samout.close()
            dict_samout[cat] = out_file_prefix + '.%s.eval.%s.bam' % (cat, slugify(str(target_region)))
    return performance, dict_samout, wr, fr


def evaluate(config, m, bam_file, filter_regions_bed, threads, out_dir):
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    # parse cr/mapper from file name
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4]  # file ends with .bam or .final.bam
    cr, mapper = tmp[tmp.find('.cr') + 3:].rsplit('.', 1)
    out_file_prefix = out_dir + '/' + config['dataset_name'] + '.cr' + cr + '.' + mapper
    cr = float(cr)
    min_mapq = config['min_mapq'] if 'min_mapq' in config else 20

    # get config params
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    chromosomes = config['chromosomes'] if 'chromosomes' in config else dict_chr2len.keys()
    # collect= features for chrom
    chrom2feat = OrderedDict()
    tid2feat = OrderedDict()
    for tid, t in m.transcripts.items():
        c = t.chromosome
        if c not in chrom2feat:
            chrom2feat[c] = []
        if tid not in tid2feat:
            tid2feat[tid] = []
        exon_len = 0
        for exon in t.exons:
            fid = tid + "_ex" + str(exon['exon_number'])
            loc = Location(dict_chr2idx[t.chromosome], exon['start'], exon['end'], strand=t.strand)
            chrom2feat[c].append((loc, (tid, 'exon', fid)))
            tid2feat[tid].append((loc, (tid, 'exon', fid)))
            exon_len += len(loc)
        intron_len = 0
        for intron in t.introns:
            loc = Location(dict_chr2idx[t.chromosome], intron['start'], intron['end'], strand=t.strand)
            fid = tid + "_in" + str(intron['exon_number'] - 1 if t.strand == '+' else intron['exon_number'])
            chrom2feat[c].append((loc, (tid, 'intron', fid)))
            tid2feat[tid].append((loc, (tid, 'intron', fid)))
            intron_len += len(loc)
    # get filter regions if any
    filter_regions = None
    if filter_regions_bed is not None:
        freg = pd.read_csv(filter_regions_bed, sep='\t', header=None,
                           names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
        filter_regions = {}
        for _, row in freg.iterrows():
            if row['chromosome'] not in dict_chr2idx.keys():
                print("Skipping unknown chrom %s in filter_regions_bed" % row['chromosome'])
                continue
            chridx = dict_chr2idx[row['chromosome']]
            if chridx not in filter_regions:
                filter_regions[chridx] = IntervalTree()
            filter_regions[chridx].add(Interval(row['start'], row['end']))

    # target regions
    if 'target_region' in config:
        target_regions = [Location.from_str(config['debug_region'], dict_chr2idx)]
    else:
        target_regions = [Location(dict_chr2idx[c], 1, c_len) for c, c_len in dict_chr2len.items() if
                          c in chromosomes] + ['unmapped']  # unmapped reads
    print("Considered target regions: ", target_regions)
    # run parallel
    performance = Counter()
    dict_samout = None
    wr, fr = 0, 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_param = {executor.submit(classify_region, reg, config, bam_file, chrom2feat, tid2feat, filter_regions,
                                        out_file_prefix): (reg) for reg in target_regions}
        pbar = tqdm.tqdm(total=len(target_regions), desc="Annotated regions")
        for future in concurrent.futures.as_completed(future_param):
            para = future_param[future]
            region_performance, region_dict_samout, region_wr, region_fr = future.result()
            wr += region_wr
            fr += region_fr
            performance.update(region_performance)
            if region_dict_samout is not None:
                if dict_samout is None:
                    dict_samout = {'tx': [], 'fx': [], 'sj': []}
                for k in dict_samout.keys():
                    dict_samout[k] += [region_dict_samout[k]]
            pbar.set_description('Classified reads: %i, filtered: %i' % (wr, fr))
            pbar.update(1)
        pbar.close()
    # write performance data
    with open(out_file_prefix + '.counts.tsv', 'w') as out_raw:
        with open(out_file_prefix + '.counts.mq%i.tsv' % min_mapq, 'w') as out_fil:
            header = '\t'.join(str(x) for x in [
                'mapper', 'conversion_rate', 'mq_fil', 'class_type', 'fid', 'true_isoform', 'cv1', 'cv2', 'se1', 'se2',
                'classification', 'count'
            ])
            print(header, file=out_raw)
            print(header, file=out_fil)
            for (
            class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, classification), count in performance.items():
                out = out_raw if mq_fil == 0 else out_fil
                print('\t'.join(str(x) for x in [
                    mapper, cr, mq_fil, class_type, fid, true_isoform, cv1, cv2, se1, se2, classification, count
                ]), file=out)
    # merge, sort and index BAM files
    if dict_samout is not None:
        for cat, bam_files in dict_samout.items():
            out_file = out_file_prefix + '.%s.eval.bam' % str(cat)
            merge_bam_files(out_file, bam_files, sort_output=True, del_in_files=True)
    print("All done. Classified %i and filtered %i reads" % (wr, fr))


def calc_mappability(genomeMappability, chrom, start, end):
    """ Calculate genomic mappability for given window """
    if genomeMappability is None:
        return np.nan
    try:
        map = flatten([[float(x[3])] * (min(end + 1, int(x[2])) - max(start, int(x[1]))) for x in
                       genomeMappability.fetch('chr' + chrom, start, end, parser=pysam.asTuple())])
        map = map + [0] * (end - start + 1 - len(map))  # pad with missing zeros
        mea_map = np.mean(map) if len(map) > 0 else 0
    except ValueError as ve:
        # print(ve)
        mea_map = np.nan
    return mea_map


def calc_conservation(genomeConservation, chrom, start, end):
    """ Calculate genomic conservation for given window """
    if genomeConservation is None:
        return np.nan
    try:
        cons = genomeConservation.stats('chr' + chrom, start, end + 1).pop()
    except RuntimeError as re:
        # print(re)
        cons = np.nan
    return cons


# extract_feature_metadata
def extract_feature_metadata(config, m, out_dir):
    """ extract metadata fox tx/fx/sj """
    genomeMappability = pysam.TabixFile(config['genome_mappability'], mode="r") if 'genome_mappability' in config else None
    genomeConservation = pyBigWig.open(config['genome_conservation']) if 'genome_conservation' in config else None
    readlen = config['readlen']
    # calculate transcript overlap
    tid2info = {}
    tiv = m.build_transcript_iv()
    with open(out_dir + '/' + config['dataset_name'] + '.tx.metadata.tsv', 'w') as out_tx:
        with open(out_dir + '/' + config['dataset_name'] + '.fx.metadata.tsv', 'w') as out_fx:
            with open(out_dir + '/' + config['dataset_name'] + '.sj.metadata.tsv', 'w') as out_sj:
                print("""tid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """A\tC\tT\tG\tmean_map\tmean_cons\texon_len\tintron_len\tn_overlapping\toverlapping_tids""",
                      file=out_tx)
                print("""tid\tfid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """A\tC\tT\tG\tmean_map\tmean_cons""",
                      file=out_fx)
                print("""tid\tfid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """don_ex_A\tdon_ex_C\tdon_ex_T\tdon_ex_G\tdon_in_A\tdon_in_C\tdon_in_T\tdon_in_G\tdon_win_map\tdon_win_cons\t""" +
                      """acc_ex_A\tacc_ex_C\tacc_ex_T\tacc_ex_G\tacc_in_A\tacc_in_C\tacc_in_T\tacc_in_G\tacc_win_map\tacc_win_cons""",
                      file=out_sj)
                for tid, t in m.transcripts.items():
                    # get RNA subseq from transcript seq. NOTE: may be shorter than 2xreadlen
                    rna_seq = reverse_complement(t.transcript_seq) if t.strand == "-" else t.transcript_seq
                    # calc mean mappability
                    tx_map = calc_mappability(genomeMappability, t.chromosome, t.start, t.end)
                    tx_cons = calc_conservation(genomeConservation, t.chromosome, t.start, t.end)
                    # get overlapping tids
                    overlapping = [x.data.tid for x in tiv[t.chromosome].overlap(t.start, t.end) if x.data.tid != tid]
                    l = len(rna_seq)
                    exon_len, intron_len = 0, 0
                    # exons
                    for exon in t.exons:
                        rnk = exon['exon_number']
                        fid = tid + "_ex" + str(rnk)
                        ex_map = calc_mappability(genomeMappability, t.chromosome, exon['start'], exon['end'])
                        ex_cons = calc_conservation(genomeConservation, t.chromosome, exon['start'], exon['end'])
                        ex_rna = rna_seq[exon['start'] - t.start:exon['end'] - t.start + 1]
                        exon_len += exon['end'] - exon['start'] + 1
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,
                            'exon',
                            rnk,
                            t.chromosome,
                            exon['start'],
                            exon['end'],
                            exon['strand'],
                            ex_rna.count("A"),
                            ex_rna.count("C"),
                            ex_rna.count("T"),
                            ex_rna.count("G"),
                            ex_map,
                            ex_cons
                        ]]), file=out_fx)
                    # introns
                    for intron in t.introns:
                        rnk = intron['exon_number'] - 1 if t.strand == '+' else intron['exon_number']
                        fid = tid + "_in" + str(rnk)
                        # tx
                        in_map = calc_mappability(genomeMappability, t.chromosome, intron['start'], intron['end'])
                        in_cons = calc_conservation(genomeConservation, t.chromosome, intron['start'], intron['end'])
                        in_rna = rna_seq[intron['start'] - t.start:intron['end'] - t.start + 1]
                        intron_len += intron['end'] - intron['start'] + 1
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,
                            'intron',
                            rnk,
                            t.chromosome,
                            intron['start'],
                            intron['end'],
                            intron['strand'],
                            in_rna.count("A"),
                            in_rna.count("C"),
                            in_rna.count("T"),
                            in_rna.count("G"),
                            in_map,
                            in_cons
                        ]]), file=out_fx)
                        # SJ
                        # note that win can be shorter than readlen if we run out of sequence!
                        don_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] \
                            if intron['strand'] == "+" else [intron['end'] - readlen, intron['end'] + readlen]
                        acc_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] \
                            if intron['strand'] == "-" else [intron['end'] - readlen, intron['end'] + readlen]
                        don_win_map = calc_mappability(genomeMappability, t.chromosome, don_win_genomic[0],
                                                       don_win_genomic[1])
                        acc_win_map = calc_mappability(genomeMappability, t.chromosome, acc_win_genomic[0],
                                                       acc_win_genomic[1])
                        don_win_cons = calc_conservation(genomeConservation, t.chromosome, don_win_genomic[0],
                                                         don_win_genomic[1])
                        acc_win_cons = calc_conservation(genomeConservation, t.chromosome, acc_win_genomic[0],
                                                         acc_win_genomic[1])
                        don_seq_rna_ex = rna_seq[max(0, t.end - intron['end'] - readlen):min(t.end - intron['end'], l)] \
                            if intron['strand'] == "-" else \
                                rna_seq[max(0, intron['start'] - readlen - t.start):min(l, intron['start'] - t.start)]
                        don_seq_rna_in = rna_seq[max(0, t.end - intron['end']):min(t.end - intron['end'] + readlen, l)] \
                            if intron['strand'] == "-" else \
                            rna_seq[max(0, intron['start'] - t.start):min(l, intron['start'] + readlen - t.start)]
                        acc_seq_rna_in = rna_seq[max(0, t.end - intron['start'] - readlen + 1):min(t.end - intron['start'] + 1,l)] \
                            if intron['strand'] == "-" else \
                            rna_seq[max(0,intron['end'] - readlen - t.start + 1):min(l,intron['end'] - t.start + 1)]
                        acc_seq_rna_ex = rna_seq[max(0, t.end - intron['start'] + 1):min(t.end - intron['start'] + readlen + 1,l)] \
                            if intron['strand'] == "-" else \
                            rna_seq[max(0,intron['end'] - t.start + 1):min(l,intron['end'] + readlen - t.start + 1)]
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,
                            'intron',
                            rnk,
                            t.chromosome,
                            intron['start'],
                            intron['end'],
                            intron['strand'],
                            # donor
                            don_seq_rna_ex.count("A"),
                            don_seq_rna_ex.count("C"),
                            don_seq_rna_ex.count("T"),
                            don_seq_rna_ex.count("G"),
                            don_seq_rna_in.count("A"),
                            don_seq_rna_in.count("C"),
                            don_seq_rna_in.count("T"),
                            don_seq_rna_in.count("G"),
                            don_win_map,
                            don_win_cons,
                            # acceptor
                            acc_seq_rna_ex.count("A"),
                            acc_seq_rna_ex.count("C"),
                            acc_seq_rna_ex.count("T"),
                            acc_seq_rna_ex.count("G"),
                            acc_seq_rna_in.count("A"),
                            acc_seq_rna_in.count("C"),
                            acc_seq_rna_in.count("T"),
                            acc_seq_rna_in.count("G"),
                            acc_win_map,
                            acc_win_cons
                        ]]), file=out_sj)
                    # tx info
                    print("\t".join([str(x) for x in [
                        tid,
                        'tx',
                        len(t.exons),
                        t.chromosome,
                        t.start,
                        t.end,
                        t.strand,
                        rna_seq.count("A"),
                        rna_seq.count("C"),
                        rna_seq.count("T"),
                        rna_seq.count("G"),
                        tx_map,
                        tx_cons,
                        exon_len,
                        intron_len,
                        len(overlapping),
                        ','.join(overlapping)
                    ]]), file=out_tx)
    print("All done.")


def extract_bam_stats(bam_file, out_dir):
    """ Extract stats per BAM file """
    stats = Counter()
    logging.info("extract_bam_stats for %s" % bam_file)
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    # parse cr/mapper from file name
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    if bam_file.endswith('.eval.bam'):
        tmp = bam_file[:-9]
        tmp, ftype = tmp.rsplit('.', 1)
    elif bam_file.endswith('.final.bam'):
        ftype = 'NA'
        tmp = bam_file[:-10]
    else:
        ftype = 'NA'
        tmp = bam_file[:-4]
    cr, mapper = tmp[tmp.find('.cr') + 3:].rsplit('.', 1)
    cr = float(cr)
    samin = pysam.AlignmentFile(bam_file, "rb")
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    cols = ['n_reads', 'n_spliced_reads', 'n_softclipped_reads', 'mean_span_reads', 'len_spliced_reads', 'n_unmapped']
    with open(out_dir + '/' + Path(bam_file).stem + '.bamstats.tsv', 'w') as out:
        print("\t".join([str(x) for x in ['chromosome', 'mapper', 'cr', 'ftype'] + cols]), file=out)
        for chrom, chrlen in dict_chr2len.items():
            rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=1, end=chrlen, max_span=None,
                               flag_filter=0)  # max_span=m.max_ilen
            stats[chrom, 'mean_span_reads'] = []
            stats[chrom, 'len_spliced_reads'] = []
            for loc, r in rit:
                stats[chrom, 'n_reads'] += 1
                stats[chrom, 'mean_span_reads'].append(r.reference_end - r.reference_start + 1)
                maxn = 0
                is_spliced = 0
                is_softclipped = 0
                for op, l in r.cigartuples:
                    if op == 3:
                        maxn = max(l, maxn)  # N
                        is_spliced = 1
                    if op == 4:
                        is_softclipped = 1
                if maxn > 0:
                    stats[chrom, 'len_spliced_reads'].append(maxn)
                stats[chrom, 'n_spliced_reads'] += is_spliced
                stats[chrom, 'n_softclipped_reads'] += is_softclipped
            stats[chrom, 'mean_span_reads'] = statistics.mean(stats[chrom, 'mean_span_reads']) if len(
                stats[chrom, 'mean_span_reads']) > 0 else 'NA'
            stats[chrom, 'len_spliced_reads'] = statistics.mean(stats[chrom, 'len_spliced_reads']) if len(
                stats[chrom, 'len_spliced_reads']) > 0 else 'NA'
        # unmapped reads
        for unmapped_read in samin.fetch(contig=None, until_eof=True):
            if unmapped_read.is_unmapped:
                stats['unmapped', 'n_unmapped'] += 1
        # write results
        for chrom in list(dict_chr2len.keys()) + ['unmapped']:
            print("\t".join([str(x) for x in [
                chrom, mapper, cr, ftype] + [stats[chrom, col] for col in cols]]), file=out)


def extract_bam_stats(bam_file, out_dir):
    """ Extract stats per BAM file """
    stats = Counter()
    logging.info("extract_bam_stats for %s" % bam_file)
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    # parse cr/mapper from file name
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    if bam_file.endswith('.eval.bam'):
        tmp = bam_file[:-9]
        tmp, ftype = tmp.rsplit('.', 1)
    elif bam_file.endswith('.final.bam'):
        ftype = 'NA'
        tmp = bam_file[:-10]
    else:
        ftype = 'NA'
        tmp = bam_file[:-4]
    cr, mapper = tmp[tmp.find('.cr') + 3:].rsplit('.', 1)
    cr = float(cr)
    samin = pysam.AlignmentFile(bam_file, "rb")
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    cols = ['n_reads', 'n_spliced_reads', 'n_softclipped_reads', 'mean_span_reads', 'len_spliced_reads', 'n_unmapped']
    with open(out_dir + '/' + Path(bam_file).stem + '.bamstats.tsv', 'w') as out:
        print("\t".join([str(x) for x in ['chromosome', 'mapper', 'cr', 'ftype'] + cols]), file=out)
        for chrom, chrlen in dict_chr2len.items():
            rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=1, end=chrlen, max_span=None,
                               flag_filter=0)  # max_span=m.max_ilen
            stats[chrom, 'mean_span_reads'] = []
            stats[chrom, 'len_spliced_reads'] = []
            for loc, r in rit:
                stats[chrom, 'n_reads'] += 1
                stats[chrom, 'mean_span_reads'].append(r.reference_end - r.reference_start + 1)
                maxn = 0
                is_spliced = 0
                is_softclipped = 0
                for op, l in r.cigartuples:
                    if op == 3:
                        maxn = max(l, maxn)  # N
                        is_spliced = 1
                    if op == 4:
                        is_softclipped = 1
                if maxn > 0:
                    stats[chrom, 'len_spliced_reads'].append(maxn)
                stats[chrom, 'n_spliced_reads'] += is_spliced
                stats[chrom, 'n_softclipped_reads'] += is_softclipped
            stats[chrom, 'mean_span_reads'] = statistics.mean(stats[chrom, 'mean_span_reads']) if len(
                stats[chrom, 'mean_span_reads']) > 0 else 'NA'
            stats[chrom, 'len_spliced_reads'] = statistics.mean(stats[chrom, 'len_spliced_reads']) if len(
                stats[chrom, 'len_spliced_reads']) > 0 else 'NA'
        # unmapped reads
        for unmapped_read in samin.fetch(contig=None, until_eof=True):
            if unmapped_read.is_unmapped:
                stats['unmapped', 'n_unmapped'] += 1
        # write results
        for chrom in list(dict_chr2len.keys()) + ['unmapped']:
            print("\t".join([str(x) for x in [
                chrom, mapper, cr, ftype] + [stats[chrom, col] for col in cols]]), file=out)


def special_count_HISAT3N_strand_stats(bam_file, out_dir):
    """ Extract HISAT3N specific stats """
    stats = Counter()
    logging.info("extract_bam_stats for %s" % bam_file)
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    # parse cr/mapper from file name
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    if bam_file.endswith('.eval.bam'):
        tmp = bam_file[:-9]
        tmp, ftype = tmp.rsplit('.', 1)
    elif bam_file.endswith('.final.bam'):
        ftype = 'NA'
        tmp = bam_file[:-10]
    else:
        ftype = 'NA'
        tmp = bam_file[:-4]
    cr, mapper = tmp[tmp.find('.cr') + 3:].rsplit('.', 1)
    cr = float(cr)
    samin = pysam.AlignmentFile(bam_file, "rb")
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    cols = ['n_reads', 'n_softclipped_reads']
    with open(out_dir + '/' + Path(bam_file).stem + '.bamstats.tsv', 'w') as out:
        print("\t".join([str(x) for x in ['chromosome', 'mapper', 'cr', 'ftype', 'matching_strand'] + cols]), file=out)
        for chrom, chrlen in dict_chr2len.items():
            rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=1, end=chrlen, max_span=None,
                               flag_filter=0)  # max_span=m.max_ilen
            for loc, r in rit:
                true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = r.query_name.split(
                    '_')
                matching_strand = ((r.is_reverse) and (true_strand == '-')) or (
                            (not r.is_reverse) and (true_strand == '+'))
                is_softclipped = 4 in [x for x, _ in r.cigartuples]
                stats[chrom, matching_strand, 'n_reads'] += 1
                stats[chrom, matching_strand, 'n_softclipped_reads'] += is_softclipped
        # write results
        for chrom in list(dict_chr2len.keys()):
            for matching_strand in [True, False]:
                print("\t".join([str(x) for x in [
                    chrom, mapper, cr, ftype, 1 if matching_strand else 0] + [stats[chrom, matching_strand, col] for col
                                                                              in cols]]), file=out)


def special_extract_sc_reads(bam_file, region, outdir):
    """ Extract soft-clipped reads """
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    out_file = outdir + '/softclipped_' + os.path.basename(bam_file)
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    samin = pysam.AlignmentFile(bam_file, "rb")
    samout = pysam.AlignmentFile(out_file, "wb", template=samin)
    reg = Location.from_str(region, dict_chr2idx)
    rit = ReadIterator(samin, dict_chr2idx, reference=dict_idx2chr[reg.chr_idx], start=reg.start, end=reg.end,
                       max_span=None, flag_filter=0)  # max_span=m.max_ilen
    for loc, r in rit:
        is_softclipped = 4 in [x for x, _ in r.cigartuples]
        if is_softclipped:
            samout.write(r)
    samout.close()
    # index
    try:
        pysam.index(out_file)  # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    print("All done.")


def filter_xc_yc_reads(bam_file, outdir):
    """ extract reads where xc/yc counts differ """
    assert os.path.exists(bam_file) and os.path.exists(bam_file + '.bai'), "Could not find bam file (or idx) for %s" % (
        bam_file)
    out_bam_file = outdir + '/xcyc_' + os.path.basename(bam_file)
    dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(bam_file)
    samin = pysam.AlignmentFile(bam_file, "rb")
    samout = pysam.AlignmentFile(out_bam_file, "wb", template=samin)
    hist = Counter()
    for chrom, chrlen in dict_chr2len.items():
        rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=1, end=chrlen, max_span=None,
                           flag_filter=0)  # max_span=m.max_ilen
        for loc, r in rit:
            if r.has_tag('yc') and r.has_tag('YC'):
                rtype = r.get_tag('YC')
                if rtype == read_colors['fn']:
                    continue  # skip FN reads
                _, _, _, _, _, _, _, _, xc, _ = r.query_name.split('_')
                yc = r.get_tag('yc')
                hist[int(xc), int(yc)] += 1
                if int(xc) != int(yc):
                    samout.write(r)
    samout.close()
    out_file = outdir + '/xcyc_' + os.path.basename(bam_file) + '.hist.tsv'
    with open(out_file, 'w') as out:
        print('bam\txc\tyc\tcount', file=out)
        for xc, yc in hist:
            print('\t'.join([str(x) for x in [os.path.basename(bam_file), xc, yc, hist[xc, yc]]]), file=out)
    # index
    try:
        pysam.index(out_bam_file)  # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    print("All done.")
