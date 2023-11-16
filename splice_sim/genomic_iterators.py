import pysam
import re
from enum import Enum
from pathlib import Path

from more_itertools import peekable
from sortedcontainers import SortedSet


class Location():
    """  A genomic location (interval). Chromosome indices define the order.
         Optional strand property.
    """
    def __init__(self, chr_idx, start, end, strand=None):
        self.chr_idx = chr_idx
        #  start <= end
        self.start = max(start, 0)
        self.end = max(end, 0)
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        self.strand = strand

    @classmethod
    def from_str(cls, loc_string, dict_chr2idx):
        """ parse from chr:start-end. start+end must be >=0 """
        reference, start, end = re.split(':|-', loc_string)
        return cls(dict_chr2idx[reference], int(start), int(end))

    def __lt__(self, other):
        if self.chr_idx != other.chr_idx:
            return self.chr_idx < other.chr_idx
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end

    def __le__(self, other):
        if self.chr_idx != other.chr_idx:
            return self.chr_idx < other.chr_idx
        return self.start <= other.start

    def __len__(self):
        return self.end - self.start + 1

    def to_dict(self):
        return {(None, 'chr_idx'): self.chr_idx, (None, 'start'): self.start, (None, 'end'): self.end}

    def __repr__(self):
        return "%i:%i-%i" % (self.chr_idx, self.start, self.end) + ("" if self.strand is None else "(%s)" % self.strand)

    def to_str(self, dict_idx2chr):
        return "%s:%i-%i" % (dict_idx2chr[self.chr_idx], self.start, self.end) + (
            "" if self.strand is None else "(%s)" % self.strand)

    def __hash__(self):
        return hash(self.__repr__())

    def set_strand(self, strand):
        self.strand = strand
        return self

    def is_stranded(self):
        return self.strand is not None

    def cs_match(self, other):
        """ True if this location is on the same chrom/strand as the passed one """
        return (self.chr_idx == other.chr_idx) and (self.strand == other.strand)

    def __eq__(self, other):
        if not self.cs_match(other):
            return False
        return (self.start == other.start) and (self.end == other.end)

    def left_match(self, other):
        if not self.cs_match(other):
            return False
        return self.start == other.start

    def right_match(self, other):
        if not self.cs_match(other):
            return False
        return self.end == other.end

    def left_pos(self):
        return Location(self.chr_idx, self.start, self.start, strand=self.strand)

    def right_pos(self):
        return Location(self.chr_idx, self.end, self.end, strand=self.strand)

    def overlaps(self, other, strand_specific=True):
        if (strand_specific) and (not self.cs_match(other)):
            return False
        return self.start <= other.end and other.start <= self.end

    @classmethod
    def merge(cls, l):
        """ Merges a list of intervals. If intervals are not on the same chromosome or if strand is not matching, None is returned """
        if l is None:
            return None
        if len(l) == 1:
            return l[0]
        merged = None
        for x in l:
            if merged is None:
                merged = x
            else:
                if not merged.cs_match(x):
                    return None
                merged.start = min(merged.start, x.start)
                merged.end = max(merged.end, x.end)
        return merged

    def is_adjacent(self, other):
        """ true if intervals are directly next to each other (not overlapping!)"""
        if not self.cs_match(other):
            return False
        a, b = (self.end + 1, other.start) if self.end < other.end else (other.end + 1, self.start)
        return a == b

    def split_by_maxwidth(self, maxwidth):
        """ Splits this into n intervals of maximum width """
        k, m = divmod(self.end - self.start + 1, maxwidth)
        ret = [
            Location(self.chr_idx, self.start + i * maxwidth, self.start + (i + 1) * maxwidth - 1, strand=self.strand)
            for i in range(k)]
        if m > 0:
            ret += [Location(self.chr_idx, self.start + k * maxwidth, self.end, strand=self.strand)]
        return ret

    def copy(self):
        """ Deep copy """
        loc = Location(self.chr_idx, self.start, self.end)
        loc.strand = self.strand
        return loc

    def to_file_str(self):
        """ returns a string "<chridx>_<start>_<end>"
            NOTE: strand is not included! Use with care if your annotation contains equal
        """
        return "%i_%i_%i" % (self.chr_idx, self.start, self.end)


# currently supported file formats (tsv supports generic tab-separated file formats).
supported_formats = ['fasta', 'sam', 'bam', 'bed', 'vcf', 'bcf', 'tsv']


def open_pysam_obj(fh, file_format=None,
                   file_extensions={
                       'fasta': ('.fa', '.fasta', '.fna'),
                       'sam': ('.sam'),
                       'bam': ('.bam'),
                       'tsv': ('.tsv', '.tsv.gz'),
                       'bed': ('.bed', '.bed.gz', '.bedgraph', '.bedgraph.gz'),
                       'vcf': ('.vcf', '.vcf.gz'),
                       'bcf': ('.bcf')
                   }):
    """ Ensures and instantiated pysam object and returns the respective file format as string.
    
    If a string or posix path was passed, a respective pysam object is instantiated and was_opened is True.
    If a pysam object was passed, it is returned and was_opened is False.
    
    Parameters
    ----------
    fh : str or pysam object
        if a string/file path is passed, the file type is detected from the file extension and a respective pysam object is instantiated.
         
    file_format : str
        Can be any <supported_formats> or NONE for auto-detection from filename (valid file extensions can be configured).
        Only useful if a str/path is passed.

    Returns
    -------
    file_handle : instance
        pysam object
    file_type : str
        detected file type 
    """
    was_opened = False
    if isinstance(fh, str) or isinstance(fh, Path):
        fh = str(fh)  # convert path to str
        if file_format is None:  # auto detect via file extension
            for ff, ext in file_extensions.items():
                if fh.endswith(ext):
                    file_format = ff
                    break
        if (file_format is None) or (file_format not in supported_formats):
            raise NotImplementedError("Unsupported file format: %s" % file_format)
        # instantiate pysam object        
        if file_format == 'fasta':
            fh = pysam.Fastafile(fh)  # @UndefinedVariable
            was_opened = True
        elif file_format == 'sam':
            fh = pysam.AlignmentFile(fh, "r")  # @UndefinedVariable
            was_opened = True
        elif file_format == 'bam':
            fh = pysam.AlignmentFile(fh, "rb")  # @UndefinedVariable
            was_opened = True
        elif file_format == 'bed':
            fh = pysam.TabixFile(fh, mode="r")  # @UndefinedVariable
        elif file_format == 'tsv':
            fh = pysam.TabixFile(fh, mode="r")  # @UndefinedVariable
        elif file_format == 'vcf':
            fh = pysam.VariantFile(fh, mode="r")  # @UndefinedVariable
        elif file_format == 'bcf':
            fh = pysam.VariantFile(fh, mode="rb")  # @UndefinedVariable
        else:
            raise NotImplementedError("Unsupported input format: %s / %s" % (file_format, str(type(fh))))
    return fh, file_format, was_opened


def get_chrom_dicts(fh, file_format=None):
    """ Extracts dicts from FASTA or BAM files for translating chromosome/chrom_idx and chromosome lengths.
    
    Parameters
    ----------
    in_file : str or pysam AlignmentFile
        A file path or 
         
    file_format : str
        @see open_pysam_obj()
        If an open pysam AlignmentFile file is passed it will not be closed.
        
    Returns
    -------
    tuple (dict, dict, dict)
        dict_chr2idx, dict_idx2chr, dict_chr2len
    
    Raises
    ------
    NotImplementedError
        if input type is not supported yet

    """
    fh, file_format, was_opened = open_pysam_obj(fh, file_format)
    try:
        if isinstance(fh, pysam.Fastafile):  # @UndefinedVariable
            chromosomes = fh.references
            dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
            dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
            dict_chr2len = {c: fh.get_reference_length(c) for c in fh.references}
        elif isinstance(fh, pysam.AlignmentFile):  # @UndefinedVariable
            chromosomes = fh.header.references
            dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
            dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
            dict_chr2len = {c: fh.header.get_reference_length(c) for c in chromosomes}
        elif isinstance(fh, pysam.TabixFile):  # @UndefinedVariable
            chromosomes = fh.contigs
            dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
            dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
            dict_chr2len = None  # no such info in tabix
        elif isinstance(fh, pysam.VariantFile):  # @UndefinedVariable
            chromosomes = list(fh.header.contigs)
            dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
            dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
            dict_chr2len = {c: fh.header.contigs.get(c).length for c in chromosomes}
        else:
            raise NotImplementedError("Unknown input file format %s" % file_format)
    finally:
        if was_opened:
            fh.close()

    return dict_chr2idx, dict_idx2chr, dict_chr2len


# ========================== iterators ==============================================


# BAM flags, @see https://broadinstitute.github.io/picard/explain-flags.html
# @abstract the read is paired in sequencing, no matter whether it is mapped in a pair
BAM_FPAIRED = 0x1
# ! @abstract the read is mapped in a proper pair 
BAM_FPROPER_PAIR = 0x2
# ! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR 
BAM_FUNMAP = 0x4
# ! @abstract the mate is unmapped 
BAM_FMUNMAP = 0x8
# ! @abstract the read is mapped to the reverse strand 
BAM_FREVERSE = 0x10
# ! @abstract the mate is mapped to the reverse strand 
BAM_FMREVERSE = 0x20
# ! @abstract this is read1 
BAM_FREAD1 = 0x40
# ! @abstract this is read2 
BAM_FREAD2 = 0x80
# ! @abstract not primary alignment 
BAM_FSECONDARY = 0x100
# ! @abstract QC failure 
BAM_FQCFAIL = 0x200
# ! @abstract optical or PCR duplicate 
BAM_FDUP = 0x400
# ! @abstract optical or PCR duplicate 
BAM_SUPPLEMENTARY = 0x800


class ReadIterator():
    """ Iterates over a BAM alignment.
        NOTE: bam_file can be a file path (which will open a new pysam.AlignmentFile) or an existing pysam.AlignmentFile object """

    def __init__(self, bam_file, dict_chr2idx, reference, start, end, min_mapping_quality=0,
                 flag_filter=BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY, max_span=None,
                 skip_pcr_dup=False):
        self.close_samfile = isinstance(bam_file, str)
        self.samfile = pysam.AlignmentFile(bam_file, "rb") if self.close_samfile else bam_file  # @UndefinedVariable
        self.start = max(0, start) if start is not None else None  # avoid negative start coordinates!
        self.end = end
        self.dict_chr2idx = dict_chr2idx
        self.reference = reference
        self.chr_idx = self.dict_chr2idx[self.reference] if self.reference is not None else None
        self.flag_filter = flag_filter
        # self.region = Location(self.chr_idx, start, end) if self.chr_idx is not None else None
        self.max_span = max_span
        self.min_mapping_quality = min_mapping_quality
        self.skip_pcr_dup = skip_pcr_dup
        self.TAG_dup = 'xx'

    def __iter__(self) -> (Location, str):
        if self.reference not in self.dict_chr2idx:
            raise StopIteration
        for r in self.samfile.fetch(contig=self.reference,
                                    start=self.start,
                                    end=self.end,
                                    until_eof=True):
            # test flags
            if r.flag & self.flag_filter:
                continue
            if r.mapping_quality < self.min_mapping_quality:
                continue
            if self.skip_pcr_dup and r.has_tag(self.TAG_dup) and r.get_tag(self.TAG_dup) == 1:
                continue
            rloc = Location(self.chr_idx if self.chr_idx is not None else self.dict_chr2idx[r.reference_name],
                            r.reference_start + 1, r.reference_end) if r.reference_name is not None else None
            rloc.set_strand('-' if r.is_reverse else '+')
            # test max_span and drop reads that span larger genomic regions
            if (self.max_span is not None) and (len(rloc) > self.max_span):
                continue
            # test location
            # if rloc.overlaps(self.region):
            yield rloc, r

    def close(self):
        if self.close_samfile:
            self.samfile.close()


class BlockStrategy(Enum):
    LEFT = 1  # same start coord
    RIGHT = 2  # same end coord
    BOTH = 3  # same start and end
    OVERLAP = 4  # overlapping


class BlockLocationIterator():
    """ returns locations and lists of values that share the same location wrt. a 
        given matching strategy (e.g., same start, same end, same coords, overlapping) """

    def __init__(self, it, strategy=BlockStrategy.LEFT, stranded=True):
        self.orgit = it
        self.it = peekable(it)
        self.strategy = strategy

    def __iter__(self) -> (Location, tuple):
        values = None
        locations = None
        for loc, value in self.it:
            mloc = loc.copy()
            values = [value]
            locations = [loc]
            if self.strategy == BlockStrategy.LEFT:
                while self.it.peek(None) and self.it.peek()[0].left_match(mloc):
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc = Location.merge((mloc, l))
            elif self.strategy == BlockStrategy.RIGHT:
                while self.it.peek(None) and self.it.peek()[0].right_match(mloc):
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc = Location.merge((mloc, l))
            elif self.strategy == BlockStrategy.BOTH:
                while self.it.peek(None) and self.it.peek()[0] == mloc:
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc = Location.merge((mloc, l))
            elif self.strategy == BlockStrategy.OVERLAP:
                while self.it.peek(None) and self.it.peek()[0].overlaps(mloc):
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc = Location.merge((mloc, l))
            yield mloc, (locations, values)

    def close(self):
        # print("close BlockLocationIterator")
        try:
            self.orgit.close()
        except AttributeError:
            pass


class BlockedIterator():
    """ gets n BlockLocationIterators. Will yield overlapping (covered) regions """

    def __init__(self, block_its):
        self.iterators = [peekable(it) for it in block_its]
        self.next_pos = SortedSet([it.peek()[0] for it in self.iterators if it.peek(None) is not None])

    def __iter__(self) -> (Location, tuple):
        mloc = None
        while len(self.next_pos) > 0:
            mloc = self.next_pos.pop(0)  # leftmost
            contributers_loc = {}
            contributers_dat = {}
            was_found = True
            while (was_found):
                was_found = False
                for idx, it in enumerate(self.iterators):
                    while it.peek(None) and it.peek()[0].overlaps(mloc):
                        aloc, dat = next(it)
                        mloc = Location.merge((mloc, aloc))
                        contributers_loc[idx] = contributers_loc[idx] + [aloc] if idx in contributers_loc else [aloc]
                        contributers_dat[idx] = contributers_dat[idx] + [dat] if idx in contributers_dat else [dat]
                        was_found = True
            yield mloc, (contributers_loc.keys(), contributers_loc, contributers_dat)
            self.next_pos = SortedSet([it.peek()[0] for it in self.iterators if it.peek(None) is not None])

    def close(self):
        for it in self.iterators:
            try:
                it.close()
            except AttributeError:
                pass


class AnnotationOverlapIterator():
    """ gets one annotation iterator (sorted) and n BlockLocationIterators. Will yield the annotation location and
    all overlapping intervals from the BlockLocationIterators.
    If the annotation iterator is a read iterator then setting check_read_alignment to true will check whether the read actually
    aligns to the annotation intervals.
     """

    def __init__(self, anno_it, block_its, check_read_alignment=False, report_splicing=False):
        self.anno_it = anno_it
        self.check_read_alignment = check_read_alignment
        self.report_splicing = report_splicing
        if check_read_alignment:
            assert isinstance(anno_it, ReadIterator), "Param check_read_alignment is defined only for ReadIterators"
        self.iterators = [peekable(it) for it in block_its]
        self.buffer = {}
        for idx, it in enumerate(self.iterators):
            self.buffer[idx] = []

    def aligns_to(self, readloc, read, loc, report_splicing=False, max_diff=5):
        """ Checks overlap of aligned intervals of the passed read with the passed location."""
        if readloc.chr_idx != loc.chr_idx:
            return False
            # Tests whether read envelops annotation
        envelops = False
        if report_splicing and (read.reference_start < loc.start) and (read.reference_end >= loc.end) and (
                'N' in read.cigarstring):
            envelops = True
        aln_blocks = read.get_blocks()
        lastb = None
        for b in aln_blocks:
            if (b[0] < loc.end and loc.start <= b[1]):
                return True
            if report_splicing and envelops and (lastb is not None):
                splice_block = (lastb[1] + 1, b[0] - 1)
                diff = abs(splice_block[0] - loc.start) + abs(splice_block[1] - loc.end)
                if diff <= max_diff:
                    return True
            lastb = b
        return False

    def __iter__(self) -> (Location, tuple):
        for aloc, a in self.anno_it:
            ret = []
            for idx, it in enumerate(self.iterators):
                # remove non-overlapping items from beginning of list
                self.buffer[idx] = [(loc, dat) for loc, dat in self.buffer[idx] if
                                    (loc.overlaps(aloc, strand_specific=False)) or (loc > aloc)]
                # skip intervals left of annotation
                while it.peek(None) and it.peek()[0] < aloc:
                    loc, dat = next(it)
                    if loc.overlaps(aloc, strand_specific=False):
                        self.buffer[idx] += [(loc, dat)]
                # add new overlapping items
                while it.peek(None) and it.peek()[0].overlaps(aloc, strand_specific=False):
                    loc, dat = next(it)
                    self.buffer[idx] += [(loc, dat)]
                # report current lists (second test for overlap to avoid bug with overhanging annotations
                # example:
                #  |-------------------|
                #    |--|
                #           <anno>
                overlapping = []
                for _, (bloc, bdat) in self.buffer[idx]:
                    for idx in range(0, len(bloc)):
                        if (not self.check_read_alignment) or (
                                self.aligns_to(aloc, a, bloc[idx], self.report_splicing)):
                            overlapping.append((bloc[idx], bdat[idx]))
                ret += [overlapping]
            yield aloc, (a, ret)

    def close(self):
        try:
            self.anno_it.close()
        except AttributeError:
            pass
        for it in self.iterators:
            try:
                it.close()
            except AttributeError:
                pass


class TabixIterator():
    """ Iterates over tabix-indexed files and returns tuples of the column values """

    def __init__(self, tabix_file, reference=None, start=None, end=None):
        self.tabix_file = tabix_file
        self.reference = reference
        self.start = start
        self.end = end

    def __iter__(self):
        f = pysam.TabixFile(self.tabix_file, mode="r")  # @UndefinedVariable
        try:
            # assume last header line has fields
            h = list(f.header)
            if len(h) > 0:
                yield tuple(h[-1].split("\t"))
            # data rows
            for row in f.fetch(reference=self.reference, start=self.start, end=self.end,
                               parser=pysam.asTuple()):  # @UndefinedVariable
                yield tuple(row)
        except:
            raise
        finally:
            f.close()
