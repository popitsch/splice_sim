from enum import Enum
from more_itertools import peekable
import pysam
from sortedcontainers import SortedSet


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


class Location():
    """  A genomic location (interval). chromosome indices define the order. """
    def __init__(self, chr_idx, start, end):
        self.chr_idx = chr_idx
        #  start <= end
        self.start = min(start, end)
        self.end = max(start, end)
    @classmethod
    def from_str(cls, loc_string, dict_chr2idx):
        """ parse from chr:start-end. m """
        reference, tmp = loc_string.split(':')
        start, end = tmp.split('-')
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
    def len(self):
        return self.end - self.start + 1
    def to_dict(self):
        return {(None, 'chr_idx'): self.chr_idx, (None, 'start'): self.start, (None, 'end'): self.end }
    def __repr__(self):
        return "%i:%i-%i" % (self.chr_idx, self.start, self.end)
    def to_str(self, dict_idx2chr):
        return "%s:%i-%i" % (dict_idx2chr[self.chr_idx], self.start, self.end)
    def __hash__(self):
        return hash(self.__repr__())
    def __eq__(self, other):
        if not isinstance(other, Location):
            return False
        return ((self.chr_idx == other.chr_idx) and (self.start == other.start) and (self.end == other.end))
    def left_match(self, other):
        if not isinstance(other, Location):
            return False
        if self.chr_idx != other.chr_idx:
            return False
        return self.start == other.start
    def right_match(self, other):
        if not isinstance(other, Location):
            return False
        if self.chr_idx != other.chr_idx:
            return False
        return self.end == other.end
    def overlaps(self, other):
        if not isinstance(other, Location):
            return False
        if self.chr_idx != other.chr_idx:
            return False
        return self.start <= other.end and other.start <= self.end
    @classmethod
    def merge(cls, a, b):
        if isinstance(a, Location) and isinstance(b, Location) and a.chr_idx == b.chr_idx:
            start = min(a.start, b.start)
            end = max(a.end, b.end) 
            return Location(a.chr_idx, start, end)          
        return None
    def is_adjacent(self, other):
        """ true if intervals are directly next to each other (not overlapping!)"""
        if isinstance(other, Location) and self.chr_idx == other.chr_idx:
            a,b = (self.end+1, other.start) if self.end < other.end else (other.end+1, self.start)
            return a == b
        return False
    def copy(self):
        return Location(self.chr_idx, self.start, self.end)
    def to_file_str(self):
        return "%i_%i_%i" % (self.chr_idx, self.start, self.end)

class TabixIterator():
    """ Iterates over tabix-indexed files and returns tuples of the column values """
    def __init__(self, tabix_file, reference=None, start=None, end=None):
        self.tabix_file = tabix_file
        self.reference = reference
        self.start = start
        self.end = end
    def __iter__(self):
        f = pysam.TabixFile(self.tabix_file, mode="r") # @UndefinedVariable
        try:
            # assume last header line has fields
            h = list(f.header)
            if len(h) > 0:
                yield tuple(h[-1].split("\t"))
            # data rows
            for row in f.fetch(reference=self.reference, start=self.start, end=self.end, parser=pysam.asTuple()): # @UndefinedVariable
                yield tuple(row) 
        except:
            raise
        finally:
            f.close()


class LocationTabixIterator(TabixIterator):
    """ Iterates over tabix-indexed files and returns location/value pairs.
        Requires a tabix-indexed file with a numeric (float) value in the 4th column and a dict that converts chromosome names
        to indices. Chromosome 'chr' prefixes will be removed automatically if the chomosome names in the passed dict do not have them.
        The iterated part of the input file can be defined with reference, start and end. 
        
        Example usage:
        it1 = LocationIterator('mydata.bedgraph.gz', dict_chr2idx, reference='chr1')
        it2 = LocationIterator('mydata.bedgraph.gz', dict_chr2idx, reference='chr2')
        for (loc1, value1), (loc2, value2) in zip(it1, it2):
            print("%s\t%f / %s\t%f" % (loc1, value1, loc2, value2))
            
    """
    def __init__(self, tabix_file, dict_chr2idx, reference=None, start=None, end=None, chr_prefix_operator=None):
        super().__init__(tabix_file, reference, start, end)
        self.dict_chr2idx = dict_chr2idx
        # chr_prefix_operator can be 'add', 'remove' or None
        self.chr_prefix_operator = chr_prefix_operator
        
    def escape_chr(self, reference):
        if not self.chr_prefix_operator:
            return(reference)
        elif self.chr_prefix_operator == 'add':
            return 'chr' + reference
        elif self.chr_prefix_operator == 'del':
            return reference[3:]         
    def unescape_chr(self, reference):
        if not self.chr_prefix_operator:
            return(reference)
        elif self.chr_prefix_operator == 'del':
            return 'chr' + reference
        elif self.chr_prefix_operator == 'add':
            return reference[3:]                     
    def chr2idx(self, reference):
        return self.dict_chr2idx[reference]       
    def get_contigs(self):
        f = pysam.TabixFile(self.tabix_file, mode="r")# @UndefinedVariable
        return f.contigs 
    def __iter__(self) -> (Location, float):
        f = pysam.TabixFile(self.tabix_file, mode="r") # @UndefinedVariable
        try:
            for row in f.fetch(reference=self.escape_chr(self.reference), start=self.start - 1, end=self.end, parser=pysam.asTuple()): # @UndefinedVariable
                reference = self.unescape_chr(row[0])
                yield Location(self.chr2idx(reference), int(row[1]) + 1, int(row[2])), float(row[3])
        except:
            raise
        finally:
            f.close()

class LocationTabixBEDIterator(TabixIterator):
    """ Iterates over tabix-indexed BED files and returns location/(name,score,strand,color) pairs.
        Requires a tabix-indexed file with a numeric (float) value in the 4th column and a dict that converts chromosome names
        to indices. Chromosome 'chr' prefixes will be removed automatically if the chomosome names in the passed dict do not have them.
        The iterated part of the input file can be defined with reference, start and end. 
        
        up/dowstream slop can be configured.
        
        Example usage:
        it1 = LocationIterator('mydata.bedgraph.gz', dict_chr2idx, reference='chr1')
        it2 = LocationIterator('mydata.bedgraph.gz', dict_chr2idx, reference='chr2')
        for (loc1, value1), (loc2, value2) in zip(it1, it2):
            print("%s\t%f / %s\t%f" % (loc1, value1, loc2, value2))
            
    """
    def __init__(self, tabix_file, dict_chr2idx, reference=None, start=None, end=None, chr_prefix_operator=None, slop_left=0, slop_right=0):
        super().__init__(tabix_file, reference, start, end)
        self.dict_chr2idx = dict_chr2idx
        # chr_prefix_operator can be 'add', 'remove' or None
        self.chr_prefix_operator = chr_prefix_operator
        self.slop_left=slop_left
        self.slop_right=slop_right       
    def escape_chr(self, reference):
        if not self.chr_prefix_operator:
            return(reference)
        elif self.chr_prefix_operator == 'add':
            return 'chr' + reference
        elif self.chr_prefix_operator == 'del':
            return reference[3:]  
    def unescape_chr(self, reference):
        if not self.chr_prefix_operator:
            return(reference)
        elif self.chr_prefix_operator == 'del':
            return 'chr' + reference
        elif self.chr_prefix_operator == 'add':
            return reference[3:]                   
    def chr2idx(self, reference):
        return self.dict_chr2idx[reference]   
    def get_contigs(self):
        f = pysam.TabixFile(self.tabix_file, mode="r")# @UndefinedVariable
        return f.contigs 
    def __iter__(self) -> (Location, (str,str,str,str)):
        f = pysam.TabixFile(self.tabix_file, mode="r") # @UndefinedVariable
        try:
            for row in f.fetch(reference=self.escape_chr(self.reference), start=self.start - 1, end=self.end, parser=pysam.asTuple()): # @UndefinedVariable
                reference = self.unescape_chr(row[0])
                name = row[3] if len(row)>3 else None
                score=int(row[4]) if len(row)>4 else None
                strand=row[5] if len(row)>5 else None
                col=row[6] if len(row)>6 else None
                yield Location(self.chr2idx(reference), int(row[1]) + 1 - self.slop_left, int(row[2])+ self.slop_right) , ( name,score,strand,col )
        except:
            raise
        finally:
            f.close()
            
class LocationTabixPerPositionIterator(LocationTabixIterator):
    """ Iterates over tabix-indexed files and returns location/value pairs but splits intervals into single genomic locations.
        chr_prefix_operator defines how to convert input reference names to reference names found in the tabix file:
            1 will add a 'chr' prefix, -1 will remove a 'chr' prefix
        Works only for non-overlapping intervals!
    """

    def __init__(self, tabix_file, dict_chr2idx, reference, start, end, chr_prefix_operator=None, 
                 seq_col=0, start_col=1, end_col=2, is_zero_based=False,
                 fill=True, fill_value=None):
        super().__init__(tabix_file, dict_chr2idx, reference, start, end, chr_prefix_operator)
        self.fill = fill
        self.fill_value = fill_value
        self.seq_col=seq_col
        self.start_col=start_col
        self.end_col=end_col
        self.off = 1 if is_zero_based else 0
    def get_references(self):
        return pysam.TabixFile(self.tabix_file, mode="r").contigs # @UndefinedVariable
    def __iter__(self) -> (Location, tuple):
        f = pysam.TabixFile(self.tabix_file, mode="r") # @UndefinedVariable
        try:
            # data rows
            pos1 = self.start
            reference = self.reference
            chr_idx = self.chr2idx(reference)
            for row in f.fetch(reference=self.escape_chr(self.reference), 
                               start=self.start - 1, 
                               end=self.end, 
                               parser=pysam.asTuple()): # @UndefinedVariable
                if self.fill:
                    while pos1 < int(row[self.start_col])+self.off:
                        yield Location(chr_idx, pos1, pos1), self.fill_value
                        pos1 += 1
                if self.start_col == self.end_col: # position based
                    yield Location(chr_idx, pos1, pos1), row[self.end_col+1:]
                else:
                    for pos1 in range(int(row[self.start_col])+self.off, int(row[self.end_col])+self.off): # interval based
                        if pos1 < self.start:  # tabix intervals might start before configured coord
                            continue
                        if pos1 > self.end:
                            continue
                        yield Location(chr_idx, pos1, pos1), row[self.end_col+1:]
                pos1 += 1
            # fill until the end
            if self.fill and chr_idx is not None:
                while pos1 <= self.end:
                    yield Location(chr_idx, pos1, pos1), self.fill_value
                    pos1 += 1
        except:
            raise
        finally:
            f.close()


    
class LocationFastaIterator():
    """Iterates over a FASTA file by genomic location"""
    def __init__(self, fasta_file, dict_chr2idx, reference, start, end):
        self.fasta_file = fasta_file
        self.start = start
        self.end = end
        self.dict_chr2idx = dict_chr2idx
        self.reference = reference
        self.chr_idx = self.dict_chr2idx[self.reference]
        self.fasta = None
    def __iter__(self) -> (Location, str):    
        self.fasta = pysam.FastaFile(self.fasta_file).fetch(reference=self.reference, start=self.start - 1, end=self.end)
        self.pos1 = 1 if self.start is None else self.start
        for base in self.fasta:
            yield Location(self.chr_idx, self.pos1, self.pos1), base
            self.pos1 = self.pos1 + 1
    def close(self):
        if self.fasta:
            self.fasta.close()

            
class LocationPileupIterator():
    """ Iterates a BAM file per position. NOTE: not peekable as iterator will progress! """
    def __init__(self, bam_file, dict_chr2idx, reference, start, end, max_depth=100000, min_base_quality=0, flag_filter=BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY):
        self.bam_file = bam_file
        self.start = start
        self.end = end
        self.dict_chr2idx = dict_chr2idx
        self.reference = reference
        self.chr_idx = self.dict_chr2idx[self.reference]
        self.min_base_quality = min_base_quality
        self.samfile = pysam.AlignmentFile(self.bam_file, "rb")   # @UndefinedVariable
        self.flag_filter = flag_filter   
        self.max_depth = max_depth
    def __iter__(self) -> (Location, str):
        self.pit = self.samfile.pileup(self.reference,
                                   self.start - 1,
                                   self.end,
                                   flag_filter=self.flag_filter,
                                   truncate=True,
                                   mark_matches=False,
                                   mark_ends=False,
                                   add_indels=False,
                                   min_base_quality=self.min_base_quality,
                                   max_depth=self.max_depth)
        for pile in self.pit:
            ref_pos1 = pile.reference_pos + 1
            yield Location(self.chr_idx, ref_pos1, ref_pos1), pile

    def close(self):
        #print("closing pileup iterator %s@%s:%i-%i" %(self.bam_file, self.reference, self.start, self.end))
        self.samfile.close()

class BlockStrategy(Enum):
    LEFT = 1  # same start coord
    RIGHT = 2  # same end coord
    BOTH = 3  # same start and end
    OVERLAP = 4  # overlapping


class BlockLocationIterator():
    """ returns locations and lists of values that share the same location wrt. a given matching strategy (e.g., same start, same end, same coords, overlapping) """
    def __init__(self, it, strategy=BlockStrategy.LEFT):
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
                    mloc=Location.merge(mloc, l)
            elif self.strategy == BlockStrategy.RIGHT:
                while self.it.peek(None) and self.it.peek()[0].right_match(mloc):
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc=Location.merge(mloc, l)
            elif self.strategy == BlockStrategy.BOTH:
                while self.it.peek(None) and self.it.peek()[0] == mloc:
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc=Location.merge(mloc, l)
            elif self.strategy == BlockStrategy.OVERLAP:
                while self.it.peek(None) and self.it.peek()[0].overlaps(mloc):
                    l, v = next(self.it)
                    locations += [l]
                    values += [v]
                    mloc=Location.merge(mloc, l)
            yield mloc, (locations, values)
    def close(self):
        #print("close BlockLocationIterator")
        try:
            self.orgit.close()
        except AttributeError:
            pass

class OverlappingLocationPerPositionIterator():
    """ returns (overlapping) locations per position """
    def __init__(self, it):
        self.orgit=it
        self.it = peekable(BlockLocationIterator(it, strategy=BlockStrategy.LEFT))
    def __iter__(self) -> (Location, tuple):
        if not self.it:
            return
        buf = []
        chr_idx = None
        while self.it or len(buf) > 0:
            if chr_idx is None: 
                next_loc, _ = self.it.peek()
                chr_idx = next_loc.chr_idx
                pos1 = next_loc.start
                progress_to = next_loc.end
            else:
                if len(buf) > 0:
                    progress_to = max([l.end for l, _ in buf])
                    pos1 += 1
                else:
                    next_loc, _ = self.it.peek()
                    chr_idx = next_loc.chr_idx
                    pos1 = next_loc.start
                    progress_to = next_loc.end
                
            for pos1 in range(pos1, progress_to + 1):
                loc = Location(chr_idx, pos1, pos1)
                # add new intervals from next block
                while self.it and self.it.peek()[0].overlaps(loc):
                    _, iitems = next(self.it) 
                    for iloc, ivalue in zip(iitems[0], iitems[1]):
                        buf += [(iloc, ivalue)]
                # remove past intervals
                buf = [x for x in buf if loc.overlaps(x[0])]
                yield loc, buf
            # remove past intervals
            buf = [x for x in buf if x[0].end > loc.end] 
    def close(self):
        #print("close OverlappingLocationPerPositionIterator")
        try:
            self.orgit.close()
        except AttributeError:
            pass             


class ReadIterator():
    """ Iterates over a BAM alignment """
    def __init__(self, bam_file, dict_chr2idx, reference, start, end, flag_filter=BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY, max_span=None):
        self.bam_file = bam_file
        self.start = max(0, start) # avoid negative start coordinates!
        self.end = end
        self.dict_chr2idx = dict_chr2idx
        self.reference = reference
        self.chr_idx = self.dict_chr2idx[self.reference]
        self.samfile = pysam.AlignmentFile(self.bam_file, "rb")   # @UndefinedVariable
        self.flag_filter = flag_filter   
        self.region = Location(self.chr_idx, start, end)
        self.max_span = max_span
    def __iter__(self) -> (Location, str):
        for r in self.samfile.fetch(contig=self.reference, 
                                    start=self.start, 
                                    end=self.end,
                                    until_eof=True):
            # test flags
            if r.flag & self.flag_filter:
                #print("skip %s %i %i" % (r, r.flag, self.flag_filter) )
                continue
            rloc = Location(self.chr_idx, r.reference_start+1, r.reference_end)
            # test max_span and drop reads that span larger genomic regions
            if self.max_span and rloc.len() > self.max_span:
                continue
            # test location
            #if rloc.overlaps(self.region):
            yield rloc, r
    def close(self):
        #print("close ReadIterator")
        self.samfile.close()  
        
class PyrangeIterator():
    """ Iterates over a pyranges/pandas dataframe """
    def __init__(self, df, dict_chr2idx, feature):
        self.df = df
        self.dict_chr2idx = dict_chr2idx
        self.feature=feature
    def __iter__(self) -> (Location, str):
        for index, row in self.df.iterrows():
            chr_idx=self.dict_chr2idx[row['Chromosome']]
            start=row['Start']
            end=row['End']
            rloc = Location(chr_idx, start, end)
            name=row[self.feature]
            yield rloc, name
    def close(self):
        pass
    
class ReadMisMatchIterator():
    """ Iterates over a BAM alignment and returns reads and a list of mismatches (requires nM flag to be set) """
    def __init__(self, bam_file, fasta_file, dict_chr2idx, reference, start, end, flag_filter=BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY, max_span=None):
        self.bam_file = bam_file
        self.fasta_file = fasta_file
        self.start = max(0, start) # avoid negative start coordinates!
        self.end = end
        self.dict_chr2idx = dict_chr2idx
        self.reference = reference
        self.chr_idx = self.dict_chr2idx[self.reference]
        self.samfile = pysam.AlignmentFile(self.bam_file, "rb")   # @UndefinedVariable
        self.flag_filter = flag_filter   
        self.region = Location(self.chr_idx, start, end)
        self.max_span = max_span

    def __iter__(self):
        region_ref = pysam.FastaFile(self.fasta_file).fetch(
            reference=self.reference, 
            start=self.start - 1, 
            end=self.end)
        region_off = self.start - 1
        region_len = self.end - self.start
        for r in self.samfile.fetch(contig=self.reference, 
                                    start=self.start, 
                                    end=self.end,
                                    until_eof=True):
            # test flags
            if r.flag & self.flag_filter:
                #print("skip %s %i %i" % (r, r.flag, self.flag_filter) )
                continue
            rloc = Location(self.chr_idx, r.reference_start+1, r.reference_end)
            # test max_span and drop reads that span larger genomic regions
            if self.max_span and rloc.len() > self.max_span:
                continue
            if not rloc.overlaps(self.region):
                continue
                
            # tuples (x,y) where x is offset in r.query_alignment_sequence and y is offset in region_ref
            try:
                rel_pos = [(x,y-region_off) for (x,y) in r.get_aligned_pairs() if ((y is not None and x is not None) and (y-region_off>=0) and (y-region_off<region_len))]
                seq = [r.query_sequence[x] for (x,y) in rel_pos]
                ref = [region_ref[y] for (x,y) in rel_pos]
            except:
                print("ERROR in read %s" %r)
                print(rel_pos)
                print(r.query_sequence)
                print(r.get_aligned_pairs())
                raise
            # calculate mismatches (ref/alt)
            mm=[(r,a) for r,a in zip(ref, seq) if r != a]
            if len(seq)>0:
                yield rloc, (r, mm)
            
    def close(self):
        #print("close ReadIterator")
        self.samfile.close()   
        
class BlockedPerPositionIterator():
    """ Wraps the passed per-position iterators and synchronizes by genomic location. NOTE: all input iterators must be sorted by location
        Example:
        it = BlockedPerPositionIterator([it1, it2, it3])
        for loc,(v1,v2,v3) in it:
            print(loc,v1,v2,v3) 
    """
    def __init__(self, iterables):
        self.iterables = iterables
        self.iterators = [iter(it) for it in iterables]
        self.pos=SortedSet()
        self.next_values=[]
        for it in self.iterators:
            value = next(it, None)
            if value:
                self.next_values+=[value]
                self.pos.add(value[0])
            else:
                self.next_values+=[None]
    def __iter__(self) -> (Location, tuple):
        while self.pos:
            loc = self.pos.pop(0)
            ret=[]
            for idx, it in enumerate(self.iterators):
                if not self.next_values[idx]:
                    ret+=[None]
                else:
                    l = self.next_values[idx][0] 
                    ret+=[self.next_values[idx][1] if l==loc else None]
            yield loc, ret
            for idx, it in enumerate(self.iterators):
                while self.next_values[idx] and self.next_values[idx][0]<=loc:
                    self.next_values[idx]=next(it, None)
                    if self.next_values[idx]==None:
                        break
                if self.next_values[idx]:
                    self.pos.add(self.next_values[idx][0])
    def close(self):
        for it in self.iterables:
            try:
                it.close()
            except AttributeError:
                pass


class AnnotationOverlapIterator():
    """ gets one annotation iterator (sorted) and n BlockLocationIterators. Will yield the annotation location and 
    all overlapping intervals from the BlockLocationIterators """
    def __init__(self, anno_it, block_its):
        self.anno_it = anno_it
        self.iterators = [ peekable(it) for it in block_its ]
        self.buffer={}
        for idx, it in enumerate(self.iterators):
            self.buffer[idx] = []
    def __iter__(self) -> (Location, tuple):
        for aloc, a in self.anno_it:
            ret=[]
            for idx, it in enumerate(self.iterators):
                # remove non-overlapping items from beginning of list 
                self.buffer[idx] = [(loc, dat) for loc,dat in self.buffer[idx] if (loc.overlaps(aloc)) or (loc > aloc) ]
                # skip reads left of annotation
                while it.peek(None) and it.peek()[0] < aloc:
                    loc, dat = next(it)
                    if loc.overlaps(aloc):
                        self.buffer[idx] += [(loc, dat)]
                # add new overlapping items
                while it.peek(None) and it.peek()[0].overlaps(aloc):
                    loc, dat = next(it)
                    self.buffer[idx] += [(loc, dat)]
                # report current lists
                ret+=[self.buffer[idx]]
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
     

class SimulatedRead:

    def __init__(self, name, chromosome, start, end, splicing, spliceSites, bamRead):
        # Name of the parsed read
        self.name = name
        # Chromosome
        self.chromosome = chromosome
        # Absolute starting position of the read alignment
        self.start = start
        # Absolute end position of the read alignment
        self.end = end
        # Splicing status
        self.splicing = splicing
        # Splicing status
        self.spliceSites = spliceSites
        # Pysam read object
        self.bamRead = bamRead

    def hasSpliceSite(self, junction):
        if self.splicing:
            for spliceSite in self.spliceSites:
                res = junction.envelop(spliceSite[0], spliceSite[1])
                if (res) :
                    iv = res.pop()
                    begin, end, data = iv
                    if begin == spliceSite[0] and end == spliceSite[1]:
                        return True
        return False

    def __repr__(self):
        return "\t".join([self.name, self.chromosome, str(self.start), str(self.end), str(self.splicing), str(self.spliceSites)])


class SimulatedReadIterator:

    def __init__(self, readIterator, flagFilter=BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_SUPPLEMENTARY):
        self._readIterator = readIterator
        self._flagFilter = flagFilter

    def __iter__(self):
        return self

    def __next__(self):

        read = self._readIterator.__next__()

        while read.flag & self._flagFilter:
            read = self._readIterator.__next__()

        name = read.query_name

        chromosome = read.reference_name

        start = read.reference_start

        end = read.reference_end

        softclippedlist = read.get_reference_positions(full_length=True)

        offsetStart = 0
        offsetEnd = 0

        while (softclippedlist[offsetStart] != start) :
            offsetStart += 1

        softclippedlist.reverse()

        # reference_end points to one past the last aligned residue
        while (softclippedlist[offsetEnd] != (end - 1)) :
            offsetEnd += 1

        start -= offsetStart
        end += offsetEnd

        spliceSites = list()
        spliced = False

        refPositions = read.get_aligned_pairs(matches_only=False, with_seq=False)
        offset = 0

        for operation in read.cigartuples:
            if operation[0] == 3:
                spliced = True
                if not refPositions[offset + operation[1]][1]:
                    if not refPositions[offset + operation[1] - 1][1]:
                        if not refPositions[offset + operation[1] + 1][1]:
                            continue
                        else :
                            spliceSite = (
                            refPositions[offset - 1][1] + 1, refPositions[offset + operation[1] + 1][1] + 1)
                            spliceSites.append(spliceSite)
                    else :
                        spliceSite = (refPositions[offset - 1][1] + 1, refPositions[offset + operation[1] - 1][1] + 1)
                        spliceSites.append(spliceSite)
                else :
                    spliceSite = (refPositions[offset - 1][1] + 1, refPositions[offset + operation[1]][1] + 1)
                    spliceSites.append(spliceSite)

            offset += operation[1]

        simulatedRead = SimulatedRead(name, chromosome, start + 1, end, spliced, spliceSites, read)

        return simulatedRead
    
    
    
def get_chrom_dicts_from_bam(bam_file):
    """ returns dicts to translate chrom/chrom_idx and chrom lengths """
    samfile = pysam.AlignmentFile(bam_file, "rb")  # @UndefinedVariable
    chromosomes = samfile.header.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
    dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
    dict_chr2len = {c: samfile.header.get_reference_length(c) for c in chromosomes}
    samfile.close()
    return dict_chr2idx, dict_idx2chr, dict_chr2len