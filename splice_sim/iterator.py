import pysam
from genomic_iterators.bam_iterators import BAM_FUNMAP,BAM_FSECONDARY,BAM_FQCFAIL,BAM_FDUP,BAM_SUPPLEMENTARY

class SimulatedRead:
    def __init__(self, name, chromosome, start, end, splicing, spliceSites, bamRead, trueTid, trueSpliced):
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
        # Pysam read object
        self.trueTid = trueTid
        # Pysam read object
        self.trueSpliced = trueSpliced
    def hasSpliceSite(self, intron_iv, max_diff=5):
        """ Checks overlap of passed intron interval and read splice-sites. 
            Overlap check is fuzzy to account for some small alignment errors.
            Total diff between endpoints must be <= max_diff
        """ 
        for ss in self.spliceSites:
            diff=abs(ss[0]-intron_iv.begin)+abs(ss[1]-intron_iv.end)
            if diff <= max_diff:
                return True
        return False
    def hasExactSpliceSite(self, junction):
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

def get_unaligned_blocks(r, min_len=10):
    """ get blocks of unaligned regions in genome coordinates. 
        Only blocks of a given minimum length are reported """
    last=None
    missing=[]
    for (x,y) in r.get_blocks():
        if last is None:
            last=(x,y)
        else:
            if last[1]==x: # special case, adjacent aligned blocks
                last=(last[0],y)
                continue
            missing+=[(last[1], x)]
        last=(x,y)
    return [(x,y) for x,y in missing if (y-x+1)>=min_len]

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
        true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = name.split('_')
        trueSpliced = "N" in true_cigar
        spliceSites=set([(x,y+1) for x,y in get_unaligned_blocks(read)])
        spliced = len(spliceSites)>0
        chromosome = read.reference_name
        start = read.reference_start
        end = read.reference_end
        simulatedRead = SimulatedRead(name, chromosome, start + 1, end, spliced, spliceSites, read, true_tid, trueSpliced)
        return simulatedRead
 