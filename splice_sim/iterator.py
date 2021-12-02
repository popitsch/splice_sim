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
        true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = name.split('_')
        trueSpliced = True if "N" in true_cigar else False
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
        simulatedRead = SimulatedRead(name, chromosome, start + 1, end, spliced, spliceSites, read, true_tid, trueSpliced)
        return simulatedRead
 