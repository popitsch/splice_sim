"""
    @author: @author: niko.popitsch@univie.ac.at

    Utility methods for splice_sim
"""
import os
import re
import unicodedata

import pysam

# Table of reverse complement bases
COMP_TABLE = {
    "A": 'T', "C": 'G', "T": 'A', "G": 'C'
}

rcmap = bytes.maketrans(b'ATCGatcgNn', b'TAGCTAGCNN')


def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    return seq[::-1].translate(rcmap)


# Colors for read classifications
colors = {"exonexonfalse": "228,26,28",
          "exonexontrue": "55,126,184",
          "exonintron": "77,175,74",
          "garbage": "255,127,0",
          "intron": "200,200,200"
          }


def parse_info(info):
    """ parse GFF3 info section """
    return {k: v for k, v in [a.split('=') for a in info.split(';') if '=' in a]}


def pad_n(seq, minlen):
    """ Add up-/downstream padding with N's to ensure a given minimum length of the passed sequence """
    ret = seq
    if (len(ret) < minlen):
        pad0 = "N" * int((minlen - len(ret)) / 2)
        pad1 = "N" * int(minlen - (len(ret) + len(pad0)))
        ret = pad0 + ret + pad1
    return (ret)


def bgzip_and_tabix(in_file, out_file=None, create_index=True, del_uncompressed=True, seq_col=0, start_col=1, end_col=1,
                    line_skip=0, zerobased=False):
    """ Will BGZIP the passed file and creates a tabix index with the given params if create_index is True"""
    if out_file is None:
        out_file = in_file + '.gz'
    pysam.tabix_compress(in_file, out_file, force=True)  # @UndefinedVariable
    if create_index:
        pysam.tabix_index(out_file, force=True, seq_col=seq_col, start_col=start_col, end_col=end_col, meta_char='#',
                          line_skip=line_skip, zerobased=zerobased)  # @UndefinedVariable
    if del_uncompressed:
        os.remove(in_file)


def localize_config(config):
    """ Recursively iterates json and adds /Volumes prefix to /group paths"""
    for k, v in config.items():
        if isinstance(v, str) and v.startswith('/groups'):
            config[k] = '/Volumes' + config[k]
        elif isinstance(v, dict):
            config[k] = localize_config(v)
        elif isinstance(v, list):
            config[k] = ['/Volumes' + x if str(x).startswith('/groups') else x for x in v]
    return config


def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-')


def merge_bam_files(out_file, bam_files, sort_output=False, del_in_files=False):
    """ merge multiple BAM files and sort + index results """
    if bam_files is None or len(bam_files) == 0:
        print("no input BAM file provided")
        return None
    samfile = pysam.AlignmentFile(bam_files[0], "rb")  # @UndefinedVariable
    out = pysam.AlignmentFile(out_file + '.unsorted.bam', "wb", template=samfile)  # @UndefinedVariable
    for f in bam_files:
        samfile = None
        try:
            samfile = pysam.AlignmentFile(f, "rb")  # @UndefinedVariable
            for read in samfile.fetch(until_eof=True):
                out.write(read)
        except Exception as e:
            print("error opening bam %s: %s" % (f, e))
        finally:
            if samfile:
                samfile.close()
    out.close()
    if sort_output:
        try:
            pysam.sort("-o", out_file, out_file + '.unsorted.bam')  # @UndefinedVariable
            os.remove(out_file + '.unsorted.bam')
            if del_in_files:
                for f in bam_files + [b + '.bai' for b in bam_files]:
                    if os.path.exists(f):
                        os.remove(f)
        except Exception as e:
            print("error sorting this bam: %s" % e)
    else:
        os.rename(out_file + '.unsorted.bam', out_file)
        if del_in_files:
            for f in bam_files + [b + '.bai' for b in bam_files]:
                os.remove(f)
    # index
    try:
        pysam.index(out_file)  # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    return out_file


def flatten(t):
    """ Flatten nested list """
    return [i for s in t for i in s]
