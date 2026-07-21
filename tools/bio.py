"""
Basic biology related functions
"""
import string
import os
import pysam
from .fileOps import opengz


class _FastaRecord:
    """Proxy for a single chromosome that supports integer indexing and slicing."""

    __slots__ = ('_fasta', '_chrom', '_upper', '_length')

    def __init__(self, fasta, chrom, upper, length):
        self._fasta = fasta
        self._chrom = chrom
        self._upper = upper
        self._length = length

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else self._length
            seq = self._fasta.fetch(self._chrom, start, stop)
        else:
            seq = self._fasta.fetch(self._chrom, key, key + 1)
        return seq.upper() if self._upper else seq

    def __len__(self):
        return self._length

    def __str__(self):
        seq = self._fasta.fetch(self._chrom)
        return seq.upper() if self._upper else seq


class FastaDict:
    """Dict-like wrapper around ``pysam.FastaFile``.

    Returns ``_FastaRecord`` proxies that support integer indexing and
    Python slicing of chromosome sequences without loading the whole FASTA
    into memory.
    """

    def __init__(self, file_path, upper=False):
        self._fasta = pysam.FastaFile(file_path)
        self._upper = upper

    def __getitem__(self, chrom):
        return _FastaRecord(self._fasta, chrom, self._upper,
                            self._fasta.get_reference_length(chrom))

    def __contains__(self, chrom):
        return chrom in self._fasta.references

    def __iter__(self):
        return iter(self._fasta.references)

    def __len__(self):
        return self._fasta.nreferences

    def keys(self):
        return list(self._fasta.references)

    def close(self):
        self._fasta.close()


def write_fasta(path_or_handle, name, seq, chunk_size=100, validate=None):
    """Writes out fasta file. if path ends in gz, will be gzipped.
    """
    if isinstance(path_or_handle, str):
        fh = opengz(path_or_handle, 'w')
    else:
        fh = path_or_handle
    if validate == 'DNA':
        valid_chars = set('ACGTUYSWKMBDHVNacgtuyswkmbdhvn.-*')
    elif validate == 'protein':
        valid_chars = set('ABCDEFGHIKLMPQSRTVWXYZUabcdefghiklmpqsrtvwxyzuNn.-*')
    else:
        valid_chars = set()
    try:
        assert any([isinstance(seq, str), isinstance(seq, str)])
    except AssertionError:
        raise RuntimeError("Sequence is not unicode or string")
    if validate is not None:
        try:
            assert all(x in valid_chars for x in seq)
        except AssertionError:
            bad_chars = {x for x in seq if x not in valid_chars}
            raise RuntimeError("Invalid FASTA character(s) seen in fasta sequence: {}".format(bad_chars))
    fh.write(">%s\n" % name)
    for i in range(0, len(seq), chunk_size):
        fh.write("%s\n" % seq[i:i+chunk_size])
    if isinstance(path_or_handle, str):
        fh.close()


def complement(seq, comp=str.maketrans("ATGCatgc", "TACGtacg")):
    """
    given a sequence, return the complement.
    """
    return str(seq).translate(comp)


def reverse_complement(seq):
    """
    Given a sequence, return the reverse complement.
    """
    return complement(seq)[::-1]


_codon_table = {
    'ATG': 'M',
    'TAA': '*', 'TAG': '*', 'TGA': '*', 'TAR': '*', 'TRA': '*',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCN': 'A',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R',
    'AGG': 'R', 'CGN': 'R', 'MGR': 'R',
    'AAT': 'N', 'AAC': 'N', 'AAY': 'N',
    'GAT': 'D', 'GAC': 'D', 'GAY': 'D',
    'TGT': 'C', 'TGC': 'C', 'TGY': 'C',
    'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q',
    'GAA': 'E', 'GAG': 'E', 'GAR': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
    'CAT': 'H', 'CAC': 'H', 'CAY': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATH': 'I',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
    'CTG': 'L', 'YTR': 'L', 'CTN': 'L',
    'AAA': 'K', 'AAG': 'K', 'AAR': 'K',
    'TTT': 'F', 'TTC': 'F', 'TTY': 'F',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S',
    'AGC': 'S', 'TCN': 'S', 'AGY': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACN': 'T',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTN': 'V',
    '': ''
    }


def codon_to_amino_acid(c):
    """
    Given a codon C, return an amino acid or ??? if codon unrecognized.
    Codons could be unrecognized due to ambiguity in IUPAC characters.
    """
    assert len(c) == 3, c
    if c is None:
        return None
    if c in _codon_table:
        return _codon_table[c]
    return 'X'


def translate_sequence(sequence):
    """
    Translates a given DNA sequence to single-letter amino acid
    space. If the sequence is not a multiple of 3 and is not a unique degenerate codon it will be truncated silently.
    """
    result = []
    sequence = sequence.upper()
    i = 0
    for i in range(0, len(sequence) - len(sequence) % 3, 3):
        result.append(codon_to_amino_acid(sequence[i: i + 3]))
    if len(sequence) % 3 == 2:
        c = codon_to_amino_acid(sequence[i + 3:] + 'N')
        if c != 'X':
            result.append(c)
    return ''.join(result)


def read_codons(seq, offset=0, skip_last=True):
    """
    Provides an iterator that reads through a sequence one codon at a time.
    """
    l = len(seq)
    if skip_last:
        l -= 3
    for i in range(offset,  l - l % 3, 3):
        yield seq[i:i + 3]


def read_codons_with_position(seq, offset=0, skip_last=True):
    """
    Provides an iterator that reads through a sequence one codon at a time,
    returning both the codon and the start position in the sequence.
    """
    l = len(seq)
    if skip_last:
        l -= 3
    for i in range(offset, l - l % 3, 3):
        yield i, seq[i:i + 3]


def get_sequence_dict(file_path, upper=True):
    """
    Returns a dict-like object of fasta records backed by pysam.
    Requires a .fai index (created by ``samtools faidx``).
    If *upper* is True, all bases will be uppercased on access.
    """
    assert os.path.exists(file_path), ('Error: FASTA file {} does not exist'.format(file_path))
    fai_path = file_path + ".fai"
    assert os.path.exists(fai_path), (
        "Error: .fai index does not exist for {}. "
        "Run 'samtools faidx {}' first.".format(file_path, file_path))
    return FastaDict(file_path, upper=upper)
