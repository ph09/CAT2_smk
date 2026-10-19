import collections
import re
import pandas as pd
from . import fileOps
from . import transcripts
from . import misc


reserved_keys = ['gene_biotype',
                 'transcript_biotype',
                 'gene_type',
                 'transcript_type',
                 'gene_name',
                 'gene_id',
                 'transcript_id',
                 'transcript_name',
                 'ID',
                 'Name',
                 'Parent']


# Reference gene IDs of the current run (versionless). Set by the consensus
# entry points from the reference annotation; inherited by forked workers.
_REF_GENE_IDS = None
# Fallback when unset: any Ensembl species gene ID (ENSG, ENSMUSG, ENSGALG...).
ENSEMBL_GENE_RE = re.compile(r"^ENS[A-Z]*G\d+")


def set_reference_gene_ids(gene_ids):
    global _REF_GENE_IDS
    _REF_GENE_IDS = {str(g).split('.')[0] for g in gene_ids}


def is_reference_gene(gene_id):
    """True if gene_id is a gene of the reference annotation (any naming scheme)."""
    g = str(gene_id)
    if _REF_GENE_IDS is not None:
        return g.split('.')[0] in _REF_GENE_IDS
    return bool(ENSEMBL_GENE_RE.match(g))


def readthrough_genes(ref_df):
    """Reference readthrough gene IDs.

    GENCODE (and converted RefSeq) tag them readthrough_gene / readthrough_transcript.
    Native Ensembl annotations of other species carry no such tag: there a gene
    counts when its description says readthrough, or, only if nothing was found
    that way, when its name is NAME1-NAME2 with both parts reference genes.
    """
    if 'ExtraTags' not in ref_df.columns:
        return set()
    tags = ref_df['ExtraTags'].fillna('')
    mask = tags.str.contains('readthrough_gene') | tags.str.contains('readthrough_transcript')
    found = set(ref_df[mask]['GeneId'].unique())
    if found:
        return found
    found = set(ref_df[tags.str.contains(r'description=[^;]*readthrough', case=False)]['GeneId'].unique())
    if found or 'GeneName' not in ref_df.columns:
        return found
    names = set(ref_df['GeneName'].dropna().astype(str).str.upper())
    for gid, name in zip(ref_df['GeneId'], ref_df['GeneName'].astype(str).str.upper()):
        parts = name.split('-')
        if len(parts) == 2 and all(p in names and p != name for p in parts):
            found.add(gid)
    return found


def check_gff3_for_disjoint_gene_ids(tx_dict, by_chrom=True):
    """Reports on gene identifiers that are disjoint in genomic space"""
    # group by name2
    tx_name2_map = transcripts.group_transcripts_by_name2(tx_dict.values())
    bad_genes = []
    for gene_id, tx_group in tx_name2_map.items():
        if by_chrom == False:
            locations = sorted((tx.interval for tx in tx_group))
            for i in range(0, len(locations) - 1):
                if not locations[i].overlap(locations[i+1]):
                    bad_genes.append(gene_id)
        else:
            if len({tx.chromosome for tx in tx_group}) != 1:
                bad_genes.append(gene_id)
    return bad_genes


def parse_gff3(annotation_attrs, annotation_gp, is_external_reference=False):
    def parse_attrs(attrs):
        r = collections.defaultdict(dict)
        for tx_id, key, value in fileOps.iter_lines(attrs):
            r[tx_id][key] = value
        return r

    attrs_dict = parse_attrs(annotation_attrs)
    tx_dict = transcripts.get_gene_pred_dict(annotation_gp)
    disjoint_genes = check_gff3_for_disjoint_gene_ids(tx_dict, by_chrom=True)
    if disjoint_genes:
        raise Exception(f"Found chromosome disjoint gene_ids {', '.join(disjoint_genes)} in the GFF3. "
                        f"This is not allowed. Use the script programs/fix_chrom_disjoint_genes to "
                        f"try and resolve this.")
    tx_name_map = {x: y.name2 for x, y in tx_dict.items()}
    results = []
    for tx_id, gene_id in tx_name_map.items():
        d = attrs_dict[tx_id]
        gene_biotype = d.get('gene_biotype', d.get('gene_type', None))
        if gene_biotype is None:
            raise Exception("Did not find a gene biotype or gene type for {} (attrs={})".format(gene_id, d))
        tx_biotype = d.get('transcript_biotype', d.get('transcript_type', None))
        if tx_biotype is None:
            raise Exception("Did not find a transcript biotype or type for {} (attrs={})".format(tx_id, d))
        gene_name = d['gene_name']
        gene_id = d['gene_id']
        tx_id = d['transcript_id']
        tx_name = d['transcript_name']
        extra_tags = ';'.join(['{}={}'.format(x, y.replace(';', '%3B').replace('=', '%3D'))
                               for x, y in d.items() if x not in reserved_keys])
        if len(extra_tags) > 0:
            try:
                misc.parse_gff_attr_line(extra_tags)
            except:
                raise Exception(f'Error parsing extra tags in input GFF3 {extra_tags}')
        if is_external_reference == True:
            # hack to fix names
            gene_id = f'exRef-{gene_id}'
            tx_id = f'exRef-{tx_id}'
        results.append([gene_id, tx_id, tx_name, gene_name, gene_biotype, tx_biotype, extra_tags])
    df = pd.DataFrame(results, columns=['GeneId', 'TranscriptId', 'TranscriptName', 'GeneName',
                                        'GeneBiotype', 'TranscriptBiotype', 'ExtraTags'])
    df = df.set_index('TranscriptId')
    return df


def convert_gff3_cmd(annotation_attrs, annotation):
    cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
            '-refseqHacks',
           '-attrsOut={}'.format(annotation_attrs), annotation, '/dev/stdout']
    return cmd
