"""
Perform name conversions on transMap/AugustusTMR transcripts.
"""

import re

# Global mapping for alignment ID to source type
alignment_source_map = {}


def set_alignment_source_map(source_map):
    """Set the global alignment source mapping"""
    global alignment_source_map
    alignment_source_map = source_map


def remove_alignment_number(aln_id, aln_re=re.compile("-[0-9]+$"), txTM_copy_re=re.compile("_[0-9]+$")):
    """
    If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    Also handles txTM copy suffixes like _7 in ENST00000567343.1_7
    :param aln_id: name string
    :param aln_re: compiled regular expression for transMap alignment numbers
    :param txTM_copy_re: compiled regular expression for txTM copy numbers
    :return: string
    """
    # Remove both types of suffixes, repeating until no more matches
    # This handles cases like ENST00000567343.1-2_7
    prev_aln_id = None
    while prev_aln_id != aln_id:
        prev_aln_id = aln_id
        # Try to remove txTM-style copy numbers (_d)
        aln_id = txTM_copy_re.split(aln_id)[0]
        # Try to remove transMap-style alignment numbers (-d)
        aln_id = aln_re.split(aln_id)[0]
    return aln_id


def remove_augustus_alignment_number(aln_id, aug_re=re.compile("^(aug(TM|TMR|MP|PB)|txTM|strg|exRef)-")):
    """
    removes the alignment numbers prepended by AugustusTM/AugustusTMR/AugustusMP/txTM/strg/exRef
    Format: aug(TM|TMR|MP|PB)-ENSMUST00000169901.2-1
    Format: txTM-ENST00000567343.1_7
    Format: strg-1234.t1
    Format: exRef-ENST00000567343.1
    :param aln_id: name string
    :param aug_re: compiled regular expression
    :return: string
    """
    return aug_re.split(aln_id)[-1]


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    :param aln_id: name string
    :return: string
    """
    return remove_alignment_number(remove_augustus_alignment_number(aln_id))


def alignment_id_to_ref_transcript_id(aln_id):
    """
    Normalize an alignment ID to the reference GenePred transcript name.

    Applies ``strip_alignment_numbers`` (mode prefixes, transMap ``-N``, txTM ``_N``)
    and strips miniprot/augMP ``rna-`` accession prefixes so IDs like
    ``augMP-rna-XM_024988618.2`` map to ``XM_024988618.2``.
    """
    ref_id = strip_alignment_numbers(aln_id)
    if ref_id.startswith('rna-'):
        ref_id = ref_id[4:]
    return ref_id


def strip_alignment_numbers_preserve_txTM_copies(aln_id):
    """
    Strip alignment numbers but preserve txTM copy numbers for gene families.
    This is used in consensus to preserve gene family expansion identified by txTM.
    :param aln_id: name string
    :return: string
    """
    # Remove Augustus prefix but preserve txTM copy numbers
    aln_id = remove_augustus_alignment_number(aln_id)
    
    # Only remove transMap-style alignment numbers (-d), not txTM copy numbers (_d)
    aln_re = re.compile("-[0-9]+$")
    return aln_re.split(aln_id)[0]


def aln_id_is_augustus(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is an Augustus transcript
    :param aln_id: name string
    :return: boolean
    """
    return True if remove_augustus_alignment_number(aln_id) != aln_id else False


def aln_id_is_transmap(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is a transMap transcript
    Includes both regular transMap alignments (ending with -d) and chained projections (ending with _cp)
    :param aln_id: name string
    :return: boolean
    """
    # Check for chained projection IDs (ending with _cp)
    if aln_id.endswith('_cp'):
        return True
    # Check for regular transMap alignment IDs (ending with -d)
    return True if remove_augustus_alignment_number(aln_id) == aln_id and remove_alignment_number(aln_id) != aln_id else False


def aln_id_is_augustus_tm(aln_id):
    return aln_id.startswith('augTM-')


def aln_id_is_augustus_tmr(aln_id):
    return aln_id.startswith('augTMR-')


def aln_id_is_augustus_mp(aln_id):
    return aln_id.startswith('augMP-')


def aln_id_is_txTM(aln_id):
    return aln_id.startswith('txTM-')

def aln_id_is_strg(aln_id):
    return aln_id.startswith('strg-')

def aln_id_is_pb(aln_id):
    return aln_id.startswith('augPB-')



def aln_id_is_exref(aln_id):
    return aln_id.startswith('exRef-')


def aln_id_is_denovo(aln_id):
    return aln_id_is_pb(aln_id) or aln_id_is_strg(aln_id)


def alignment_type(aln_id):
    """returns what type of alignment this ID is"""
    if aln_id_is_augustus_tmr(aln_id):
        return 'augTMR'
    elif aln_id_is_augustus_mp(aln_id):
        return 'augMP'
    elif aln_id_is_augustus_tm(aln_id):
        return 'augTM'
    elif aln_id_is_pb(aln_id):
        return 'augPB'
    elif aln_id_is_exref(aln_id):
        return 'exRef'
    elif aln_id_is_txTM(aln_id):
        return 'txTM'
    elif aln_id_is_strg(aln_id):
        return 'strg'
    elif aln_id_is_transmap(aln_id):
        return 'transMap'
    else:
        # Try to use the source mapping if available
        global alignment_source_map
        if aln_id in alignment_source_map:
            return alignment_source_map[aln_id]
        
        # For plain transcript IDs that don't have prefixes, we can't determine the type
        # from the ID alone. This could be transMap or txTM data without prefixes.
        # Return a generic identifier that won't cause the join to fail.
        return 'unknown'
