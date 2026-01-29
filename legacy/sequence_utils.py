from Bio.Seq import Seq


def set_aa_seq(svfusion): # translates the svfusion transcript to amino acids (to be fed into NetMHCpan)
    """
    :function: translate the nucelotide sequence in svfusion to amino acid sequence
    :param: svfusion: a SVFusion class
    :return: the corresponding amino acid sequence
    :NOTE: if there is a stop codon, the process will stop before translating all nucleotides
    """
    dna_seq = Seq(trim_to_3x(svfusion.nt_sequence))
    mrna_seq = dna_seq.transcribe()
    aa_seq = mrna_seq.translate(to_stop=True)
    return str(aa_seq)


def set_nt_seq(svfusion):
    """
    :function: given a svfusion, paste the nucleotides in cc_1, cc_2 and 3utr
    :param svfusion: a SVFusion class
    :return: the nucleotide sequence from start codon to the end of 3utr
    """
    return svfusion.nt_sequence_cds + svfusion.nt_sequence_3utr


def trim_to_3x(nt_sequence):
    """
    :function: trim a sequence to be divisible by 3
    :return: trimmed sequence
    """
    remainder = len(nt_sequence) % 3
    if remainder:
        return nt_sequence[: -remainder]
    else:
        return nt_sequence


def generate_neoepitopes(svfusion, window_range): 
    """
    :function: cut the WT and MUT protein sequence by window_range, then get the MUT specific peptides
    :param svfusion: a SVFusion class
    :param window_range: a list, specifying the range of window size, e.g. [8,9,10,11]
    :return: a list of neopeptides for svfusion
    """
    mut_peptides = []
    for window in window_range:
        mut_peptides = mut_peptides + cut_sequence(svfusion.aa_sequence, window)
    wt_peptides = []
    for window in window_range:
        wt_peptides = wt_peptides + cut_sequence(svfusion.cc_1.transcript.protein_sequence, window)
        wt_peptides = wt_peptides + cut_sequence(svfusion.cc_2.transcript.protein_sequence, window)
    return list(set(mut_peptides)-set(wt_peptides)) # so the whole time, we take out wildtype peptides from the mutant peptides for the given genes

def junction_indices(svfusion):
    if svfusion.cc_1.part == '5':
        left_nt = svfusion.cc_1.nt_sequence + svfusion.nt_sequence_ins
    else:
        left_nt = svfusion.cc_2.nt_sequence + svfusion.nt_sequence_ins
        
    junction_nt_index = len(left_nt)
    junction_aa_index = junction_nt_index // 3
    
    return junction_nt_index, junction_aa_index

def find_peptide_occurrences(aa_sequence, peptide):
    starts = []
    if not peptide:
        return starts
    i = aa_sequence.find(peptide)
    while i != -1:
        starts.append(i)
        i = aa_sequence.find(peptide, i + 1)
        
    return starts

def generate_neoepitopes_with_meta(svfusion, window_range):
    """Generate neoepitopes plus per-peptide metadata for breakpoint/junction analysis.

    This keeps the original NeoSV logic (mutant peptides not found in either WT partner protein),
    but also records where each neoepitope sits in the fused protein sequence and whether it spans
    the estimated junction.

    Returns
    -------
    (neoepitopes, meta)
        neoepitopes: list[str]
        meta: dict[str, dict] keyed by peptide sequence
    """
    neoepitopes = generate_neoepitopes(svfusion, window_range)
    j_nt, j_aa = junction_indices(svfusion)

    meta = {}
    for pep in neoepitopes:
        starts = find_peptide_occurrences(svfusion.aa_sequence, pep)
        # Choose the first occurrence for simple tabular output, but keep all for completeness.
        first_start = starts[0] if starts else None
        first_end = (first_start + len(pep) - 1) if first_start is not None else None
        spans = any(start < j_aa < (start + len(pep)) for start in starts) # does any occurrence span junction? returns T/F boolean
        meta[pep] = {
            'junction_nt': j_nt,
            'junction_aa': j_aa,
            'pep_starts_aa': starts,
            'num_pep_occurrences': len(starts),
            'pep_start_aa': first_start,
            'pep_end_aa': first_end,
            'spans_junction': spans,
        }
    return neoepitopes, meta
    

def cut_sequence(aa_sequence, window):
    """
    :param aa_sequence: amino acid sequence
    :param window: the size of window, window should be smaller than len(aa_sequence)
    :return: all possible slices with length = window in this sequence
    """
    if len(aa_sequence) < window:
        return []
    else:
        return [aa_sequence[i: i+window] for i in range(len(aa_sequence)-window+1)]
