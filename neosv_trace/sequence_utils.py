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

def compress_genomic_blocks(origins, svfusion):
    """Compress per-nt origin dicts into human-readable genomic blocks.

    Returns (origin_type, summary_str, feature_str)

    - If any INS nt is present -> ('INS','INS','INS')
    - If empty -> ('EMPTY','','')
    - Otherwise, compress consecutive genomic positions into blocks.
      Blocks are tracked separately for cc_side=1 vs cc_side=2 so that exon numbering
      is computed against the correct transcript.
    """
    if not origins:
        return 'EMPTY', '', ''
    if any(o.get('type') == 'INS' for o in origins):
        return 'INS', 'INS', 'INS'

    blocks = []
    cur = None
    for o in origins:
        chrom = o.get('chrom')
        pos = int(o.get('pos'))
        feat = o.get('feature', '')
        cc_side = int(o.get('cc_side', 0) or 0)

        if cur is None:
            cur = {'chrom': chrom, 'start': pos, 'end': pos, 'feat': feat, 'cc_side': cc_side}
            continue

        adjacent = (
            chrom == cur['chrom']
            and cc_side == cur['cc_side']
            and feat == cur['feat']
            and abs(pos - cur['end']) == 1
        )
        if adjacent:
            cur['end'] = pos
        else:
            blocks.append(cur)
            cur = {'chrom': chrom, 'start': pos, 'end': pos, 'feat': feat, 'cc_side': cc_side}

    if cur is not None:
        blocks.append(cur)

    def _tx_for_side(side: int):
        if side == 1:
            return getattr(svfusion.cc_1, "transcript", None)
        if side == 2:
            return getattr(svfusion.cc_2, "transcript", None)
        return None

    def fmt_block(b):
        tx = _tx_for_side(int(b.get('cc_side', 0) or 0))
        s, e = int(b['start']), int(b['end'])
        # Print in ascending genomic coordinates for readability
        lo, hi = (s, e) if s <= e else (e, s)

        start_exon_num, start_exon_id = exon_info_at_pos(tx, lo)
        end_exon_num, end_exon_id = exon_info_at_pos(tx, hi)

        exon_part = ""
        if start_exon_num is not None and end_exon_num is not None:
            if start_exon_num == end_exon_num:
                exon_part = f",E{start_exon_num}"
            else:
                exon_part = f",E{start_exon_num}-E{end_exon_num}"
        # Keep it compact; exon_id often None depending on pyensembl version
        return f"{b['chrom']}:{lo}-{hi}({b['feat']}{exon_part})"

    summary = ';'.join(fmt_block(b) for b in blocks)
    feature = ';'.join(sorted(set(b['feat'] for b in blocks if b.get('feat'))))
    return 'GENOMIC', summary, feature


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
        
        # Build per-peptide genomic origin summary (based on CDS nucleotides only)
        if not hasattr(svfusion, 'cds_origin') or svfusion.cds_origin is None:
            try:
                svfusion.build_cds_origin()
            except Exception:
                svfusion.cds_origin = []

        pep_origin_type = ''
        pep_origin = ''
        pep_feature = ''
        if first_start is not None:
            nt0 = first_start * 3
            nt1 = nt0 + len(pep) * 3  # exclusive
            origins = svfusion.cds_origin[nt0:nt1] if svfusion.cds_origin else []
            pep_origin_type, pep_origin, pep_feature = compress_genomic_blocks(origins, svfusion)

        
        meta[pep] = {
            'junction_nt': j_nt,
            'junction_aa': j_aa,
            'pep_starts_aa': starts,
            'num_pep_occurrences': len(starts),
            'pep_start_aa': first_start,
            'pep_end_aa': first_end,
            'spans_junction': spans,
            'pep_origin_type': pep_origin_type,
            'pep_origin': pep_origin,
            'pep_feature': pep_feature,
        }
    return neoepitopes, meta

def compute_peptide_meta(svfusion, pep: str) -> dict:
    """Compute metadata for a single peptide *on demand*.

    This avoids building large per-nucleotide origin vectors for fusions/peptides
    that never pass MHC filtering.
    """
    j_nt, j_aa = junction_indices(svfusion)
    transcript = svfusion.cc_1.transcript if svfusion.cc_1.part == '5' else svfusion.cc_2.transcript #FIXME this is a bit hacky, but we need the transcript to get exon info for the peptide origin summary. We could consider storing the transcript directly in the SVFusion CC objects in the future to avoid this.


    starts = find_peptide_occurrences(svfusion.aa_sequence, pep)
    first_start = starts[0] if starts else None
    first_end = (first_start + len(pep) - 1) if first_start is not None else None
    spans = any(start < j_aa < (start + len(pep)) for start in starts)

    pep_origin_type = ''
    pep_origin = ''
    pep_feature = ''
    if first_start is not None:
        # Only now build the per-nt CDS origin vector (can be large)
        if not hasattr(svfusion, 'cds_origin') or svfusion.cds_origin is None:
            try:
                svfusion.build_cds_origin()
            except Exception:
                svfusion.cds_origin = []

        nt0 = first_start * 3
        nt1 = nt0 + len(pep) * 3  # exclusive
        origins = svfusion.cds_origin[nt0:nt1] if svfusion.cds_origin else []
        pep_origin_type, pep_origin, pep_feature = compress_genomic_blocks(origins, svfusion)

    return {
        'junction_nt': j_nt,
        'junction_aa': j_aa,
        'pep_starts_aa': starts,
        'num_pep_occurrences': len(starts),
        'pep_start_aa': first_start,
        'pep_end_aa': first_end,
        'spans_junction': spans,
        'pep_origin_type': pep_origin_type,
        'pep_origin': pep_origin,
        'pep_feature': pep_feature,
    }
    

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


def exon_info_at_pos(transcript, pos: int):
    """Return (exon_number_1based, exon_id) for a genomic position in a transcript.

    Exon number is in transcript 5'->3' order:
      + strand: increasing genomic coordinate
      - strand: decreasing genomic coordinate

    Returns (None, None) if transcript is None or pos is not in any exon.
    """
    if transcript is None:
        return (None, None)

    exons = list(getattr(transcript, "exons", []) or [])
    # sort exons in transcript order (5' -> 3')
    if getattr(transcript, "strand", "+") == "+":
        exons.sort(key=lambda e: (e.start, e.end))
    else:
        exons.sort(key=lambda e: (e.start, e.end), reverse=True)

    for i, exon in enumerate(exons, start=1):
        if exon.start <= pos <= exon.end:
            exon_id = getattr(exon, "id", None)
            return (i, exon_id)

    return (None, None)


def exon_info_for_block(transcript, start: int, end: int):
    """
    Docstring for exon_info_for_block
    
    :param transcript: Description
    :param start: Description
    :type start: int
    :param end: Description
    :type end: int
    """

    exons = list(transcript.exons)
    if transcript.strand == "+": exons.sort(key=lambda e: (e.start, e.end))
    else: exons.sort(key=lambda e: (e.start, e.end), reverse=True)

    hits = []

    for i, exon in enumerate(exons, start=1):
        if not (end < exon.start or start > exon.end):
            hits.append({"exon_number": i, "exon_id": getattr(exon, "id", None), "start": exon.start, "end": exon.end}) 
            
    return hits