import re


def write_annot(filepath, sveffects, prefix):
    """
    :param filepath: outfile for annotation
    :param sveffects: a list of sveffect classes
    :return: None
    """
    header = '\t'.join(['prefix', 'sv_id', 'chrom1', 'pos1', 'function1', 'gene1', 'transcript_id1', 'strand1', 'transcript_retain1',
                        'chrom2', 'pos2', 'function2', 'gene2', 'transcript_id2', 'strand2', 'transcript_retain2',
                        'svpattern', 'svtype', 'fusion'])
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for sveffect in sveffects:
            sv_id = getattr(sveffect.sv, "sv_id", None)
            sv_id = "None" if sv_id is None else str(sv_id)
            row = [prefix, sv_id] + sveffect.output()
            f.write('\t'.join('' if x is None else str(x) for x in row) + '\n')
            
def write_fusion(filepath, svfusions, dict_neo, prefix):
    """
    Write fusion-derived neoantigens with provenance columns.

    dict_neo:
        {neoepitope: {allele: (affinity, BA_rank, EL_rank, PASS/FILTER)}}
    """
    with open(filepath, 'w') as f:
        header_cols = [
            'prefix',
            # Base SV/fusion columns (match SVFusion.output())
            'sv_id',
            'chrom1', 'pos1', 'gene1', 'transcript_id1',
            'chrom2', 'pos2', 'gene2', 'transcript_id2',
            'svpattern', 'svtype', 'frameshift',
            # Provenance / localization columns
            'junction_nt', 'junction_aa',
            'num_pep_occurrences', 'pep_start_aa', 'pep_end_aa', 'spans_junction',
            'pep_origin_type', 'pep_origin', 'pep_feature',
            # NetMHC columns
            'neoantigen', 'allele', 'affinity', 'BA_rank', 'EL_rank'
        ]
        f.write('\t'.join(header_cols) + '\n')

        seen = set()


        # Compute peptide metadata lazily to avoid storing per-nucleotide origin
        # vectors for peptides that never pass MHC filtering.
        from .sequence_utils import compute_peptide_meta

        for svfusion in svfusions:
            base = svfusion.output()  # includes sv_id already as col 1

            for neoepitope in svfusion.neoepitopes:
                if neoepitope not in dict_neo:
                    continue

                # Only compute metadata if at least one allele passes MHC filtering.
                any_pass = any(v[3] == 'PASS' for v in dict_neo[neoepitope].values())
                if not any_pass:
                    continue

                # Lazily compute/cache meta now that we know we'll write at least one row.
                if not hasattr(svfusion, 'neoepitope_meta') or svfusion.neoepitope_meta is None:
                    svfusion.neoepitope_meta = {}
                if neoepitope not in svfusion.neoepitope_meta:
                    svfusion.neoepitope_meta[neoepitope] = compute_peptide_meta(svfusion, neoepitope)
                meta = svfusion.neoepitope_meta.get(neoepitope, {}) or {}

                for allele, vals in dict_neo[neoepitope].items():
                    if vals[3] != 'PASS':
                        continue

                    affinity = vals[0]
                    ba_rank = vals[1]
                    el_rank = vals[2]

                    row = (
                        [prefix]
                        + base
                        + [
                            meta.get('junction_nt'),
                            meta.get('junction_aa'),
                            meta.get('num_pep_occurrences'),
                            meta.get('pep_start_aa'),
                            meta.get('pep_end_aa'),
                            meta.get('spans_junction', False),
                            meta.get('pep_origin_type', ''),
                            meta.get('pep_origin', ''),
                            meta.get('pep_feature', ''),
                            neoepitope,
                            allele,
                            affinity,
                            ba_rank,
                            el_rank,
                        ]
                    )

                    # event = sv_event_key(svfusion.sv)
                    dedup_key = (
                        canonical_sv_id(getattr(svfusion.sv, "sv_id", None)) or getattr(svfusion.sv, "sv_id", None),
                        neoepitope,        # peptide string
                        allele,            # HLA allele string
                    )

                    if dedup_key in seen:
                        continue
                    seen.add(dedup_key)

                    # stringify at the edge
                    f.write('\t'.join('' if x is None else str(x) for x in row) + '\n')

            # Free potentially-large per-nucleotide origin arrays once this SVFusion is done.
            if hasattr(svfusion, 'cds_origin'):
                svfusion.cds_origin = None

def canonical_sv_id(sv_id: str | None) -> str | None:
    if not sv_id:
        return None
    return re.sub(r"_[12]$", "", sv_id)

def sv_event_key(sv) -> tuple:
    """
    Canonical identity for an SV event for output dedup.
    """
    sv_id = canonical_sv_id(getattr(sv, "sv_id", None)) or getattr(sv, "sv_id", None)

    a = getattr(sv, "chrom1", None)
    b = getattr(sv, "chrom2", None)
    lo, hi = sorted([a, b])

    # include svtype/pattern to avoid accidental collisions
    return (sv_id, lo, hi, getattr(sv, "svtype", None), getattr(sv, "pattern", None))

# def write_annot(filepath, sveffects, prefix):
#     """
#     :param filepath: outfile for annotation
#     :param sveffects: a list of sveffect classes
#     :return: None
#     """
#     header = '\t'.join(['prefix', 'sv_id', 'chrom1', 'pos1', 'function1', 'gene1', 'transcript_id1', 'strand1', 'transcript_retain1',
#                         'chrom2', 'pos2', 'function2', 'gene2', 'transcript_id2', 'strand2', 'transcript_retain2',
#                         'svpattern', 'svtype', 'fusion'])
#     with open(filepath, 'w') as f:
#         f.write(header + '\n')
#         for sveffect in sveffects:
#             sv_id = getattr(sveffect.sv, "sv_id", None)
#             sv_id = "None" if sv_id is None else str(sv_id)
#             row = [prefix, sv_id] + sveffect.output()
#             f.write('\t'.join('' if x is None else str(x) for x in row) + '\n')
            
# def write_fusion(filepath, svfusions, dict_neo, prefix):
#     """
#     Write fusion-derived neoantigens with provenance columns.

#     dict_neo:
#         {neoepitope: {allele: (affinity, BA_rank, EL_rank, PASS/FILTER)}}
#     """
#     with open(filepath, 'w') as f:
#         header_cols = [
#             'prefix',
#             # Base SV/fusion columns (match SVFusion.output())
#             'sv_id',
#             'chrom1', 'pos1', 'gene1', 'transcript_id1',
#             'chrom2', 'pos2', 'gene2', 'transcript_id2',
#             'svpattern', 'svtype', 'frameshift',
#             # Provenance columns
#             # 'focus_gene', 'bp1_dist_to_focus', 'bp2_dist_to_focus', - removed 9/2/26
#             'junction_nt', 'junction_aa',
#             'num_pep_occurrences', 'pep_start_aa', 'pep_end_aa', 'spans_junction',
#             # NetMHC columns
#             'neoantigen', 'allele', 'affinity', 'BA_rank', 'EL_rank'
#         ]
#         f.write('\t'.join(header_cols) + '\n')

#         for svfusion in svfusions:
#             base = svfusion.output()  # includes sv_id already as col 1

#             # focus_info = getattr(svfusion, 'focus_info', None) or {}
#             # focus_gene = focus_info.get('gene', '')

#             # bp1_dist = getattr(svfusion, 'bp1_dist_to_focus', None)
#             # bp2_dist = getattr(svfusion, 'bp2_dist_to_focus', None)

#             for neoepitope in svfusion.neoepitopes:
#                 if neoepitope not in dict_neo:
#                     continue

#                 meta = getattr(svfusion, 'neoepitope_meta', {}).get(neoepitope, {}) or {}

#                 for allele, vals in dict_neo[neoepitope].items():
#                     if vals[3] != 'PASS':
#                         continue

#                     affinity = vals[0]
#                     ba_rank = vals[1]
#                     el_rank = vals[2]

#                     row = (
#                         [prefix]
#                         + base
#                         + [
#                             # focus_gene,
#                             # bp1_dist,
#                             # bp2_dist,
#                             meta.get('junction_nt'),
#                             meta.get('junction_aa'),
#                             meta.get('num_pep_occurrences'),
#                             meta.get('pep_start_aa'),
#                             meta.get('pep_end_aa'),
#                             meta.get('spans_junction', False),
#                             neoepitope,
#                             allele,
#                             affinity,
#                             ba_rank,
#                             el_rank,
#                         ]
#                     )

#                     # stringify at the edge
#                     f.write('\t'.join('' if x is None else str(x) for x in row) + '\n')




# # def write_fusion(filepath, svfusions, dict_neo, prefix):
# #     """
# #     :param filepath: outfile for fusion derived neoantigen
# #     :param svfusions: a list of svfusion classes
# #     :param dict_neo: a dictionary for neoepitopes from netmhcpan (key: neoepitope, value: {allele: (affinity, BA_rank, EL_rank, pass/fail)})
# #     :param prefix: a string prefix for the sample/patient ID
# #     :return: None
# #     """
# #     with open(filepath, 'w') as f:
# #         header = '\t'.join(['prefix', 'sv_id', 'chrom1', 'pos1', 'gene1', 'transcript_id1',
# #                             'chrom2', 'pos2', 'gene2', 'transcript_id2',
# #                             'svpattern', 'svtype', 'frameshift',
# #                             'neoantigen', 'allele', 'affinity', 'BA_rank', 'EL_rank'])
# #         f.write(header + '\n')
# #         for svfusion in svfusions:
# #             sv_id = getattr(svfusion.sv, "sv_id", None)
# #             sv_id = "None" if sv_id is None else str(sv_id)
            
# #             for neoepitope in svfusion.neoepitopes:
# #                 if neoepitope not in dict_neo:
# #                     continue
                
# #                 for allele, vals in dict_neo[neoepitope].items():
# #                     if vals[3] == 'PASS':
# #                         affinity = str(vals[0]) # binding affinity from neoepitope dictionary
# #                         ba_rank = str(vals[1])
# #                         el_rank = str(vals[2])
# #                         meta = getattr(svfusion, 'neoepitope_meta', {}).get(neoepitope, {})
# #                         focus_info = meta.get('focus_info', None) or {}
# #                         focus_gene_name = focus_info.get('gene', '')
# #                         bp1_dist_to_focus = getattr(svfusion, 'bp1_dist_to_focus', None)
# #                         bp2_dist_to_focus = getattr(svfusion, 'bp2_dist_to_focus', None)
                        
# #                         row = svfusion.output() + [
# #                             str(focus_gene_name),
# #                             '' if bp1_dist_to_focus is None else str(bp1_dist_to_focus),
# #                             '' if bp2_dist_to_focus is None else str(bp2_dist_to_focus),
# #                             '' if meta.get('junction_nt') is None else str(meta.get('junction_nt')),
# #                             '' if meta.get('junction_aa') is None else str(meta.get('junction_aa')),
# #                             '' if meta.get('num_pep_occurrences') is None else str(meta.get('num_pep_occurrences')),
# #                             '' if meta.get('pep_start_aa') is None else str(meta.get('pep_start_aa')),
# #                             '' if meta.get('pep_end_aa') is None else str(meta.get('pep_end_aa')),
# #                             str(meta.get('spans_junction', False)),
# #                             neoepitope, allele, affinity, ba_rank, el_rank
# #                         ]
# #                         f.write('\t'.join([prefix, sv_id] + row) + '\n')
