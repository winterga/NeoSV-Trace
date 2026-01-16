def write_annot(filepath, sveffects, prefix):
    """
    :param filepath: outfile for annotation
    :param sveffects: a list of sveffect classes
    :return: None
    """
    header = '\t'.join(['prefix', 'chrom1', 'pos1', 'function1', 'gene1', 'transcript_id1', 'strand1', 'transcript_retain1',
                        'chrom2', 'pos2', 'function2', 'gene2', 'transcript_id2', 'strand2', 'transcript_retain2',
                        'svpattern', 'svtype', 'fusion'])
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for sveffect in sveffects:
            sv_id = getattr(sveffect.sv, "sv_id", None)
            sv_id = "None" if sv_id is None else str(sv_id)
            f.write('\t'.join([prefix, sv_id] + sveffect.output()) + '\n')


def write_fusion(filepath, svfusions, dict_neo, prefix):
    """
    :param filepath: outfile for fusion derived neoantigen
    :param svfusions: a list of svfusion classes
    :param dict_neo: a dictionary for neoepitopes from netmhcpan (key: neoepitope, value: {allele: (affinity, BA_rank, EL_rank, pass/fail)})
    :param prefix: a string prefix for the sample/patient ID
    :return: None
    """
    with open(filepath, 'w') as f:
        header = '\t'.join(['prefix', 'sv_id', 'chrom1', 'pos1', 'gene1', 'transcript_id1',
                            'chrom2', 'pos2', 'gene2', 'transcript_id2',
                            'svpattern', 'svtype', 'frameshift',
                            'neoantigen', 'allele', 'affinity', 'BA_rank', 'EL_rank'])
        f.write(header + '\n')
        for svfusion in svfusions:
            sv_id = getattr(svfusion.sv, "sv_id", None)
            sv_id = "None" if sv_id is None else str(sv_id)
            
            for neoepitope in svfusion.neoepitopes:
                if neoepitope not in dict_neo:
                    continue
                
                for allele, vals in dict_neo[neoepitope].items():
                    if vals[3] == 'PASS':
                        affinity = str(vals[0])
                        ba_rank = str(vals[1])
                        el_rank = str(vals[2])
                        f.write('\t'.join([prefix, sv_id] + svfusion.output() + [neoepitope, allele, affinity, ba_rank, el_rank]) + '\n')