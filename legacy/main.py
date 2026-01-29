import os
import sys
from .args import create_arg_parser
from .input import vcf_load, bedpe_load, hla_load, ensembl_load, get_window_range
from .sv_utils import sv_pattern_infer_vcf, sv_pattern_infer_bedpe, remove_duplicate, dist_to_interval
from .annotation_utils import sv_to_sveffect
from .fusion_utils import sv_to_svfusion
<<<<<<< HEAD
from .sequence_utils import set_nt_seq, set_aa_seq, generate_neoepitopes
from .mhc import mhc_predict_pep_prep, mhc_predict_run, mhc_predict_reload, mhc_filter
=======
from .sequence_utils import generate_neoepitopes_with_meta, set_nt_seq, set_aa_seq, generate_neoepitopes
from .netMHC import netmhc_pep_prep, netmhc_run, netmhc_reload, netmhc_filter
>>>>>>> trace-test1
from .output import write_annot, write_fusion


def main():
    # generate the args
    args = create_arg_parser()
    # load a pyensembl Genome class
    ensembl = ensembl_load(args.release, args.gtffile, args.cdnafile, args.cachedir)
    
    # get coordinates for the focus gene (default: FHIT) to support breakpoint-distance output columns

    focus_gene_id = args.focus_gene_id
    focus_info = None
    
    try: 
        gene = ensembl.gene_by_id(focus_gene_id)
        gene_name = ensembl.gene_name_of_gene_id(focus_gene_id)
        focus_info = {
            'gene': gene_name,
            'chrom': str(gene.contig).replace('chr', ''),
            'start': int(gene.start),
            'end': int(gene.end)
        }
        
    except Exception:
        focus_info = None
        print(f"Warning: Could not find gene ID {focus_gene_id} in the reference. Breakpoint-distance columns will be empty.")

    if args.svfile.endswith('.vcf'):
        # load vcf file and transform SVs into VariantCallingFormat class
        sv_vcfs = vcf_load(args.svfile)
        # transform SV from VariantCallingFormat class to StructuralVariant class
        svs = [sv_pattern_infer_vcf(sv_vcf) for sv_vcf in sv_vcfs]
    elif args.svfile.endswith('.bedpe'):
        # load bedpe file SVs
        sv_beds = bedpe_load(args.svfile)
        # transform SV from BEDPE class to StructuralVariant class
        svs = [sv_pattern_infer_bedpe(sv_bed) for sv_bed in sv_beds]
    else:
        sys.exit('The input file must end with .vcf or .bedpe. Other format is not allowed.')
    # remove duplicated SVs
    svs = remove_duplicate(svs)
    
    # Get patient/sample ID from svfile
    neoprefix = os.path.basename(args.svfile).split('.')[0]

    # annotate SVs and write to file_anno
    sv_effects = [sv_to_sveffect(sv, ensembl, args.complete) for sv in svs]
    sv_effects_flat = [sv_effect_unit for sv_effect in sv_effects for sv_effect_unit in sv_effect]
    file_anno = os.path.join(args.outdir, args.prefix + '.anno.txt')
    write_annot(file_anno, sv_effects_flat, prefix=neoprefix) # FIXME needs to incorporate adding sv_id so that it can be traced back from neoantigen file

    if not args.anno:
        # load the hla file and join the alleles by ,
        hla_alleles = hla_load(args.hlafile)
        # specify the lengths of neoepitopes you want to predict
        window_range = get_window_range(args.window)
        # transform StructuralVariant class to SVFusion class
        sv_fusions = [sv_to_svfusion(sv, ensembl) for sv in svs]
        # remove those empty SVFusions
        sv_fusions = [sv_fusion for sv_fusion in sv_fusions if not sv_fusion.is_empty()]

        # predict neoepitopes for each SVFusion
        for sv_fusion in sv_fusions:
            sv_fusion.nt_sequence = set_nt_seq(sv_fusion)
            sv_fusion.aa_sequence = set_aa_seq(sv_fusion)
<<<<<<< HEAD
            sv_fusion.neoepitopes = generate_neoepitopes(sv_fusion, window_range)

        # predict binding affinity using MHC predictor
        file_mhc_in = os.path.join(args.outdir, args.prefix + '.net.in.txt')
        file_mhc_out = os.path.join(args.outdir, args.prefix + '.net.out.txt')
=======
            sv_fusion.neoepitopes, sv_fusion.neoepitope_meta = generate_neoepitopes_with_meta(sv_fusion, window_range)
            sv_fusion.focus_info = focus_info
            
            if focus_info is not None:
                # distance is 0 if breakpoint lies inside the gene interval; otherwise min distance to edge
                
                sv_fusion.bp1_dist_to_focus = dist_to_interval(focus_info=focus_info, chrom=sv_fusion.sv.chrom1, pos=sv_fusion.sv.pos1)
                sv_fusion.bp2_dist_to_focus = dist_to_interval(focus_info=focus_info, chrom=sv_fusion.sv.chrom2, pos=sv_fusion.sv.pos2)
            else:
                sv_fusion.bp1_dist_to_focus = None
                sv_fusion.bp2_dist_to_focus = None
        
        # predict binding affinity using netMHC
        file_netmhc_in = os.path.join(args.outdir, args.prefix + '.net.in.txt')
        file_netmhc_out = os.path.join(args.outdir, args.prefix + '.net.out.txt')
>>>>>>> trace-test1
        if sv_fusions:
            mhc_predict_pep_prep(file_mhc_in, sv_fusions)
            mhc_predict_run(args.netmhc, file_mhc_in, hla_alleles, file_mhc_out)
        else:
            open(file_mhc_in, 'w').close()
            open(file_mhc_out, 'w').close()

        # reload and filter the neoepitopes based on netMHC result
        dict_epitope = mhc_predict_reload(args.netmhc, file_mhc_out)
        dict_epitope = mhc_filter(dict_epitope, args.aff_cutoff, args.ba_rank_cutoff, args.el_rank_cutoff)

        # generate the final outfile
        file_fusion = os.path.join(args.outdir, args.prefix + '.neoantigen.txt')
        write_fusion(file_fusion, sv_fusions, dict_epitope, prefix=neoprefix)


if __name__ == '__main__':
    main()
