import subprocess
import re


def netmhc_pep_prep(filepath, svfusions):
    """
    :param svfusions: a list of svfusion classes
    :param filepath: the output file
    :return: None
    """
    print("Preparing peptides for netMHC.")
    with open(filepath, 'w') as f:
        for svfusion in svfusions:
            for neoepitope in svfusion.neoepitopes:
                f.write(neoepitope + '\n')


def netmhc_run(netmhcpath, peppath, alleles, outpath):
    """
    :param netmhcpath: absolute path of the netmhc execution file
    :param peppath: input peptides file
    :param alleles: HLA alleles separated by ,
    :param outpath: outfile
    :return: None
    """
    print("Running netMHCpan.")
    cmd = netmhcpath \
        + ' -a ' + alleles \
        + ' -f ' + peppath \
        + ' -inptype 1' \
        + ' -BA' \
        + ' > ' + outpath
    p = subprocess.Popen(cmd, shell=True)
    p.wait()


def netmhc_reload(resultpath):
    """
    :param resultpath: the result of netmhc
    :return: a dictionary {neoepitope: {allele: [affinity, rank, FILTER]}}
    """
    pep_dic = {}
    with open(resultpath, 'r') as f:
        for line in f:
            if re.match(r'^\s+1\s+H', line):
                tmpline = re.split(r'\s+', line)
                allele = tmpline[2]
                pep = tmpline[3]
                affinity = float(tmpline[16])
                ba_rank = float(tmpline[15])
                el_rank = float(tmpline[13])
                if pep not in pep_dic:
                    pep_dic[pep] = {}
                pep_dic[pep][allele] = [affinity, ba_rank, el_rank, 'FILTER']
    return pep_dic


# ------------------
# mhcflurry support
# ------------------

def mhcflurry_run(mhcflurry_path, peppath, alleles, outpath, chunk_size=500):
    '''
    Run mhcflurry-predict on the given peptides and alleles, writing CSV output to outpath.

    Parameters
    ----------
    mhcflurry_path : str
        Path to mhcflurry-predict executable (or command on PATH).
    peppath : str
        File with one peptide per line.
    alleles : str
        Comma-separated HLA alleles, e.g. 'HLA-A02:01,HLA-B07:02'.
    outpath : str
        Output CSV file path.
    chunk_size : int
        Number of peptides to send per CLI call to avoid OS argv limits.

    '''
    print("Running MHCflurry.")
    # MHCflurry expects alleles separated by spaces, not commas like netMHCpan
    allele_list = [a.strip() for a in alleles.split(',')]

    with open(peppath, 'r') as pep_file:
        peptides = [line.strip() for line in pep_file if line.strip()]

    if not peptides:
        open(outpath, 'w').close()
        return
    
    header_written = False
    with open(outpath, 'w') as outf:
        for i in range(0, len(peptides), chunk_size):
            chunk = peptides[i:i+chunk_size]

            # mhcflurry-predict prints output to stdout
            cmd = [mhcflurry_path, '--alleles'] + allele_list + ['--peptides'] + chunk
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            out, err = p.communicate()

            if p.returncode != 0:
                raise RuntimeError(f"mhcflurry-predict failed: (exit {p.returncode}).\nSTDERR:\n{err}")
            lines = [line for line in out.splitlines() if line.strip()]

            if not lines:
                continue
            if not header_written:
                outf.write(lines[0] + "\n")
                start_idx = 1
                header_written = True
            else:
                if lines[0].lower().startswith('allele,peptide,'):
                    start_idx = 1
                else:
                    start_idx = 0
                
            for line in lines[start_idx:]:
                outf.write(line + "\n")

def mhcflurry_reload(resultpath):
    '''
    Load MHCflurry CSV output into a dictionary.

    Parameters
    ----------
    resultpath : str
        Path to MHCflurry CSV output file.

    Returns
    -------
    dict
        A dictionary {neoepitope: {allele: [affinity, rank, FILTER]}}.

    '''
    pep_dic = {}
    with open(resultpath, 'r') as f:
        first = f.read(1)
        if not first:
            return pep_dic # empty file
        f.seek(0)

        header = f.readline().strip().split(',')
        idx = {name: i for i, name in enumerate(header)}
        required_cols = ['allele', 'peptide', 'mhcflurry_affinity', 'mhcflurry_affinity_percentile', 'mhcflurry_presentation_percentile']
        for col in required_cols:
            if col not in idx:
                raise ValueError(f"Unexpected mhcflurry output: missing column '{col}'. Columns={header}")
        presentation_col = 'mhcflurry_presentation_percentile'
        has_presentation = presentation_col in idx

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            allele = parts[idx['allele']]
            peptide = parts[idx['peptide']]
            affinity = float(parts[idx['mhcflurry_affinity']])
            ba_rank = float(parts[idx['mhcflurry_affinity_percentile']])
            if has_presentation:
                el_rank = float(parts[idx[presentation_col]])
            else:
                el_rank = 100.0 # if not provided, set to high value so it won't pass filtering
            
            if peptide not in pep_dic:
                pep_dic[peptide] = {}
            pep_dic[peptide][allele] = [affinity, ba_rank, el_rank, 'FILTER']
        return pep_dic
    
def mhc_predict_pep_prep(filepath, svfusions):
    # predicctor-agnostic peptide prep
    netmhc_pep_prep(filepath, svfusions)

def mhc_predict_run(predictor, predictor_path, peppath, alleles, outpath):
    if predictor == 'netmhcpan':
        netmhc_run(predictor_path, peppath, alleles, outpath)
    elif predictor == 'mhcflurry':
        mhcflurry_run(predictor_path, peppath, alleles, outpath)
    else:
        raise ValueError(f"Unsupported MHC predictor: {predictor}")

def mhc_predict_reload(predictor, resultpath):
    if predictor == 'netmhcpan':
        return netmhc_reload(resultpath)
    elif predictor == 'mhcflurry':
        return mhcflurry_reload(resultpath)
    else:
        raise ValueError(f"Unsupported MHC predictor: {predictor}")

def mhc_filter(pep_dic, aff_thre, ba_rank_thre, el_rank_thre):
    """
    :param pep_dic: the dictionary for neoepitopes generated by netmhc_load
    :param aff_thre: threshold for binding affinity
    :param ba_rank_thre: threshold for BA rank percentile
    :param el_rank_thre: threshold for EL rank percentile
    :return: a dictionary in which passed neoepitopes are labeled PASS
    """
    print("Filtering MHC predictor results")
    for pep in pep_dic:
        for allele in pep_dic[pep]:
            if pep_dic[pep][allele][0] <= aff_thre and pep_dic[pep][allele][1] <= ba_rank_thre and pep_dic[pep][allele][2] <= el_rank_thre:
                pep_dic[pep][allele][3] = 'PASS'
    return pep_dic
