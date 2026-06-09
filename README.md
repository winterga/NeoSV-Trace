# NeoSV-Trace

**Predict neoantigens arising from structural variants in cancer genomes.**
*NOTE:* This repository is highly based off of [NeoSV](https://github.com/ysbioinfo/NeoSV). NeoSV-Trace is an extension of NeoSV, so much of the logic described here is directly from NeoSV.

---

## Table of Contents

1. [Background](#background)
2. [How It Works](#how-it-works)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Input Formats](#input-formats)
6. [Full CLI Reference](#full-cli-reference)
7. [Output Files](#output-files)
8. [MHC Predictors](#mhc-predictors)
9. [Ensembl Reference](#ensembl-reference)
10. [Citation / License](#citation--license)

---

## Background

### Structural Variants in Cancer

Structural variants (SVs) are large-scale rearrangements of the genome — deletions, duplications, inversions, and translocations — that move, flip, or fuse genomic segments. Unlike single-nucleotide variants (SNVs), SVs can join sequences from entirely different chromosomes or genomic loci. In cancer, SVs are common and often drive tumour biology: the BCR-ABL fusion in chronic myeloid leukaemia is an archetypal example of an SV producing a pathological gene fusion.

### Neoantigens and Cancer Immunotherapy

When a tumour cell acquires a somatic mutation, the resulting altered protein is broken down into short peptides (typically 8–11 amino acids) by the proteasome. These peptides are loaded onto MHC (major histocompatibility complex) molecules — called HLA molecules in humans — and displayed on the cell surface. If a peptide derives from a tumour-specific mutation, it may be recognised by cytotoxic T cells as foreign, triggering an immune response. Such tumour-specific peptides are called **neoantigens**. Neoantigens are a central focus of cancer immunotherapy: they underlie the efficacy of checkpoint inhibitors and are the basis of personalised cancer vaccines.

### MHC Binding and HLA Typing

Whether a peptide becomes a neoantigen depends critically on whether it binds the patient's HLA molecules. HLA genes are among the most polymorphic in the human genome — each person carries a unique combination of alleles (e.g. `HLA-A*02:01`), and each allele has distinct preferences for which peptide sequences it will bind. A peptide that binds tightly (low IC50) and ranks well among random peptides (low percentile rank) is a strong candidate neoantigen. This means neoantigen prediction is patient-specific: you need both the tumour's mutations and the patient's HLA type.

### Why SVs Need Their Own Tool

Most neoantigen prediction pipelines are built around SNVs and small indels, which alter a single protein at a predictable location. SVs are fundamentally different: a translocation can join the coding sequences of two genes on different chromosomes, creating a fusion protein with a completely novel amino acid sequence at the junction. This junction sequence does not exist in the normal proteome and may be highly immunogenic. However, correctly predicting the peptide sequence at an SV junction requires knowing the precise genomic coordinates of both breakpoints, the transcript structures involved, the reading frame, and whether the junction is in-frame or causes a frameshift. NeoSV-Trace automates this entire process — from raw SV calls to filtered, patient-specific neoantigen candidates.

---

## How It Works

NeoSV-Trace runs a sequential pipeline from SV calls to predicted neoantigens:

```
┌─────────────────────────┐
│   SV calls              │
│   (.vcf or .bedpe)      │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│   SV annotation         │  ← PyEnsembl (Ensembl reference)
│   Breakpoint genomic    │
│   context, gene,        │
│   transcript, function  │
│   → <prefix>.anno.txt   │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│   Fusion transcript     │  ← SV pattern inference (DEL/DUP/INV/TRA)
│   construction          │     reading frame, strand, retained exons
│   Nucleotide + AA seqs  │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│   Candidate neopeptides │  ← Sliding window over fusion AA sequence
│   8–11-mers spanning    │     (length range configurable)
│   or near SV junction   │
│   → <prefix>.all_       │
│     neopeptides.txt     │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│   MHC binding prediction│  ← netMHCpan or MHCflurry
│   IC50, BA rank,        │     patient HLA alleles
│   EL rank per allele    │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│   Filtered neoantigens  │  ← IC50 ≤ 500 nM AND BA rank ≤ 2% AND
│   → <prefix>.           │     EL rank ≤ 2% (all thresholds configurable)
│     neoantigen.txt      │
└─────────────────────────┘
```

**Step-by-step:**

1. **Load reference** — PyEnsembl loads the chosen Ensembl release, providing transcript structures, exon boundaries, and coding sequences.
2. **Parse SV calls** — VCF breakend (BND) records or BEDPE rows are read and converted into `StructuralVariant` objects, with duplicate SVs removed.
3. **Infer SV patterns** — Each SV is classified (deletion, duplication, inversion, translocation) and the orientation of the two breakpoints is resolved.
4. **Annotate breakpoints** — Each breakpoint is mapped to a gene and transcript; its genomic function (exon, intron, UTR) and retained transcript segments are recorded. Results are written to `<prefix>.anno.txt`.
5. **Build fusion transcripts** — For SVs that create gene fusions, the nucleotide sequence of the fusion transcript is assembled and translated to an amino acid sequence.
6. **Generate candidate peptides** — A sliding window generates all k-mers of the specified lengths. All candidates are written to `<prefix>.all_neopeptides.txt`.
7. **Run MHC predictor** — Unique peptides are submitted to netMHCpan or MHCflurry against the patient's HLA alleles.
8. **Filter and annotate** — Peptides passing the binding thresholds are annotated with their genomic origin and written to `<prefix>.neoantigen.txt`.

---

## Installation

### Python package

```bash
pip install git+https://github.com/macintyrelab/NeoSV-Trace.git
```

**Python >= 3.8 required.** Python dependencies (`pandas`, `numpy`, `pyensembl`, `biopython`) are installed automatically.

### External MHC predictor (required)

NeoSV-Trace does not bundle an MHC binding predictor. You must install at least one:

**Option A — netMHCpan** (recommended for most users)

netMHCpan is free for academic use and must be downloaded from the [DTU Health Tech server](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/). Follow their installation instructions and ensure the `netMHCpan` executable is on your `PATH`, or supply the path with `--netmhc-path`.

**Option B — MHCflurry**

MHCflurry is open-source and pip-installable:

```bash
pip install mhcflurry
mhcflurry-downloads fetch
```

Ensure `mhcflurry-predict` is on your `PATH`, or supply the directory with `--mhcflurry-path`.

### Ensembl reference cache

NeoSV-Trace uses PyEnsembl to look up transcript structures. Download the reference for your genome build before running:

```bash
# hg19 / GRCh37 (Ensembl release 75 — the default)
pyensembl install --release 75 --species human

# hg38 / GRCh38 (Ensembl release 95 or 115)
pyensembl install --release 95 --species human
```

See [Ensembl Reference](#ensembl-reference) for a full list of supported releases.

---

## Quick Start

```bash
neosv-trace \
  --sv-file    tumor_svs.bedpe \
  --hla-file   patient_hla.txt \
  --out        results/ \
  --prefix     sample01 \
  --release    75
```

This uses the default settings: hg19 reference (Ensembl release 75), netMHCpan predictor, 8–11-mer peptides, IC50 ≤ 500 nM, BA rank ≤ 2%, EL rank ≤ 2%.

To run annotation only (no HLA file needed):

```bash
neosv-trace \
  --sv-file    tumor_svs.vcf \
  --out        results/ \
  --prefix     sample01 \
  --release    75 \
  --anno-only
```

---

## Input Formats

### SV calls — VCF format

A VCF (Variant Call Format) file is the standard output of most SV callers (DELLY, MANTA, SVABA, etc.). NeoSV-Trace reads **breakend (BND) records**, which encode paired genomic positions where DNA sequence continuity is disrupted.

A BND ALT field encodes which strand the two segments join on. For example:

- `A]4:56789]` — the reference base `A` at this position is followed by the sequence at chromosome 4, position 56789, read in the forward direction.
- `]4:56789]A` — the sequence arriving from chromosome 4:56789 (reverse complemented) precedes the reference base `A`.

Non-BND records (e.g. `<DEL>` symbolic alleles) are skipped with a warning.

**Required columns (tab-separated):**

| Column index | Field | Notes |
|---|---|---|
| 0 | CHROM | `chr` prefix stripped automatically |
| 1 | POS | 1-based breakpoint position |
| 2 | ID | Used as `sv_id` |
| 3 | REF | Reference base |
| 4 | ALT | BND notation (see above) |

Lines beginning with `#` are treated as header/comment lines and skipped.

**Example:**

```
##fileformat=VCFv4.1
#CHROM  POS     ID   REF  ALT           QUAL  FILTER  INFO
1       12345   sv1  A    A]4:56789]    .     .       .
4       56789   sv1  T    T]1:12345]    .     .       .
```

---

### SV calls — BEDPE format

BEDPE (Bed Paired End) is a tab-separated format used by tools such as LUMPY and common in population-scale SV datasets. Each row describes one SV using two genomic intervals — one for each breakpoint.

**Required columns (looked up by name, order-independent):**

| Column name | Description |
|---|---|
| `chrom1` | Chromosome of the first breakpoint |
| `start1` | Start position of the first breakpoint interval |
| `chrom2` | Chromosome of the second breakpoint |
| `start2` | Start position of the second breakpoint interval |
| `strand1` | Orientation at breakpoint 1 (`+` or `-`) |
| `strand2` | Orientation at breakpoint 2 (`+` or `-`) |
| `sv_id` | (Optional) identifier for the SV; used in output for traceability |

Additional columns (`end1`, `end2`, `score`, `svtype`, etc.) are accepted and ignored. The file must have a header row.

**Example:**

```
chrom1  start1     end1       chrom2  start2     end2       sv_id   score  strand1  strand2  svtype
1       79457189   79457190   11      93161656   93161657   SV_1    .      -        +        TRA
2       205434182  205434183  2       205554216  205554217  SV_2    .      +        -        DEL
```

---

### HLA alleles file

HLA typing identifies which HLA alleles a patient carries. Each HLA molecule is encoded by a highly polymorphic gene, and the 4-digit (two-field) resolution — e.g. `HLA-A*02:01` — is the minimum needed for accurate binding prediction. The first two digits identify the allele group; the next two specify the exact protein sequence.

The HLA file should contain one allele per line. Most people carry two alleles each for HLA-A, HLA-B, and HLA-C (six alleles total), but the file can contain any number. The asterisk (`*`) is stripped automatically, so both formats below are accepted.

**Example (`patient_hla.txt`):**

```
HLA-A*02:01
HLA-A*24:02
HLA-B*07:02
HLA-B*35:01
HLA-C*07:02
HLA-C*04:01
```

---

## Full CLI Reference

```
neosv-trace [options]
```

| Flag | Dest | Type | Default | Required | Description |
|---|---|---|---|---|---|
| `-sf` / `--sv-file` | `svfile` | string | — | **Yes** | SV input file. Format auto-detected from extension: `.vcf` or `.bedpe`. |
| `-hf` / `--hla-file` | `hlafile` | string | — | **Yes** (unless `--anno-only`) | HLA alleles file. One allele per line, 4-digit resolution. |
| `-mp` / `--mhc-predictor` | `mhc_predictor` | choice | `netmhcpan` | No | MHC binding predictor to use. Choices: `netmhcpan`, `mhcflurry`. |
| `-np` / `--netmhc-path` | `netmhc` | string | auto (PATH) | No | Absolute path to the `netMHCpan` executable. Not needed if on PATH. |
| `-fp` / `--mhcflurry-path` | `mhcflurry` | string | auto (PATH) | No | Absolute path to the MHCflurry installation directory. Not needed if `mhcflurry-predict` is on PATH. |
| `-o` / `--out` | `outdir` | string | — | **Yes** | Output directory. Created if it does not exist. |
| `-p` / `--prefix` | `prefix` | string | `sample` | No | Prefix for all output filenames. |
| `-r` / `--release` | `release` | string | `75` | No | Ensembl release number (54–115) or `custom`. hg18=54, hg19=75, hg38=95 or 115. |
| `-gf` / `--gtf-file` | `gtffile` | string | — | Only with `--release custom` | GTF annotation file for a custom reference genome. |
| `-cf` / `--cdna-file` | `cdnafile` | string | — | Only with `--release custom` | cDNA FASTA file for a custom reference genome. |
| `-pd` / `--pyensembl-cache-dir` | `cachedir` | string | platform default | No | Directory where PyEnsembl stores its index files. |
| `-l` / `--epitope-lengths` | `window` | string | `8-11` | No | Range of peptide lengths to generate, as `min-max` (e.g. `8-11` generates 8, 9, 10, 11-mers). |
| `-ic` / `--ic50-cutoff` | `aff_cutoff` | float | `500` | No | IC50 binding affinity cutoff in nM. Peptides with affinity above this value are excluded. |
| `-brc` / `--ba-ranking-cutoff` | `ba_rank_cutoff` | float | `2.0` | No | Binding affinity rank percentile cutoff. Peptides ranked above this percentile are excluded. |
| `-erc` / `--el-ranking-cutoff` | `el_rank_cutoff` | float | `2.0` | No | Eluted ligand rank percentile cutoff. Peptides ranked above this percentile are excluded. |
| `-ct` / `--complete-transcript` | `complete` | flag | `False` | No | If set, only complete (fully annotated) transcripts are used for SV annotation. |
| `--anno-only` | `anno` | flag | `False` | No | Run annotation only. Skips HLA loading and MHC prediction. `--hla-file` is not required. |

**Validation notes:**

- `--release` must be an integer between 54 and 115 inclusive, or the string `custom`.
- `--release custom` requires both `--gtf-file` and `--cdna-file`.
- If `--mhc-predictor netmhcpan`, the tool `netMHCpan` or `netMHCpan-4.1` must be on PATH or `--netmhc-path` must be provided.
- If `--mhc-predictor mhcflurry`, `mhcflurry-predict` must be on PATH or `--mhcflurry-path` must be provided.

---

## Output Files

All output files are tab-separated. The `prefix` column in each file reflects the stem of the input SV filename (not the `-p` / `--prefix` flag, which controls output filenames only).

---

### `<prefix>.anno.txt` — SV Annotation

Written after breakpoint annotation. Contains one row per SV per transcript pair considered. This file is always produced, even with `--anno-only`.

| Column | Description |
|---|---|
| `prefix` | Stem of the input SV filename |
| `sv_id` | SV identifier from the input file |
| `chrom1` | Chromosome of breakpoint 1 |
| `pos1` | Position of breakpoint 1 |
| `function1` | Genomic function at breakpoint 1 (e.g. `Exon`, `Intron`, `UTR`) |
| `gene1` | Gene name at breakpoint 1 |
| `transcript_id1` | Ensembl transcript ID at breakpoint 1 |
| `strand1` | Strand of the transcript at breakpoint 1 |
| `transcript_retain1` | Portion of transcript 1 retained in the fusion (e.g. `E1-i1` = exon 1 through intron 1) |
| `chrom2` | Chromosome of breakpoint 2 |
| `pos2` | Position of breakpoint 2 |
| `function2` | Genomic function at breakpoint 2 |
| `gene2` | Gene name at breakpoint 2 |
| `transcript_id2` | Ensembl transcript ID at breakpoint 2 |
| `strand2` | Strand of the transcript at breakpoint 2 |
| `transcript_retain2` | Portion of transcript 2 retained in the fusion |
| `svpattern` | Integer code for the resolved SV orientation pattern |
| `svtype` | SV class (`TRA`, `DEL`, `DUP`, `INV`) |
| `fusion` | Whether the SV produces a gene fusion (`Fusion` or `No fusion`) |

**Example row:**

```
POG-CA  None  1  12174267  Intron  TNFRSF1B  ENST00000376259  +  E1-i1  4  140362305  Intron  SCOC  ENST00000394205  +  i2-E5  1  TRA  Fusion
```

---

### `<prefix>.all_neopeptides.txt` — All Candidate Neopeptides

Written before MHC filtering. Contains all k-mer peptides generated from fusion sequences, regardless of predicted binding. Useful for inspecting the full candidate set or running alternative binding analyses.

| Column | Description |
|---|---|
| `prefix` | Stem of the input SV filename |
| `sv_id` | SV identifier |
| `chrom1` | Chromosome of breakpoint 1 |
| `pos1` | Position of breakpoint 1 |
| `gene1` | Gene at breakpoint 1 |
| `transcript_id1` | Transcript at breakpoint 1 |
| `chrom2` | Chromosome of breakpoint 2 |
| `pos2` | Position of breakpoint 2 |
| `gene2` | Gene at breakpoint 2 |
| `transcript_id2` | Transcript at breakpoint 2 |
| `svpattern` | Integer SV orientation pattern code |
| `svtype` | SV class |
| `frameshift` | Reading frame consequence (`In-frame` or `Frameshift`) |
| `junction_nt` | Nucleotide index of the SV junction in the assembled fusion sequence |
| `junction_aa` | Amino acid index of the SV junction in the fusion protein |
| `spans_junction` | `True` if this peptide spans the SV junction; `False` if it falls entirely within one gene's sequence |
| `neopeptide` | Peptide amino acid sequence |
| `pep_length` | Length of the peptide in amino acids |

**Example row:**

```
23561  BEDPE_2_205434182_2_205554216_1  2  205434182  PARD3B  ENST00000406610  2  205554216  PARD3B  ENST00000406610  1  DEL  In-frame  120  40  True  LQRYLKTREKK  11
```

---

### `<prefix>.neoantigen.txt` — Filtered Neoantigens

The primary output. Contains only peptides that passed all MHC binding filters, annotated with their genomic origin and per-allele binding scores. Rows are deduplicated by `(sv_id, peptide, allele)`.

| Column | Description |
|---|---|
| `prefix` | Stem of the input SV filename |
| `sv_id` | SV identifier |
| `chrom1` | Chromosome of breakpoint 1 |
| `pos1` | Position of breakpoint 1 |
| `gene1` | Gene at breakpoint 1 |
| `transcript_id1` | Transcript at breakpoint 1 |
| `chrom2` | Chromosome of breakpoint 2 |
| `pos2` | Position of breakpoint 2 |
| `gene2` | Gene at breakpoint 2 |
| `transcript_id2` | Transcript at breakpoint 2 |
| `svpattern` | Integer SV orientation pattern code |
| `svtype` | SV class |
| `frameshift` | Reading frame consequence (`In-frame` or `Frameshift`) |
| `junction_nt` | Nucleotide index of the SV junction |
| `junction_aa` | Amino acid index of the SV junction |
| `num_pep_occurrences` | Number of times this peptide sequence occurs in the fusion protein |
| `pep_start_aa` | 0-based start position of the peptide in the fusion amino acid sequence |
| `pep_end_aa` | 0-based end position of the peptide in the fusion amino acid sequence |
| `spans_junction` | Whether the peptide spans the SV junction |
| `pep_origin_type` | Classification of which component of the fusion the peptide originates from |
| `pep_origin` | Detailed genomic origin annotation |
| `pep_feature` | Genomic feature annotation for the peptide origin |
| `neoantigen` | Peptide amino acid sequence |
| `allele` | HLA allele (e.g. `HLA-A02:01`) |
| `affinity` | Predicted IC50 binding affinity in nM (lower = stronger binding) |
| `BA_rank` | Binding affinity rank percentile (lower = stronger relative to random peptides) |
| `EL_rank` | Eluted ligand rank percentile (lower = more likely to be naturally presented) |

**Example row:**

```
POG-CA  sv_11  11  128782118  FLI1  ENST00000344954  22  29288739  EWSR1  ENST00000629659  3  TRA  In-frame  796  265  1  263  271  True  ...  QQNPSYDSV  HLA-A02:01  917.14  1.686  2.335
```

---

### Intermediate files

These files are written to the output directory and can be useful for debugging or reuse:

| File | Description |
|---|---|
| `<prefix>.net.in.txt` | One unique peptide per line; input submitted to the MHC predictor |
| `<prefix>.net.out.txt` | Raw output from the MHC predictor (netMHCpan text or MHCflurry CSV) |
| `<prefix>.net.out.txt.input.csv` | (MHCflurry only) Cartesian product of all alleles × all peptides sent to `mhcflurry-predict` |
| `<prefix>.net.out.txt.log` | (MHCflurry only) Subprocess stdout/stderr from the MHCflurry call |

---

## MHC Predictors

NeoSV-Trace supports two external MHC binding predictors. Both predict IC50 affinity and rank percentile scores; the choice affects accuracy, speed, and licensing.

### netMHCpan

- **Algorithm:** Pan-specific neural network trained on MHC ligand data from IEDB and other sources. The current version (4.1) is widely used and benchmarked.
- **Scores:** Binding affinity (BA) IC50 in nM, BA rank percentile, and eluted ligand (EL) rank percentile. EL rank is generally a better predictor of natural presentation than BA alone.
- **License:** Free for academic, non-commercial use. Requires manual download and registration from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/).
- **Invocation:** `netMHCpan -a <alleles> -f <pepfile> -inptype 1 -BA`

### MHCflurry

- **Algorithm:** Neural network predictor; open-source and pip-installable.
- **Scores:** Affinity IC50, affinity percentile (BA rank), and presentation percentile (EL rank). If the presentation percentile is absent from the output, EL rank is set to 100.0 and the peptide will not pass the EL rank filter.
- **License:** Apache 2.0 (open source). No registration required.
- **Invocation:** `mhcflurry-predict <input.csv> --out <output.csv>`

### Which predictor to choose?

| Consideration | netMHCpan | MHCflurry |
|---|---|---|
| License | Academic/non-commercial only | Open source (Apache 2.0) |
| Installation | Manual download required | `pip install mhcflurry` |
| Benchmark performance | Extensively published | Strong, slightly less benchmarked |
| EL rank support | Yes | Yes (if presentation model downloaded) |

For most research use cases, **netMHCpan is recommended** because it has the most extensive independent benchmarking. MHCflurry is a good choice when open-source licensing is required or in automated pipelines.

### Filter thresholds

All three filters must be satisfied simultaneously (AND logic):

| Filter | Flag | Default | Meaning |
|---|---|---|---|
| IC50 affinity | `--ic50-cutoff` | 500 nM | Strong binders conventionally defined as < 500 nM |
| BA rank | `--ba-ranking-cutoff` | 2.0% | Top 2% of random peptides by binding affinity |
| EL rank | `--el-ranking-cutoff` | 2.0% | Top 2% of random peptides by eluted ligand score |

The BA and EL rank cutoffs are generally more robust than the raw IC50 because they account for allele-specific score distributions. Loosening any cutoff (e.g. `--el-ranking-cutoff 10`) will increase sensitivity at the cost of more false positives.

---

## Ensembl Reference

NeoSV-Trace uses [PyEnsembl](https://github.com/openvax/pyensembl) to retrieve transcript structures, exon boundaries, and coding sequences. You must download the appropriate Ensembl release for your genome build before running.

### Supported releases and genome builds

| Genome build | Ensembl release | Download command |
|---|---|---|
| hg18 / GRCh36 | 54 | `pyensembl install --release 54 --species human` |
| hg19 / GRCh37 | **75** (default) | `pyensembl install --release 75 --species human` |
| hg38 / GRCh38 | 95 | `pyensembl install --release 95 --species human` |
| hg38 / GRCh38 | 115 | `pyensembl install --release 115 --species human` |
| Any | 54–115 | `pyensembl install --release <N> --species human` |

Any integer release between 54 and 115 is accepted. **Release 75 is the default** because it corresponds to GRCh37/hg19, which remains the most common coordinate system in published cancer SV datasets.

### Custom reference genome

If your SV calls are in coordinates not covered by a standard Ensembl release (e.g. a non-human organism or a patched assembly), you can supply a custom GTF annotation and cDNA FASTA:

```bash
neosv-trace \
  --sv-file    tumor_svs.bedpe \
  --hla-file   patient_hla.txt \
  --out        results/ \
  --release    custom \
  --gtf-file   /path/to/annotation.gtf \
  --cdna-file  /path/to/cdna.fa
```

### PyEnsembl cache location

By default, PyEnsembl stores its index files in a platform-specific cache directory. To specify a custom location (e.g. on a shared filesystem in an HPC environment):

```bash
neosv-trace \
  --sv-file    tumor_svs.bedpe \
  --hla-file   patient_hla.txt \
  --out        results/ \
  --pyensembl-cache-dir /scratch/shared/pyensembl_cache
```

The same directory must be used when pre-downloading the reference:

```bash
PYENSEMBL_CACHE_DIR=/scratch/shared/pyensembl_cache \
  pyensembl install --release 75 --species human
```

---

## Citation / License

If you use NeoSV-Trace in your research, please cite:

> Wintergerst G, *et al.* NeoSV-Trace: prediction of neoantigens arising from structural variants in cancer genomes. *(in preparation)*

and

> Shi, Y., Jing, B. & Xi, R. Comprehensive analysis of neoantigens derived from structural variation across whole genomes from 2528 tumors. Genome Biol 24, 169 (2023)

**License:** See `LICENSE` file in this repository.

**Author:** Greyson Wintergerst — [gwintergerst@ext.cnio.es](mailto:gwintergerst@ext.cnio.es)