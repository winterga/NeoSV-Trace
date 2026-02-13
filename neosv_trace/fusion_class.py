from Bio.Seq import Seq


class CDS(object):
    def __init__(self, start, end, intact, startcodon, stopcodon):
        self.start = int(start)
        self.end = int(end)
        self.intact = intact
        self.startcodon = startcodon
        self.stopcodon = stopcodon

    def __str__(self):
        return "%s(start = %d, end = %d, intact = %r, startcodon = %r, stopcodon = %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact,
            self.startcodon,
            self.stopcodon
        )

    def __repr__(self):
        return "%s(%d, %d, %r, %r, %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact,
            self.startcodon,
            self.stopcodon
        )


class CDSCollection(object):
    def __init__(self, transcript, cdslist, part):
        self.transcript = transcript
        self.cdslist = cdslist
        self.part = part

    def __str__(self):
        return "%s(transcript_id = %s, transcript_name = %s, part= %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.part
        )

    def __repr__(self):
        return "%s(%s, %s, %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.part
        )

    @property
    def cds_range(self):
        if self.part == '5':
            start = 0
            end = start + self.cut_length - 1
        else:
            end = self.comp_length - 1
            # pyensembl coding_sequence_position_ranges does not include the stop codon
            # but coding_sequence includes the stop codon, so we should subtract 3
            # additionally for 3' part transcript - https://github.com/openvax/pyensembl/issues/176
            start = end - self.cut_length + 1 - 3
        return start, end

    @property
    def transcript_id(self):
        """
        :return: transcript id of this collection of exons
        """
        return self.transcript.transcript_id

    @property
    def strand(self):
        return self.transcript.strand

    @property
    def gene_name(self):
        """
        :return: gene name of this collection of exons
        """
        return self.transcript.gene.gene_name

    @property
    def transcript_name(self):
        """
        :return: transcript name of this collection of exons
        """
        return self.transcript.transcript_name

    @property
    def nt_sequence(self):
        seq_start, seq_end = self.cds_range
        return self.transcript.coding_sequence[seq_start: seq_end+1]

    @property
    def cut_length(self):
        return sum([(cds.end-cds.start+1) for cds in self.cdslist])

    @property
    def comp_length(self):
        return len(self.transcript.coding_sequence)

# An object to store the fusion of two transcripts, including the SV information, the CDS collections of 
# both sides, and the final nt and aa sequence of the fusion transcript.
class SVFusion(object):
    def __init__(self, sv, cdscollection_1, cdscollection_2):
        self.sv = sv
        self.cc_1 = cdscollection_1
        self.cc_2 = cdscollection_2
        self.nt_sequence = None
        self.aa_sequence = None
        self.neoepitopes = None

    def __str__(self):
        return "%s(sv = %r, cc_1 = %r, cc_2 = %r, nt_sequence = %s, aa_sequence = %s)" % (
            self.__class__.__name__,
            self.sv,
            self.cc_1,
            self.cc_2,
            self.nt_sequence,
            self.aa_sequence
        )

    def __repr__(self):
        return "%s(%r, %r, %r, %s, %s)" % (
            self.__class__.__name__,
            self.sv,
            self.cc_1,
            self.cc_2,
            self.nt_sequence,
            self.aa_sequence
        )

    def __eq__(self, other):
        return self.sv.sorted_coord == other.sv.sorted_coord

    def __hash__(self):
        return hash(self.sv.sorted_coord)

    def is_empty(self):
        if not self.cc_1 or not self.cc_2:
            return True
        
    def get_transcript(self):
        if self.cc_1.strand == '+': return self.cc_1.transcript
        else: return self.cc_2.transcript
        

    def _cdscollection_positions(self, cc, cc_side: int):
        """Return a list of genomic positions for the truncated CDSCollection in transcript coding order.

        Each entry is a dict with keys: chrom, pos, feature, gene, transcript_id.
        """
        if cc is None:
            return []
        chrom = str(cc.transcript.contig).replace('chr', '')
        gene = cc.gene_name
        tid = cc.transcript_id

        # Ensure exon/CDS pieces are ordered 5'->3' for the transcript
        cdslist = list(cc.cdslist or [])
        if cc.strand == '+':
            cdslist.sort(key=lambda c: (c.start, c.end))
        else:
            # on - strand, 5'->3' is descending genomic coordinate
            cdslist.sort(key=lambda c: (c.start, c.end), reverse=True)

        out = []
        for cds in cdslist:
            if cc.strand == '+':
                rng = range(int(cds.start), int(cds.end) + 1)
            else:
                rng = range(int(cds.end), int(cds.start) - 1, -1)
            for p in rng:
                out.append({
                    'type': 'GENOMIC',
                    'chrom': chrom,
                    'pos': int(p),
                    'feature': 'CDS',  # provenance: drawn from transcript CDS pieces
                    'gene': gene,
                    'transcript_id': tid,
                    'strand': cc.strand,
                    'cc_side': int(cc_side),
                })
        return out

    def build_cds_origin(self):
        """Build per-nucleotide origin labels for nt_sequence_cds.

        - Native transcript-derived nucleotides are mapped back to genomic positions.
        - Inserted junction sequence nucleotides are labeled as type='INS'.

        Result is stored on `self.cds_origin` (list, length == len(self.nt_sequence_cds)).
        """
        left_cc, right_cc = (self.cc_1, self.cc_2) if self.cc_1.part == '5' else (self.cc_2, self.cc_1)

        left = self._cdscollection_positions(left_cc, 1 if left_cc is self.cc_1 else 2)
        right = self._cdscollection_positions(right_cc, 1 if right_cc is self.cc_1 else 2)

        ins = [{'type': 'INS'} for _ in range(len(self.nt_sequence_ins))]

        self.cds_origin = left + ins + right
        return self.cds_origin

    @property
    def nt_sequence_ins(self):
        # the insertion sequence is for the forward strand,
        # so should be adjusted by the direction of final transcript
        if self.cc_1.strand == '+':
            return self.sv.insertion
        else:
            return str(Seq(self.sv.insertion).reverse_complement())

    @property
    def nt_sequence_cds(self):
        if self.cc_1.part == '5':
            return self.cc_1.nt_sequence + self.nt_sequence_ins + self.cc_2.nt_sequence
        else:
            return self.cc_2.nt_sequence + self.nt_sequence_ins + self.cc_1.nt_sequence

    @property
    def nt_sequence_3utr(self):
        if self.cc_1.part == '5':
            return self.cc_2.transcript.three_prime_utr_sequence
        else:
            return self.cc_1.transcript.three_prime_utr_sequence

    @property
    def frame_effect(self):
        # fixed a bug: the 5' cds collection could be empty, or less than 3 aa,
        # in this case, traditional start codon is lost and the program will search for the next start codon automatically.
        # but such prediction is of low reliability, we should annotate it and remove these fusions when necessary.
        if self.cc_1.part == '5' and not self.cc_1.nt_sequence.startswith('ATG'):
            return 'Start-loss'
        if self.cc_2.part == '5' and not self.cc_2.nt_sequence.startswith('ATG'):
            return 'Start-loss'
        if len(self.nt_sequence_cds) == 3*(len(self.aa_sequence)+1):
            return 'In-frame'
        elif len(self.nt_sequence_cds) > 3*(len(self.aa_sequence)+1):
            return 'Stop-gain'
        else:
            return 'Stop-loss'

    def output(self):
        tran_id_1 = self.cc_1.transcript_id
        gene_name_1 = self.cc_1.gene_name
        tran_id_2 = self.cc_2.transcript_id
        gene_name_2 = self.cc_2.gene_name
        chrom_1, chrom_2 = self.sv.chrom1, self.sv.chrom2
        pos_1, pos_2 = self.sv.pos1, self.sv.pos2
        
        sv_id = '' if self.sv.sv_id is None else str(self.sv.sv_id)
        
        return [sv_id, chrom_1, str(pos_1), gene_name_1, tran_id_1,
                chrom_2, str(pos_2), gene_name_2, tran_id_2, 
                str(self.sv.pattern), self.sv.svtype, self.frame_effect]