class StructuralVariant(object):
    """Container for a single structural variant / BND-style breakpoint pair.

    Core fields are intentionally minimal so both VCF(BND) and BEDPE inputs can map here.

    Attributes
    ----------
    chrom1, pos1, chrom2, pos2
        Genomic breakpoints (1-based positions; consistent with VCF POS).
    insertion
        Inserted sequence at the junction (may be empty; BEDPE loses this).
    pattern
        Integer 1-4 encoding the four VCF BND bracket patterns.
    sv_id
        Optional identifier (e.g., VCF ID, BEDPE name).
    strand1, strand2
        Optional BEDPE-style strand orientation (+/-). If missing we infer a
        consistent pair from `pattern`.
    """

    def __init__(self, chrom1, pos1, chrom2, pos2, insertion, pattern, sv_id=None, strand1=None, strand2=None):
        self.chrom1 = str(chrom1).replace('chr', '')
        self.pos1 = int(pos1)
        self.chrom2 = str(chrom2).replace('chr', '')
        self.pos2 = int(pos2)
        self.insertion = insertion
        self.pattern = int(pattern)
        self.sv_id = sv_id
        self._strand1 = strand1
        self._strand2 = strand2

    def __str__(self):
        return (
            "%s(chrom1=%s, pos1=%d, strand1=%s, chrom2=%s, pos2=%d, strand2=%s, insertion='%s', pattern=%d, sv_id=%r)"
            % (
                self.__class__.__name__,
                self.chrom1,
                self.pos1,
                self.strand1,
                self.chrom2,
                self.pos2,
                self.strand2,
                self.insertion,
                self.pattern,
                self.sv_id,
            )
        )

    def __repr__(self):
        return (
            "%s(%s, %d, %s, %d, '%s', %d, sv_id=%r, strand1=%r, strand2=%r)"
            % (
                self.__class__.__name__,
                self.chrom1,
                self.pos1,
                self.chrom2,
                self.pos2,
                self.insertion,
                self.pattern,
                self.sv_id,
                self._strand1,
                self._strand2,
            )
        )

    def __eq__(self, other):
        return self.sorted_coord == other.sorted_coord

    def __hash__(self):
        return hash(self.sorted_coord)

    @property
    def sorted_coord(self):
        """Sorted breakpoints (chr_pos strings) used for dedup identity."""
        bp1 = self.chrom1 + '_' + str(self.pos1)
        bp2 = self.chrom2 + '_' + str(self.pos2)
        return tuple(sorted([bp1, bp2]))

    @property
    def strand1(self):
        """Breakpoint 1 orientation in BEDPE-style (+/-).

        If missing, we infer from `pattern` using a mapping consistent with
        sv_utils.sv_pattern_infer_bedpe():

          pattern 1 -> (+, -)
          pattern 2 -> (+, +)
          pattern 3 -> (-, +)
          pattern 4 -> (-, -)
        """
        if self._strand1 in ['+', '-']:
            return self._strand1
        return '+' if self.pattern in (1, 2) else '-'

    @property
    def strand2(self):
        """Breakpoint 2 orientation in BEDPE-style (+/-)."""
        if self._strand2 in ['+', '-']:
            return self._strand2
        if self.pattern == 1:
            return '-'
        if self.pattern == 2:
            return '+'
        if self.pattern == 3:
            return '+'
        return '-'

    @property
    def svtype(self):
        """High-level SV category derived from chroms/positions/pattern."""
        if self.chrom1 != self.chrom2:
            return 'TRA'
        # same chromosome
        if self.pos1 < self.pos2:
            if self.pattern == 1:
                return 'DEL'
            if self.pattern == 2:
                return 'h2hINV'
            if self.pattern == 3:
                return 'DUP'
            return 't2tINV'
        # reversed ordering
        if self.pattern == 1:
            return 'DUP'
        if self.pattern == 2:
            return 'h2hINV'
        if self.pattern == 3:
            return 'DEL'
        return 't2tINV'
