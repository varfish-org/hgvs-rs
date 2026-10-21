# Differences from biocommons/hgvs

hgvs-rs is a port of [biocommons/hgvs][bc] and aims to give the same results. This list covers the differences named in the pull requests, issues and code comments of hgvs-rs. It is not the result of a full comparison.

Each entry names the situation, what biocommons/hgvs does, and what hgvs-rs does and why. The entries under [Proposed fixes for biocommons/hgvs](#proposed-fixes-for-biocommonshgvs) are deviations that hgvs-rs makes on purpose and that are to be proposed upstream. Links to the code of biocommons/hgvs point to commit [13289ab][bc-commit]. The examples on mock transcripts were run with both libraries.

## Normalization

### Intronic variants

- **Situation:** `assembly::Mapper` with `normalize` returns a c. or n. variant with an intron offset, e.g. from `g_to_c()`.
- **biocommons/hgvs:** The normalizer raises `HGVSUnsupportedOperationError` for an intronic variant, and `AssemblyMapper` returns the variant unnormalized ([normalizer.py][bc-norm-intron], [assemblymapper.py][bc-am-norm]).
- **hgvs-rs:** It maps the variant to the genome, normalizes it there within its intron, and maps it back. The variant thus follows the 3' rule on the transcript, as in VariantValidator ([#257]).

## Genome data

### Genome alignment in `c_to_p`

- **Situation:** `c_to_p` of the assembly mapper on a transcript without an alignment to the assembly.
- **biocommons/hgvs:** `AssemblyMapper.c_to_p()` looks up the genome alignment for every variant and raises `HGVSDataNotAvailableError` if there is none ([assemblymapper.py][bc-am-c-to-p]).
- **hgvs-rs:** `assembly::Mapper::c_to_p()` looks up the alignment only for a variant with an intron offset, and only with `genome_seq_available`, a setting that biocommons/hgvs does not have. Only such a variant needs bases from the genome. As in biocommons/hgvs, the genome gives the bases of a dup such as `NM_003777.3:c.4110+1dup` ([#301]). Other variants work without the alignment. With `shift_over_boundary` and the preference `Intron` or `Exon`, it also looks up the alignment for an insertion or dup in an exon or intron, as only these shift.

## Protein prediction (`c_to_p`)

### Uncertain protein changes

- **Situation:** a protein change that `c_to_p` predicts, e.g. for `NM_000059.3:c.7791A>G`.
- **biocommons/hgvs:** By default, it marks an inferred protein change as uncertain (`inferred_p_is_uncertain`) and gives `NP_000050.2:p.(Lys2597=)` ([defaults.ini][bc-uncertain], [test][bc-am-test]).
- **hgvs-rs:** It always gives the certain form, `NP_000050.2:p.Lys2597=`. It has no global configuration such as `inferred_p_is_uncertain` (comment in `AltSeqToHgvsp::create_variant()` in [src/mapper/altseq.rs](src/mapper/altseq.rs)).

### CDS with an incomplete 3' end

- **Situation:** a CDS whose length is not a multiple of 3, e.g. of a GENCODE transcript with the tag `cds_end_NF`.
- **biocommons/hgvs:** `RefTranscriptData` raises `NotImplementedError`, so `c_to_p` fails for every variant on the transcript ([reftranscriptdata.py][bc-cds-len]). A mitochondrial CDS that ends with `T` or `TA` is the exception, as the poly(A) tail completes its stop codon.
- **hgvs-rs:** It translates the CDS up to its last complete codon, and a change after that codon gives `p.?`. It still rejects such a CDS if a stop codon is in frame, because then the frame changes, e.g. by ribosomal slippage in OAZ1 ([#299]). The mitochondrial exception is the same as in biocommons/hgvs.

### Variants after a CDS without stop codon

- **Situation:** the CDS has no stop codon, and a variant ends with a position after its last base c.N, e.g. `c.N_*1insA`, `c.*1G>A` or `c.(N-1)_*1del`.
- **biocommons/hgvs:** It counts the last codon as a stop codon. A variant with a c.* start, and an insertion or dup with a c.* end, is a 3' UTR variant and gives `p.?` ([`get_variant_region()`][bc-region]). Other variants that end with a c.* position get an extension. If the CDS ends with `GCG AAA`, followed by `GGGGG`, `c.26_*1del` gives `p.(Ter9ArgextTer?)`. A dup that ends at c.N fails with `IndexError` if bases follow c.N.
- **hgvs-rs:** Without a stop codon, the CDS goes on after c.N, so there is no 3' UTR. These variants are CDS variants, and `c_to_p` gives `p.?`, as the amino acids after c.N are unknown. The same holds for a variant that can go after c.N, e.g. `c.Ndup`. The example gives `p.?`.

### Alt protein of a CDS without stop codon

- **Situation:** the CDS has no stop codon, and a variant in it forms a new stop codon, or a stop codon follows c.N. E.g., `c.15_16insGGGTAA` turns `TGG GAA` at c.13 to c.18 into `TGG GGG TAA GAA`.
- **biocommons/hgvs:** It counts the last codon as a stop codon. `AltTranscriptData` ends the alt protein with the last stop codon up to `cds_stop`, or else with the next one after it ([altseqbuilder.py][bc-alt-cut]). The example gives `p.(Ter6GlyextTer1)`.
- **hgvs-rs:** The alt protein ends at c.N, as moved by the variant, because the amino acids after c.N are unknown. Otherwise, the alt protein would not line up with the reference protein. The example gives `p.Trp5_Glu6insGlyTer`, as with a stop codon after c.18. A frameshift that ends after c.N gives `fsTer?`. For the same reason, a dup moves `cds_stop` with its copy, as an insertion does.

### Deletion with more reference bases than positions

- **Situation:** a deletion whose reference is longer than its interval, e.g. `NM_000059.3:c.7790delAAG`.
- **biocommons/hgvs:** `c_to_p` and the normalizer raise `HGVSInvalidVariantError` ([tests][bc-am-test], [normalizer tests][bc-norm-test]).
- **hgvs-rs:** It does not check this yet. `c_to_p` gives `NP_000050.2:p.Glu2598LysfsTer50`, and the normalizer gives `NM_000059.3:c.7791delA` ([#96]). This difference is a gap, not a decision.

### Deletion of the stop codon

- **Situation:** an in-frame deletion runs from the CDS through the stop codon into the 3' UTR. E.g., `c.26_*1del` turns `GCG AAA TAG AGG AGG TAG` into `GCG AGG AGG TAG`.
- **biocommons/hgvs:** It treats a deletion as a frameshift if its start is a c. position and its end a c.* position ([altseqbuilder.py][bc-del-stop], [biocommons/hgvs#775][bc-775]). HGVS describes such deletions as frameshifts. The example gives `p.(Lys9ArgfsTer3)`.
- **hgvs-rs:** It decides the frame by the length of the deletion only. The example gives `p.Lys9delinsArgArg`, which hides the lost stop codon. This is a known difference, not yet fixed.

### Deletion of the stop codon up to the transcript end

- **Situation:** a deletion removes the stop codon and every base after it, so no codon takes the place of the stop codon. E.g., `c.28_*3del` on a transcript that ends at c.*3.
- **biocommons/hgvs:** `c_to_p` raises `IndexError`. `AltSeqToHgvsp` searches for the first changed amino acid without a bound, but the alt protein is the reference protein without its stop ([altseq_to_hgvsp.py][bc-fs-start]). [biocommons/hgvs#529][bc-529] reported a similar `IndexError`.
- **hgvs-rs:** It gives `p.Ter10extTer?`. This extension lacks the new amino acid, as the alt transcript has none. This is a known difference, not yet fixed.

### Dup of the first amino acids

- **Situation:** an in-frame insertion repeats the amino acids from Met1 on, right after them. E.g., `c.6_7insATGAAG` turns `ATG AAA GCG` into `ATG AAA ATG AAG GCG`, so `MKA` becomes `MKMKA`.
- **biocommons/hgvs:** It describes the change as a dup and gives `p.(Met1_Lys2dup)` ([altseq_to_hgvsp.py][bc-ins-dup]).
- **hgvs-rs:** It gives `p.Lys2_Ala3insMetLys`. `check_if_ins_is_dup()` in [src/mapper/altseq.rs](src/mapper/altseq.rs) skips a dup that starts at the first amino acid. HGVS describes such an insertion as a dup. This is a known difference, not yet fixed.

### Stop codon changed into another stop codon

- **Situation:** a variant changes the stop codon into another stop codon, e.g. `c.30G>A` turns `TAG` into `TAA`.
- **biocommons/hgvs:** It checks for a stop codon at the variant position before it checks for a silent change ([altseq_to_hgvsp.py][bc-stop-at-variant]). So it writes the change as a substitution, `p.(Ter10Ter)`.
- **hgvs-rs:** It writes `p.Ter10=`, as for any other silent change. HGVS writes a silent change with `=`, e.g. `p.Cys188=`. `AltSeqToHgvsp::create_variant()` in [src/mapper/altseq.rs](src/mapper/altseq.rs) turns such a substitution into `=` ([#218]).

### New stop codon in front of the old one

- **Situation:** an in-frame insertion forms a new stop codon right in front of the old one, so the protein stays the same. E.g., if the CDS ends with `AAA TAG`, `c.28_29insAAT` gives `AAA TAA TAG`.
- **biocommons/hgvs:** The insertion moves `cds_stop` ([altseqbuilder.py][bc-delins-cds-stop]), so the alt protein ends with both stop codons. `c_to_p` reports the second one ([altseq_to_hgvsp.py][bc-stop-at-variant]) and gives `p.(Ter11Ter)`, a position after the end of the protein. A dup keeps `cds_stop` ([altseqbuilder.py][bc-dup-cds-stop]), so the same change `c.26_28dup` gives `p.(Ter10Ter)`. [biocommons/hgvs#895][bc-895] reports this, and [biocommons/hgvs#896][bc-896] proposes a fix.
- **hgvs-rs:** It gives `p.=` for a change after the end of the reference protein (`AltSeqToHgvsp::convert_to_hgvs_variant()` in [src/mapper/altseq.rs](src/mapper/altseq.rs)), so the insertion gives `p.=`. The dup gives `p.Ter10=`, see [Stop codon changed into another stop codon](#stop-codon-changed-into-another-stop-codon). E.g., `NM_001253909.2:c.416_417insGTG` turns the stop codon `TAA` into `TAG TGA` and gives `p.=`.

## Proposed fixes for biocommons/hgvs

Unless an entry says otherwise, the examples use a transcript whose 5' UTR ends with `AAAATCAAAC`, whose CDS is `ATG AAA GCG AAA GCG TTT CGC GCG AAA TAG` (`MKAKAFRAK*`), and whose 3' UTR starts with `GGGGGCCCCC`, in a single exon. So c.-6 is a T, c.-1 a C, c.30 a G and c.*1 a G.

### Variants that leave the CDS intact

- **Situation:** an insertion, dup, deletion, delins or inversion whose alt transcript sequence keeps the start codon or the stop codon where it is, because the same change fits in front of c.1 or after the stop codon. Every placement of such a change gives the same alt sequence, so they cannot differ in their effect.
- **biocommons/hgvs:** It decides the region by the positions in the description ([`get_variant_region()`][bc-region]). A description that touches c.1 or the stop codon counts as exonic, and the protein is built from it literally. `c.-1_1insGCC`, `c.1dup`, `c.-1_1del` and `c.-6_1inv` give `p.Met1?`, and `c.1_3dup` gives `p.(Met1dup)`, although the ATG stays intact. `c.30dup` and `c.28_30dup` give `p.(=)`, and `c.30del` and `c.29_30insG` give `p.(Ter10Ter)`, although the stop codon is unchanged. Likewise, `NM_022051.2:c.-1_1insGCC` gives `NP_071334.1:p.Met1?` (`INITMET02` in [ext.tsv][bc-ext]), and `NM_001051.2:c.1257dupG`, where c.1257 is the last base of the stop codon, gives `NP_001042.1:p.(=)` ([test][bc-dup-cds-end]). With `ins_at_boundary_is_intronic` set to `False` and c.1 at an exon start, `c.1-1_1insC` gives `p.Met1?`. A dup from far into the 5' UTR, e.g. `c.-21_2dup` on a 5' UTR of 25 bases, leaves no ATG at the old start index and gives `p.Met1?`, while the intact CDS just moves 23 bases downstream.
- **hgvs-rs (proposed fix):** It checks the alt transcript sequence. If it ends with the reference from c.1 on, or starts with the reference up to the last base of the CDS, the variant is a 5' or 3' UTR variant, and `c_to_p` gives `p.?`, as for any UTR variant. All of the examples give `p.?`, and every description of the same change gives the same result. Without a stop codon, there is no 3' UTR, see [Variants after a CDS without stop codon](#variants-after-a-cds-without-stop-codon).

### Positions after a completed mitochondrial stop codon

- **Situation:** a variant with a c.* position on a mitochondrial transcript whose stop codon the poly(A) tail completes, e.g. `c.16_*1del` on a CDS `ATG GCC AAG CTG TGG T` that is followed by `CCCC`.
- **biocommons/hgvs:** `RefTranscriptData` puts the completing `A`s into the transcript sequence after the CDS ([reftranscriptdata.py][bc-cds-len]), but `cds_stop` stays, and `AltSeqBuilder` counts the c.* positions from `cds_stop` ([altseqbuilder.py][bc-setup]). So c.*1 is the first completing `A`, not the first base after the CDS. The example deletes `TA` instead of `TC` and gives `p.(Ter6ThrextTer?)`.
- **hgvs-rs (proposed fix):** The completing `A`s go into the translated bases only, so the c.* positions count the bases of the transcript. The example gives `p.Ter6ProextTer?`.

## Shifts over exon boundaries (`shift_over_boundary`)

### Dup rewritten as an insertion

- **Situation:** a shift moves a dup fully 3' or 5' on the genome and rewrites it as an insertion of its bases, e.g. `NM_000314.4:c.209+1_209+2dup` (PTEN, plus strand) with `Exon`.
- **biocommons/hgvs:** Its docstring says that the rewrite shifts the variant "farther than would normally be possible". The code inserts before the duplicated bases for a 3' shift on the plus strand and for a 5' shift on the minus strand, else after them ([`_var_g_shift_with_rewrite()`][bc-shift-rewrite]). The shifts work on the genome, where the strand plays no role, so on the plus strand the insertion sits at the near end of the shift. The 3' shift of the example becomes `c.209_209+1insGT` and the 5' shift `c.209+1_209+2insTG`. Both are intronic, so the example gives `p.?`. Likewise, `NM_000314.4:c.802dup` with `Exon` becomes `c.802-1_802insG` and gives `p.?`, and `NM_000314.4:c.803_810dup` with `Intron` never reaches the intron and gives `p.(Phe271ThrfsTer8)`.
- **hgvs-rs:** It inserts at the end that the shift reached, on both strands, as the docstring describes. The examples become `c.208_209insTG` (`p.Cys71PhefsTer29`), `c.802_803insG` (`p.Asp268GlyfsTer30`) and `c.802-1_802insGACAAAAT` (`p.?`). On the minus strand, both rules give the same insertion, e.g. `NM_004985.4:c.111+1dup` (KRAS) with `Exon` gives `p.Asp38GlyfsTer10` in both libraries. These results come from running biocommons/hgvs on the test data of hgvs-rs.

### `Default` skips the shifts

- **Situation:** `shift_over_boundary` with the preference `Default` on an insertion or dup in an exon or intron.
- **biocommons/hgvs:** It computes the shifts but never keeps one, so the protein is the same as without the option. Without `alt_ac`, `VariantMapper.c_to_p()` still raises `HGVSUnsupportedOperationError` ([variantmapper.py][bc-shift-loop]).
- **hgvs-rs:** It skips the shifts. So `Default` needs no genome, and `variant::Mapper::c_to_p()` does not fail.

### Errors during a shift

- **Situation:** a shift fails.
- **biocommons/hgvs:** It skips the shift on `HGVSInvalidVariantError`, `HGVSInvalidIntervalError` or `HGVSUnsupportedOperationError` and raises any other error ([`_var_c_shifts()`][bc-shift-errors]).
- **hgvs-rs:** Its errors do not map one to one to these. It skips the shift on any error except missing data, and it returns `DataError` and the `ReplaceReferenceFailed` of the normalizer.

### Translation table of a shift

- **Situation:** the protein comes from a shifted variant on a transcript with a translation table other than the standard one.
- **biocommons/hgvs:** It builds the `AltSeqBuilder` of the shifted variant without the translation table, so it translates the shifted sequence with the standard table ([variantmapper.py][bc-shift-loop]).
- **hgvs-rs:** It uses the reference data of the transcript, including its translation table.

### Insertions and dups that leave the CDS intact

- **Situation:** `Intron` or `Exon` on an insertion or dup that fits in front of c.1 or after the stop codon, e.g. `c.1dup` if c.-1 is not an A. See [Variants that leave the CDS intact](#variants-that-leave-the-cds-intact).
- **biocommons/hgvs:** It counts the variant as exonic or intronic by its description and shifts it ([`get_variant_region()`][bc-region]).
- **hgvs-rs:** The variant is a 5' or 3' UTR variant, which lies neither in an exon nor in an intron. `Intron` and `Exon` keep it, which gives `p.?`. They also skip a shift that ends in such a placement.

[bc]: https://github.com/biocommons/hgvs
[bc-commit]: https://github.com/biocommons/hgvs/tree/13289ab3b9db3c751918906c543991c74917ff61
[bc-norm-intron]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/normalizer.py#L99-L101
[bc-am-norm]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/assemblymapper.py#L251-L262
[bc-am-c-to-p]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/assemblymapper.py#L169-L206
[bc-uncertain]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/_data/defaults.ini#L5
[bc-am-test]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_assemblymapper.py#L39-L43
[bc-norm-test]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_normalizer.py#L187-L188
[bc-cds-len]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/reftranscriptdata.py#L28-L39
[bc-setup]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L376-L413
[bc-region]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L205-L225
[bc-ext]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/data/ext.tsv#L12
[bc-dup-cds-end]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_variantmapper.py#L113-L117
[bc-alt-cut]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L58-L70
[bc-del-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L272-L276
[bc-775]: https://github.com/biocommons/hgvs/pull/775
[bc-fs-start]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L83-L86
[bc-529]: https://github.com/biocommons/hgvs/issues/529
[bc-ins-dup]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L296-L315
[bc-stop-at-variant]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L193-L203
[bc-delins-cds-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L255-L256
[bc-dup-cds-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L302-L333
[bc-895]: https://github.com/biocommons/hgvs/issues/895
[bc-896]: https://github.com/biocommons/hgvs/pull/896
[bc-shift-loop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/variantmapper.py#L508-L534
[bc-shift-errors]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/variantmapper.py#L971-L992
[bc-shift-rewrite]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/variantmapper.py#L1013-L1045
[#96]: https://github.com/varfish-org/hgvs-rs/pull/96
[#218]: https://github.com/varfish-org/hgvs-rs/pull/218
[#257]: https://github.com/varfish-org/hgvs-rs/pull/257
[#299]: https://github.com/varfish-org/hgvs-rs/pull/299
[#301]: https://github.com/varfish-org/hgvs-rs/pull/301
