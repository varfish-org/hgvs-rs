# Differences from biocommons/hgvs

hgvs-rs is a port of [biocommons/hgvs][bc] and aims to give the same results. This list covers the differences named in the pull requests, issues and code comments of hgvs-rs. It is not the result of a full comparison.

Each entry names the situation, what biocommons/hgvs does, and what hgvs-rs does and why. Links to the code of biocommons/hgvs point to commit [13289ab][bc-commit].

## Normalization

### Intronic variants

- **Situation:** `assembly::Mapper` with `normalize` returns a c. or n. variant with an intron offset, e.g. from `g_to_c()`.
- **biocommons/hgvs:** The normalizer raises `HGVSUnsupportedOperationError` for an intronic variant, and `AssemblyMapper` returns the variant unnormalized ([normalizer.py][bc-norm-intron], [assemblymapper.py][bc-am-norm]).
- **hgvs-rs:** It maps the variant to the genome, normalizes it there within its intron, and maps it back. The variant thus follows the 3' rule on the transcript, as in VariantValidator ([#257]).

## Genome data

### Genome alignment in `c_to_p`

- **Situation:** `c_to_p` of the assembly mapper on a transcript without an alignment to the assembly.
- **biocommons/hgvs:** `AssemblyMapper.c_to_p()` looks up the genome alignment for every variant and raises `HGVSDataNotAvailableError` if there is none ([assemblymapper.py][bc-am-c-to-p]).
- **hgvs-rs:** `assembly::Mapper::c_to_p()` looks up the alignment only for a variant with an intron offset, and only with `genome_seq_available`, a setting that biocommons/hgvs does not have. Only such a variant needs bases from the genome. As in biocommons/hgvs, the genome gives the bases of a dup such as `NM_003777.3:c.4110+1dup` ([#301]). Other variants work without the alignment.

## Protein prediction (`c_to_p`)

### Uncertain protein changes

- **Situation:** a protein change that `c_to_p` predicts, e.g. for `NM_000059.3:c.7791A>G`.
- **biocommons/hgvs:** By default, it marks an inferred protein change as uncertain (`inferred_p_is_uncertain`) and gives `NP_000050.2:p.(Lys2597=)` ([defaults.ini][bc-uncertain], [test][bc-am-test]).
- **hgvs-rs:** It always gives the certain form, `NP_000050.2:p.Lys2597=`. It has no global configuration such as `inferred_p_is_uncertain` (comment in `AltSeqToHgvsp::create_variant()` in [src/mapper/altseq.rs](src/mapper/altseq.rs)).

### CDS with an incomplete 3' end

- **Situation:** a CDS whose length is not a multiple of 3, e.g. of a GENCODE transcript with the tag `cds_end_NF`.
- **biocommons/hgvs:** `RefTranscriptData` raises `NotImplementedError`, so `c_to_p` fails for every variant on the transcript ([reftranscriptdata.py][bc-cds-len]). A mitochondrial CDS that ends with `T` or `TA` is the exception, as the poly(A) tail completes its stop codon.
- **hgvs-rs:** It translates the CDS up to its last complete codon, and a change after that codon gives `p.?`. It still rejects such a CDS if a stop codon is in frame, because then the frame changes, e.g. by ribosomal slippage in OAZ1 ([#299]).

### Reference with more than one stop codon

- **Situation:** the translation of the reference CDS has more than one stop codon.
- **biocommons/hgvs:** `c_to_p` gives `p.?` for every variant on the transcript ([altseqbuilder.py][bc-multi-stop]). A request for `p.Met1?` instead was closed as working as designed ([biocommons/hgvs#651][bc-651]).
- **hgvs-rs:** It still predicts the protein for a variant that changes no base from the first stop codon on, so that a change of Met1 gets a result ([#96]). A variant that changes the first stop codon or a base after it gives `p.?`. The check places the changed bases as far 3' as possible, so every description of a change gives the same result. E.g., if the CDS starts with `ATG GGG TAG`, both `c.2_6del` and `c.3_7del` give `p.?`. The prediction ends with the first stop codon. A frameshift in front of it changes none of its bases, so `c.3_4insG` gives the extension `p.Ter3ValextTer?`, as with one stop codon. A dup moves `cds_stop` with its copy, as an insertion does, so the alt protein keeps the stop codons after the copy.

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

### Stop codon in the inserted bases

- **Situation:** a frameshifting insertion or dup puts a stop codon into the reading frame within its own bases. E.g., `c.21_25dup` turns `CGC GTA GAA TAG` into `CGC GTA GCG TAG AAT AG`.
- **biocommons/hgvs:** `_insert_stop_seq_end()` finds the stop codon in the inserted bases, so the variant counts as in frame ([altseqbuilder.py][bc-insert-stop], [biocommons/hgvs#830][bc-830]). The example gives `p.(Glu9Ala)`.
- **hgvs-rs:** It decides the frame by the length of the insertion only. The example gives `p.Glu9AlafsTer2`. This is a known difference, not yet fixed.

### Silent change by a frameshifting insertion

- **Situation:** a frameshifting insertion or dup forms a new stop codon in front of the old one and keeps the protein. E.g., `c.22_26dup` turns `GTA GAA TAG` into `GTA GAG TAG ATA G`.
- **biocommons/hgvs:** `_insert_stop_seq_end()` finds the stop codon in the inserted bases ([altseqbuilder.py][bc-insert-stop]), and `c_to_p` gives `p.(Glu9=)`.
- **hgvs-rs:** It gives `p.Glu9GlufsTer?`. `AltSeqToHgvsp` keeps the frameshift flag for the silent change and checks this flag first. This is a known difference, not yet fixed.

### Dup of the first amino acids

- **Situation:** an in-frame insertion repeats the amino acids from Met1 on, right after them. E.g., `c.6_7insATGAAG` turns `ATG AAA GCG` into `ATG AAA ATG AAG GCG`, so `MKA` becomes `MKMKA`.
- **biocommons/hgvs:** It describes the change as a dup and gives `p.(Met1_Lys2dup)` ([altseq_to_hgvsp.py][bc-ins-dup]).
- **hgvs-rs:** It gives `p.Lys2_Ala3insMetLys`. `check_if_ins_is_dup()` in [src/mapper/altseq.rs](src/mapper/altseq.rs) skips a dup that starts at the first amino acid. HGVS describes such an insertion as a dup. This is a known difference, not yet fixed.

### Variants that leave the CDS intact

- **Situation:** an insertion, dup, deletion, delins or inversion that changes no base of the CDS, because it can go before c.1 or after the stop codon. Examples are `NM_022051.2:c.-1_1insGCC` and `NM_001051.2:c.1257dupG`, where c.1257 is the last base of the stop codon. If c.-6 is a T and c.1 an A, `c.-6_1inv` keeps c.1 and is the same change as `c.-5_-1inv`.
- **biocommons/hgvs:** It decides the region by the positions in the description ([`get_variant_region()`][bc-region]). Both examples count as exonic. `c_to_p` gives `NP_071334.1:p.Met1?` (`INITMET02` in [ext.tsv][bc-ext]) and `NP_001042.1:p.(=)` ([test][bc-dup-cds-end]). It cuts an inversion at c.1, so `c.-6_1inv` gives `p.Met1?`, while `c.-5_-1inv` gives `p.?`.
- **hgvs-rs:** It decides the region by the transcript sequence. If a placement before c.1 or after the stop codon gives the same sequence, the variant is a 5' or 3' UTR variant, and `c_to_p` gives `p.?` for both examples. So every description of the same change gives the same result, and no result places the variant in a start or stop codon that it does not change.

### Stop codon changed into another stop codon

- **Situation:** a variant changes the stop codon into another stop codon, e.g. `c.30G>A` turns `TAG` into `TAA`.
- **biocommons/hgvs:** It checks for a stop codon at the variant position before it checks for a silent change ([altseq_to_hgvsp.py][bc-stop-at-variant]). So it writes the change as a substitution, `p.(Ter10Ter)`.
- **hgvs-rs:** It writes `p.Ter10=`, as for any other silent change. HGVS writes a silent change with `=`, e.g. `p.Cys188=`. `AltSeqToHgvsp::create_variant()` in [src/mapper/altseq.rs](src/mapper/altseq.rs) turns such a substitution into `=` ([#218]).

### Silent change of the stop codon

- **Situation:** an insertion or dup changes the stop codon or forms a new one in front of it, so the protein stays the same. E.g., `NM_001253909.2:c.416_417insGTG` turns the stop codon `TAA` into `TAG TGA`, and `c.417A>G` turns it into `TAG`.
- **biocommons/hgvs:** The insertion moves `cds_stop` ([altseqbuilder.py][bc-delins-cds-stop]), so the alt protein ends with two stop codons. `c_to_p` reports the second one ([altseq_to_hgvsp.py][bc-stop-at-variant]) and gives `NP_001240838.1:p.(Ter140Ter)`, a position after the end of the protein. For `c.417A>G`, it gives `p.(Ter139Ter)`. A dup keeps `cds_stop` ([altseqbuilder.py][bc-dup-cds-stop]), so the result also depends on the description. If the CDS ends with `AAA TAG`, `c.26_28dup` gives `p.(Ter10Ter)`, and the same change `c.28_29insAAT` gives `p.(Ter11Ter)`. [biocommons/hgvs#895][bc-895] reports these problems, and [biocommons/hgvs#896][bc-896] proposes a fix.
- **hgvs-rs:** Translation ends at the first stop codon, so `c_to_p` counts such an alt protein as unchanged. It gives `p.Ter139=` for both variants on `NM_001253909.2`, and `p.Ter10=` for every placement of the dup. The position is the codon of the first changed base, as the 3' rule picks it. A dup moves `cds_stop` with its copy, as an insertion does.

### Variants after a CDS without stop codon

- **Situation:** the CDS has no stop codon, and a variant changes bases after its last base c.N, e.g., `c.N_*1insA`, `c.*1G>A` or `c.Ndup`.
- **biocommons/hgvs:** It decides the region by the positions in the description. A variant with a c.* end counts as a 3' UTR variant and gives `p.?`. A dup that ends at c.N fails with `IndexError`.
- **hgvs-rs:** Without a stop codon, the CDS goes on after c.N, so there is no 3' UTR. These variants change the CDS, and `c_to_p` gives `p.?`, as the amino acids after c.N are unknown.

### Alt protein of a CDS without stop codon

- **Situation:** the CDS has no stop codon, and a variant in it forms a new stop codon, or a stop codon follows c.N. E.g., `c.15_16insGGGTAA` turns `TGG GAA` at c.13 to c.18 into `TGG GGG TAA GAA`.
- **biocommons/hgvs:** It counts the last codon as a stop codon. `AltTranscriptData` ends the alt protein with the last stop codon up to `cds_stop`, or else with the next one after it ([altseqbuilder.py][bc-alt-cut]). The example gives `p.(Ter6GlyextTer1)`.
- **hgvs-rs:** The alt protein ends at c.N, as moved by the variant, because the amino acids after c.N are unknown. Otherwise, the alt protein would not line up with the reference protein. The example gives `p.Trp5_Glu6insGlyTer`, as with a stop codon after c.18. A frameshift that ends after c.N gives `fsTer?`.

[bc]: https://github.com/biocommons/hgvs
[bc-commit]: https://github.com/biocommons/hgvs/tree/13289ab3b9db3c751918906c543991c74917ff61
[bc-norm-intron]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/normalizer.py#L99-L101
[bc-am-norm]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/assemblymapper.py#L251-L262
[bc-am-c-to-p]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/assemblymapper.py#L169-L206
[bc-uncertain]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/_data/defaults.ini#L5
[bc-am-test]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_assemblymapper.py#L39-L43
[bc-norm-test]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_normalizer.py#L187-L188
[bc-cds-len]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/reftranscriptdata.py#L28-L39
[bc-multi-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L107-L108
[bc-651]: https://github.com/biocommons/hgvs/issues/651
[bc-del-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L272-L276
[bc-775]: https://github.com/biocommons/hgvs/pull/775
[bc-fs-start]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L83-L86
[bc-529]: https://github.com/biocommons/hgvs/issues/529
[bc-insert-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L361-L374
[bc-830]: https://github.com/biocommons/hgvs/pull/830
[bc-ins-dup]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L296-L315
[bc-region]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L205-L225
[bc-ext]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/data/ext.tsv#L12
[bc-dup-cds-end]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/tests/test_hgvs_variantmapper.py#L113-L117
[bc-delins-cds-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L255-L256
[bc-dup-cds-stop]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L302-L333
[bc-alt-cut]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseqbuilder.py#L58-L70
[bc-stop-at-variant]: https://github.com/biocommons/hgvs/blob/13289ab3b9db3c751918906c543991c74917ff61/src/hgvs/utils/altseq_to_hgvsp.py#L193-L203
[bc-895]: https://github.com/biocommons/hgvs/issues/895
[bc-896]: https://github.com/biocommons/hgvs/pull/896
[#96]: https://github.com/varfish-org/hgvs-rs/pull/96
[#218]: https://github.com/varfish-org/hgvs-rs/pull/218
[#257]: https://github.com/varfish-org/hgvs-rs/pull/257
[#299]: https://github.com/varfish-org/hgvs-rs/pull/299
[#301]: https://github.com/varfish-org/hgvs-rs/pull/301
