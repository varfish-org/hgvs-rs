# Changelog

## [0.17.1](https://github.com/varfish-org/hgvs-rs/compare/v0.17.0...v0.17.1) (2024-07-26)


### Bug Fixes

* add selenoprotein biotype ([#197](https://github.com/varfish-org/hgvs-rs/issues/197)) ([e4a43b6](https://github.com/varfish-org/hgvs-rs/commit/e4a43b6e0027b819b243b2926bbdf026c0c77cd3))

## [0.17.0](https://github.com/varfish-org/hgvs-rs/compare/v0.16.3...v0.17.0) (2024-07-25)


### Features

* add variant region condition for dups ending in 3' UTR ([#192](https://github.com/varfish-org/hgvs-rs/issues/192)) ([#195](https://github.com/varfish-org/hgvs-rs/issues/195)) ([362909c](https://github.com/varfish-org/hgvs-rs/commit/362909ca0a7887cda6588f7da5671cb01ee4500f))

## [0.16.3](https://github.com/varfish-org/hgvs-rs/compare/v0.16.2...v0.16.3) (2024-07-16)


### Bug Fixes

* BioType naming for Ig* types, derive Hash for BioType ([#189](https://github.com/varfish-org/hgvs-rs/issues/189)) ([62519c7](https://github.com/varfish-org/hgvs-rs/commit/62519c7c6344b465fb65964aa52d503688b091a2))

## [0.16.2](https://github.com/varfish-org/hgvs-rs/compare/v0.16.1...v0.16.2) (2024-07-14)


### Bug Fixes

* add biotype vaultRNA_primary_transcript ([#180](https://github.com/varfish-org/hgvs-rs/issues/180)) ([d59783c](https://github.com/varfish-org/hgvs-rs/commit/d59783c51c30feeb36c389b2c8242e79ef9b5ad4))
* add biotype vaultRNA_primary_transcript ([#183](https://github.com/varfish-org/hgvs-rs/issues/183)) ([4345377](https://github.com/varfish-org/hgvs-rs/commit/4345377031220063a63e34be8b1af90c567e94f2))
* add partial field to cdot json representation struct ([#179](https://github.com/varfish-org/hgvs-rs/issues/179)) ([113a06a](https://github.com/varfish-org/hgvs-rs/commit/113a06a4646e6418b739da2423f953eadaac4ccc))
* add Serialize to data::cdot::json representations ([#176](https://github.com/varfish-org/hgvs-rs/issues/176)) ([2de1b8a](https://github.com/varfish-org/hgvs-rs/commit/2de1b8ab03d7671099b04e9bdfe58251f80e94f9))

## [0.16.1](https://github.com/varfish-org/hgvs-rs/compare/v0.16.0...v0.16.1) (2024-06-10)


### Bug Fixes

* use array instead of Vec for Codon representation ([#169](https://github.com/varfish-org/hgvs-rs/issues/169)) ([e455168](https://github.com/varfish-org/hgvs-rs/commit/e455168bf72b36f7e1c662bfa5dc3c10d2e448b0))


### Performance Improvements

* cache mappers, avoid allocations, codons are arrays ([#172](https://github.com/varfish-org/hgvs-rs/issues/172)) ([34acaf4](https://github.com/varfish-org/hgvs-rs/commit/34acaf4e5151de7186ce06ef42850ae6be716e18))

## [0.16.0](https://github.com/varfish-org/hgvs-rs/compare/v0.15.0...v0.16.0) (2024-03-05)


### Features

* adding vertebrate mitochondrial code ([#160](https://github.com/varfish-org/hgvs-rs/issues/160)) ([#161](https://github.com/varfish-org/hgvs-rs/issues/161)) ([eb6739b](https://github.com/varfish-org/hgvs-rs/commit/eb6739b4e19d00beb4bae60bbced6fa16f33cd06))

## [0.15.0](https://github.com/varfish-org/hgvs-rs/compare/v0.14.1...v0.15.0) (2024-02-08)


### Miscellaneous Chores

* bumping dependencies ([#157](https://github.com/varfish-org/hgvs-rs/issues/157)) ([6e75cd5](https://github.com/varfish-org/hgvs-rs/commit/6e75cd5dac6b5a885e682fdf57046497a8a45373))

## [0.14.1](https://github.com/varfish-org/hgvs-rs/compare/v0.14.0...v0.14.1) (2023-11-19)


### Bug Fixes

* allow import of BioType guide_RNA ([#148](https://github.com/varfish-org/hgvs-rs/issues/148)) ([0e1240f](https://github.com/varfish-org/hgvs-rs/commit/0e1240f3522d5dbf2b53e60ddc2df2d6f81d3b3a))

## [0.14.0](https://github.com/varfish-org/hgvs-rs/compare/v0.13.2...v0.14.0) (2023-11-19)


### Features

* adding support for selenoproteins ([#145](https://github.com/varfish-org/hgvs-rs/issues/145)) ([#146](https://github.com/varfish-org/hgvs-rs/issues/146)) ([c5e21e2](https://github.com/varfish-org/hgvs-rs/commit/c5e21e28aa9de2bd4738e1f52cd96930fb1a3e48))

## [0.13.2](https://github.com/varfish-org/hgvs-rs/compare/v0.13.1...v0.13.2) (2023-11-08)


### Bug Fixes

* add more biotypes from cdot JSON for vep110 release ([#141](https://github.com/varfish-org/hgvs-rs/issues/141)) ([4cb7cb6](https://github.com/varfish-org/hgvs-rs/commit/4cb7cb6bd4456d7341260e3f753ef58aa32aab48))

## [0.13.1](https://github.com/varfish-org/hgvs-rs/compare/v0.13.0...v0.13.1) (2023-11-08)


### Bug Fixes

* adding back missing Gene::biotype member ([#138](https://github.com/varfish-org/hgvs-rs/issues/138)) ([84db645](https://github.com/varfish-org/hgvs-rs/commit/84db6458b0e3551913488662a74dd738c1279308))

## [0.13.0](https://github.com/varfish-org/hgvs-rs/compare/v0.12.0...v0.13.0) (2023-11-08)


### Features

* make compatible with cdot v0.2.21 ([#136](https://github.com/varfish-org/hgvs-rs/issues/136)) ([c677b33](https://github.com/varfish-org/hgvs-rs/commit/c677b33d000117914359412e0f44c1317a82cc7c))

## [0.12.0](https://github.com/varfish-org/hgvs-rs/compare/v0.11.0...v0.12.0) (2023-10-21)


### Features

* move out assemblies to biocommons_bioutils crate ([#134](https://github.com/varfish-org/hgvs-rs/issues/134)) ([#135](https://github.com/varfish-org/hgvs-rs/issues/135)) ([d06bd28](https://github.com/varfish-org/hgvs-rs/commit/d06bd28d17884df1999669a378a516d317ab28e0))


### Bug Fixes

* problem with annotating stop_retained insertions ([#131](https://github.com/varfish-org/hgvs-rs/issues/131)) ([#132](https://github.com/varfish-org/hgvs-rs/issues/132)) ([0f42d20](https://github.com/varfish-org/hgvs-rs/commit/0f42d20fa89bfc3db429ce6bfe53984dd03ba29c))

## [0.11.0](https://github.com/varfish-org/hgvs-rs/compare/v0.10.1...v0.11.0) (2023-09-18)


### Features

* various dependency updates, bump version to 0.11.0 ([#128](https://github.com/varfish-org/hgvs-rs/issues/128)) ([09ed016](https://github.com/varfish-org/hgvs-rs/commit/09ed0167624ea67354eef78e698f9b0439050b11))

## [0.10.1](https://github.com/varfish-org/hgvs-rs/compare/v0.10.0...v0.10.1) (2023-07-04)


### Bug Fixes

* properly configure dependabot for noodles ([#119](https://github.com/varfish-org/hgvs-rs/issues/119)) ([b539921](https://github.com/varfish-org/hgvs-rs/commit/b5399215294e555df06bc78c7ffca4d71acb4c75))

## [0.10.0](https://github.com/varfish-org/hgvs-rs/compare/v0.9.0...v0.10.0) (2023-07-04)


### Miscellaneous Chores

* bump version to 0.10.0 ([7bf4cbd](https://github.com/varfish-org/hgvs-rs/commit/7bf4cbde1c4f09756943b339a7d28d05fc0b7e24))

## [0.9.0](https://github.com/varfish-org/hgvs-rs/compare/v0.8.2...v0.9.0) (2023-06-12)


### Features

* express existing thread safety (as most things are immutable) ([#115](https://github.com/varfish-org/hgvs-rs/issues/115)) ([0d6f241](https://github.com/varfish-org/hgvs-rs/commit/0d6f24177cfafa62796e9b5da069f8d1ca807aad))

## [0.8.2](https://github.com/varfish-org/hgvs-rs/compare/v0.8.1...v0.8.2) (2023-06-08)


### Code Refactoring

* replace linked-hash-map by indexmap ([#111](https://github.com/varfish-org/hgvs-rs/issues/111)) ([#112](https://github.com/varfish-org/hgvs-rs/issues/112)) ([f6d3ab4](https://github.com/varfish-org/hgvs-rs/commit/f6d3ab47dc79daab2ad837c5c739976061991926))

## [0.8.1](https://github.com/varfish-org/hgvs-rs/compare/v0.8.0...v0.8.1) (2023-05-23)


### Bug Fixes

* bump dependencies ([#108](https://github.com/varfish-org/hgvs-rs/issues/108)) ([af75b48](https://github.com/varfish-org/hgvs-rs/commit/af75b48b3f189010a9e6e27d8cd6a29477cb5198))

## [0.8.0](https://github.com/varfish-org/hgvs-rs/compare/v0.7.0...v0.8.0) (2023-05-23)


### Features

* losen dependencies ([#106](https://github.com/varfish-org/hgvs-rs/issues/106)) ([c654507](https://github.com/varfish-org/hgvs-rs/commit/c6545077e5e0bad33d4e0168cf455c3dcc1c928a))

## [0.7.0](https://github.com/varfish-org/hgvs-rs/compare/v0.6.2...v0.7.0) (2023-04-24)


### Features

* proper error handling with enums and thiserror ([#69](https://github.com/varfish-org/hgvs-rs/issues/69)) ([#103](https://github.com/varfish-org/hgvs-rs/issues/103)) ([add8248](https://github.com/varfish-org/hgvs-rs/commit/add8248cabe25d1f19993e7041b1e97279e18bdd))

## [0.6.2](https://github.com/varfish-org/hgvs-rs/compare/v0.6.1...v0.6.2) (2023-04-18)


### Bug Fixes

* issue with non-dup insertion at start of protein ([#99](https://github.com/varfish-org/hgvs-rs/issues/99)) ([#100](https://github.com/varfish-org/hgvs-rs/issues/100)) ([bc5b5cf](https://github.com/varfish-org/hgvs-rs/commit/bc5b5cf4cd77e56b96be879eb8c939d29e0c58e0))

## [0.6.1](https://github.com/varfish-org/hgvs-rs/compare/v0.6.0...v0.6.1) (2023-04-06)


### Bug Fixes

* cases where dup/inv goes beyond CDS ([#89](https://github.com/varfish-org/hgvs-rs/issues/89)) ([5d951b1](https://github.com/varfish-org/hgvs-rs/commit/5d951b131a1295fc6e83be2676834e5b1ea34244))
* only warn on combining RefAlt with whole gene deletion ([#92](https://github.com/varfish-org/hgvs-rs/issues/92)) ([bff6c72](https://github.com/varfish-org/hgvs-rs/commit/bff6c725ab10d56acb1f01809e7f4926f4a08d2e))
* out of bound panic in case of long deletions ([#91](https://github.com/varfish-org/hgvs-rs/issues/91)) ([f0bcaa6](https://github.com/varfish-org/hgvs-rs/commit/f0bcaa6a7aeaf70212cb7d6b8e49b5ace1aeeb92))
* out of bounds issue for protein sequence creation ([#93](https://github.com/varfish-org/hgvs-rs/issues/93)) ([a40a5f5](https://github.com/varfish-org/hgvs-rs/commit/a40a5f547406feefc9f952b209838914721857fd))
* out of bounds issue on 5'-to-3' shifting beyond CDS ([#94](https://github.com/varfish-org/hgvs-rs/issues/94)) ([4fccdb8](https://github.com/varfish-org/hgvs-rs/commit/4fccdb826fea15539b59630b2f9d742271aac3c4))
* problem with variants in multi-stop codon txs ([#95](https://github.com/varfish-org/hgvs-rs/issues/95)) ([#96](https://github.com/varfish-org/hgvs-rs/issues/96)) ([f25658c](https://github.com/varfish-org/hgvs-rs/commit/f25658c7ada6483970c46d0cf43bf53fd7a2b791))

## [0.6.0](https://github.com/varfish-org/hgvs-rs/compare/v0.5.2...v0.6.0) (2023-04-05)


### Features

* make some tables visible in `hgvs::sequences` ([#87](https://github.com/varfish-org/hgvs-rs/issues/87)) ([d81bf8c](https://github.com/varfish-org/hgvs-rs/commit/d81bf8c5ee8548471972828b2985fe99be94eccc))

## [0.5.2](https://github.com/varfish-org/hgvs-rs/compare/v0.5.1...v0.5.2) (2023-04-04)


### Performance Improvements

* further speeding up translate_cds code ([#83](https://github.com/varfish-org/hgvs-rs/issues/83)) ([#85](https://github.com/varfish-org/hgvs-rs/issues/85)) ([60b071d](https://github.com/varfish-org/hgvs-rs/commit/60b071db524df11b1f8e9633074df7c3213fb8ed))

## [0.5.1](https://github.com/varfish-org/hgvs-rs/compare/v0.5.0...v0.5.1) (2023-04-03)


### Performance Improvements

* tune translate_cds implementation ([#80](https://github.com/varfish-org/hgvs-rs/issues/80)) ([#81](https://github.com/varfish-org/hgvs-rs/issues/81)) ([a608ba1](https://github.com/varfish-org/hgvs-rs/commit/a608ba1a62892b9b49e85c06d02665db94c4e4ad))

## [0.5.0](https://github.com/varfish-org/hgvs-rs/compare/v0.4.0...v0.5.0) (2023-03-31)


### Features

* allow configuring that there is no genome sequence ([#65](https://github.com/varfish-org/hgvs-rs/issues/65)) ([cd5b7bb](https://github.com/varfish-org/hgvs-rs/commit/cd5b7bb0f04369b34f4d8857d558915ae2ddccbb))
* replacing usages of unwrap() with expect or Result ([#70](https://github.com/varfish-org/hgvs-rs/issues/70)) ([#73](https://github.com/varfish-org/hgvs-rs/issues/73)) ([94d6f88](https://github.com/varfish-org/hgvs-rs/commit/94d6f88f44f1f574dad7ca3015c2dba2852b869e))


### Bug Fixes

* case of missing stop codon (in particular ENSEMBL) ([#72](https://github.com/varfish-org/hgvs-rs/issues/72)) ([a44e28c](https://github.com/varfish-org/hgvs-rs/commit/a44e28c3ea57bea66cf9354c17b0f93e9a920636))
* fixing wrong validation assumption about insertions ([#64](https://github.com/varfish-org/hgvs-rs/issues/64)) ([f58ff53](https://github.com/varfish-org/hgvs-rs/commit/f58ff53ad017c65968afef64ad50c532be64f605))
* issue with transcripts missing stop codon ([#67](https://github.com/varfish-org/hgvs-rs/issues/67)) ([f70bb91](https://github.com/varfish-org/hgvs-rs/commit/f70bb91bba07f41110f32bc7da62f0cbc873e85d))
* return error if problem in normalization ([#71](https://github.com/varfish-org/hgvs-rs/issues/71)) ([b656d2a](https://github.com/varfish-org/hgvs-rs/commit/b656d2ae87e0cfbbb08a7ae6b055f2d9ce13f1bc))
* return error in normalization instead of unwrap() ([#68](https://github.com/varfish-org/hgvs-rs/issues/68)) ([8144db6](https://github.com/varfish-org/hgvs-rs/commit/8144db60b74bb5fa2ff625cfef73be2297ddfd0f))

## [0.4.0](https://github.com/varfish-org/hgvs-rs/compare/v0.3.1...v0.4.0) (2023-03-27)


### Features

* allow for disabling renormalization ([#53](https://github.com/varfish-org/hgvs-rs/issues/53)) ([b2f5dbe](https://github.com/varfish-org/hgvs-rs/commit/b2f5dbeb9904fd70494a030eab37bcf1edb8845a))
* improved normalizer configuration ([#52](https://github.com/varfish-org/hgvs-rs/issues/52)) ([2db5e92](https://github.com/varfish-org/hgvs-rs/commit/2db5e92316930e76d4b58d47ccc1d2ba2c1010fe))
* make sequences module public ([#54](https://github.com/varfish-org/hgvs-rs/issues/54)) ([cfbb134](https://github.com/varfish-org/hgvs-rs/commit/cfbb134b38d853273d734db7bdb7ce1beafa9ad8))
* more comprehensive clone for field-less enums ([#58](https://github.com/varfish-org/hgvs-rs/issues/58)) ([3ada6b2](https://github.com/varfish-org/hgvs-rs/commit/3ada6b2a21e24fea2feeb8fb89888710072d4b52))


### Bug Fixes

* cases of empty sequence ([#61](https://github.com/varfish-org/hgvs-rs/issues/61)) ([0a6a109](https://github.com/varfish-org/hgvs-rs/commit/0a6a1094f7387daa92c7bbcbb1c08e2c0c3620aa))
* fixing projection issues with cdot::json::Provider ([#55](https://github.com/varfish-org/hgvs-rs/issues/55)) ([edd0de9](https://github.com/varfish-org/hgvs-rs/commit/edd0de9e12e781a3fbeb8fc49d623b8d66484e8d))
* fixing wrong validation assumption about insertions ([#62](https://github.com/varfish-org/hgvs-rs/issues/62)) ([bf027dc](https://github.com/varfish-org/hgvs-rs/commit/bf027dc3cdef0b07af4cc234f79499b9146426cf))
* interpret replace_reference condition ([#49](https://github.com/varfish-org/hgvs-rs/issues/49)) ([8326a55](https://github.com/varfish-org/hgvs-rs/commit/8326a558eeee75defe8c466ec996ba02f81116cd))
* issue with out of bounds deletion ([#59](https://github.com/varfish-org/hgvs-rs/issues/59)) ([475398d](https://github.com/varfish-org/hgvs-rs/commit/475398d0e831afa43abc53327abb14135df40aea))
* remove ref length validation ([#57](https://github.com/varfish-org/hgvs-rs/issues/57)) ([763b7a3](https://github.com/varfish-org/hgvs-rs/commit/763b7a38f69f1c1b87472fa655099a988c424294))
* truncate AA insertion at stop codon ([#60](https://github.com/varfish-org/hgvs-rs/issues/60)) ([54231b9](https://github.com/varfish-org/hgvs-rs/commit/54231b983b1cbc4480631bde05ab55e2a668ed72))

## [0.3.1](https://github.com/varfish-org/hgvs-rs/compare/v0.3.0...v0.3.1) (2023-03-13)


### Bug Fixes

* deriving PartialEq and Eq for enums ([#47](https://github.com/varfish-org/hgvs-rs/issues/47)) ([abd946e](https://github.com/varfish-org/hgvs-rs/commit/abd946e0b37444222ff4f30da99eb61d67ac1a3d))

## [0.3.0](https://github.com/varfish-org/hgvs-rs/compare/v0.2.0...v0.3.0) (2023-03-13)


### Features

* implement cdot data provider ([#43](https://github.com/varfish-org/hgvs-rs/issues/43)) ([#44](https://github.com/varfish-org/hgvs-rs/issues/44)) ([3a8ed9d](https://github.com/varfish-org/hgvs-rs/commit/3a8ed9d49c1c34bb7295091afe82b7011d6826ef))

## [0.2.0](https://github.com/varfish-org/hgvs-rs/compare/v0.1.1...v0.2.0) (2023-03-06)


### Features

* implement PartialEq for Assembly ([#40](https://github.com/varfish-org/hgvs-rs/issues/40)) ([8808211](https://github.com/varfish-org/hgvs-rs/commit/8808211ed3f26c187f4d2787c23e680bbcdf38c3))


### Bug Fixes

* git lfs untrack 'src/static_data/**/*.json*' ([#41](https://github.com/varfish-org/hgvs-rs/issues/41)) ([fa1af68](https://github.com/varfish-org/hgvs-rs/commit/fa1af68c13b76bf1e9ba329159d3cf29b5893620))

## [0.1.1](https://github.com/varfish-org/hgvs-rs/compare/v0.1.0...v0.1.1) (2023-03-03)


### Bug Fixes

* release 0.1.1 on crates.io ([#38](https://github.com/varfish-org/hgvs-rs/issues/38)) ([e318abf](https://github.com/varfish-org/hgvs-rs/commit/e318abf6368b1f0b7160ad25f0880706a92fc662))

## 0.1.0 (2023-03-03)


### Features

* implement AlignmentMapper ([#14](https://github.com/varfish-org/hgvs-rs/issues/14)) ([#15](https://github.com/varfish-org/hgvs-rs/issues/15)) ([07b573d](https://github.com/varfish-org/hgvs-rs/commit/07b573df79601dab3bbb933258693f5afd55c55c))
* implement Display for HgvsVariant and related types ([#4](https://github.com/varfish-org/hgvs-rs/issues/4)) ([#6](https://github.com/varfish-org/hgvs-rs/issues/6)) ([8e81536](https://github.com/varfish-org/hgvs-rs/commit/8e815366b19639b932a837c916384e663524043a))
* implement parsing of HGVS variants ([#2](https://github.com/varfish-org/hgvs-rs/issues/2)) ([#3](https://github.com/varfish-org/hgvs-rs/issues/3)) ([dbcfd05](https://github.com/varfish-org/hgvs-rs/commit/dbcfd059802459d6a5ca595b560c955de2d5f4ac))
* implement VariantMapper ([#12](https://github.com/varfish-org/hgvs-rs/issues/12)) ([#13](https://github.com/varfish-org/hgvs-rs/issues/13)) ([44b10de](https://github.com/varfish-org/hgvs-rs/commit/44b10de9bc312814dfdb6ca160dedd5322570e0c))
* implementing AssemblyMapper ([#7](https://github.com/varfish-org/hgvs-rs/issues/7)) ([#35](https://github.com/varfish-org/hgvs-rs/issues/35)) ([9f3e0f3](https://github.com/varfish-org/hgvs-rs/commit/9f3e0f3198b5b78b394692d3a99fedc65091c467))
* port over access to the UTA data structures ([#10](https://github.com/varfish-org/hgvs-rs/issues/10)) ([#11](https://github.com/varfish-org/hgvs-rs/issues/11)) ([3e71e32](https://github.com/varfish-org/hgvs-rs/commit/3e71e3285215eed80e152d520864686d543af2b0))
* port over assembly info from bioutils ([#8](https://github.com/varfish-org/hgvs-rs/issues/8)) ([#9](https://github.com/varfish-org/hgvs-rs/issues/9)) ([d583556](https://github.com/varfish-org/hgvs-rs/commit/d5835565c358f2b132e94f5496a77c01a1b96096))
* port over variant normalizer ([#19](https://github.com/varfish-org/hgvs-rs/issues/19)) ([#20](https://github.com/varfish-org/hgvs-rs/issues/20)) ([8593068](https://github.com/varfish-org/hgvs-rs/commit/8593068608a045e753565cd8e17decb0f429b26e))
* porting over test_hgvs_variantmapper_gcp ([#21](https://github.com/varfish-org/hgvs-rs/issues/21)) ([#27](https://github.com/varfish-org/hgvs-rs/issues/27)) ([7f925c8](https://github.com/varfish-org/hgvs-rs/commit/7f925c845e8c2188d7cf06cb0beea75fa38c2ec3))
* porting over test_variantmapper_cp_real ([#21](https://github.com/varfish-org/hgvs-rs/issues/21)) ([#26](https://github.com/varfish-org/hgvs-rs/issues/26)) ([6873324](https://github.com/varfish-org/hgvs-rs/commit/687332422d317cb173733491080b460773b6b2f2))
* start to port over test_hgvs_grammer_full.py ([#21](https://github.com/varfish-org/hgvs-rs/issues/21)) ([#31](https://github.com/varfish-org/hgvs-rs/issues/31)) ([9d3e3b5](https://github.com/varfish-org/hgvs-rs/commit/9d3e3b532922095278abdb5ae57108f0c5067109))


### Bug Fixes

* variant mapper g_to_t used wrong logic ([#32](https://github.com/varfish-org/hgvs-rs/issues/32)) ([040e2d3](https://github.com/varfish-org/hgvs-rs/commit/040e2d3cd9ec77b2eecc755c4a8fc39a83f43101))

## Changelog
