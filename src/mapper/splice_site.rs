//! Placements of a variant at an exon boundary, for
//! `ShiftOverBoundaryPreference::KeepSpliceSite`.
//!
//! An insertion or duplication in a repeat fits at several equivalent places on the genome, its
//! placements.  Positions are 1-based genome positions.

use crate::mapper::alignment;
use crate::mapper::cigar::CigarOp;

/// The region that a placement lies in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Region {
    /// In an exon.  An insertion at the edge of an exon counts as in the exon.
    Exon,
    /// In an intron, keeping its essential splice-site bases.
    Intron,
}

/// A placement of a variant on the genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Placement {
    /// An insertion between `g.x` and `g.x+1`.
    Insertion(i32),
    /// The bases from `g.start` to `g.end`, e.g. the bases that a duplication copies.
    Bases(i32, i32),
}

/// All equivalent placements of a variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Placements {
    /// Insertions between `g.x` and `g.x+1` for `x` in `lo..=hi`.
    Insertions { lo: i32, hi: i32 },
}

/// The exons of a transcript on the genome.
///
/// The essential splice-site bases are the two bases at each end of an intron.  On the plus
/// strand, the start of an intron has `+1` and `+2` of the donor, and its end has `-2` and `-1`
/// of the acceptor.  The minus strand swaps donor and acceptor, but the bases are the same.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Exons {
    /// `(start, end)` of each exon, inclusive, in genome order.
    spans: Vec<(i32, i32)>,
    /// Whether the transcript is on the plus strand.
    plus_strand: bool,
}

impl Exons {
    #[cfg(test)]
    fn new(spans: Vec<(i32, i32)>, plus_strand: bool) -> Self {
        Self { spans, plus_strand }
    }

    /// The exons of the alignment of `mapper`, i.e., its parts between introns.
    pub(crate) fn from_alignment(mapper: &alignment::Mapper) -> Self {
        let cigar = &mapper.cigar_mapper;
        let mut spans = Vec::new();
        let mut start = None;
        for (i, op) in cigar.cigar_op.iter().enumerate() {
            if *op == CigarOp::Skip {
                if let Some(start) = start.take() {
                    spans.push((start, cigar.ref_pos[i]));
                }
            } else if start.is_none() {
                start = Some(cigar.ref_pos[i]);
            }
        }
        if let (Some(start), Some(end)) = (start, cigar.ref_pos.last()) {
            spans.push((start, *end));
        }
        Self {
            // From 0-based, right-open positions relative to the alignment start.
            spans: spans
                .into_iter()
                .filter(|(start, end)| start < end)
                .map(|(start, end)| (start + mapper.gc_offset + 1, end + mapper.gc_offset))
                .collect(),
            plus_strand: mapper.strand != -1,
        }
    }

    /// The introns as `(start, end)`, inclusive.
    fn introns(&self) -> impl Iterator<Item = (i32, i32)> + '_ {
        self.spans
            .windows(2)
            .map(|pair| (pair[0].1 + 1, pair[1].0 - 1))
    }

    /// Insertion points between two bases of an exon, as inclusive ranges.
    fn exon_points(&self) -> Vec<(i32, i32)> {
        self.spans
            .iter()
            .map(|(start, end)| (*start, end - 1))
            .collect()
    }

    /// Insertion points at the edge of an exon next to an intron, as inclusive ranges.
    fn edge_points(&self) -> Vec<(i32, i32)> {
        self.introns()
            .flat_map(|(start, end)| [(start - 1, start - 1), (end, end)])
            .collect()
    }

    /// Insertion points in an intron that keep its essential splice-site bases, as inclusive
    /// ranges.  They exclude the points between the two bases at each end.
    fn intron_points(&self) -> Vec<(i32, i32)> {
        self.introns()
            .map(|(start, end)| (start + 1, end - 2))
            .collect()
    }

    /// The region of `placement`, if it lies fully in an exon or keeps the splice site in an
    /// intron.
    pub(crate) fn region_of(&self, placement: Placement) -> Option<Region> {
        let contains = |ranges: &[(i32, i32)], start: i32, end: i32| {
            ranges.iter().any(|(lo, hi)| *lo <= start && end <= *hi)
        };
        match placement {
            Placement::Insertion(x) => {
                if contains(&self.intron_points(), x, x) {
                    Some(Region::Intron)
                } else if contains(&self.exon_points(), x, x) || contains(&self.edge_points(), x, x)
                {
                    Some(Region::Exon)
                } else {
                    None
                }
            }
            Placement::Bases(start, end) => {
                let cores: Vec<_> = self
                    .introns()
                    .map(|(start, end)| (start + 2, end - 2))
                    .collect();
                if contains(&cores, start, end) {
                    Some(Region::Intron)
                } else if contains(&self.spans, start, end) {
                    Some(Region::Exon)
                } else {
                    None
                }
            }
        }
    }

    /// Choose the placement that `c_to_p` uses for the protein, given the placement `given` of
    /// the variant and all its placements `all`.  `None` means that every placement changes an
    /// essential splice-site base.
    ///
    /// A placement in an intron that keeps the splice site goes first, then one in an exon.
    /// `given` wins in its region, else the most 3' placement on the transcript.  In an exon, an
    /// insertion between two exon bases goes before one at the edge.
    pub(crate) fn choose(&self, given: Placement, all: Placements) -> Option<(Region, Placement)> {
        if self.region_of(given) == Some(Region::Intron) {
            return Some((Region::Intron, given));
        }
        if let Some(placement) = self.most_3_prime(all, &self.intron_points()) {
            return Some((Region::Intron, placement));
        }
        if self.region_of(given) == Some(Region::Exon) {
            return Some((Region::Exon, given));
        }
        self.most_3_prime(all, &self.exon_points())
            .or_else(|| self.most_3_prime(all, &self.edge_points()))
            .map(|placement| (Region::Exon, placement))
    }

    /// The most 3' placement of `all` on the transcript within `ranges`.
    fn most_3_prime(&self, all: Placements, ranges: &[(i32, i32)]) -> Option<Placement> {
        let Placements::Insertions { lo, hi } = all;
        let candidates = ranges
            .iter()
            .filter_map(|(start, end)| {
                let (start, end) = (lo.max(*start), hi.min(*end));
                (start <= end).then_some((start, end))
            })
            .map(|(start, end)| if self.plus_strand { end } else { start });
        let x = if self.plus_strand {
            candidates.max()
        } else {
            candidates.min()
        }?;
        Some(Placement::Insertion(x))
    }
}

/// `seq` for a placement `by` positions further 5' on the genome, i.e., rotated by `by`.
///
/// An insertion of `seq` between `g.x` and `g.x+1` equals the insertion of
/// `rotated(seq, by)` between `g.(x-by)` and `g.(x-by+1)` if both lie in the same repeat.
pub(crate) fn rotated(seq: &str, by: i32) -> String {
    let Ok(len) = i32::try_from(seq.len()) else {
        return seq.to_string();
    };
    if len == 0 || !seq.is_ascii() {
        return seq.to_string();
    }
    let split = (len - by.rem_euclid(len)) as usize;
    format!("{}{}", &seq[split..], &seq[..split])
}

#[cfg(test)]
mod test {
    use super::{rotated, Exons, Placement, Placements, Region};

    /// Exons `g.100_200` and `g.301_400` with the intron `g.201_300` between them.
    fn exons(plus_strand: bool) -> Exons {
        Exons::new(vec![(100, 200), (301, 400)], plus_strand)
    }

    #[test]
    fn region_of_insertion() {
        // The same on both strands.  On the plus strand, `g.201_202` are `+1` and `+2` and
        // `g.299_300` are `-2` and `-1`.  On the minus strand, `g.299_300` are `+2` and `+1`.
        for plus_strand in [true, false] {
            let exons = exons(plus_strand);
            for (x, expected) in [
                (150, Some(Region::Exon)),
                // At the exon edge: before `+1` and after `-1`.
                (200, Some(Region::Exon)),
                (300, Some(Region::Exon)),
                // Between the two essential bases.
                (201, None),
                (299, None),
                // In the intron next to the essential bases.
                (202, Some(Region::Intron)),
                (298, Some(Region::Intron)),
                (250, Some(Region::Intron)),
            ] {
                assert_eq!(exons.region_of(Placement::Insertion(x)), expected, "{x}");
            }
        }
    }

    #[test]
    fn region_of_bases() {
        let exons = exons(true);
        for (start, end, expected) in [
            (150, 200, Some(Region::Exon)),
            (301, 301, Some(Region::Exon)),
            (200, 201, None),
            (201, 201, None),
            (202, 202, None),
            (203, 298, Some(Region::Intron)),
            (299, 299, None),
            (300, 301, None),
        ] {
            assert_eq!(
                exons.region_of(Placement::Bases(start, end)),
                expected,
                "{start}-{end}"
            );
        }
    }

    #[test]
    fn choose_insertion() {
        let plus = exons(true);
        let minus = exons(false);
        // In the intron after the essential bases: intronic, whatever the given placement.
        let all = Placements::Insertions { lo: 199, hi: 203 };
        for exons in [&plus, &minus] {
            assert_eq!(
                exons.choose(Placement::Insertion(201), all),
                Some((
                    Region::Intron,
                    Placement::Insertion(if exons.plus_strand { 203 } else { 202 })
                ))
            );
        }
        // Only up to the essential bases: in the exon, before the edge.
        let all = Placements::Insertions { lo: 198, hi: 201 };
        assert_eq!(
            plus.choose(Placement::Insertion(201), all),
            Some((Region::Exon, Placement::Insertion(199)))
        );
        assert_eq!(
            minus.choose(Placement::Insertion(201), all),
            Some((Region::Exon, Placement::Insertion(198)))
        );
        // The given placement wins in its region.
        assert_eq!(
            plus.choose(Placement::Insertion(198), all),
            Some((Region::Exon, Placement::Insertion(198)))
        );
        // A duplication of an exon base counts as in the exon.
        assert_eq!(
            plus.choose(Placement::Bases(200, 200), all),
            Some((Region::Exon, Placement::Bases(200, 200)))
        );
        // Only the edge keeps the splice site.
        let all = Placements::Insertions { lo: 299, hi: 300 };
        for exons in [&plus, &minus] {
            assert_eq!(
                exons.choose(Placement::Insertion(299), all),
                Some((Region::Exon, Placement::Insertion(300)))
            );
        }
        // Only between the essential bases.
        let all = Placements::Insertions { lo: 201, hi: 201 };
        assert_eq!(plus.choose(Placement::Insertion(201), all), None);
    }

    #[test]
    fn rotated_seq() {
        assert_eq!(rotated("CGG", 0), "CGG");
        assert_eq!(rotated("CGG", 1), "GCG");
        assert_eq!(rotated("CGG", 2), "GGC");
        assert_eq!(rotated("CGG", 4), "GCG");
        assert_eq!(rotated("CGG", -1), "GGC");
        assert_eq!(rotated("", 3), "");
    }
}
