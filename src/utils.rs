pub(crate) fn trim_common_prefixes(reference: &str, alternative: &str) -> (usize, String, String) {
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.to_string(), alternative.to_string());
    }

    let mut trim = 0;
    while trim < reference.len() && trim < alternative.len() {
        if reference.chars().nth(trim) != alternative.chars().nth(trim) {
            break;
        }

        trim += 1;
    }

    (
        trim,
        reference[trim..].to_string(),
        alternative[trim..].to_string(),
    )
}

pub(crate) fn trim_common_suffixes(reference: &str, alternative: &str) -> (usize, String, String) {
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.to_string(), alternative.to_string());
    }

    let mut trim = 0;
    let mut i_r = reference.len();
    let mut i_a = alternative.len();
    let mut pad = 0;
    while trim < reference.len() && trim < alternative.len() {
        trim += 1;
        assert!(i_r > 0);
        assert!(i_a > 0);
        i_r -= 1;
        i_a -= 1;

        if reference.chars().nth(i_r) != alternative.chars().nth(i_a) {
            pad = 1;
            break;
        }
    }

    (
        trim - pad,
        reference[..(i_r + pad)].to_string(),
        alternative[..(i_a + pad)].to_string(),
    )
}

/// Reverse complementing shortcut.
pub(crate) fn revcomp(seq: &str) -> String {
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes()))
        .unwrap()
        .to_string()
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::utils::{revcomp, trim_common_prefixes, trim_common_suffixes};

    #[test]
    fn suffix_trimming() {
        assert_eq!(
            trim_common_suffixes("", ""),
            (0, "".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_suffixes("", "C"),
            (0, "".to_string(), "C".to_string())
        );
        assert_eq!(
            trim_common_suffixes("C", ""),
            (0, "C".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_suffixes("A", "AA"),
            (1, "".to_string(), "A".to_string())
        );
        assert_eq!(
            trim_common_suffixes("AT", "AG"),
            (0, "AT".to_string(), "AG".to_string())
        );
        assert_eq!(
            trim_common_suffixes("ATCG", "AGCG"),
            (2, "AT".to_string(), "AG".to_string())
        );
    }

    #[test]
    fn prefix_trimming() {
        assert_eq!(
            trim_common_prefixes("", ""),
            (0, "".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_prefixes("", "C"),
            (0, "".to_string(), "C".to_string())
        );
        assert_eq!(
            trim_common_prefixes("C", ""),
            (0, "C".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_prefixes("TA", "GA"),
            (0, "TA".to_string(), "GA".to_string())
        );
        assert_eq!(
            trim_common_prefixes("CGTA", "CGGA"),
            (2, "TA".to_string(), "GA".to_string())
        );
    }

    #[test]
    fn revcomp_cases() {
        assert_eq!(revcomp(""), "");
        assert_eq!(revcomp("A"), "T");
        assert_eq!(revcomp("AG"), "CT");
        assert_eq!(revcomp("CGAG"), "CTCG");
    }
}
