/// Canonical k-mer support for treating a k-mer and its reverse complement as identical.
///
/// This is useful for double-stranded DNA analysis where strand orientation doesn't matter.

/// Compute the reverse complement of a packed k-mer.
///
/// The reverse complement operation:
/// 1. Reverses the order of nucleotides
/// 2. Complements each nucleotide: A<->T, C<->G
///
/// The complement mapping uses XOR with 0b11:
/// - A (00) XOR 11 = T (11)
/// - C (01) XOR 11 = G (10)
/// - G (10) XOR 11 = C (01)
/// - T (11) XOR 11 = A (00)
///
/// # Arguments
/// * `packed` - The packed k-mer as u64
/// * `k` - The length of the k-mer (must be in range 1..=32)
///
/// # Panics
/// In debug builds, panics if k is 0 or > 32 (exceeds u64 capacity).
/// In release builds with invalid k, produces undefined behavior.
///
/// # Performance
/// This iterative implementation runs in O(k) time. For k <= 32, this is typically
/// fast enough (~32 loop iterations max) and provides correct behavior for all k values.
///
/// # Examples
/// ```
/// use biopython_kmer_counter::packing::{pack_kmer, unpack_kmer};
/// use biopython_kmer_counter::canonical::reverse_complement_packed;
///
/// // ACG -> CGT (reverse complement)
/// let acg = pack_kmer(b"ACG").unwrap();
/// let rc = reverse_complement_packed(acg, 3);
/// assert_eq!(unpack_kmer(rc, 3), "CGT");
/// ```
#[inline]
pub fn reverse_complement_packed(packed: u64, k: usize) -> u64 {
    // Validate k is within valid range for 2-bit packing in u64
    // u64 has 64 bits, each nucleotide uses 2 bits, so max k = 32
    debug_assert!(k > 0 && k <= 32, "k must be in range 1..=32, got {}", k);

    // Use simple iterative approach for correctness
    // Each nucleotide is 2 bits, complement mapping: A(00)<->T(11), C(01)<->G(10)
    let mut result = 0u64;

    for i in 0..k {
        // Extract nucleotide at position i from the right (LSB side)
        let nuc = (packed >> (i * 2)) & 0b11;

        // Complement: XOR with 0b11 (flips both bits)
        // A(00) -> T(11), C(01) -> G(10), G(10) -> C(01), T(11) -> A(00)
        let comp_nuc = nuc ^ 0b11;

        // Place at position (k-1-i) in result (reverse order)
        result |= comp_nuc << ((k - 1 - i) * 2);
    }

    result
}

/// Get the canonical representation of a k-mer.
///
/// The canonical k-mer is the lexicographically smaller of the k-mer and its reverse complement.
/// This ensures that a k-mer and its reverse complement are treated as the same entity.
///
/// # Arguments
/// * `packed` - The packed k-mer as u64
/// * `k` - The length of the k-mer (must be in range 1..=32)
///
/// # Returns
/// The packed representation of the canonical k-mer (min of forward and reverse complement).
///
/// # Panics
/// In debug builds, panics if k is 0 or > 32 (exceeds u64 capacity).
/// In release builds with invalid k, produces undefined behavior.
///
/// # Examples
/// ```
/// use biopython_kmer_counter::packing::pack_kmer;
/// use biopython_kmer_counter::canonical::canonical_kmer;
///
/// let acgt = pack_kmer(b"ACGT").unwrap();
/// let acgt_rc = pack_kmer(b"ACGT").unwrap(); // ACGT is its own RC (palindrome)
/// assert_eq!(canonical_kmer(acgt, 4), canonical_kmer(acgt_rc, 4));
/// ```
#[inline]
pub fn canonical_kmer(packed: u64, k: usize) -> u64 {
    debug_assert!(k > 0 && k <= 32, "k must be in range 1..=32, got {}", k);
    let rc = reverse_complement_packed(packed, k);
    packed.min(rc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packing::{pack_kmer, unpack_kmer};

    #[test]
    fn test_reverse_complement_single() {
        // A -> T
        let a = pack_kmer(b"A").unwrap();
        let t = pack_kmer(b"T").unwrap();
        assert_eq!(reverse_complement_packed(a, 1), t);

        // C -> G
        let c = pack_kmer(b"C").unwrap();
        let g = pack_kmer(b"G").unwrap();
        assert_eq!(reverse_complement_packed(c, 1), g);

        // G -> C
        assert_eq!(reverse_complement_packed(g, 1), c);

        // T -> A
        assert_eq!(reverse_complement_packed(t, 1), a);
    }

    #[test]
    fn test_reverse_complement_multiple() {
        // AC -> GT
        let ac = pack_kmer(b"AC").unwrap();
        let gt = pack_kmer(b"GT").unwrap();
        assert_eq!(reverse_complement_packed(ac, 2), gt);

        // ACGT -> ACGT (palindrome)
        let acgt = pack_kmer(b"ACGT").unwrap();
        assert_eq!(reverse_complement_packed(acgt, 4), acgt);

        // AAAA -> TTTT
        let aaaa = pack_kmer(b"AAAA").unwrap();
        let tttt = pack_kmer(b"TTTT").unwrap();
        assert_eq!(reverse_complement_packed(aaaa, 4), tttt);

        // ACGTACGT -> ACGTACGT (palindrome)
        let acgtacgt = pack_kmer(b"ACGTACGT").unwrap();
        assert_eq!(reverse_complement_packed(acgtacgt, 8), acgtacgt);
    }

    #[test]
    fn test_reverse_complement_verify_string() {
        let test_cases = vec![
            ("A", "T"),
            ("C", "G"),
            ("G", "C"),
            ("T", "A"),
            ("AC", "GT"),
            ("GT", "AC"),
            ("ACG", "CGT"),
            ("ACGT", "ACGT"), // palindrome
            ("AAAA", "TTTT"),
            ("CGCG", "CGCG"), // palindrome
        ];

        for (forward, expected_rc) in test_cases {
            let packed = pack_kmer(forward.as_bytes()).unwrap();
            let rc_packed = reverse_complement_packed(packed, forward.len());
            let rc_string = unpack_kmer(rc_packed, forward.len());
            assert_eq!(
                rc_string, expected_rc,
                "RC of {} should be {}, got {}",
                forward, expected_rc, rc_string
            );
        }
    }

    #[test]
    fn test_reverse_complement_double() {
        // RC of RC should give original
        let test_kmers = vec!["A", "C", "AC", "GT", "ACGT", "ACGTACGT"];

        for kmer in test_kmers {
            let packed = pack_kmer(kmer.as_bytes()).unwrap();
            let k = kmer.len();
            let rc1 = reverse_complement_packed(packed, k);
            let rc2 = reverse_complement_packed(rc1, k);
            assert_eq!(
                packed, rc2,
                "RC of RC of {} should equal original",
                kmer
            );
        }
    }

    #[test]
    fn test_canonical_kmer() {
        // For non-palindromes, canonical should be the smaller one
        let ac = pack_kmer(b"AC").unwrap();
        let gt = pack_kmer(b"GT").unwrap();
        assert_eq!(canonical_kmer(ac, 2), canonical_kmer(gt, 2));
        assert_eq!(canonical_kmer(ac, 2), ac); // AC < GT lexicographically

        // For palindromes, canonical should equal original
        let acgt = pack_kmer(b"ACGT").unwrap();
        assert_eq!(canonical_kmer(acgt, 4), acgt);
    }

    #[test]
    fn test_canonical_symmetry() {
        // A k-mer and its RC should have the same canonical form
        let test_kmers = vec!["AC", "GT", "ACG", "CGT", "AAAC", "GTTT"];

        for kmer in test_kmers {
            let packed = pack_kmer(kmer.as_bytes()).unwrap();
            let k = kmer.len();
            let rc_packed = reverse_complement_packed(packed, k);

            let canonical1 = canonical_kmer(packed, k);
            let canonical2 = canonical_kmer(rc_packed, k);

            assert_eq!(
                canonical1, canonical2,
                "Canonical form of {} and its RC should be equal",
                kmer
            );
        }
    }

    #[test]
    fn test_canonical_ordering() {
        // Test that canonical returns the lexicographically smaller k-mer
        let test_cases = vec![
            ("AAA", "TTT", "AAA"), // AAA < TTT
            ("AAC", "GTT", "AAC"), // AAC < GTT
            ("ACG", "CGT", "ACG"), // ACG < CGT
            ("ACGT", "ACGT", "ACGT"), // Palindrome
        ];

        for (kmer1, kmer2, expected_canonical) in test_cases {
            let packed1 = pack_kmer(kmer1.as_bytes()).unwrap();
            let packed2 = pack_kmer(kmer2.as_bytes()).unwrap();
            let k = kmer1.len();

            let canonical1 = canonical_kmer(packed1, k);
            let canonical2 = canonical_kmer(packed2, k);

            let canonical_str = unpack_kmer(canonical1, k);
            assert_eq!(
                canonical_str, expected_canonical,
                "Canonical of {} and {} should be {}",
                kmer1, kmer2, expected_canonical
            );
            assert_eq!(canonical1, canonical2, "Canonical forms should match");
        }
    }
}
