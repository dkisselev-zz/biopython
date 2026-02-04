/// 2-bit packing for DNA k-mers.
///
/// Encoding scheme:
/// - A = 00 (0)
/// - C = 01 (1)
/// - G = 10 (2)
/// - T = 11 (3)
///
/// This allows packing up to 32 nucleotides in a single u64.

/// Convert a nucleotide character to its 2-bit representation.
///
/// Returns None if the character is not A, C, G, or T.
#[inline]
pub fn encode_nucleotide(nuc: u8) -> Option<u8> {
    match nuc {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None, // N or other ambiguous base
    }
}

/// Convert a 2-bit representation back to a nucleotide character.
#[inline]
fn decode_nucleotide(bits: u8) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!(),
    }
}

/// Pack a DNA k-mer into a u64 using 2-bit encoding.
///
/// Returns None if:
/// - The k-mer contains non-ACGT characters (N, ambiguous bases)
/// - The k-mer length exceeds 32 (can't fit in u64)
///
/// # Examples
/// ```
/// use biopython_kmer_counter::packing::pack_kmer;
/// assert_eq!(pack_kmer(b"A"), Some(0b00));
/// assert_eq!(pack_kmer(b"C"), Some(0b01));
/// assert_eq!(pack_kmer(b"AC"), Some(0b0001));
/// assert_eq!(pack_kmer(b"ACGT"), Some(0b00011011));
/// assert_eq!(pack_kmer(b"ACGN"), None); // Contains N
/// ```
#[inline]
pub fn pack_kmer(kmer: &[u8]) -> Option<u64> {
    if kmer.len() > 32 {
        return None;
    }

    let mut packed: u64 = 0;

    for &nuc in kmer {
        let bits = encode_nucleotide(nuc)?;
        packed = (packed << 2) | (bits as u64);
    }

    Some(packed)
}

/// Unpack a u64 back to a DNA k-mer string.
///
/// # Arguments
/// * `packed` - The packed k-mer as u64
/// * `k` - The length of the k-mer (must be in range 1..=32)
///
/// # Panics
/// In debug builds, panics if k is 0 or > 32 (exceeds u64 capacity).
/// In release builds with invalid k, produces undefined behavior.
///
/// # Examples
/// ```
/// use biopython_kmer_counter::packing::unpack_kmer;
/// assert_eq!(unpack_kmer(0b00, 1), "A");
/// assert_eq!(unpack_kmer(0b01, 1), "C");
/// assert_eq!(unpack_kmer(0b0001, 2), "AC");
/// assert_eq!(unpack_kmer(0b00011011, 4), "ACGT");
/// ```
#[inline]
pub fn unpack_kmer(mut packed: u64, k: usize) -> String {
    // Validate k is within valid range for 2-bit packing in u64
    // u64 has 64 bits, each nucleotide uses 2 bits, so max k = 32
    debug_assert!(k > 0 && k <= 32, "k must be in range 1..=32, got {}", k);

    let mut result = Vec::with_capacity(k);

    for _ in 0..k {
        let nuc = decode_nucleotide((packed & 0b11) as u8);
        result.push(nuc);
        packed >>= 2;
    }

    result.reverse();
    String::from_utf8(result).expect("Invalid UTF-8 in unpacked k-mer")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_nucleotide() {
        assert_eq!(encode_nucleotide(b'A'), Some(0));
        assert_eq!(encode_nucleotide(b'C'), Some(1));
        assert_eq!(encode_nucleotide(b'G'), Some(2));
        assert_eq!(encode_nucleotide(b'T'), Some(3));
        assert_eq!(encode_nucleotide(b'N'), None);
        assert_eq!(encode_nucleotide(b'X'), None);
    }

    #[test]
    fn test_case_insensitive() {
        assert_eq!(encode_nucleotide(b'a'), Some(0));
        assert_eq!(encode_nucleotide(b'c'), Some(1));
        assert_eq!(encode_nucleotide(b'g'), Some(2));
        assert_eq!(encode_nucleotide(b't'), Some(3));
    }

    #[test]
    fn test_decode_nucleotide() {
        assert_eq!(decode_nucleotide(0), b'A');
        assert_eq!(decode_nucleotide(1), b'C');
        assert_eq!(decode_nucleotide(2), b'G');
        assert_eq!(decode_nucleotide(3), b'T');
    }

    #[test]
    fn test_pack_single_nucleotide() {
        assert_eq!(pack_kmer(b"A"), Some(0b00));
        assert_eq!(pack_kmer(b"C"), Some(0b01));
        assert_eq!(pack_kmer(b"G"), Some(0b10));
        assert_eq!(pack_kmer(b"T"), Some(0b11));
    }

    #[test]
    fn test_pack_multiple_nucleotides() {
        assert_eq!(pack_kmer(b"AC"), Some(0b0001));
        assert_eq!(pack_kmer(b"GT"), Some(0b1011));
        assert_eq!(pack_kmer(b"ACGT"), Some(0b00011011));
        assert_eq!(pack_kmer(b"AAAA"), Some(0b00000000));
        assert_eq!(pack_kmer(b"TTTT"), Some(0b11111111));
    }

    #[test]
    fn test_pack_with_ambiguous() {
        assert_eq!(pack_kmer(b"ACGN"), None);
        assert_eq!(pack_kmer(b"NACGT"), None);
        assert_eq!(pack_kmer(b"ACXGT"), None);
    }

    #[test]
    fn test_pack_too_long() {
        let long_kmer = b"A".repeat(33);
        assert_eq!(pack_kmer(&long_kmer), None);
    }

    #[test]
    fn test_unpack_single_nucleotide() {
        assert_eq!(unpack_kmer(0b00, 1), "A");
        assert_eq!(unpack_kmer(0b01, 1), "C");
        assert_eq!(unpack_kmer(0b10, 1), "G");
        assert_eq!(unpack_kmer(0b11, 1), "T");
    }

    #[test]
    fn test_unpack_multiple_nucleotides() {
        assert_eq!(unpack_kmer(0b0001, 2), "AC");
        assert_eq!(unpack_kmer(0b1011, 2), "GT");
        assert_eq!(unpack_kmer(0b00011011, 4), "ACGT");
    }

    #[test]
    fn test_pack_unpack_roundtrip() {
        let test_cases = vec![
            "A", "C", "G", "T",
            "AC", "GT", "ACGT",
            "ACGTACGT",
            "AAAAAAAAAA",
            "TTTTTTTTTT",
            "ACGTACGTACGTACGT",
        ];

        for kmer in test_cases {
            let k = kmer.len();
            let packed = pack_kmer(kmer.as_bytes()).expect("Pack failed");
            let unpacked = unpack_kmer(packed, k);
            assert_eq!(unpacked, kmer, "Roundtrip failed for {}", kmer);
        }
    }

    #[test]
    fn test_pack_max_length() {
        // Test k=32 (maximum that fits in u64)
        let kmer_32 = "A".repeat(32);
        assert!(pack_kmer(kmer_32.as_bytes()).is_some());

        // Verify unpack
        let packed = pack_kmer(kmer_32.as_bytes()).unwrap();
        let unpacked = unpack_kmer(packed, 32);
        assert_eq!(unpacked, kmer_32);
    }
}
