/// Parallel k-mer counting using Rayon for multi-threading.

use ahash::AHashMap;
use rayon::prelude::*;

use crate::canonical::canonical_kmer;
use crate::packing::{encode_nucleotide, pack_kmer};

/// Count k-mers in a single sequence using sliding window bit-shift approach.
///
/// This implementation maintains a running packed k-mer and shifts it left by 2 bits
/// for each new nucleotide, achieving O(n) complexity instead of O(n*k).
///
/// # Arguments
/// * `sequence` - The DNA sequence as bytes
/// * `k` - The k-mer length
/// * `canonical` - Whether to use canonical k-mers
///
/// # Returns
/// HashMap of packed k-mers to their counts
///
/// # Performance
/// - O(n) complexity vs O(n*k) for pack_kmer approach
/// - Maintains sliding window buffer for constant-time updates
/// - Resets buffer on invalid characters (N or ambiguous bases)
fn count_kmers_in_sequence(
    sequence: &[u8],
    k: usize,
    canonical: bool,
) -> AHashMap<u64, usize> {
    let mut counts = AHashMap::new();

    if sequence.len() < k {
        return counts;
    }

    // Validate k is in valid range for 2-bit packing
    debug_assert!(k > 0 && k <= 32, "k must be in range 1..=32, got {}", k);

    // Mask to keep only k*2 bits (k nucleotides * 2 bits each)
    let mask = (1u64 << (k * 2)) - 1;

    // Sliding window buffer and valid nucleotide counter
    let mut current_packed = 0u64;
    let mut valid_count = 0usize; // Number of consecutive valid nucleotides seen

    for &nuc in sequence.iter() {
        // Try to encode the nucleotide
        if let Some(bits) = encode_nucleotide(nuc) {
            // Shift left by 2 bits, add new nucleotide, and mask to k*2 bits
            current_packed = ((current_packed << 2) | (bits as u64)) & mask;
            valid_count += 1;

            // Once we have k valid nucleotides, we have a complete k-mer
            if valid_count >= k {
                let mut packed = current_packed;

                // Apply canonical transformation if requested
                if canonical {
                    packed = canonical_kmer(packed, k);
                }

                *counts.entry(packed).or_insert(0) += 1;
            }
        } else {
            // Invalid character (N or ambiguous base) - reset sliding window
            current_packed = 0;
            valid_count = 0;
        }
    }

    counts
}

/// Generic function that counts k-mers from any parallel iterator of sequences.
///
/// This generic implementation allows streaming sequences without loading everything into RAM.
/// It accepts any parallel iterator that yields byte slices, enabling zero-copy processing.
///
/// # Design Rationale
/// By accepting `impl ParallelIterator<Item = &[u8]>`, this function can process:
/// - Slices from a Vec (current use case)
/// - Memory-mapped files
/// - Streaming decompressed data
/// - Any other parallel iterator of byte sequences
///
/// This enables future Python API extensions for true streaming from Bio.SeqIO
/// without loading entire FASTA/FASTQ files into RAM.
///
/// # Arguments
/// * `sequences` - Any parallel iterator yielding `&[u8]` byte slices
/// * `k` - The k-mer length
/// * `canonical` - Whether to use canonical k-mers
///
/// # Returns
/// HashMap of packed k-mers to their total counts
///
/// # Performance
/// - Works with streaming data without full materialization
/// - Zero-copy: processes byte slices directly without allocation
/// - Thread-local accumulators reduce allocation pressure
/// - Parallel merge eliminates sequential bottleneck
pub fn count_kmers_parallel_generic<'a, I>(
    sequences: I,
    k: usize,
    canonical: bool,
) -> AHashMap<u64, usize>
where
    I: ParallelIterator<Item = &'a [u8]>,
{
    sequences
        // Fold: Each thread accumulates counts into a thread-local HashMap
        .fold(
            || AHashMap::new(), // Initialize thread-local HashMap
            |mut acc, seq| {
                // Count k-mers in this sequence
                let seq_counts = count_kmers_in_sequence(seq, k, canonical);

                // Merge into thread-local accumulator
                for (kmer, count) in seq_counts {
                    *acc.entry(kmer).or_insert(0) += count;
                }

                acc
            },
        )
        // Reduce: Merge thread-local HashMaps in parallel
        .reduce(
            || AHashMap::new(), // Identity for reduction
            |mut map1, map2| {
                // Merge map2 into map1
                for (kmer, count) in map2 {
                    *map1.entry(kmer).or_insert(0) += count;
                }
                map1
            },
        )
}

/// Count k-mers across multiple sequences in parallel using fold/reduce pattern.
///
/// This is the public API that accepts a slice of Vec<u8>. Internally delegates to
/// the generic streaming implementation for maximum flexibility.
///
/// # Arguments
/// * `sequences` - Slice of DNA sequences as Vec<u8>
/// * `k` - The k-mer length
/// * `canonical` - Whether to use canonical k-mers
///
/// # Returns
/// HashMap of packed k-mers to their total counts across all sequences
///
/// # Performance
/// - Parallelizes across CPU cores using Rayon
/// - Uses ahash for fast hashing of DNA k-mers
/// - Thread-local accumulators reduce allocation pressure
/// - Parallel merge eliminates sequential bottleneck
/// - Zero-copy: sequences are processed as slices without duplication
pub fn count_kmers_parallel(
    sequences: &[Vec<u8>],
    k: usize,
    canonical: bool,
) -> AHashMap<u64, usize> {
    // Convert Vec<Vec<u8>> to parallel iterator of &[u8] slices (zero-copy)
    count_kmers_parallel_generic(sequences.par_iter().map(|v| v.as_slice()), k, canonical)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::packing::unpack_kmer;

    #[test]
    fn test_count_kmers_in_sequence() {
        let seq = b"ACGTACGT";
        let counts = count_kmers_in_sequence(seq, 3, false);

        // Verify counts
        let acg_packed = pack_kmer(b"ACG").unwrap();
        let cgt_packed = pack_kmer(b"CGT").unwrap();
        let gta_packed = pack_kmer(b"GTA").unwrap();
        let tac_packed = pack_kmer(b"TAC").unwrap();

        assert_eq!(counts.get(&acg_packed), Some(&2));
        assert_eq!(counts.get(&cgt_packed), Some(&2));
        assert_eq!(counts.get(&gta_packed), Some(&2));
        assert_eq!(counts.get(&tac_packed), Some(&1));
    }

    #[test]
    fn test_count_kmers_empty_sequence() {
        let seq = b"";
        let counts = count_kmers_in_sequence(seq, 3, false);
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_count_kmers_shorter_than_k() {
        let seq = b"AC";
        let counts = count_kmers_in_sequence(seq, 10, false);
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_count_kmers_with_ambiguous() {
        let seq = b"ACNGT";
        let counts = count_kmers_in_sequence(seq, 3, false);

        // ACN, CNT, NGT all contain N and should be skipped
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_count_kmers_canonical() {
        let seq = b"ACGTACGT";
        let counts = count_kmers_in_sequence(seq, 3, true);

        // With canonical, ACG and CGT should be counted together
        // ACG and CGT are reverse complements: ACG <-> CGT
        // The canonical (lexicographically smaller) is ACG

        let acg_packed = pack_kmer(b"ACG").unwrap();
        let canonical_acg = canonical_kmer(acg_packed, 3);

        // ACG appears 2 times, CGT appears 2 times
        // With canonical, they should be counted together as 4
        assert_eq!(counts.get(&canonical_acg), Some(&4));
    }

    #[test]
    fn test_count_kmers_parallel_single_sequence() {
        let sequences = vec![b"ACGTACGT".to_vec()];
        let counts = count_kmers_parallel(&sequences, 3, false);

        let acg_packed = pack_kmer(b"ACG").unwrap();
        let cgt_packed = pack_kmer(b"CGT").unwrap();

        assert_eq!(counts.get(&acg_packed), Some(&2));
        assert_eq!(counts.get(&cgt_packed), Some(&2));
    }

    #[test]
    fn test_count_kmers_parallel_multiple_sequences() {
        let sequences = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"TGCA".to_vec(),
        ];
        let counts = count_kmers_parallel(&sequences, 3, false);

        let acg_packed = pack_kmer(b"ACG").unwrap();
        let cgt_packed = pack_kmer(b"CGT").unwrap();
        let tgc_packed = pack_kmer(b"TGC").unwrap();
        let gca_packed = pack_kmer(b"GCA").unwrap();

        assert_eq!(counts.get(&acg_packed), Some(&2)); // Appears in first two sequences
        assert_eq!(counts.get(&cgt_packed), Some(&2)); // Appears in first two sequences
        assert_eq!(counts.get(&tgc_packed), Some(&1)); // Appears in third sequence
        assert_eq!(counts.get(&gca_packed), Some(&1)); // Appears in third sequence
    }

    #[test]
    fn test_count_kmers_parallel_canonical() {
        let sequences = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
        ];
        let counts = count_kmers_parallel(&sequences, 3, true);

        // ACG and CGT are reverse complements
        let acg_packed = pack_kmer(b"ACG").unwrap();
        let canonical_acg = canonical_kmer(acg_packed, 3);

        // ACG appears 2 times, CGT appears 2 times across both sequences
        // Canonical should count them together: 4 total
        assert_eq!(counts.get(&canonical_acg), Some(&4));
    }


    #[test]
    fn test_large_sequence() {
        // Generate a longer sequence
        let sequence: Vec<u8> = b"ACGT".iter().cycle().take(1000).copied().collect();
        let counts = count_kmers_in_sequence(&sequence, 21, false);

        // Should have some k-mers counted
        assert!(counts.len() > 0);

        // Verify total count (should be sequence.len() - k + 1)
        let total: usize = counts.values().sum();
        assert_eq!(total, sequence.len() - 21 + 1);
    }

    #[test]
    fn test_verify_unpacking() {
        // Verify that we can unpack the counted k-mers back to strings
        let seq = b"ACGTACGT";
        let counts = count_kmers_in_sequence(seq, 3, false);

        for (packed, _count) in counts.iter() {
            let unpacked = unpack_kmer(*packed, 3);
            // Verify it's a valid 3-mer
            assert_eq!(unpacked.len(), 3);
            assert!(unpacked.chars().all(|c| "ACGT".contains(c)));
        }
    }
}
