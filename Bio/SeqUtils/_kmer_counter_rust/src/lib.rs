/// Biopython k-mer counter Rust extension.
///
/// This module provides high-performance k-mer counting using:
/// - 2-bit packing for memory efficiency
/// - Parallel processing with Rayon
/// - Fast hashing with ahash
/// - Stateful accumulation for batch processing
/// - Zero-copy byte slice processing (no .to_uppercase() allocation)
/// - Case-insensitive encoding (handles both 'A'/'a', 'C'/'c', etc.)
///
/// Exposed to Python via PyO3 bindings.

use ahash::AHashMap;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

mod canonical;
mod counting;
mod packing;

use counting::count_kmers_parallel_generic;
use packing::unpack_kmer;

/// Stateful k-mer counter that can accumulate counts over multiple batches.
///
/// This class allows incremental counting, which is useful for:
/// - Streaming large FASTA/FASTQ files
/// - Processing data in chunks to control memory usage
/// - Accumulating counts from multiple sources
///
/// # Examples
/// ```python
/// from Bio.SeqUtils._kmer_counter_rust import KmerCounter
///
/// # Create counter
/// counter = KmerCounter(k=3, canonical=False)
///
/// # Add sequences in batches
/// counter.update(["ACGTACGT"])
/// counter.update(["TGCATGCA"])
///
/// # Get final counts
/// counts = counter.get_counts()
/// ```
#[pyclass]
struct KmerCounter {
    /// K-mer length
    k: usize,
    /// Whether to use canonical k-mers (treat k-mer and RC as identical)
    canonical: bool,
    /// Internal state: packed k-mers to counts
    internal_map: AHashMap<u64, usize>,
}

#[pymethods]
impl KmerCounter {
    /// Create a new k-mer counter.
    ///
    /// # Arguments
    /// * `k` - K-mer length (must be > 0 and < 32)
    /// * `canonical` - Whether to use canonical k-mers
    ///
    /// # Raises
    /// * `ValueError` - If k is 0 or >= 32
    ///
    /// # Examples
    /// ```python
    /// counter = KmerCounter(k=21, canonical=True)
    /// ```
    #[new]
    fn new(k: usize, canonical: bool) -> PyResult<Self> {
        // Validate k
        if k == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "k must be positive",
            ));
        }

        if k > 31 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "k must be < 32 for Rust implementation (use Python fallback for larger k)",
            ));
        }

        Ok(KmerCounter {
            k,
            canonical,
            internal_map: AHashMap::new(),
        })
    }

    /// Update the counter with a batch of sequences.
    ///
    /// This method can be called multiple times to accumulate counts
    /// from different batches of sequences. Counts are merged into the
    /// internal state.
    ///
    /// # Arguments
    /// * `sequences` - List of DNA sequences as strings (case-insensitive)
    ///
    /// # Examples
    /// ```python
    /// counter = KmerCounter(k=3, canonical=False)
    /// counter.update(["ACGTACGT"])
    /// counter.update(["TGCATGCA"])  # Accumulates with previous batch
    /// ```
    fn update(&mut self, sequences: Vec<String>) -> PyResult<()> {
        // Convert strings to byte slices (zero-copy view of existing strings)
        // Note: encode_nucleotide in packing module handles both uppercase and lowercase
        // This avoids allocating new uppercased strings
        let seq_bytes: Vec<&[u8]> = sequences
            .iter()
            .map(|s| s.as_bytes())
            .collect();

        // Count k-mers in this batch using generic iterator (zero-copy)
        let batch_counts = count_kmers_parallel_generic(
            seq_bytes.par_iter().copied(),
            self.k,
            self.canonical,
        );

        // Merge batch counts into internal state
        for (packed_kmer, count) in batch_counts {
            *self.internal_map.entry(packed_kmer).or_insert(0) += count;
        }

        Ok(())
    }

    /// Get the accumulated k-mer counts.
    ///
    /// Returns a dictionary mapping k-mer strings to their total counts
    /// across all batches that have been processed with `update()`.
    ///
    /// # Returns
    /// Dictionary of k-mer strings to counts
    ///
    /// # Examples
    /// ```python
    /// counter = KmerCounter(k=3, canonical=False)
    /// counter.update(["ACGTACGT"])
    /// counts = counter.get_counts()
    /// print(counts["ACG"])  # 2
    /// ```
    fn get_counts(&self) -> PyResult<HashMap<String, usize>> {
        let mut result = HashMap::new();

        for (&packed_kmer, &count) in self.internal_map.iter() {
            let kmer_string = unpack_kmer(packed_kmer, self.k);
            result.insert(kmer_string, count);
        }

        Ok(result)
    }

    /// Clear all accumulated counts and reset to empty state.
    ///
    /// Useful for reusing the same counter object with different data.
    ///
    /// # Examples
    /// ```python
    /// counter = KmerCounter(k=3, canonical=False)
    /// counter.update(["ACGTACGT"])
    /// counter.clear()
    /// counts = counter.get_counts()  # Empty dict
    /// ```
    fn clear(&mut self) -> PyResult<()> {
        self.internal_map.clear();
        Ok(())
    }

    /// Return the number of unique k-mers counted so far.
    ///
    /// # Examples
    /// ```python
    /// counter = KmerCounter(k=3, canonical=False)
    /// counter.update(["ACGTACGT"])
    /// print(len(counter))  # Number of unique k-mers
    /// ```
    fn __len__(&self) -> PyResult<usize> {
        Ok(self.internal_map.len())
    }

    /// String representation of the counter.
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "KmerCounter(k={}, canonical={}, unique_kmers={})",
            self.k,
            self.canonical,
            self.internal_map.len()
        ))
    }
}

/// Count k-mers in DNA sequences (convenience function for backward compatibility).
///
/// This is a convenience function that creates a KmerCounter, updates it once,
/// and returns the counts. For batch processing, use the KmerCounter class directly.
///
/// # Arguments
/// * `sequences` - List of DNA sequences as strings (case-insensitive)
/// * `k` - K-mer length (must be < 32)
/// * `canonical` - Whether to use canonical k-mers (treat k-mer and RC as same)
///
/// # Returns
/// Dictionary mapping k-mer strings to their counts
///
/// # Examples
/// ```python
/// from Bio.SeqUtils._kmer_counter_rust import count_kmers_rust
/// counts = count_kmers_rust(["ACGTACGT"], k=3, canonical=False)
/// assert counts["ACG"] == 2
/// ```
#[pyfunction]
fn count_kmers_rust(
    sequences: Vec<String>,
    k: usize,
    canonical: bool,
) -> PyResult<HashMap<String, usize>> {
    // Use the KmerCounter class for implementation
    let mut counter = KmerCounter::new(k, canonical)?;
    counter.update(sequences)?;
    counter.get_counts()
}

/// Python module definition.
#[pymodule]
fn _kmer_counter_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<KmerCounter>()?;
    m.add_function(wrap_pyfunction!(count_kmers_rust, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_counter_new() {
        let counter = KmerCounter::new(3, false).unwrap();
        assert_eq!(counter.k, 3);
        assert_eq!(counter.canonical, false);
        assert_eq!(counter.internal_map.len(), 0);
    }

    #[test]
    fn test_kmer_counter_invalid_k() {
        let result = KmerCounter::new(0, false);
        assert!(result.is_err());

        let result = KmerCounter::new(32, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_kmer_counter_update() {
        let mut counter = KmerCounter::new(3, false).unwrap();
        counter.update(vec!["ACGTACGT".to_string()]).unwrap();

        let counts = counter.get_counts().unwrap();
        assert_eq!(counts.get("ACG"), Some(&2));
        assert_eq!(counts.get("CGT"), Some(&2));
    }

    #[test]
    fn test_kmer_counter_accumulate() {
        let mut counter = KmerCounter::new(3, false).unwrap();

        // First batch
        counter.update(vec!["ACGT".to_string()]).unwrap();

        // Second batch - should accumulate
        counter.update(vec!["ACGT".to_string()]).unwrap();

        let counts = counter.get_counts().unwrap();
        assert_eq!(counts.get("ACG"), Some(&2)); // Count doubled
        assert_eq!(counts.get("CGT"), Some(&2)); // Count doubled
    }

    #[test]
    fn test_kmer_counter_clear() {
        let mut counter = KmerCounter::new(3, false).unwrap();
        counter.update(vec!["ACGTACGT".to_string()]).unwrap();

        assert_eq!(counter.__len__().unwrap(), 4); // 4 unique k-mers

        counter.clear().unwrap();
        assert_eq!(counter.__len__().unwrap(), 0);

        let counts = counter.get_counts().unwrap();
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_kmer_counter_canonical() {
        let mut counter = KmerCounter::new(3, true).unwrap();
        counter.update(vec!["ACGTACGT".to_string()]).unwrap();

        let counts = counter.get_counts().unwrap();
        // ACG and CGT are reverse complements, counted together
        assert_eq!(counts.get("ACG"), Some(&4));
    }

    // Backward compatibility tests for count_kmers_rust function
    #[test]
    fn test_count_kmers_rust_basic() {
        let sequences = vec!["ACGTACGT".to_string()];
        let counts = count_kmers_rust(sequences, 3, false).unwrap();

        assert_eq!(counts.get("ACG"), Some(&2));
        assert_eq!(counts.get("CGT"), Some(&2));
        assert_eq!(counts.get("GTA"), Some(&2));
        assert_eq!(counts.get("TAC"), Some(&1));
    }

    #[test]
    fn test_count_kmers_rust_canonical() {
        let sequences = vec!["ACGTACGT".to_string()];
        let counts = count_kmers_rust(sequences, 3, true).unwrap();

        // ACG and CGT are reverse complements
        // With canonical, they should be counted together
        // The canonical form will be ACG (lexicographically smaller)
        assert_eq!(counts.get("ACG"), Some(&4));
    }

    #[test]
    fn test_count_kmers_rust_multiple_sequences() {
        let sequences = vec!["ACGT".to_string(), "ACGT".to_string()];
        let counts = count_kmers_rust(sequences, 3, false).unwrap();

        assert_eq!(counts.get("ACG"), Some(&2));
        assert_eq!(counts.get("CGT"), Some(&2));
    }

    #[test]
    fn test_count_kmers_rust_case_insensitive() {
        let sequences_lower = vec!["acgt".to_string()];
        let sequences_upper = vec!["ACGT".to_string()];

        let counts_lower = count_kmers_rust(sequences_lower, 3, false).unwrap();
        let counts_upper = count_kmers_rust(sequences_upper, 3, false).unwrap();

        assert_eq!(counts_lower, counts_upper);
    }

    #[test]
    fn test_count_kmers_rust_with_ambiguous() {
        let sequences = vec!["ACNGT".to_string()];
        let counts = count_kmers_rust(sequences, 3, false).unwrap();

        // All k-mers contain N, so should be empty
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_count_kmers_rust_invalid_k() {
        let sequences = vec!["ACGT".to_string()];

        // k = 0 should error
        let result = count_kmers_rust(sequences.clone(), 0, false);
        assert!(result.is_err());

        // k > 31 should error
        let result = count_kmers_rust(sequences, 32, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_count_kmers_rust_empty_sequence() {
        let sequences = vec!["".to_string()];
        let counts = count_kmers_rust(sequences, 3, false).unwrap();
        assert_eq!(counts.len(), 0);
    }

    #[test]
    fn test_count_kmers_rust_sequence_shorter_than_k() {
        let sequences = vec!["AC".to_string()];
        let counts = count_kmers_rust(sequences, 10, false).unwrap();
        assert_eq!(counts.len(), 0);
    }
}
