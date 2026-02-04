use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// Encode a single ASCII byte to a 2-bit value.
/// A=0, C=1, G=2, T=3.  Returns None for any other byte.
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Decode a 2-bit-packed k-mer back to an ASCII string.
fn decode_kmer(kmer: u64, k: u32) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut buf = vec![0u8; k as usize];
    let mut v = kmer;
    for i in (0..k as usize).rev() {
        buf[i] = bases[(v & 0x3) as usize];
        v >>= 2;
    }
    // buf contains only bytes from `bases` (A/C/G/T) — always valid UTF-8.
    String::from_utf8(buf).unwrap()
}

/// Accumulate k-mers from one sequence directly into an existing map.
/// Both the forward k-mer and its reverse complement are maintained as
/// rolling 2-bit-packed integers; no per-position O(k) recomputation.
/// Any non-ACGT byte resets both windows; k-mers that would span such a
/// byte are never emitted.  No allocation occurs inside this function.
fn count_into(counts: &mut FxHashMap<u64, u64>, seq: &[u8], k: usize, use_canonical: bool) {
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let rc_shift: u32 = 2 * (k as u32 - 1); // MSB bit-pair position in the RC word
    let mut kmer: u64 = 0;
    let mut rc: u64 = 0;
    let mut valid_run: usize = 0;

    for &byte in seq {
        if let Some(b) = encode_base(byte) {
            // Forward: new base enters at the LSB.
            kmer = ((kmer << 2) | b) & mask;
            // Reverse complement: complement (b ^ 0x3) enters at the MSB;
            // the previous LSB pair drops off the right end naturally.
            rc = (rc >> 2) | ((b ^ 0x3) << rc_shift);
            valid_run += 1;
            if valid_run >= k {
                let key = if use_canonical { kmer.min(rc) } else { kmer };
                *counts.entry(key).or_insert(0) += 1;
            }
        } else {
            valid_run = 0;
            kmer = 0;
            rc = 0;
        }
    }
}

/// Merge two FxHashMaps; values for shared keys are summed.
fn merge(mut a: FxHashMap<u64, u64>, b: FxHashMap<u64, u64>) -> FxHashMap<u64, u64> {
    for (key, val) in b {
        *a.entry(key).or_insert(0) += val;
    }
    a
}

/// Count k-mers across a batch of sequences in parallel using Rayon.
/// `fold` allocates exactly one FxHashMap per thread; every sequence assigned
/// to that thread is accumulated directly into it via `count_into`.
/// `reduce` merges the per-thread maps with a single pass at the end.
/// The caller owns the thread pool and reuses it across chunks.
fn count_batch(pool: &rayon::ThreadPool, seqs: &[&[u8]], k: usize, use_canonical: bool) -> FxHashMap<u64, u64> {
    pool.install(|| {
        seqs.par_iter()
            .fold(|| FxHashMap::default(), |mut acc, seq| {
                count_into(&mut acc, seq, k, use_canonical);
                acc
            })
            .reduce(|| FxHashMap::default(), merge)
    })
}

/// Number of sequences copied into Rust memory per chunk.
/// 4096 × 150 bp ≈ 600 KB — comfortably within L2 cache on most CPUs.
const CHUNK_SIZE: usize = 4096;

/// Python-facing entry point.  The input list is borrowed directly from
/// Python (no upfront copy).  Sequences are extracted in chunks of
/// ``CHUNK_SIZE``: each chunk is copied into owned Rust memory via a single
/// ``memcpy`` per sequence (through the ``PyBytes`` buffer protocol), the GIL
/// is released, and Rayon counts k-mers in parallel across the chunk.  Peak
/// Rust-side heap usage is therefore bounded to one chunk regardless of
/// total input length.  The thread pool is created once and reused for every
/// chunk.
///
/// Arguments
/// ---------
/// sequences : list[bytes]
///     Each element is a byte-string of one nucleotide sequence.
/// k : int
///     K-mer length.  Must be 1 ≤ k ≤ 32.
/// canonical : bool
///     If true, each k-mer is replaced by min(kmer, reverse_complement(kmer)).
/// threads : int
///     Number of Rayon worker threads.
///
/// Returns
/// -------
/// dict[str, int]
///     Mapping from ACGT k-mer string to occurrence count.
#[pyfunction]
fn count_kmers<'py>(
    py: Python<'py>,
    sequences: &Bound<'py, PyList>,
    k: usize,
    canonical: bool,
    threads: usize,
) -> PyResult<Bound<'py, PyDict>> {
    if k == 0 || k > 32 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be between 1 and 32 inclusive",
        ));
    }

    // Build the thread pool once; reused across every chunk.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build Rayon thread pool");

    let len = sequences.len();
    let mut counts: FxHashMap<u64, u64> = FxHashMap::default();

    for chunk_start in (0..len).step_by(CHUNK_SIZE) {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(len);

        // Copy this chunk into owned Rust memory (GIL held).
        // extract::<&[u8]> reads zero-copy from the PyBytes internal buffer;
        // to_vec() is the single memcpy that produces the owned data the
        // allow_threads boundary requires.
        let chunk: Vec<Vec<u8>> = (chunk_start..chunk_end)
            .map(|i| {
                let item = sequences.get_item(i)?;
                let slice: &[u8] = item.extract()?;
                Ok(slice.to_vec())
            })
            .collect::<PyResult<_>>()?;

        // Release the GIL; count k-mers across the chunk in parallel.
        let chunk_counts = py.allow_threads(|| {
            let refs: Vec<&[u8]> = chunk.iter().map(|v| v.as_slice()).collect();
            count_batch(&pool, &refs, k, canonical)
        });

        // Merge this chunk's counts into the running total.
        counts = merge(counts, chunk_counts);
    }

    // Build the Python dict from the merged counts.
    let dict = PyDict::new_bound(py);
    let k_u32 = k as u32;
    for (kmer_u64, count) in &counts {
        let kmer_str = decode_kmer(*kmer_u64, k_u32);
        dict.set_item(&kmer_str, *count)?;
    }
    Ok(dict)
}

#[pymodule]
fn _kmer_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(pyo3::wrap_pyfunction!(count_kmers, m)?)?;
    Ok(())
}
