//! # tcrdist_rs
//!
//! A library for computing the distances between immune cell receptors.
mod distance;
mod match_table;
use crate::distance as _distance;

#[cfg(all(feature = "pyo3"))]
use pyo3::prelude::*;
use triple_accel;

pub use match_table::{
    amino_distances, cdr1_distances, cdr2_distances, gene_distance, phmc_distances,
};

/// Compute the Hamming distance between two strings.
///
/// The strings being compared must have equal lengths or the program will terminate.
/// Moreover, the strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// s1 : str
///     A string.
/// s2 : str
///     Another string
///
/// Returns
/// -------
/// int
///     The Hamming distance between the two strings.
///
/// Examples
/// --------
/// >>> s1 = "abcdefg"
/// >>> s2 = "abddefg"
/// >>> assert(hamming(s1, s2) == 1)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn hamming(s1: &str, s2: &str) -> PyResult<u32> {
    Ok(triple_accel::hamming(&s1.as_bytes(), &s2.as_bytes()))
}

/// Compute the Hamming distance matrix on an iterable of strings.
///
/// This returns the upper right triangle of the distance matrix.
/// The strings being compared must have equal lengths or the program will terminate.
/// Moreover, the strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs : iterable of str
///     The strings to be compared. They must all be the same length and have
///     an appropriate representation in bytes.
///
/// Returns
/// -------
/// list of int
///     The Hamming distances among the strings.
///
/// Examples
/// --------
/// >>> seqs = ["abc", "abd", "fcd"]
/// >>> assert(hamming_matrix(seqs, parallel=False) == [1, 3, 2])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, parallel=false))]
fn hamming_matrix(seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_matrix(&seqs, parallel, "hamming"))
}

/// Compute the Hamming distance between one string and many others.
///
/// The strings being compared must have equal lengths or the program will terminate.
/// Moreover, the strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seq : str
///     The string against which all others will be compared.
/// seqs : iterable of str
///     The other strings being compared.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Hamming distances among the strings.
///
/// Examples
/// --------
/// >>> seq = "abb"
/// >>> seqs = ["abc", "abd", "fcd"]
/// >>> assert(hamming_one_to_many(seq, seqs, parallel=False) == [1, 1, 3])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, parallel=false))]
fn hamming_one_to_many(seq: &str, seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_one_to_many(
        &seq, &seqs, parallel, "hamming",
    ))
}

/// Compute the Hamming distance between many strings and many others.
///
/// The strings being compared must have equal lengths or the program will terminate.
/// Moreover, the strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs1 : iterable of str
///     The first iterable of strings.
/// seqs2 : iterable of str
///     The other iterable of strings.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Hamming distances among the strings.
///
/// Examples
/// --------
/// >>> seqs1 = ["abb", "abc"]
/// >>> seqs2 = ["abc", "abd", "fcd"]
/// >>> assert(hamming_many_to_many(seqs1, seqs2, parallel=False) == [1, 1, 3, 0, 1, 3])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn hamming_many_to_many(seqs1: Vec<&str>, seqs2: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_many_to_many(
        &seqs1, &seqs2, parallel, "hamming",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, threshold, parallel=false))]
fn hamming_neighbor_matrix(
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_matrix(
        &seqs, threshold, parallel, "hamming",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, threshold, parallel=false))]
fn hamming_neighbor_one_to_many(
    seq: &str,
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<usize>> {
    Ok(_distance::str_neighbor_one_to_many(
        seq, &seqs, threshold, parallel, "hamming",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, threshold, parallel=false))]
fn hamming_neighbor_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_many_to_many(
        &seqs1, &seqs2, threshold, parallel, "hamming",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn hamming_bin_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    parallel: bool,
) -> PyResult<Vec<u32>> {
    Ok(_distance::str_bin_many_to_many(
        &seqs1, &seqs2, parallel, "hamming",
    ))
}

/// Compute the Levenshtein distance between two strings.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// s1 : str
///     A string.
/// s2 : str
///     Another string
///
/// Returns
/// -------
/// int
///     The Levenshtein distance between the two strings.
///
/// Examples
/// --------
/// >>> s1 = "abcdefg"
/// >>> s2 = "abdcd defgggg"
/// >>> assert(levenshtein(s1, s2) == 6)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn levenshtein(s1: &str, s2: &str) -> PyResult<u32> {
    Ok(triple_accel::levenshtein(&s1.as_bytes(), &s2.as_bytes()))
}

/// Compute the Levenshtein distance matrix on an iterable of strings.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs : iterable str
///     The strings to be compared. The must have an appropriate representation
///     as byte strings.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seqs = ["Monday", "Tuesday", "Wednesday", "Thursday"]
/// >>> assert(levenshtein_matrix(seqs, parallel=False) == [4, 5, 5, 4, 2, 5])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, parallel=false))]
fn levenshtein_matrix(seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_matrix(&seqs, parallel, "levenshtein"))
}

/// Compute the Levenshtein distance between one string and many others.
///
/// This returns the upper right triangle of the distance matrix.
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seq : str
///     The string against which all others will be compared.
/// seqs : iterable of str
///     The other strings being compared.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seq = "Sunday"
/// >>> seqs = ["Monday", "Tuesday", "Wednesday", "Thursday"]
/// >>> assert(levenshtein_one_to_many(seq, seqs, parallel=False) == [2, 3, 5, 4])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, parallel=false))]
fn levenshtein_one_to_many(seq: &str, seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_one_to_many(
        &seq,
        &seqs,
        parallel,
        "levenshtein",
    ))
}

/// Compute the Levenshtein distance between many strings and many others.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs1 : iterable of str
///     The first iterable of strings.
/// seqs : iterable of str
///     The second iterable of strings.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seqs1 = ["Sunday", "Saturday"]
/// >>> seqs2 = ["Monday", "Tuesday", "Wednesday"]
/// >>> assert(levenshtein_many_to_many(seqs1, seqs2, parallel=False) == [2, 3, 5, 5, 5, 6])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn levenshtein_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    parallel: bool,
) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_many_to_many(
        &seqs1,
        &seqs2,
        parallel,
        "levenshtein",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, threshold, parallel=false))]
fn levenshtein_neighbor_matrix(
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_matrix(
        &seqs,
        threshold,
        parallel,
        "levenshtein",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, threshold, parallel=false))]
fn levenshtein_neighbor_one_to_many(
    seq: &str,
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<usize>> {
    Ok(_distance::str_neighbor_one_to_many(
        seq,
        &seqs,
        threshold,
        parallel,
        "levenshtein",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, threshold, parallel=false))]
fn levenshtein_neighbor_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_many_to_many(
        &seqs1,
        &seqs2,
        threshold,
        parallel,
        "levenshtein",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn levenshtein_bin_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    parallel: bool,
) -> PyResult<Vec<u32>> {
    Ok(_distance::str_bin_many_to_many(
        &seqs1,
        &seqs2,
        parallel,
        "levenshtein",
    ))
}

/// Compute the Levenshtein distance between two strings using an exponential search.
///
/// The strings must be representable as byte strings. This uses an exponential
/// search to estimate the number of edits. It will be more efficient than
/// levenshtein_distance when the number of edits is small.
///
/// Parameters
/// ----------
/// s1 : str
///     A string.
/// s2 : str
///     Another string
///
/// Returns
/// -------
/// int
///     The Levenshtein distance between the two strings.
///
/// Examples
/// --------
/// >>> s1 = "abcdefg"
/// >>> s2 = "abdcd defgggg"
/// >>> assert(levenshtein_exp(s1, s2) == 6)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn levenshtein_exp(s1: &str, s2: &str) -> PyResult<u32> {
    Ok(triple_accel::levenshtein_exp(
        &s1.as_bytes(),
        &s2.as_bytes(),
    ))
}

/// Compute the Levenshtein distance matrix on an iterable of strings using exponential search.
///
/// The strings must be representable as byte strings. This uses an exponential
/// search to estimate the number of edits. It will be more efficient than
/// levenshtein_exp_matrix when the number of edits is small.
///
/// Parameters
/// ----------
/// seqs : iterable str
///     The strings to be compared. The must have an appropriate representation
///     as byte strings.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seqs = ["Monday", "Tuesday", "Wednesday", "Thursday"]
/// >>> assert(levenshtein_exp_matrix(seqs, parallel=False) == [4, 5, 5, 4, 2, 5])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, parallel=false))]
fn levenshtein_exp_matrix(seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_matrix(
        &seqs,
        parallel,
        "levenshtein_exp",
    ))
}

/// Compute the Levenshtein distance between one string and many others using exponential search.
///
/// This returns the upper right triangle of the distance matrix.
/// The strings must be representable as byte strings. This uses an exponential
/// search to estimate the number of edits. It will be more efficient than
/// levenshtein_one_to_many when the number of edits is small.
///
/// Parameters
/// ----------
/// seq : str
///     The string against which all others will be compared.
/// seqs : iterable of str
///     The other strings being compared.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seq = "Sunday"
/// >>> seqs = ["Monday", "Tuesday", "Wednesday", "Thursday"]
/// >>> assert(levenshtein_one_to_many(seq, seqs, parallel=False) == [2, 3, 5, 4])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, parallel=false))]
fn levenshtein_exp_one_to_many(seq: &str, seqs: Vec<&str>, parallel: bool) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_one_to_many(
        &seq,
        &seqs,
        parallel,
        "levenshtein_exp",
    ))
}

/// Compute the Levenshtein distance between many strings and many others.
///
/// The strings must be representable as byte strings. This uses an exponential
/// search to estimate the number of edits. It will be more efficient than
/// levenshtein_many_to_many when the number of edits is small.
///
/// Parameters
/// ----------
/// seqs1 : iterable of str
///     The first iterable of strings.
/// seqs : iterable of str
///     The second iterable of strings.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The Levenshtein distances among the strings.
///
/// Examples
/// --------
/// >>> seqs1 = ["Sunday", "Saturday"]
/// >>> seqs2 = ["Monday", "Tuesday", "Wednesday"]
/// >>> assert(levenshtein_many_to_many(seqs1, seqs2, parallel=False) == [2, 3, 5, 5, 5, 6])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn levenshtein_exp_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    parallel: bool,
) -> PyResult<Vec<u32>> {
    Ok(_distance::str_cmp_many_to_many(
        &seqs1,
        &seqs2,
        parallel,
        "levenshtein_exp",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, threshold, parallel=false))]
fn levenshtein_exp_neighbor_matrix(
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_matrix(
        &seqs,
        threshold,
        parallel,
        "levenshtein_exp",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, threshold, parallel=false))]
fn levenshtein_exp_neighbor_one_to_many(
    seq: &str,
    seqs: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<usize>> {
    Ok(_distance::str_neighbor_one_to_many(
        seq,
        &seqs,
        threshold,
        parallel,
        "levenshtein_exp",
    ))
}
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, threshold, parallel=false))]
fn levenshtein_exp_neighbor_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    threshold: u32,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::str_neighbor_many_to_many(
        &seqs1,
        &seqs2,
        threshold,
        parallel,
        "levenshtein_exp",
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, parallel=false))]
fn levenshtein_exp_bin_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    parallel: bool,
) -> PyResult<Vec<u32>> {
    Ok(_distance::str_bin_many_to_many(
        &seqs1,
        &seqs2,
        parallel,
        "levenshtein_exp",
    ))
}

/// Compute the distance between two amino acids using BLOSUM62 substitution penalities.
///
/// This function is invariant to case. I.e., lowercase and uppercase residues
/// can be compared.
///
/// Parameters
/// ----------
/// s1 : str
///     An amino acid residue.
/// s2 : str
///     Another amino acid residue.
///
/// Returns
/// -------
/// int
///     The distance between the two residues.
///
/// Examples
/// --------
/// >>> s1 = "C"
/// >>> s2 = "a"
/// >>> assert(amino_acid_distance(s1, s2) == 4)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn amino_acid_distance(s1: &str, s2: &str) -> PyResult<u16> {
    Ok(amino_distances(&s1.as_bytes()[0], &s2.as_bytes()[0]))
}

/// Compute the distance between V genes using precomputed tcrdists.
///
/// Parameters
/// ----------
/// s1 : str
///     A V gene.
/// s2 : str
///     A V gene.
///
/// Returns
/// -------
/// int
///     The distance between two V genes.
///
/// Examples
/// --------
/// >>> s1 = "C"
/// >>> s2 = "a"
/// >>> assert(amino_acid_distance(s1, s2) == 4)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn v_gene_distance(s1: &str, s2: &str) -> PyResult<u16> {
    Ok(gene_distance(&s1.as_bytes(), &s2.as_bytes()))
}

/// Compute the pMHC distance between V alleles using precomputed tcrdists.
///
/// tcrdists were precomputed using ntrim = ctrim = 0, dist_weight = 1,
/// gap_penalty = 4, and fixed_gappos = True.
///
/// Parameters
/// ----------
/// s1 : str
///     A V allele.
/// s2 : str
///     Another V allele.
///
/// Returns
/// -------
/// int
///     The distance between the V alleles' pMHCs.
///
/// Examples
/// --------
/// s1 = "TRBV2*01"
/// s2 = "TRBV6-2*01"
/// assert(phmc_distance(s1, s2) == 16)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn phmc_distance(s1: &str, s2: &str) -> PyResult<u16> {
    Ok(phmc_distances(s1.as_bytes(), s2.as_bytes()))
}

/// Compute the CDR1 distance between V alleles using precomputed tcrdists.
///
/// tcrdists were precomputed using ntrim = ctrim = 0, dist_weight = 1,
/// gap_penalty = 4, and fixed_gappos = True.
///
/// Parameters
/// ----------
/// s1 : str
///     A V allele.
/// s2 : str
///     Another V allele.
///
/// Returns
/// -------
/// int
///     The distance between the V alleles' CDR1s.
///
/// Examples
/// --------
/// s1 = "TRBV2*01"
/// s2 = "TRBV6-2*01"
/// assert(cdr1_distance(s1, s2) == 8)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn cdr1_distance(s1: &str, s2: &str) -> PyResult<u16> {
    Ok(cdr1_distances(s1.as_bytes(), s2.as_bytes()))
}

/// Compute the CDR2 distance between V alleles using precomputed tcrdists.
///
/// tcrdists were precomputed using ntrim = ctrim = 0, dist_weight = 1,
/// gap_penalty = 4, and fixed_gappos = True.
///
/// Parameters
/// ----------
/// s1 : str
///     A V allele.
/// s2 : str
///     Another V allele.
///
/// Returns
/// -------
/// int
///     The distance between the V alleles' CDR2s.
///
/// Examples
/// --------
/// >>> s1 = "TRBV2*01"
/// >>> s2 = "TRBV6-2*01"
/// >>> assert(cdr2_distance(s1, s2) == 24)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
fn cdr2_distance(s1: &str, s2: &str) -> PyResult<u16> {
    Ok(cdr2_distances(s1.as_bytes(), s2.as_bytes()))
}

/// Compute the tcrdist between two strings.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// s1 : str
///     A string. Ideally, this should be a string of amino acid residues.
/// s2 : str
///     A string. Ideally, this should be a string of amino acid residues.
/// dist_weight : int, default 1
///     A weight applied to the mismatch distances. This weight is not applied
///     to the gap penalties.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
///
/// Returns
/// -------
/// int
///     The tcrdist between two strings.
///
/// Examples
/// --------
/// >>> s1 = "CASRTGTVYEQYF"
/// >>> s2 = "CASSTLDRVYNSPLHF"
/// >>> dist_weight = 1
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist(s1, s2, dist_weight, gap_penalty, ntrim, ctrim, fixed_gappos)
/// >>> assert(dist == 40)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (s1, s2, dist_weight=1, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false))]
fn tcrdist(
    s1: &str,
    s2: &str,
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
) -> PyResult<u16> {
    Ok(_distance::tcrdist(
        s1.as_bytes(),
        s2.as_bytes(),
        dist_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
    ))
}

/// Compute the tcrdist matrix on an iterable of strings.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs : iterable of str
///     Iterable of strings.
/// dist_weight : int, default 1
///     A weight applied to the mismatch distances. This weight is not applied
///     to the gap penalties.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the strings.
///
/// Examples
/// --------
/// >>> seqs = ["CASRTGTVYEQYF", "CASSTLDRVYNSPLHF", "CASSESGGQVDTQYF", "CASSPTGPTDTQYF"]
/// >>> dist_weight = 1
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_matrix(seqs, dist_weight, gap_penalty, ntrim, ctrim, fixed_gappos,
/// parallel=False)
/// >>> assert(dist == [40, 28, 28, 40, 40, 19])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, dist_weight=1, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_matrix(
    seqs: Vec<&str>,
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_matrix(
        &seqs,
        dist_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between one string and many others.
///
/// This returns the upper right triangle of the distance matrix.
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seq : str
///     The string against which all others will be compared.
/// seqs : str
///     The other strings being compared.
/// dist_weight : int, default 1
///     A weight applied to the mismatch distances. This weight is not applied
///     to the gap penalties.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the strings.
///
/// Examples
/// --------
/// >>> seq = "CASSYPIEGGRAFTGELFF"
/// >>> seqs = ["CASRTGTVYEQYF", "CASSTLDRVYNSPLHF", "CASSESGGQVDTQYF", "CASSPTGPTDTQYF"]
/// >>> dist_weight = 1
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_one_to_many(seq, seqs, dist_weight, gap_penalty, ntrim, ctrim, fixed_gappos,
/// parallel=False)
/// >>> assert(dist == [52, 41, 52, 48])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, dist_weight=1, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_one_to_many(
    seq: &str,
    seqs: Vec<&str>,
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_one_to_many(
        seq,
        &seqs,
        dist_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between many strings and many others.
///
/// The strings must be representable as byte strings.
///
/// Parameters
/// ----------
/// seqs1 : iterable of str
///     The first iterable of strings.
/// seqs2 : iterable of str
///     The other iterable of strings.
/// dist_weight : int, default 1
///     A weight applied to the mismatch distances. This weight is not applied
///     to the gap penalties.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the strings.
///
/// Examples
/// --------
/// >>> seqs1 = ["CASSPTGPTDTQYF", "CASSYPIEGGRAFTGELFF"]
/// >>> seqs2 = ["CASRTGTVYEQYF", "CASSTLDRVYNSPLHF", "CASSESGGQVDTQYF"]
/// >>> dist_weight = 1
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_many_to_many(seqs1, seqs2, dist_weight, gap_penalty, ntrim, ctrim, fixed_gappos, parallel=False)
/// >>> assert(dist == [52, 41, 52, 48])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, dist_weight=1, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_many_to_many(
    seqs1: Vec<&str>,
    seqs2: Vec<&str>,
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_many_to_many(
        &seqs1,
        &seqs2,
        dist_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between two CDR3-V allele pairs.
///
/// This incorporates differences between the pMHC, CDR1, CDR2, and CDR3.
///
/// Parameters
/// ----------
/// s1 : iterable of str
///     An iterable of the CDR3 amino acid sequence and V allele.
/// s2 : iterable of str
///     An iterable of the CDR3 amino acid sequence and V allele.
/// phmc_weight : int, default 1
///     How much the difference in pMHCs contributes to the distance.
/// cdr1_weight : int, default 1
///     How much the difference in CDR1s contributes to the distance.
/// cdr2_weight : int, default 1
///     How much the difference in CDR2s contributes to the distance.
/// cdr3_weight : int, default 3
///     How much the difference in CDR3s contributes to the distance.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
///
/// Returns
/// -------
/// int
///     The tcrdist between two CDR3-V allele pairs.
///
/// Examples
/// --------
/// >>> s1 = ["CASRTGTVYEQYF", "TRBV2*01"]
/// >>> s2 = ["CASSTLDRVYNSPLHF", "TRBV6-2*01"]
/// >>> phmc_weight = 1
/// >>> cdr1_weight = 1
/// >>> cdr2_weight = 1
/// >>> cdr3_weight = 3
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_allele(s1, s2, phmc_weight, cdr1_weight, cdr2_weight, cdr3_weight, gap_penalty, ntrim, ctrim, fixed_gappos)
/// >>> assert(dist == 168)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (s1, s2, phmc_weight=1, cdr1_weight=1, cdr2_weight=1, cdr3_weight=3, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false))]
fn tcrdist_allele(
    s1: [&str; 2],
    s2: [&str; 2],
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
) -> PyResult<u16> {
    Ok(_distance::tcrdist_allele(
        s1,
        s2,
        phmc_weight,
        cdr1_weight,
        cdr2_weight,
        cdr3_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
    ))
}

/// Compute the tcrdist matrix on an iterable of CDR3-V allele pairs.
///
/// This returns the upper right triangle of the distance matrix.
/// This incorporates differences between the pMHC, CDR1, CDR2, and CDR3.
///
/// Parameters
/// ----------
/// seqs : iterable of iterable of str
///     An iterable containing iterables of CDR3 amino acid sequences and V alleles.
/// phmc_weight : int, default 1
///     How much the difference in pMHCs contributes to the distance.
/// cdr1_weight : int, default 1
///     How much the difference in CDR1s contributes to the distance.
/// cdr2_weight : int, default 1
///     How much the difference in CDR2s contributes to the distance.
/// cdr3_weight : int, default 3
///     How much the difference in CDR3s contributes to the distance.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, default False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the CDR3-V allele pairs.
///
/// Examples
/// --------
/// >>> seqs = [["CASRTGTVYEQYF", "TRBV2*01"], ["CASSTLDRVYNSPLHF", "TRBV6-2*01"], ["CASSESGGQVDTQYF", "TRBV6-4*01"], ["CASSPTGPTDTQYF", "TRBV18*01"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5*01"]]
/// >>> phmc_weight = 1
/// >>> cdr1_weight = 1
/// >>> cdr2_weight = 1
/// >>> cdr3_weight = 3
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_allele_matrix(seqs, phmc_weight, cdr1_weight, cdr2_weight, cdr3_weight, gap_penalty, ntrim, ctrim, fixed_gappos, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203, 163, 169, 148, 116, 189, 198])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, phmc_weight=1, cdr1_weight=1, cdr2_weight=1, cdr3_weight=3, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_allele_matrix(
    seqs: Vec<[&str; 2]>,
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_allele_matrix(
        &seqs,
        phmc_weight,
        cdr1_weight,
        cdr2_weight,
        cdr3_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between one CDR3-V allele pair and many others.
///
/// This incorporates differences between the pMHC, CDR1, CDR2, and CDR3.
///
/// Parameters
/// ----------
/// s1 : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V allele.
/// seqs : iterable of iterable of str
///     An iterable of iterables containing pairs of CDR3 amino acid sequences and V alleles.
/// phmc_weight : int, default 1
///     How much the difference in pMHCs contributes to the distance.
/// cdr1_weight : int, default 1
///     How much the difference in CDR1s contributes to the distance.
/// cdr2_weight : int, default 1
///     How much the difference in CDR2s contributes to the distance.
/// cdr3_weight : int, default 3
///     How much the difference in CDR3s contributes to the distance.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the CDR3-V allele sequences.
///
/// Examples
/// --------
/// >>> seq = ["CASRTGTVYEQYF", "TRBV2*01"]
/// >>> seqs = [["CASSTLDRVYNSPLHF", "TRBV6-2*01"], ["CASSESGGQVDTQYF", "TRBV6-4*01"], ["CASSPTGPTDTQYF", "TRBV18*01"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5*01"]]
/// >>> phmc_weight = 1
/// >>> cdr1_weight = 1
/// >>> cdr2_weight = 1
/// >>> cdr3_weight = 3
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_allele_one_to_many(seq, seqs, phmc_weight, cdr1_weight, cdr2_weight, cdr3_weight, gap_penalty, ntrim, ctrim, fixed_gappos, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, phmc_weight=1, cdr1_weight=1, cdr2_weight=1, cdr3_weight=3, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_allele_one_to_many(
    seq: [&str; 2],
    seqs: Vec<[&str; 2]>,
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_allele_one_to_many(
        seq,
        &seqs,
        phmc_weight,
        cdr1_weight,
        cdr2_weight,
        cdr3_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between many CDR3-V allele pairs and many others.
///
/// This incorporates differences between the pMHC, CDR1, CDR2, and CDR3.
///
/// Parameters
/// ----------
/// seqs1 : iterable of iterable of str
///     An iterable of iterables containing pairs of CDR3 amino acid sequences and V alleles.
/// seqs2 : iterable of iterable of str
///     An iterable of iterables containing pairs of CDR3 amino acid sequences and V alleles.
/// phmc_weight : int, default 1
///     How much the difference in pMHCs contributes to the distance.
/// cdr1_weight : int, default 1
///     How much the difference in CDR1s contributes to the distance.
/// cdr2_weight : int, default 1
///     How much the difference in CDR2s contributes to the distance.
/// cdr3_weight : int, default 3
///     How much the difference in CDR3s contributes to the distance.
/// gap_penalty : int, default 4
///     The penalty given to the difference in length of the strings.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// fixed_gappos : bool, False
///     If True, insert gaps at a fixed position after the cysteine residue
///     starting the CDR3 (typically position 6). If False, find the "optimal"
///     position for inserting the gaps to make up the difference in length.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the CDR3-V allele pairs.
///
/// Examples
/// --------
/// >>> seqs1 = [["CASRTGTVYEQYF", "TRBV2*01"], ["CASSYSEEPSSPLHF", "TRBV6-6*01"]]
/// >>> seqs2 = [["CASSTLDRVYNSPLHF", "TRBV6-2*01"], ["CASSESGGQVDTQYF", "TRBV6-4*01"], ["CASSPTGPTDTQYF", "TRBV18*01"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5*01"]]
/// >>> phmc_weight = 1
/// >>> cdr1_weight = 1
/// >>> cdr2_weight = 1
/// >>> cdr3_weight = 3
/// >>> gap_penalty = 4
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> fixed_gappos = False
/// >>> dist = tcrdist_allele_many_to_many(seqs1, seqs2, phmc_weight, cdr1_weight, cdr2_weight, cdr3_weight, gap_penalty, ntrim, ctrim, fixed_gappos, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203, 104, 125, 143, 121])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, phmc_weight=1, cdr1_weight=1, cdr2_weight=1, cdr3_weight=3, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=false, parallel=false))]
fn tcrdist_allele_many_to_many(
    seqs1: Vec<[&str; 2]>,
    seqs2: Vec<[&str; 2]>,
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_allele_many_to_many(
        &seqs1,
        &seqs2,
        phmc_weight,
        cdr1_weight,
        cdr2_weight,
        cdr3_weight,
        gap_penalty,
        ntrim,
        ctrim,
        fixed_gappos,
        parallel,
    ))
}

/// Compute the tcrdist between two CDR3-V gene pairs.
///
/// Parameters
/// ----------
/// s1 : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V gene pair.
/// s2 : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V gene pair.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
///
/// Returns
/// -------
/// int
///     The tcrdist between two CDR3-V gene pairs.
///
/// Examples
/// --------
/// >>> s1 = ["CASRTGTVYEQYF", "TRBV2"]
/// >>> s2 = ["CASSTLDRVYNSPLHF", "TRBV6-2"]
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> dist = tcrdist_gene(s1, s2, ntrim, ctrim)
/// >>> assert(dist == 168)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (s1, s2, ntrim=3, ctrim=2))]
fn tcrdist_gene(s1: [&str; 2], s2: [&str; 2], ntrim: usize, ctrim: usize) -> PyResult<u16> {
    Ok(_distance::tcrdist_gene(s1, s2, ntrim, ctrim))
}

/// Compute the tcrdist matrix on an iterable of CDR3-V gene pairs.
///
/// This returns the upper right triangle of the distance matrix.
///
/// Parameters
/// ----------
/// seqs : iterable of iterable of str
///     An iterable containing iterables of CDR3 amino acid sequence and V gene pairs.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdist among the CDR3-V gene pairs.
///
/// Examples
/// --------
/// >>> seqs = [["CASRTGTVYEQYF", "TRBV2"], ["CASSTLDRVYNSPLHF", "TRBV6-2"], ["CASSESGGQVDTQYF", "TRBV6-4"], ["CASSPTGPTDTQYF", "TRBV18"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5"]]
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> dist = tcrdist_gene_matrix(seqs, ntrim, ctrim, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203, 163, 169, 148, 116, 189, 198])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_matrix(
    seqs: Vec<[&str; 2]>,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_gene_matrix(
        &seqs, ntrim, ctrim, parallel,
    ))
}

/// Compute the tcrdist between one CDR3-V gene pair and many others.
///
/// Parameters
/// ----------
/// seq : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V allele pair.
/// seqs : iterable of iterable of str
///     An iterable containing iterables of CDR3 amino acid sequence and V gene pairs.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the CDR3-V gene pairs.
///
/// Examples
/// --------
/// >>> seq = ["CASRTGTVYEQYF", "TRBV2"]
/// >>> seqs = [["CASSTLDRVYNSPLHF", "TRBV6-2"], ["CASSESGGQVDTQYF", "TRBV6-4"], ["CASSPTGPTDTQYF", "TRBV18"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5"]]
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> dist = tcrdist_gene_one_to_many(seq, seqs, ntrim, ctrim, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_one_to_many(
    seq: [&str; 2],
    seqs: Vec<[&str; 2]>,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_gene_one_to_many(
        seq, &seqs, ntrim, ctrim, parallel,
    ))
}

/// Compute the tcrdist between many CDR3-V gene pairs and many others.
///
/// Parameters
/// ----------
/// seqs1 : iterable of iterable of str
///     An iterable containing iterables of CDR3 amino acid sequence and V gene pairs.
/// seqs2 : iterable of iterable of str
///     An iterable containing iterables of CDR3 amino acid sequence and V gene pairs.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
/// parallel: bool, default False
///     Bool to specify if computation should be parallelized.
///
/// Returns
/// -------
/// list of int
///     The tcrdists among the CDR3-V gene pairs.
///
/// Examples
/// --------
/// >>> seqs1 = [["CASRTGTVYEQYF", "TRBV2"], ["CASSYSEEPSSPLHF", "TRBV6-6"]]
/// >>> seqs2 = [["CASSTLDRVYNSPLHF", "TRBV6-2"], ["CASSESGGQVDTQYF", "TRBV6-4"], ["CASSPTGPTDTQYF", "TRBV18"], ["CASSYPIEGGRAFTGELFF", "TRBV6-5"]]
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> dist = tcrdist_gene_many_to_many(seqs1, seqs2, ntrim, ctrim, parallel=False)
/// >>> assert(dist == [168, 142, 134, 203, 104, 125, 143, 121])
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_many_to_many(
    seqs1: Vec<[&str; 2]>,
    seqs2: Vec<[&str; 2]>,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<u16>> {
    Ok(_distance::tcrdist_gene_many_to_many(
        &seqs1, &seqs2, ntrim, ctrim, parallel,
    ))
}

/// Compute whether two CDR3-V gene pairs are neighbors with tcrdist_gene.
///
/// This function is quicker than using the tcrdist_gene function since it computes
/// whether the V genes are within the distance threshold and whether the difference
/// in lengths won't incur a penalty larger than the distance threshold. With these
/// two checks, many unnecessary calculations are avoided.
///
/// Parameters
/// ----------
/// s1 : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V gene pair.
/// s2 : iterable of str
///     An iterable containing a CDR3 amino acid sequence and V gene pair.
/// threshold: int
///     The distance threshold that will be used to call sequences neighbors.
/// ntrim : int, default 3
///     The position at which the distance calculation will begin.
///     This parameter must be >= 0.
/// ctrim : int, default 2
///     The position, counted from the end, at which the calculation will end.
///     This parameter must be >= 0.
///
/// Returns
/// -------
/// bool
///     Whether the two CDR3-V gene pairs are have tcrdist within the threshold.
///
/// Examples
/// --------
/// >>> s1 = ["CASRTGTVYEQYF", "TRBV2"]
/// >>> s2 = ["CASSTLDRVYNSPLHF", "TRBV6-2"]
/// >>> threshold = 20
/// >>> ntrim = 3
/// >>> ctrim = 2
/// >>> are_neighbors = tcrdist_gene(s1, s2, threshold, ntrim, ctrim)
/// >>> assert(are_neighbors == False)
#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (s1, s2, threshold, ntrim=3, ctrim=2))]
fn tcrdist_gene_neighbor(
    s1: [&str; 2],
    s2: [&str; 2],
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
) -> PyResult<bool> {
    Ok(_distance::tcrdist_gene_neighbor(
        s1, s2, threshold, ntrim, ctrim,
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs, threshold, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_neighbor_matrix(
    seqs: Vec<[&str; 2]>,
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::tcrdist_gene_neighbor_matrix(
        &seqs, threshold, ntrim, ctrim, parallel,
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seq, seqs, threshold, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_neighbor_one_to_many(
    seq: [&str; 2],
    seqs: Vec<[&str; 2]>,
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<usize>> {
    Ok(_distance::tcrdist_gene_neighbor_one_to_many(
        seq, &seqs, threshold, ntrim, ctrim, parallel,
    ))
}

#[cfg(all(feature = "pyo3"))]
#[pyfunction]
#[pyo3(signature = (seqs1, seqs2, threshold, ntrim=3, ctrim=2, parallel=false))]
fn tcrdist_gene_neighbor_many_to_many(
    seqs1: Vec<[&str; 2]>,
    seqs2: Vec<[&str; 2]>,
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> PyResult<Vec<[usize; 2]>> {
    Ok(_distance::tcrdist_gene_neighbor_many_to_many(
        &seqs1, &seqs2, threshold, ntrim, ctrim, parallel,
    ))
}

#[pymodule]
#[pyo3(name = "tcrdist_rs")]
pub fn tcrdist_rs(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hamming, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(hamming_neighbor_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_neighbor_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_neighbor_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(levenshtein_exp, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(levenshtein_exp_neighbor_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_neighbor_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_neighbor_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(levenshtein, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(levenshtein_neighbor_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_neighbor_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_neighbor_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(hamming_bin_many_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_exp_bin_many_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(levenshtein_bin_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(amino_acid_distance, m)?)?;
    m.add_function(wrap_pyfunction!(v_gene_distance, m)?)?;
    m.add_function(wrap_pyfunction!(phmc_distance, m)?)?;
    m.add_function(wrap_pyfunction!(cdr1_distance, m)?)?;
    m.add_function(wrap_pyfunction!(cdr2_distance, m)?)?;

    m.add_function(wrap_pyfunction!(tcrdist, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(tcrdist_allele, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_allele_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_allele_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_allele_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(tcrdist_gene, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_many_to_many, m)?)?;

    m.add_function(wrap_pyfunction!(tcrdist_gene_neighbor, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_neighbor_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_neighbor_one_to_many, m)?)?;
    m.add_function(wrap_pyfunction!(tcrdist_gene_neighbor_many_to_many, m)?)?;

    Ok(())
}
