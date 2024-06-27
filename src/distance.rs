use crate::match_table;

use once_cell::sync::Lazy;
use rayon::prelude::*;
use std::cmp;
use std::collections::HashMap;
use triple_accel::*;

pub static POOL: Lazy<rayon::ThreadPool> = Lazy::new(|| {
    rayon::ThreadPoolBuilder::new()
        .num_threads(
            std::thread::available_parallelism()
                .unwrap_or(std::num::NonZeroUsize::new(1).unwrap())
                .get(),
        )
        .thread_name(move |i| format!("{}-{}", "tcrdist", i))
        .build()
        .expect("could not spawn threads")
});

/// Compute the tcrdist between two byte strings.
///
/// * `s1` - A byte string.
/// * `s2` - Another byte string.
/// * `dist_weight` - The weight applied to the mismatches and not the gaps.
/// * `gap_penalty` - The penalty applied to the difference in the strings.
/// * `ntrim` - The position at which the distance computation will begin.
/// * `ctrim` - The position, counted from the end, at which the calculation will end.
/// * `fixed_gappos` - If true, insert gaps at a fixed position ater the cysteine residue starting
/// the CDR3 (typically after position 6). If false, find the "optimal" position for inserting the
/// gaps to make up the difference in length.
///
/// # Examples
///
/// ```
/// let s1 = b"CASRTGTVYEQYF";
/// let s2 = b"CASSTLDRVYNSPLHF";
///
/// let dist_weight = 1;
/// let gap_penalty = 4;
/// let ntrim = 3;
/// let ctrim = 2;
/// let fixed_gappos = false;
///
/// let dist = tcrdist(s1, s2, dist_weight, gap_penalty, ntrim, ctrim,
///                    fixed_gappos);
/// assert_eq!(dist, 40);
/// ```
pub fn tcrdist(
    s1: &[u8],
    s2: &[u8],
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
) -> u16 {
    let s1_len: usize = s1.len();
    let s2_len: usize = s2.len();
    if s1_len == s2_len {
        let mut distance: u16 = 0;
        for idx in ntrim..s1_len - ctrim {
            distance += match_table::amino_distances(&s1[idx], &s2[idx]);
        }
        return distance;
    }

    let len_diff: usize;
    let short_len: usize;

    if s1_len > s2_len {
        short_len = s2_len;
        len_diff = s1_len - s2_len;
    } else {
        short_len = s1_len;
        len_diff = s2_len - s1_len;
    };

    let i_short_len: i8 = short_len as i8;

    let min_gappos: i8;
    let max_gappos: i8;

    if fixed_gappos {
        min_gappos = cmp::min(6, (i_short_len + 1) / 2);
        max_gappos = min_gappos;
    } else {
        if 10 > short_len {
            let num: i8 = if i_short_len % 2 != 0 { 1 } else { 0 };
            max_gappos = i_short_len / 2 + num;
            min_gappos = max_gappos - num;
        } else {
            max_gappos = i_short_len - 5;
            min_gappos = 5;
        }
    }

    let u_min_gappos: usize = min_gappos as usize;
    let u_max_gappos: usize = max_gappos as usize;
    let mut min_dist: u16 = u16::MAX;
    for gappos in u_min_gappos..u_max_gappos + 1 {
        let mut tmp_dist: u16 = 0;

        for idx in ntrim..gappos {
            tmp_dist += match_table::amino_distances(&s1[idx], &s2[idx]);
        }

        for c_i in ctrim..short_len - gappos {
            let idx1: usize = s1_len - 1 - c_i;
            let idx2: usize = s2_len - 1 - c_i;
            tmp_dist += match_table::amino_distances(&s1[idx1], &s2[idx2]);
        }
        if tmp_dist < min_dist {
            min_dist = tmp_dist;
        }
        if min_dist == 0 {
            break;
        }
    }
    return min_dist * dist_weight + (len_diff as u16) * gap_penalty;
}

/// Compute the tcrdist matrix on an iterable of strings.
pub fn tcrdist_matrix(
    seqs: &[&str],
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let seqs_len: usize = seqs.len();
        let num_combinations: usize = seqs_len * (seqs_len - 1) / 2;
        let mut dists: Vec<u16> = vec![0; num_combinations];
        let mut counter: usize = 0;

        for (i, &s1) in seqs.iter().enumerate() {
            for &s2 in seqs[i + 1..].iter() {
                dists[counter] = tcrdist(
                    s1.as_bytes(),
                    s2.as_bytes(),
                    dist_weight,
                    gap_penalty,
                    ntrim,
                    ctrim,
                    fixed_gappos,
                );
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(i, &s1)| {
                    seqs[i + 1..]
                        .iter()
                        .map(|&s2| {
                            tcrdist(
                                s1.as_bytes(),
                                s2.as_bytes(),
                                dist_weight,
                                gap_penalty,
                                ntrim,
                                ctrim,
                                fixed_gappos,
                            )
                        })
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute the tcrdist between one string and many other strings.
pub fn tcrdist_one_to_many(
    seq: &str,
    seqs: &[&str],
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    let seq_bytes = seq.as_bytes();
    if parallel == false {
        seqs.iter()
            .map(|&s| {
                tcrdist(
                    seq_bytes,
                    s.as_bytes(),
                    dist_weight,
                    gap_penalty,
                    ntrim,
                    ctrim,
                    fixed_gappos,
                )
            })
            .collect()
    } else {
        seqs.par_iter()
            .map(|&s| {
                tcrdist(
                    seq_bytes,
                    s.as_bytes(),
                    dist_weight,
                    gap_penalty,
                    ntrim,
                    ctrim,
                    fixed_gappos,
                )
            })
            .collect()
    }
}

/// Compute the tcrdist between many strings and many other strings.
pub fn tcrdist_many_to_many(
    seqs1: &[&str],
    seqs2: &[&str],
    dist_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let mut dists: Vec<u16> = vec![0; seqs1.len() * seqs2.len()];
        let mut counter: usize = 0;

        for &s1 in seqs1.iter() {
            let s1_bytes = s1.as_bytes();
            for &s2 in seqs2.iter() {
                dists[counter] = tcrdist(
                    s1_bytes,
                    s2.as_bytes(),
                    dist_weight,
                    gap_penalty,
                    ntrim,
                    ctrim,
                    fixed_gappos,
                );
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .flat_map(|&s1| {
                    let s1_bytes = s1.as_bytes();
                    seqs2
                        .iter()
                        .map(|&s2| {
                            tcrdist(
                                s1_bytes,
                                s2.as_bytes(),
                                dist_weight,
                                gap_penalty,
                                ntrim,
                                ctrim,
                                fixed_gappos,
                            )
                        })
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute the distance between V alleles which are written as byte strings.
///
/// This function is memoized to speed up V alleles distance computations further.
///
/// TODO Memoize this function so only 1 lookup compared to 3.
fn v_allele_distance(
    v_allele_1: &[u8],
    v_allele_2: &[u8],
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
) -> u16 {
    if v_allele_1 == v_allele_2 {
        0
    } else {
        match_table::phmc_distances(v_allele_1, v_allele_2) * phmc_weight
            + match_table::cdr1_distances(v_allele_1, v_allele_2) * cdr1_weight
            + match_table::cdr2_distances(v_allele_1, v_allele_2) * cdr2_weight
    }
}

/// Compute the full tcrdist between two CDR3-V allele pairs.
///
/// This incorporates differences between the pMHC, CDR1, CDR2, and CDR3.
///
/// * `s1` - An array of a CDR3 amino acid sequence and V gene allele.
/// * `s2` - Another array of a CDR3 amino acid sequence and V gene allele.
/// * `phmc_weight` - The weight applied to the pMHC distance.
/// * `cdr1_weight` - The weight applied to the CDR1 distance.
/// * `cdr2_weight` - The weight applied to the CDR2 distance.
/// * `cdr3_weight` - The weight applied to the CDR3 distance.
/// * `gap_penalty` - The penalty applied to the difference in the strings.
/// * `ntrim` - The position at which the distance computation will begin.
/// * `ctrim` - The position, counted from the end, at which the calculation will end.
/// * `fixed_gappos` - If true, insert gaps at a fixed position ater the cysteine residue starting
/// the CDR3 (typically after position 6). If false, find the "optimal" position for inserting the
/// gaps to make up the difference in length.
///
/// # Examples
/// ```
/// let s1 = ["CASRTGTVYEQYF", "TRBV2*01"];
/// let s2 = ["CASSTLDRVYNSPLHF", "TRBV6-2*01"];
///
/// let phmc_weight = 1;
/// let cdr1_weight = 1;
/// let cdr2_weight = 1;
/// let cdr3_weight = 3;
/// let gap_penalty = 4;
/// let ntrim = 3;
/// let ctrim = 2;
/// let fixed_gappos = false;
///
/// dist = tcrdist_allele(s1, s2, phmc_weight, cdr1_weight, cdr2_weight, cdr3_weight,
///                       gap_penalty, ntrim, ctrim, fixed_gappos);
///
/// assert_eq!(dist, 168);
/// ```
pub fn tcrdist_allele(
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
) -> u16 {
    v_allele_distance(
        s1[1].as_bytes(),
        s2[1].as_bytes(),
        phmc_weight,
        cdr1_weight,
        cdr2_weight,
    ) + cdr3_weight
        * tcrdist(
            s1[0].as_bytes(),
            s2[0].as_bytes(),
            1,
            gap_penalty,
            ntrim,
            ctrim,
            fixed_gappos,
        )
}

/// Compute the full tcrdist matrix on an iterable of CDR3-V allele arrays.
pub fn tcrdist_allele_matrix(
    seqs: &[[&str; 2]],
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let seqs_len: usize = seqs.len();
        let num_combinations: usize = seqs_len * (seqs_len - 1) / 2;
        let mut dists: Vec<u16> = vec![0; num_combinations];
        let mut counter: usize = 0;

        for (i, &s1) in seqs.iter().enumerate() {
            for &s2 in seqs[i + 1..].iter() {
                dists[counter] = tcrdist_allele(
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
                );
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(i, &s1)| {
                    seqs[i + 1..]
                        .iter()
                        .map(|&s2| {
                            tcrdist_allele(
                                s1,
                                s2,
                                phmc_weight,
                                cdr1_weight,
                                cdr3_weight,
                                cdr3_weight,
                                gap_penalty,
                                ntrim,
                                ctrim,
                                fixed_gappos,
                            )
                        })
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute the full tcrdist between a CDR3-V allele array and many others.
pub fn tcrdist_allele_one_to_many(
    seq: [&str; 2],
    seqs: &[[&str; 2]],
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        seqs.iter()
            .map(|&s| {
                tcrdist_allele(
                    seq,
                    s,
                    phmc_weight,
                    cdr1_weight,
                    cdr2_weight,
                    cdr3_weight,
                    gap_penalty,
                    ntrim,
                    ctrim,
                    fixed_gappos,
                )
            })
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .map(|&s| {
                    tcrdist_allele(
                        seq,
                        s,
                        phmc_weight,
                        cdr1_weight,
                        cdr2_weight,
                        cdr3_weight,
                        gap_penalty,
                        ntrim,
                        ctrim,
                        fixed_gappos,
                    )
                })
                .collect::<Vec<u16>>()
        })
    }
}

/// Compute the full tcrdist between many CDR3-V allele arrays and many others.
pub fn tcrdist_allele_many_to_many(
    seqs1: &[[&str; 2]],
    seqs2: &[[&str; 2]],
    phmc_weight: u16,
    cdr1_weight: u16,
    cdr2_weight: u16,
    cdr3_weight: u16,
    gap_penalty: u16,
    ntrim: usize,
    ctrim: usize,
    fixed_gappos: bool,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let mut dists: Vec<u16> = vec![0; seqs1.len() * seqs2.len()];
        let mut counter: usize = 0;

        for &s1 in seqs1.iter() {
            for &s2 in seqs2.iter() {
                dists[counter] = tcrdist_allele(
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
                );
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .flat_map(|&s1| {
                    seqs2
                        .iter()
                        .map(|&s2| {
                            tcrdist_allele(
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
                            )
                        })
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute the full tcrdist between two CDR3-V gene pairs.
///
/// * `s1` - An array of a CDR3 amino acid sequence and V gene.
/// * `s2` - Another array of a CDR3 amino acid sequence and V gene.
/// * `ntrim` - The position at which the distance computation will begin.
/// * `ctrim` - The position, counted from the end, at which the calculation will end.
///
/// # Examples
/// ```
/// let s1 = ["CASRTGTVYEQYF", "TRBV2"];
/// let s2 = ["CASSTLDRVYNSPLHF", "TRBV6-2"];
///
/// let ntrim = 3;
/// let ctrim = 2;
///
/// dist = tcrdist_allele(s1, s2, ntrim, ctrim);
///
/// assert_eq!(dist, 168);
/// ```
pub fn tcrdist_gene(s1: [&str; 2], s2: [&str; 2], ntrim: usize, ctrim: usize) -> u16 {
    match_table::gene_distance(s1[1].as_bytes(), s2[1].as_bytes())
        + 3 * tcrdist(
            s1[0].as_bytes(),
            s2[0].as_bytes(),
            1,
            4,
            ntrim,
            ctrim,
            false,
        )
}

/// Compute the full tcrdist matrix on an iterable of CDR3-V gene arrays.
pub fn tcrdist_gene_matrix(
    seqs: &[[&str; 2]],
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let seqs_len: usize = seqs.len();
        let num_combinations: usize = seqs_len * (seqs_len - 1) / 2;
        let mut dists: Vec<u16> = vec![0; num_combinations];
        let mut counter: usize = 0;

        for (i, &s1) in seqs.iter().enumerate() {
            for &s2 in seqs[i + 1..].iter() {
                dists[counter] = tcrdist_gene(s1, s2, ntrim, ctrim);
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(i, &s1)| {
                    seqs[i + 1..]
                        .iter()
                        .map(|&s2| tcrdist_gene(s1, s2, ntrim, ctrim))
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute the full tcrdist between a CDR3-V gene array and many others.
pub fn tcrdist_gene_one_to_many(
    seq: [&str; 2],
    seqs: &[[&str; 2]],
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        seqs.iter()
            .map(|&s| tcrdist_gene(seq, s, ntrim, ctrim))
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .map(|&s| tcrdist_gene(seq, s, ntrim, ctrim))
                .collect::<Vec<u16>>()
        })
    }
}

/// Compute the full tcrdist between many CDR3-V gene arrays and many others.
pub fn tcrdist_gene_many_to_many(
    seqs1: &[[&str; 2]],
    seqs2: &[[&str; 2]],
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<u16> {
    if parallel == false {
        let seqs1_len: usize = seqs1.len();
        let seqs2_len: usize = seqs2.len();
        let mut dists: Vec<u16> = vec![0; seqs1_len * seqs2_len];
        let mut counter: usize = 0;

        for &s1 in seqs1.iter() {
            for &s2 in seqs2.iter() {
                dists[counter] = tcrdist_gene(s1, s2, ntrim, ctrim);
                counter += 1;
            }
        }
        dists
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .flat_map(|&s1| {
                    seqs2
                        .iter()
                        .map(|&s2| tcrdist_gene(s1, s2, ntrim, ctrim))
                        .collect::<Vec<u16>>()
                })
                .collect()
        })
    }
}

/// Compute whether two CDR3-V gene arrays are tcrdist-gene neighbors.
pub fn tcrdist_gene_neighbor(
    s1: [&str; 2],
    s2: [&str; 2],
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
) -> bool {
    let s1_bytes: &[u8] = s1[0].as_bytes();
    let s2_bytes: &[u8] = s2[0].as_bytes();
    let s1_len: usize = s1_bytes.len();
    let s2_len: usize = s2_bytes.len();
    let len_diff: u16 = if s1_len > s2_len {
        (s1_len - s2_len) as u16
    } else {
        (s2_len - s1_len) as u16
    };

    // Gap penalty is 12. Stop computation if lengths are too different.
    if len_diff * 12 > threshold {
        return false;
    }

    // Stop computation if V gene distance and length difference are too different.
    let v_gene_dist = match_table::gene_distance(s1[1].as_bytes(), s2[1].as_bytes());
    if v_gene_dist + len_diff > threshold {
        return false;
    }

    (v_gene_dist + 3 * tcrdist(s1_bytes, s2_bytes, 1, 4, ntrim, ctrim, false)) <= threshold
}

pub fn tcrdist_gene_neighbor_matrix(
    seqs: &[[&str; 2]],
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<[usize; 2]> {
    if parallel == false {
        seqs.iter()
            .enumerate()
            .flat_map(|(idx, &s1)| {
                seqs[idx + 1..]
                    .iter()
                    .enumerate()
                    .filter(move |(_, &s2)| tcrdist_gene_neighbor(s1, s2, threshold, ntrim, ctrim))
                    .map(move |(jdx, _)| [idx, idx + 1 + jdx])
            })
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(idx, &s1)| {
                    seqs[idx + 1..]
                        .iter()
                        .enumerate()
                        .filter(|(_, &s2)| tcrdist_gene_neighbor(s1, s2, threshold, ntrim, ctrim))
                        .map(move |(jdx, _)| [idx, idx + 1 + jdx])
                        .collect::<Vec<[usize; 2]>>()
                })
                .collect()
        })
    }
}

pub fn tcrdist_gene_neighbor_one_to_many(
    seq: [&str; 2],
    seqs: &[[&str; 2]],
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<usize> {
    if parallel == false {
        seqs.iter()
            .enumerate()
            .filter(|(_, &s)| tcrdist_gene_neighbor(seq, s, threshold, ntrim, ctrim))
            .map(|(idx, _)| idx)
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .filter(|(_, &s)| tcrdist_gene_neighbor(seq, s, threshold, ntrim, ctrim))
                .map(|(idx, _)| idx)
                .collect::<Vec<usize>>()
        })
    }
}

pub fn tcrdist_gene_neighbor_many_to_many(
    seqs1: &[[&str; 2]],
    seqs2: &[[&str; 2]],
    threshold: u16,
    ntrim: usize,
    ctrim: usize,
    parallel: bool,
) -> Vec<[usize; 2]> {
    if parallel == false {
        seqs1
            .iter()
            .enumerate()
            .flat_map(|(idx, &s1)| {
                seqs2
                    .iter()
                    .enumerate()
                    .filter(|(_, &s2)| tcrdist_gene_neighbor(s1, s2, threshold, ntrim, ctrim))
                    .map(|(jdx, _)| [idx, jdx])
                    .collect::<Vec<[usize; 2]>>()
            })
            .collect()
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .enumerate()
                .flat_map(|(idx, &s1)| {
                    seqs2
                        .iter()
                        .enumerate()
                        .filter(|(_, &s2)| tcrdist_gene_neighbor(s1, s2, threshold, ntrim, ctrim))
                        .map(|(jdx, _)| [idx, jdx])
                        .collect::<Vec<[usize; 2]>>()
                })
                .collect()
        })
    }
}

pub fn map_metric(metric_name: &str) -> fn(&[u8], &[u8]) -> u32 {
    let metric: Result<fn(&[u8], &[u8]) -> u32, &'static str> = match metric_name {
        "hamming" => Ok(hamming),
        "levenshtein" => Ok(levenshtein),
        "levenshtein_exp" => Ok(levenshtein_exp),
        _ => Err("The given metric is not an acceptable option. Try hamming, levenshtein, or levenshtein_exp.")
    };
    metric.unwrap()
}

pub fn str_cmp_matrix(seqs: &[&str], parallel: bool, metric: &str) -> Vec<u32> {
    let metric_fn = map_metric(metric);
    if parallel == false {
        let seqs_len: usize = seqs.len();
        let num_combinations: usize = seqs_len * (seqs_len - 1) / 2;
        let mut dists: Vec<u32> = vec![0; num_combinations];
        let mut counter: usize = 0;

        for (i, &s1) in seqs.iter().enumerate() {
            for &s2 in seqs[i + 1..].iter() {
                dists[counter] = metric_fn(s1.as_bytes(), s2.as_bytes());
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(i, &s1)| {
                    seqs[i + 1..]
                        .iter()
                        .map(|&s2| metric_fn(s1.as_bytes(), s2.as_bytes()))
                        .collect::<Vec<u32>>()
                })
                .collect()
        })
    }
}

pub fn str_cmp_one_to_many(seq: &str, seqs: &[&str], parallel: bool, metric: &str) -> Vec<u32> {
    let metric_fn = map_metric(metric);
    let seq_bytes: &[u8] = seq.as_bytes();
    if parallel == false {
        seqs.iter()
            .map(|&s| metric_fn(seq_bytes, s.as_bytes()))
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .map(|&s| metric_fn(seq_bytes, s.as_bytes()))
                .collect::<Vec<u32>>()
        })
    }
}

pub fn str_cmp_many_to_many(
    seqs1: &[&str],
    seqs2: &[&str],
    parallel: bool,
    metric: &str,
) -> Vec<u32> {
    let metric_fn = map_metric(metric);
    if parallel == false {
        let mut dists: Vec<u32> = vec![0; seqs1.len() * seqs2.len()];
        let mut counter: usize = 0;

        for &s1 in seqs1.iter() {
            for &s2 in seqs2.iter() {
                dists[counter] = metric_fn(s1.as_bytes(), s2.as_bytes());
                counter += 1;
            }
        }

        dists
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .flat_map(|&s1| {
                    seqs2
                        .iter()
                        .map(|&s2| metric_fn(s1.as_bytes(), s2.as_bytes()))
                        .collect::<Vec<u32>>()
                })
                .collect()
        })
    }
}

pub fn str_neighbor_matrix(
    seqs: &[&str],
    threshold: u32,
    parallel: bool,
    metric: &str,
) -> Vec<[usize; 2]> {
    let metric_fn = map_metric(metric);
    if parallel == false {
        seqs.iter()
            .enumerate()
            .flat_map(|(idx, &s1)| {
                seqs[idx + 1..]
                    .iter()
                    .enumerate()
                    .filter(move |(_, &s2)| metric_fn(s1.as_bytes(), s2.as_bytes()) <= threshold)
                    .map(move |(jdx, _)| [idx, idx + 1 + jdx])
            })
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .flat_map(|(idx, &s1)| {
                    seqs[idx + 1..]
                        .iter()
                        .enumerate()
                        .filter(move |(_, &s2)| {
                            metric_fn(s1.as_bytes(), s2.as_bytes()) <= threshold
                        })
                        .map(move |(jdx, _)| [idx, idx + 1 + jdx])
                        .collect::<Vec<[usize; 2]>>()
                })
                .collect()
        })
    }
}

pub fn str_neighbor_one_to_many(
    seq: &str,
    seqs: &[&str],
    threshold: u32,
    parallel: bool,
    metric: &str,
) -> Vec<usize> {
    let metric_fn = map_metric(metric);
    let seq_bytes = seq.as_bytes();
    if parallel == false {
        seqs.iter()
            .enumerate()
            .filter(|(_, &s)| metric_fn(seq_bytes, s.as_bytes()) <= threshold)
            .map(|(idx, _)| idx)
            .collect()
    } else {
        POOL.install(|| {
            seqs.par_iter()
                .enumerate()
                .filter(|(_, &s)| metric_fn(seq_bytes, s.as_bytes()) <= threshold)
                .map(|(idx, _)| idx)
                .collect::<Vec<usize>>()
        })
    }
}

pub fn str_neighbor_many_to_many(
    seqs1: &[&str],
    seqs2: &[&str],
    threshold: u32,
    parallel: bool,
    metric: &str,
) -> Vec<[usize; 2]> {
    let metric_fn = map_metric(metric);
    if parallel == false {
        seqs1
            .iter()
            .enumerate()
            .flat_map(|(idx, &s1)| {
                seqs2
                    .iter()
                    .enumerate()
                    .filter(|(_, &s2)| metric_fn(s1.as_bytes(), s2.as_bytes()) <= threshold)
                    .map(|(jdx, _)| [idx, jdx])
                    .collect::<Vec<[usize; 2]>>()
            })
            .collect()
    } else {
        POOL.install(|| {
            seqs1
                .par_iter()
                .enumerate()
                .flat_map(|(idx, &s1)| {
                    seqs2
                        .iter()
                        .enumerate()
                        .filter(|(_, &s2)| metric_fn(s1.as_bytes(), s2.as_bytes()) <= threshold)
                        .map(|(jdx, _)| [idx, jdx])
                        .collect::<Vec<[usize; 2]>>()
                })
                .collect()
        })
    }
}

pub fn str_bin_many_to_many(
    seqs1: &[&str],
    seqs2: &[&str],
    parallel: bool,
    metric: &str,
) -> Vec<u32> {
    let metric_fn = map_metric(metric);
    if parallel == false {
        let counter = seqs1.iter().fold(HashMap::new(), |mut acc, &s1| {
            seqs2.iter().for_each(|&s2| {
                *acc.entry(metric_fn(s1.as_bytes(), s2.as_bytes()))
                    .or_insert(0) += 1;
            });
            acc
        });
        let max_distance: usize = *counter
            .iter()
            .max_by(|a, b| a.0.cmp(&b.0))
            .map(|(k, _v)| k)
            .unwrap() as usize;

        let mut bincount = vec![0; max_distance + 1];
        for (key, value) in counter.iter() {
            bincount[*key as usize] = *value;
        }
        bincount
    } else {
        let counter = POOL.install(|| {
            seqs1
                .par_iter()
                .fold(
                    || HashMap::new(),
                    |mut acc, &s1| {
                        for s2 in seqs2.iter() {
                            *acc.entry(metric_fn(s1.as_bytes(), s2.as_bytes()))
                                .or_insert(0) += 1;
                        }
                        acc
                    },
                )
                .reduce(
                    || HashMap::new(),
                    |m1, m2| {
                        m2.iter().fold(m1, |mut acc, (k, vs)| {
                            *acc.entry(k.clone()).or_insert(0) += vs;
                            acc
                        })
                    },
                )
        });
        let max_distance: usize = *counter
            .iter()
            .max_by(|a, b| a.0.cmp(&b.0))
            .map(|(k, _v)| k)
            .unwrap() as usize;

        let mut bincount = vec![0; max_distance + 1];
        for (key, value) in counter.iter() {
            bincount[*key as usize] = *value;
        }
        bincount
    }
}
