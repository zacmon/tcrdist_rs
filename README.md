# tcrdist_rs
This is a Rust implementation of TCRdist. Previous versions have been written in [python/numba](https://tcrdist3.readthedocs.io/en/latest/index.html) and [C++](https://github.com/phbradley/conga/tree/master/tcrdist_cpp). This package also includes Python bindings to the [triple_accel](https://github.com/Daniel-Liu-c0deb0t/triple_accel) package for fast, SIMD Hamming and Levenshtein distances.

Each distance function has the following functions. We take `tcrdist_gene` as an example.

Distances can be computed with
- `tcrdist_gene(seq1, seq2, ...)`
- `tcrdist_gene_matrix(seqs, ...)` (computes the upper right triangle of the distance matrix among sequences)
- `tcrdist_gene_one_to_many(seq, seqs, ...)`
- `tcrdist_gene_many_to_many(seqs1, seqs2, ...)`
- `tcrdist_gene_pairwise(seqs1, seqs2, ...)`
  
where `...` implies other parameters, see the docstrings (e.g., `help(tcrdist_gene)`).
The above functions yield the distances as a list.

Neighbors can be found with
- `tcrdist_gene_neighbor(seqs1, seqs2, threshold, ...)`
- `tcrdist_gene_neighbor_matrix(seqs1, seqs2, threshold, ...)`
- `tcrdist_gene_neighbor_one_to_many(seqs1, seqs2, threshold, ...)`
- `tcrdist_gene_neighbor_many_to_many(seqs1, seqs2, threshold, ...)`
- `tcrdist_gene_neighbor_pairwise(seqs1, seqs2, threshold, ...)`
  
Whereas `tcrdist_gene_neighbor` returns a bool, all the other functions will return a list of lists, with each list containing the indices of the neighbors and the distance between the neighbors.

All the distances currently supported:

- `hamming`
- `levenshtein`
- `levenshtein_exp` (possibly faster Levenshtein distance which uses an exponential search)
- `tcrdist`
- `tcrdist_allele` (tcrdist computed using amino acid CDR3-V allele pairs)
- `tcrdist_gene` (tcrdist computed using amino acid CDR3-V gene pairs)

Presently, the non-TCRdist functions also have the following functions:
- `hamming_bin_many_to_many(seqs1, seqs2, parallel)`
- `levenshtein_bin_many_to_many(seqs1, seqs2, parallel)`
- `levenshtein_exp_bin_many_to_many(seqs1, seqs2, parallel)`
  
which compute all pairwise distances between `seqs1` and `seqs2` and bin them.
These are useful for characterizing distributions of [coincidence](https://www.pnas.org/doi/full/10.1073/pnas.2213264120) across large datasets.

## Installation
Given an environment with Python >= 3.7, this package can be installed with
```bash
pip install tcrdist_rs
```

## Development environment setup
Install Rust if necessary.
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Download the `maturin` library, assuming you are in an environment which has Python >= 3.7.

```bash
pip install maturin
```

To compile tcrdist, run

```bash
maturin develop --release --profile release
```

in the cloned repository.

## Example

We use the same [settings](https://tcrdist3.readthedocs.io/en/latest/tcrdistances.html#i-want-complete-control) as Python TCRdist. Notably, the keyword `parallel` enables the use of all CPU cores when computing many distances.
However, be mindful of enabling it and use it only when there are enough sequences such that the overhead of parallelizing is less costly than the actual computation.

```python
import tcrdist_rs

tcrs = [['CASRTGTVYEQYF', 'TRBV2*01'], ['CASSTLDRVYNSPLHF', 'TRBV6-2*01'],
        ['CASSESGGQVDTQYF', 'TRBV6-4*01'], ['CASSPTGPTDTQYF', 'TRBV18*01'],
        ['CASSYPIEGGRAFTGELFF', 'TRBV6-5*01']]

phmc_weight = 1
cdr1_weight = 1
cdr2_weight = 1
cdr3_weight = 3
gap_penalty = 4
ntrim = 3
ctrim = 2
fixed_gappos = False
parallel = False
distances = tcrdist_rs.tcrdist_allele_matrix(tcrs, phmc_weight, cdr1_weight,
                                             cdr2_weight, cdr3_weight,
                                             gap_penalty, ntrim, ctrim, fixed_gappos,
                                             parallel)
```

TCRdist can also be computed at the level of genes.

```python
import tcrdist_rs


tcrs = [['CASRTGTVYEQYF', 'TRBV2'], ['CASSTLDRVYNSPLHF', 'TRBV6-2'],
        ['CASSESGGQVDTQYF', 'TRBV6-4'], ['CASSPTGPTDTQYF', 'TRBV18'],
        ['CASSYPIEGGRAFTGELFF', 'TRBV6-5']]


distances = tcrdist_rs.tcrdist_gene_matrix(tcrs, parallel=False)
```

## Outlook
- Vectorize/simdify computations for amino acid lookup?
- Improve lookup table performance (hash map vs. match vs. array)
- Precompute V allele lookup tables from custom databases as opposed to hard code (or memoize V allele dists as they are used in the input sequences)

## References
- Dash et al. (2017) "Quantifiable predictive features define epitope-specific T cell receptor repertoires." Nature 547(89-93). https://doi.org/10.1038/nature22383
- Mayer-Blackwell et al. (2021) "TCR meta-clonotypes for biomarker discovery with tcrdist3 enabled identification of public, HLA-restricted clusters of SARS-CoV-2 TCRs." eLife 10:68605. https://doi.org/10.7554/eLife.68605
