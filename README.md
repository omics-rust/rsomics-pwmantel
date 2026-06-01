# rsomics-pwmantel

Pairwise Mantel test across N distance matrices.

For every pair `(i < j)` of the input distance matrices, runs a Mantel test —
the Pearson (default) or Spearman correlation between the upper triangles of the
two matrices, with a permutation test for significance — and collates one row
per pair into a results table. Drop-in compatible with
`skbio.stats.distance.pwmantel`.

```
rsomics-pwmantel dm1.tsv dm2.tsv [dm3.tsv ...] \
    [--method pearson|spearman] [--permutations 999] \
    [--alternative two-sided|greater|less] [--seed S] [-o results.tsv]
```

Each input is an lsmat-format distance matrix (a blank top-left corner, a
tab-separated id header, then one `id<TAB>values…` row per sample). The matrices
must share ids; each pair's second matrix is reordered onto the first's id order
before correlating, so ids need not be in the same order. Matrices must be at
least 3×3.

The output is a TSV with one row per pair and columns `dm1`, `dm2`,
`statistic`, `p-value`, `n`, `method`, `permutations`, `alternative`. Pairs are
emitted in input order (`itertools.combinations` semantics: first matrix outer,
second inner). The `dm1`/`dm2` labels are the input filename stems.

## Statistic vs p-value

The **correlation statistic** is deterministic and reproduces scikit-bio's
`pwmantel()` to floating-point tolerance for both methods (tested value-exact to
1e-9 against a committed skbio-captured golden and, where the venv is present, a
live skbio oracle; ~5e-13 on a 400×400 fixture). Spearman is Pearson on the
average-ranked distances.

The **p-value** is a permutation Monte-Carlo estimate: the rows and columns of
the first matrix in each pair are permuted `--permutations` times and the
proportion of permuted statistics at least as extreme as the observed one (with
the `(count+1)/(perms+1)` correction) is reported. The permutations come from
this crate's own seeded RNG (SplitMix64 + Lemire-bounded Fisher-Yates),
reproducible across runs and thread counts for a given `--seed`, but **not** a
bit-for-bit reproduction of numpy's PCG64 stream — so the p-value is an estimate
that converges to skbio's as permutations grow, not an identical draw. Each pair
gets an independent draw by folding its index into the seed, mirroring how skbio
threads one rng through the pairwise loop.

## Origin

This crate is an independent Rust reimplementation of
`skbio.stats.distance.pwmantel`, informed by its BSD-3-licensed source (the
`itertools.combinations` pairing over the labelled matrices, the per-pair
reorder onto the first matrix's ids, the inline standardize-and-dot Mantel
kernel, the upper-triangle condensed form, and the `(count+1)/(perms+1)`
p-value) and by the method's primary reference:

- Mantel, N. (1967). "The detection of disease clustering and a generalized
  regression approach." *Cancer Research* 27(2): 209–220. PMID: 6018555.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio <https://scikit-bio.org> (BSD-3-Clause).
