# Compare Alleles and Return Phase

R port of GraphLD's `merge_alleles()` helper. Alleles are compared
case-insensitively.

## Usage

``` r
ldgm_merge_alleles(anc_alleles, deriv_alleles, ref_alleles, alt_alleles)
```

## Arguments

- anc_alleles:

  Ancestral alleles from the LDGM SNP list.

- deriv_alleles:

  Derived alleles from the LDGM SNP list.

- ref_alleles:

  Reference alleles from summary statistics.

- alt_alleles:

  Alternative alleles from summary statistics.

## Value

Numeric phase vector: `1` for exact matches, `-1` for swapped matches,
`0` for mismatches, and `NA_real_` when both summary-stat alleles are
missing.
