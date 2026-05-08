# Compute BLUP Weights for One LDGM Block

Initial single-block R port of GraphLD's BLUP kernel. It computes
`sqrt(sample_size) * sigmasq * solve(P + sample_size * sigmasq * I, P %*% z)`.

## Usage

``` r
ldgm_blup_block(precision, z, sample_size, sigmasq)
```

## Arguments

- precision:

  Sparse LDGM precision matrix `P`.

- z:

  Numeric vector of Z scores for the block.

- sample_size:

  GWAS sample size.

- sigmasq:

  Per-variant heritability variance.

## Value

Numeric vector of BLUP weights.
