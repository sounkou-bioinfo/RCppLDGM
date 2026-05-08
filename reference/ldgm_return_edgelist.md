# Return an LDGM Edge List

Normalizes an LDGM graph edge list and rounds weights using the upstream
Python `ldgm.return_edgelist()` convention of four decimal places by
default.

## Usage

``` r
ldgm_return_edgelist(graph, digits = 4L)
```

## Arguments

- graph:

  A data frame with columns `from`, `to`, and `weight`.

- digits:

  Number of decimal places used to round edge weights.

## Value

A data frame with columns `from`, `to`, and `weight`.
