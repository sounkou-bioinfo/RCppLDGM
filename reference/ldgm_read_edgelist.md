# Read an LDGM Edge List File

Reads a comma-separated, headerless LDGM edge-list file. GraphLD/LDGM
files store node ids as zero-based integers; this reader converts them
to ordinary one-based R ids in the returned data frame.

## Usage

``` r
ldgm_read_edgelist(path)
```

## Arguments

- path:

  Path to the edge-list file.

## Value

A validated edge-list data frame with one-based `from`/`to` ids.
