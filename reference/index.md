# Package index

## LDGM construction

- [`ldgm_tree_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables.md)
  : Bundle Canonical Tree-Diff Tables for LDGM Construction
- [`ldgm_tree_tables_from_tskit()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables_from_tskit.md)
  : Build LDGM Tree Tables from tskit Inputs
- [`ldgm_brick_ts()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_ts.md)
  : Brick an LDGM Tree-Table Bundle
- [`ldgm_make_ldgm()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_ldgm.md)
  : Make an LDGM from a Tree-Table Bundle
- [`ldgm_make_ldgm_from_tree_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_ldgm_from_tree_tables.md)
  : Make an LDGM from Canonical Tree-Diff Tables
- [`ldgm_make_ldgm_from_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_ldgm_from_tables.md)
  : Make an LDGM from Canonical Bricked Tables
- [`ldgm_brick_edges_from_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_edges_from_tables.md)
  : Brick Tree-Sequence Edge Tables from Canonical Tree-Diff Tables
- [`ldgm_brick_graph_inputs_from_edges()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_graph_inputs_from_edges.md)
  : Build Brick-Haplotype Input Tables from Bricked Edges
- [`ldgm_mutations_to_bricks()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_mutations_to_bricks.md)
  : Map Mutations to Bricks from Canonical Tables
- [`ldgm_brick_haplo_graph()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_haplo_graph.md)
  : Construct a Brick-Haplotype Graph from Canonical Tables
- [`ldgm_make_snplist()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_snplist.md)
  : Make an LDGM SNP List from Canonical Tables

## Graph operations

- [`ldgm_edge_list()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_edge_list.md)
  : Validate and Normalize an LDGM Edge List
- [`ldgm_remove_node()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_remove_node.md)
  : Remove a Node from a Directed Weighted LDGM Graph
- [`ldgm_reduce_graph()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reduce_graph.md)
  : Reduce a Brick-Haplotype Graph to an LDGM SNP Graph
- [`ldgm_return_edgelist()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_return_edgelist.md)
  : Return an LDGM Edge List

## Precision operators

- [`ldgm_sparse_precision()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_sparse_precision.md)
  : Build a Sparse LDGM Precision Matrix from an Edge List
- [`ldgm_precision()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision.md)
  : Create an LDGM Precision Object
- [`ldgm_precision_select()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_select.md)
  : Select a Schur-Complement View of an LDGM Precision Object
- [`ldgm_precision_matrix()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_matrix.md)
  : Extract an LDGM Precision Matrix
- [`ldgm_precision_multiply()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_multiply.md)
  : Multiply by an LDGM Precision Matrix
- [`ldgm_precision_solve()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_solve.md)
  : Solve an LDGM Precision Linear System
- [`ldgm_precision_scale()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_scale.md)
  : Scale an LDGM Precision Matrix
- [`ldgm_precision_update()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_update.md)
  : Update an LDGM Precision Matrix Diagonal
- [`ldgm_precision_update_element()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_update_element.md)
  : Update One LDGM Precision Diagonal Element
- [`ldgm_precision_logdet()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_precision_logdet.md)
  : Log Determinant of an LDGM Precision Matrix
- [`ldgm_inverse_diagonal()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_inverse_diagonal.md)
  : Diagonal of the Inverse Precision Matrix
- [`ldgm_gaussian_likelihood()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_gaussian_likelihood.md)
  : Gaussian Log-Likelihood for GraphLD Precision-Premultiplied
  Statistics
- [`ldgm_gaussian_likelihood_gradient()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_gaussian_likelihood_gradient.md)
  : Gaussian Log-Likelihood Gradient
- [`ldgm_gaussian_likelihood_hessian()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_gaussian_likelihood_hessian.md)
  : Gaussian Log-Likelihood Hessian Approximation
- [`ldgm_variant_solve()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_variant_solve.md)
  : Solve Variant-Level Right-Hand Sides with Duplicate Precision
  Indices
- [`ldgm_variant_info()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_variant_info.md)
  : Extract LDGM Variant Information

## GraphLD workflows

- [`ldgm_load_ldgm()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_load_ldgm.md)
  : Load a GraphLD-Style LDGM Block
- [`ldgm_read_edgelist()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_read_edgelist.md)
  : Read an LDGM Edge List File
- [`ldgm_read_snplist()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_read_snplist.md)
  : Read an LDGM SNP List File
- [`ldgm_merge_alleles()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_merge_alleles.md)
  : Compare Alleles and Return Phase
- [`ldgm_merge_snplists()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_merge_snplists.md)
  : Merge an LDGM Precision Object with Summary Statistics
- [`ldgm_partition_variants()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_partition_variants.md)
  : Partition Variants by LDGM Metadata Blocks
- [`ldgm_blup_block()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_blup_block.md)
  : Compute BLUP Weights for One LDGM Block
- [`ldgm_run_blup()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_run_blup.md)
  : Compute BLUP Weights Across LDGM Blocks
- [`ldgm_run_clump()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_run_clump.md)
  : Perform LD Clumping Across LDGM Blocks
- [`ldgm_reml_link()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reml_link.md)
  : Softplus Link Used by GraphREML Heritability Models
- [`ldgm_reml_block()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reml_block.md)
  : Compute One GraphREML Block Likelihood Slice
- [`ldgm_run_reml()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_run_reml.md)
  : Run a Serial GraphREML-Style Optimizer

## Score-test HDF5

- [`ldgm_hdf5_filter_info()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_hdf5_filter_info.md)
  : Inspect HDF5 Compression Filters
- [`ldgm_write_score_test_hdf5()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_write_score_test_hdf5.md)
  : Write GraphLD-Style Score-Test HDF5 Output
- [`ldgm_read_score_test_hdf5()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_read_score_test_hdf5.md)
  : Read GraphLD-Style Score-Test HDF5 Output

## Native runtime

- [`ldgm_openmp_info()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_openmp_info.md)
  : Inspect Native OpenMP Support
- [`ldgm_set_openmp_threads()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_set_openmp_threads.md)
  : Set Native OpenMP Threads
