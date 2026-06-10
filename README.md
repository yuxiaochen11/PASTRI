# PASTRI

PASTRI infers a directed cell-state transition matrix from a reconstructed
lineage tree and terminal-cell annotations. It is designed for lineage trees
derived from CRISPR/Cas9 barcodes, somatic mutations, or similar heritable
marks.

The method does not infer transitions from expression pseudotime. Instead, it
uses the topology of the lineage tree: terminal cells whose latest common
ancestor (LCA) is close to the leaves are treated as observations from a short
evolutionary separation, while pairs sharing an older LCA represent a longer
separation. PASTRI combines these pairwise relationships across depth values,
then applies user-defined upper bounds to encode allowed or forbidden
directions.

## Algorithm overview

PASTRI runs in four stages.

### 1. Convert the lineage tree into pairwise LCA features

`calculate_lca_depths()` reads a Newick tree and, for every internal node,
counts:

- `distance_to_root`: number of edges from the node to the root;
- `distance_to_farthest_leaf`: maximum number of descendant edges to a leaf.

Let \(D_{\max}\) be the largest `distance_to_farthest_leaf` in the tree. For
each internal node, the implementation calculates

\[
\mathrm{normalized\_depth}
= \frac{1}{2}\left(
\frac{d_{\mathrm{root}}}{D_{\max}}
+ 1
- \frac{d_{\mathrm{leaf}}}{D_{\max}}
\right)
\]

and

\[
\mathrm{lca\_normalized\_height}
= 1 - \mathrm{normalized\_depth}.
\]

Every unordered pair of terminal nodes is then assigned the metrics of its
LCA. The result contains `node1`, `node2`, `subtree_root`, `lca_depth`,
`lca_depth_norm`, and `lca_normalized_height`.

### 2. Estimate two-cell correlations at each depth

Terminal nodes are matched to `cell_info$nodeLabel` and assigned cell types.
For each selected depth value \(u\), PASTRI counts unordered cell-type pairs
and constructs a symmetric two-cell correlation matrix \(C(u)\). Cell-type
frequencies among the selected terminal cells are used to standardize this
matrix.

Only depth values satisfying `u <= fi_depth` are retained. Therefore,
`fi_depth` is a depth cutoff, not an optimization tolerance.

### 3. Recover one transition matrix per depth

For each retained \(C(u)\), PASTRI:

1. standardizes the correlation matrix by cell-type frequencies;
2. performs an eigendecomposition;
3. truncates negative eigenvalues to zero;
4. rescales eigenvalues by the depth-dependent power `1 / (2 * u)`;
5. removes negative or undefined entries;
6. normalizes every column to sum to one.

In all returned matrices, **columns are start cell types and rows are end cell
types**. Thus `matrix["C5", "C4"]` is the inferred transition strength from
`C4` to `C5`.

### 4. Apply constraints and obtain the consensus matrix

The mean of the depth-specific matrices initializes an `L-BFGS-B`
optimization. PASTRI minimizes the sum of Frobenius distances between the
candidate matrix \(T\) and all depth-specific matrices:

\[
\min_T \sum_u \lVert T - T(u) \rVert_F.
\]

`Bound_Matrix` supplies the element-wise upper bounds:

- `0` forbids a transition;
- `Inf` leaves a transition unconstrained above;
- a finite positive value sets a numerical upper bound.

The bound matrix must have the same dimensions and cell-type order as
`cell_type_list`. After optimization, negative/undefined entries are removed
and columns are normalized again.

## Input requirements

### Lineage tree

- Newick format;
- terminal labels must be unique;
- terminal labels must match `cell_info$nodeLabel`.

### Cell metadata

`cell_info` must contain:

| Column | Meaning |
| --- | --- |
| `nodeLabel` | label matching a terminal node in the Newick tree |
| `celltype` | cell-state annotation |
| `cellNum` | number of cells represented by the node; only `cellNum == 1` is used |

Each value in `cell_type_list` should occur among the retained terminal cells.

## Installation

```r
install.packages("remotes")
remotes::install_github("yuxiaochen11/PASTRI")
```

For a reproducible Conda/R/Jupyter setup and solutions to common installation
errors, see [Installation and troubleshooting](articles/troubleshooting.html).

## Quick start

```r
library(PASTRI)

# 1. Load the example lineage tree.
tree_file <- system.file("extdata", "HSC.nwk", package = "PASTRI")

# 2. Calculate pairwise LCA features.
node_pair_depth <- calculate_lca_depths(
  tree.path = tree_file,
  mc.cores = 1
)

# 3. Load terminal-cell annotations.
cell_file <- system.file("extdata", "HSC_nodeInfos.csv", package = "PASTRI")
cell_info <- read.csv(cell_file, stringsAsFactors = FALSE)

# The order defines the row/column order of every matrix.
cell_types <- c("C4", "C5", "C6", "C7", "C8")

# 4. Encode an irreversible progression C4 -> C5 -> C6 -> C7 -> C8.
# Rows are destinations and columns are sources. Entries above the diagonal
# represent backward transitions and are forbidden here.
bound_matrix <- matrix(
  Inf,
  nrow = length(cell_types),
  ncol = length(cell_types),
  dimnames = list(cell_types, cell_types)
)
bound_matrix[upper.tri(bound_matrix)] <- 0

# 5. Infer the constrained transition matrix.
result <- get_optimal_transition_matrix(
  node_pair_depth = node_pair_depth,
  cell_info = cell_info,
  Sel_u = "lca_normalized_height",
  fi_depth = 0.2,
  Bound_Matrix = bound_matrix,
  cell_type_list = cell_types,
  mc.cores = 1
)

# Column = start state; row = end state.
result$optimal_matrix

# Long-form representation of the same matrix.
head(result$optimal_norm_df_dataframe)
```

## Main functions

| Function | Purpose |
| --- | --- |
| `calculate_lca_depths()` | calculate pairwise LCA depth/height features from a Newick tree |
| `calculate_transition_matrix()` | build depth-specific matrices and optimize their constrained consensus |
| `get_optimal_transition_matrix()` | recommended high-level interface for filtering cells and formatting results |
| `calculate_transition_matrix_error()` | Frobenius-distance objective used by the optimizer |

## Output

`get_optimal_transition_matrix()` returns:

- `optimal_matrix`: constrained transition matrix, with destinations in rows
  and sources in columns;
- `optimal_norm_df_dataframe`: long-form table containing `start_cell`,
  `end_cell`, `norm_optimal_T`, and the evaluated depth label.

## Practical notes

- On Windows, PASTRI automatically uses one core because
  `parallel::mclapply()` does not support fork-based multicore execution.
- The number of terminal-node pairs grows as \(n(n-1)/2\); large trees can
  require substantial memory and runtime.
- At least one observed depth must satisfy `u <= fi_depth`. If none does,
  increase `fi_depth` or inspect the selected depth column.
- `cell_type_list` order and `Bound_Matrix` order must be identical.
- The output represents model-derived transition strengths under the supplied
  lineage and constraints; causal interpretation requires independent
  biological evidence.
