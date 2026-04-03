# PASTRI

PASTRI infers the transition probability matrix between cell states along a cell
lineage tree (e.g. CRISPR lineage tracing / somatic mutation phylogeny).  
It uses LCA depth between pairs of terminal cells, applies biologically
irreversible constraints, and solves for an optimal directed transition matrix
across specified cell types.

## Installation

```r
# install.packages("devtools")
devtools::install_github("yuxiaochen11/PASTRI")
```


## Quick Start

```
library(PASTRI)

# 1. Load example lineage tree shipped with the package
tree_file <- system.file("extdata", "HSC.nwk", package = "PASTRI")

# 2. Compute pairwise LCA depth / height features
node_depth_df <- calculate_lca_depths(tree.path = tree_file)

- If you want to manually specify the depth of an internal node (i.e., the root of a subtree), you can add this information to the node_depth_df data frame.

# 3. Load example cell metadata
cell_file <- system.file("extdata", "HSC_nodeInfos.csv", package = "PASTRI")
cell_info  <- read.csv(cell_file, stringsAsFactors = FALSE)

# 4. Infer optimal constrained transition matrix
res <- get_optimal_transition_matrix(
  node_pair_depth = node_depth_df,
  cell_info       = cell_info,
  Sel_u           = "lca_normalized_height",
  fi_depth        = 0.2,
  Bound_Matrix    = matrix(c(Inf,0,0,0,0,
                             Inf,Inf,0,0,0,
                             Inf,Inf,Inf,0,0,
                             Inf,Inf,Inf,Inf,0,
                             Inf,Inf,Inf,Inf,Inf), nrow=5, byrow=TRUE),
  cell_type_list  = c("C4","C5","C6","C7","C8"),
  mc.cores        = 1
)

res$optimal_matrix
```
