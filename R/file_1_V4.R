#' @title Calculate the Depth of Latest Common Ancestor (LCA) for Cell Lineage Tree Nodes.
#'
#' @description
#' This function calculates the depth of the latest common ancestor (LCA) for pairs of leaf nodes in a cell lineage tree.
#' It computes the distance from each node to the root (distance_to_root), the distance to the farthest leaf (distance_to_farthest_leaf),
#' and normalized versions of these quantities. It then enumerates all pairs of leaf nodes, finds their LCA,
#' and reports the depth/height metrics of that LCA.
#'
#' @param tree.path The input path of the cell lineage tree in Newick format
#'   (for example, a path obtained by
#'   \code{system.file("extdata", "HSC.nwk", package = "PASTRI")}).
#' @param mc.cores Integer. Desired number of parallel workers for internal
#'   computations. On Windows, values greater than 1 are not supported by
#'   \code{parallel::mclapply} and will be silently coerced to 1.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \strong{node1}: The label of the first leaf node in the pair.
#'   \item \strong{node2}: The label of the second leaf node in the pair.
#'   \item \strong{subtree_root}: The internal node ID corresponding to the latest common ancestor (LCA).
#'   \item \strong{lca_depth}: The distance from the LCA to the farthest leaf in that subtree.
#'   \item \strong{lca_depth_norm}: The normalized distance of the LCA to the farthest leaf.
#'   \item \strong{lca_normalized_height}: The normalized height (1 - normalized depth) of the LCA.
#' }
#'
#' @examples
#' \dontrun{
#' tree_file <- system.file("extdata", "HSC.nwk", package = "PASTRI")
#' res <- calculate_lca_depths(tree.path = tree_file, mc.cores = 2)
#' head(res)
#' }
#'
#' @export
calculate_lca_depths <- function(tree.path, mc.cores = parallel::detectCores()) {
  
  # --- make mc.cores safe on Windows ---------------------------------------
  # parallel::mclapply() uses fork and does not support mc.cores > 1 on Windows.
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1L
  } else {
    # just in case someone passed something weird like 0 or negative
    mc.cores <- max(1L, as.integer(mc.cores))
  }
  
  # --- read tree -----------------------------------------------------------
  tree <- read.tree(tree.path)
  
  # fortify() (from ggtree / ggplot2 ecosystem) turns 'phylo' into a data.frame
  tree_df <- fortify(tree)
  
  # all unique internal nodes that appear as 'parent'
  parent_nodes <- unique(tree_df$parent)
  
  # --- compute depth metrics for each internal node / subtree -------------
  subtree_depths <- parallel::mclapply(
    X = seq_along(parent_nodes),
    mc.cores = mc.cores,
    FUN = function(j) {
      
      current_parent <- parent_nodes[j]
      
      ## distance_to_root: climb parent pointers to the root
      node_index <- which(tree_df$node == current_parent)
      distance_to_root <- 0
      repeat {
        # if node is its own parent we are at the root
        if (tree_df$node[node_index] == tree_df$parent[node_index]) {
          break
        } else {
          parent_node <- tree_df$parent[node_index]
          node_index <- which(tree_df$node == parent_node)
          distance_to_root <- distance_to_root + 1
        }
      }
      
      ## distance_to_farthest_leaf: walk downward breadth-wise until no children
      daughter_nodes <- tree_df %>%
        dplyr::filter(parent %in% current_parent & (parent != node))
      distance_to_farthest_leaf <- 0
      repeat {
        if (nrow(daughter_nodes) == 0) break
        daughter_nodes <- tree_df %>%
          dplyr::filter((parent %in% daughter_nodes$node) & (parent != node))
        distance_to_farthest_leaf <- distance_to_farthest_leaf + 1
      }
      
      data.frame(
        subtree_root              = current_parent,
        distance_to_root          = distance_to_root,
        distance_to_farthest_leaf = distance_to_farthest_leaf
      )
    }
  ) %>% plyr::rbind.fill()
  
  # --- normalize depths ----------------------------------------------------
  max_far <- max(subtree_depths$distance_to_farthest_leaf)
  subtree_depths$max_distance_to_farthest_leaf <- max_far
  
  subtree_depths$distance_to_farthest_leaf_norm <-
    subtree_depths$distance_to_farthest_leaf / max_far
  
  subtree_depths$normalized_depth <- (
    (subtree_depths$distance_to_root / max_far) +
      (1 - subtree_depths$distance_to_farthest_leaf / max_far)
  ) / 2
  
  subtree_depths$lca_normalized_height <- 1 - subtree_depths$normalized_depth
  
  # --- enumerate all leaf node pairs --------------------------------------
  leaf_nodes <- tree$tip.label
  
  node_pairs_combinations <- combn(unique(leaf_nodes), 2)
  node_pairs_df <- data.frame(
    node1 = node_pairs_combinations[1, ],
    node2 = node_pairs_combinations[2, ],
    stringsAsFactors = FALSE
  )
  
  # --- for each leaf pair, compute LCA and extract depth metrics ----------
  node_pair_distances <- parallel::mclapply(
    X = seq(nrow(node_pairs_df)),
    mc.cores = mc.cores,
    FUN = function(n) {
      
      node1 <- node_pairs_df$node1[n]
      node2 <- node_pairs_df$node2[n]
      
      # latest common ancestor in 'phylo' tree
      latest_common_ancestor <- getMRCA(tree, c(node1, node2))
      
      # look up that ancestor's metrics
      lca_row <- dplyr::filter(subtree_depths,
                               subtree_root == latest_common_ancestor)
      
      data.frame(
        node1                  = as.character(node1),
        node2                  = as.character(node2),
        subtree_root           = latest_common_ancestor,
        lca_depth              = lca_row$distance_to_farthest_leaf,
        lca_depth_norm         = lca_row$distance_to_farthest_leaf_norm,
        lca_normalized_height  = lca_row$lca_normalized_height,
        stringsAsFactors       = FALSE
      )
    }
  ) %>% dplyr::bind_rows()
  
  return(node_pair_distances)
}
