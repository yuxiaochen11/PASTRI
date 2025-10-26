#' @title Get Optimal Transition Matrix
#'
#' @description
#' Compute the optimal cell transition matrix by integrating node-pair depth
#' information with per-cell annotations. The function:
#' (1) filters to terminal (single-leaf) cells of selected cell types,
#' (2) calls \code{calculate_transition_matrix()} to infer the best
#'     transition matrix under structural constraints,
#' (3) flattens that matrix into a long data frame of start/end cell-type
#'     pairs with normalized transition strengths, and
#' (4) labels results with the evaluated depth cutoff.
#'
#' This version is Windows-safe: all internal parallelism is automatically
#' coerced to single-core on Windows so that vignette building (pkgdown)
#' will not fail with \code{'mc.cores' > 1 is not supported on Windows}.
#'
#' @param node_pair_depth A data frame containing node pair depth information
#'   (e.g. from \code{calculate_lca_depths()}).
#' @param cell_info A data frame with per-cell info; must contain at least
#'   \code{cellNum} and \code{celltype}.
#'   Only rows with \code{cellNum == 1} are retained (terminal leaves),
#'   and then filtered to the provided \code{cell_type_list}.
#' @param Sel_u Character string choosing which depth/height measure to use
#'   from \code{node_pair_depth} (e.g. "lca_normalized_height").
#' @param fi_depth Numeric. Depth cutoff / threshold passed to
#'   \code{calculate_transition_matrix()}.
#' @param Bound_Matrix Numeric matrix encoding constraints (e.g. irreversible
#'   transitions). Diagonal/upper-triangular structure often encodes prior
#'   directionality.
#' @param cell_type_list Character vector of cell types of interest, in the
#'   order that defines the transition structure.
#' @param mc.cores Integer. Desired parallel cores. On Windows this will be
#'   silently coerced to 1. On non-Windows it will be coerced to at least 1.
#'
#' @return A list with:
#' \itemize{
#'   \item \strong{optimal_matrix}: the inferred optimal transition matrix
#'         (cell type x cell type).
#'   \item \strong{optimal_norm_df_dataframe}: a long-form data frame with
#'         columns \code{type_combined_new}, \code{start_cell}, \code{end_cell},
#'         \code{norm_optimal_T}, and \code{depth} ("d_<fi_depth>").
#' }
#'
#' @examples
#' \dontrun{
#' res <- get_optimal_transition_matrix(
#'   node_pair_depth   = node_pair_depth,
#'   cell_info         = cell_info,
#'   Sel_u             = "lca_normalized_height",
#'   fi_depth          = 0.2,
#'   Bound_Matrix      = Bound_Matrix,
#'   cell_type_list    = c("C4","C5","C6","C7","C8"),
#'   mc.cores          = 8
#' )
#' str(res$optimal_matrix)
#' head(res$optimal_norm_df_dataframe)
#' }
#'
#' @export
get_optimal_transition_matrix <- function(
    node_pair_depth,
    cell_info,
    Sel_u,
    fi_depth,
    Bound_Matrix,
    cell_type_list,
    mc.cores = parallel::detectCores()
) {
  ## ----------------------------------------------------------------------
  ## 0. make mc.cores safe on Windows
  ##    (pkgdown builds vignette on Windows-like env where mclapply>1 dies)
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1L
  } else {
    mc.cores <- max(1L, as.integer(mc.cores))
  }

  ## helper to safely do parallel::mclapply or fallback to lapply
  safe_mclapply <- function(X, FUN, ...) {
    if (.Platform$OS.type == "windows" || mc.cores <= 1L) {
      out <- lapply(X, FUN, ...)
    } else {
      out <- parallel::mclapply(X, FUN, mc.cores = mc.cores, ...)
    }
    out
  }

  ## ----------------------------------------------------------------------
  ## 1. subset to terminal / single-leaf cells of desired types
  single_cell_info <- dplyr::filter(cell_info, .data$cellNum == 1)
  terminal_cells   <- dplyr::filter(single_cell_info,
                                    .data$celltype %in% cell_type_list)

  ## ----------------------------------------------------------------------
  ## 2. compute optimal transition matrix under constraints
  ## NOTE: we forward mc.cores (as 'use_cores') to your calculate_transition_matrix
  ## so that downstream also respects Windows safety.
  optimal_trans_matrix_obj <- calculate_transition_matrix(
    node_pair_depth   = node_pair_depth,
    node_info         = terminal_cells,
    sel_u             = Sel_u,
    use_cores         = mc.cores,      # <-- was hard-coded 80
    fi_depth          = fi_depth,
    bound_matrix      = Bound_Matrix,
    cell_type_list    = cell_type_list
  )

  optimal_matrix <- optimal_trans_matrix_obj$optimal_trans_matrix

  ## ----------------------------------------------------------------------
  ## 3. flatten matrix into long-form df
  combinations <- expand.grid(
    start_cell = colnames(optimal_matrix),
    end_cell   = rownames(optimal_matrix),
    stringsAsFactors = FALSE
  )

  optimal_norm_list <- safe_mclapply(
    X = seq_len(nrow(combinations)),
    FUN = function(i) {
      start_cell_i <- combinations$start_cell[i]
      end_cell_i   <- combinations$end_cell[i]

      value_i <- optimal_matrix[end_cell_i, start_cell_i]

      dplyr::tibble(
        type_combined_new = paste(start_cell_i, end_cell_i, sep = "_"),
        start_cell        = start_cell_i,
        end_cell          = end_cell_i,
        norm_optimal_T    = value_i
      )
    }
  )

  optimal_norm_df <- dplyr::bind_rows(optimal_norm_list)

  ## annotate depth label
  optimal_norm_df <- optimal_norm_df %>%
    dplyr::mutate(depth = paste0("d_", fi_depth))

  ## ----------------------------------------------------------------------
  ## 4. return structured result
  return(list(
    optimal_matrix             = optimal_matrix,
    optimal_norm_df_dataframe  = optimal_norm_df
  ))
}
