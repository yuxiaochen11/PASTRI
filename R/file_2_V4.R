#' @title Calculate Transition Matrix Error
#'
#' @description
#' This function computes the total error between an initial guess transition matrix and a list
#' of observed transition matrices. The error is defined as the cumulative Frobenius norm
#' difference between the guess and each observed matrix.
#'
#' @param initial_vector A numeric vector representing the flattened initial guess for the transition matrix.
#' @param observed_matrices A list of observed transition matrices.
#' @param matrix_size The size (number of rows/columns) of the transition matrix.
#' @param type_list A character vector of cell types used to set the row and column names of the matrix.
#'
#' @return A numeric value representing the total error.
#' @export
#'
#' @examples
#' # For example, if you have a flattened 5x5 initial guess vector and two 5x5 observed matrices:
#' initial_vec <- runif(25)  # generate a random initial guess for a 5x5 matrix
#' observed_list <- list(matrix(runif(25), nrow = 5), matrix(runif(25), nrow = 5))
#' cell_types <- c("C4", "C5", "C6", "C7", "C8")
#' error_value <- calculate_transition_matrix_error(initial_vec, observed_list, 5, cell_types)
calculate_transition_matrix_error <- function(initial_vector, observed_matrices, matrix_size, type_list) {
  # Convert the initial vector to a transition matrix guess
  transition_matrix_guess <- matrix(initial_vector, nrow = matrix_size, byrow = FALSE)
  colnames(transition_matrix_guess) <- type_list
  rownames(transition_matrix_guess) <- type_list

  total_error <- 0.0

  # Loop over each observed matrix and accumulate the Frobenius norm error
  for (i in seq_along(observed_matrices)) {
    obs_matrix <- as.matrix(observed_matrices[[i]])
    if (!is.matrix(obs_matrix)) {
      stop("Observed matrix at index ", i, " is not a matrix!")
    }
    total_error <- total_error + norm(transition_matrix_guess - obs_matrix, type = "F")
  }
  return(total_error)
}

#' @title Calculate Transition Matrix from Two-Cell Correlation
#'
#' @description
#' This function calculates cell transition matrices based on two-cell correlation functions C(u).
#' It first computes the frequency of each cell type (using only single-leaf nodes), then builds
#' a standardized cell type combination to avoid duplicate counts. Using the node pair depth data
#' and selected depth measure (e.g., "lca_normalized_height" or "lca_depth"), the function computes
#' the correlation frequencies for each depth. For each depth (up to a maximum threshold fi_depth),
#' a corresponding transition matrix is computed using eigenvalue decomposition, where the eigenvalues
#' are adjusted according to the depth. Finally, the average transition matrix from all depths is
#' used as the initial guess for an optimization procedure (via the L-BFGS-B method) to derive the
#' optimal transition matrix.
#'
#' @param node_pair_depth A data frame containing node pair depth information and related metrics
#' (e.g., normalized LCA height).
#' @param node_info A data frame containing detailed information for each node, including columns
#' 'celltype' (cell type) and 'nodeLabel' (node label).
#' @param sel_u A string indicating the depth column to be used (e.g., "lca_normalized_height" or "lca_depth").
#' @param use_cores An integer specifying the number of cores to use for parallel computation (default is 80).
#' @param fi_depth A numeric value indicating the maximum depth threshold for computing the transition matrix.
#' @param bound_matrix A constrained matrix
#' @param cell_type_list A list of selected cell types

#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{all_trans_matrix_list}: A list of transition matrices computed for different depths.
#'   \item \strong{optimal_trans_matrix}: The optimal transition matrix obtained from the optimization procedure.
#' }
#' @export
#'
#' @examples
#' # Suppose you have node pair depth data (node_pair_depth) and node information (node_info), and you select
#' # "lca_normalized_height" as the depth measure with fi_depth set to 0.2:
#' result <- calculate_transition_matrix(node_pair_depth, node_info, sel_u = "lca_normalized_height", use_cores = 80, fi_depth = 0.2, bound_matrix, cell_type_list)
#' # To view the optimal transition matrix:
#' print(result$optimal_trans_matrix)
calculate_transition_matrix <- function(node_pair_depth, node_info, sel_u, use_cores, fi_depth, bound_matrix, cell_type_list) {
  options(digits = 8)

  #### 1. Compute the frequency of each cell type (only for single-leaf nodes)
  node_info$celltype <- factor(node_info$celltype, levels = cell_type_list)
  freq_celltype <- node_info %>%
    filter(celltype %in% cell_type_list) %>%
    `[`("celltype") %>%
    table() %>%
    as.data.frame() %>%
    mutate(freq = Freq / sum(Freq))

  #### 2. Compute two-cell correlation functions C(u)
  # 2.1 Retain and round the selected depth measure in node_pair_depth
  node_pair_depth <- node_pair_depth %>%
    mutate(lca_normalized_height_raw = lca_normalized_height) %>%
    mutate(lca_normalized_height = round(lca_normalized_height, 8))

  # 2.2 Add cell type information for each node in the node pair data
  node_pair_depth$node1_type <- node_info$celltype[match(node_pair_depth$node1, node_info$nodeLabel)]
  node_pair_depth$node2_type <- node_info$celltype[match(node_pair_depth$node2, node_info$nodeLabel)]
  node_pair_depth <- node_pair_depth %>%
    filter(node1_type %in% cell_type_list & node2_type %in% cell_type_list)

  # 2.3 Construct standardized cell type combinations to avoid duplicate counts
  unique_type_combn <- t(combn(unique(cell_type_list), 2)) %>% as.data.frame(stringsAsFactors = FALSE)
  colnames(unique_type_combn) <- c("type_1", "type_2")
  unique_type_combn <- bind_rows(
    data.frame(type_1 = cell_type_list, type_2 = cell_type_list, stringsAsFactors = FALSE),
    unique_type_combn
  )
  unique_type_combn$ref_comb <- paste(unique_type_combn$type_1, unique_type_combn$type_2, sep = "_")
  unique_type_combn$reverse_comb <- paste(unique_type_combn$type_2, unique_type_combn$type_1, sep = "_")

  combn_standard_df <- data.frame(
    ref_comb = c(unique_type_combn$ref_comb, unique_type_combn$ref_comb),
    present_comb = c(unique_type_combn$ref_comb, unique_type_combn$reverse_comb),
    stringsAsFactors = FALSE
  )
  rm(unique_type_combn)

  # 2.4 Standardize the cell type combination in node_pair_depth
  node_pair_depth$type_combined <- paste(node_pair_depth$node1_type, node_pair_depth$node2_type, sep = "_")
  node_pair_depth$type_combined <- combn_standard_df$ref_comb[match(node_pair_depth$type_combined, combn_standard_df$present_comb)]
  rm(combn_standard_df)

  # 2.5 Create a cell type frequency matrix
  freq_celltype_matrix <- matrix(unlist(freq_celltype$freq),
                                 nrow = nrow(freq_celltype), ncol = nrow(freq_celltype), byrow = TRUE)
  colnames(freq_celltype_matrix) <- freq_celltype$celltype
  freq_celltype_matrix <- freq_celltype_matrix[, cell_type_list]

  # 2.6 Compute the two-cell correlation functions C(u) using the selected depth measure (sel_u)
  depth_values <- sort(unique(unlist(node_pair_depth[, sel_u])))
  Cu <- mclapply(seq_along(depth_values), function(d) {
    current_depth <- depth_values[d]
    cu_subset <- node_pair_depth %>% filter(!!sym(sel_u) == current_depth)
    type_comb_counts <- table(cu_subset$type_combined) %>% as.data.frame()
    colnames(type_comb_counts) <- c("type_combined", "combined_count")
    type_comb_counts$Cu_freq <- type_comb_counts$combined_count / sum(type_comb_counts$combined_count)
    type_comb_counts$u <- current_depth
    return(type_comb_counts)
  }, mc.cores = use_cores) %>% bind_rows()

  # Split the combined cell types into two separate columns
  Cu <- Cu %>% mutate(temp = type_combined) %>% separate(temp, into = c("cell_1", "cell_2"), sep = "_")

  #### 3. Compute the transition matrix for each depth
  all_trans_matrix_list <- list()
  depth_values <- unique(Cu$u)
  # Only consider depths less than or equal to fi_depth
  depth_values <- depth_values[depth_values <= fi_depth]

  lapply(seq_along(depth_values), function(depth_index) {
    current_depth <- depth_values[depth_index]
    cu_subset <- Cu %>% filter(u == current_depth)

    # Build an initial C(u) count matrix with rows and columns corresponding to cell_type_list
    cu_matrix <- matrix(0, nrow = length(cell_type_list), ncol = length(cell_type_list))
    colnames(cu_matrix) <- cell_type_list
    rownames(cu_matrix) <- cell_type_list

    lapply(seq_len(nrow(cu_subset)), function(s) {
      cell_1 <- cu_subset$cell_1[s]
      cell_2 <- cu_subset$cell_2[s]
      count <- cu_subset$combined_count[s]
      cu_matrix[cell_1, cell_2] <<- cu_matrix[cell_1, cell_2] + count
      cu_matrix[cell_2, cell_1] <<- cu_matrix[cell_2, cell_1] + count
      return(NULL)
    })

    # Normalize the C(u) matrix
    cu_matrix <- cu_matrix / sum(cu_matrix)

    # Calculate an intermediate matrix Tmp based on cu_matrix and the cell type frequency matrix
    Tmp <- cu_matrix / (t(sqrt(freq_celltype_matrix)) %*% sqrt(freq_celltype_matrix))
    Tmp <- replace(Tmp, is.nan(Tmp), 0)

    # Perform eigen decomposition on Tmp
    eig_decomp <- eigen(Tmp)
    eigen_values <- eig_decomp$values
    # Replace negative eigenvalues with 0
    eigen_values[eigen_values < 0] <- 0
    eigen_vectors <- eig_decomp$vectors

    # Compute the transition matrix by adjusting eigenvalues based on the current depth
    trans_matrix <- eigen_vectors %*% diag(eigen_values^(1/(2 * current_depth))) %*% t(eigen_vectors) *
      (t(sqrt(freq_celltype_matrix)) %*% (1/sqrt(freq_celltype_matrix)))

    # Use only the real part and set negative or NaN values to 0
    trans_matrix <- Re(trans_matrix)
    trans_matrix[trans_matrix < 0 | is.na(trans_matrix)] <- 0
    trans_matrix[is.nan(trans_matrix)] <- 0

    # Normalize each column of the transition matrix so that each column sums to 1
    trans_matrix <- apply(trans_matrix, 2, function(col) { col / sum(col) })
    trans_matrix[trans_matrix < 0 | is.na(trans_matrix)] <- 0
    trans_matrix[is.nan(trans_matrix)] <- 0

    colnames(trans_matrix) <- cell_type_list
    rownames(trans_matrix) <- cell_type_list

    # Store the transition matrix for the current depth in the list (name format: "depth_current_depth")
    all_trans_matrix_list[[paste0("depth_", current_depth)]] <<- trans_matrix
  })

  #### 4. Optimize to obtain the optimal transition matrix
  # 4.1 Compute the average of all transition matrices as the initial guess
  sum_matrix <- Reduce("+", all_trans_matrix_list)
  mean_matrix <- sum_matrix / length(all_trans_matrix_list)
  mean_matrix_flat <- as.vector(mean_matrix)

  #
  upper_bound <- as.vector(bound_matrix)

  # 4.2 Perform optimization using the L-BFGS-B method
  optim_result <- optim(mean_matrix_flat,
                        calculate_transition_matrix_error,
                        observed_matrices = all_trans_matrix_list,
                        matrix_size = nrow(mean_matrix),
                        type_list = cell_type_list,
                        upper = upper_bound,
                        method = "L-BFGS-B")

  optimal_matrix <- matrix(optim_result$par, nrow = nrow(mean_matrix), byrow = FALSE)
  colnames(optimal_matrix) <- cell_type_list
  rownames(optimal_matrix) <- cell_type_list
  optimal_matrix <- apply(optimal_matrix, 2, function(col) { col / sum(col) })
  optimal_matrix <- replace(optimal_matrix, is.na(optimal_matrix), 0)

  #### Return the list of all transition matrices and the optimized transition matrix
  return(list(all_trans_matrix_list = all_trans_matrix_list,
              optimal_trans_matrix = optimal_matrix))
}
