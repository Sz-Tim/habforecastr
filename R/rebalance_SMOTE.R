
#' Generate Synthetic Samples Using Modified SMOTE
#'
#' This function generates synthetic samples for the minority class using a modified version of SMOTE
#' from the smotefamily package. It then combines the original and synthetic data.
#'
#' @param og_df A data frame containing the original dataset.
#' @param K The number of nearest neighbors to use for generating synthetic samples.
#' @param dup_size The desired size of the synthetic dataset.
#'
#' @return A data frame containing the original and synthetic data.
#' @import dplyr
#' @importFrom lubridate ymd
#' @importFrom smotefamily knearest
#' @export
#'
#' @examples
#' \dontrun{
#' og_df <- data.frame(
#'   siteid = rep(1:10, each = 10),
#'   lon = runif(100, -180, 180),
#'   lat = runif(100, -90, 90),
#'   obsid = 1:100,
#'   y = rnorm(100),
#'   date = Sys.Date() - 1:100,
#'   alert = sample(c("A1", "A2"), 100, replace = TRUE),
#'   alert1 = sample(1:2, 100, replace = TRUE),
#'   alert2 = sample(1:2, 100, replace = TRUE),
#'   tl = sample(1:5, 100, replace = TRUE)
#' )
#' K <- 5
#' dup_size <- 100
#' result <- generate_synthetic_samples(og_df, K, dup_size)
#' }
rebalance_smote <- function (og_df, K = 5, dup_size = 0) {
  # Extract unique site coordinates
  site_coords <- og_df |>
    select(siteid, lon, lat) |>
    slice_head(n = 1, by = siteid)

  # Prepare the dataset for k-nearest neighbors (KNN) calculation
  knn_df <- og_df |>
    select(-obsid, -y, -date, -alert, -lon, -lat) |>
    mutate(alert1 = as.numeric(alert1),
           alert2 = as.numeric(alert2),
           tl = as.numeric(tl))

  # Define the target variable
  target <- og_df$alert

  # Number of columns in the KNN dataframe
  ncD <- ncol(knn_df)

  # Count the number of instances in each class
  n_target <- table(target)

  # Create the positive (minority) and negative (majority) sets
  P_set <- subset(knn_df, target == "A1")[sample(min(n_target)), ]
  N_set <- subset(knn_df, target != "A1")

  # Sizes of the positive and negative sets
  sizeP <- nrow(P_set)
  sizeN <- nrow(N_set)

  # Determine the number of nearest neighbors (K)
  K <- min(sizeP - 1, K)

  # Calculate the K-nearest neighbors for the positive set
  knear <- knearest(P_set, P_set, K)

  # Calculate the number of synthetic samples to generate
  sum_dup <- n_dup_max(sizeP + sizeN, sizeP, sizeN, dup_size)

  # Initialize an empty dataframe for synthetic data
  syn_dat <- NULL

  # Column index for the site ID
  site_col <- which(names(knn_df) == "siteid")

  # Generate synthetic samples
  for (i in 1:sizeP) {
    if (is.matrix(knear)) {
      pair_idx <- knear[i, ceiling(runif(sum_dup) * K)]
    } else {
      pair_idx <- rep(knear[i], sum_dup)
    }
    g <- runif(sum_dup)
    P_i <- matrix(unlist(P_set[i, ]), sum_dup, ncD, byrow = TRUE)
    Q_i <- as.matrix(P_set[pair_idx, ])
    syn_i <- P_i + g * (Q_i - P_i)
    site_binom <- rbinom(nrow(syn_i), 1, prob = g)
    syn_i[, 3] <- P_i[, 3] * (1 - site_binom) + Q_i[, 3] * site_binom
    syn_dat <- rbind(syn_dat, syn_i)
  }

  # Convert synthetic data to a dataframe and merge with site coordinates
  syn_df <- as.data.frame(syn_dat) |>
    set_names(names(knn_df)) |>
    left_join(site_coords) |>
    mutate(year = round(year),
           yday = round(yday),
           date = ymd(paste(year, "-01-01")) + yday - 1,
           alert = factor("A1", levels = levels(og_df$alert)),
           alert1 = round(alert1) |> factor(levels = c(1, 2), labels = levels(og_df$alert1)),
           alert2 = round(alert2) |> factor(levels = c(1, 2), labels = levels(og_df$alert2)),
           tl = round(tl) |> factor(levels = seq_along(levels(og_df$tl)), labels = levels(og_df$tl), ordered = TRUE),
           y = og_df$y[1],
           obsid = 10000 + row_number())

  # Combine the original and synthetic data
  return(bind_rows(og_df, syn_df))
}




#' Find k-Nearest Neighbors
#'
#' This function finds the k-nearest neighbors for each point in the dataset using the FNN package.
#'
#' @param D A data frame or matrix representing the dataset.
#' @param P A data frame or matrix representing the points for which to find the nearest neighbors.
#' @param n_clust The number of nearest neighbors to find.
#'
#' @return A matrix where each row contains the indices of the k-nearest neighbors for each point in P.
#' @importFrom FNN knnx.index
#' @export
#'
#' @examples
#' \dontrun{
#' D <- matrix(runif(100), ncol = 2)
#' P <- matrix(runif(20), ncol = 2)
#' n_clust <- 3
#' neighbors <- knearest(D, P, n_clust)
#' }
knearest <- function(D, P, n_clust) {
  # Check if the FNN package is installed, and install it if necessary
  if (!requireNamespace("FNN", quietly = TRUE)) {
    install.packages("FNN", quiet = TRUE)
  }

  # Find the k-nearest neighbors using the FNN package
  knD <- FNN::knnx.index(D, P, k = (n_clust + 1), algorithm = "kd_tree")

  # Remove self-references (where a point is its own nearest neighbor)
  knD <- knD * (knD != row(knD))

  # Identify rows where the first neighbor is non-zero
  que <- which(knD[, 1] > 0)

  # Replace zero entries with the first neighbor's index
  for (i in que) {
    knD[i, which(knD[i, ] == 0)] <- knD[i, 1]
    knD[i, 1] <- 0
  }

  # Return the indices of the k-nearest neighbors
  return(knD[, 2:(n_clust + 1)])
}




#' Calculate Maximum Number of Duplicates
#'
#' This function calculates the maximum number of duplicates to generate based on the input sizes and duplication size.
#'
#' @param size_input The total size of the input dataset.
#' @param size_P The size of the positive (minority) class.
#' @param size_N The size of the negative (majority) class.
#' @param dup_size The desired duplication size. Default is 0.
#'
#' @return The maximum number of duplicates to generate.
#' @export
#'
#' @examples
#' \dontrun{
#' size_input <- 1000
#' size_P <- 50
#' size_N <- 950
#' dup_size <- 0
#' max_dup <- n_dup_max(size_input, size_P, size_N, dup_size)
#' }
n_dup_max <- function(size_input, size_P, size_N, dup_size = 0) {
  # Check if dup_size is a vector with more than one element
  if (is.vector(dup_size) && length(dup_size) > 1) {
    # If any element of dup_size is 0, calculate sizeM based on the input sizes
    if (length(which(dup_size == 0)) > 0) {
      sizeM <- floor((2 * size_N - size_input) / size_P)
    }
    # If no element of dup_size is 0, set sizeM to the maximum value in dup_size
    if (length(which(dup_size == 0)) == 0) {
      sizeM <- max(dup_size)
    }
  }

  # Check if dup_size is not a vector or has only one element
  if (!is.vector(dup_size) || length(dup_size) == 1) {
    # If dup_size is 0, calculate sizeM based on the input sizes
    if (dup_size == 0) {
      sizeM <- floor((2 * size_N - size_input) / size_P)
    }
    # If dup_size is not 0, set sizeM to dup_size
    if (dup_size != 0) {
      sizeM <- dup_size
    }
  }

  # Return the calculated maximum number of duplicates
  return(sizeM)
}
