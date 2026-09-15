rm(list = ls())
#set.seed()
# setwd("")

############################################################
# 1. Packages and function files
############################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(pROC)
  library(mclust)
  library(sparcl)
  library(MASS)
  library(glmnet)
  library(huge)
  library(expm)
  library(mnormt)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(DensParcorr)
  library(sparsegl)
  library(igraph)
})

source("SW function.R")
source("bridge_functions.R")
source("Evaluation_functions.R")

############################################################
# 2. Parameter settings
############################################################

l <- 500
p <- 80  
q <- 50
N <- 150
times <- 100
train_ratio <- 0.7

N_tr <- as.integer(round(train_ratio * N))
N_ts <- N - N_tr

############################################################
# 3. Parallel 
############################################################

cl <- parallel::makeCluster(
  max(1L, min(8L, parallel::detectCores() - 1L))
)
doParallel::registerDoParallel(cl)

clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, {
  library(Matrix)
  library(sparcl)
  library(MASS)
  library(glmnet)
  library(huge)
  library(expm)
  library(mnormt)
  library(foreach)
  library(DensParcorr)
  library(sparsegl)

  source("SW function.R")
  source("bridge_functions.R")
  source("Evaluation_functions.R")
  NULL
})

############################################################
# 4. Generate temporal covariance matrices and spatial networks for both groups
############################################################

temporal_obj <- build_temporal_cov(q = q, sigmaT = "band")

spatial_obj <- build_spatial_base(
  p = p,
  spatial_type = "smallworld",
  hub_g = 5,
  sw_m = 10,
  sw_banded_n = 6,
  sw_source = "SW function.R"
)

delta_image <- which(spatial_obj$delta_vec != 0)

subject_precision <- generate_subject_precision(
  N = N,
  p = p,
  Theta1 = spatial_obj$Theta1,
  omega1.total = spatial_obj$omega1.total,
  omega2.total = spatial_obj$omega2.total,
  noise_sd = 0.2
)

############################################################
# 5. Result containers
############################################################

bridge_results <- vector("list", times)
module_metrics_list <- vector("list", times)
screen_hit_rate <- rep(NA_real_, times)

############################################################
# 6.  Independent simulation repetitions
# Perform the outer train/test split before standardization or screening
# The test set is standardized using the means and standard deviations estimated from the outer training set. 
# PCA loadings are also estimated only from the training set. The test data are projected onto the training-derived PCA loadings.
############################################################

for (r in seq_len(times)) {
  set.seed(r)
  cat("BRIDGE repetition", r, "of", times, "\n")
  
  ## Generate brain imaging data
  X_AD <- simulate_group_X(
    omega_list = subject_precision$omega1_N,
    sigmaT = temporal_obj$sigmaT1,
    N = N,
    p = p
  )
  Omega_AD <- estimate_group_precision(X_list = X_AD, N = N)

  X_HC <- simulate_group_X(
    omega_list = subject_precision$omega2_N,
    sigmaT = temporal_obj$sigmaT2,
    N = N,
    p = p
  )
  Omega_HC <- estimate_group_precision(X_list = X_HC, N = N)

  B_AD_raw <- vectorize_fisher_features(Omega_AD, N = N, p = p)$X_vec
  B_HC_raw <- vectorize_fisher_features(Omega_HC, N = N, p = p)$X_vec


  train_id <- sample(seq_len(N), N_tr, replace = FALSE)
  test_id <- setdiff(seq_len(N), train_id)

  B_train_raw <- rbind(
    B_AD_raw[train_id, , drop = FALSE],
    B_HC_raw[train_id, , drop = FALSE]
  )
  B_test_raw <- rbind(
    B_AD_raw[test_id, , drop = FALSE],
    B_HC_raw[test_id, , drop = FALSE]
  )

  diagnosis_train <- factor(
    c(rep("AD", N_tr), rep("HC", N_tr)),
    levels = c("HC", "AD")
  )
  diagnosis_test <- factor(
    c(rep("AD", N_ts), rep("HC", N_ts)),
    levels = c("HC", "AD")
  )
  
  screen_obj <- screen_image_train_test(
    B_train_raw = B_train_raw,
    B_test_raw = B_test_raw,
    y_train = as.numeric(diagnosis_train == "AD"),
    top_k = 500
  )

  B_train <- screen_obj$B_train
  B_test <- screen_obj$B_test
  order_image_number <- screen_obj$selected_index

  screen_hit_rate[r] <-
    length(intersect(order_image_number, delta_image)) /
    length(delta_image)


  # Construct theta matrix and generate gene expression data
  pB <- ncol(B_train)
  pG <- l
  p_total <- pB + pG

  if (pB != 500L || pG != 500L) {
    stop("This example requires 500 screened imaging features and 500 genes.")
  }

  theta <- generate_theta_blocks(
    pB = pB,
    pG = pG,
    block_setting = "block10"
  )

  G_train_raw <- B_train %*% theta + matrix(
    rnorm(nrow(B_train) * pG),
    nrow = nrow(B_train),
    ncol = pG
  )
  G_test_raw <- B_test %*% theta + matrix(
    rnorm(nrow(B_test) * pG),
    nrow = nrow(B_test),
    ncol = pG
  )

  gene_obj <- standardize_gene_train_test(G_train_raw, G_test_raw)
  G_train_std <- gene_obj$G_train
  G_test_std <- gene_obj$G_test

  true_image <- get_delta_vec_number(
    delta_image = delta_image,
    order_image_number = order_image_number
  )

  if (length(true_image) > 0L) {
    true_gene_raw <- unique(which(
      theta[true_image, , drop = FALSE] != 0,
      arr.ind = TRUE
    )[, 2])
  } else {
    true_gene_raw <- integer(0)
  }

  true_active <- sort(unique(c(true_image, pB + true_gene_raw)))


  # multivariate sparse regression to estimate theta and iterative biclustering to identify modules

  invisible(capture.output(
    discovery_obj <- tryCatch(
      run_module_discovery_once(
        B_train = B_train,
        G_train = G_train_raw,
        S = 10,
        alpha = 1,
        nfolds = 5
      ),
      error = function(e) {
        list(
          error = conditionMessage(e),
          theta_hat = matrix(0, nrow = pB, ncol = pG),
          bicluster_module = matrix(0, nrow = pB, ncol = pG),
          lay_reduce = matrix(0, nrow = pB, ncol = pG),
          module_number = vector("list", 10),
          module_number_gene = vector("list", 10),
          lay_list = vector("list", 10),
          theta_list = vector("list", 10)
        )
      }
    )
  ))

  discovery_obj$module_metrics <- compute_module_discovery_metrics(
    theta_true = theta,
    bicluster_module = discovery_obj$bicluster_module
  )
  module_metrics_list[[r]] <- discovery_obj$module_metrics

  # Module-wise PCA and BRIDGE bootstrap classification

  method_data <- build_bridge_method_data(
    B_train = B_train,
    G_train = G_train_std,
    B_test = B_test,
    G_test = G_test_std,
    discovery_obj = discovery_obj,
    var_explained = 0.90
  )

  full_fit <- tryCatch(
    run_bridge_bootstrap_method(
      method_data = method_data,
      y_train = diagnosis_train,
      y_test = diagnosis_test,
      boot_strap = 100,
      true_active = true_active,
      p_total = p_total,
      tune_nfolds = 10,
      seed = 1000 + r
    ),
    error = function(e) {
      list(
        method = "BRIDGE",
        error = conditionMessage(e),
        train_acc = NA_real_,
        test_acc = NA_real_,
        train_auc = NA_real_,
        test_auc = NA_real_,
        selection_eval = NULL
      )
    }
  )


  bridge_results[[r]] <- list(
    method = full_fit$method,
    error = full_fit$error,
    train_acc = full_fit$train_acc,
    test_acc = full_fit$test_acc,
    train_auc = full_fit$train_auc,
    test_auc = full_fit$test_auc,
    selection_eval = full_fit$selection_eval
  )

  invisible(gc())
}

parallel::stopCluster(cl)

############################################################
# 7. Results
############################################################

bridge_summary <- summarize_method_results_auc(bridge_results)

classification_selection_summary <- bridge_summary$overall[
  bridge_summary$overall$metric %in% c("test_acc", "test_auc", "TPR", "FPR"),
  ,
  drop = FALSE
]

module_metrics_df <- compute_module_metrics_df(module_metrics_list)
module_metric_names <- setdiff(names(module_metrics_df), "rep")

module_summary <- do.call(rbind, lapply(module_metric_names, function(metric) {
  values <- module_metrics_df[[metric]]
  values <- values[is.finite(values)]

  data.frame(
    metric = metric,
    mean = if (length(values) > 0L) mean(values) else NA_real_,
    sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
    row.names = NULL
  )
}))

cat("\nBRIDGE classification and feature-selection summary\n")
print(classification_selection_summary, row.names = FALSE)

cat("\nBRIDGE module-recovery summary\n")
print(module_summary, row.names = FALSE)

cat("\nMean imaging-screening hit rate:",
    round(mean(screen_hit_rate, na.rm = TRUE), 4), "\n")
