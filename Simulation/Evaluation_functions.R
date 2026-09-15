############################################################
# Evaluation_functions.R
# Evaluation functions: AUC, module recovery, feature selection, repetition summaries, and plotting.
############################################################

############################################################
# safe_auc()
# Purpose: Compute binary AUC robustly.
############################################################
safe_auc <- function(labels, scores, positive = "AD") {
  labels <- as.factor(labels)
  scores <- as.numeric(scores)
  
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]
  scores <- scores[ok]
  
  if (length(labels) == 0) return(NA_real_)
  if (length(unique(labels)) < 2) return(NA_real_)
  
  y <- as.integer(labels == positive)
  if (length(unique(y)) < 2) return(NA_real_)
  
  out <- tryCatch({
    roc_obj <- pROC::roc(
      response = y,
      predictor = scores,
      direction = "<",
      quiet = TRUE
    )
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
  
  out
}

############################################################
# get_best_median_worst_indices()
# Purpose: Locate representative best, median, and worst simulation repetitions.
############################################################
get_best_median_worst_indices <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  ok <- which(is.finite(x))
  if (length(ok) == 0) {
    return(list(best = NA_integer_, median = NA_integer_, worst = NA_integer_))
  }
  
  x_ok <- x[ok]
  if (higher_better) {
    best_idx <- ok[which.max(x_ok)]
    worst_idx <- ok[which.min(x_ok)]
  } else {
    best_idx <- ok[which.min(x_ok)]
    worst_idx <- ok[which.max(x_ok)]
  }
  
  med_idx <- ok[which.min(abs(x_ok - median(x_ok, na.rm = TRUE)))]
  
  list(
    best = as.integer(best_idx),
    median = as.integer(med_idx),
    worst = as.integer(worst_idx)
  )
}

############################################################
# get_scalar_or_na()
# Purpose: Safely extract one scalar metric from a method-result list.
############################################################
get_scalar_or_na <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!name %in% names(x)) return(NA_real_)
  val <- x[[name]]
  if (length(val) == 0) return(NA_real_)
  as.numeric(val[1])
}

############################################################
# matrix_component_labels()
# Purpose: Label four-neighbor connected components in a binary module matrix.
############################################################
matrix_component_labels <- function(mask) {
  mask <- as.matrix(mask != 0)
  nr <- nrow(mask)
  nc <- ncol(mask)
  lab <- matrix(0L, nrow = nr, ncol = nc)
  cur_lab <- 0L
  
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!mask[i, j] || lab[i, j] != 0L) next
      
      cur_lab <- cur_lab + 1L
      queue_i <- c(i)
      queue_j <- c(j)
      head_idx <- 1L
      lab[i, j] <- cur_lab
      
      while (head_idx <= length(queue_i)) {
        x <- queue_i[head_idx]
        y <- queue_j[head_idx]
        head_idx <- head_idx + 1L
        
        neigh <- rbind(
          c(x - 1L, y),
          c(x + 1L, y),
          c(x, y - 1L),
          c(x, y + 1L)
        )
        
        for (k in seq_len(nrow(neigh))) {
          nx <- neigh[k, 1]
          ny <- neigh[k, 2]
          if (nx >= 1L && nx <= nr && ny >= 1L && ny <= nc) {
            if (mask[nx, ny] && lab[nx, ny] == 0L) {
              lab[nx, ny] <- cur_lab
              queue_i <- c(queue_i, nx)
              queue_j <- c(queue_j, ny)
            }
          }
        }
      }
    }
  }
  
  lab
}

############################################################
# compute_module_discovery_metrics()
# Purpose: Compare the nonzero truth and estimated module masks and compute recovery metrics.
############################################################
compute_module_discovery_metrics <- function(theta_true, bicluster_module) {
  true_mask <- as.matrix(theta_true != 0)
  est_mask  <- as.matrix(bicluster_module != 0)
  
  if (!all(dim(true_mask) == dim(est_mask))) {
    stop("theta_true and bicluster_module must have the same dimension.")
  }
  
  tp <- sum(true_mask & est_mask)
  fp <- sum(!true_mask & est_mask)
  fn <- sum(true_mask & !est_mask)
  tn <- sum(!true_mask & !est_mask)
  
  precision <- ifelse(tp + fp == 0, NA_real_, tp / (tp + fp))
  recall    <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
  f1        <- ifelse(is.na(precision) || is.na(recall) || (precision + recall == 0),
                      NA_real_,
                      2 * precision * recall / (precision + recall))
  jaccard   <- ifelse(tp + fp + fn == 0, NA_real_, tp / (tp + fp + fn))
  specificity <- ifelse(tn + fp == 0, NA_real_, tn / (tn + fp))
  fpr <- ifelse(tn + fp == 0, NA_real_, fp / (tn + fp))
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  
  true_lab <- matrix_component_labels(true_mask)
  est_lab  <- matrix_component_labels(est_mask)
  idx_union <- which(true_mask | est_mask)
  
  ari <- NA_real_
  if (length(idx_union) >= 2) {
    ari <- tryCatch({
      mclust::adjustedRandIndex(true_lab[idx_union], est_lab[idx_union])
    }, error = function(e) NA_real_)
  }
  
  list(
    tp = tp, fp = fp, fn = fn, tn = tn,
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = jaccard,
    specificity = specificity,
    fpr = fpr,
    accuracy = accuracy,
    ari = ari
  )
}

############################################################
# compute_module_metrics_df()
# Purpose: Arrange per-repetition module metric lists into a data frame.
############################################################
compute_module_metrics_df <- function(module_metrics_list) {
  data.frame(
    rep = seq_along(module_metrics_list),
    ARI = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "ari")),
    Jaccard = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "jaccard")),
    F1 = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "f1")),
    Precision = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "precision")),
    Recall = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "recall")),
    FPR = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "fpr")),
    Accuracy = sapply(module_metrics_list, function(x) get_scalar_or_na(x, "accuracy"))
  )
}

#####################Visualization#######################################

.plot_heatmap_cont <- function(mat, main = "", zlim = NULL,
                               xlab = "Brain functional Connections",
                               ylab = "Genes") {
  mat <- as.matrix(mat)
  if (is.null(zlim)) {
    zmax <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(zmax) || zmax == 0) zmax <- 1
    zlim <- c(-zmax, zmax)
  }
  
  par(bg = "white")
  image(
    x = seq_len(nrow(mat)),
    y = seq_len(ncol(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(160),
    zlim = zlim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    useRaster = TRUE
  )
  box(col = "grey35", lwd = 0.8)
  mtext(xlab, side = 1, line = 2.0, cex = 0.95)
  mtext(ylab, side = 2, line = 2.0, cex = 0.95)
}


.plot_heatmap_bin <- function(mat, main = "",
                              xlab = "Brain functional Connections",
                              ylab = "Genes") {
  mat <- as.matrix(mat)
  par(bg = "white")
  image(
    x = seq_len(nrow(mat)),
    y = seq_len(ncol(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = c("white", "#e41a1c"),
    breaks = c(-0.5, 0.5, 1.5),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    useRaster = TRUE
  )
  box(col = "grey35", lwd = 0.8)
  mtext(xlab, side = 1, line = 2.0, cex = 0.95)
  mtext(ylab, side = 2, line = 2.0, cex = 0.95)
}

############################################################
# save_triptych_heatmap()。
# Purpose: Save side-by-side heatmaps of true theta, estimated theta, and identified modules.
############################################################
save_triptych_heatmap <- function(theta_true, theta_hat, bicluster_module, file, main_prefix = "") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  png(file, width = 2100, height = 760, res = 170, bg = "white")
  par(mfrow = c(1, 3), mar = c(3.2, 3.2, 3.2, 1.0), oma = c(0, 0, 0, 0))
  
  zmax <- max(abs(c(theta_true, theta_hat)), na.rm = TRUE)
  if (!is.finite(zmax) || zmax == 0) zmax <- 1
  
  .plot_heatmap_cont(theta_true, paste0(main_prefix, "True Theta"), zlim = c(-zmax, zmax))
  .plot_heatmap_cont(theta_hat, paste0(main_prefix, "Estimated Theta Hat"), zlim = c(-zmax, zmax))
  .plot_heatmap_bin((bicluster_module != 0) * 1, paste0(main_prefix, "Identified Modules"))
  
  dev.off()
}

############################################################
# compute_frequency_matrix()
# Purpose: Compute how often every imaging-gene cell is selected as a module across repetitions.
############################################################
compute_frequency_matrix <- function(bicluster_module_list) {
  dims <- lapply(bicluster_module_list, dim)
  nr <- unique(sapply(dims, `[`, 1))
  nc <- unique(sapply(dims, `[`, 2))
  
  if (length(nr) != 1 || length(nc) != 1) {
    stop("All bicluster_module matrices must have the same dimension.")
  }
  
  arr <- simplify2array(lapply(bicluster_module_list, function(x) (x != 0) * 1))
  if (length(dim(arr)) == 2) return(arr)
  apply(arr, c(1, 2), mean, na.rm = TRUE)
}

############################################################
# save_frequency_heatmap()
# Purpose: Save the module-selection frequency matrix as a heatmap.
############################################################
save_frequency_heatmap <- function(freq_mat, file, main = "Module Frequency Heatmap") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  png(file, width = 900, height = 780, res = 170, bg = "white")
  par(mar = c(3.2, 3.2, 3.2, 1.0))
  image(
    x = seq_len(nrow(freq_mat)),
    y = seq_len(ncol(freq_mat)),
    z = t(freq_mat[nrow(freq_mat):1, , drop = FALSE]),
    col = colorRampPalette(c("white", "#fee08b", "#f46d43", "#d73027"))(140),
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    useRaster = TRUE
  )
  box(col = "grey35", lwd = 0.8)
  mtext("Brain functional Connections", side = 1, line = 2.0, cex = 0.95)
  mtext("Genes", side = 2, line = 2.0, cex = 0.95)
  dev.off()
}

############################################################
# save_module_metrics_boxplot()
# Purpose: Save boxplots of module-recovery metrics across repetitions.
############################################################
save_module_metrics_boxplot <- function(module_metrics_df, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  png(file, width = 1300, height = 760, res = 170, bg = "white")
  par(mar = c(7, 4, 3, 1))
  boxplot(
    module_metrics_df[, c("ARI", "Jaccard", "F1", "Precision", "Recall", "FPR")],
    las = 2,
    col = "grey90",
    border = "grey35",
    main = "Module Recovery Metrics Across Repetitions",
    ylab = "Metric Value"
  )
  dev.off()
}

############################################################
# compute_tpr_fpr_from_selected()
# Purpose: Compute feature-selection TPR and FPR.
############################################################
compute_tpr_fpr_from_selected <- function(selected_all, p_total, true_active) {
  selected_all <- sort(unique(selected_all))
  selected_all <- selected_all[selected_all >= 1 & selected_all <= p_total]
  
  true_inactive <- setdiff(seq_len(p_total), true_active)
  
  tpr <- if (length(true_active) == 0) {
    NA_real_
  } else {
    length(intersect(true_active, selected_all)) / length(true_active)
  }
  
  fpr <- if (length(true_inactive) == 0) {
    NA_real_
  } else {
    length(intersect(true_inactive, selected_all)) / length(true_inactive)
  }
  
  list(tpr = tpr, fpr = fpr)
}

############################################################
# bridge_selection_from_bootstrap()
# Purpose: Calculate the complete TPR/FPR curve over bootstrap-frequency
# cutoffs. The common threshold is chosen later from the across-repetition
# mean curve, not separately within each repetition.
############################################################
bridge_selection_from_bootstrap <- function(sel_mat, true_active, p_total) {
  B <- nrow(sel_mat)
  freq <- colSums(sel_mat, na.rm = TRUE)
  
  curve_df <- data.frame(
    threshold_count = seq_len(B),
    selected_size = NA_integer_,
    TPR = NA_real_,
    FPR = NA_real_,
    Youden = NA_real_
  )
  
  for (k in seq_len(B)) {
    selected_all <- which(freq >= k)
    met <- compute_tpr_fpr_from_selected(selected_all, p_total, true_active)
    curve_df$selected_size[k] <- length(selected_all)
    curve_df$TPR[k] <- met$tpr
    curve_df$FPR[k] <- met$fpr
    curve_df$Youden[k] <- met$tpr - met$fpr
  }
  
  list(
    selection_frequency = freq / B,
    curve = curve_df,
    threshold_rule = "selected_later_from_across_repetition_mean_curve",
    best_threshold_count = NA_integer_,
    selected_all = integer(0),
    tpr = NA_real_,
    fpr = NA_real_
  )
}
#######################Summarize results#####################################

summarize_method_results_auc <- function(result_list) {
  per_rep <- data.frame(
    rep = seq_along(result_list),
    train_acc = sapply(result_list, function(x) get_scalar_or_na(x, "train_acc")),
    test_acc  = sapply(result_list, function(x) get_scalar_or_na(x, "test_acc")),
    train_auc = sapply(result_list, function(x) get_scalar_or_na(x, "train_auc")),
    test_auc  = sapply(result_list, function(x) get_scalar_or_na(x, "test_auc")),
    tpr = sapply(result_list, function(x) get_scalar_or_na(x, "tpr")),
    fpr = sapply(result_list, function(x) get_scalar_or_na(x, "fpr"))
  )
  

  #  maximize the Youden index of the averaged curve;
  valid_curve_id <- which(vapply(result_list, function(x) {
    !is.null(x) &&
      !is.null(x$selection_eval) &&
      !is.null(x$selection_eval$curve) &&
      nrow(x$selection_eval$curve) > 0
  }, logical(1)))
  
  selected_threshold_count <- NA_integer_
  selected_threshold_frequency <- NA_real_
  selected_sets_at_common_threshold <- vector("list", length(result_list))
  
  if (length(valid_curve_id) > 0) {
    common_thresholds <- Reduce(
      intersect,
      lapply(valid_curve_id, function(r) {
        result_list[[r]]$selection_eval$curve$threshold_count
      })
    )
    common_thresholds <- sort(unique(as.integer(common_thresholds)))
    
    if (length(common_thresholds) > 0) {
      tpr_curve <- matrix(
        NA_real_,
        nrow = length(result_list),
        ncol = length(common_thresholds)
      )
      fpr_curve <- matrix(
        NA_real_,
        nrow = length(result_list),
        ncol = length(common_thresholds)
      )
      
      for (r in valid_curve_id) {
        z <- result_list[[r]]$selection_eval$curve
        pos <- match(common_thresholds, z$threshold_count)
        tpr_curve[r, ] <- as.numeric(z$TPR[pos])
        fpr_curve[r, ] <- as.numeric(z$FPR[pos])
      }
      
      mean_tpr_curve <- colMeans(tpr_curve, na.rm = TRUE)
      mean_fpr_curve <- colMeans(fpr_curve, na.rm = TRUE)
      mean_youden_curve <- mean_tpr_curve - mean_fpr_curve
      
      ok <- which(is.finite(mean_youden_curve))
      if (length(ok) > 0) {
        best_position <- ok[which.max(mean_youden_curve[ok])]
        selected_threshold_count <- common_thresholds[best_position]
        selected_threshold_frequency <-
          selected_threshold_count / max(common_thresholds)
        

        per_rep$tpr <- tpr_curve[, best_position]
        per_rep$fpr <- fpr_curve[, best_position]
        
        for (r in valid_curve_id) {
          selection_frequency <-
            result_list[[r]]$selection_eval$selection_frequency
          bootstrap_number <-
            nrow(result_list[[r]]$selection_eval$curve)
          selected_sets_at_common_threshold[[r]] <- which(
            selection_frequency * bootstrap_number >=
              selected_threshold_count - sqrt(.Machine$double.eps)
          )
        }
      }
    }
  }
  
  overall <- data.frame(
    metric = c("train_acc", "test_acc", "train_auc", "test_auc", "TPR", "FPR"),
    mean = c(
      mean(per_rep$train_acc, na.rm = TRUE),
      mean(per_rep$test_acc, na.rm = TRUE),
      mean(per_rep$train_auc, na.rm = TRUE),
      mean(per_rep$test_auc, na.rm = TRUE),
      mean(per_rep$tpr, na.rm = TRUE),
      mean(per_rep$fpr, na.rm = TRUE)
    ),
    sd = c(
      sd(per_rep$train_acc, na.rm = TRUE),
      sd(per_rep$test_acc, na.rm = TRUE),
      sd(per_rep$train_auc, na.rm = TRUE),
      sd(per_rep$test_auc, na.rm = TRUE),
      sd(per_rep$tpr, na.rm = TRUE),
      sd(per_rep$fpr, na.rm = TRUE)
    )
  )
  
  list(
    per_rep = per_rep,
    overall = overall,
    feature_threshold_rule = "Youden index of the across-repetition mean TPR/FPR curve",
    feature_threshold_count = selected_threshold_count,
    feature_threshold_frequency = selected_threshold_frequency,
    selected_sets_at_common_threshold = selected_sets_at_common_threshold
  )
}


make_method_summary_long <- function(method_results_named_list) {
  out <- list()
  for (nm in names(method_results_named_list)) {
    tmp <- summarize_method_results_auc(method_results_named_list[[nm]])$overall
    tmp$method <- nm
    out[[nm]] <- tmp
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
