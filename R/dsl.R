#' Estimating Regression using the DSL framework
#' @param model A regression model \code{dsl} currently supports \code{lm} (linear regression), \code{logit} (logistic regression), \code{poisson} (Poisson regression for count outcomes, estimated as Poisson Pseudo Maximum Likelihood), \code{felm} (fixed-effects regression), \code{fepois} (fixed-effects Poisson regression), \code{negbin} (negative binomial regression for overdispersed count outcomes), and \code{fenegbin} (fixed-effects negative binomial regression).
#' @param formula A formula used in the specified regression model. For count models (\code{poisson}, \code{fepois}, \code{negbin}, and \code{fenegbin}), it can include an offset, e.g., \code{Y ~ X + offset(log(exposure))}, where the offset cannot contain variables in \code{predicted_var}.
#' @param predicted_var A vector of column names in the data that correspond to variables that need to be predicted. Predicted variables should be numeric. For count models, a predicted independent variable can also be categorical (factor, character, or logical) or binary (numeric 0/1). \code{dsl} then predicts the probability of each level and uses the expectation of the moment conditions over the predicted levels.
#' @param prediction A vector of column names in the data that correspond to predictions of \code{predicted_var}.
#' @param data A data frame. The class should be \code{data.frame}.
#' @param cluster A column name in the data that indicates the level at which cluster standard errors are calculated. Default is \code{NULL}. (\code{fepois} and \code{fenegbin} with predicted independent variables) When \code{cluster} is specified, standard errors are clustered by the smallest clusters that contain both the groups of the fixed effects and the clusters in \code{cluster} (\code{summary()} reports it when they differ from \code{cluster}).
#' @param labeled (Optional) A column name in the data that indicates which observation is labeled. It should be a vector of 1 (labeled) and 0 (non-labeled). When \code{NULL}, the function assumes that observations that have \code{NA} in \code{predicted_var} are non-labeled and other observations are labeled.
#' @param sample_prob (Optional) A column name in the data that correspond to the sampling probability for labeling a particular observation. When \code{NULL}, the function assumes random sampling with equal probabilities.
#' @param fixed_effect (Used when \code{model = "felm"}, \code{model = "fepois"}, or \code{model = "fenegbin"}) A type of fixed effects regression you run. \code{oneway} (one-way fixed effects) or \code{twoways} (two-way fixed effects).
#' @param index (Used when \code{model = "felm"}, \code{model = "fepois"}, or \code{model = "fenegbin"}) A vector of column names specifying fixed effects. When \code{fixed_effect = oneway}, it has one element. When \code{fixed_effect = twoways}, it has two elements, e.g., \code{index = c("state", "year")}.
#' @param sl_method A name of a supervised machine learning model used internally to predict \code{predicted_var} by fine-tuning \code{prediction} or using predictors (specified in \code{feature}) when \code{prediction = NULL}. Users can run \code{available_method()} to see available supervised machine learning methods. Default is \code{grf} (generalized random forest).
#' @param feature A vector of column names in the data that correspond to predictors used to fit a supervised machine learning (specified in \code{sl_method}). When independent variables are predicted, including the outcome and the other independent variables in \code{feature} usually makes the estimates much more precise.
#' @param family (Used when making predictions with \code{SuperLearner} methods) \code{gaussian} (default) or \code{binomial} (for binary \code{predicted_var}). For count variables, use \code{gaussian}. Predictions are kept within the range of the labeled values of \code{predicted_var}.
#' @param cross_fit The fold of cross-fitting. Default is \code{5}.
#' @param sample_split The number of sampling-splitting. Default is \code{10}. (Count models) When the moment conditions cannot be solved in some sample splits, these sample splits are dropped with a warning.
#' @param seed Numeric \code{seed} used internally. Default is \code{1234}.
#' @rawNamespace import(Matrix, except = summary)
#' @importFrom grf regression_forest
#' @importFrom estimatr lm_robust
#' @importFrom matrixcalc is.positive.definite
#' @importFrom arm model.matrixBayes
#' @importFrom stats .getXlevels as.formula ave delete.response glm lm median model.frame model.matrix model.offset model.response optim predict sd pnorm qnorm relevel terms uniroot var
#' @importFrom utils capture.output
#' @importFrom graphics points
#' @import tidyverse
#' @import SuperLearner
#' @return \code{dsl} returns an object of \code{dsl} class.
#'  \itemize{
#'    \item \code{coefficients}: Estimated coefficients.
#'    \item \code{standard_errors}: Estimated standard errors.
#'    \item \code{vcov}: Estimated variance-covariance matrix.
#'    \item \code{RMSE}: Root mean squared error in the internal prediction step.
#'    \item \code{theta}: (\code{model = "negbin"} or \code{"fenegbin"}) Estimated dispersion parameter of the negative binomial distribution, where \code{Var(Y | X) = mu + mu^2/theta}.
#'    \item \code{internal}: Outputs used only for the internal use.
#'  }
#' @export

dsl <- function(model = "lm",
                formula,
                predicted_var,
                prediction = NULL,
                data,
                cluster = NULL,
                labeled = NULL,
                sample_prob = NULL,
                index = NULL,
                fixed_effect = "oneway",
                sl_method = "grf",
                feature = NULL,
                family = "gaussian",
                cross_fit = 5,
                sample_split = 10,
                seed = 1234){

  # ##################
  # Setup
  # ##################
  optim_method <-  "L-BFGS-B"
  lambda <-  0

  # data.frame
  class(data) <- "data.frame"

  if((model %in% c("lm", "logit", "poisson", "felm", "fepois", "negbin", "fenegbin")) == FALSE){
    stop(" `model` should be either `lm`, `logit`, `poisson`, `felm`, `fepois`, `negbin`, or `fenegbin` ")
  }

  if(is.null(prediction) & is.null(feature)){
    stop(" Please specify at least one of `prediction` and `feature` ")
  }

  if(is.function(family) == TRUE){
    family <- family()
  }
  if(inherits(family, "family") == TRUE){
    family <- family$family
  }
  if(isTRUE(is.character(family) && length(family) == 1 && family %in% c("gaussian", "binomial")) == FALSE){
    stop(" `family` should be either `gaussian` or `binomial`. For count variables, please use `gaussian`. ")
  }
  count_model <- model %in% c("poisson", "fepois", "negbin", "fenegbin")

  # Offset (count models only)
  terms_formula <- terms(formula)
  if(is.null(attr(terms_formula, "offset")) == FALSE){
    if(count_model == FALSE){
      stop(" `offset()` in `formula` is supported only for count models (`poisson`, `fepois`, `negbin`, and `fenegbin`). ")
    }
    offset_var <- unlist(lapply(attr(terms_formula, "offset"), function(j) all.vars(attr(terms_formula, "variables")[[j + 1]])))
    if(any(offset_var %in% predicted_var)){
      stop(" `offset()` in `formula` cannot contain variables in `predicted_var`. ")
    }
  }

  # Model-Specific Check
  if(model %in% c("felm", "fepois", "fenegbin")){
    if(is.null(index) == TRUE){
      stop(paste0("Please specify `index` for `model = ", model, "` "))
    }
    if(fixed_effect == "twoways"){
      if(length(index) != 2){
        stop("Please specify two indices for the two-way fixed effects estimator.")
      }
    }else if(fixed_effect == "oneway"){
      index_use <- index[1]
    }
  }
  if(model %in% c("fepois", "fenegbin")){
    if((fixed_effect %in% c("oneway", "twoways")) == FALSE){
      stop(" `fixed_effect` should be either `oneway` or `twoways` ")
    }
    if(attr(terms(formula), "intercept") == 0){
      stop(paste0(" `formula` should include an intercept for `model = ", model, "` (the intercept is absorbed by fixed effects). "))
    }
  }

  # labeled
  if(is.null(labeled) == TRUE){
    label_index <- 1 - as.numeric(apply(data[, predicted_var, drop = FALSE], 1, function(x) any(is.na(x))))
    data$labeled___ <- label_index
    labeled <- "labeled___"

    if(mean(data[,labeled]) == 1){
      warnings(" Please use `NA` to indicate observations that do not have labels. ")
    }
  }

  # sample_prob
  if(is.null(sample_prob) == TRUE){
    sample_prob0 <- sum(data[, labeled], na.rm = TRUE)/nrow(data)
    data$sample_prob___ <- sample_prob0
    sample_prob <- "sample_prob___"
    equal_prob <- TRUE
  }else{
    if(var(data[, sample_prob]) == 0){
      equal_prob <- TRUE
    }else{
      equal_prob <- FALSE
    }
  }

  # Setup feature
  covariates_use <- unique(c(prediction, feature))
  covariates_use <- setdiff(covariates_use, predicted_var)

  # Handle Methods
  if(any(sl_method != "grf")){
    all_method <- available_method(print_out = FALSE)
    if(any((sl_method %in% all_method) == FALSE)){
      stop(" `sl_method` is not supported in `dsl`. Please run `available_method() to see what `sl_method` is supported. ")
    }
  }

  # Handle cluster standard errors
  if(is.null(cluster) == TRUE){
    data$cluster___ <- seq(1:nrow(data))
    clustered <- FALSE
  }else{
    data$cluster___ <- as.numeric(as.factor(data[, cluster]))
    clustered <- TRUE
  }

  # Handling NAs
  if(any(is.na(data[data[, labeled] == 1, predicted_var]))){
    stop("NA in `predicted_var` of labeled data")
  }
  keep_var <- setdiff(unique(c("cluster___", sample_prob, labeled, covariates_use, all.vars(formula))), predicted_var)
  if(model %in% c("fepois", "fenegbin")){
    keep_var <- unique(c(keep_var, index))
  }
  remove <- apply(data[, keep_var, drop = FALSE], 1, function(x) any(is.na(x)))
  data <- data[remove == FALSE, , drop = FALSE]

  # compute the total num of samples and labeled data
  num_expert <- sum(data[, labeled], na.rm = TRUE)
  num_data   <- nrow(data)

  if(any(is.na(data[data[, labeled] == 0, predicted_var]) == FALSE)){
    stop("Some `predicted_var` in non-labeled data are not NA. Please check the data. ")
  }

  # Checking values
  if(all(sort(unique(data[, labeled])) %in% c(0,1)) == FALSE){
    stop(" `labeled` in `data` should be 0 or 1. Please check the data. ")
  }
  if(all(data[, sample_prob] > 0 & data[, sample_prob] <= 1) == FALSE){
    stop(" `sample_prob` in `data` should be greater than 0 and equal to or smaller than 1. Please check the data. ")
  }
  if(model %in% c("poisson", "fepois", "negbin", "fenegbin")){
    mf_check <- model.frame(formula, data = data, na.action = "na.pass")
    if(any(model.response(mf_check) < 0, na.rm = TRUE)){
      stop(paste0(" The outcome should be non-negative for `model = ", model, "`. Please check the data. "))
    }
    rm(mf_check)
  }

  # Discrete predicted covariates.
  # (count models) Factor, character, logical, or numeric 0/1 covariates in `predicted_var` are predicted as probabilities of their levels,
  # and the prediction part of the moments is the expectation over the predicted distribution (see `count_stack`).
  # (other models) Predicted covariates should be numeric.
  rhs_var <- all.vars(delete.response(terms_formula))
  pred_levels <- NULL
  for(v in intersect(predicted_var, rhs_var)){
    x_lab <- data[data[, labeled] == 1, v]
    if(is.factor(data[, v]) | is.character(data[, v]) | is.logical(data[, v])){
      if(count_model == FALSE){
        stop(paste0(" `", v, "` in `predicted_var` should be numeric for `model = ", model, "`. Please code a categorical variable as numeric 0/1 indicator variables. "))
      }
      if(is.factor(data[, v]) == TRUE){
        lev <- levels(data[, v])[levels(data[, v]) %in% as.character(x_lab)]
      }else{
        lev <- sort(unique(as.character(x_lab)))
      }
      if(length(lev) < 2){
        stop(paste0(" `", v, "` in `predicted_var` takes only one value in the labeled data. "))
      }
      data[, v] <- factor(as.character(data[, v]), levels = lev)
      pred_levels[[v]] <- lev
    }else if(count_model == TRUE & all(x_lab %in% c(0, 1)) & length(unique(x_lab)) == 2){
      pred_levels[[v]] <- c(0, 1)
    }
  }
  if(is.null(pred_levels) == FALSE){
    num_scenario <- prod(sapply(pred_levels, length))
    if(num_scenario > 100){
      stop(" Discrete predicted covariates have too many combinations of levels (more than 100). ")
    }
  }

  # Character predictors are converted to factors with levels in the whole data (so that every fold has the same columns)
  for(v in covariates_use){
    if(is.character(data[, v]) == TRUE){
      data[, v] <- factor(data[, v])
    }
  }

  # Standardize sample_prob
  mean_prob <- sum(data[, labeled], na.rm = TRUE)/nrow(data)
  data[, sample_prob] <- data[, sample_prob]*(mean_prob/mean(data[, sample_prob], na.rm = TRUE))


  # ##################
  # Data Preparation
  # ##################
  # Unique Cluster
  # Units of cross-fitting folds: clusters (observations without `cluster`).
  # (fepois, fenegbin with predicted covariates) The pairwise moment combines documents of the same group. Folds keep groups together
  # so that predictions do not depend on the labels of other documents in the same group.
  data$fold___ <- data$cluster___
  if(model %in% c("fepois", "fenegbin")){
    if(count_pairwise(formula, predicted_var) == TRUE){
      cluster_fold <- NULL
      if(clustered == TRUE){
        cluster_fold <- data$cluster___
      }
      index_fold <- count_fe_index(data, index, fixed_effect, TRUE, cluster_fold)
      data$fold___ <- count_components(data[, index_fold], cluster_fold)
    }
  }
  uniq_cluster <- sort(unique(data$fold___))
  if(length(uniq_cluster) < cross_fit){
    stop(paste0(" `cross_fit` (", cross_fit, ") should not be larger than the number of clusters (", length(uniq_cluster), ") used for cross-fitting. "))
  }

  # Number of Columns we analyze
  ncol_X <- ncol(model.matrix(formula, data = data))
  if(model == "felm"){
    ncol_X <- ncol_X - 1   # Remove Intercept

    ## including fixed effects
    if(fixed_effect != "twoways"){
      ncol_X_exp <- ncol_X + length(unique(data[, index_use]))
    }else if(fixed_effect == "twoways"){
      ncol_X_exp <- ncol_X + length(unique(data[, index[1]])) + length(unique(data[, index[2]])) - 1
    }
  }else if(model %in% c("fepois", "fenegbin")){
    ncol_X <- ncol_X - 1   # Remove Intercept

    ## (twoways) including fixed effects of the index with fewer levels
    if(fixed_effect == "twoways"){
      ncol_X_exp <- ncol_X + min(length(unique(data[, index[1]])), length(unique(data[, index[2]]))) - 1
    }else{
      ncol_X_exp <- ncol_X
    }
  }else{
    ncol_X_exp <- ncol_X
  }
  if(model %in% c("negbin", "fenegbin")){
    ncol_X_exp <- ncol_X_exp + 1 # including log(theta)
  }

  # matrices
  dm_point_l <- dm_se_l <- matrix(NA, ncol = ncol_X, nrow = sample_split)
  dm_vcov_l <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  dm_vcov_main_1_l <- dm_vcov_main_23_l <- array(NA, dim = c(ncol_X_exp, ncol_X_exp, sample_split))
  dm_vcov_J_l <- array(NA, dim = c(ncol_X_exp, ncol_X_exp, sample_split))
  dm_vcov_D_l <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  RMSE_cv <- matrix(NA, nrow = length(predicted_var), ncol = sample_split)
  dm_theta_l <- rep(NA, sample_split)
  split_fail <- rep(FALSE, sample_split)
  message_fail <- c()

  # ######################
  # Sample Splitting
  # ######################

  set.seed(seed)

  # Draw fold assignments for all sample splits before the loop.
  # fit_model() calls set.seed(seed), so drawing them inside the loop would give the same folds in every split.
  id_base0 <- rep(1:cross_fit, each = floor(length(uniq_cluster)/cross_fit))
  extra <- length(uniq_cluster) - length(id_base0)
  if(extra > 0){
    id_base0 <- c(id_base0, 1:extra)
  }
  id_base_all <- lapply(1:sample_split, function(x) sample(id_base0, size = length(id_base0), replace = FALSE))

  cat("Cross-Fitting: ")
  for(ss_use in 1:sample_split){
    cat(paste0(ss_use,"/", sample_split, ".."))

    id_base <- id_base_all[[ss_use]]

    dm_predicted_var <- list()
    dm_prob <- list()
    id_final <- c()
    RMSE_cv0 <- matrix(NA, nrow = length(predicted_var), ncol = cross_fit)
    # ############
    # Cross-Fit
    # ############
    for(i_use in 1:cross_fit){
      train_id <- uniq_cluster[id_base != i_use]
      test_id  <- setdiff(uniq_cluster, train_id)

      train_which <- unlist(sapply(train_id, function(x) which(data$fold___ == x)))
      test_which  <- unlist(sapply(test_id, function(x) which(data$fold___ == x)))

      data_train     <- data[train_which, , drop = FALSE]
      data_test      <- data[test_which, , drop = FALSE]

      # Fit the model E( predicted_var | covariates_use, R = 1) in the training sample

      for(v_use in 1:length(predicted_var)){
        v <- predicted_var[v_use]
        if(v %in% names(pred_levels)){
          # (count models) discrete predicted covariates: predict the probability of each level
          # (with two levels, only the probability of the second level)
          lev <- pred_levels[[v]]
          prob_test <- matrix(NA, nrow = nrow(data_test), ncol = length(lev))
          level_fit <- if(length(lev) == 2) 2 else seq_along(lev)
          RMSE_k <- c()
          for(k in level_fit){
            data_train$indicator___ <- as.numeric(data_train[, v] == lev[k])
            data_test$indicator___  <- as.numeric(data_test[, v] == lev[k])
            fit_train <- fit_model(outcome = "indicator___", labeled = labeled, prediction = prediction,
                                   covariates = covariates_use,
                                   data = data_train[data_train[, labeled] == 1, , drop = FALSE],
                                   seed = seed, method = sl_method, family = family)
            out_test <- fit_test(fit_out = fit_train,
                                 outcome = "indicator___", labeled = labeled, covariates = covariates_use,
                                 method = sl_method, family = family,
                                 data = data_test, seed = seed)
            prob_test[, k] <- pmin(pmax(out_test$Y_hat, 0), 1)
            RMSE_k <- c(RMSE_k, out_test$RMSE)
          }
          if(length(lev) == 2){
            prob_test[, 1] <- 1 - prob_test[, 2]
          }
          prob_sum <- rowSums(prob_test)
          prob_test[prob_sum == 0, ] <- 1
          prob_test <- prob_test/rowSums(prob_test)

          # placeholder for the column of the predicted covariate (the moments use prob_test)
          if(is.factor(data[, v]) == TRUE){
            predicted_var_test <- rep(lev[1], nrow(data_test))  # converted to a factor below
          }else{
            predicted_var_test <- prob_test[, 2]
          }
          if(i_use == 1){
            dm_prob[[v]] <- prob_test
          }else{
            dm_prob[[v]] <- rbind(prob_test, dm_prob[[v]])
          }
          RMSE_test <- mean(RMSE_k)
        }else{
          fit_train <- fit_model(outcome = v, labeled = labeled, prediction = prediction,
                                 covariates = covariates_use,
                                 data = data_train[data_train[, labeled] == 1, , drop = FALSE],
                                 seed = seed, method = sl_method, family = family)

          # Predict E(predicted_var | covariates_use, R = 1) in the test data
          out_test <- fit_test(fit_out = fit_train,
                               outcome = v, labeled = labeled, covariates = covariates_use,
                               method = sl_method, family = family,
                               data = data_test, seed = seed)

          # Predictions are kept within the range of the labeled outcome in the training data
          # (some learners, e.g., with the log link, extrapolate explosively)
          range_train <- range(data_train[data_train[, labeled] == 1, v])
          predicted_var_test <- pmin(pmax(out_test$Y_hat, range_train[1]), range_train[2])
          RMSE_test <- out_test$RMSE
        }

        # Store results
        if(i_use == 1){
          dm_predicted_var[[v_use]] <- predicted_var_test
        }else{
          dm_predicted_var[[v_use]] <- c(predicted_var_test, dm_predicted_var[[v_use]])
        }

        RMSE_cv0[v_use, i_use] <- RMSE_test
      }
      rm(data_train)
      rm(data_test)
      # id
      id_final <- c(test_which, id_final)
    }
    rownames(RMSE_cv0) <- predicted_var

    # data we use for estimation
    data_pred <- data_orig <- data[id_final, , drop = FALSE]

    for(v_use in 1:length(predicted_var)){
      if(is.factor(data[, predicted_var[v_use]]) == TRUE){
        data_pred[, predicted_var[v_use]] <- factor(dm_predicted_var[[v_use]], levels = levels(data[, predicted_var[v_use]]))
      }else{
        data_pred[, predicted_var[v_use]] <- dm_predicted_var[[v_use]]
      }
    }
    RMSE_cv[1:length(predicted_var), ss_use] <- apply(RMSE_cv0, 1, mean, na.rm = TRUE)

    # ##############
    # DSL Estimator
    # ##############
    # (count models) When the moment conditions cannot be solved in this sample split, we drop the split (see below)
    fit_dm <- tryCatch(dsl_general_moment_est(model = model,
                                     formula = formula,
                                     labeled = labeled,
                                     sample_prob = sample_prob,
                                     predicted_var = predicted_var,
                                     data_orig = data_orig,
                                     data_pred = data_pred,
                                     index = index,
                                     fixed_effect = fixed_effect,
                                     clustered = clustered,
                                     optim_method = optim_method,
                                     lambda = lambda,
                                     pred_levels = pred_levels,
                                     pred_prob = dm_prob),
                       dsl_convergence = function(e) e)
    if(inherits(fit_dm, "dsl_convergence") == TRUE){
      split_fail[ss_use] <- TRUE
      message_fail <- unique(c(message_fail, conditionMessage(fit_dm)))
      next
    }
    coef_name <- names(fit_dm$coefficients)
    cluster_note <- fit_dm$cluster_note

    # (fepois, fenegbin with predicted covariates and `cluster`) the index whose fixed effects are estimated as parameters
    # can be the one with more levels (see `dsl_general_moment_est`)
    if(ncol(fit_dm$J) != ncol_X_exp){
      ncol_X_exp <- ncol(fit_dm$J)
      dm_vcov_main_1_l <- dm_vcov_main_23_l <- dm_vcov_J_l <- array(NA, dim = c(ncol_X_exp, ncol_X_exp, sample_split))
    }

    dm_point_l[ss_use, 1:ncol_X] <- fit_dm$coefficients
    dm_se_l[ss_use, 1:ncol_X]    <- fit_dm$standard_errors
    dm_vcov_l[1:ncol_X, 1:ncol_X, ss_use] <- as.matrix(fit_dm$vcov)

    # dm_vcov_meat_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat)
    dm_vcov_main_1_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat_decomp$main_1)
    dm_vcov_main_23_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat_decomp$main_23)
    dm_vcov_J_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$J)
    dm_vcov_D_l[1:ncol_X, 1:ncol_X, ss_use] <- fit_dm$D
    dm_theta_l[ss_use] <- fit_dm$theta
  }
  rm(data)

  # Sample splits whose moment conditions could not be solved are dropped
  if(all(split_fail) == TRUE){
    stop(paste(message_fail, collapse = " | "), call. = FALSE)
  }
  if(any(split_fail) == TRUE){
    warning(paste0(sum(split_fail), " of ", sample_split, " sample splits were dropped because the estimation failed:", paste(message_fail, collapse = " | ")), call. = FALSE)
    keep <- (split_fail == FALSE)
    RMSE_cv <- RMSE_cv[, keep, drop = FALSE]
    dm_point_l <- dm_point_l[keep, , drop = FALSE]
    dm_se_l <- dm_se_l[keep, , drop = FALSE]
    dm_vcov_l <- dm_vcov_l[, , keep, drop = FALSE]
    dm_vcov_main_1_l <- dm_vcov_main_1_l[, , keep, drop = FALSE]
    dm_vcov_main_23_l <- dm_vcov_main_23_l[, , keep, drop = FALSE]
    dm_vcov_J_l <- dm_vcov_J_l[, , keep, drop = FALSE]
    dm_vcov_D_l <- dm_vcov_D_l[, , keep, drop = FALSE]
    dm_theta_l <- dm_theta_l[keep]
    sample_split <- sum(keep)
  }

  # ################################
  # Combine across sampling-splitting (Take into account the variation across sample splitting)
  # ################################

  # point estimates + standard errors
  dm_point <- apply(dm_point_l, 2, median)
  dm_se_0  <- dm_se_l^2 + (t(t(dm_point_l) - dm_point))^2
  dm_se    <- sqrt(apply(dm_se_0, 2, median))
  names(dm_point) <- names(dm_se) <- coef_name

  # vcov (Take into account the variation across sample splitting)
  dm_vcov_mat <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  for(ss_use in 1:sample_split){
    dm_point_use <- dm_point_l[ss_use, 1:ncol_X] - dm_point
    dm_vcov_adj  <- matrix(dm_point_use, nrow = ncol_X, ncol =1) %*% matrix(dm_point_use, nrow = 1, ncol = ncol_X)
    dm_vcov_mat[1:ncol_X, 1:ncol_X, ss_use] <- dm_vcov_l[1:ncol_X, 1:ncol_X, ss_use] + dm_vcov_adj
  }
  dm_vcov0 <- apply(dm_vcov_mat, MARGIN = 1:2, median)
  dm_vcov  <- (dm_vcov0 + t(dm_vcov0))/2 # make it symmetric
  if(is.positive.definite(dm_vcov) == FALSE){
    dm_vcov <- nearPD(dm_vcov, keepDiag = TRUE)
  }

  internal <- list("formula" = formula, "model" = model,
                   "predicted_var" = predicted_var, "prediction" = prediction,
                   "num_expert" = num_expert, "equal_prob" = equal_prob, "num_data" = num_data,
                   "cluster" = cluster, "cluster_note" = cluster_note,
                   "index" = index, "fixed_effect" = fixed_effect,
                   "sl_method"= sl_method, "feature" = feature, "cross_fit" = cross_fit, "sample_split" = sample_split,
                   "vcov_main_1" = dm_vcov_main_1_l,
                   "vcov_main_23" = dm_vcov_main_23_l,
                   "vcov_J" = dm_vcov_J_l,
                   "vcov_D" = dm_vcov_D_l)

  out <- list("coefficients" = dm_point,
              "standard_errors" = dm_se,
              "vcov" = dm_vcov,
              "RMSE" = RMSE_cv,
              "internal" = internal)
  if(model %in% c("negbin", "fenegbin")){
    out$theta <- median(dm_theta_l)
  }

  class(out)  <- c(class(out), "dsl")
  return(out)
}
