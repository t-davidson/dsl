dsl_general_moment_est <- function(model, formula, labeled, sample_prob, predicted_var, data_orig, data_pred, index = NULL, fixed_effect = NULL,
                                   clustered, optim_method = "L-BFGS-B",
                                   lambda = 0.00001, pred_levels = NULL, pred_prob = NULL){

  # Prepare data
  labeled_ind <- data_orig[, labeled]
  sample_prob_use <- data_orig[, sample_prob]
  count_model <- model %in% c("poisson", "fepois", "negbin", "fenegbin")

  mf_orig <- model.frame(formula, data = data_orig, na.action = "na.pass")
  X_orig <- model.matrix(formula, data = mf_orig)
  Y_orig <- as.numeric(model.response(mf_orig))
  mf_pred <- model.frame(formula, data = data_pred, na.action = "na.pass")
  X_pred <- model.matrix(formula, data = mf_pred)
  Y_pred <- as.numeric(model.response(mf_pred))

  # (count models) offset of the linear predictor (it cannot contain predicted variables)
  offset <- NULL
  if(count_model == TRUE){
    offset <- model.offset(mf_pred)
    if(is.null(offset) == TRUE){
      offset <- rep(0, nrow(X_pred))
    }
    offset <- as.numeric(offset)
  }
  rm(mf_orig)
  rm(mf_pred)
  col_name_keep <- colnames(X_orig)

  # (count models) Discrete predicted covariates. The prediction part of the moments is the expectation over the predicted
  # distribution of the covariates (see count_stack). X_scen[[k]] is the design matrix when the predicted covariates take
  # the k-th combination of their levels, and prob_scen[, k] is its predicted probability.
  # X_pred is the expected design matrix (used for standardization and starting values).
  X_scen <- NULL
  if(is.null(pred_levels) == FALSE){
    scen <- count_scenario(formula, data_pred, pred_levels, pred_prob)
    X_scen <- scen$X
    prob_scen <- scen$prob
    rm(scen)
    X_pred <- Reduce("+", lapply(seq_along(X_scen), function(k) X_scen[[k]] * prob_scen[, k]))
    # Columns that do not depend on the predicted covariates (e.g., the intercept) are copied exactly (sum_k prob_k = 1 only up to rounding)
    for(j in seq_len(ncol(X_pred))){
      if(all(sapply(X_scen, function(X) identical(X[, j], X_scen[[1]][, j]))) == TRUE){
        X_pred[, j] <- X_scen[[1]][, j]
      }
    }
  }


  # Standardize to make estimation easier
  X_orig_use <- X_pred_use <- matrix(1, nrow = nrow(X_orig), ncol = ncol(X_orig))
  scale_cov  <- which(apply(X_pred, 2, function(x) sd(x, na.rm = TRUE)) > 0)

  if(all(X_orig[,1] == 1) == TRUE){
    with_intercept <- TRUE
    mean_X <- rep(0, ncol(X_pred))
    sd_X   <- rep(1, ncol(X_pred))
    for(j in scale_cov){
      X_pred_use[, j] <- scale(X_pred[, j])
      X_orig_use[, j] <- (X_orig[, j] - mean(X_pred[, j], na.rm = TRUE))/sd(X_pred[, j], na.rm = TRUE)
      mean_X[j] <- mean(X_pred[, j], na.rm = TRUE)
      sd_X[j] <- sd(X_pred[, j], na.rm = TRUE)
    }
    mean_X <- mean_X[-1] # remove the first element
    sd_X   <- sd_X[-1] # remove the first element
  }else{
    with_intercept <- FALSE

    sd_X <- rep(1, ncol(X_pred))
    for(j in scale_cov){
      X_pred_use[, j] <- X_pred[, j]/sd(X_pred[, j], na.rm = TRUE)
      X_orig_use[, j] <- X_orig[, j]/sd(X_pred[, j], na.rm = TRUE)
      sd_X[j] <- sd(X_pred[, j], na.rm = TRUE)
    }
  }
  colnames(X_orig_use) <- colnames(X_orig)
  colnames(X_pred_use) <- colnames(X_pred)

  # (count models) standardize the design matrices of discrete predicted covariates in the same way, and stack them
  info <- NULL
  if(count_model == TRUE){
    info <- list("offset" = offset, "pred" = NULL)
    if(is.null(X_scen) == FALSE){
      X_scen_use <- lapply(X_scen, function(X){
        X_use <- matrix(1, nrow = nrow(X), ncol = ncol(X))
        for(j in scale_cov){
          if(with_intercept == TRUE){
            X_use[, j] <- (X[, j] - mean(X_pred[, j], na.rm = TRUE))/sd(X_pred[, j], na.rm = TRUE)
          }else{
            X_use[, j] <- X[, j]/sd(X_pred[, j], na.rm = TRUE)
          }
        }
        colnames(X_use) <- colnames(X)
        return(X_use)
      })
      if(model %in% c("fepois", "fenegbin")){
        X_scen_use <- lapply(X_scen_use, function(X) X[, -1, drop = FALSE])  # intercept is absorbed by fixed effects
      }
      info$pred <- list("X" = do.call("rbind", X_scen_use), "prob" = as.numeric(prob_scen),
                        "id" = rep(seq_len(nrow(X_pred)), length(X_scen_use)))
      rm(X_scen_use); rm(X_scen)
    }
  }

  # OLD
  # X_orig_use <- X_pred_use <- matrix(1, nrow = nrow(X_orig), ncol = ncol(X_orig))
  # # scale_cov  <- which(apply(X_pred, 2, function(x) sd(x, na.rm = TRUE)) > 0)
  # scale_cov  <- seq(from = 2, to = ncol(X_orig))
  # mean_X <- sd_X <- c()
  # for(j in scale_cov){
  #   X_pred_use[, j] <- scale(X_pred[, j])
  #   X_orig_use[, j] <- (X_orig[, j] - mean(X_pred[, j], na.rm = TRUE))/sd(X_pred[, j], na.rm = TRUE)
  #   mean_X[j-1] <- mean(X_pred[, j], na.rm = TRUE)
  #   sd_X[j-1] <- sd(X_pred[, j], na.rm = TRUE)
  # }
  # colnames(X_orig_use) <- colnames(X_orig)
  # colnames(X_pred_use) <- colnames(X_pred)

  rm(X_orig); rm(X_pred)

  # ##################################
  # Point Estimation
  # ##################################
  if(model == "felm"){
    # For Fixed Effects, we remove intercept
    X_orig_use <- X_orig_use[, -1, drop = FALSE]
    X_pred_use <- X_pred_use[, -1, drop = FALSE]

    # Fixed Effects
    adj_Y0 <- as.numeric(labeled_ind/sample_prob_use)*as.numeric(Y_pred - Y_orig)
    adj_Y0[labeled_ind == 0] <- 0
    adj_Y <- as.numeric(Y_pred) - adj_Y0
    adj_X0 <- (X_pred_use - X_orig_use) * as.numeric(labeled_ind/sample_prob_use)
    adj_X0[labeled_ind == 0, ] <- 0
    adj_X <- X_pred_use - adj_X0

    if(fixed_effect != "twoways"){
      fixed_effect_use <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[1]) ## This part is correct.

      fe_Y <- fixed_effect_use$adj_Y_avg_exp
      fe_X <- fixed_effect_use$adj_X_avg_exp
    }else if(fixed_effect == "twoways"){
      fixed_effect_use_1 <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[1])
      fixed_effect_use_2 <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[2])

      fe_Y <- fixed_effect_use_1$adj_Y_avg_exp + fixed_effect_use_2$adj_Y_avg_exp - mean(adj_Y, na.rm = TRUE)
      fe_X <- t(t(fixed_effect_use_1$adj_X_avg_exp + fixed_effect_use_2$adj_X_avg_exp) - apply(adj_X, 2, mean, na.rm = TRUE))
    }
  }else if(model %in% c("fepois", "fenegbin")){
    # For Fixed Effects, we remove intercept
    X_orig_use <- X_orig_use[, -1, drop = FALSE]
    X_pred_use <- X_pred_use[, -1, drop = FALSE]
    fe_Y <- fe_X <- NULL

    # (pairwise) When a covariate is predicted, we use the pairwise moment (see fepois)
    rhs_var  <- all.vars(delete.response(terms(formula)))
    pairwise <- any(predicted_var %in% rhs_var)

    # Fixed effects of the concentrated index are concentrated out (see fepois).
    # (twoways) Fixed effects of the other index are estimated with par (as kappa_t in fepois).
    # The concentrated index is the index with more levels. (pairwise) Standard errors are clustered at least at the level of
    # the concentrated index, so we prefer the index giving more independent clusters (e.g., the index nested in `cluster`).
    if(fixed_effect == "twoways"){
      num_level <- sapply(index, function(x) length(unique(data_orig[, x])))
      if(pairwise == TRUE & clustered == TRUE){
        num_comp  <- sapply(index, function(x) max(count_components(data_orig[, x], data_orig$cluster___)))
        index_use <- index[order(-num_comp, -num_level)[1]]
      }else{
        index_use <- index[which.max(num_level)]
      }
      index_dummy <- setdiff(index, index_use)
    }else{
      index_use   <- index[1]
      index_dummy <- NULL
    }
    fe_index <- match(data_orig[, index_use], unique(data_orig[, index_use])) # 1, ..., G (no unused levels)

    # Covariates should vary within fixed effects
    for(index_check in c(index_use, index_dummy)){
      no_var <- apply(X_pred_use, 2, function(x) all(tapply(x, data_orig[, index_check], function(z) max(z) - min(z)) == 0))
      if(any(no_var)){
        stop(paste0(" `", paste(colnames(X_pred_use)[no_var], collapse = "`, `"), "` does not vary within `", index_check,
                    "` and cannot be estimated with fixed effects of `", index_check, "`. "))
      }
    }

    fe_dummy <- kappa_ini <- NULL
    if(is.null(index_dummy) == FALSE){
      # The reference level (kappa_t = 1) is the level with the largest DSL-corrected sum of Y.
      # kappa_t starts at the ratio of the DSL-corrected means of Y (per unit of exp(offset)) to that of the reference level.
      Y_dr <- Y_pred + ifelse(labeled_ind == 1, Y_orig - Y_pred, 0) * as.numeric(labeled_ind/sample_prob_use)
      index_dummy_use <- factor(as.character(data_orig[, index_dummy])) # drop unused levels
      Y_dr_sum <- tapply(Y_dr, index_dummy_use, sum)
      rate <- Y_dr_sum/tapply(exp(offset), index_dummy_use, sum)
      index_dummy_use <- relevel(index_dummy_use, ref = names(Y_dr_sum)[which.max(Y_dr_sum)])
      fe_dummy <- model.matrix(~ index_dummy_use)[, -1, drop = FALSE]
      kappa_ini <- as.numeric(rate[levels(index_dummy_use)[-1]]/rate[levels(index_dummy_use)[1]])
      kappa_ini <- pmax(kappa_ini, 0.01)
    }

    # (pairwise) groups are weighted by 1/(number of observations), and standard errors are clustered by groups
    # (or by the smallest clusters that contain both groups and `cluster`)
    info$index <- fe_index
    info$dummy <- fe_dummy
    info$pairwise <- pairwise
    if(pairwise == TRUE){
      info$omega <- 1/as.numeric(table(fe_index))
      cluster_comp <- count_components(fe_index, if(clustered == TRUE) data_orig$cluster___ else NULL)
      if(max(cluster_comp) < 2){
        stop(paste0(" With predicted covariates, `", model, "` clusters standard errors at least at the level of `", index_use,
                    "`. Each cluster in `cluster` contains observations of all levels of `", index_use, "`, so standard errors cannot be estimated. "))
      }
    }
  }else{
    fe_Y <- fe_X <- NULL
  }

  inioptim <- rep(0, ncol(X_orig_use))
  if(model %in% c("fepois", "fenegbin")){
    inioptim <- c(inioptim, kappa_ini) # (twoways) kappa
  }

  if(count_model == TRUE){
    # (negbin, fenegbin) We first estimate the Poisson model, which gives the starting values
    model_pois <- ifelse(model %in% c("poisson", "negbin"), "poisson", "fepois")
    if(model_pois == "poisson"){
      # (poisson) Start the intercept at log of the DSL estimate of E(Y) (per unit of exp(offset)),
      # which solves the moment conditions when all slopes are zero
      Y_dr <- Y_pred + ifelse(labeled_ind == 1, Y_orig - Y_pred, 0) * as.numeric(labeled_ind/sample_prob_use)
      if(with_intercept == TRUE & mean(Y_dr) > 0){
        if(all(offset == 0)){
          inioptim[1] <- log(mean(Y_dr))
        }else{
          inioptim[1] <- log(sum(Y_dr)/sum(exp(offset)))
        }
      }
    }

    # (poisson, fepois) exp() makes `optim` unstable (overflow). We solve the same moment conditions with Newton-Raphson
    est0 <- dsl_general_newton(par = inioptim,
                               labeled_ind = labeled_ind,
                               sample_prob_use = sample_prob_use,
                               Y_orig = Y_orig,
                               X_orig = X_orig_use,
                               Y_pred = Y_pred,
                               X_pred = X_pred_use,
                               model = model_pois,
                               info = info)

    if(model %in% c("negbin", "fenegbin")){
      r_pi <- as.numeric(labeled_ind/sample_prob_use)
      Y_orig_0 <- ifelse(labeled_ind == 1, Y_orig, 0)
      if(model == "negbin"){
        X_orig_0 <- X_orig_use
        X_orig_0[labeled_ind == 0, ] <- 0
        mu_orig <- as.numeric(exp(X_orig_0 %*% est0 + offset))
      }else{
        # (fenegbin) mu_pilot = scale * exp(X par_X + offset) with the fepois estimates of par_X and one overall scale (see fenegbin)
        par_X_pilot <- est0[1:ncol(X_orig_use)]
        st <- count_stack(X_pred_use, Y_pred, info)
        X_orig_0 <- X_orig_use
        X_orig_0[labeled_ind == 0, ] <- 0
        exp_pilot_pred <- count_collapse(as.numeric(exp(st$X %*% par_X_pilot + st$offset)), st)
        exp_pilot_orig <- as.numeric(exp(X_orig_0 %*% par_X_pilot + offset))
        Y_dr <- Y_pred + (Y_orig_0 - Y_pred) * r_pi
        scale_pilot <- sum(Y_dr)/sum(exp_pilot_pred + (exp_pilot_orig - exp_pilot_pred) * r_pi)
        if(is.finite(scale_pilot) == FALSE | scale_pilot <= 0){
          dsl_stop_convergence(" `fenegbin` could not be estimated: the DSL-corrected sum of the outcome is not positive. ")
        }
        info$pilot <- list("par_X" = par_X_pilot, "scale" = scale_pilot)
        rm(st)

        mu_orig <- fenegbin_dsl_mu_orig(est0, labeled_ind, sample_prob_use, Y_orig, X_orig_use, Y_pred, X_pred_use, info)$mu_orig
        mu_orig <- pmax(mu_orig, 0)
      }
      # Starting value of theta: method of moments with labeled observations, E[(Y - mu)^2 - mu] = E[mu^2]/theta
      mu_orig[labeled_ind == 0] <- 0
      var_excess <- sum(r_pi * ((Y_orig_0 - mu_orig)^2 - mu_orig))
      mu_sq <- sum(r_pi * mu_orig^2)
      theta_ini <- ifelse(var_excess > 0, min(mu_sq/var_excess, 100), 100)

      # As in MASS::glm.nb, we alternate between the moments for the mean given theta and the moment for theta given the mean
      est0 <- dsl_negbin_solve(par = c(est0, log(theta_ini)),
                               labeled_ind = labeled_ind,
                               sample_prob_use = sample_prob_use,
                               Y_orig = Y_orig,
                               X_orig = X_orig_use,
                               Y_pred = Y_pred,
                               X_pred = X_pred_use,
                               model = model,
                               info = info)
    }
  }else{
    est0 <- optim(par = inioptim,
                  fn = dsl_general_moment,
                  labeled_ind = labeled_ind,
                  sample_prob_use = sample_prob_use,
                  Y_orig = Y_orig,
                  X_orig = X_orig_use,
                  Y_pred = Y_pred,
                  X_pred = X_pred_use,
                  fe_Y = fe_Y,
                  fe_X = fe_X,
                  lambda = lambda,
                  model = model,
                  method = optim_method,
                  hessian = FALSE,
                  control = list(maxit = 5000))$par
  }
  name_par <- colnames(X_orig_use)
  if(model %in% c("fepois", "fenegbin")){
    name_par <- c(name_par, colnames(info$dummy))
  }
  if(model %in% c("negbin", "fenegbin")){
    name_par <- c(name_par, "log(theta)")
    theta <- as.numeric(exp(est0[length(est0)]))
  }else{
    theta <- NA
  }
  names(est0) <- name_par

  # ###################################################
  # (felm) we compute estimates for other paramters
  # ###################################################
  if(model == "felm"){
    if(fixed_effect == "oneway"){
      index_use <- index[1]
      adj_data <- fixed_effect_use$adj_data
      est_fe0   <- as.numeric(adj_data[, "adj_Y_avg"]) - as.matrix(adj_data[, colnames(X_orig_use), drop = FALSE])%*% est0
      formula_use <- paste0("~factor(",index_use, ")")
      X_fe <- model.matrixBayes(as.formula(formula_use), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe) <- gsub(paste0("factor\\(",index_use, "\\)"), "", colnames(X_fe))
      est_fe <- est_fe0[match(colnames(X_fe), rownames(adj_data)), ]

      X_orig_use_exp <- cbind(X_orig_use, X_fe)
      X_pred_use_exp <- cbind(X_pred_use, X_fe)
      est0_exp <- c(est0, est_fe)
      rm(X_fe)
    }else if(fixed_effect == "twoways"){
      adj_data_1 <- fixed_effect_use_1$adj_data
      adj_data_2 <- fixed_effect_use_2$adj_data

      est_fe_base <- mean(adj_Y, na.rm = TRUE) - sum(est0* apply(adj_X, 2, mean, na.rm = TRUE))

      # alpha
      est_fe1_base   <- as.numeric(adj_data_1[, "adj_Y_avg"]) - as.matrix(adj_data_1[, colnames(X_orig_use), drop = FALSE])%*% est0 - est_fe_base
      formula_use1 <- paste0("~factor(",index[1], ")")
      X_fe1 <- model.matrixBayes(as.formula(formula_use1), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe1) <- gsub(paste0("factor\\(",index[1], "\\)"), "", colnames(X_fe1))
      est_fe1_base <- est_fe1_base[match(colnames(X_fe1), rownames(adj_data_1)),]

      # gamma
      est_fe2_base   <- as.numeric(adj_data_2[, "adj_Y_avg"]) - as.matrix(adj_data_2[, colnames(X_orig_use), drop = FALSE])%*% est0
      formula_use2 <- paste0("~factor(",index[2], ")")
      X_fe2 <- model.matrixBayes(as.formula(formula_use2), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe2) <- gsub(paste0("factor\\(",index[2], "\\)"), "", colnames(X_fe2))
      est_fe2_base <- est_fe2_base[match(colnames(X_fe2), rownames(adj_data_2)),]

      # we set gamma_0 = 0
      alpha_mean <- est_fe2_base[1]
      est_fe_1   <- est_fe1_base + alpha_mean
      est_fe_2   <- est_fe2_base - alpha_mean

      X_orig_use_exp <- cbind(X_orig_use, X_fe1, X_fe2[,-1, drop = FALSE])
      X_pred_use_exp <- cbind(X_pred_use, X_fe1, X_fe2[,-1, drop = FALSE])
      est0_exp <- c(est0, est_fe_1, est_fe_2[-1])
      rm(X_fe1); rm(X_fe2)
    }
  }else{
    est0_exp <- est0
    X_orig_use_exp <- X_orig_use
    X_pred_use_exp <- X_pred_use
    if(model %in% c("fepois", "negbin", "fenegbin")){
      est0 <- est0[1:ncol(X_orig_use)] # remove (twoways) fixed effects of the other index and (negbin, fenegbin) log(theta)
    }
  }
  rm(X_orig_use); rm(X_pred_use)

  # ################################
  # Variance Estimation
  # ################################
  n <- nrow(X_orig_use_exp)

  # Meat
  # (fepois, fenegbin with predicted covariates) clustered by groups of the concentrated index (or the smallest clusters containing them)
  clustered_meat <- clustered
  cluster_meat <- data_orig$cluster___
  if(isTRUE(info$pairwise)){
    clustered_meat <- TRUE
    cluster_meat <- cluster_comp
  }
  Meat_decomp <- dsl_general_moment_base_decomp(par = est0_exp,
                                                labeled_ind = labeled_ind,
                                                sample_prob_use = sample_prob_use,
                                                Y_orig = Y_orig,
                                                X_orig = X_orig_use_exp,
                                                Y_pred = Y_pred,
                                                X_pred = X_pred_use_exp,
                                                model  = model,
                                                clustered = clustered_meat,
                                                cluster = cluster_meat,
                                                info = info)
  # Meat   <- dsl_general_Meat(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model,
  #                            clustered, data_orig$cluster___)

  Meat <- Meat_decomp$main_1 + Meat_decomp$main_23

  # Jacobian
  J <- dsl_general_Jacobian(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model, info)

  # Variance
  s_J <- solve(J)
  s_J_use <- s_J[1:length(est0), , drop = FALSE]
  rm(s_J)

  V0 <- (s_J_use %*% Meat %*% t(s_J_use))/n

  # Variance of the target paramter
  if(model %in% c("felm", "fepois", "fenegbin")){ # with intercept
    # Rescale-back the coefficients and V
    D_2 <- cbind(diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X))) # we don't need intercept
    D   <- rbind(D_2)
  }else{
    # Rescale-back the coefficients and V
    if(with_intercept == TRUE){
      D_1 <- c(1, - mean_X/sd_X)
      D_2 <- cbind(0, diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X)))
      D   <- rbind(D_1, D_2)
    }else{
      D_2 <- cbind(diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X)))
      D   <- rbind(D_2)
    }
  }

  # Estimates and Variance
  est <- as.numeric(D %*% est0)
  V   <- D %*% V0 %*% t(D)

  if(model %in% c("felm", "fepois", "fenegbin")){
    col_name_keep <- col_name_keep[-1]
  }

  names(est) <- col_name_keep
  se <- sqrt(diag(V))
  names(se) <- col_name_keep

  # (fepois, fenegbin with predicted covariates) level of clustering when it differs from `cluster`
  cluster_note <- NULL
  if(isTRUE(info$pairwise)){
    if(clustered == FALSE){
      cluster_note <- paste0("`", index_use, "`")
    }else if(max(cluster_comp) < length(unique(data_orig$cluster___))){
      cluster_note <- paste0("the smallest clusters containing both `", index_use, "` and `cluster`")
    }
  }

  out <- list("coefficients" = est, "standard_errors" = se, "vcov" = V,
              "Meat" = Meat, "Meat_decomp" = Meat_decomp, "J" = J, "D" = D, "vcov0" = V0, "theta" = theta,
              "cluster_note" = cluster_note)
  return(out)
}
# ##############
# General
# ##############
dsl_general_moment <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X, model, lambda){

  if(model == "lm"){
    m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "logit"){
    m_dr   <- logit_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "felm"){
    m_dr   <- felm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X)
  }

  g <- apply(m_dr, 2, mean) # m_n(theta) in equation (6)
  g_out <- sum(g^2) + lambda*mean(par^2) # penalty

  return(g_out)
}

# ###############
# Newton-Raphson
# ###############
# Error signaled when the moment conditions of a count model cannot be solved in a sample split.
# dsl() drops such sample splits (with a warning), and stops only when every sample split fails.
dsl_stop_convergence <- function(message){
  cond <- structure(class = c("dsl_convergence", "error", "condition"), list(message = message, call = NULL))
  stop(cond)
}

dsl_convergence_message <- function(model, info){
  base <- paste0(" The DSL moment conditions of `", model, "` could not be solved (Newton-Raphson did not converge). ")
  if(model == "fepois" & isTRUE(info$pairwise)){
    reason <- paste0("With predicted covariates, this can happen when their predictions are poor. ",
                     "Including the outcome and the other covariates in `feature` usually improves the predictions. ")
  }else if(model == "fepois"){
    reason <- paste0("This can happen when the DSL-corrected sum of the outcome or of the fitted values is not positive in some groups ",
                     "(e.g., groups with few labeled observations and large prediction errors). ")
  }else{
    reason <- paste0("This can happen when predictions are poor or when some covariates perfectly predict zero outcomes (separation). ")
  }
  return(paste0(base, reason))
}

# Solves m_n(theta) = 0 by Newton-Raphson, using the analytical Jacobian (also used for the variance).
# Steps are halved until the objective in `dsl_general_moment`, sum(m_n(theta)^2), decreases.
# The Newton step is always a descent direction of this objective.
# (count models) info: see count_stack. (negbin, fepois) When theta is given, only the moments for the mean are solved given theta.
dsl_general_newton <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model,
                               info = NULL, theta = NULL, tol = 1e-8, maxit = 100){

  moment_mean <- function(par){
    if(model == "poisson"){
      m_dr <- poisson_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    }else if(model == "fepois"){
      m_dr <- fepois_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
    }else if(model == "negbin"){
      m_dr <- negbin_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
    }
    return(apply(m_dr, 2, mean)) # m_n(theta) in equation (6)
  }

  g   <- moment_mean(par)
  obj <- sum(g^2)
  if(is.finite(obj) == FALSE){
    dsl_stop_convergence(paste0(dsl_convergence_message(model, info), "(The moment conditions are not defined at the starting values.) "))
  }
  converged <- FALSE
  for(iter in 1:maxit){
    J <- dsl_general_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, info, theta)
    step <- tryCatch(as.numeric(solve(J, g)), error = function(e) NULL) # J is -d m_n(theta)/d theta
    if(is.null(step) == TRUE || all(is.finite(step)) == FALSE){
      break
    }
    if(max(abs(step)) < tol){
      converged <- TRUE
      break
    }

    # step-halving
    step_size <- 1
    improved  <- FALSE
    while(step_size > 1e-10){
      par_new <- par + step_size*step
      g_new   <- moment_mean(par_new)
      obj_new <- sum(g_new^2)
      if(is.finite(obj_new) & obj_new < obj){
        improved <- TRUE
        break
      }
      step_size <- step_size/2
    }
    if(improved == FALSE){
      # The objective cannot decrease further. This is a solution when the remaining step is negligible (numerical precision)
      converged <- max(abs(step)) < 1e-6
      break
    }
    par <- par_new; g <- g_new; obj <- obj_new
  }

  if(converged == FALSE){
    dsl_stop_convergence(dsl_convergence_message(model, info))
  }
  return(par)
}

# ##################################
# Estimation of negbin and fenegbin
# ##################################
# As in MASS::glm.nb, we alternate between (1) the moments for the mean parameters given theta (Newton-Raphson) and
# (2) the moment for log(theta) given the mean parameters (one-dimensional root finding, bracketed in theta in (1e-4, 1e6)),
# until both converge. The moment for log(theta) uses only labeled observations (see negbin).
# (negbin) par = (par_X, log(theta)). (fenegbin) par = (par_X, kappa, log(theta)); the moments for the mean are those of fepois
# with the negative binomial weights (see fenegbin).
dsl_negbin_solve <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, info, tol = 1e-8, maxit = 200){

  model_mean <- ifelse(model == "negbin", "negbin", "fepois")
  model_pois <- ifelse(model == "negbin", "poisson", "fepois")
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig_0 <- X_orig
  X_orig_0[labeled_ind == 0, ] <- 0
  Y_orig_0 <- ifelse(labeled_ind == 1, Y_orig, 0)

  # mu of labeled observations and their weights in the moment for log(theta)
  theta_parts <- function(par_mean){
    if(model == "negbin"){
      out <- list("mu_orig" = as.numeric(exp(X_orig_0 %*% par_mean + info$offset)), "w" = r_pi)
    }else{
      out <- fenegbin_dsl_mu_orig(par_mean, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    }
    return(out)
  }

  par_mean  <- par[-length(par)]
  log_theta <- par[length(par)]
  converged <- FALSE
  for(iter in 1:maxit){
    par_mean_new <- dsl_general_newton(par = par_mean, labeled_ind = labeled_ind, sample_prob_use = sample_prob_use,
                                       Y_orig = Y_orig, X_orig = X_orig, Y_pred = Y_pred, X_pred = X_pred,
                                       model = model_mean, info = info, theta = exp(log_theta))

    # log(theta) in (log(1e-4), log(1e6))
    th <- theta_parts(par_mean_new)
    moment_theta <- function(log_theta) negbin_dsl_theta_moment(log_theta, Y_orig_0, th$mu_orig, th$w)
    range_theta <- log(c(1e-4, 1e6))
    moment_range <- sapply(range_theta, moment_theta)
    if(any(is.finite(moment_range) == FALSE)){
      dsl_stop_convergence(paste0(" The moment condition of the dispersion parameter theta of `", model, "` is not defined. "))
    }
    if(moment_range[2] > 0){
      # The likelihood of the labeled observations increases up to theta = 1e6 (the Poisson limit)
      dsl_stop_convergence(paste0(" The dispersion parameter theta of `", model, "` diverges to infinity: the labeled observations are not overdispersed ",
                                  "relative to the fitted means, and the negative binomial model reduces to the Poisson model. Please use `model = ", model_pois, "`. "))
    }
    if(moment_range[1] < 0){
      dsl_stop_convergence(paste0(" The dispersion parameter theta of `", model, "` goes to zero: the labeled observations are far more dispersed ",
                                  "than the negative binomial model allows given the fitted means. Please use `model = ", model_pois, "`, which does not model the variance. "))
    }
    log_theta_new <- uniroot(moment_theta, interval = range_theta,
                             f.lower = moment_range[1], f.upper = moment_range[2], tol = 1e-10)$root

    change <- max(abs(c(par_mean_new - par_mean, log_theta_new - log_theta)))
    par_mean  <- par_mean_new
    log_theta <- log_theta_new
    if(change < 1e-7){
      converged <- TRUE
      break
    }
  }

  if(converged == FALSE){
    dsl_stop_convergence(paste0(" The estimation of `", model, "` did not converge (the mean parameters and theta keep changing). Please use `model = ", model_pois, "`. "))
  }
  return(c(par_mean, log_theta))
}

# (fepois, fenegbin with predicted covariates) Smallest clusters such that every group and every cluster in `cluster` is contained
# in one of them (connected components). Without `cluster`, the groups themselves. Returns 1, 2, ...
count_components <- function(group, cluster = NULL){
  comp <- match(group, unique(group))
  if(is.null(cluster) == FALSE){
    repeat{
      comp_new <- ave(comp, cluster, FUN = min)
      comp_new <- ave(comp_new, group, FUN = min)
      if(all(comp_new == comp)){
        break
      }
      comp <- comp_new
    }
    comp <- match(comp, unique(comp))
  }
  return(comp)
}

# (count models) Design matrices for every combination of the levels of the discrete predicted covariates and their probabilities.
# pred_levels: list of levels of each discrete predicted covariate. pred_prob: list of n x (number of levels) matrices of predicted probabilities.
# The probability of a combination is the product of the probabilities of its levels.
count_scenario <- function(formula, data_pred, pred_levels, pred_prob){
  grid <- expand.grid(lapply(pred_levels, seq_along), KEEP.OUT.ATTRS = FALSE)
  X_list <- list()
  prob <- matrix(1, nrow = nrow(data_pred), ncol = nrow(grid))
  for(k in 1:nrow(grid)){
    data_k <- data_pred
    for(v in names(pred_levels)){
      level_k <- pred_levels[[v]][grid[k, v]]
      if(is.factor(data_pred[, v]) == TRUE){
        data_k[, v] <- factor(rep(level_k, nrow(data_pred)), levels = levels(data_pred[, v]))
      }else{
        data_k[, v] <- rep(level_k, nrow(data_pred))
      }
      prob[, k] <- prob[, k] * pred_prob[[v]][, grid[k, v]]
    }
    mf_k <- model.frame(formula, data = data_k, na.action = "na.pass")
    X_list[[k]] <- model.matrix(formula, data = mf_k)
  }
  out <- list("X" = X_list, "prob" = prob)
  return(out)
}

# dsl_general_Meat <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, clustered, cluster){
#
#   n <- nrow(X_orig)
#
#   if(model == "lm"){
#     m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }else if(model == "logit"){
#     m_dr   <- logit_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }else if(model == "felm"){ # we can use the same function as lm.
#     m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }
#
#   if(clustered == FALSE){
#     Meat   <- (t(m_dr)%*% m_dr)/n
#   }else if(clustered == TRUE){
#     uniq_cluster <- sort(unique(cluster))
#     Meat_0 <- matrix(0, nrow = ncol(m_dr), ncol = ncol(m_dr))
#     for(z in 1:length(uniq_cluster)){
#       m_dr_cl <- m_dr[cluster == uniq_cluster[z], ,drop = FALSE]
#       matrix_ind <- matrix(1, nrow = nrow(m_dr_cl), ncol = nrow(m_dr_cl))
#       Meat_0 <- Meat_0 + t(m_dr_cl) %*% matrix_ind %*% m_dr_cl
#     }
#     Meat <- Meat_0/n
#   }
#
#   return(Meat)
# }

dsl_general_moment_base_decomp <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, clustered, cluster, info = NULL){

  if(model == "lm"){
    m_orig <- lm_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- lm_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "logit"){
    m_orig <- logit_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- logit_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "poisson"){
    m_orig <- poisson_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    m_pred <- poisson_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)

  }else if(model == "fepois"){
    m_orig <- fepois_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    m_pred <- fepois_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)

  }else if(model == "negbin"){
    m_orig <- negbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    m_pred <- negbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)

  }else if(model == "fenegbin"){
    m_orig <- fenegbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
    m_pred <- fenegbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)

  }else if(model == "felm"){ # we can use the same function as "lm"
    m_orig <- lm_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- lm_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

    m_orig <- Matrix(m_orig, sparse = TRUE)
    m_pred <- Matrix(m_pred, sparse = TRUE)
  }

  n <- nrow(X_orig)

  # Clustering
  if(clustered == FALSE){
    # diag_1 <- Matrix(diag(labeled_ind/(sample_prob_use^2), ncol = length(labeled_ind), nrow = length(labeled_ind)), sparse = TRUE)
    # diag_2 <- Matrix(diag(labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind)), sparse = TRUE)

    diag_1 <- Diagonal(x = labeled_ind/(sample_prob_use^2))
    diag_2 <- Diagonal(x = labeled_ind/sample_prob_use)

    main_1 <- (t(m_orig - m_pred) %*% diag_1 %*% (m_orig - m_pred))/n
    main_2 <- (t(m_pred)%*% m_pred)/n
    main_30 <- t(m_pred)%*% diag_2 %*% (m_orig - m_pred)
    main_3 <- (main_30 + t(main_30))/n
  }else if(clustered == TRUE){
    m_orig_m_pred <- (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
    if(model %in% c("poisson", "fepois", "negbin", "fenegbin")){
      # sums within clusters
      m_pred_sum <- rowsum(m_pred, cluster)
      m_orig_m_pred_sum <- rowsum(m_orig_m_pred, cluster)
    }else{
      uniq_cluster <- sort(unique(cluster))
      m_pred_sum <- matrix(NA, ncol = ncol(m_pred), nrow = length(uniq_cluster))
      m_orig_m_pred_sum <- matrix(NA, ncol = ncol(m_orig_m_pred), nrow = length(uniq_cluster))
      for(z in 1:length(uniq_cluster)){
        m_pred_sum[z,1:ncol(m_pred)] <- colSums(m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
        m_orig_m_pred_sum[z, 1:ncol(m_orig_m_pred)] <- colSums(m_orig_m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
      }
    }
    m_pred_sum <- Matrix(m_pred_sum, sparse = TRUE)
    m_orig_m_pred_sum <- Matrix(m_orig_m_pred_sum, sparse = TRUE)
    main_1  <- (t(m_orig_m_pred_sum)%*%m_orig_m_pred_sum)/n
    main_2  <- (t(m_pred_sum)%*%m_pred_sum)/n
    main_3_c0 <- (t(m_pred_sum)%*%m_orig_m_pred_sum)/n
    main_3  <- main_3_c0 + t(main_3_c0)
    # }else{
    # previous loop version
    #   uniq_cluster <- sort(unique(cluster))
    #   main_1_b <- main_2_b <- main_3_b <- matrix(0, nrow = ncol(m_pred), ncol = ncol(m_pred))
    #   for(z in 1:length(uniq_cluster)){
    #     m_orig_cl <- m_orig[cluster == uniq_cluster[z], ,drop = FALSE]
    #     m_pred_cl <- m_pred[cluster == uniq_cluster[z], ,drop = FALSE]
    #
    #     matrix_ind <- matrix(1, nrow = nrow(m_orig_cl), ncol = nrow(m_orig_cl))
    #
    #     m_orig_m_pred_cl <- (m_orig_cl - m_pred_cl) * as.numeric(labeled_ind/sample_prob_use)[cluster == uniq_cluster[z]]
    #
    #     if(model == "felm"){
    #       m_orig_m_pred_cl <- Matrix(m_orig_m_pred_cl, sparse = TRUE)
    #       matrix_ind <- Matrix(matrix_ind, sparse = TRUE)
    #     }
    #
    #     main_1_b <- main_1_b + (t(m_orig_m_pred_cl) %*% matrix_ind %*% (m_orig_m_pred_cl))
    #     main_2_b <- main_2_b + (t(m_pred_cl)%*% matrix_ind %*% m_pred_cl)
    #     main_3_b <- main_3_b + t(m_pred_cl)%*% matrix_ind %*% m_orig_m_pred_cl + t(m_orig_m_pred_cl)%*% matrix_ind %*% m_pred_cl
    #   }
    #   main_1 <- main_1_b/n
    #   main_2 <- main_2_b/n
    #   main_3 <- main_3_b/n
    # }
  }

  out   <- list("main_1" = main_1, "main_23" = main_2 + main_3)
  return(out)
}

dsl_general_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, info = NULL, theta = NULL){

  if(model == "lm"){
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }else if(model == "logit"){
    J   <- logit_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "poisson"){
    J   <- poisson_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  }else if(model == "fepois"){
    J   <- fepois_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  }else if(model == "negbin"){
    J   <- negbin_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  }else if(model == "fenegbin"){
    J   <- fenegbin_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  }else if(model == "felm"){
    # we can use the same function as "lm"
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }
  return(J)
}
