dsl_general_moment_est <- function(model, formula, labeled, sample_prob, predicted_var, data_orig, data_pred, index = NULL, fixed_effect = NULL,
                                   clustered, optim_method = "L-BFGS-B",
                                   lambda = 0.00001){

  # Prepare data
  labeled_ind <- data_orig[, labeled]
  sample_prob_use <- data_orig[, sample_prob]

  mf_orig <- model.frame(formula, data = data_orig, na.action = "na.pass")
  X_orig <- model.matrix(formula, data = mf_orig)
  Y_orig <- as.numeric(model.response(mf_orig))
  mf_pred <- model.frame(formula, data = data_pred, na.action = "na.pass")
  X_pred <- model.matrix(formula, data = mf_pred)
  Y_pred <- as.numeric(model.response(mf_pred))
  rm(mf_orig)
  rm(mf_pred)
  col_name_keep <- colnames(X_orig)


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

    # Fixed effects of the index with more levels are concentrated out (see `fepois_dsl_fe`).
    # (twoways) Fixed effects of the other index are estimated with par (as kappa_t in `fepois_dsl_fe`).
    if(fixed_effect == "twoways"){
      num_level   <- sapply(index, function(x) length(unique(data_orig[, x])))
      index_use   <- index[which.max(num_level)]
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

    fe_dummy <- NULL
    if(is.null(index_dummy) == FALSE){
      # The reference level (kappa_t = 1) is the level with the largest DSL-corrected sum of Y
      Y_dr <- Y_pred + ifelse(labeled_ind == 1, Y_orig - Y_pred, 0) * as.numeric(labeled_ind/sample_prob_use)
      index_dummy_use <- factor(as.character(data_orig[, index_dummy])) # drop unused levels
      Y_dr_sum <- tapply(Y_dr, index_dummy_use, sum)
      index_dummy_use <- relevel(index_dummy_use, ref = names(Y_dr_sum)[which.max(Y_dr_sum)])
      fe_dummy <- model.matrix(~ index_dummy_use)[, -1, drop = FALSE]
    }
    fe_info <- list("index" = fe_index, "dummy" = fe_dummy)
  }else{
    fe_Y <- fe_X <- NULL
  }

  inioptim <- rep(0, ncol(X_orig_use))
  if((model %in% c("fepois", "fenegbin")) == FALSE){
    fe_info <- NULL
  }else if(is.null(fe_info$dummy) == FALSE){
    inioptim <- c(inioptim, rep(1, ncol(fe_info$dummy))) # (fepois, fenegbin, twoways) kappa starts at 1
  }

  if(model %in% c("poisson", "fepois", "negbin", "fenegbin")){
    # (negbin, fenegbin) We first estimate the Poisson model, which gives the starting values
    model_pois <- ifelse(model %in% c("poisson", "negbin"), "poisson", "fepois")
    if(model_pois == "poisson"){
      # (poisson) Start the intercept at log of the DSL estimate of E(Y), which solves the moment conditions when all slopes are zero
      Y_dr <- Y_pred + ifelse(labeled_ind == 1, Y_orig - Y_pred, 0) * as.numeric(labeled_ind/sample_prob_use)
      if(with_intercept == TRUE & mean(Y_dr) > 0){
        inioptim[1] <- log(mean(Y_dr))
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
                               fe_info = fe_info)

    if(model %in% c("negbin", "fenegbin")){
      # Starting value of theta: method of moments with labeled observations, E[(Y - mu)^2 - mu] = E[mu^2]/theta
      if(model == "negbin"){
        mu_orig <- as.numeric(exp(X_orig_use %*% est0))
      }else{
        fe_pois <- fepois_dsl_fe(est0, labeled_ind, sample_prob_use, Y_orig, X_orig_use, Y_pred, X_pred_use, fe_info)
        mu_orig <- pmax(fe_pois$exp_fe * fe_pois$exp_orig, 0)

        # (fenegbin) mu_pilot for the weights: fepois estimates of (par_X, kappa) and fixed effects estimated only with predictions
        exp_fe_pilot <- as.numeric(rowsum(Y_pred, fe_info$index))/as.numeric(rowsum(fe_pois$exp_pred, fe_info$index))
        exp_fe_pilot[is.finite(exp_fe_pilot) == FALSE] <- 0
        fe_info$mu_pilot_pred <- pmax(exp_fe_pilot[fe_info$index] * fe_pois$exp_pred, 0)
        fe_info$mu_pilot_orig <- pmax(exp_fe_pilot[fe_info$index] * fe_pois$exp_orig, 0)
        rm(fe_pois)
      }
      r_pi <- as.numeric(labeled_ind/sample_prob_use)
      Y_orig_0 <- ifelse(labeled_ind == 1, Y_orig, 0)
      mu_orig[labeled_ind == 0] <- 0
      var_excess <- sum(r_pi * ((Y_orig_0 - mu_orig)^2 - mu_orig))
      mu_sq <- sum(r_pi * mu_orig^2)
      theta_ini <- ifelse(var_excess > 0, min(mu_sq/var_excess, 100), 100)

      if(model == "negbin"){
        est0 <- dsl_general_newton(par = c(est0, log(theta_ini)),
                                   labeled_ind = labeled_ind,
                                   sample_prob_use = sample_prob_use,
                                   Y_orig = Y_orig,
                                   X_orig = X_orig_use,
                                   Y_pred = Y_pred,
                                   X_pred = X_pred_use,
                                   model = model)
      }else{
        est0 <- dsl_fenegbin_solve(par = c(est0, log(theta_ini)),
                                   labeled_ind = labeled_ind,
                                   sample_prob_use = sample_prob_use,
                                   Y_orig = Y_orig,
                                   X_orig = X_orig_use,
                                   Y_pred = Y_pred,
                                   X_pred = X_pred_use,
                                   fe_info = fe_info)
      }
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
    name_par <- c(name_par, colnames(fe_info$dummy))
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
  Meat_decomp <- dsl_general_moment_base_decomp(par = est0_exp,
                                                labeled_ind = labeled_ind,
                                                sample_prob_use = sample_prob_use,
                                                Y_orig = Y_orig,
                                                X_orig = X_orig_use_exp,
                                                Y_pred = Y_pred,
                                                X_pred = X_pred_use_exp,
                                                model  = model,
                                                clustered = clustered,
                                                cluster = data_orig$cluster___,
                                                fe_info = fe_info)
  # Meat   <- dsl_general_Meat(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model,
  #                            clustered, data_orig$cluster___)

  Meat <- Meat_decomp$main_1 + Meat_decomp$main_23

  # Jacobian
  J <- dsl_general_Jacobian(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model, fe_info)

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

  out <- list("coefficients" = est, "standard_errors" = se, "vcov" = V,
              "Meat" = Meat, "Meat_decomp" = Meat_decomp, "J" = J, "D" = D, "vcov0" = V0, "theta" = theta)
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
  }else if(model == "poisson"){
    m_dr   <- poisson_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
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
# Solves m_n(theta) = 0 by Newton-Raphson, using the analytical Jacobian (also used for the variance).
# Steps are halved until the objective in `dsl_general_moment`, sum(m_n(theta)^2), decreases.
# The Newton step is always a descent direction of this objective.
dsl_general_newton <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model,
                               fe_info = NULL, theta = NULL, tol = 1e-8, maxit = 100){

  moment_mean <- function(par){
    if(model == "poisson"){
      m_dr <- poisson_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    }else if(model == "fepois"){
      m_dr <- fepois_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
    }else if(model == "negbin"){
      m_dr <- negbin_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    }
    return(apply(m_dr, 2, mean)) # m_n(theta) in equation (6)
  }

  g   <- moment_mean(par)
  obj <- sum(g^2)
  if(is.finite(obj) == FALSE){
    stop(" The moment conditions are not defined at the starting values. ")
  }
  converged <- FALSE
  for(iter in 1:maxit){
    J <- dsl_general_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, fe_info, theta)
    step <- as.numeric(solve(J, g)) # J is -d m_n(theta)/d theta
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
      break
    }
    par <- par_new; g <- g_new; obj <- obj_new

    # (negbin) the last element of par is log(theta)
    if(model == "negbin"){
      if(exp(par[length(par)]) > 1e6){
        stop(" The dispersion parameter theta diverges to infinity: the outcome does not appear to be overdispersed, and the negative binomial model reduces to the Poisson model. Please use `model = poisson`. ")
      }
    }
  }

  if(converged == FALSE){
    warning(" Newton-Raphson did not converge. Estimates may be unreliable. ")
  }
  return(par)
}

# ##########################
# Estimation of fenegbin
# ##########################
# As in MASS::glm.nb, we alternate between (1) the moments for (par_X, kappa) given theta (fepois with negative binomial weights; Newton-Raphson)
# and (2) the moment for log(theta) given (par_X, kappa) (one-dimensional root finding), until both converge.
dsl_fenegbin_solve <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, tol = 1e-8, maxit = 100){

  moment_theta <- function(log_theta, par_mean){
    m_dr <- fenegbin_dsl_moment_base(c(par_mean, log_theta), labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
    return(mean(m_dr[, ncol(m_dr)]))
  }

  message_fail <- paste0(" `fenegbin` could not be estimated: the dispersion parameter theta is not identified. This happens when there are ",
                         "too few labeled observations per fixed effect. Please use `model = fepois`, which is valid for overdispersed outcomes. ")

  par_mean  <- par[-length(par)]
  log_theta <- par[length(par)]
  converged <- FALSE
  for(iter in 1:maxit){
    par_mean_new <- tryCatch(dsl_general_newton(par = par_mean, labeled_ind = labeled_ind, sample_prob_use = sample_prob_use,
                                                Y_orig = Y_orig, X_orig = X_orig, Y_pred = Y_pred, X_pred = X_pred,
                                                model = "fepois", fe_info = fe_info, theta = exp(log_theta)),
                             error = function(e) stop(message_fail))

    # log(theta) in (log(1e-4), log(1e6))
    range_theta <- log(c(1e-4, 1e6))
    moment_range <- sapply(range_theta, moment_theta, par_mean = par_mean_new)
    if(any(is.finite(moment_range) == FALSE)){
      stop(message_fail)
    }
    if(moment_range[2] > 0){
      stop(" The dispersion parameter theta diverges to infinity: the outcome does not appear to be overdispersed, and the negative binomial model reduces to the Poisson model. Please use `model = fepois`. ")
    }
    if(moment_range[1] < 0){
      stop(message_fail)
    }
    log_theta_new <- uniroot(moment_theta, interval = range_theta, par_mean = par_mean_new,
                             f.lower = moment_range[1], f.upper = moment_range[2], tol = tol)$root

    change <- max(abs(c(par_mean_new - par_mean, log_theta_new - log_theta)))
    par_mean  <- par_mean_new
    log_theta <- log_theta_new
    if(change < 1e-6){
      converged <- TRUE
      break
    }
  }

  if(converged == FALSE){
    warning(" The estimation of fenegbin did not converge. Estimates may be unreliable. ")
  }
  return(c(par_mean, log_theta))
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

dsl_general_moment_base_decomp <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, clustered, cluster, fe_info = NULL){

  if(model == "lm"){
    m_orig <- lm_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- lm_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "logit"){
    m_orig <- logit_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- logit_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "poisson"){
    m_orig <- poisson_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- poisson_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "fepois"){
    m_orig <- fepois_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
    m_pred <- fepois_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)

  }else if(model == "negbin"){
    m_orig <- negbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
    m_pred <- negbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)

  }else if(model == "fenegbin"){
    m_orig <- fenegbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
    m_pred <- fenegbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)

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
    uniq_cluster <- sort(unique(cluster))
    m_orig_m_pred <- (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
    m_pred_sum <- matrix(NA, ncol = ncol(m_pred), nrow = length(uniq_cluster))
    m_orig_m_pred_sum <- matrix(NA, ncol = ncol(m_orig_m_pred), nrow = length(uniq_cluster))
    for(z in 1:length(uniq_cluster)){
      m_pred_sum[z,1:ncol(m_pred)] <- colSums(m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
      m_orig_m_pred_sum[z, 1:ncol(m_orig_m_pred)] <- colSums(m_orig_m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
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

dsl_general_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, fe_info = NULL, theta = NULL){

  if(model == "lm"){
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }else if(model == "logit"){
    J   <- logit_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "poisson"){
    J   <- poisson_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "fepois"){
    J   <- fepois_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  }else if(model == "negbin"){
    J   <- negbin_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "fenegbin"){
    J   <- fenegbin_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  }else if(model == "felm"){
    # we can use the same function as "lm"
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }
  return(J)
}
