# ######
# lm
# ######
lm_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_orig <- X_orig * as.numeric(Y_orig - X_orig %*% par)

  m_orig[labeled_ind == 0, ] <- 0

  m_pred <- X_pred * as.numeric(Y_pred - X_pred %*% par)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

lm_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_orig <- X_orig * as.numeric(Y_orig - X_orig %*% par)
  m_orig[labeled_ind == 0, ] <- 0
  return(m_orig)
}

lm_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_pred <- X_pred * as.numeric(Y_pred - X_pred %*% par)
  return(m_pred)
}

lm_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model){
  X_orig[labeled_ind == 0, ] <- 0

  if(model == "felm"){
    X_orig <- Matrix(X_orig, sparse = TRUE)
    X_pred <- Matrix(X_pred, sparse = TRUE)
  }
  diag_1 <- Diagonal(x = labeled_ind/sample_prob_use)
  diag_2 <- Diagonal(x = 1 - labeled_ind/sample_prob_use)

  # diag_1 <- diag(labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind))
  # diag_2 <- diag(1 - labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind))

  J <- (t(X_orig) %*% diag_1 %*% X_orig + t(X_pred) %*% diag_2 %*% X_pred)/nrow(X_orig)
  return(J)
}

# #####
# felm
# #####
felm_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X){
  ## We just need to compute the alpha and gamma first.
  fe_use <- as.numeric(fe_Y) - as.matrix(fe_X) %*% as.numeric(par)

  m_orig <- X_orig * as.numeric(Y_orig - fe_use - X_orig %*% par)

  m_orig[labeled_ind == 0, ] <- 0

  m_pred <- X_pred * as.numeric(Y_pred - fe_use - X_pred %*% par)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

demean_dsl <- function(data_base, adj_Y, adj_X, index_use){
  data_base$uniq_id___ <- seq(1:nrow(data_base))
  adj_Y_avg  <- tapply(adj_Y, data_base[, index_use], mean, na.rm = TRUE)
  adj_X_avg0 <- lapply(1:ncol(adj_X), FUN = function(x) tapply(adj_X[, x], data_base[, index_use], mean, na.rm = TRUE))
  adj_X_avg  <- do.call("cbind", adj_X_avg0)
  colnames(adj_X_avg) <- colnames(adj_X)
  adj_data   <- data.frame(cbind(adj_Y_avg, adj_X_avg))
  adj_data$index___ <- rownames(adj_data)
  fixed_effect_use <- merge(data_base[, c(index_use, "uniq_id___"), drop = FALSE], adj_data, by.x = index_use, by.y = "index___", all.x = TRUE, sort = FALSE)
  fixed_effect_use <- fixed_effect_use[order(fixed_effect_use$uniq_id___, decreasing = FALSE), , drop = FALSE]
  rm(data_base)
  out <- list("adj_Y_avg_exp" = fixed_effect_use[, "adj_Y_avg"], "adj_X_avg_exp" = fixed_effect_use[, colnames(adj_X), drop = FALSE],
              "adj_data" = adj_data)
  return(out)
}

# #####
# logit
# #####
logit_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit <- 1/(1+exp(-X_orig %*% par))
  m_orig <- X_orig * as.numeric(Y_orig - inv_logit)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  inv_logit_pred <- 1/(1+exp(-X_pred %*% par))
  m_pred <- X_pred * as.numeric(Y_pred - inv_logit_pred)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

logit_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit <- 1/(1+exp(-X_orig %*% par))
  m_orig <- X_orig * as.numeric(Y_orig - inv_logit)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

logit_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit_pred <- 1/(1+exp(-X_pred %*% par))
  m_pred <- X_pred * as.numeric(Y_pred - inv_logit_pred)
  return(m_pred)
}

logit_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){

  #
  # inv_logit_pred <- 1/(1+exp(-X_pred%*%par))
  # inv_logit_pred2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  inv_logit_pred <- exp(-X_pred%*%par)
  inv_logit_pred2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  diag_pred2 <- Diagonal(x = inv_logit_pred2)
  diag_pred2_R <- Diagonal(x = inv_logit_pred2*labeled_ind/sample_prob_use)

  grad_pred <- (t(X_pred) %*% diag_pred2 %*% X_pred)/nrow(X_pred)
  grad_pred_R <- (t(X_pred) %*% diag_pred2_R %*% X_pred)/nrow(X_pred)

  # original
  # inv_logit_orig  <- 1/(1+exp(-X_orig%*%par))
  # # inv_logit_orig2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  # inv_logit_orig2 <- inv_logit_orig/(1 + inv_logit_orig)^2
  inv_logit_orig  <- exp(-X_orig%*%par)
  inv_logit_orig2 <- inv_logit_orig/(1 + inv_logit_orig)^2
  inv_logit_orig2[labeled_ind == 0] <- 0  # r/pi * Y
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  diag_orig2_R <- Diagonal(x = inv_logit_orig2 *labeled_ind/sample_prob_use)
  grad_orig <- (t(X_orig) %*% diag_orig2_R %*% X_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig - grad_pred_R

  out <- grad_main
  return(out)
}

# #######
# poisson
# #######
# Poisson regression with the log link. The moment is the Poisson score X * (Y - exp(X par)).
# As the variance is estimated with the sandwich formula, this is Poisson Pseudo Maximum Likelihood (PPML):
# it only requires E(Y | X) = exp(X par), and Y does not need to be an integer or follow a Poisson distribution.
poisson_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  mu_orig <- exp(X_orig %*% par)
  m_orig <- X_orig * as.numeric(Y_orig - mu_orig)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  mu_pred <- exp(X_pred %*% par)
  m_pred <- X_pred * as.numeric(Y_pred - mu_pred)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

poisson_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  mu_orig <- exp(X_orig %*% par)
  m_orig <- X_orig * as.numeric(Y_orig - mu_orig)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

poisson_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  mu_pred <- exp(X_pred %*% par)
  m_pred <- X_pred * as.numeric(Y_pred - mu_pred)
  return(m_pred)
}

poisson_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){

  # prediction (d/d par of exp(X par) is exp(X par) * X)
  mu_pred <- as.numeric(exp(X_pred %*% par))
  diag_pred2 <- Diagonal(x = mu_pred)
  diag_pred2_R <- Diagonal(x = mu_pred*labeled_ind/sample_prob_use)

  grad_pred <- (t(X_pred) %*% diag_pred2 %*% X_pred)/nrow(X_pred)
  grad_pred_R <- (t(X_pred) %*% diag_pred2_R %*% X_pred)/nrow(X_pred)

  # original
  mu_orig <- as.numeric(exp(X_orig %*% par))
  mu_orig[labeled_ind == 0] <- 0  # r/pi * Y
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  diag_orig2_R <- Diagonal(x = mu_orig*labeled_ind/sample_prob_use)
  grad_orig <- (t(X_orig) %*% diag_orig2_R %*% X_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig - grad_pred_R

  out <- grad_main
  return(out)
}

# ######
# fepois
# ######
# Poisson regression with fixed effects (PPML with fixed effects). The conditional mean is exp_fe_g * kappa_t * exp(X par_X),
# where exp_fe_g = exp(alpha_g) is the fixed effect of group g of the index with more levels, and
# kappa_t = exp(gamma_t) is the fixed effect of the other index (twoways only; kappa_t = 1 for the reference level).
# exp_fe_g is concentrated out. For a given par, exp_fe_g solves the DSL moment of the fixed effect,
#   sum_{i in g} [Y_dr_i - exp_fe_g * w_dr_i] = 0,
# where Y_dr and w_dr are the DSL versions of Y and kappa_t * exp(X par_X). Thus, exp_fe_g = sum_g Y_dr / sum_g w_dr.
# Fixed effects are estimated on the linear scale (exp_fe_g and kappa_t rather than alpha_g and gamma_t). They exist even when
# the DSL-corrected sum of Y within a group is negative (e.g., a small group with a labeled observation with a large prediction error).
# Substituting exp_fe_g gives the moments for par = (par_X, kappa) with the design Z = [X, dummy variables of the other index]
# centered at center_g = sum_g DSL(Z * kappa_t * exp(X par_X)) / sum_g w_dr. This accounts for the estimation of exp_fe_g (as demeaning in felm).
# (fenegbin) When theta is given, every term is multiplied by the negative binomial weight theta/(theta + mu_pilot) (see fenegbin).
fepois_dsl_fe <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta = NULL){
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0

  par_X <- par[1:ncol(X_pred)]
  exp_pred_X <- as.numeric(exp(X_pred %*% par_X))
  exp_orig_X <- as.numeric(exp(X_orig %*% par_X))
  exp_orig_X[labeled_ind == 0] <- 0

  if(is.null(fe_info$dummy) == TRUE){
    kappa  <- 1
    Z_pred <- X_pred
    Z_orig <- X_orig
  }else{
    kappa  <- as.numeric(1 + fe_info$dummy %*% (par[-(1:ncol(X_pred))] - 1))
    Z_pred <- cbind(X_pred, fe_info$dummy)
    Z_orig <- cbind(X_orig, fe_info$dummy)
  }
  exp_pred <- kappa * exp_pred_X
  exp_orig <- kappa * exp_orig_X

  # derivative of kappa_t * exp(X par_X) with respect to par
  if(is.null(fe_info$dummy) == TRUE){
    d_pred <- X_pred * exp_pred
    d_orig <- X_orig * exp_orig
  }else{
    d_pred <- cbind(X_pred * exp_pred, fe_info$dummy * exp_pred_X)
    d_orig <- cbind(X_orig * exp_orig, fe_info$dummy * exp_orig_X)
  }

  # (fenegbin) weights theta/(theta + mu_pilot). mu_pilot is fixed during estimation (see fenegbin)
  if(is.null(theta) == TRUE){
    weight_pred <- weight_orig <- 1
  }else{
    weight_pred <- theta/(theta + fe_info$mu_pilot_pred)
    weight_orig <- theta/(theta + fe_info$mu_pilot_orig)
  }

  Y_dr <- Y_pred * weight_pred + (Y_orig * weight_orig - Y_pred * weight_pred) * r_pi
  w_dr <- exp_pred * weight_pred + (exp_orig * weight_orig - exp_pred * weight_pred) * r_pi
  v_dr <- Z_pred * (exp_pred * weight_pred) + (Z_orig * (exp_orig * weight_orig) - Z_pred * (exp_pred * weight_pred)) * r_pi

  # fe_info$index takes values 1, ..., G. So, the g-th row of rowsum() corresponds to group g
  w_sum <- as.numeric(rowsum(w_dr, fe_info$index))
  w_sum[w_sum <= 0] <- NA  # exp_fe_g is not defined
  exp_fe <- (as.numeric(rowsum(Y_dr, fe_info$index))/w_sum)[fe_info$index]
  center <- (rowsum(v_dr, fe_info$index)/w_sum)[fe_info$index, , drop = FALSE]

  out <- list("exp_fe" = exp_fe, "center" = center, "exp_orig" = exp_orig, "exp_pred" = exp_pred,
              "Z_orig" = Z_orig, "Z_pred" = Z_pred, "d_orig" = d_orig, "d_pred" = d_pred,
              "weight_orig" = weight_orig, "weight_pred" = weight_pred)
  return(out)
}

fepois_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta = NULL){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)

  m_orig <- (fe$Z_orig - fe$center) * as.numeric(fe$weight_orig * (Y_orig - fe$exp_fe * fe$exp_orig))
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  m_pred <- (fe$Z_pred - fe$center) * as.numeric(fe$weight_pred * (Y_pred - fe$exp_fe * fe$exp_pred))
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

fepois_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta = NULL){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  m_orig <- (fe$Z_orig - fe$center) * as.numeric(fe$weight_orig * (Y_orig - fe$exp_fe * fe$exp_orig))
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

fepois_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta = NULL){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  m_pred <- (fe$Z_pred - fe$center) * as.numeric(fe$weight_pred * (Y_pred - fe$exp_fe * fe$exp_pred))
  return(m_pred)
}

fepois_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta = NULL){
  # Jacobian of the concentrated moment. Centering Z at center_g accounts for the dependence of exp_fe_g on par.
  # (fenegbin) Weights are treated as fixed; their derivatives multiply (Y - mu) and have mean zero.
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)

  # prediction
  diag_pred2 <- Diagonal(x = fe$exp_fe * fe$weight_pred * (1 - labeled_ind/sample_prob_use))
  grad_pred <- (t(fe$Z_pred - fe$center) %*% diag_pred2 %*% fe$d_pred)/nrow(X_pred)

  # original
  diag_orig2_R <- Diagonal(x = fe$exp_fe * fe$weight_orig * labeled_ind/sample_prob_use)
  grad_orig <- (t(fe$Z_orig - fe$center) %*% diag_orig2_R %*% fe$d_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig

  out <- grad_main
  return(out)
}

# ######
# negbin
# ######
# Negative binomial (NB2) regression: Var(Y | X) = mu + mu^2/theta with mu = exp(X par_X). par = (par_X, log(theta)).
# Moments are the maximum likelihood scores:
#   (1) mean:       X * (Y - mu) * theta/(theta + mu)  (Poisson score weighted by theta/(theta + mu))
#   (2) dispersion: d loglik / d log(theta)
# The score for the mean is valid when E(Y | X) = mu even if Y does not follow the negative binomial distribution.
# Predictions carry no information about the dispersion of Y (Y_pred is an estimate of the conditional mean), so the score for log(theta)
# uses only labeled observations weighted by 1/pi (m_pred = 0 in the DSL moment). theta only affects the weights and not the validity of the coefficients.

# Scores for one observation given Y, mu, and theta (and their derivatives)
negbin_dsl_parts <- function(Y, mu, theta){
  mu_pos <- pmax(mu, 0)  # (fenegbin) fixed effects on the linear scale can make mu negative. Then, weights are those of Poisson.
  Y_pos  <- pmax(Y, 0)   # digamma(Y + theta) requires Y + theta > 0 (only matters when predictions are negative)

  resid_w <- (Y - mu) * theta/(theta + mu_pos)                              # moment for the mean
  a <- ifelse(mu > 0, theta * (theta + Y)/(theta + mu_pos)^2, 1)             # - d resid_w / d mu
  resid_w_phi <- (Y - mu) * theta * mu_pos/(theta + mu_pos)^2                # d resid_w / d log(theta)

  s_theta <- digamma(Y_pos + theta) - digamma(theta) + log(theta/(theta + mu_pos)) + (mu_pos - Y_pos)/(theta + mu_pos)
  s_theta_theta <- trigamma(Y_pos + theta) - trigamma(theta) + 1/theta - 1/(theta + mu_pos) - (mu_pos - Y_pos)/(theta + mu_pos)^2
  s_phi <- theta * s_theta                                                   # moment for log(theta)
  s_phi_mu <- theta * (Y_pos - mu_pos)/(theta + mu_pos)^2 * as.numeric(mu > 0)  # d s_phi / d mu
  s_phi_phi <- theta * s_theta + theta^2 * s_theta_theta                     # d s_phi / d log(theta)

  out <- list("resid_w" = resid_w, "a" = a, "resid_w_phi" = resid_w_phi,
              "s_phi" = s_phi, "s_phi_mu" = s_phi_mu, "s_phi_phi" = s_phi_phi)
  return(out)
}

# Moments (orig and pred) and Jacobian (J = - d m_n/d par)
negbin_dsl_block <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  n <- nrow(X_pred)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0

  theta <- exp(par[length(par)])
  par_X <- par[1:ncol(X_pred)]
  mu_pred <- as.numeric(exp(X_pred %*% par_X))
  mu_orig <- as.numeric(exp(X_orig %*% par_X))

  s_pred <- negbin_dsl_parts(Y_pred, mu_pred, theta)
  s_orig <- negbin_dsl_parts(Y_orig, mu_orig, theta)

  m_pred <- cbind(X_pred * s_pred$resid_w, 0)  # the score for log(theta) uses only labeled observations
  m_orig <- cbind(X_orig * s_orig$resid_w, s_orig$s_phi)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  # m_dr = (1 - r/pi) * m_pred + r/pi * m_orig
  w_pred <- 1 - r_pi
  w_orig <- r_pi
  J_1 <- t(X_pred) %*% (X_pred * (s_pred$a * mu_pred * w_pred)) + t(X_orig) %*% (X_orig * (s_orig$a * mu_orig * w_orig))
  J_2 <- - colSums(X_pred * (s_pred$resid_w_phi * w_pred)) - colSums(X_orig * (s_orig$resid_w_phi * w_orig))
  J_3 <- - colSums(X_orig * (s_orig$s_phi_mu * mu_orig * w_orig))
  J_4 <- - sum(s_orig$s_phi_phi * w_orig)
  J <- rbind(cbind(J_1, J_2), c(J_3, J_4))/n

  out <- list("m_pred" = m_pred, "m_orig" = m_orig, "J" = J, "mu_pred" = mu_pred, "mu_orig" = mu_orig)
  return(out)
}

negbin_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  blk <- negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  m_dr <- blk$m_pred + (blk$m_orig - blk$m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

negbin_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)$m_orig)
}

negbin_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)$m_pred)
}

negbin_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)$J)
}

# ########
# fenegbin
# ########
# Negative binomial (NB2) regression with fixed effects. par = (par_X, kappa, log(theta)) (see fepois for kappa).
# The moments for (par_X, kappa) are those of fepois weighted by the negative binomial weight theta/(theta + mu_pilot), and
# the moment for log(theta) is the maximum likelihood score (see negbin) with mu = exp_fe_g * kappa_t * exp(X par_X).
# Maximum likelihood uses the weights theta/(theta + mu) with the estimated fixed effect exp_fe_g. With DSL, exp_fe_g depends heavily on
# the few labeled observations of the group (weighted by 1/pi), and weights that depend on them bias the coefficients
# (incidental parameter problem). Thus, mu_pilot = exp_fe_pilot_g * kappa_t * exp(X par_X) uses (par_X, kappa) from the fepois estimates and
# the fixed effect estimated only with predictions (exp_fe_pilot_g = sum_g Y_pred / sum_g kappa_t * exp(X par_X)), and it is fixed during estimation.
# The moments for (par_X, kappa) are then linear in Y, as in fepois, and remain valid for any weights.
# As in negbin, the score for log(theta) uses only labeled observations weighted by 1/pi.
# Cross-derivatives between the mean parameters and theta multiply (Y - mu) and have mean zero, so the Jacobian omits them.
fenegbin_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  fe <- fepois_dsl_fe(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  m_orig <- (fe$Z_orig - fe$center) * as.numeric(fe$weight_orig * (Y_orig - fe$exp_fe * fe$exp_orig))
  s_orig <- negbin_dsl_parts(Y_orig, fe$exp_fe * fe$exp_orig, theta)
  m_orig <- cbind(m_orig, s_orig$s_phi)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

fenegbin_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  fe <- fepois_dsl_fe(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  m_pred <- (fe$Z_pred - fe$center) * as.numeric(fe$weight_pred * (Y_pred - fe$exp_fe * fe$exp_pred))
  m_pred <- cbind(m_pred, 0)  # the score for log(theta) uses only labeled observations
  return(m_pred)
}

fenegbin_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  m_orig <- fenegbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_pred <- fenegbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

fenegbin_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  J_fe <- fepois_dsl_Jacobian(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)

  fe <- fepois_dsl_fe(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info, theta)
  Y_orig[labeled_ind == 0] <- 0
  s_orig <- negbin_dsl_parts(Y_orig, fe$exp_fe * fe$exp_orig, theta)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  J_theta <- - sum(s_orig$s_phi_phi * r_pi)/nrow(X_pred)

  J <- rbind(cbind(as.matrix(J_fe), 0), c(rep(0, ncol(J_fe)), J_theta))
  return(J)
}
