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
fepois_dsl_fe <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
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

  Y_dr <- Y_pred + (Y_orig - Y_pred) * r_pi
  w_dr <- exp_pred + (exp_orig - exp_pred) * r_pi
  v_dr <- Z_pred * exp_pred + (Z_orig * exp_orig - Z_pred * exp_pred) * r_pi

  # fe_info$index takes values 1, ..., G. So, the g-th row of rowsum() corresponds to group g
  w_sum <- as.numeric(rowsum(w_dr, fe_info$index))
  w_sum[w_sum <= 0] <- NA  # exp_fe_g is not defined
  exp_fe <- (as.numeric(rowsum(Y_dr, fe_info$index))/w_sum)[fe_info$index]
  center <- (rowsum(v_dr, fe_info$index)/w_sum)[fe_info$index, , drop = FALSE]

  out <- list("exp_fe" = exp_fe, "center" = center, "exp_orig" = exp_orig, "exp_pred" = exp_pred,
              "Z_orig" = Z_orig, "Z_pred" = Z_pred, "d_orig" = d_orig, "d_pred" = d_pred)
  return(out)
}

fepois_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)

  m_orig <- (fe$Z_orig - fe$center) * as.numeric(Y_orig - fe$exp_fe * fe$exp_orig)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  m_pred <- (fe$Z_pred - fe$center) * as.numeric(Y_pred - fe$exp_fe * fe$exp_pred)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

fepois_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_orig <- (fe$Z_orig - fe$center) * as.numeric(Y_orig - fe$exp_fe * fe$exp_orig)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

fepois_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_pred <- (fe$Z_pred - fe$center) * as.numeric(Y_pred - fe$exp_fe * fe$exp_pred)
  return(m_pred)
}

fepois_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info){
  # Jacobian of the concentrated moment. Centering Z at center_g accounts for the dependence of exp_fe_g on par.
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)

  # prediction
  diag_pred2 <- Diagonal(x = fe$exp_fe * (1 - labeled_ind/sample_prob_use))
  grad_pred <- (t(fe$Z_pred - fe$center) %*% diag_pred2 %*% fe$d_pred)/nrow(X_pred)

  # original
  diag_orig2_R <- Diagonal(x = fe$exp_fe * labeled_ind/sample_prob_use)
  grad_orig <- (t(fe$Z_orig - fe$center) %*% diag_orig2_R %*% fe$d_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig

  out <- grad_main
  return(out)
}

# ###############
# negbin, fenegbin
# ###############
# Negative binomial (NB2) regression: Var(Y | X) = mu + mu^2/theta, with mu = exp(X par_X) (negbin) or
# mu = exp_fe_g * kappa_t * exp(X par_X) (fenegbin; fixed effects on the linear scale, as in fepois).
# par = (par_X, kappa, log(theta)). Moments are the maximum likelihood scores:
#   (1) mean:       Z * (Y - mu) * theta/(theta + mu)  (Poisson score weighted by theta/(theta + mu))
#   (2) dispersion: d loglik / d log(theta)
# The score for the mean is valid when E(Y | X) = mu even if Y does not follow the negative binomial distribution.
# (fenegbin) The fixed effects solve the moment of each group, sum_{i in g} (Y - mu) * theta/(theta + mu) = 0. They have no closed form,
# so they are estimated jointly with par (see `dsl_negbin_newton`), and the moments for par account for their estimation.
# When exp_fe_g <= 0 (a group whose DSL-corrected sum of Y is negative; see fepois), mu <= 0 and we use Poisson weights (theta/(theta + 0) = 1).

# Scores for one observation given Y, mu, and theta (and their derivatives)
negbin_dsl_parts <- function(Y, mu, theta){
  mu_pos <- pmax(mu, 0)
  Y_pos  <- pmax(Y, 0)  # digamma(Y + theta) requires Y + theta > 0 (only matters when predictions are negative)

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

# Moments and Jacobian blocks. M1: moments for par. M2: moments for the fixed effects (fenegbin).
# A = - d M1/d par, B = - d M1/d exp_fe, C = - d M2/d par, D = - d M2/d exp_fe (diagonal). All are averaged over observations.
negbin_dsl_block <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info = NULL){
  n <- nrow(X_pred)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0

  theta <- exp(par[length(par)])
  par_X <- par[1:ncol(X_pred)]
  exp_pred_X <- as.numeric(exp(X_pred %*% par_X))
  exp_orig_X <- as.numeric(exp(X_orig %*% par_X))

  if(is.null(fe_info) == TRUE){
    exp_fe <- kappa <- 1
    Z_pred <- X_pred
    Z_orig <- X_orig
  }else{
    exp_fe <- fe_info$exp_fe[fe_info$index]
    if(is.null(fe_info$dummy) == TRUE){
      kappa  <- 1
      Z_pred <- X_pred
      Z_orig <- X_orig
    }else{
      kappa  <- as.numeric(1 + fe_info$dummy %*% (par[(ncol(X_pred) + 1):(length(par) - 1)] - 1))
      Z_pred <- cbind(X_pred, fe_info$dummy)
      Z_orig <- cbind(X_orig, fe_info$dummy)
    }
  }
  mu_pred <- exp_fe * kappa * exp_pred_X
  mu_orig <- exp_fe * kappa * exp_orig_X

  # derivatives of mu with respect to (par_X, kappa) and exp_fe
  if(is.null(fe_info$dummy) == TRUE){
    dmu_pred <- X_pred * mu_pred
    dmu_orig <- X_orig * mu_orig
  }else{
    dmu_pred <- cbind(X_pred * mu_pred, fe_info$dummy * exp_fe * exp_pred_X)
    dmu_orig <- cbind(X_orig * mu_orig, fe_info$dummy * exp_fe * exp_orig_X)
  }
  dmu_fe_pred <- kappa * exp_pred_X
  dmu_fe_orig <- kappa * exp_orig_X

  s_pred <- negbin_dsl_parts(Y_pred, mu_pred, theta)
  s_orig <- negbin_dsl_parts(Y_orig, mu_orig, theta)

  M1_pred <- cbind(Z_pred * s_pred$resid_w, s_pred$s_phi)
  M1_orig <- cbind(Z_orig * s_orig$resid_w, s_orig$s_phi)
  M1_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  M2_pred <- s_pred$resid_w
  M2_orig <- s_orig$resid_w
  M2_orig[labeled_ind == 0] <- 0

  # m_dr = (1 - r/pi) * m_pred + r/pi * m_orig
  w_pred <- 1 - r_pi
  w_orig <- r_pi
  A_1 <- t(Z_pred) %*% (dmu_pred * (s_pred$a * w_pred)) + t(Z_orig) %*% (dmu_orig * (s_orig$a * w_orig))
  A_2 <- - colSums(Z_pred * (s_pred$resid_w_phi * w_pred)) - colSums(Z_orig * (s_orig$resid_w_phi * w_orig))
  A_3 <- - colSums(dmu_pred * (s_pred$s_phi_mu * w_pred)) - colSums(dmu_orig * (s_orig$s_phi_mu * w_orig))
  A_4 <- - sum(s_pred$s_phi_phi * w_pred) - sum(s_orig$s_phi_phi * w_orig)
  A <- rbind(cbind(A_1, A_2), c(A_3, A_4))/n

  out <- list("M1_pred" = M1_pred, "M1_orig" = M1_orig, "M2_pred" = M2_pred, "M2_orig" = M2_orig, "A" = A, "r_pi" = r_pi,
              "mu_pred" = mu_pred, "mu_orig" = mu_orig, "Y_orig" = Y_orig)

  if(is.null(fe_info) == FALSE){
    B_obs <- cbind(Z_pred * (s_pred$a * dmu_fe_pred * w_pred) + Z_orig * (s_orig$a * dmu_fe_orig * w_orig),
                   - s_pred$s_phi_mu * dmu_fe_pred * w_pred - s_orig$s_phi_mu * dmu_fe_orig * w_orig)
    C_obs <- cbind(dmu_pred * (s_pred$a * w_pred) + dmu_orig * (s_orig$a * w_orig),
                   - s_pred$resid_w_phi * w_pred - s_orig$resid_w_phi * w_orig)
    out$B <- rowsum(B_obs, fe_info$index)/n
    out$C <- rowsum(C_obs, fe_info$index)/n
    out$D <- as.numeric(rowsum(s_pred$a * dmu_fe_pred * w_pred + s_orig$a * dmu_fe_orig * w_orig, fe_info$index))/n
  }
  return(out)
}

# Moments for par. (fenegbin) They account for the estimation of the fixed effects: m_1 - (B/D) * m_2 (see fepois).
negbin_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info = NULL){
  blk <- negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_orig <- blk$M1_orig
  if(is.null(fe_info) == FALSE){
    m_orig <- m_orig - (blk$B/blk$D)[fe_info$index, , drop = FALSE] * blk$M2_orig
  }
  return(m_orig)
}

negbin_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info = NULL){
  blk <- negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  m_pred <- blk$M1_pred
  if(is.null(fe_info) == FALSE){
    m_pred <- m_pred - (blk$B/blk$D)[fe_info$index, , drop = FALSE] * blk$M2_pred
  }
  return(m_pred)
}

negbin_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info = NULL){
  blk <- negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_info)
  J <- blk$A
  if(is.null(fe_info) == FALSE){
    J <- J - t(blk$B/blk$D) %*% blk$C # Schur complement (the fixed effects are concentrated out)
  }
  return(J)
}
