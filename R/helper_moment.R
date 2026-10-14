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

# #############################
# Count models: shared helpers
# #############################
# Count models (poisson, fepois, negbin, fenegbin) take `info`, a list with
#   offset: offset of the linear predictor (length n; 0 when the formula has no offset).
#   pred:   (discrete predicted covariates) list(X, prob, id). X stacks the design matrices X^(k) of all observations for every
#           combination k of the levels of the predicted discrete covariates, prob is the predicted probability of the combination,
#           and id is the observation of each row. The prediction part of each moment is then sum_k prob_k * m(Y_pred, X^(k)),
#           the expectation of the moment over the predicted distribution of the covariates.
#           NULL when no discrete covariate is predicted. The prediction part is then m(Y_pred, X_pred).
#   (fepois, fenegbin) index, dummy, pairwise, omega, and pilot (see fepois and fenegbin).
count_stack <- function(X_pred, Y_pred, info){
  if(is.null(info$pred) == TRUE){
    return(list("X" = X_pred, "Y" = Y_pred, "offset" = info$offset, "prob" = 1, "id" = NULL))
  }
  K <- length(info$pred$id)/nrow(X_pred)
  out <- list("X" = info$pred$X, "Y" = rep(Y_pred, K), "offset" = rep(info$offset, K), "prob" = info$pred$prob, "id" = info$pred$id)
  return(out)
}

# For each observation, sum_k prob_k * val_k over the rows of the stacked prediction design
count_collapse <- function(val, st){
  if(is.null(st$id) == TRUE){
    return(val)
  }
  out <- rowsum(val * st$prob, st$id)
  if(is.null(dim(val)) == TRUE){
    out <- as.numeric(out)
  }
  return(out)
}

# Values of observations (e.g., r/pi) for the rows of the stacked prediction design
count_expand <- function(val, st){
  if(is.null(st$id) == TRUE){
    return(val)
  }
  return(val[st$id])
}

# #######
# poisson
# #######
# Poisson regression with the log link. The moment is the Poisson score X * (Y - exp(X par + offset)).
# As the variance is estimated with the sandwich formula, this is Poisson Pseudo Maximum Likelihood (PPML):
# it only requires E(Y | X) = exp(X par + offset), and Y does not need to be an integer or follow a Poisson distribution.
poisson_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  m_orig <- poisson_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  m_pred <- poisson_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

poisson_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  mu_orig <- exp(X_orig %*% par + info$offset)
  m_orig <- X_orig * as.numeric(Y_orig - mu_orig)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

poisson_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  st <- count_stack(X_pred, Y_pred, info)
  mu_pred <- exp(st$X %*% par + st$offset)
  m_pred <- st$X * as.numeric(st$Y - mu_pred)
  return(count_collapse(m_pred, st))
}

poisson_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  r_pi <- as.numeric(labeled_ind/sample_prob_use)

  # prediction (d/d par of exp(X par) is exp(X par) * X). The weight is (1 - r/pi) as m_dr = (1 - r/pi) * m_pred + r/pi * m_orig
  st <- count_stack(X_pred, Y_pred, info)
  mu_pred <- as.numeric(exp(st$X %*% par + st$offset))
  grad_pred <- t(st$X) %*% (st$X * (mu_pred * st$prob * (1 - count_expand(r_pi, st))))

  # original
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  mu_orig <- as.numeric(exp(X_orig %*% par + info$offset))
  grad_orig <- t(X_orig) %*% (X_orig * (mu_orig * r_pi))

  out <- (grad_pred + grad_orig)/nrow(X_orig)
  return(out)
}

# ######
# fepois
# ######
# Poisson regression with fixed effects (PPML with fixed effects). The conditional mean is exp_fe_g * kappa_t * exp(X par_X + offset),
# where exp_fe_g = exp(alpha_g) is the fixed effect of group g of the concentrated index, and
# kappa_t = exp(gamma_t) is the fixed effect of the other index (twoways only; kappa_t = 1 for the reference level).
# Fixed effects are on the linear scale (exp_fe_g and kappa_t rather than alpha_g and gamma_t). par = (par_X, kappa).
#
# (1) Only the outcome is predicted (info$pairwise = FALSE). exp_fe_g is concentrated out. For a given par, exp_fe_g solves the
#   DSL moment of the fixed effect, sum_{i in g} [Y_dr_i - exp_fe_g * w_dr_i] = 0, where Y_dr and w_dr are the DSL versions of Y and
#   kappa_t * exp(X par_X + offset). Thus, exp_fe_g = sum_g Y_dr / sum_g w_dr. It exists even when the DSL-corrected sum of Y within
#   a group is negative (e.g., a small group with a labeled observation with a large prediction error).
#   Substituting exp_fe_g gives the moments for par with the design Z = [X, dummy variables of the other index]
#   centered at center_g = sum_g DSL(Z * kappa_t * exp(X par_X)) / sum_g w_dr. This accounts for the estimation of exp_fe_g (as demeaning in felm).
#   As X is observed, the moments are linear in Y_dr, and they are exactly unbiased for any group size.
#
# (2) A covariate is predicted (info$pairwise = TRUE). exp_fe_g and center_g would then be ratios of DSL sums over the few labeled
#   observations of each group, which biases the coefficients (incidental parameter problem; the bias does not vanish as the number
#   of groups grows). Instead, we use the pairwise (conditional likelihood) form of the moment, which does not involve exp_fe_g:
#     M_g = sum_{t != s in g} (Z_t - Z_s) Y_t f_s,   f_s = kappa_s * exp(X_s par_X + offset_s).
#   E[(Z_t - Z_s) Y_t f_s] = exp_fe_g (Z_t - Z_s) f_t f_s, so M_g has mean zero (the terms of (t, s) and (s, t) cancel).
#   Writing M_g = sum_{t != s} [A_t B_s - Y_t C_s] with A = Z Y, B = f, C = Z f, every term is a product of quantities of two
#   different observations. We replace each quantity by its DSL version (e.g., A_dr = A_pred + r/pi * (A_orig - A_pred)). As the labeling
#   of different observations is independent, the product of the DSL versions is unbiased for the product, and M_g remains exactly
#   unbiased for any group size. Groups are weighted by omega_g = 1/(number of observations in g). The weights must not depend on
#   predicted variables. When nothing is predicted, M_g/sum_g f equals the moment in (1), and omega_g = 1/T_g gives nearly the same weights.
#   Moments of different observations within a group are dependent, so standard errors are clustered at least at the group level.
#
# (fenegbin) When theta is given, every quantity is multiplied by the negative binomial weight theta/(theta + mu_pilot) (see fenegbin).
fenegbin_weight <- function(X, offset, info, theta){
  mu_pilot <- info$pilot$scale * as.numeric(exp(X %*% info$pilot$par_X + offset))
  return(theta/(theta + mu_pilot))
}

fepois_dsl_fe <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0

  par_X <- par[1:ncol(X_pred)]
  exp_pred_X <- as.numeric(exp(X_pred %*% par_X + info$offset))
  exp_orig_X <- as.numeric(exp(X_orig %*% par_X + info$offset))
  exp_orig_X[labeled_ind == 0] <- 0

  if(is.null(info$dummy) == TRUE){
    kappa  <- 1
    Z_pred <- X_pred
    Z_orig <- X_orig
  }else{
    kappa  <- as.numeric(1 + info$dummy %*% (par[-(1:ncol(X_pred))] - 1))
    Z_pred <- cbind(X_pred, info$dummy)
    Z_orig <- cbind(X_orig, info$dummy)
  }
  exp_pred <- kappa * exp_pred_X
  exp_orig <- kappa * exp_orig_X

  # derivative of kappa_t * exp(X par_X) with respect to par
  if(is.null(info$dummy) == TRUE){
    d_pred <- X_pred * exp_pred
    d_orig <- X_orig * exp_orig
  }else{
    d_pred <- cbind(X_pred * exp_pred, info$dummy * exp_pred_X)
    d_orig <- cbind(X_orig * exp_orig, info$dummy * exp_orig_X)
  }

  # (fenegbin) weights theta/(theta + mu_pilot). mu_pilot is fixed during estimation (see fenegbin)
  if(is.null(theta) == TRUE){
    weight_pred <- weight_orig <- 1
  }else{
    weight_pred <- fenegbin_weight(X_pred, info$offset, info, theta)
    weight_orig <- fenegbin_weight(X_orig, info$offset, info, theta)
  }

  Y_dr <- Y_pred * weight_pred + (Y_orig * weight_orig - Y_pred * weight_pred) * r_pi
  w_dr <- exp_pred * weight_pred + (exp_orig * weight_orig - exp_pred * weight_pred) * r_pi
  v_dr <- Z_pred * (exp_pred * weight_pred) + (Z_orig * (exp_orig * weight_orig) - Z_pred * (exp_pred * weight_pred)) * r_pi

  # info$index takes values 1, ..., G. So, the g-th row of rowsum() corresponds to group g
  w_sum <- as.numeric(rowsum(w_dr, info$index))
  w_sum[w_sum <= 0] <- NA  # exp_fe_g is not defined
  exp_fe <- (as.numeric(rowsum(Y_dr, info$index))/w_sum)[info$index]
  center <- (rowsum(v_dr, info$index)/w_sum)[info$index, , drop = FALSE]

  out <- list("exp_fe" = exp_fe, "center" = center, "exp_orig" = exp_orig, "exp_pred" = exp_pred,
              "Z_orig" = Z_orig, "Z_pred" = Z_pred, "d_orig" = d_orig, "d_pred" = d_pred,
              "weight_orig" = weight_orig, "weight_pred" = weight_pred)
  return(out)
}

# Pairwise moment (2). Returns m_pred and m_orig for each observation such that, within group g,
#   sum_t [m_pred_t + r_t/pi_t * (m_orig_t - m_pred_t)] = omega_g * M_g,
# and (when jacobian = TRUE) J = - d mean(m_dr)/d par.
# For two different observations t and s, the product of DSL versions U_dr_t * V_dr_s decomposes exactly as
#   U_pred_t * V_pred_s + r_t/pi_t * (U_orig_t - U_pred_t) * V_dr_s + r_s/pi_s * (V_orig_s - V_pred_s) * U_pred_t.
# m_pred collects the first terms and m_orig - m_pred the others, so the meat (with clusters containing groups) is exact.
fepois_dsl_pair <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL, jacobian = FALSE, meat = FALSE){
  n <- nrow(X_orig)
  p <- ncol(X_orig)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0
  st <- count_stack(X_pred, Y_pred, info)
  dummy_st <- info$dummy
  if(is.null(info$dummy) == FALSE & is.null(st$id) == FALSE){
    dummy_st <- info$dummy[st$id, , drop = FALSE]
  }

  # f = kappa_t * exp(X par_X + offset), Z = [X, dummy], D = d f/d par, and (fenegbin) weights v
  fe_parts <- function(X, offset, dummy){
    e <- as.numeric(exp(X %*% par[1:p] + offset))
    if(is.null(dummy) == TRUE){
      out <- list("f" = e, "Z" = X, "D" = X * e)
    }else{
      kappa <- as.numeric(1 + dummy %*% (par[-(1:p)] - 1))
      out <- list("f" = kappa * e, "Z" = cbind(X, dummy), "D" = cbind(X * (kappa * e), dummy * e))
    }
    if(is.null(theta) == TRUE){
      out$v <- 1
    }else{
      out$v <- fenegbin_weight(X, offset, info, theta)
    }
    return(out)
  }
  po <- fe_parts(X_orig, info$offset, info$dummy)
  ps <- fe_parts(st$X, st$offset, dummy_st)

  # Quantities A = v Z Y, B = v f, C = v Z f, and Y_w = v Y: prediction (expectation over the predicted distribution) and original
  a_pred <- count_collapse(ps$Z * (ps$v * st$Y), st)
  b_pred <- count_collapse(ps$v * ps$f, st)
  c_pred <- count_collapse(ps$Z * (ps$v * ps$f), st)
  y_pred <- count_collapse(ps$v * st$Y, st)
  a_orig <- po$Z * (po$v * Y_orig)
  b_orig <- po$v * po$f
  c_orig <- po$Z * (po$v * po$f)
  y_orig <- po$v * Y_orig

  # DSL versions
  A_dr <- a_pred + (a_orig - a_pred) * r_pi
  B_dr <- b_pred + (b_orig - b_pred) * r_pi
  C_dr <- c_pred + (c_orig - c_pred) * r_pi
  Y_dr <- y_pred + (y_orig - y_pred) * r_pi

  # Sums within groups (info$index takes values 1, ..., G)
  g <- info$index
  omega <- info$omega[g]
  sum_g <- function(x){
    out <- rowsum(x, g)
    if(is.null(dim(x)) == TRUE){
      out <- as.numeric(out)
    }
    return(out)
  }
  SA_dr <- sum_g(A_dr); SB_dr <- sum_g(B_dr); SC_dr <- sum_g(C_dr); SY_dr <- sum_g(Y_dr)
  SA_pred <- sum_g(a_pred); SB_pred <- sum_g(b_pred); SC_pred <- sum_g(c_pred); SY_pred <- sum_g(y_pred)

  # M_g = sum_{t != s} [A_t B_s - Y_t C_s]
  m_pred <- (a_pred * (SB_pred[g] - b_pred) - (SC_pred[g, , drop = FALSE] - c_pred) * y_pred) * omega
  m_corr <- (a_orig - a_pred) * (SB_dr[g] - B_dr) + (SA_pred[g, , drop = FALSE] - a_pred) * (b_orig - b_pred) -
    (SC_dr[g, , drop = FALSE] - C_dr) * (y_orig - y_pred) - (c_orig - c_pred) * (SY_pred[g] - y_pred)
  m_orig <- m_pred + m_corr * omega
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  colnames(m_pred) <- colnames(m_orig) <- NULL

  out <- list("m_pred" = m_pred, "m_orig" = m_orig)

  if(jacobian == TRUE){
    # d M_g/d par = SA_g SDb_g' - sum_t A_t Db_t' + sum_t (Y_t - SY_g) Dc_t, where Db and Dc are the DSL versions of
    # d B/d par = v D and d C/d par = v Z D' (the weights v are fixed)
    Db_pred <- count_collapse(ps$D * ps$v, st)
    Db_dr <- Db_pred + (po$D * po$v - Db_pred) * r_pi
    SDb_dr <- rowsum(Db_dr, g)
    w <- omega * (Y_dr - SY_dr[g])
    grad_1 <- t(SA_dr * info$omega) %*% SDb_dr
    grad_2 <- t(A_dr * omega) %*% Db_dr
    grad_3 <- t(ps$Z * (count_expand(w * (1 - r_pi), st) * st$prob * ps$v)) %*% ps$D + t(po$Z * (w * r_pi * po$v)) %*% po$D
    J <- (grad_2 - grad_1 - grad_3)/n
    dimnames(J) <- NULL
    out$J <- J
  }

  if(meat == TRUE){
    # Meat when observations are independent (no `cluster`), from the variance of the U-statistic (see fepois_dsl_pair_meat).
    # main_1 is the part from the DSL corrections (difference from the prediction-only moment), as in other models.
    meat_dr   <- fepois_dsl_pair_meat(A_dr, B_dr, C_dr, Y_dr, g, info$omega)
    psi_corr  <- meat_dr$psi - fepois_dsl_pair_meat(a_pred, b_pred, c_pred, y_pred, g, info$omega, psi_only = TRUE)$psi
    out$main_1  <- crossprod(psi_corr * omega)/n
    out$main_23 <- meat_dr$meat/n - out$main_1
  }
  return(out)
}

# Variance of sum_g omega_g * M_g when observations are independent (also within groups).
# M_g = sum_{t < s} h_ts is a U-statistic with h_ts = P_ts + P_st and P_ts = A_t B_s - Y_t C_s. With psi_t = sum_{s != t} h_ts,
#   E[sum_t psi_t psi_t' - sum_{t < s} h_ts h_ts'] = Var(M_g),
# because pairs of pairs that share no observation are independent. Every term is computed from sums within groups (no loop over pairs).
fepois_dsl_pair_meat <- function(A, B, C, Y, g, omega, psi_only = FALSE){
  w  <- omega[g]
  SA <- rowsum(A, g); SB <- as.numeric(rowsum(B, g)); SC <- rowsum(C, g); SY <- as.numeric(rowsum(Y, g))
  psi <- A * (SB[g] - B) - (SC[g, , drop = FALSE] - C) * Y + (SA[g, , drop = FALSE] - A) * B - C * (SY[g] - Y)
  if(psi_only == TRUE){
    return(list("psi" = psi))
  }

  w2 <- w^2
  # sum_g omega_g^2 sum_{t != s} U_t V_t' c_s
  type_1 <- function(U, V, c){
    Sc <- as.numeric(rowsum(c, g))
    return(crossprod(U * (w2 * (Sc[g] - c)), V))
  }
  # sum_g omega_g^2 sum_{t != s} U_t V_s'
  type_2 <- function(U, V){
    return(crossprod(rowsum(U, g) * omega^2, rowsum(V, g)) - crossprod(U * w2, V))
  }
  T_2 <- type_2(A * Y, C * B)
  half_hh <- type_1(A, A, B^2) - T_2 - t(T_2) + type_1(C, C, Y^2) +
    type_2(A * B, A * B) - type_1(A, C, B * Y) - type_1(C, A, B * Y) + type_2(C * Y, C * Y)

  meat <- crossprod(psi * w) - half_hh
  dimnames(meat) <- NULL
  return(list("meat" = meat, "psi" = psi))
}

fepois_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  m_orig <- fepois_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_pred <- fepois_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

fepois_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  if(isTRUE(info$pairwise)){
    return(fepois_dsl_pair(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)$m_orig)
  }
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_orig <- (fe$Z_orig - fe$center) * as.numeric(fe$weight_orig * (Y_orig - fe$exp_fe * fe$exp_orig))
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

fepois_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  if(isTRUE(info$pairwise)){
    return(fepois_dsl_pair(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)$m_pred)
  }
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_pred <- (fe$Z_pred - fe$center) * as.numeric(fe$weight_pred * (Y_pred - fe$exp_fe * fe$exp_pred))
  return(m_pred)
}

fepois_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  if(isTRUE(info$pairwise)){
    return(fepois_dsl_pair(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta, jacobian = TRUE)$J)
  }
  # Jacobian of the concentrated moment. Centering Z at center_g accounts for the dependence of exp_fe_g on par.
  # (fenegbin) The weights are fixed.
  fe <- fepois_dsl_fe(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)

  # prediction
  diag_pred2 <- Diagonal(x = fe$exp_fe * fe$weight_pred * (1 - labeled_ind/sample_prob_use))
  grad_pred <- (t(fe$Z_pred - fe$center) %*% diag_pred2 %*% fe$d_pred)/nrow(X_pred)

  # original
  diag_orig2_R <- Diagonal(x = fe$exp_fe * fe$weight_orig * labeled_ind/sample_prob_use)
  grad_orig <- (t(fe$Z_orig - fe$center) %*% diag_orig2_R %*% fe$d_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig

  out <- as.matrix(grad_main)
  return(out)
}

# ######
# negbin
# ######
# Negative binomial (NB2) regression: Var(Y | X) = mu + mu^2/theta with mu = exp(X par_X + offset). par = (par_X, log(theta)).
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

# Moments (orig and pred) and Jacobian (J = - d m_n/d par).
# When theta is given, par = par_X and the block contains only the moments for the mean (used when solving for par_X given theta).
negbin_dsl_block <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  n <- nrow(X_pred)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0

  theta_given <- is.null(theta) == FALSE
  if(theta_given == FALSE){
    theta <- exp(par[length(par)])
  }
  par_X <- par[1:ncol(X_pred)]
  st <- count_stack(X_pred, Y_pred, info)
  mu_pred <- as.numeric(exp(st$X %*% par_X + st$offset))
  mu_orig <- as.numeric(exp(X_orig %*% par_X + info$offset))

  s_pred <- negbin_dsl_parts(st$Y, mu_pred, theta)
  s_orig <- negbin_dsl_parts(Y_orig, mu_orig, theta)

  m_pred <- count_collapse(st$X * s_pred$resid_w, st)
  m_orig <- X_orig * s_orig$resid_w

  # m_dr = (1 - r/pi) * m_pred + r/pi * m_orig
  w_pred <- st$prob * (1 - count_expand(r_pi, st))
  w_orig <- r_pi
  J_1 <- t(st$X) %*% (st$X * (s_pred$a * mu_pred * w_pred)) + t(X_orig) %*% (X_orig * (s_orig$a * mu_orig * w_orig))

  if(theta_given == TRUE){
    m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
    J <- J_1/n
  }else{
    m_pred <- cbind(m_pred, 0)  # the score for log(theta) uses only labeled observations
    m_orig <- cbind(m_orig, s_orig$s_phi)
    m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

    J_2 <- - colSums(st$X * (s_pred$resid_w_phi * w_pred)) - colSums(X_orig * (s_orig$resid_w_phi * w_orig))
    J_3 <- - colSums(X_orig * (s_orig$s_phi_mu * mu_orig * w_orig))
    J_4 <- - sum(s_orig$s_phi_phi * w_orig)
    J <- rbind(cbind(J_1, J_2), c(J_3, J_4))/n
  }
  dimnames(J) <- NULL

  out <- list("m_pred" = m_pred, "m_orig" = m_orig, "J" = J, "mu_orig" = mu_orig)
  return(out)
}

negbin_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  blk <- negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_dr <- blk$m_pred + (blk$m_orig - blk$m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

negbin_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)$m_orig)
}

negbin_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)$m_pred)
}

negbin_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta = NULL){
  return(negbin_dsl_block(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)$J)
}

# Mean moment of the score for log(theta) (labeled observations weighted by w = r/pi), given mu of labeled observations
negbin_dsl_theta_moment <- function(log_theta, Y_orig, mu_orig, w){
  s <- negbin_dsl_parts(Y_orig, mu_orig, exp(log_theta))
  return(sum(w * s$s_phi)/length(w))
}

# ########
# fenegbin
# ########
# Negative binomial (NB2) regression with fixed effects. par = (par_X, kappa, log(theta)) (see fepois for kappa).
# The moments for (par_X, kappa) are those of fepois with every observation weighted by the negative binomial weight
# theta/(theta + mu_pilot), and the moment for log(theta) is the maximum likelihood score (see negbin).
# Maximum likelihood would use the weights theta/(theta + mu) with each group's estimated fixed effect. The weights of a group then depend
# on the outcomes (or predictions) of the same group, which biases the coefficients in short panels (incidental parameter problem).
# Thus, mu_pilot = scale * exp(X par_X_pilot + offset) uses the fepois estimates par_X_pilot and one overall scale
# (sum Y_dr / sum DSL(exp(X par_X_pilot + offset))), without fixed effects, and it is fixed during estimation.
# As the weights are functions of the covariates only, the moments for (par_X, kappa) remain unbiased, as in fepois (including the pairwise
# moment when covariates are predicted). The weights only affect efficiency.
# The score for log(theta) uses only labeled observations weighted by 1/pi, with mu = exp_fe_g * kappa_t * exp(X par_X + offset) and
# exp_fe_g = sum_g Y_dr / sum_g DSL(kappa_t * exp(X par_X + offset)). theta is an auxiliary parameter for the weights;
# with few labeled observations per group, it is biased downward (the estimated fixed effects add noise to the residuals).
# Cross-derivatives between the mean parameters and theta multiply (Y - mu) and have mean zero, so the Jacobian omits them.

# mu of labeled observations for the score of log(theta), and w = r/pi (0 for groups whose fixed effect is not defined)
fenegbin_dsl_mu_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  p <- ncol(X_orig)
  r_pi <- as.numeric(labeled_ind/sample_prob_use)
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  Y_orig[labeled_ind == 0] <- 0
  st <- count_stack(X_pred, Y_pred, info)
  kappa_fun <- function(dummy){
    if(is.null(dummy) == TRUE){
      return(1)
    }
    return(as.numeric(1 + dummy %*% (par[-(1:p)] - 1)))
  }
  dummy_st <- info$dummy
  if(is.null(info$dummy) == FALSE & is.null(st$id) == FALSE){
    dummy_st <- info$dummy[st$id, , drop = FALSE]
  }
  f_pred <- count_collapse(kappa_fun(dummy_st) * as.numeric(exp(st$X %*% par[1:p] + st$offset)), st)
  f_orig <- kappa_fun(info$dummy) * as.numeric(exp(X_orig %*% par[1:p] + info$offset))
  f_dr <- f_pred + (f_orig - f_pred) * r_pi
  Y_dr <- Y_pred + (Y_orig - Y_pred) * r_pi

  sum_f <- as.numeric(rowsum(f_dr, info$index))
  exp_fe <- as.numeric(rowsum(Y_dr, info$index))/sum_f
  valid <- (sum_f > 0) & (exp_fe >= 0)
  exp_fe[valid == FALSE] <- 0

  out <- list("mu_orig" = exp_fe[info$index] * f_orig, "w" = r_pi * valid[info$index], "Y_orig" = Y_orig)
  return(out)
}

fenegbin_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  m_orig <- fepois_dsl_moment_orig(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  th <- fenegbin_dsl_mu_orig(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  s_orig <- negbin_dsl_parts(th$Y_orig, th$mu_orig, theta)
  m_orig <- cbind(m_orig, s_orig$s_phi * as.numeric(th$w > 0))
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

fenegbin_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  m_pred <- fepois_dsl_moment_pred(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)
  m_pred <- cbind(m_pred, 0)  # the score for log(theta) uses only labeled observations
  return(m_pred)
}

fenegbin_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  m_orig <- fenegbin_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  m_pred <- fenegbin_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

fenegbin_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info){
  theta <- exp(par[length(par)])
  par_fe <- par[-length(par)]
  J_fe <- fepois_dsl_Jacobian(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info, theta)

  th <- fenegbin_dsl_mu_orig(par_fe, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, info)
  s_orig <- negbin_dsl_parts(th$Y_orig, th$mu_orig, theta)
  J_theta <- - sum(s_orig$s_phi_phi * th$w)/nrow(X_pred)

  J <- rbind(cbind(as.matrix(J_fe), 0), c(rep(0, ncol(J_fe)), J_theta))
  dimnames(J) <- NULL
  return(J)
}
