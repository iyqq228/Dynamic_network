# --- Line-graph S2: retain the original make_S2 (sparse S2 for hyperparameter updates) ---
{
  ## ---- thin wrappers to C++ ----
  vg_beta_cpp <- function(b, delta, scan, ties_method){
    cpp_value_grad_beta(
      beta  = b,
      delta = delta,
      X     = scan$X_rows,
      id    = as.integer(scan$id_row),
      start_sorted = scan$start_sorted,
      stop_sorted  = scan$stop_sorted,
      rows_start   = as.integer(scan$rows_start),
      rows_stop    = as.integer(scan$rows_stop),
      times        = scan$times,
      events_by_time = scan$events_rowidx_by_time,
      ties_method  = ties_method
    )
  }
  vg_delta_cpp <- function(delta, beta, scan, ties_method){
    cpp_value_grad_delta(
      delta = delta, beta = beta, X = scan$X_rows,
      id    = as.integer(scan$id_row),
      start_sorted = scan$start_sorted, stop_sorted = scan$stop_sorted,
      rows_start = as.integer(scan$rows_start), rows_stop = as.integer(scan$rows_stop),
      times = scan$times, events_by_time = scan$events_rowidx_by_time,
      ties_method = ties_method
    )
  }
  cumLambda_cpp <- function(scan, beta, delta, method){
    cpp_compute_cumLambda(
      X = scan$X_rows, beta = beta, delta = delta,
      id = as.integer(scan$id_row),
      start_sorted = scan$start_sorted, stop_sorted = scan$stop_sorted,
      rows_start = as.integer(scan$rows_start), rows_stop = as.integer(scan$rows_stop),
      times = scan$times, events_by_time = scan$events_rowidx_by_time,
      method = method
    )
  }
  K_diag_full_cpp <- function(scan, beta, delta, cumLambda_by_time){
    cpp_K_diag_full(
      X = scan$X_rows, beta = beta, delta = delta,
      id = as.integer(scan$id_row), times = scan$times,
      rows_start = as.integer(scan$rows_start), rows_stop = as.integer(scan$rows_stop),
      start_sorted = scan$start_sorted, stop_sorted = scan$stop_sorted,
      cumLambda_by_time = cumLambda_by_time
    )
  }
  K_matvec_cpp <- function(x, beta, delta, scan, ties_method){
    cpp_K_matvec(
      x = x, beta = beta, delta = delta, X = scan$X_rows,
      id = as.integer(scan$id_row),
      start_sorted = scan$start_sorted, stop_sorted = scan$stop_sorted,
      rows_start = as.integer(scan$rows_start), rows_stop = as.integer(scan$rows_stop),
      times = scan$times, events_by_time = scan$events_rowidx_by_time,
      ties_method = ties_method
    )
  }
  
  ## ======== blocks ========
  
  beta_block_update_tv_fast2_cpp <- function(beta0, delta, scan, maxit=50,
                                             ties_method=c("breslow","efron")){
    ties_method <- match.arg(ties_method)
    fn <- function(b) -vg_beta_cpp(b, delta, scan, ties_method)$value
    gr <- function(b) -vg_beta_cpp(b, delta, scan, ties_method)$grad
    opt <- optim(beta0, fn, gr, method="L-BFGS-B",
                 control=list(maxit=maxit, factr=1e7, pgtol=1e-8))
    list(beta = opt$par, value = -opt$value, convergence = opt$convergence)
  }
  
  approx_diag_K <- function(scan, beta, delta, ties_method=c("breslow","efron")){
    # Continue using the R approximation as the preconditioner; replace it with the full version if needed
    ties_method <- match.arg(ties_method)
    X  <- scan$X_rows; id <- scan$id_row
    a_row <- exp(as.vector(X %*% beta)); d_id <- exp(delta)
    ss <- scan$start_sorted; es <- scan$stop_sorted
    rs <- scan$rows_start;    re <- scan$rows_stop
    times <- scan$times; nR <- nrow(X); m <- length(delta)
    mt_by_time <- vapply(scan$events_rowidx_by_time, length, 0L)
    grp <- scan$events_rowidx_by_time
    Wid <- numeric(m); diagK <- numeric(m)
    si <- 1L; ei <- 1L; S0 <- 0.0
    for (k in seq_along(times)) {
      t <- times[k]
      while (si <= nR && ss[si] <= t) { r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0+wr; Wid[eid] <- Wid[eid]+wr; si <- si+1L }
      while (ei <= nR && es[ei] <  t) { r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0-wr; Wid[eid] <- Wid[eid]-wr; ei <- ei+1L }
      if (!(S0 > 0) || !is.finite(S0)) next
      mk <- mt_by_time[k]; if (mk == 0L) next
      if (ties_method == "breslow") {
        w <- Wid / S0; diagK <- diagK + mk * (w * (1 - w))
      } else {
        ev_rows <- grp[[k]]
        wevent <- numeric(m)
        if (length(ev_rows)) {
          wr_ev <- a_row[ev_rows] * d_id[id[ev_rows]]
          acc   <- rowsum(wr_ev, group = id[ev_rows], reorder = FALSE)
          wevent[as.integer(rownames(acc))] <- as.numeric(acc)
        }
        Ek <- sum(wevent)
        if (Ek <= 0 || !is.finite(Ek)) {
          w <- Wid / S0; diagK <- diagK + mk * (w * (1 - w))
        } else {
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk; denom <- S0 - alpha * Ek
            if (denom <= 0 || !is.finite(denom)) next
            w_ell <- (Wid - alpha * wevent) / denom
            diagK <- diagK + w_ell * (1 - w_ell)
          }
        }
      }
    }
    pmax(diagK, 0)
  }
  
  delta_block_update_fista_tv_fast2_cpp <- function(
    delta0, beta, scan, ops, sigma2, rho,
    max_iter=10, tol_grad_rms=1e-4,
    backtrack=TRUE, step_init=NULL, bt_factor=0.5, ls_armijo=1e-4,
    grad_clip_rms=1e3,
    chol_threshold=25000L, pcg_tol=1e-6, pcg_maxit=2000,
    ties_method=c("breslow","efron")
  ){
    ties_method <- match.arg(ties_method)
    Minv <- make_M_solver(ops, rho, chol_threshold, pcg_tol, pcg_maxit)
    
    K_diag <- approx_diag_K_cpp(scan, beta, delta0, ties_method)
    epsK <- 1e-12
    Hdiag <- (1 / sigma2) + pmax(as.numeric(K_diag), epsK)
    invHdiag <- 1 / Hdiag
    step <- if (is.null(step_init)) 1/mean(Hdiag) else min(step_init, 1/mean(Hdiag))
    step <- step * 20
    eval_f_g <- function(delta){
      core <- vg_delta_cpp(delta, beta, scan, ties_method)
      Minv_delta <- Minv(delta)
      f <- core$value - 0.5 * sum(delta * Minv_delta) / sigma2
      g <- core$grad_delta - Minv_delta / sigma2
      g_rms <- sqrt(mean(g^2))
      if (g_rms > grad_clip_rms) g <- g * (grad_clip_rms / g_rms)
      list(f=f, g=g)
    }
    
    y <- delta0; delta <- delta0; t_k <- 1
    eg0 <- eval_f_g(delta0); f_prev <- eg0$f
   
    for (it in seq_len(max_iter)){
      print(it)
      eg_y <- eval_f_g(y); g <- eg_y$g
      g_rms <- sqrt(mean(g^2))
      if (!is.finite(g_rms) || g_rms < tol_grad_rms) break
      
      d <- g * invHdiag
      s <- step; f_y <- eg_y$f; gain_lin <- sum(g * d)
      delta_try <- y + s * d
      eg_try <- eval_f_g(delta_try); f_try <- eg_try$f
      
      if (backtrack){
        thresh <- f_y + ls_armijo * s * gain_lin
        while (!is.finite(f_try) || f_try < thresh){
          s <- s * bt_factor
          if (s < 1e-16) return(list(delta=delta, value=f_prev, iters=it-1L, step=step))
          delta_try <- y + s*d
          eg_try <- eval_f_g(delta_try); f_try <- eg_try$f
          thresh <- f_y + ls_armijo * s * gain_lin
        }
      }
      
      t_k1 <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
      if (!is.finite(f_try) || f_try < f_prev){
        t_k1 <- 1.0
        y_anchor <- delta
        s2 <- min(s * bt_factor, step)
        delta_try <- y_anchor + s2 * d
        eg_try <- eval_f_g(delta_try); f_try <- eg_try$f
        s <- s2
      }
      
      step <- s
      y_new <- delta_try + ((t_k - 1)/t_k1) * (delta_try - delta)
      delta <- delta_try; y <- y_new; t_k <- t_k1; f_prev <- f_try
    }
    list(delta=delta, value=f_prev, iters=it-1L, step=step)
  }
  
  ## ======== Hyperparameter update (the original ratio_fix version is unchanged) ========
  update_sigma_rho_ratio_fix <- function(
    delta, n, S2, rho_old, sigma2_old,
    K_matvec=NULL, K_diag=NULL,
    damping=0.3,
    chol_threshold=5000L, n_probe=50, seed=1L,
    pcg_tol=1e-6, pcg_maxit=1000,
    sigma2_min=1e-6, sigma2_max=1e+3
  ){
    stopifnot(inherits(S2, "sparseMatrix"))
    m <- length(delta)
    if (is.null(K_diag)) K_diag <- rep(0, m)
    K_diag <- as.numeric(K_diag)
    
    degL <- as.numeric(Matrix::rowSums(S2))
    use_chol <- (m <= chol_threshold)
    if (use_chol){
      M <- Diagonal(m) + rho_old * S2
      cf <- tryCatch(Cholesky(M, LDL=FALSE, Imult=0, perm=TRUE),
                     error=function(e) Cholesky(M + Diagonal(m,1e-10), LDL=FALSE, Imult=0, perm=TRUE))
      A_apply <- function(x) as.numeric(solve(cf, x))
      A_apply_many <- function(B) as.matrix(solve(cf, B))
    } else {
      M_mv <- function(x) x + as.numeric(rho_old * (S2 %*% x))
      Mdiag <- pmax(1 + rho_old * degL, 1e-12)
      A_apply <- function(b) pcg_solve(M_mv, b, Mdiag=Mdiag, tol=pcg_tol, maxit=pcg_maxit)
      A_apply_many <- function(B) apply(B, 2L, A_apply)
    }
    
    q0 <- sum(delta * A_apply(delta))
    Oinv_vec <- function(x) A_apply(x) / sigma2_old
    H_mv <- function(x) Oinv_vec(x) + K_diag * x
    Hdiag  <- (1 / sigma2_old) + K_diag
    Hsolve <- function(b) pcg_solve(H_mv, b, Mdiag=pmax(Hdiag,1e-12),
                                    tol=min(1e-8, pcg_tol), maxit=max(3000, pcg_maxit))
    
    set.seed(seed); S <- max(10L, as.integer(n_probe))
    V <- matrix(sample(c(-1,1), m*S, replace=TRUE), m, S)
    
    Av <- A_apply_many(V)
    HinvAv <- apply(Av, 2L, Hsolve)
    t0 <- mean(colSums(V * HinvAv))
    
    sigma2_corr <- (q0+t0)/m
    sigma2_corr <- max(sigma2_min, min(sigma2_corr, sigma2_max))
    sigma2_new  <- (1 - damping) * sigma2_old + damping * sigma2_corr
    
    a <- A_apply(delta); Sa <- as.numeric(S2 %*% a); q1 <- sum(a * Sa)
    AS2V <- A_apply_many(as.matrix(S2 %*% V)); t1 <- mean(colSums(V * AS2V))
    S2Av <- as.matrix(S2 %*% Av); AS2Av <- A_apply_many(S2Av)
    Hinv_AS2Av <- apply(AS2Av, 2L, Hsolve); t2 <- mean(colSums(V * Hinv_AS2Av))
    S2_AS2V <- as.matrix(S2 %*% AS2V); AS2AS2V <- A_apply_many(S2_AS2V); t_SS <- mean(colSums(V * AS2AS2V))
    AaSa <- A_apply(Sa); d_tmp <- as.numeric(S2 %*% AaSa); Ad_tmp <- A_apply(d_tmp); u_SS_delta <- sum(delta * Ad_tmp)
    S2_AS2Av2 <- as.matrix(S2 %*% AS2Av); AS2AS2Av2 <- A_apply_many(S2_AS2Av2)
    Hinv_AS2AS2Av <- apply(AS2AS2Av2, 2L, Hsolve); t_SS_H <- mean(colSums(V * Hinv_AS2AS2Av))
    
    g_rho <- 0.5 * ( (q1 + t2) / sigma2_old - t1 )
    h_rho <- 0.5 * ( t_SS - (2 / sigma2_old) * u_SS_delta - (2 / sigma2_old) * t_SS_H )
    
    h_safe <- if (!is.finite(h_rho) || h_rho <= 0) max(abs(h_rho), 1e-8) else h_rho
    step_cap <- 0.10 * (0.5 - 0.0)
    delta_rho <- g_rho / h_safe
    delta_rho <- max(-step_cap, min(step_cap, delta_rho))
    
    rho_hat <- rho_old + delta_rho
    rho_lo <- 0.0; rho_hi <- 0.5; eps <- 1e-6
    rho_hat <- max(rho_lo + eps, min(rho_hat, rho_hi - eps))
    rho_new <- (1 - damping) * rho_old + damping * rho_hat
    
    list(sigma2=sigma2_new, rho=rho_new,
         stats=list(q0=q0, t0=t0, q1=q1, t1=t1, t2=t2, g_rho=g_rho, h_rho=h_rho))
  }
  
  ## ======== Top-level fit (C++-accelerated) ========
  
  fit_edge_frailty_ppl_block_pcg_tv_fast_cpp <- function(
    data_sim,
    zi1_pattern="^Z1i_", zi2_pattern="^Z2i_",
    max_outer=10, tol=1e-3, beta_maxit=50,
    damping=0.2, verbose=TRUE,
    ties_method=c("breslow","efron"),
    chol_exact_threshold=5000
  ){
    ties_method <- match.arg(ties_method)
    stopifnot(is.list(data_sim), "df" %in% names(data_sim))
    df <- data_sim$df
    n  <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
    
    ops  <- make_linegraph_ops(n, pairs = data_sim$pairs)
    scan <- build_tv_scanline_index(df, zedge_pat="^Zedge_", zi1_pat=zi1_pattern, zi2_pat=zi2_pattern)
    
    X_rows <- scan$X_rows
    
    p_edge <- length(grep("^Zedge_", colnames(X_rows)))
    gamma_names <- colnames(X_rows)[(p_edge+1):ncol(X_rows)]
    p <- ncol(X_rows); m <- max(scan$id_row)
    
    beta  <- rep(0, p)
    delta <- rep(0, m)
    sigma2 <- 0.9
    rho    <- 0.2
    
    S2 <- make_S2(n, pairs = data_sim$pairs)
    history <- vector("list", max_outer)
    
    for (it in 1:max_outer){
      # -- beta update (C++)
      
      # res_b   <- beta_block_update_tv_fast2_cpp(beta, delta, scan, maxit=beta_maxit, ties_method=ties_method)
      # 
      # beta_new <- res_b$beta
      vdelta <- rep(0, m)  # Replace this line if diag(Hdelta^{-1}) is already available
      
    
      res_beta <- beta_update_via_coxph_fit(
        X_rows = scan$X_rows,
        df = df,
        scan =scan,
        delta  = delta,
        vdelta = vdelta,
        ties   = ties_method ,
        strata = rep(1L, nrow(df))   # Alternatively, pass a stratification variable such as df$type
      )
      beta_new <- res_beta$beta

      # -- delta update (FISTA + C++ value/grad)
      iterrrr <-10
      if(it == 1)
      {
        iterrrr <- 20
      }
      res_d <- delta_block_update_fista_tv_fast2_cpp(
        delta, beta_new, scan, ops, sigma2, rho,
        max_iter=iterrrr, tol_grad_rms=1e-4, step_init=5e-1,
        chol_threshold=chol_exact_threshold, pcg_tol=5e-7,
        ties_method=ties_method
      )
      delta_new <- res_d$delta
      
      # -- hyper-params (use full K for stability)
      # lam   <- cumLambda_cpp(scan, beta_new, delta_new, method=ties_method)
      # Kdiag <- K_diag_full_cpp(scan, beta_new, delta_new, lam$cumLambda_by_time)
      Kdiag <- compute_K_diag_partial_cpp(scan,beta_new,delta_new,ties_method = ties_method)
      up <- update_sigma_rho_ratio_fix(
        delta = delta_new, n = n, S2 = S2,
        rho_old = rho, sigma2_old = sigma2,
        K_matvec = function(x) K_matvec_cpp(x, beta_new, delta_new, scan, ties_method),
        K_diag   = Kdiag, damping = damping,
        chol_threshold = chol_exact_threshold, pcg_tol = 5e-7
      )
      sigma2_new <- up$sigma2
      rho_new    <- up$rho
      
      d_beta  <- sqrt(sum((beta_new  - beta )^2))
      d_delta <- sqrt(sum((delta_new - delta)^2))
      d_hyper <- sqrt(sum((c(sigma2_new, rho_new) - c(sigma2, rho))^2))
      
      beta  <- beta_new
      delta <- delta_new
      sigma2 <- sigma2_new
      rho    <- rho_new
      
      beta_edge <- beta[seq_len(p_edge)]
      gamma_vec <- setNames(beta[(p_edge+1):p], gamma_names)
      
      history[[it]] <- list(iter=it, beta_edge=beta_edge, gamma=gamma_vec,
                            sigma=sqrt(sigma2), rho=rho,
                            d_beta=d_beta, d_delta=d_delta, d_hyper=d_hyper)
      
      if (isTRUE(verbose)){
        cat(sprintf("[TV-fast %02d] sigma=%.4g rho=%.4f | dβ=%.2e dδ=%.2e dϑ=%.2e\n",
                    it, sqrt(sigma2), rho, d_beta, d_delta, d_hyper))
        cat("  beta_edge =", paste(round(beta_edge,4), collapse=", "),
            "| gamma =", paste(paste0(names(gamma_vec),"=",round(gamma_vec,4)), collapse=", "), "\n")
      }
      if (max(d_beta, d_delta, d_hyper) < tol) break
    }
    
    list(
      beta_edge = beta_edge,
      gamma     = gamma_vec,
      beta_all  = beta,
      sigma     = sqrt(sigma2),
      rho       = rho,
      delta     = delta,
      trace     = history,
      X_names   = colnames(X_rows)
    )
  }
}


Rcpp::sourceCpp("real_data/RCPP/edge_frailty_core.cpp")
Rcpp::sourceCpp("real_data/RCPP/approx_diag_K.cpp")
# Continue using the existing build_tv_scanline_index() and update_sigma_rho_ratio_fix()
