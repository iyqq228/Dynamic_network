

### Production functions


######## Complete workflow
## A new, faster implementation
{

  sample_edge_delta_undirected <- function(n, sigma2, rho,
                                           method = c("node-factor","sparse-chol")) {
    method <- match.arg(method)
    if (n < 3) stop("n must be >= 3 for undirected dyads.")
    rho_lo <- -1 / (2 * (n - 2))
    rho_hi <-  1 / 2
    if (rho < rho_lo || rho > rho_hi) {
      warning(sprintf("rho=%.4f clipped into [%.4f, %.4f] to ensure PSD.", rho, rho_lo, rho_hi))
      rho <- max(rho_lo, min(rho, rho_hi))
    }
    pairs <- as.matrix(t(combn(n, 2)))
    m <- nrow(pairs)
    
    if (method == "node-factor") {
      if (rho < 0 || rho > 0.5) {
        stop("node-factor 方法仅适用 rho ∈ [0, 0.5]。需要负相关请用 'sparse-chol'。")
      }
      tau2 <- max(rho * sigma2, 0)
      s2   <- max(sigma2 - 2 * tau2, 0)
      u   <- rnorm(n, 0, sqrt(tau2))
      eps <- rnorm(m, 0, sqrt(s2))
      delta <- u[pairs[,1]] + u[pairs[,2]] + eps
      return(list(delta = delta, pairs = pairs))
    }
    

    i_idx <- pairs[,1]; j_idx <- pairs[,2]
    ii <- c(i_idx, j_idx); jj <- rep.int(seq_len(m), 2L)
    B <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, m))
    S2 <- crossprod(B) - Diagonal(m, 2)
    Omega <- sigma2 * (Diagonal(m) + rho * S2)
    

    L <- NULL; ok <- TRUE
    L <- tryCatch(chol(Omega), error = function(e) { ok <<- FALSE; NULL })
    if (!ok) {
      jit <- 1e-10 * mean(diag(Omega))
      L <- tryCatch(chol(Omega + Diagonal(m, jit)),
                    error = function(e) NULL)
      if (is.null(L)) {
        jit <- 1e-6 * mean(diag(Omega))
        L <- chol(Omega + Diagonal(m, jit))
        warning("Applied extra jitter to stabilize sparse Cholesky.")
      } else {
        warning("Applied small jitter to stabilize sparse Cholesky.")
      }
    }
    z <- rnorm(m)
    delta <- as.numeric(t(L) %*% z)
    list(delta = delta, pairs = pairs)
  }
  

  .vapply_f <- function(x, f) {
    out <- tryCatch(f(x), error = function(e) NULL)
    if (is.null(out) || length(out) != length(x)) {
      return(vapply(x, f, numeric(1)))
    }
    out
  }
  

  make_z_ij_vec <- function(p_edge_eff, prob = 0.5, mu = 0, sd = 1) {
    stopifnot(p_edge_eff >= 1)
    if (p_edge_eff == 1) {
      return(rbinom(1, 1, prob))
    }

    c(rnorm(p_edge_eff - 1, mean = mu, sd = sd), rbinom(1, 1, prob)-prob)
  }
  
  simulate_joint_network_data_tv <- function(
    n = 50,
    T_max = 10,
    beta = 1,
    eta  = 1,
    h0   = 0.1,
    seed = 123,
    alpha_true = NULL,
    edge_frailty = FALSE,
    Z_i = NULL,
    Z_edge_gen = c("avg","rand"),
    p_node = NULL,

    sigma_edge = 0.01,
    rho_edge   = 0,
    edge_method = c("node-factor","sparse-chol"),

    tv_breaks   = 4,
    tv_edge_fun = NULL,
    tv_node_fun = NULL,
    tv_edge_amp = 0.4,
    tv_node_amp = 0.4,

    ties = FALSE,
    ties_bins = 100
  ){
    set.seed(seed)
    

    .vapply_f <- function(x, f) {
      out <- tryCatch(f(x), error = function(e) NULL)
      if (!is.null(out) && is.numeric(out) && length(out) == length(x)) return(out)
      vapply(x, f, numeric(1))
    }
    

    resolve_ties_times <- function(ev_times, t_censor, dt, left_only = TRUE) {
      if (!length(ev_times)) return(ev_times)
      ev_times <- sort(ev_times)
      if (!anyDuplicated(ev_times)) return(ev_times)
      eps <- max(dt * 1e-4, 1e-6)
      r <- rle(ev_times)
      idx0 <- cumsum(c(1L, head(r$lengths, -1L)))
      for (b in which(r$lengths > 1L)) {
        i0 <- idx0[b]; L <- r$lengths[b]; base <- ev_times[i0]
        if (left_only) {
          offs <- eps * seq_len(L-1L)
          placed <- pmax(base - rev(offs), 0)
          prev <- if (i0 == 1L) 0 else ev_times[i0 - 1L]
          placed <- pmax(placed, prev + 1e-12)
          ev_times[i0:(i0+L-2L)] <- placed
          ev_times[i0+L-1L] <- base
        } else {
          kL <- floor((L-1L)/2); kR <- (L-1L) - kL
          lefts  <- base - eps * seq_len(kL)
          rights <- base + eps * seq_len(kR)
          arr <- sort(c(lefts, base, rights))
          arr[arr < 0]        <- 0
          arr[arr > t_censor] <- t_censor
          ev_times[i0:(i0+L-1L)] <- arr
        }
      }
      sort(ev_times)
    }

    
    Z_edge_gen  <- match.arg(Z_edge_gen)
    edge_method <- match.arg(edge_method)
    
    beta_vec <- as.numeric(beta);  p_edge_eff <- length(beta_vec)
    

    use_node_terms <- !is.null(eta)
    
    if (use_node_terms) {
      eta_vec <- as.numeric(eta)

      if (is.null(Z_i)) {
        q <- if (is.null(p_node)) length(eta_vec) else as.integer(p_node)
        if (q <= 0) stop("eta 非 NULL 时，节点维度 q 必须 > 0。可通过 p_node 指定。")
        Z_i <- matrix(rnorm(n * q, 0, 0.8), nrow = n, ncol = q)
      } else {
        q <- ncol(Z_i)
        if (q != length(eta_vec))
          stop("Z_i 的列数必须与 eta 的长度一致。")
      }
      node_effc <- as.vector(Z_i %*% eta_vec)
      z1i_cols  <- paste0("Z1i_",  seq_len(q))
      z2i_cols  <- paste0("Z2i_",  seq_len(q))
    } else {

      q <- 0
      Z_i <- NULL
      node_effc <- rep(0, n)
      z1i_cols <- z2i_cols <- character(0)
    }

    
    z_cols <- paste0("Zedge_", seq_len(p_edge_eff))
    

    mu_i    <- if (!is.null(alpha_true)) as.numeric(alpha_true) else rep(0, n)
    gamma_i <- rep(0, n)
    

    pairs_mat <- as.matrix(t(combn(n, 2)))
    m_edges   <- nrow(pairs_mat)
    

    delta_vec <- NULL
    if (edge_frailty) {
      samp <- sample_edge_delta_undirected(
        n = n, sigma2 = sigma_edge^2, rho = rho_edge, method = edge_method
      )
      delta_vec <- samp$delta
    }
    

    lambda0_fun <- function(t) 2 + 0.5 * sin(t)
    lambda0_max <- 10
    

    if (is.null(tv_edge_fun)) {
      a_e <- tv_edge_amp
      tv_edge_fun <- function(t) { ifelse(t > 1, 1, 0.5) }


      f_e_min <- 0.1; f_e_max <- 2
    } else { f_e_min <- 0.1; f_e_max <- 2 }
    
    if (use_node_terms) {
      if (is.null(tv_node_fun)) {
        a_n <- tv_node_amp
        tv_node_fun <- function(t) { ifelse(t > 1, 1, 0.5) }
        f_n_min <- 1 - abs(a_n); f_n_max <- 1 + abs(a_n)
      } else { f_n_min <- 0.1; f_n_max <- 2 }
      if (f_n_min <= 0) stop("时变节点缩放必须为正；请调小 tv_node_amp 或自定义正函数。")
    } else {

      f_n_min <- f_n_max <- 1
    }
    if (f_e_min <= 0) stop("时变边缩放必须为正；请调小 tv_edge_amp 或自定义正函数。")
    

    rate_node <- h0 * exp(gamma_i)
    rate_node[rate_node <= 0] <- .Machine$double.eps
    death_time <- rexp(n, rate = rate_node) + T_max/2
    

    if (ties) {
      dt <- T_max / ties_bins
      snap_to_grid <- function(x) pmin(T_max, pmax(0, round(x / dt) * dt))
    } else {
      snap_to_grid <- function(x) x
      dt <- T_max
    }
    

    ccc_i <- matrix(rnorm(n * 2, 0, 0.8), nrow = n, ncol = 2)
    
    rows_list <- vector("list", m_edges)
    for (e_id in seq_len(m_edges)) {
      i <- pairs_mat[e_id, 1]; j <- pairs_mat[e_id, 2]
      t_censor <- min(death_time[i], death_time[j], T_max)
      if (t_censor <= 0) next
      

      if (Z_edge_gen == "avg") {
        z_ij_vec <- (ccc_i[i, ] + ccc_i[j, ]) / 2





      } else {
        z_ij_vec <- make_z_ij_vec(p_edge_eff)
      }
     
      z_edge_lin_base <- sum(beta_vec * z_ij_vec)
      
      node_pair_base  <- if (use_node_terms) (node_effc[i] + node_effc[j]) else 0
      delta_ij <- if (edge_frailty) delta_vec[e_id] else 0
      

      edge_sup <- if (z_edge_lin_base >= 0) f_e_max * z_edge_lin_base else f_e_min * z_edge_lin_base
      if (use_node_terms) {
        node_sup <- if (node_pair_base  >= 0) f_n_max * node_pair_base  else f_n_min * node_pair_base
      } else {
        node_sup <- 0
      }
      rate_max <- lambda0_max * exp(edge_sup + node_sup + mu_i[i] + mu_i[j] + delta_ij)
      if (!is.finite(rate_max) || rate_max <= 0) next
      

      K <- rpois(1L, lambda = rate_max * t_censor)
      ev_times <- numeric(0)
      if (K > 0) {
        t_prop <- sort(runif(K, 0, t_censor))
        fe <- .vapply_f(t_prop, tv_edge_fun)
        if (use_node_terms) fn <- .vapply_f(t_prop, tv_node_fun)
        
        lin_t <- fe * z_edge_lin_base + mu_i[i] + mu_i[j] + delta_ij
        if (use_node_terms) lin_t <- lin_t + fn * node_pair_base
        
        lam_t <- exp(lin_t) * .vapply_f(t_prop, lambda0_fun)
        u <- runif(K)
        keep <- (u <= pmin(pmax(lam_t / rate_max, 0), 1))
        if (any(keep)) ev_times <- t_prop[keep]
      }
      

      if (length(ev_times)) ev_times <- snap_to_grid(ev_times)
      if (length(ev_times) && ties) {
        ev_times <- resolve_ties_times(ev_times, t_censor, dt, left_only = TRUE)
      }
      

      if (length(tv_breaks) == 1L) {
        Kb <- max(as.integer(tv_breaks), 1L)
        base_grid <- seq(0, t_censor, length.out = Kb + 1L)
      } else {
        probs <- sort(pmin(pmax(as.numeric(tv_breaks), 0), 1))
        base_grid <- sort(unique(c(0, probs * t_censor, t_censor)))
      }
      cuts <- sort(unique(c(base_grid, ev_times)))
      if (length(cuts) < 2L) cuts <- c(0, t_censor)
      
      s <- cuts[-length(cuts)]
      e <- cuts[-1L]
      n_int <- length(s)
      t_mid <- 0.5 * (s + e)
      
      fe_mid <- .vapply_f(t_mid, tv_edge_fun)
      if (use_node_terms) fn_mid <- .vapply_f(t_mid, tv_node_fun)
      

      if (length(ev_times)) {
        M <- abs(outer(e, ev_times, "-")) < 1e-12
        event_vec <- as.integer(rowSums(M) > 0L)
      } else {
        event_vec <- integer(n_int)
      }
      

      zedge_mat <- if (p_edge_eff > 0) {
        out <- matrix(rep(z_ij_vec, each = n_int), nrow = n_int) * fe_mid
        colnames(out) <- z_cols
        out
      } else NULL
      
      if (use_node_terms) {
        z1_mat <- matrix(rep(Z_i[i, ], each = n_int), nrow = n_int) * fn_mid
        z2_mat <- matrix(rep(Z_i[j, ], each = n_int), nrow = n_int) * fn_mid
        colnames(z1_mat) <- z1i_cols
        colnames(z2_mat) <- z2i_cols
      } else {
        z1_mat <- z2_mat <- NULL
      }
      
      const_df <- data.frame(
        start = s, stop = e, event = event_vec,
        id = rep.int(e_id, n_int),
        node1 = rep.int(i, n_int),
        node2 = rep.int(j, n_int),
        type  = ifelse(seq_len(n_int) == n_int, "death", "interval"),
        alpha_i = rep.int(mu_i[i], n_int),
        alpha_j = rep.int(mu_i[j], n_int),
        delta   = rep.int(delta_ij, n_int),
        check.names = FALSE
      )
      
      df_edge <- cbind(
        const_df,
        if (!is.null(zedge_mat)) as.data.frame(zedge_mat, check.names = FALSE) else NULL,
        if (!is.null(z1_mat))    as.data.frame(z1_mat,    check.names = FALSE) else NULL,
        if (!is.null(z2_mat))    as.data.frame(z2_mat,    check.names = FALSE) else NULL
      )
      rows_list[[e_id]] <- df_edge
    }
    

    df <- do.call(rbind, rows_list)
    if (!is.null(df) && nrow(df)) {
      df <- df[order(df$stop, df$start, df$id), , drop = FALSE]
      rownames(df) <- NULL
      df$time  <- df$stop
      df$death <- as.integer(df$type == "death")
    } else {
      df <- data.frame()
    }
    
    list(
      df          = df,
      mu_i        = mu_i,
      gamma_i     = gamma_i,
      Z_i         = if (use_node_terms) Z_i else NULL,
      delta_vec   = if (edge_frailty) delta_vec else NULL,
      pairs       = pairs_mat
    )
  }
  
  
  
  
  
}

##### Estimation part I
{
  suppressPackageStartupMessages({ library(Matrix) })
  library(survival)
  make_S2 <- function(n, pairs) {

    m <- nrow(pairs)
    ii <- c(pairs[,1], pairs[,2])
    jj <- rep.int(seq_len(m), 2L)
    B  <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, m))
    crossprod(B) - Diagonal(m, 2)
  }
  

  pcg_solve <- function(matvec, b, Mdiag = NULL, tol = 1e-6, maxit = 2000) {
    n <- length(b); x <- numeric(n)
    r <- b - matvec(x)
    z <- if (is.null(Mdiag)) r else r / Mdiag
    p <- z; rz_old <- sum(r * z)
    if (sqrt(rz_old) < tol) return(x)
    for (k in 1:maxit) {
      Ap <- matvec(p)
      denom <- sum(p * Ap); if (abs(denom) < 1e-30) break
      alpha <- rz_old / denom
      x <- x + alpha * p
      r <- r - alpha * Ap
      z <- if (is.null(Mdiag)) r else r / Mdiag
      rz_new <- sum(r * z)
      if (sqrt(rz_new) < tol) break
      beta <- rz_new / rz_old
      p <- z + beta * p
      rz_old <- rz_new
    }
    x
  }
  

  slq_logdet <- function(matvec, n, n_vec = 20, m_lanczos = 30, seed = 1L) {
    set.seed(seed)
    ests <- numeric(n_vec)
    for (k in 1:n_vec) {
      v <- sample(c(-1, 1), n, replace = TRUE) / sqrt(n)
      al <- numeric(m_lanczos); be <- numeric(m_lanczos - 1L)
      q_old <- rep(0, n); q <- v; beta <- 0; jlen <- 0L
      for (j in 1:m_lanczos) {
        z <- matvec(q); if (j > 1) z <- z - beta * q_old
        alpha <- sum(q * z); z <- z - alpha * q; beta <- sqrt(sum(z * z))
        al[j] <- alpha; if (j < m_lanczos) be[j] <- beta
        q_old <- q; jlen <- j
        if (beta < 1e-14) break
        q <- z / beta
      }
      Tj <- diag(al[1:jlen])
      if (jlen > 1L) {
        idx <- 1:(jlen - 1L)
        Tj[cbind(idx, idx + 1L)] <- be[1:(jlen - 1L)]
        Tj[cbind(idx + 1L, idx)] <- be[1:(jlen - 1L)]
      }
      ev <- eigen(Tj, symmetric = TRUE)
      logT11 <- sum((ev$vectors[1, ]^2) * log(pmax(ev$values, .Machine$double.eps)))
      ests[k] <- logT11
    }
    n * mean(ests)
  }
  

  
  make_linegraph_ops <- function(n, pairs) {
    m <- nrow(pairs)
    ii <- c(pairs[,1], pairs[,2])
    jj <- rep.int(seq_len(m), 2L)
    B  <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, m))
    S2_matvec <- function(x) as.numeric(crossprod(B, B %*% x) - 2 * x)
    S2_matmat <- function(X) as.matrix(crossprod(B, B %*% X) - 2 * X)
    M_matvec  <- function(rho) {
      if (abs(rho) < .Machine$double.eps) return(function(x) x)
      function(x) x + rho * (as.numeric(crossprod(B, B %*% x)) - 2 * x)
    }
    list(B = B, S2_matvec = S2_matvec, S2_matmat = S2_matmat, M_matvec = M_matvec)
  }
  
  make_probe_matrix <- function(m, n_probe = 16, seed = 1L) {
    set.seed(seed)
    V <- matrix(sample(c(-1L,1L), m * n_probe, replace = TRUE), m, n_probe)
    storage.mode(V) <- "double"; V
  }
  


  make_M_solver <- function(ops, rho, chol_threshold = 5000L,
                            pcg_tol = 1e-6, pcg_maxit = 2000) {
    m <- ncol(ops$B)
    if (m <= chol_threshold) {
      M <- Diagonal(m) + rho * (crossprod(ops$B) - Diagonal(m, 2))
      ok <- TRUE
      cf <- tryCatch(Cholesky(M, LDL = FALSE, Imult = 0, perm = TRUE),
                     error = function(e) { ok <<- FALSE; NULL })
      if (!ok) cf <- Cholesky(M + Diagonal(m, 1e-10), LDL = FALSE, Imult = 0.0, perm = TRUE)
      solvefun <- function(b) as.numeric(solve(cf, b))
    } else {
      mv <- ops$M_matvec(rho)




      Mdiag <- NULL
      solvefun <- function(b) pcg_solve(mv, b, Mdiag = Mdiag, tol = pcg_tol, maxit = pcg_maxit)
    }
    solvefun
  }
  


  approx_diag_K <- function(scan, beta, delta, ties_method = c("breslow","efron")) {
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
      while (si <= nR && ss[si] < t) {
        r <- rs[si]; eid <- id[r]; wr <- a_row[r] * d_id[eid]
        S0 <- S0 + wr; Wid[eid] <- Wid[eid] + wr; si <- si + 1L
      }
      while (ei <= nR && es[ei] < t) {
        r <- re[ei]; eid <- id[r]; wr <- a_row[r] * d_id[eid]
        S0 <- S0 - wr; Wid[eid] <- Wid[eid] - wr; ei <- ei + 1L
      }
      if (!(S0 > 0) || !is.finite(S0)) next
      
      mk <- mt_by_time[k]
      if (mk == 0L) next
      
      if (ties_method == "breslow") {
        w <- Wid / S0
        diagK <- diagK + mk * (w * (1 - w))
      } else {
        ev_rows <- grp[[k]]
        ev_eids <- id[ev_rows]
        wevent <- numeric(m)
        if (length(ev_rows)) {
          wr_ev <- a_row[ev_rows] * d_id[id[ev_rows]]
          acc   <- rowsum(wr_ev, group = ev_eids, reorder = FALSE)
          wevent[as.integer(rownames(acc))] <- as.numeric(acc)
        }
        Ek <- sum(wevent)
        if (Ek <= 0 || !is.finite(Ek)) {
          w <- Wid / S0
          diagK <- diagK + mk * (w * (1 - w))
        } else {
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (denom <= 0 || !is.finite(denom)) next
            w_ell <- (Wid - alpha * wevent) / denom
            diagK <- diagK + w_ell * (1 - w_ell)
          }
        }
      }
    }
    
    
    pmax(diagK, 0)
  }
  
  








































  
  

  compute_K_diag_partial <- function(scan, beta, delta, ties_method = c("breslow","efron")) {
    ties_method <- match.arg(ties_method)
    

    X   <- scan$X_rows
    id  <- as.integer(scan$id_row)
    m   <- length(delta)
    nR  <- nrow(X)
    tt  <- scan$times
    K   <- length(tt)
    

    uniq_id <- sort(unique(id))
    if (!identical(uniq_id, seq_len(length(uniq_id)))) {
      map <- match(id, uniq_id)
      id  <- as.integer(map)

      m   <- max(id)

      delta <- delta[match(uniq_id, uniq_id)]
    }
    

    a_row <- exp(as.vector(X %*% beta))
    d_id  <- exp(delta)
    

    ss <- scan$start_sorted
    es <- scan$stop_sorted
    rs <- scan$rows_start
    re <- scan$rows_stop
    grp <- scan$events_rowidx_by_time
    mt_by_time <- vapply(grp, length, 0L)
    

    Wid   <- numeric(m)
    S0    <- 0.0
    si    <- 1L
    ei    <- 1L
    diagK <- numeric(m)
    

    for (k in seq_len(K)) {
      t <- tt[k]
      

      while (si <= nR && ss[si] < t) {
        r   <- rs[si]
        eid <- id[r]
        wr  <- a_row[r] * d_id[eid]
        S0      <- S0 + wr
        Wid[eid] <- Wid[eid] + wr
        si <- si + 1L
      }
      

      while (ei <= nR && es[ei] < t) {
        r   <- re[ei]
        eid <- id[r]
        wr  <- a_row[r] * d_id[eid]
        S0      <- S0 - wr
        Wid[eid] <- Wid[eid] - wr
        ei <- ei + 1L
      }
      
      mk <- mt_by_time[k]
      if (mk == 0L || !(S0 > 0) || !is.finite(S0)) next
      
      if (ties_method == "breslow") {

        w <- Wid / S0

        diagK <- diagK + mk * (w * (1 - w))
        
      } else {

        ev_rows <- grp[[k]]

        if (length(ev_rows)) {
          wr_ev <- a_row[ev_rows] * d_id[id[ev_rows]]
          acc   <- rowsum(wr_ev, group = id[ev_rows], reorder = FALSE)
          wevent <- numeric(m)
          wevent[as.integer(rownames(acc))] <- as.numeric(acc)
        } else {
          wevent <- numeric(m)
        }
        Ek <- sum(wevent)
        
        if (!(Ek > 0) || !is.finite(Ek)) {

          w <- Wid / S0
          diagK <- diagK + mk * (w * (1 - w))
        } else {

          for (ell in 0:(mk - 1L)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (!(denom > 0) || !is.finite(denom)) next
            w_ell <- (Wid - alpha * wevent) / denom
            diagK <- diagK + (w_ell * (1 - w_ell))
          }
        }
      }
    }
    

    diagK[diagK < 0] <- 0
    diagK
  }
  
  
  

  build_tv_scanline_index <- function(df,
                                      zedge_pat="^Zedge_",
                                      zi1_pat="^Z1i_",
                                      zi2_pat="^Z2i_",
                                      scan_node = FALSE) {
    need <- c("id","start","stop","event")
    if (!all(need %in% names(df))) stop("df 需要列: id,start,stop,event")
    


    id_row <- as.integer(df$id)
    m <- max(id_row)
    

    zcols <- grep(zedge_pat, names(df), value = TRUE)
    if (!length(zcols)) stop("未找到 Zedge_*")
    

    zi1 <- zi2 <- character(0)
    common <- character(0)
    if (isTRUE(scan_node)) {
      zi1 <- grep(zi1_pat, names(df), value = TRUE)
      zi2 <- grep(zi2_pat, names(df), value = TRUE)
      if (!length(zi1) || !length(zi2)) stop("未找到 Z1i_*/Z2i_*")
      suf <- function(v, pat) sub(pat, "", v, perl=TRUE)
      s1 <- suf(zi1, zi1_pat); s2 <- suf(zi2, zi2_pat)
      common <- intersect(s1, s2)
      if (!length(common)) stop("Z1i_* 与 Z2i_* 无共同后缀")
      zi1 <- zi1[match(common, s1)]
      zi2 <- zi2[match(common, s2)]
    }
    

    if (isTRUE(scan_node)) {
      X_rows <- cbind(as.matrix(df[, zcols, drop=FALSE]),
                      as.matrix(df[, zi1,  drop=FALSE]) + as.matrix(df[, zi2, drop=FALSE]))
      colnames(X_rows) <- c(zcols, paste0("Xi_plus_Xj_", common))
    } else {
      X_rows <- as.matrix(df[, zcols, drop=FALSE])
      colnames(X_rows) <- zcols
    }
    storage.mode(X_rows) <- "double"
    

    idx_event <- which(df$event == 1L)
    if (!length(idx_event)) stop("没有事件行（event==1）")
    times <- sort(unique(df$stop[idx_event]))
    grp   <- split(idx_event, f = match(df$stop[idx_event], times))
    

    ord_start <- order(df$start, df$id)
    ord_stop  <- order(df$stop,  df$id)
    
    nR <- nrow(df)
    start_by_row <- numeric(nR); stop_by_row <- numeric(nR)
    start_by_row[ord_start] <- df$start[ord_start]
    stop_by_row [ord_stop ] <- df$stop [ord_stop ]
    

    k_enter <- findInterval(start_by_row, times) + 1L
    k_leave <- findInterval(stop_by_row,  times) + 1L
    
    K <- length(times)
    idx_enter <- which(k_enter >= 1L & k_enter <= K)
    idx_leave <- which(k_leave >= 1L & k_leave <= K)
    
    rows_enter_by_time <- split(idx_enter, k_enter[idx_enter])
    rows_leave_by_time <- split(idx_leave, k_leave[idx_leave])
    

    list(
      times = times,
      events_rowidx_by_time = grp,
      start_sorted = df$start[ord_start],
      stop_sorted  = df$stop[ord_stop],
      rows_start   = ord_start,
      rows_stop    = ord_stop,
      id_row = id_row,
      X_rows = X_rows,
      mt_by_time = vapply(grp, length, 0L),
      rows_enter_by_time = rows_enter_by_time,
      rows_leave_by_time = rows_leave_by_time,

      used_node_terms = isTRUE(scan_node),
      colnames_edge = zcols,
      colnames_node = if (isTRUE(scan_node)) paste0("Xi_plus_Xj_", common) else character(0)
    )
  }
  

  make_scan_kernel <- function(scan, beta, ties_method = c("breslow","efron")) {
    ties_method <- match.arg(ties_method)
    X  <- scan$X_rows
    id <- as.integer(scan$id_row)
    a_row <- exp(as.vector(X %*% beta))
    grp <- scan$events_rowidx_by_time
    

    xbeta_all <- as.vector(X %*% beta)
    xbeta_sum <- vapply(grp, function(idx) sum(xbeta_all[idx]), 0.0)
    
    K  <- length(scan$times)
    mt <- if (!is.null(scan$mt_by_time)) scan$mt_by_time else vapply(grp, length, 0L)
    

    value_grad_delta_fastR <- function(delta) {
      d_id <- exp(delta)
      m    <- length(delta)
      Wid  <- numeric(m)
      S0   <- 0.0
      val  <- 0.0
      g    <- numeric(m)
      
      for (k in seq_len(K)) {

        add_rows <- scan$rows_enter_by_time[[as.character(k)]]
        if (!is.null(add_rows) && length(add_rows)) {
          eids <- id[add_rows]
          wr   <- a_row[add_rows] * d_id[eids]
          S0   <- S0 + sum(wr)
          acc  <- rowsum(wr, group = eids, reorder = FALSE)
          idx  <- as.integer(rownames(acc))
          Wid[idx] <- Wid[idx] + as.numeric(acc)
        }
        

        rm_rows <- scan$rows_leave_by_time[[as.character(k)]]
        if (!is.null(rm_rows) && length(rm_rows)) {
          eids <- id[rm_rows]
          wr   <- a_row[rm_rows] * d_id[eids]
          S0   <- S0 - sum(wr)
          acc  <- rowsum(wr, group = eids, reorder = FALSE)
          idx  <- as.integer(rownames(acc))
          Wid[idx] <- Wid[idx] - as.numeric(acc)
        }
        
        mk <- mt[k]
        if (mk == 0L || !(S0 > 0) || !is.finite(S0)) next
        
        ev_rows <- grp[[k]]
        ev_eids <- id[ev_rows]
        

        val <- val + xbeta_sum[k] + sum(delta[ev_eids])
        

        g <- g + tabulate(ev_eids, nbins = m)
        

        if (ties_method == "breslow" || mk == 1L) {
          val <- val - mk * log(S0)
          g   <- g - (mk / S0) * Wid
        } else {

          wr_ev <- a_row[ev_rows] * d_id[ev_eids]
          if (length(wr_ev)) {
            acc  <- rowsum(wr_ev, group = ev_eids, reorder = FALSE)
            we_idx <- as.integer(rownames(acc))
            we_val <- as.numeric(acc)
            Ek     <- sum(we_val)
          } else {
            we_idx <- integer(0); we_val <- numeric(0); Ek <- 0
          }
          
          if (!(Ek > 0) || !is.finite(Ek)) {

            val <- val - mk * log(S0)
            g   <- g - (mk / S0) * Wid
          } else {
            alpha_vec <- (0:(mk-1)) / mk
            denom_vec <- S0 - alpha_vec * Ek
            good <- is.finite(denom_vec) & (denom_vec > 0)
            if (!any(good)) next
            denom_vec <- denom_vec[good]
            alpha_vec <- alpha_vec[good]
            

            val <- val - sum(log(denom_vec))

            s1 <- sum(1 / denom_vec)
            s2 <- sum(alpha_vec / denom_vec)
            g  <- g - s1 * Wid
            if (length(we_idx)) g[we_idx] <- g[we_idx] + s2 * we_val
          }
        }
      }
      list(value = val, grad_delta = g)
    }
    
    list(
      value_grad_delta = value_grad_delta_fastR,

      K_matvec = function(x, delta) {
        d_id <- exp(delta)
        ss <- scan$start_sorted; es <- scan$stop_sorted
        rs <- scan$rows_start;    re <- scan$rows_stop
        times <- scan$times; nR <- nrow(X); m <- length(delta)
        out <- numeric(m); Wid <- numeric(m); si <- 1L; ei <- 1L; S0 <- 0.0
        mt_by_time <- vapply(grp, length, 0L)
        for (k in seq_along(times)) {
          t <- times[k]
          while (si <= nR && ss[si] < t) { r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0+wr; Wid[eid] <- Wid[eid]+wr; si <- si+1L }
          while (ei <= nR && es[ei] <  t) { r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0-wr; Wid[eid] <- Wid[eid]-wr; ei <- ei+1L }
          mk <- mt_by_time[k]
          if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
          if (ties_method == "breslow") {
            w <- Wid / S0; dot <- sum(w * x); out <- out + mk * w * (x - dot)
          } else {
            ev_rows <- grp[[k]]
            wevent <- numeric(m)
            if (length(ev_rows)) {
              wr_ev <- a_row[ev_rows] * d_id[id[ev_rows]]
              acc <- rowsum(wr_ev, group = id[ev_rows], reorder = FALSE)
              wevent[as.integer(rownames(acc))] <- as.numeric(acc)
            }
            Ek <- sum(wevent)
            if (!(Ek > 0) || !is.finite(Ek)) {
              w <- Wid / S0; dot <- sum(w * x); out <- out + mk * w * (x - dot)
            } else {
              for (ell in 0:(mk-1)) {
                alpha <- ell / mk; denom <- S0 - alpha * Ek
                if (denom <= 0 || !is.finite(denom)) next
                w_ell <- (Wid - alpha * wevent) / denom
                dot <- sum(w_ell * x)
                out <- out + w_ell * (x - dot)
              }
            }
          }
        }
        out
      },
      K_matmul = function(Xin, delta) {

        Xmat <- if (is.matrix(Xin)) Xin else cbind(Xin)
        b <- ncol(Xmat)
        d_id <- exp(delta)
        ss <- scan$start_sorted; es <- scan$stop_sorted
        rs <- scan$rows_start;    re <- scan$rows_stop
        times <- scan$times; nR <- nrow(X); m <- length(delta)
        out <- matrix(0.0, m, b)
        Wid <- numeric(m); si <- 1L; ei <- 1L; S0 <- 0.0
        mt_by_time <- vapply(grp, length, 0L)
        for (k in seq_along(times)) {
          t <- times[k]
          while (si <= nR && ss[si] < t) { r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0+wr; Wid[eid] <- Wid[eid]+wr; si <- si+1L }
          while (ei <= nR && es[ei] <  t) { r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]; S0 <- S0-wr; Wid[eid] <- Wid[eid]-wr; ei <- ei+1L }
          mk <- mt_by_time[k]
          if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
          if (ties_method == "breslow") {
            w <- Wid / S0; WX <- Xmat * w; dots <- crossprod(w, Xmat); out <- out + mk * (WX - (w %*% dots))
          } else {
            ev_rows <- grp[[k]]
            wevent <- numeric(m)
            if (length(ev_rows)) {
              wr_ev <- a_row[ev_rows] * d_id[id[ev_rows]]
              acc <- rowsum(wr_ev, group = id[ev_rows], reorder = FALSE)
              wevent[as.integer(rownames(acc))] <- as.numeric(acc)
            }
            Ek <- sum(wevent)
            if (!(Ek > 0) || !is.finite(Ek)) {
              w <- Wid / S0; WX <- Xmat * w; dots <- crossprod(w, Xmat); out <- out + mk * (WX - (w %*% dots))
            } else {
              for (ell in 0:(mk-1)) {
                alpha <- ell / mk; denom <- S0 - alpha * Ek
                if (denom <= 0 || !is.finite(denom)) next
                w <- (Wid - alpha * wevent) / denom
                WX <- Xmat * w; dots <- crossprod(w, Xmat); out <- out + (WX - (w %*% dots))
              }
            }
          }
        }
        if (is.matrix(Xin)) out else as.vector(out)
      }
    )
  }
  

  
  compute_cumLambda <- function(scan, beta, delta, method = c("breslow","efron")) {
    method <- match.arg(method)
    X  <- scan$X_rows
    id <- as.integer(scan$id_row)
    m  <- length(delta)
    stopifnot(max(id) <= m)
    
    a_row <- exp(as.vector(X %*% beta))
    d_id  <- exp(delta)
    
    ss <- scan$start_sorted; es <- scan$stop_sorted
    rs <- scan$rows_start;    re <- scan$rows_stop
    times <- scan$times; nR <- nrow(X)
    
    mt_by_time <- vapply(scan$events_rowidx_by_time, length, 0L)
    
    Wid <- numeric(m)
    si <- 1L; ei <- 1L; S0 <- 0.0
    dLambda <- numeric(length(times))
    S0_by_time <- numeric(length(times))
    event_risk_by_time <- numeric(length(times))
    
    for (k in seq_along(times)) {
      t <- times[k]
      while (si <= nR && ss[si] < t) {
        r   <- rs[si]; eid <- id[r]
        wr  <- a_row[r] * d_id[eid]
        S0  <- S0 + wr
        Wid[eid] <- Wid[eid] + wr
        si  <- si + 1L
      }
      while (ei <= nR && es[ei] < t) {
        r   <- re[ei]; eid <- id[r]
        wr  <- a_row[r] * d_id[eid]
        S0  <- S0 - wr
        Wid[eid] <- Wid[eid] - wr
        ei  <- ei + 1L
      }
      
      S0_by_time[k] <- S0
      mk <- mt_by_time[k]
      if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
      
      if (method == "efron") {
        ev_rows <- scan$events_rowidx_by_time[[k]]
        if (length(ev_rows)) {
          event_risk_by_time[k] <- sum(a_row[ev_rows] * d_id[id[ev_rows]])
        } else {
          event_risk_by_time[k] <- 0
        }
        Ek <- event_risk_by_time[k]
        dLambda[k] <- sum(1 / (S0 - (0:(mk-1)) * (Ek / mk)))
      } else {
        dLambda[k] <- mk / S0
        event_risk_by_time[k] <- NA_real_
      }
    }
    
    list(
      times = times,
      mt_by_time = mt_by_time,
      S0_by_time = S0_by_time,
      dLambda_by_time = dLambda,
      cumLambda_by_time = cumsum(dLambda),
      event_risk_by_time = event_risk_by_time
    )
  }
  













































































  
}



##### Estimation part II
{
  beta_update_via_coxph_fit <- function(X_rows, df,
                                        scan,
                                        delta, vdelta = NULL,
                                        beta0 = NULL,
                                        sigma2 ,
                                        ties = c("efron","breslow"),
                                        strata = NULL,
                                        weights = NULL,
                                        maxit = 50) {
    
    
    row_df <- cbind(df[,c("start","stop","event")], scan$X_rows)
    id_row <- scan$id_row
    uids <- sort(unique(id_row))
    if (!identical(uids, seq_len(length(uids)))) id_row <- match(id_row, uids)
    m <- max(id_row); if (!exists("vdelta")) vdelta <- rep(0, m)
    row_df$off <- delta[id_row]
    


    form <- as.formula(
      paste("Surv(start, stop, event) ~",
            paste(colnames(scan$X_rows), collapse = " + "),
            "+ offset(off)")
    )
    fit_off <- coxph(form, data=row_df, ties="efron",
                     model=FALSE, x=FALSE, y=FALSE,
                     control = coxph.control(timefix = FALSE))
    beta_off <- coef(fit_off)
    
    list(
      beta   = as.numeric(beta_off)
    )
  }
  
  
  






























































































































































  {
    {
      
      delta_block_update_fista_tv_fast2 <- function(
    delta0, beta, scan, ops, sigma2, rho,
    max_iter = 20, tol_grad_rms = 1e-4,
    backtrack = TRUE, step_init = NULL, bt_factor = 0.5, ls_armijo = 1e-4,
    grad_clip_rms = 1e3,
    chol_threshold = 25000L, pcg_tol = 1e-6, pcg_maxit = 2000,
    ties_method = c("breslow","efron")
      ){

        ties_method <- match.arg(ties_method)
        m <- length(delta0)
        if (!is.finite(sigma2) || sigma2 <= 0) stop("sigma2 must be > 0")
        if (!is.numeric(rho)) stop("rho must be numeric")
        if (length(beta) == 0) stop("beta must be provided")
        if (length(delta0) == 0) stop("delta0 must be non-empty")
        

        kern <- make_scan_kernel(scan, beta, ties_method = ties_method)
        Minv <- make_M_solver(ops, rho, chol_threshold, pcg_tol, pcg_maxit)
        

        K_diag <- compute_K_diag_partial(scan, beta, delta0, ties_method = ties_method)
        K_diag <- as.numeric(K_diag)
        epsK   <- 1e-12
        


        A_diag_ub <- rep(1, m)
        if (!is.null(ops$S2)) {

          degL <- tryCatch(as.numeric(Matrix::rowSums(ops$S2)), error = function(e) NULL)
          if (!is.null(degL)) {

            A_diag_ub <- 1 / pmax(1 + rho * (-abs(degL)), 1e-12)
          }
        }
        Hdiag   <- pmax(K_diag, epsK) + (A_diag_ub / sigma2)
        invHdiag <- 1 / Hdiag
        

        Hscale <- median(Hdiag)
        if (!is.finite(Hscale) || Hscale <= 0) Hscale <- mean(Hdiag)
        if (!is.finite(Hscale) || Hscale <= 0) Hscale <- 1
        step <- if (is.null(step_init)) 1 / Hscale else (step_init / Hscale)
        step <- min(max(step, 1e-6), 1e6)
        

        eval_f_g <- function(delta) {
          core <- kern$value_grad_delta(delta)

          Ad <- Minv(delta)

          f <- core$value - 0.5 * sum(delta * Ad) / sigma2

          g <- core$grad_delta - Ad / sigma2
          

          g_rms <- sqrt(mean(g^2))
          if (g_rms > grad_clip_rms) g <- g * (grad_clip_rms / g_rms)
          list(f = f, g = g)
        }
        

        y      <- delta0
        delta  <- delta0
        t_k    <- 1.0
        
        eg0    <- eval_f_g(delta0)
        f_prev <- eg0$f
        if (!is.finite(f_prev)) stop("Initial objective is not finite. Check inputs.")

        g0_rms <- sqrt(mean(eg0$g^2))
        if (!is.finite(g0_rms) || g0_rms < tol_grad_rms) {
          return(list(delta = delta0, value = f_prev, iters = 0L, step = step, grad_rms = g0_rms))
        }
        
        it_done <- 0L
        for (it in seq_len(max_iter)) {

          eg_y <- eval_f_g(y)
          g    <- eg_y$g
          g_rms_y <- sqrt(mean(g^2))
          if (!is.finite(g_rms_y)) break
          
          d      <- g * invHdiag
          f_y    <- eg_y$f
          gain_lin <- sum(g * d)
          

          s <- step
          eg_try <- NULL; f_try <- -Inf; delta_try <- NULL
          if (backtrack) {
            repeat {
              delta_try <- y + s * d
              eg_try    <- eval_f_g(delta_try)
              f_try     <- eg_try$f

              if (is.finite(f_try) && f_try >= f_y + ls_armijo * s * gain_lin) break
              s <- s * bt_factor
              if (s < 1e-16) {
                return(list(delta = delta, value = f_prev, iters = it - 1L, step = step, grad_rms = g_rms_y))
              }
            }
          } else {
            delta_try <- y + s * d
            eg_try    <- eval_f_g(delta_try)
            f_try     <- eg_try$f
            if (!is.finite(f_try)) {

              return(list(delta = delta, value = f_prev, iters = it - 1L, step = step, grad_rms = g_rms_y))
            }
          }
          

          if (f_try < f_prev) {

            t_k1 <- 1.0
            y    <- delta
            
            eg_y <- eval_f_g(y)
            g    <- eg_y$g
            g_rms_y <- sqrt(mean(g^2))
            if (!is.finite(g_rms_y)) break
            
            d        <- g * invHdiag
            f_y      <- eg_y$f
            gain_lin <- sum(g * d)
            
            s <- min(s * bt_factor, step)
            repeat {
              delta_try <- y + s * d
              eg_try    <- eval_f_g(delta_try)
              f_try     <- eg_try$f
              if (is.finite(f_try) && f_try >= f_y + ls_armijo * s * gain_lin) break
              s <- s * bt_factor
              if (s < 1e-16) {
                return(list(delta = delta, value = f_prev, iters = it - 1L, step = step, grad_rms = g_rms_y))
              }
            }
          } else {

            t_k1 <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
          }
          

          step <- s
          y_new <- delta_try + ((t_k - 1) / t_k1) * (delta_try - delta)
          
          g_rms_accept <- sqrt(mean(eg_try$g^2))
          delta  <- delta_try
          y      <- y_new
          t_k    <- t_k1
          f_prev <- f_try
          it_done <- it
          
          if (!is.finite(g_rms_accept) || g_rms_accept < tol_grad_rms) break
        }
        
        list(delta = delta, value = f_prev, iters = it_done, step = step,
             grad_rms = if (exists("g_rms_accept")) g_rms_accept else NA_real_)
      }
      
    }
  }
  
  

  update_sigma_rho_ratio_fix <- function(
    delta, n, S2, rho_old, sigma2_old,

    K_matvec = NULL, K_diag = NULL,
    damping  = 0.3,
    chol_threshold = 5000L, n_probe = 50, seed = 1L,
    pcg_tol = 1e-6, pcg_maxit = 1000,
    sigma2_min = 1e-6, sigma2_max = 1e+3
  ){
    stopifnot(inherits(S2, "sparseMatrix"))
    m <- length(delta)
    if (is.null(K_diag)) K_diag <- rep(0, m)
    K_diag <- as.numeric(K_diag)
    

    degL <- as.numeric(Matrix::rowSums(S2))
    cf <- NULL
    use_chol <- (m <= chol_threshold)
    if (use_chol) {
      M <- Diagonal(m) + rho_old * S2
      cf <- tryCatch(
        Cholesky(M, LDL = FALSE, Imult = 0, perm = TRUE),
        error = function(e) Cholesky(M + Diagonal(m, 1e-10), LDL = FALSE, Imult = 0, perm = TRUE)
      )
      A_apply      <- function(x) as.numeric(solve(cf, x))
      A_apply_many <- function(B) as.matrix(solve(cf, B))
    } else {
      M_mv  <- function(x) x + as.numeric(rho_old * (S2 %*% x))
      Mdiag <- pmax(1 + rho_old * degL, 1e-12)
      A_apply <- function(b) pcg_solve(M_mv, b, Mdiag = Mdiag, tol = pcg_tol, maxit = pcg_maxit)
      A_apply_many <- function(B) apply(B, 2L, A_apply)
    }
    

    q0 <- sum(delta * A_apply(delta))
    

    Oinv_vec <- function(x) A_apply(x) / sigma2_old
    H_mv     <- function(x) Oinv_vec(x) + K_diag * x
    

    Hdiag  <- (1 / sigma2_old) + K_diag
    Hsolve <- function(b) pcg_solve(H_mv, b, Mdiag = pmax(Hdiag, 1e-12),
                                    tol = min(1e-8, pcg_tol), maxit = max(3000, pcg_maxit))
    

    set.seed(seed)
    S <- max(10L, as.integer(n_probe))
    V <- matrix(sample(c(-1,1), m * S, replace = TRUE), m, S)
    
    Av      <- A_apply_many(V)
    HinvAv  <- apply(Av, 2L, Hsolve)
    t0      <- mean(colSums(V * HinvAv))
    
    print((q0+t0)/q0)
    
    
    sigma2_corr <- (q0+t0 ) / m
    sigma2_corr <- max(sigma2_min, min(sigma2_corr, sigma2_max))
    sigma2_new  <- (1 - damping) * sigma2_old + damping * sigma2_corr
    


    a   <- A_apply(delta)
    Sa  <- as.numeric(S2 %*% a)
    q1  <- sum(a * Sa)
    

    AS2V <- A_apply_many(as.matrix(S2 %*% V))
    t1   <- mean(colSums(V * AS2V))
    

    S2Av        <- as.matrix(S2 %*% Av)
    AS2Av       <- A_apply_many(S2Av)
    Hinv_AS2Av  <- apply(AS2Av, 2L, Hsolve)
    t2          <- mean(colSums(V * Hinv_AS2Av))
    

    S2_AS2V   <- as.matrix(S2 %*% AS2V)
    AS2AS2V   <- A_apply_many(S2_AS2V)
    t_SS      <- mean(colSums(V * AS2AS2V))
    

    AaSa       <- A_apply(Sa)
    d_tmp      <- as.numeric(S2 %*% AaSa)
    Ad_tmp     <- A_apply(d_tmp)
    u_SS_delta <- sum(delta * Ad_tmp)
    

    S2_AS2Av      <- as.matrix(S2 %*% AS2Av)
    AS2AS2Av      <- A_apply_many(S2_AS2Av)
    Hinv_AS2AS2Av <- apply(AS2AS2Av, 2L, Hsolve)
    t_SS_H        <- mean(colSums(V * Hinv_AS2AS2Av))
    

    g_rho <- 0.5 * ( (q1 + t2) / sigma2_old - t1 )
    h_rho <- 0.5 * ( t_SS - (2 / sigma2_old) * u_SS_delta - (2 / sigma2_old) * t_SS_H )
    

    h_safe    <- if (!is.finite(h_rho) || h_rho <= 0) max(abs(h_rho), 1e-8) else h_rho
    step_cap  <- 0.10 * (0.5 - 0.0)
    delta_rho <-  g_rho / h_safe
    delta_rho <- max(-step_cap, min(step_cap, delta_rho))
    
    rho_hat <- rho_old + delta_rho
    

    rho_lo <- 0.0
    rho_hi <- 0.5
    eps    <- 1e-8
    rho_hat <- max(rho_lo + eps, min(rho_hat, rho_hi - eps))
    rho_new <- (1 - damping) * rho_old + damping * rho_hat
    
    list(
      sigma2 = sigma2_new,
      rho    = rho_new,
      stats  = list(
        q0 = q0, t0 = t0, sigma2_corr = sigma2_corr,
        q1 = q1, t1 = t1, t2 = t2,
        t_SS = t_SS, u_SS_delta = u_SS_delta, t_SS_H = t_SS_H,
        g_rho = g_rho, h_rho = h_rho, delta_rho = delta_rho,
        rho_hat = rho_hat
      )
    )
  }
  
  
  


  fit_edge_frailty_ppl_block_pcg_tv_fast <- function(
    data_sim,
    zi1_pattern = "^Z1i_", zi2_pattern = "^Z2i_",
    max_outer = 10, tol = 5e-4, beta_maxit = 100,
    damping = 0.2, verbose = TRUE,
    sigma_fix0 = FALSE,
    rho_fix0 = FALSE,
    beta_ini = NULL,
    sigma2_ini =NULL,rho_ini=NULL,
    ties_method = c("breslow","efron"),
    chol_exact_threshold = 5000
  ) {
    stopifnot(is.list(data_sim), "df" %in% names(data_sim))
    df <- data_sim$df

    n  <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
    
    ops <- make_linegraph_ops(n, pairs = data_sim$pairs)
    scan <- build_tv_scanline_index(df,
                                    zedge_pat="^Zedge_",
                                    zi1_pat = zi1_pattern, zi2_pat = zi2_pattern)
    X_rows <- scan$X_rows
    p_edge <- length(grep("^Zedge_", colnames(X_rows)))
    gamma_names <- colnames(X_rows)[(p_edge+1):ncol(X_rows)]
    p <- ncol(X_rows); m <- max(scan$id_row)
    
    beta  <- rep(0, p)
    delta <- rep(0, m)
    if(is.null(sigma2_ini))
    {
      sigma2 <-1
    }
    else
    {
      sigma2 <- sigma2_ini
    }
    if(is.null(rho_ini))
    {
      rho <- 0.01
    }
    else
    {
      rho <- rho_ini
    }
    
    
    S2 <- make_S2(n, pairs = data_sim$pairs)
    history <- vector("list", max_outer)
    sig_dif <- 0.1
    for (it in 1:max_outer) {




      
      vdelta <- rep(0, m)
      res_beta <- beta_update_via_coxph_fit(
        X_rows = scan$X_rows,
        df = data_sim$df,
        scan =scan,
        sigma2 =sigma2,
        delta  = delta,
        vdelta = vdelta,
        ties   = ties_method ,
        strata = rep(1L, nrow(df))
      )
      beta_new <- res_beta$beta
      
      
      t_max <- max(scan$times)
      
      


       step_ddd <- t_max/10
 

      max_iter11 <- 5
      if(it==1)
      {
        max_iter11 <- 10
      }

      res_d <- delta_block_update_fista_tv_fast2(
        delta, beta_new, scan, ops, sigma2, rho,
        step_init = step_ddd,
        max_iter =  max_iter11, tol_grad_rms = 1e-4,
        chol_threshold = chol_exact_threshold, pcg_tol = 5e-7,
        ties_method = ties_method
      )
      delta_new <- res_d$delta

      







      
      if(sigma_fix0)
      {
        delta_new <- rep(0,length(res_d$delta))
      }
      







      
      
      









      

      
      {
        kern  <- make_scan_kernel(scan, beta_new,ties_method = ties_method)
        lam <- compute_cumLambda(scan, beta_new, delta_new, method = ties_method)
        


        Kdiag <- compute_K_diag_partial(
          scan, beta_new, delta_new,ties_method = ties_method)
        


        Hdiag <- (1 / sigma2) + pmax(Kdiag, 1e-12)
        



        



        
        up <- update_sigma_rho_ratio_fix(
          delta = delta_new, n = n, S2 = S2,
          rho_old = rho, sigma2_old = sigma2,
          K_matvec = function(x) kern$K_matvec(x, delta_new),
          K_diag   = Kdiag
        )
        
        sig_dif <- abs(up$sigma2-sigma2)
        
        sigma2_new <- up$sigma2
        rho_new    <- up$rho
        if(rho_fix0)
        {
          rho_new <- 0
        }
      }
      
      {















        
      }
      
      
      d_beta  <- sqrt(sum((beta_new  - beta )^2))
      d_delta <- sqrt(sum((delta_new - delta)^2))
      d_hyper <- sqrt(sum((c(sigma2_new, rho_new) - c(sigma2, rho))^2))
      
      beta  <- beta_new
      delta <- delta_new
      sigma2 <- sigma2_new
      rho    <- rho_new
      
      beta_edge <- beta[seq_len(p_edge)]
      gamma_vec <- setNames(beta[(p_edge+1):p], gamma_names)
      
      history[[it]] <- list(
        iter = it,
        beta_edge = beta_edge,
        gamma = gamma_vec,
        sigma = sqrt(sigma2),
        rho   = rho,
        d_beta = d_beta, d_delta = d_delta, d_hyper = d_hyper
      )
      
      if (isTRUE(verbose)) {
        cat(sprintf("[TV-fast %02d] sigma=%.4g rho=%.4f | dβ=%.2e dδ=%.2e dϑ=%.2e\n",
                    it, sqrt(sigma2), rho, d_beta, d_delta, d_hyper))
        cat("  beta_edge =", paste(round(beta_edge, 4), collapse=", "),
            "| gamma =", paste(paste0(names(gamma_vec), "=", round(gamma_vec, 4)),
                               collapse=", "), "\n")
      }

      if (max( d_hyper) < tol) break
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



########## Regression parameter inference
{
  
  build_theta_delta_ops <- function(scan, beta, delta, ties_method = c("breslow","efron")) {
    ties_method <- match.arg(ties_method)
    X  <- scan$X_rows
    id <- as.integer(scan$id_row)
    p  <- ncol(X); m <- length(delta)
    
    a_row <- exp(as.vector(X %*% beta))
    d_id  <- exp(delta)
    
    ss <- scan$start_sorted; es <- scan$stop_sorted
    rs <- scan$rows_start;    re <- scan$rows_stop
    times <- scan$times; nR <- nrow(X)
    
    mt_by_time <- vapply(scan$events_rowidx_by_time, length, 0L)
    grp        <- scan$events_rowidx_by_time
    

    Htt_matrix <- function() {
      Wid <- numeric(m)
      Wxe <- matrix(0.0, m, p)
      S0  <- 0.0
      Wxx <- matrix(0.0, p, p)
      si <- 1L; ei <- 1L
      Htt <- matrix(0.0, p, p)
      
      for (k in seq_along(times)) {
        t <- times[k]

        while (si <= nR && ss[si] < t) {
          r <- rs[si]; eid <- id[r]
          wr <- a_row[r] * d_id[eid]
          xr <- X[r, , drop = FALSE]
          S0 <- S0 + wr
          Wid[eid] <- Wid[eid] + wr
          Wxe[eid, ] <- Wxe[eid, ] + wr * xr
          Wxx <- Wxx + wr * crossprod(xr)
          si <- si + 1L
        }

        while (ei <= nR && es[ei] < t) {
          r <- re[ei]; eid <- id[r]
          wr <- a_row[r] * d_id[eid]
          xr <- X[r, , drop = FALSE]
          S0 <- S0 - wr
          Wid[eid] <- Wid[eid] - wr
          Wxe[eid, ] <- Wxe[eid, ] - wr * xr
          Wxx <- Wxx - wr * crossprod(xr)
          ei <- ei + 1L
        }
        
        mk <- mt_by_time[k]
        if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
        
        if (ties_method == "breslow") {
          barx <- colSums(Wxe) / S0
          Htt  <- Htt + mk * ( Wxx / S0 - tcrossprod(barx, barx) )
        } else {

          ev_rows <- grp[[k]]
          if (length(ev_rows)) {
            eids <- id[ev_rows]
            wrev <- a_row[ev_rows] * d_id[eids]
            Xev  <- X[ev_rows, , drop = FALSE]
            Ek   <- sum(wrev)
            WXevent_total  <- colSums(Xev * wrev)
            WXXevent_total <- crossprod(Xev, wrev * Xev)
          } else {
            Ek <- 0
            WXevent_total  <- rep(0.0, p)
            WXXevent_total <- matrix(0.0, p, p)
          }
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (!(denom > 0) || !is.finite(denom)) next
            Wxx_ell <- Wxx - alpha * WXXevent_total
            barx    <- (colSums(Wxe) - alpha * WXevent_total) / denom
            Htt     <- Htt + ( Wxx_ell / denom - tcrossprod(barx, barx) )
          }
        }
      }
      Htt
    }
    

    Htd_matvec <- function(u_delta) {
      Wid <- numeric(m); Wxe <- matrix(0.0, m, p); S0 <- 0.0
      si <- 1L; ei <- 1L
      out <- numeric(p)
      
      for (k in seq_along(times)) {
        t <- times[k]
        while (si <= nR && ss[si] < t) {
          r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 + wr
          Wid[eid] <- Wid[eid] + wr
          Wxe[eid, ] <- Wxe[eid, ] + wr * X[r, ]
          si <- si + 1L
        }
        while (ei <= nR && es[ei] < t) {
          r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 - wr
          Wid[eid] <- Wid[eid] - wr
          Wxe[eid, ] <- Wxe[eid, ] - wr * X[r, ]
          ei <- ei + 1L
        }
        mk <- mt_by_time[k]
        if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
        
        if (ties_method == "breslow") {
          barx <- colSums(Wxe) / S0

          s_vec <- colSums(Wxe * as.numeric(u_delta)) / S0
          s_scl <- sum(Wid * as.numeric(u_delta)) / S0
          out <- out + mk * (s_vec - s_scl * barx)
        } else {
          ev_rows <- grp[[k]]
          wevent  <- numeric(m); WXevent <- matrix(0.0, m, p)
          if (length(ev_rows)) {
            eids <- id[ev_rows]; wrev <- a_row[ev_rows] * d_id[eids]
            acc1 <- rowsum(wrev, group = eids, reorder = FALSE)
            wevent[as.integer(rownames(acc1))] <- as.numeric(acc1)
            WX_tmp <- wrev * X[ev_rows, , drop = FALSE]
            sp <- split(seq_along(eids), f = eids)
            for (kk in names(sp)) {
              idx <- sp[[kk]]
              WXevent[as.integer(kk), ] <- colSums(WX_tmp[idx, , drop = FALSE])
            }
          }
          Ek <- sum(wevent)
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (!(denom > 0) || !is.finite(denom)) next
            Wid_ell <- Wid - alpha * wevent
            Wxe_ell <- Wxe - alpha * WXevent
            barx <- colSums(Wxe_ell) / denom
            s_vec <- colSums(Wxe_ell * as.numeric(u_delta)) / denom
            s_scl <- sum(Wid_ell * as.numeric(u_delta)) / denom
            out <- out + (s_vec - s_scl * barx)
          }
        }
      }
      out
    }
    

    Hdt_matvec <- function(v_theta) {
      Wid <- numeric(m); Wxe <- matrix(0.0, m, p); S0 <- 0.0
      si <- 1L; ei <- 1L
      out <- numeric(m)
      
      for (k in seq_along(times)) {
        t <- times[k]
        while (si <= nR && ss[si] < t) {
          r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 + wr
          Wid[eid] <- Wid[eid] + wr
          Wxe[eid, ] <- Wxe[eid, ] + wr * X[r, ]
          si <- si + 1L
        }
        while (ei <= nR && es[ei] < t) {
          r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 - wr
          Wid[eid] <- Wid[eid] - wr
          Wxe[eid, ] <- Wxe[eid, ] - wr * X[r, ]
          ei <- ei + 1L
        }
        mk <- mt_by_time[k]
        if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
        
        if (ties_method == "breslow") {
          barx <- colSums(Wxe) / S0
          s1   <- as.numeric((Wxe %*% v_theta) / S0)
          s2   <- as.numeric((barx %*% v_theta))
          out  <- out + mk * (s1 - s2 * (Wid / S0))
        } else {
          ev_rows <- grp[[k]]
          wevent  <- numeric(m); WXevent <- matrix(0.0, m, p)
          if (length(ev_rows)) {
            eids <- id[ev_rows]; wrev <- a_row[ev_rows] * d_id[eids]
            acc1 <- rowsum(wrev, group = eids, reorder = FALSE)
            wevent[as.integer(rownames(acc1))] <- as.numeric(acc1)
            WX_tmp <- wrev * X[ev_rows, , drop = FALSE]
            sp <- split(seq_along(eids), f = eids)
            for (kk in names(sp)) {
              idx <- sp[[kk]]
              WXevent[as.integer(kk), ] <- colSums(WX_tmp[idx, , drop = FALSE])
            }
          }
          Ek <- sum(wevent)
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (!(denom > 0) || !is.finite(denom)) next
            Wid_ell <- Wid - alpha * wevent
            Wxe_ell <- Wxe - alpha * WXevent
            barx <- colSums(Wxe_ell) / denom
            s1   <- as.numeric((Wxe_ell %*% v_theta) / denom)
            s2   <- as.numeric((barx %*% v_theta))
            out  <- out + (s1 - s2 * (Wid_ell / denom))
          }
        }
      }
      out
    }
    
    list(Htt_matrix = Htt_matrix,
         Htd_matvec = Htd_matvec,
         Hdt_matvec = Hdt_matvec)
  }
  



  edge_joint_scores <- function(scan, beta, delta, ties_method = c("breslow","efron"), S2 = NULL, sigma2 = NULL, rho = 0) {
    ties_method <- match.arg(ties_method)
    X  <- scan$X_rows
    id <- as.integer(scan$id_row)
    p  <- ncol(X); m <- length(delta)
    
    a_row <- exp(as.vector(X %*% beta))
    d_id  <- exp(delta)
    
    ss <- scan$start_sorted; es <- scan$stop_sorted
    rs <- scan$rows_start;    re <- scan$rows_stop
    times <- scan$times; nR <- nrow(X)
    
    mt_by_time <- vapply(scan$events_rowidx_by_time, length, 0L)
    grp        <- scan$events_rowidx_by_time
    

    Utheta_by_edge <- matrix(0.0, m, p)
    Udelta_by_edge <- numeric(m)
    
    Wid <- numeric(m)
    Wxe <- matrix(0.0, m, p)
    S0  <- 0.0
    si <- 1L; ei <- 1L
    
    for (k in seq_along(times)) {
      t <- times[k]
      while (si <= nR && ss[si] < t) {
        r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
        S0 <- S0 + wr
        Wid[eid] <- Wid[eid] + wr
        Wxe[eid, ] <- Wxe[eid, ] + wr * X[r, ]
        si <- si + 1L
      }
      while (ei <= nR && es[ei] < t) {
        r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
        S0 <- S0 - wr
        Wid[eid] <- Wid[eid] - wr
        Wxe[eid, ] <- Wxe[eid, ] - wr * X[r, ]
        ei <- ei + 1L
      }
      mk <- mt_by_time[k]
      if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
      

      ev_rows <- grp[[k]]
      if (length(ev_rows)) {
        eids <- id[ev_rows]



        sp <- split(seq_along(eids), f = eids)
        for (kk in names(sp)) {
          idx <- sp[[kk]]
          Utheta_by_edge[as.integer(kk), ] <- Utheta_by_edge[as.integer(kk), ] +
            colSums(X[ev_rows[idx], , drop = FALSE])
        }

        tab <- table(eids)
        Udelta_by_edge[as.integer(names(tab))] <- Udelta_by_edge[as.integer(names(tab))] + as.numeric(tab)
      }
      
      if (ties_method == "breslow") {

        Utheta_by_edge <- Utheta_by_edge - (mk / S0) * Wxe
        Udelta_by_edge <- Udelta_by_edge - (mk / S0) * Wid
      } else {

        wevent  <- numeric(m); WXevent <- matrix(0.0, m, p)
        if (length(ev_rows)) {
          eids <- id[ev_rows]; wrev <- a_row[ev_rows] * d_id[eids]
          acc1 <- rowsum(wrev, group = eids, reorder = FALSE)
          wevent[as.integer(rownames(acc1))] <- as.numeric(acc1)
          WX_tmp <- wrev * X[ev_rows, , drop = FALSE]
          sp2 <- split(seq_along(eids), f = eids)
          for (kk in names(sp2)) {
            idx <- sp2[[kk]]
            WXevent[as.integer(kk), ] <- colSums(WX_tmp[idx, , drop = FALSE])
          }
        }
        Ek <- sum(wevent)
        for (ell in 0:(mk-1)) {
          alpha <- ell / mk
          denom <- S0 - alpha * Ek
          if (!(denom > 0) || !is.finite(denom)) next
          Utheta_by_edge <- Utheta_by_edge - ( (Wxe - alpha * WXevent) / denom )
          Udelta_by_edge <- Udelta_by_edge - ( (Wid - alpha * wevent)  / denom )
        }
      }
    }
    
    list(U_theta_by_edge = Utheta_by_edge,
         U_delta_by_edge = Udelta_by_edge)
  }
  
  

  {
    
    edge_joint_scores <- function(
    scan, beta, delta,
    ties_method = c("breslow","efron"),

    Q_matvec = NULL,

    S2 = NULL, sigma2 = NULL, rho = 0
    ){
      ties_method <- match.arg(ties_method)
      X  <- scan$X_rows
      id <- as.integer(scan$id_row)
      p  <- ncol(X); m <- length(delta)
      
      a_row <- exp(as.vector(X %*% beta))
      d_id  <- exp(delta)
      
      ss <- scan$start_sorted; es <- scan$stop_sorted
      rs <- scan$rows_start;    re <- scan$rows_stop
      times <- scan$times; nR <- nrow(X)
      
      mt_by_time <- vapply(scan$events_rowidx_by_time, length, 0L)
      grp        <- scan$events_rowidx_by_time
      
      Utheta_by_edge     <- matrix(0.0, m, p)
      Udelta_raw_by_edge <- numeric(m)
      
      Wid <- numeric(m)
      Wxe <- matrix(0.0, m, p)
      S0  <- 0.0
      si <- 1L; ei <- 1L
      
      for (k in seq_along(times)) {
        t <- times[k]

        while (si <= nR && ss[si] < t) {
          r <- rs[si]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 + wr
          Wid[eid] <- Wid[eid] + wr
          Wxe[eid, ] <- Wxe[eid, ] + wr * X[r, ]
          si <- si + 1L
        }
        while (ei <= nR && es[ei] < t) {
          r <- re[ei]; eid <- id[r]; wr <- a_row[r]*d_id[eid]
          S0 <- S0 - wr
          Wid[eid] <- Wid[eid] - wr
          Wxe[eid, ] <- Wxe[eid, ] - wr * X[r, ]
          ei <- ei + 1L
        }
        
        mk <- mt_by_time[k]
        if (!(S0 > 0) || !is.finite(S0) || mk == 0L) next
        

        ev_rows <- grp[[k]]
        if (length(ev_rows)) {
          eids <- id[ev_rows]
          sp <- split(seq_along(eids), f = eids)
          for (kk in names(sp)) {
            idx <- sp[[kk]]
            Utheta_by_edge[as.integer(kk), ] <-
              Utheta_by_edge[as.integer(kk), ] +
              colSums(X[ev_rows[idx], , drop = FALSE])
          }
          tab <- table(eids)
          Udelta_raw_by_edge[as.integer(names(tab))] <-
            Udelta_raw_by_edge[as.integer(names(tab))] + as.numeric(tab)
        }
        

        if (ties_method == "breslow") {
          Utheta_by_edge     <- Utheta_by_edge - (mk / S0) * Wxe
          Udelta_raw_by_edge <- Udelta_raw_by_edge - (mk / S0) * Wid
        } else {
          wevent  <- numeric(m); WXevent <- matrix(0.0, m, p)
          if (length(ev_rows)) {
            eids <- id[ev_rows]; wrev <- a_row[ev_rows] * d_id[eids]
            acc1 <- rowsum(wrev, group = eids, reorder = FALSE)
            wevent[as.integer(rownames(acc1))] <- as.numeric(acc1)
            WX_tmp <- wrev * X[ev_rows, , drop = FALSE]
            sp2 <- split(seq_along(eids), f = eids)
            for (kk in names(sp2)) {
              idx <- sp2[[kk]]
              WXevent[as.integer(kk), ] <- colSums(WX_tmp[idx, , drop = FALSE])
            }
          }
          Ek <- sum(wevent)
          for (ell in 0:(mk-1)) {
            alpha <- ell / mk
            denom <- S0 - alpha * Ek
            if (!(denom > 0) || !is.finite(denom)) next
            Utheta_by_edge     <- Utheta_by_edge - ( (Wxe  - alpha * WXevent) / denom )
            Udelta_raw_by_edge <- Udelta_raw_by_edge - ( (Wid  - alpha * wevent)  / denom )
          }
        }
      }
      

      if (!is.null(Q_matvec)) {
        Qdelta <- as.numeric(Q_matvec(delta))
        
      } else if (!is.null(S2) && !is.null(sigma2)) {
        if (!inherits(S2, "sparseMatrix")) {
          S2 <- Matrix::Matrix(S2, sparse = TRUE)
        }

        if (!inherits(S2, "symmetricMatrix")) {
          S2 <- Matrix::forceSymmetric(S2, uplo = "U")
        }
        if (!inherits(S2, "dsCMatrix")) {
          S2 <- methods::as(S2, "dsCMatrix")
        }
        if (nrow(S2) != m) stop("nrow(S2) must equal length(delta)")
        if (!(is.numeric(sigma2) && length(sigma2) == 1L && sigma2 > 0)) {
          stop("sigma2 must be a positive scalar.")
        }
        

        A <- Matrix::Diagonal(m) + rho * S2
        


















      b <-  as.matrix(delta, ncol = 1)
      y  <- Matrix::solve(A, b)
     
      Qdelta <- as.numeric(y) / sigma2
        
        
      Udelta_pen_by_edge <- Udelta_raw_by_edge - Qdelta
      
      list(
        U_theta_by_edge      = Utheta_by_edge,
        U_delta_raw_by_edge  = Udelta_raw_by_edge,
        U_delta_by_edge  = Udelta_pen_by_edge,
        Qdelta               = Qdelta
      )
    }
    
    
  }
  
  
  
  

  {







    cox_frailty_theta_vcov <- function(
    scan, pairs,
    beta_hat,
    delta_hat, sigma2, rho,
    ties_method = c("breslow","efron"),
    sand = FALSE,
    cluster = c("dyad","gray","edge","node"),
    beta_col_pattern  = "^Zedge_",
    gamma_col_pattern = "^Xi_plus_Xj_",
    ridge = 1e-8, chol_threshold = 5000L,
    pcg_tol = 1e-6, pcg_maxit = 2000,beta_adjust=TRUE,S2 =NULL
    ){
      ties_method <- match.arg(ties_method)
      cluster     <- match.arg(cluster)

      apply_BTB <- function(X) {
        X <- if (is.matrix(X)) X else cbind(X)
        if (!is.null(ops$B)) {
          as.matrix(crossprod(ops$B, ops$B %*% X))
        } else if (!is.null(ops$S2_matmat) && is.function(ops$S2_matmat)) {
          ops$S2_matmat(X) + 2 * X
        } else {
          stop("Neither ops$B nor ops$S2_matmat function is available.")
        }
      }
      
      X  <- scan$X_rows
      p  <- ncol(X)
      m  <- length(delta_hat)
      theta_hat <- beta_hat
      

      beta_idx  <- grep(beta_col_pattern,  colnames(X), value = FALSE)
      gamma_idx <- grep(gamma_col_pattern, colnames(X), value = FALSE)







      n_node <- max(pairs)
      ops    <- make_linegraph_ops(n_node, pairs)
      Msolve <- make_M_solver(ops, rho, chol_threshold, pcg_tol, pcg_maxit)
      OmegaInv_apply <- function(x) (1/sigma2) * Msolve(x)
      

      sk <- make_scan_kernel(scan, beta = theta_hat, ties_method = ties_method)
      K_apply <- function(x) sk$K_matvec(x, delta = delta_hat)
      

      Hdd_matvec <- function(x) K_apply(x) + OmegaInv_apply(x)
      
      Mdiag_pcg  <- approx_diag_K(scan, theta_hat, delta_hat, ties_method) + (1/sigma2)
      Hdd_solve  <- function(rhs) pcg_solve(Hdd_matvec, rhs, Mdiag = Mdiag_pcg,
                                            tol = pcg_tol, maxit = pcg_maxit)
      
      
      Hdd_solve_mat <- function(R){ R <- if (is.matrix(R)) R else cbind(R)
      out <- matrix(0.0, nrow(R), ncol(R))
      for (j in seq_len(ncol(R))) out[, j] <- Hdd_solve(R[, j])
      out
      }
      

      Htt_fallback <- function() {
        H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
        H_ops$Htt_matrix()
      }
      Htd_matmul_fallback <- function(U){
        H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
        U <- if (is.matrix(U)) U else cbind(U)
        do.call(cbind, lapply(seq_len(ncol(U)), function(j) H_ops$Htd_matvec(U[, j])))
      }
      Hdt_matmul_fallback <- function(V){
        H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
        V <- if (is.matrix(V)) V else cbind(V)
        do.call(cbind, lapply(seq_len(ncol(V)), function(j) H_ops$Hdt_matvec(V[, j])))
      }
      if (exists("build_theta_delta_ops_matmul") && is.function(build_theta_delta_ops_matmul)) {
        TDops <- build_theta_delta_ops_matmul(scan, theta_hat, delta_hat, ties_method)
        I_tt        <- TDops$Htt_matrix() + diag(ridge, p)
        Htd_matmul  <- TDops$Htd_matmul
        Hdt_matmul  <- TDops$Hdt_matmul
      } else {
        I_tt        <- Htt_fallback() + diag(ridge, p)
        Htd_matmul  <- Htd_matmul_fallback
        Hdt_matmul  <- Hdt_matmul_fallback
      }
      

      Ip <- diag(1, p)
      Q  <- Hdt_matmul(Ip)
      Z  <- Hdd_solve_mat(Q)
      S  <- I_tt - Htd_matmul(Z)
      cf_S <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(cf_S)) { S <- S + diag(ridge, p); cf_S <- chol(S) }
      S_inv <- chol2inv(cf_S)



      if (!sand) {
        Var_theta <- S_inv
      } else {

        if (cluster == "gray") {


          Mid <- matrix(0.0, p, p)
          for (j in 1:p) {
            ej <- rep(0.0, p); ej[j] <- 1
            col1 <- I_tt[, j]
            z1   <- Hdd_solve(Hdt_matmul(ej))
            t12  <- Htd_matmul(z1)
            z2   <- Hdd_solve(K_apply(z1))
            t3   <- Htd_matmul(z2)
            Mid[, j] <- col1 - 2 * t12 + t3
          }
        } else {
          

          es <- edge_joint_scores(scan, beta_hat, delta_hat, ties_method ,S2 =S2,sigma2 = sigma2,rho=rho)
          Utheta_e <- as.matrix(es$U_theta_by_edge)
          udelta   <- as.numeric(es$U_delta_by_edge)
    



          
          if (cluster == "edge") {
            M_tt <- crossprod(Utheta_e)
            M_dt <- Utheta_e * udelta
            Mdd_apply <- function(X) (udelta^2) * (if (is.matrix(X)) X else cbind(X))
          } else if (cluster == "node") {
            BTB_Uθ <- apply_BTB(Utheta_e)
            M_tt   <- crossprod(Utheta_e, BTB_Uθ)
            M_dt   <- BTB_Uθ * udelta
            Mdd_apply <- function(X) {
              X  <- if (is.matrix(X)) X else cbind(X)
              DX <- X * udelta
              udelta * apply_BTB(DX)
            }
          } else {
            c_v <- if (n_node > 1) n_node/(n_node - 1) else 1
            c_e <- if (m      > 1) m     /(m      - 1) else 1
            
            BTB_Uθ     <- apply_BTB(Utheta_e)
            M_nodes_tt <- crossprod(Utheta_e, BTB_Uθ)
            M_edges_tt <- crossprod(Utheta_e)
            M_tt <- 2 * (c_v * M_nodes_tt) - (c_e * M_edges_tt)
            
            M_nodes_dt <- BTB_Uθ * udelta
            M_edges_dt <- Utheta_e * udelta
            M_dt <- 2 * (c_v * M_nodes_dt) - (c_e * M_edges_dt)
            
            Mdd_apply <- function(X) {
              X  <- if (is.matrix(X)) X else cbind(X)
              DX <- X * udelta
              BTB_DX <- apply_BTB(DX)
              2 * (udelta * BTB_DX) - (udelta^2) * X
            }
          }
          
          

          Mid <- M_tt
          Y   <- Hdd_solve_mat(M_dt)
          Mid <- Mid - Htd_matmul(Y)
          Mid <- Mid - crossprod(M_dt, Z)
          T2  <- Hdd_solve_mat(Mdd_apply(Z))
          Mid <- Mid + Htd_matmul(T2)
        }
        
        Var_theta <- S_inv %*% Mid %*% S_inv 
      }
      Var_theta_unad <- S_inv

      if(beta_adjust)
      {
        M_tt <- crossprod(Utheta_e)
        M_dt <- Utheta_e * udelta
        Mdd_apply <- function(X) (udelta^2) * (if (is.matrix(X)) X else cbind(X))

        Mid <- M_tt
        Y   <- Hdd_solve_mat(M_dt)
        Mid <- Mid - Htd_matmul(Y)
        Mid <- Mid - crossprod(M_dt, Z)
        T2  <- Hdd_solve_mat(Mdd_apply(Z))
        Mid <- Mid + Htd_matmul(T2)
        Var_theta_unad <- S_inv %*% Mid %*% S_inv
      }
      
      
      take <- function(M, r, c) M[r, c, drop = FALSE]
      Var_beta  <- take(Var_theta, beta_idx,  beta_idx)
      Var_gamma <- take(Var_theta, gamma_idx, gamma_idx)
      Cov_bg    <- take(Var_theta, beta_idx,  gamma_idx)
      
      num_n <- max(pairs)
      se <- sqrt(pmax(diag(Var_theta), 0)) 
        
      

      se[beta_idx] <- sqrt(pmax(diag(Var_theta_unad)[beta_idx], 0))
      names(se) <- colnames(X)
      ci <- cbind(lower = theta_hat - 1.96 * se,
                  upper = theta_hat + 1.96 * se)
      rownames(ci) <- colnames(X)
      
      list(
        Var_theta = Var_theta,
        Var_beta  = Var_beta,
        Var_gamma = Var_gamma,
        Cov_beta_gamma = Cov_bg,
        se = se,
        ci = ci,
        S_inv= S_inv,
        index = list(beta = beta_idx, gamma = gamma_idx),
        names = list(theta = colnames(X),
                     beta  = colnames(X)[beta_idx],
                     gamma = colnames(X)[gamma_idx]),
        info = list(
          sand = sand, cluster = if (sand) cluster else "none",
          ties_method = ties_method,
          notes = if (!sand) "Model-based: [H^{-1}]_{θθ}" else
            switch(cluster,
                   gray = "Sandwich: [H^{-1}]_{θθ} * Gray I * [H^{-1}]_{θθ}",
                   edge = "Sandwich: edge clusters (Σ_e U_eU_e^T)",
                   node = "Sandwich: node clusters (Σ_v U_vU_v^T)",
                   dyad = "Sandwich: dyadic two-way (2Σ_v - Σ_e)"))
      )
    }
    
  }
}
  






























































































































































































































































































































































































































































































































































































































  
 

}



## Variance inference
{

  compute_edge_M_and_K_strict <- function(scan, beta,
                                          ties_method = c("breslow","efron")) {
    ties_method <- match.arg(ties_method)
    
    X  <- scan$X_rows
    id <- as.integer(scan$id_row)
    m  <- max(id)
    a_row <- exp(as.vector(X %*% beta))
    
    times <- scan$times
    grp   <- scan$events_rowidx_by_time
    
    
    mt_by_time <- vapply(grp, length, 0L)
    
  
    
    ss <- scan$start_sorted; es <- scan$stop_sorted
    rs <- scan$rows_start;    re <- scan$rows_stop
    nR <- nrow(X)
    
    S_vec <- numeric(m)
    Kdiag <- numeric(m)
    Wid   <- numeric(m)
    S0    <- 0.0
    si <- 1L; ei <- 1L
    
    for (k in seq_along(times)) {
      t <- times[k]
      

      while (si <= nR && ss[si] < t) {
        r <- rs[si]; eid <- id[r]; wr <- a_row[r]
        S0 <- S0 + wr
        Wid[eid] <- Wid[eid] + wr
        si <- si + 1L
      }
      

      while (ei <= nR && es[ei] < t) {
        r <- re[ei]; eid <- id[r]; wr <- a_row[r]
        S0 <- S0 - wr
        Wid[eid] <- Wid[eid] - wr
        ei <- ei + 1L
      }
      
      mk <- mt_by_time[k]
      if (!(S0 > 0) || !is.finite(S0)) next
      
      ev_rows <- grp[[k]]
      if (mk > 0L && length(ev_rows)) {
        ev_eids <- id[ev_rows]
        S_vec[ev_eids] <- S_vec[ev_eids] + 1
      }
      
      if (mk > 0L) {
        if (ties_method == "breslow") {
          p <- Wid / S0
          S_vec <- S_vec - mk * p
          Kdiag <- Kdiag + mk *  ( (p ))
        } else {

          Evec <- numeric(m)
          if (length(ev_rows)) {
            for (r in ev_rows) Evec[id[r]] <- Evec[id[r]] + a_row[r]
          }
          E_tot <- sum(Evec)
          
          for (l in 0:(mk - 1L)) {
            frac <- l / mk
            S0_l  <- S0  - frac * E_tot
            if (!(S0_l > 0) || !is.finite(S0_l)) next
            Wid_l <- Wid - frac * Evec
            p_l   <- Wid_l / S0_l
            S_vec <- S_vec - p_l
            Kdiag <- Kdiag + ((p_l))
          }
        }
      }
    }
    
    Kdiag <- pmax(Kdiag, 0)
    list(M = S_vec, K_diag = Kdiag)
  }
  
  
  
  {

 

      {

        sigma2_Z_from_SK <- function(S_vec, K_diag,
                                     sandwich = TRUE,
                                     eps = 1e-16) {
          stopifnot(length(S_vec) == length(K_diag))
          Kd <- pmax(as.numeric(K_diag), eps)
          invK   <- 1 / Kd
          invK2  <- invK * invK
          Q   <- sum((S_vec * invK)^2)
          t1  <- sum(invK)
          t2  <- sum(invK2)
          

          if(sandwich)
          {
            B <-  sum(S_vec*invK2 *S_vec )
            v1 <- sqrt(2 * t2)
            v2 <- sqrt(B)
            
          }
          Z   <- (Q - t1) / v2
          p1  <- pnorm(Z, lower.tail = FALSE)
          p2  <- 2 * pnorm(abs(Z), lower.tail = FALSE)
          
          
          list(
            Z = Z, p_one_sided = p1, p_two_sided = p2,
            Q = Q, U = 0.5 * (Q - t1),
            tr_Kinv = t1, tr_Kinv2 = t2,
            method = "Z-test for H0: sigma^2 = 0 (unprojected)",
            alternative = "sigma^2 > 0"
          )
        }
        

        sigma2_Z_from_scan <- function(scan, beta, ties_method = c("breslow","efron"),
                                       eps = 1e-16) {
          ties_method <- match.arg(ties_method)
          sk <- compute_edge_M_and_K_strict(scan, beta, ties_method = ties_method)
          out <- sigma2_Z_from_SK(sk$M, sk$K_diag, eps = eps)
          out$ties_method <- ties_method
          out
        }
        
      }

     
  }

  
  score_test_sigma2_boundary <- function(
    scan,
    beta_hat_null,
    delta=NULL ,
    full = FALSE,
    ties_method = c("breslow","efron"),
    eps = 1e-12,
    prefer = c("davies","imhof")
  ){
    ties_method <- match.arg(ties_method)
    prefer <- match.arg(prefer)
    
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
      stop("Package 'CompQuadForm' is required for exact p-values (Davies/Imhof).")
    }
    

    beta_tilde <- beta_hat_null
    

    MK   <- compute_edge_M_and_K_strict(scan, beta_tilde ,ties_method="breslow")
    Svec <- as.numeric(MK$M)
    
    Kd   <- as.numeric(MK$K_diag)
    
    if(full)
    {

      estF0 <- estimate_baseline_cumhaz(scan, beta, ties_method = "efron")
      

      out <- compute_edge_M_and_K_full(
        scan, beta,
        F0_fun = estF0$F0_fun,
        delta  = NULL,
        return_mu_by_row = TRUE
      )
      
      Svec <- as.numeric(out$M)
      Kd   <- as.numeric(out$K_diag)
    }




    ok <- which(is.finite(Svec) & is.finite(Kd) & (Kd > eps))
    if (length(ok) == 0L) {
      return(list(
        statistic = list(Q0 = NA_real_),
        eigen = list(lambda = numeric(0L)),
        pval = list(davies = NA_real_, imhof = NA_real_),
        beta_tilde = beta_tilde,
        method = paste0("Exact σ^2 boundary score (", ties_method, ")"),
        note = "No valid components (non-finite or tiny K_diag)."
      ))
    }
    Svec <- Svec[ok]
    Kd   <- Kd[ok]
    

    Q0 <- sum( (Svec / Kd)^2 )
    

    lambda <- 1 / Kd
    

    p_davies <- NA_real_
    p_imhof  <- NA_real_
    

    dav <- try(
      CompQuadForm::davies(q = Q0, lambda = lambda, acc = 1e-12, lim = 1e6),
      silent = TRUE
    )
    if (!inherits(dav, "try-error") && is.list(dav) && is.numeric(dav$Qq)) {
      
      p_davies <- min(max(dav$Qq, 0), 1)
    }
    

    imh <- try(
      CompQuadForm::imhof(q = Q0, lambda = lambda),
      silent = TRUE
    )
    if (!inherits(imh, "try-error") && is.list(imh) && is.numeric(imh$Qq)) {
      p_imhof <- min(max(imh$Qq, 0), 1)
    }
    
    

    main_p <- if (prefer == "davies") p_davies else p_imhof
    
    list(
      statistic = list(Q0 = Q0),
      eigen = list(lambda = lambda, K_diag = Kd),
      pval = list(main = main_p, davies = p_davies, imhof = p_imhof),
      beta_tilde = beta_tilde,
      method = paste0("Exact σ^2 boundary score via ", toupper(prefer),
                      " (", ties_method, "); eigen(M)=1/K_diag")
    )
  }
  

  




  {
    wald_ci_sigma2_rho <- function(
    delta, n, S2,pairs = NULL,
    sigma2_hat, rho_hat,
    K_diag = NULL,
    conf_level = 0.95,
    chol_threshold = 5000L,
    sandwich =FALSE,
    n_mc = 128L, n_subspace = 24L,
    seed = 1L,
    exact_trace = TRUE,
    pcg_tol = 1e-6, pcg_maxit = 2000,
    use_nonfrozen_hessian = TRUE,
    nonfrozen_cap = 10,
    info_floor_rel = 1e-12,
    try_chol_first = TRUE,
    rho_bounds = c(-1/(2*(n-2)), 0.5),
    return_internals = TRUE
    ){
      stopifnot(inherits(S2, "sparseMatrix"))
      m <- length(delta)
      if (is.null(K_diag)) K_diag <- rep(0, m)
      K_diag <- as.numeric(K_diag)
      rho_lo <- rho_bounds[1]; rho_hi <- rho_bounds[2]
      
      if (!exists("pcg_solve")) {
        pcg_solve <- function(matvec, b, Mdiag = NULL, tol = 1e-6, maxit = 2000) {
          n <- length(b); x <- numeric(n)
          r <- b - matvec(x)
          z <- if (is.null(Mdiag)) r else r / Mdiag
          p <- z; rz_old <- sum(r * z)
          if (sqrt(rz_old) < tol) return(x)
          for (k in 1:maxit) {
            Ap <- matvec(p)
            denom <- sum(p * Ap); if (abs(denom) < 1e-30) break
            alpha <- rz_old / denom
            x <- x + alpha * p
            r <- r - alpha * Ap
            if (sqrt(sum(r*r)) < tol) break
            z <- if (is.null(Mdiag)) r else r / Mdiag
            rz_new <- sum(r * z)
            beta <- rz_new / rz_old
            p <- z + beta * p
            rz_old <- rz_new
          }
          x
        }
      }
      

      trace_hutchpp <- function(apply_op_many, n, s = 32L, k = 8L, seed = NULL) {
        if (!is.null(seed)) set.seed(seed)
        s <- max(4L, as.integer(s)); k <- max(0L, as.integer(k))
        G <- matrix(rnorm(n * (s + k)), n, s + k)
        Y <- apply_op_many(G)
        tr1 <- 0
        if (k > 0L) {
          Q <- qr.Q(qr(Y[, 1:k, drop = FALSE]))
          AQ <- apply_op_many(Q)
          tr1 <- sum(Q * AQ)

          G2 <- G[, (k+1):(s+k), drop = FALSE]
          G2 <- G2 - Q %*% (crossprod(Q, G2))
        } else {
          G2 <- G
        }

        AG2 <- apply_op_many(G2)
        tr2 <- mean(colSums(G2 * AG2))
        tr1 + tr2
      }
      

      degL <- as.numeric(Matrix::rowSums(S2))
      use_chol <- FALSE
      if (isTRUE(try_chol_first)) {
        M_try <- Diagonal(m) + rho_hat * S2
        cf_try <- try(Cholesky(M_try, LDL = FALSE, Imult = 0, perm = TRUE), silent = TRUE)
        if (!inherits(cf_try, "try-error")) {
          use_chol <- TRUE
          cf <- cf_try
        }
      }
      if (!use_chol) {
        use_chol <- (m <= chol_threshold)
        if (use_chol) {
          M <- Diagonal(m) + rho_hat * S2
          cf <- tryCatch(
            Cholesky(M, LDL = FALSE, Imult = 0, perm = TRUE),
            error = function(e) Cholesky(M + Diagonal(m, 1e-10), LDL = FALSE, Imult = 0, perm = TRUE)
          )
        }
      }
      
      if (use_chol) {
        A_apply_many <- function(B) as.matrix(solve(cf, B))
        A_apply      <- function(x) as.numeric(solve(cf, x))
      } else {
        M_mv  <- function(x) x + as.numeric(rho_hat * (S2 %*% x))
        Mdiag <- pmax(1 + rho_hat * degL, 1e-12)
        A_apply <- function(b) pcg_solve(M_mv, b, Mdiag = Mdiag, tol = pcg_tol, maxit = pcg_maxit)
        A_apply_many <- function(B) {

          out <- matrix(0.0, nrow(B), ncol(B))
          for (j in seq_len(ncol(B))) out[, j] <- A_apply(B[, j])
          out
        }
      }
      
      sigma2 <- sigma2_hat
      Oinv_apply <- function(x) A_apply(x) / sigma2
      H_mv       <- function(x) Oinv_apply(x) + K_diag * x
      Hdiag      <- (1 / sigma2) + K_diag
      Hsolve     <- function(b) pcg_solve(H_mv, b, Mdiag = pmax(Hdiag, 1e-12),
                                          tol = min(1e-8, pcg_tol), maxit = max(3000, pcg_maxit))
      Hsolve_many <- function(B) {
        out <- matrix(0.0, nrow(B), ncol(B))
        for (j in seq_len(ncol(B))) out[, j] <- Hsolve(B[, j])
        out
      }
      S2_many <- function(B) as.matrix(S2 %*% B)
      

      a   <- A_apply(delta)
      Sa  <- as.numeric(S2 %*% a)
      q0  <- sum(delta * a)
      q1  <- sum(a * Sa)
      

      op_t0   <- function(B) Hsolve_many(A_apply_many(B))
      op_t2   <- function(B) Hsolve_many(A_apply_many(S2_many(A_apply_many(B))))
      op_SS   <- function(B) A_apply_many(S2_many(A_apply_many(S2_many(B))))
      op_SS_H <- function(B) Hsolve_many(A_apply_many(S2_many(A_apply_many(S2_many(A_apply_many(B))))))
      

      op_r0 <- function(B) Hsolve_many(A_apply_many(Hsolve_many(A_apply_many(B))))
      op_r1 <- function(B) Hsolve_many(A_apply_many(Hsolve_many(A_apply_many(S2_many(A_apply_many(B))))))
      op_r2 <- function(B) {
        AB  <- A_apply_many(B)
        y1  <- Hsolve_many(A_apply_many(S2_many(AB)))
        Hsolve_many(A_apply_many(S2_many(A_apply_many(y1))))
      }
      

      set.seed(seed)

      
      
      t0   <- trace_hutchpp(op_t0,   m, s = n_mc, k = n_subspace)
      t2   <- trace_hutchpp(op_t2,   m, s = n_mc, k = n_subspace)
      t_SS <- trace_hutchpp(op_SS,   m, s = n_mc, k = n_subspace)
      t_SS_H <- trace_hutchpp(op_SS_H, m, s = n_mc, k = n_subspace)
      
      compute_traces_exact <- function(S2, sigma2, rho, K_diag) {
        stopifnot(inherits(S2, "sparseMatrix"))
        m <- nrow(S2)
        

        M  <- Diagonal(m) + rho * S2
        cfM <- Cholesky(M, LDL = FALSE, Imult = 0, perm = TRUE)
        A_full <- as.matrix(solve(cfM, Diagonal(m)))
        

        H_full <- (1 / sigma2) * A_full + diag(K_diag, m)
        R <- chol(H_full)
        
        solveH <- function(B) backsolve(R, forwardsolve(t(R), B))
        


        t0 <- sum(diag(solveH(A_full)))
        

        AS2A <- A_full %*% (as.matrix(S2) %*% A_full)
        t2   <- sum(diag(solveH(AS2A)))
        

        T_SS <- AS2A %*% as.matrix(S2)
        t_SS <- sum(diag(T_SS))
        

        T_SSH <- T_SS %*% A_full
        t_SS_H <- sum(diag(solveH(T_SSH)))
        
        list(t0 = t0, t2 = t2, t_SS = t_SS, t_SS_H = t_SS_H)
      }
      if(exact_trace)
      {
        res_trace <- compute_traces_exact(S2 , sigma2_hat ,rho_hat ,K_diag)
        
        t0 <- res_trace[[1]]
        t2<- res_trace[[2]]
        t_SS<- res_trace[[3]]
        t_SS_H <-  res_trace[[4]]
      }

      
      r0 <- r1 <- r2 <- 0
      if (isTRUE(use_nonfrozen_hessian)) {
        r0_raw <- trace_hutchpp(op_r0, m, s = n_mc, k = n_subspace)
        r1_raw <- trace_hutchpp(op_r1, m, s = n_mc, k = n_subspace)
        r2_raw <- trace_hutchpp(op_r2, m, s = n_mc, k = n_subspace)

        main_tau   <- abs((q0 + t0) / (2 * sigma2))
        main_cross <- abs((q1 + t2) / (2 * sigma2))
        main_rho   <- abs( 0.5 * ( t_SS - (2/sigma2) * (sum(delta * A_apply(as.numeric(S2 %*% A_apply(as.numeric(S2 %*% a))))) ) - (2/sigma2) * t_SS_H ) )


        
        if(exact_trace)
        {


          M   <- Diagonal(m) + rho_hat * S2
          cfM <- Cholesky(M, LDL = FALSE, Imult = 0, perm = TRUE)
          A_full <- as.matrix(solve(cfM, Diagonal(m)))
          
          H_full <- (1 / sigma2_hat) * A_full + diag(K_diag, m)
          R <- chol(H_full)
          
          solveH <- function(B) backsolve(R, forwardsolve(t(R), B))
          

          HA <- solveH(A_full)
          r0_raw <- sum(HA * t(HA))
          

          AS2A <- A_full %*% (S2 %*% A_full)

          T1 <- solveH(AS2A)
          T2 <- A_full %*% T1
          T3 <- solveH(T2)
          r1_raw <- sum(diag(T3))
          r2_raw <- sum(T1 * t(T1))
          
          
        }
        
        r0 <- r0_raw 
        r1 <- r1_raw 
        r2 <- r2_raw 
      }
      

      d2Q_tau2   <- -(q0 + t0)/(2 * sigma2) + (r0)/(2 * sigma2^2)




      AaSa       <- A_apply(as.numeric(S2 %*% a))
      d_tmp      <- as.numeric(S2 %*% AaSa)
      Ad_tmp     <- A_apply(d_tmp)
      u_SS_delta <- sum(delta * Ad_tmp)
      
      d2Q_rho2   <-  0.5 * ( t_SS - (2/sigma2) * u_SS_delta - (2/sigma2) * t_SS_H + (r2)/(sigma2^2) )
      d2Q_taurho <- -(q1 + t2)/(2 * sigma2) + (r1)/(2 * sigma2^2)
      
      HESS <- matrix(c(d2Q_tau2, d2Q_taurho,
                       d2Q_taurho, d2Q_rho2), 2, 2, byrow = TRUE)
      J <- -(HESS); J <- 0.5 * (J + t(J))
      

      lam <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
      if (any(!is.finite(lam))) stop("Information matrix has non-finite eigenvalues (trace estimates too noisy or PCG too loose).")
      floor_val <- info_floor_rel * mean(diag(J))
      if (min(lam) < floor_val) {
        J <- J + diag(floor_val - min(lam) + 1e-6, 2)
      }
      
      Vtheta <- tryCatch(solve(J), error = function(e) MASS::ginv(J, tol = 1e-12))
      se_tau <- sqrt(max(Vtheta[1,1], 0))
      se_rho <- sqrt(max(Vtheta[2,2], 0))

      G_orig <- matrix(c(sigma2_hat, 0,
                         0,          1), 2, 2, byrow = TRUE)
      V_orig <- G_orig %*% Vtheta %*% t(G_orig)
      se_sigma2 <- sqrt(max(V_orig[1,1], 0))
      

      kappa_hat <- sigma2_hat * rho_hat
      G_kappa <- matrix(c(sigma2_hat, 0,
                          kappa_hat,  sigma2_hat), 2, 2, byrow = TRUE)
      V_kappa <- G_kappa %*% Vtheta %*% t(G_kappa)
      se_kappa <- sqrt(max(V_kappa[2,2], 0))
      

      alpha <- 1 - conf_level; z <- qnorm(1 - alpha/2)
      
      tau_hat <- log(sigma2_hat)
      ci_tau  <- c(tau_hat - z * se_tau, tau_hat + z * se_tau)
      ci_sigma2 <- exp(ci_tau)

      

      x   <- (rho_hat - rho_lo) / (rho_hi - rho_lo)
      x   <- min(max(x, 1e-10), 1 - 1e-10)
      psi_hat <- qlogis(x)
      dpsi_drho <- 1 / ((rho_hi - rho_lo) * x * (1 - x))
      se_psi <- abs(dpsi_drho) * se_rho
      ci_psi <- c(psi_hat - z * se_psi, psi_hat + z * se_psi)
      inv_map <- function(psi) rho_lo + (rho_hi - rho_lo) * plogis(psi)
      ci_rho <- pmin(rho_hi, pmax(rho_lo, inv_map(ci_psi)))


      ci_kappa <- c(kappa_hat - z * se_kappa, kappa_hat + z * se_kappa)
      
 
      out <- list(
        estimate      = c(sigma2 = sigma2_hat, rho = rho_hat, kappa = kappa_hat),
        se            = c(se_sigma2 = se_sigma2, se_rho = se_rho, se_kappa = se_kappa),
        ci_sigma2     = ci_sigma2,
        ci_rho        = ci_rho,
        ci_kappa      = ci_kappa,
        conf_level    = conf_level,
        cov_tau_rho   = Vtheta,
        cov_original  = V_orig,
        cov_sigma2_kappa = V_kappa,
        notes = character(0)
      )
      

      condJ <- max(eigen(J, symmetric = TRUE, only.values = TRUE)$values) /
        max(min(eigen(J, symmetric = TRUE, only.values = TRUE)$values), 1e-16)
      note_vec <- c()
      if (condJ > 1e6 || sigma2_hat < 1e-3 || (rho_hi - rho_hat) < 0.05*(rho_hi - rho_lo)) {
        note_vec <- c(note_vec,
                      "Possible weak identification: small sigma^2 and/or rho near boundary; prefer kappa= sigma^2*rho or profile-LR for rho."
        )
      }
      if (isTRUE(use_nonfrozen_hessian) && (abs(r0) + abs(r1) + abs(r2)) < 1e-12) {
        note_vec <- c(note_vec, "Nonfrozen Hessian corrections were negligible or fully damped.")
      }
      out$notes <- note_vec
      
      if (return_internals) {
        out$internals <- list(
          q0 = q0, q1 = q1, t0 = t0, t2 = t2, t_SS = t_SS,
          u_SS_delta = u_SS_delta, t_SS_H = t_SS_H,
          r0 = r0, r1 = r1, r2 = r2,
          HESS = HESS, J = J, condJ = condJ,
          rho_bounds = c(rho_lo, rho_hi),
          n_mc = n_mc, n_subspace = n_subspace
        )
      }
      out
    }
  }
  

}


