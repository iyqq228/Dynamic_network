## ---- Parallel dependencies ---
library(survival)
library(Matrix)
source("function_FN.R", encoding = "UTF-8")
Sys.time()
{
library(foreach)
library(doSNOW)
library(doRNG)
library(data.table)


set.seed(20250917)
  params <- expand.grid(
    replicate = 1:20,
    n         = c(20),
    T_max     = c(5),
    sigma_edge= c(0.5),
    rho_edge  = c(0.1),
    beta1     = c(0.8),
    beta2     = c(-0.5),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )


num_cores <- 20
cl <- makeCluster(num_cores)
registerDoSNOW(cl)


n_jobs <- nrow(params)
pb <- txtProgressBar(min = 0, max = n_jobs, style = 3)
progress <- function(m) setTxtProgressBar(pb, m)
opts <- list(progress = progress)

p_max <- 2
beta_cols <- paste0("beta_", seq_len(p_max))
ii <- 2


res_list <- foreach(ii = seq_len(n_jobs),
                    .packages = c("data.table","Matrix","survival"),
                    .options.snow = list(
                      preschedule = FALSE,
                      progress    = progress
                    )

) %dorng% {
                                  

                                  pa <- params[ii, ]

                                  num_nodes <- params[ii,]$n
                                  T_max <- params[ii,]$T_max
                                  gamma_effect <- 0.5
                                  data_sim <- simulate_joint_network_data_tv(
                                    n = num_nodes, T_max = T_max,
                                    alpha_true=rep(0, num_nodes),
                                    beta = c(0.8, -0.5),
                                    eta  = gamma_effect*c(0),
                                    edge_frailty=TRUE, 
                                    h0=0.1,
                                    sigma_edge=params[ii,]$sigma_edge,
                                    rho_edge=params[ii,]$rho_edge,
                                    Z_edge_gen = "rand",
                                    tv_breaks = 1,
                                    tv_edge_amp = 0.2,
                                    tv_node_amp = 0.3,
                                    seed = ii*10,
                                    edge_method = "node-factor",
                                    ties = TRUE,
                                    ties_bins = T_max*100
                                  )
                                  

                                   
                                  fit <- fit_edge_frailty_ppl_block_pcg_tv_fast(
                                    data_sim = data_sim,
                                    max_outer = 50,
                                    beta_maxit = 50,
                                    sigma_fix0 = FALSE,
                                    sigma2_ini = (params[ii,]$sigma_edge)^2 ,
                                    rho_ini = params[ii,]$rho_edg,
                                    tol = 5e-4,
                                    damping = 0.99,
                                    chol_exact_threshold = 5000,
                                    ties_method = "efron"
                                  )
                                  df<-data_sim$df
                            
                    
                                  {

                                    scan <- build_tv_scanline_index(data_sim$df, zedge_pat="^Zedge_", zi1_pat="^Z1i_")
                                    n    <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
                                    ops  <- make_linegraph_ops(n, data_sim$pairs)
                                    S2   <- make_S2(n, data_sim$pairs)
               
                                    beta_hat   <- fit$beta_all
                                    delta_hat  <- fit$delta
                                    sigma2_hat <- fit$sigma^2
                                    rho_hat    <- fit$rho
                                    
                

                                    ties_method  <- "efron"
                                    kern  <- make_scan_kernel(scan, fit$beta_all,ties_method = ties_method)
                                    lam <- compute_cumLambda(scan, fit$beta_all,delta = fit$delta, method = ties_method)

                                    Kdiag <- compute_K_diag_partial(
                                      scan, fit$beta_all, fit$delta,ties_method = ties_method)

                                    ci_out <- wald_ci_sigma2_rho(
                                      delta = fit$delta, n =params[ii,]$n, S2 = S2,
                                      pairs = data_sim$pairs,
                                      sigma2_hat =  fit$sigma^2, rho_hat = fit$rho,
                                      K_diag = Kdiag, chol_threshold = 5000
                                    )


                                    c_sigma <- ci_out$ci_sigma2
                                    c_rho <- ci_out$ci_rho
                                    estsd_sigma2 <- ci_out$se[1]
                                    ci_out$se[1]

                                    {
                                      inf <- cox_frailty_theta_vcov (
                                        scan, pairs = data_sim$pairs,
                                        beta_hat =  fit$beta_all, delta_hat = fit$delta,
                                        sand = TRUE,beta_adjust = TRUE, cluster = "edge",
                                        sigma2 = (fit$sigma)^2, rho = fit$rho,
                                        ties_method = "efron",
                                        beta_col_pattern  = "^Zedge_",
                                        gamma_col_pattern = "^Xi_plus_Xj_",
                                        S2 =S2
                                      )
                                      sqrt(diag(inf$S_inv))

                                      inf$se
                                      inf$ci
                                    }




                                    {
                                      par_map <- c(
                                        beta1  = "Zedge_1",
                                        beta2  = "Zedge_2"
                                      )
                                      

                                      get_se   <- function(nm) if (!is.null(inf$se) && nm %in% names(inf$se)) unname(inf$se[[nm]]) else NA_real_
                                      get_low  <- function(nm) if (!is.null(inf$ci) && nm %in% rownames(inf$ci)) unname(inf$ci[nm, "lower"]) else NA_real_
                                      get_high <- function(nm) if (!is.null(inf$ci) && nm %in% rownames(inf$ci)) unname(inf$ci[nm, "upper"]) else NA_real_
                                      

                                      se_ci_vec <- unlist(lapply(names(par_map), function(k) {
                                        nm <- par_map[[k]]
                                        setNames(
                                          c(get_se(nm), get_low(nm), get_high(nm)),
                                          paste0(c(k, k, k), c("_se", "_low", "_high"))
                                        )
                                      }), use.names = TRUE)
                                      
                                    }
                                  }
                                  



                                  

                                  bhat <- unname(fit$beta_all)


                                  
                                

                                  return(c(
                                    n                = params[ii,]$n,
                                    sigma            = params[ii,]$sigma_edge^2,
                                    rho              = params[ii,]$rho_edge,
                                    beta1_true  = params[ii,]$beta1,
                                    beta2_true  = params[ii,]$beta2,
                                    beta1            = bhat[1],
                                    beta2            = bhat[2],
                                    T_max            = params[ii,]$T_max,
                                    sigma            = as.numeric(fit$sigma)^2,
                                    sigma_low        = c_sigma[1],
                                    sigma_high       = c_sigma[2],
                                    rho_low          = c_rho[1],
                                    rho_high         = c_rho[2],
                                    estsd_sigma2     = ci_out$se[1],
                                    estsd_rho        = ci_out$se[2],
                                    rho              = as.numeric(fit$rho),
                                    se_ci_vec
                                  ))
                                  






                                  
                                
                                }

close(pb)


stopCluster(cl)

}
Sys.time()


# res_mat <- do.call(rbind, res_list)
res_mat_all <- data.frame(do.call(rbind, res_list))
###### Summarize the first result
{
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
  })
  
  summarise_res_v3 <- function(df, truths = NULL, truth_tbl = NULL, alpha = 0.05) {
    z <- qnorm(1 - alpha/2)
    keys <- c("n","sigma","rho","T_max")
    

    if (!is.null(truth_tbl)) {
      stopifnot(all(keys %in% names(truth_tbl)))
      df <- df %>% left_join(truth_tbl, by = keys)
    }
    

    gcol <- function(d, nm, na = NA_real_) {
      if (nm %in% names(d)) d[[nm]] else rep(na, nrow(d))
    }
    

    make_rows <- function(param, est, low, high, se_col, truth_col, truth_const = NULL) {
      se_vec <- if (!is.null(se_col) && se_col %in% names(df)) gcol(df, se_col) else {

        (gcol(df, high) - gcol(df, low)) / (2*z)
      }
      tru_vec <- if (!is.null(truth_col) && truth_col %in% names(df)) gcol(df, truth_col) else {
        if (!is.null(truth_const)) rep(truth_const, nrow(df)) else rep(NA_real_, nrow(df))
      }
      tibble(
        n = df$n, sigma = df$sigma, rho = df$rho, T_max = df$T_max,
        param = param,
        est   = gcol(df, est),
        low   = gcol(df, low),
        high  = gcol(df, high),
        se    = se_vec,
        true  = tru_vec
      )
    }
    

    truth_from_list <- function(name) {
      if (!is.null(truths) && name %in% names(truths)) truths[[name]] else NULL
    }
    
    long <- bind_rows(

      make_rows("sigma", "sigma.1", "sigma_low", "sigma_high",
                "estsd_sigma2.se_sigma2", truth_col = "sigma"),
      make_rows("rho",   "rho.1",   "rho_low",   "rho_high",
                "estsd_rho.se_rho", truth_col = "rho"),
      

      make_rows("beta1", "beta1", "beta1_low", "beta1_high",
                "beta1_se", truth_col = "beta1_true", truth_const = truth_from_list("beta1")),
      make_rows("beta2", "beta2", "beta2_low", "beta2_high",
                "beta2_se", truth_col = "beta2_true", truth_const = truth_from_list("beta2"))
    )
    

    summary_long <- long %>%
      group_by(n, sigma, rho, T_max, param) %>%
      summarise(
        true        = suppressWarnings(unique(true[!is.na(true)])[1]),
        est_mean    = mean(est, na.rm = TRUE),

        bias        = if (all(is.na(true))) NA_real_ else mean(est - true, na.rm = TRUE),

        bias_sigma  = if (all(is.na(true))) NA_real_
        else mean(sqrt(pmax(est,  0)) - sqrt(pmax(true, 0)), na.rm = TRUE),
        emp_sd      = sd(est, na.rm = TRUE),
        avg_est_se  = mean(se, na.rm = TRUE),
        cover_95    = if (all(is.na(true))) NA_real_ else mean(true >= low & true <= high, na.rm = TRUE),
        mean_ci_len = mean(high - low, na.rm = TRUE),
        n_eff       = sum(!is.na(est)),
        .groups = "drop"
      ) %>%

      mutate(bias = ifelse(param == "sigma", bias_sigma, bias)) %>%
      select(-bias_sigma) %>%
      arrange(n, sigma, rho, T_max, param)
    
    summary_wide <- summary_long %>%
      select(n, sigma, rho, T_max, param, bias, emp_sd, avg_est_se, cover_95) %>%
      pivot_wider(names_from = param, values_from = c(bias, emp_sd, avg_est_se, cover_95))
    
    
    list(long = summary_long, wide = summary_wide)
  }

  res <- summarise_res_v3(res_mat_all)
  

  head(res$wide)
  rr <- res$wide

  res$long
  
}

rr1<-rr
library(dplyr)
res_mat_all <- res_mat_all %>% mutate(
  sigma=sigma^2,
  sigma.1 = sigma.1^2
)

res_mat <- res_mat_all 

res_mat <- data.frame(res_mat)
sd((res_mat[,c("sigma.1")])^2)
sd(res_mat[,c("rho.1")])
mean(res_mat[,c("estsd_sigma2.se_sigma2")])
mean(res_mat[,c("estsd_rho.se_rho")])
library(dplyr)
results_summary <- res_mat %>%
  dplyr::group_by(n, sigma, rho) %>%
  dplyr::summarise(
    beta1 = mean(beta1, na.rm = TRUE),
    beta2   = mean(beta2,   na.rm = TRUE),
    sigma1 = mean(sigma.1,   na.rm = TRUE),
    rho1 = mean(rho.1,   na.rm = TRUE),
    .groups = "drop"
  )



############ Test results:
{
  
  #######First test the Type I error
  Sys.time()
  {
    library(foreach)
    library(doSNOW)     # Display the progress bar with doSNOW
    library(doRNG)      # Reproducible parallel execution
    library(data.table) # Fast aggregation with rbindlist; can be replaced by dplyr/tibble
    
    ## ---- Parameter grid (adjust as needed) ----
    set.seed(20250917)
    params <- expand.grid(
      replicate = 1:50,          # Number of replicates
      n         = c(20),          # Number of nodesc(20,30,40,60)
      T_max     = c(5,10),
      sigma_power = c(0.1),   ####Sigma values for the power setting: c(0.1, 0.15)
      # rho_sigma_power = c(0,0.2),  ######Rho values used to test power for sigma
      sigma_edge= c(0.5),       ####Sigma value used under the rho setting
      rho_power  = c(0.02,0.03),     ####Set to 0.1 when calculating power for sigma
      beta1     = c(0),  
      beta2     = c(0),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    
    ## ---- Parallel cluster ----
    num_cores <- 5
    cl <- makeCluster(num_cores)  # Windows/cross-platform
    registerDoSNOW(cl)
    
    
    ## ---- Progress bar ----
    n_jobs <- nrow(params)
    pb <- txtProgressBar(min = 0, max = n_jobs, style = 3)
    progress <- function(m) setTxtProgressBar(pb, m)
    opts <- list(progress = progress)
    
    p_max <- 2
    beta_cols <- paste0("beta_", seq_len(p_max))
    ii <- 6
    ## ---- Main parallel loop ----
    # doRNG ensures reproducibility (the global seed is set above)
    res_list <- foreach(ii = seq_len(n_jobs),
                        .packages = c("data.table","Matrix","survival"),
                        .options.snow = list(
                          preschedule = FALSE,   # Disable prescheduling and enable dynamic scheduling
                          progress    = progress # Progress bar
                        )
                        # .errorhandling = "pass"  # Return errors without interrupting execution
    ) %dorng% {
      
      # Extract the current parameters
      pa <- params[ii, ]
      
      
      num_nodes <- params[ii,]$n
      T_max <- params[ii,]$T_max
      
      ####Type I error for sigma
      {
      data_sim <- simulate_joint_network_data_tv(
        n = num_nodes, T_max = T_max,
        alpha_true=rep(0, num_nodes),
        beta = c(0.5, 0.5),        # Two edge covariates
        eta  = c(0, 0),   # Two-dimensional node covariates; their magnitude affects the sigma estimate
        edge_frailty=FALSE,
        h0=0.1,
        sigma_edge=0 ,
        rho_edge=0 ,
        Z_edge_gen = "rand",
        tv_breaks = 1,              # Five segments per edge (then merge in event times)
        tv_edge_amp = 0.2,          # Time-varying amplitude
        tv_node_amp = 0.3,
        seed = ii*100,
        edge_method = "node-factor",
        ties = TRUE,
        ties_bins = T_max*200
      )

        ## Construct the scan index and graph operators
        scan <- build_tv_scanline_index(data_sim$df, zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_")
        n    <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
        ops  <- make_linegraph_ops(n, data_sim$pairs)
        S2   <- make_S2(n, data_sim$pairs)
        res_beta <- beta_update_via_coxph_fit(
          X_rows = scan$X_rows,
          df = data_sim$df,
          scan =scan,
          delta  =rep(0,length(data_sim$pairs)),
          ties   = "efron" ,
          strata = rep(1L, nrow(data_sim$df))   
        )
       
        res_sigma <- score_test_sigma2_boundary(scan,
                                                beta_hat_null =res_beta$beta,   
                                                ties_method="efron")

        p1_I <- res_sigma$pval$imhof

        p1_I_ad <-   p1_I

      }

      {
        data_sim <- simulate_joint_network_data_tv(
          n = num_nodes, T_max = T_max,
          alpha_true=rep(0, num_nodes),
          beta = c(0.5, 0.5),        # Two edge covariates
          eta  = c(0, 0),   # Three-dimensional node covariates
          edge_frailty=TRUE,
          h0=0.1,
          sigma_edge=params[ii,]$sigma_power,
          rho_edge= 0.2  ,
          Z_edge_gen = "rand",
          tv_breaks = 1,              # Five segments per edge (then merge in event times)
          tv_edge_amp = 0.2,          # Time-varying amplitude
          tv_node_amp = 0.3,
          seed = ii*10,
          edge_method = "node-factor",
          ties = TRUE,
          ties_bins = T_max*100
        )

        ## Construct the scan index and graph operators
        scan <- build_tv_scanline_index(data_sim$df, zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_")
        n    <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
        ops  <- make_linegraph_ops(n, data_sim$pairs)
        S2   <- make_S2(n, data_sim$pairs)
        res_beta <- beta_update_via_coxph_fit(
          X_rows = scan$X_rows,
          df = data_sim$df,
          scan =scan,
          delta  =rep(0,length(data_sim$pairs)),
          ties   = "efron" ,
          strata = rep(1L, nrow(data_sim$df))   # Alternatively, pass a stratification variable such as df$type
        )
      
        kk <- compute_edge_M_and_K_strict(scan,
                                    beta = res_beta$beta,
                                    ties_method="efron")

        res_sigma <- score_test_sigma2_boundary(scan,
                                                beta_hat_null = res_beta$beta,
                                                ties_method="efron")
        
        p1_II <- res_sigma$pval$imhof

        p1_II_ad <-   p1_I
        res_sigma$pval
      }

      
      #######2. Test rho
      #####Test the Type I error rate when rho = 0
      {
        
        data_sim <- simulate_joint_network_data_tv(
          n = num_nodes, T_max = T_max,
          alpha_true=rep(0, num_nodes),
          beta = c(0.5, 0.5),        # Two edge covariates
          eta  = c(0, 0),   # Three-dimensional node covariates; why does setting this to 0 cause problems?
          edge_frailty=TRUE,
          h0=0.1,
          sigma_edge=params[ii,]$sigma_edge,
          rho_edge=0 ,
          Z_edge_gen = "rand",
          tv_breaks = 1,              # Five segments per edge (then merge in event times)
          tv_edge_amp = 0.2,          # Time-varying amplitude
          tv_node_amp = 0.3,
          seed = ii*1000,
          edge_method = "node-factor",
          ties = TRUE,
          ties_bins = T_max*100
        )

        fit_1 <- fit_edge_frailty_ppl_block_pcg_tv_fast(
          data_sim = data_sim,           # data_sim must contain df, pairs, and Z_i
          zi1_pattern = "^Z1i_",
          zi2_pattern = "^Z2i_",
          beta_ini = fit$beta_edge,
          max_outer = 2,
          beta_maxit = 50,
          sigma2_ini = params[ii,]$sigma_edge^2,
          rho_ini = 0,
          tol = 5e-5,
          damping = 0.95,
          chol_exact_threshold = 5000,
          ties_method = "efron"
        )
        fit <- fit_1
        
        
        rho_hat_null <- fit_1$rho
        scan <- build_tv_scanline_index(data_sim$df, zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_")
        n    <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
        ops  <- make_linegraph_ops(n, data_sim$pairs)
        S2   <- make_S2(n, data_sim$pairs)
 
   
        res_rho <- score_test_rho0_theory(
          scan, S2,
          beta_hat = fit$beta_all,
          delta_hat = fit$delta,
          sigma2_hat = fit$sigma^2,
          ties_method = "efron",
          Kdiag_fn = "full",
          info = "observed",           # Use the observed curvature
          pairs = data_sim$pairs,
          emp_VAR = FALSE,
          orthogonalize = FALSE
        )
        res_rho$pval
        p3_I <- res_rho$pval
        

        
        
        {
          
          ## 3a) Interval for sigma^2 (closed-form approximation)
          ties_method  <- "efron"
          kern  <- make_scan_kernel(scan, fit_1$beta_all,ties_method = ties_method)
          Kdiag <- compute_K_diag_partial(
            scan, fit_1$beta_all, fit_1$delta,ties_method = ties_method)
          
          ci_out <- wald_ci_sigma2_rho(
            delta = fit_1$delta, n =params[ii,]$n, S2 = S2,
            pairs = data_sim$pairs,
            sigma2_hat =  fit_1$sigma^2, rho_hat = fit_1$rho,
            K_diag = Kdiag, chol_threshold = 5000
          )
          
          
          c_sigma <- ci_out$ci_sigma2   # 95% CI for sigma^2
          c_rho <- ci_out$ci_rho      # 95% CI for rho (projected onto the feasible region)
          p2_I<- 1
          if(min( c_rho)>0.000001)
          {
            p2_I <- 0.001
          }
          
          ##########One-sided test
          {
            rho_hat <- as.numeric(ci_out$estimate["rho"])
            se_rho  <- as.numeric(ci_out$se["se_rho"])
            
            if (!is.finite(se_rho) || se_rho <= 0) {
              p2_I <- NA_real_
            } else if (rho_hat <= 0) {
              p2_I <- 1
            } else {
              Z_rho <- rho_hat / se_rho
              p2_I <- 1 - pnorm(Z_rho)
            }
          }
          
          
          
        }
        
        
      }
      
      
      ####Test power for rho
      {
        data_sim <- simulate_joint_network_data_tv(
          n = num_nodes, T_max = T_max,
          alpha_true=rep(0, num_nodes),
          beta = c(0.5, 0.5),        # Two edge covariates
          eta  = c(0,0),   # Three-dimensional node covariates
          edge_frailty=TRUE,
          h0=0.1,
          sigma_edge=params[ii,]$sigma_edge,
          rho_edge=params[ii,]$rho_power ,
          Z_edge_gen = "rand",
          tv_breaks = 1,              # Five segments per edge (then merge in event times)
          tv_edge_amp = 0.2,          # Time-varying amplitude
          tv_node_amp = 0.3,
          seed = ii*1000,
          edge_method = "node-factor",
          ties = TRUE,
          ties_bins = T_max*100
        )
        
        
        fit_1 <- fit_edge_frailty_ppl_block_pcg_tv_fast(
          data_sim = data_sim,           # data_sim must contain df, pairs, and Z_i
          zi1_pattern = "^Z1i_",
          zi2_pattern = "^Z2i_",
          max_outer = 25,
          beta_maxit = 50,
          sigma2_ini = params[ii,]$sigma_edge^2,
          rho_ini = 0.001,
          tol = 1e-4,
          damping = 0.9,
          chol_exact_threshold = 5000,
          ties_method = "efron"
        )
        
        fit <- fit_1
        scan <- build_tv_scanline_index(data_sim$df, zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_")
        n    <- if (!is.null(data_sim$Z_i)) nrow(data_sim$Z_i) else max(data_sim$pairs)
        ops  <- make_linegraph_ops(n, data_sim$pairs)
        S2   <- make_S2(n, data_sim$pairs)

        res_rho <- score_test_rho0_theory(
          scan, S2,
          beta_hat = fit$beta_all,
          delta_hat = fit$delta,
          sigma2_hat = fit$sigma^2,
          ties_method = "efron",
          Kdiag_fn = "full",
          info = "observed",           # Use the observed curvature
          pairs = data_sim$pairs,
          emp_VAR = FALSE,
          orthogonalize = FALSE
        )
        res_rho$pval
        p3_II <- res_rho$pval
        
 
        ####Wald
        {
          
          ## 3a) Interval for sigma^2 (closed-form approximation)
          ties_method  <- "efron"
          kern  <- make_scan_kernel(scan, fit_1$beta_all,ties_method = ties_method)
          Kdiag <- compute_K_diag_partial(
            scan, fit_1$beta_all, fit_1$delta,ties_method = ties_method)
          
          ci_out <- wald_ci_sigma2_rho(
            delta = fit_1$delta, n =params[ii,]$n, S2 = S2,
            pairs = data_sim$pairs,
            sigma2_hat =  fit_1$sigma^2, rho_hat = fit_1$rho,
            K_diag = Kdiag, chol_threshold = 5000
          )
          
          
          c_sigma <- ci_out$ci_sigma2   # 95% CI for sigma^2
          c_rho <- ci_out$ci_rho      # 95% CI for rho (projected onto the feasible region)
          p2_II<- 1
          if(min( c_rho)>0.0001)
          {
            p2_II <- 1e-6
          }
          
          ##########One-sided test
          {
            rho_hat <- as.numeric(ci_out$estimate["rho"])
            se_rho  <- as.numeric(ci_out$se["se_rho"])
            
            if (!is.finite(se_rho) || se_rho <= 0) {
              p2_II <- NA_real_
            } else if (rho_hat <= 0) {
              p2_II <- 1
            } else {
              Z_rho <- rho_hat / se_rho
              p2_II <- 1 - pnorm(Z_rho)
            }
          }
          
        }
      }
      
      

      
      # 3) Organize the output
      bhat <- unname(fit$beta_all)
 

      return(c(n=params[ii,]$n,
               sigma=params[ii,]$sigma_power,
               T_max = params[ii,]$T_max,
               sigma_rhoset = params[ii,]$sigma_edge,
               rho=params[ii,]$rho_power,
               beta1=bhat[1],
               beta2=bhat[2],
               gamma1 = bhat[3],
               gamma2 = bhat[4],
               sigma=as.numeric(fit$sigma),
               rho.1 = as.numeric(rho_hat_null),
               sigma_p_I =p1_I,
               sigma_p_II = p1_II,
               sigma_p_I_ad =p1_I_ad,
               sigma_p_II_ad = p1_II_ad,
               rho_p_I = p2_I,
               rho_p_noeff_I =p3_I,
               rho_p_II = p2_II,
               rho_p_noeff_II =p3_II,
               rho=as.numeric(fit$rho)))
      

    }
    
    close(pb)
    
    ## ---- Cleanup and aggregation ----
    stopCluster(cl)
    

  }
  Sys.time()
  
  
  
  }
