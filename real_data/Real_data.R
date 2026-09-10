library(readxl)

time_data <- read_excel("time1.xlsx")   ########event time
cov_data <- read_excel("cov.xlsx")   ########Covariate
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

##### 1. Data process
{
  
  
{
  calculate_avg <- function(value1, value2) {

    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    
    if (is.na(value1) & is.na(value2)) {
      return(mean(cov_data$body_size, na.rm = TRUE))
    } else if (is.na(value1)) {
      return(value2) 
    } else if (is.na(value2)) {
      return(value1)
    } else {
      return(mean(c(value1, value2)))
    }
  }
  
  
  

  interaction_types <- c("Q", "N-N", "N-C", "N-F", "C-C", "C-F", "F-F")
  

  covariate_dfs <- list()
  
 
  calculate_avg <- function(value1, value2) {
    if (is.na(value1) & is.na(value2)) {
      return(mean(cov_data$body_size, na.rm = TRUE)) 
    } else if (is.na(value1)) {
      return(value2) 
    } else if (is.na(value2)) {
      return(value1) 
    } else {
      return(mean(c(value1, value2))) 
    }
  }
  

  for (day in c(1,15,25,40)) {
  
    df <- data.frame(matrix(NA, nrow = nrow(time_data), ncol = 8))
  
    if (day <= 11) {
      period_column <- "group_period1"
    } else if (day <= 21) {
      period_column <- "group_period2"
    } else if (day <= 31) {
      period_column <- "group_period3"
    } else {
      period_column <- "group_period4"
    }
    
   
    Ant_pair <- time_data$`Ant Pair`
    for (i in 1:nrow(time_data)) {
      pair <- unlist(strsplit(as.character(Ant_pair[i]), "-"))
      ant1 <- gsub("Ant", "", pair[1])  
      ant2 <- gsub("Ant", "", pair[2])  
      
    
      group_ant1 <- cov_data[which(cov_data$tag_id == ant1), period_column]
      group_ant2 <- cov_data[which(cov_data$tag_id == ant2), period_column]
      
     
      if ("Q" %in% c(group_ant1, group_ant2)) {
        interaction_column <- which(interaction_types == "Q")
      } else {
     
        order <- c("N", "C", "F")
    
        interaction <- ifelse(match(group_ant1, order) < match(group_ant2, order), 
                              paste(group_ant1, group_ant2, sep = "-"), 
                              paste(group_ant2, group_ant1, sep = "-"))
        
        interaction_column <- which(interaction_types == interaction)
      }
      df[i, interaction_column] <- 1

      body_size_1 <- cov_data[which(cov_data$tag_id == ant1), "body_size"][[1]]
      body_size_2 <- cov_data[which(cov_data$tag_id == ant2), "body_size"][[1]]
      
      age_1 <- cov_data[which(cov_data$tag_id == ant1), "age"][[1]]
      age_2 <- cov_data[which(cov_data$tag_id == ant2), "age"][[1]]
      
      body_size_avg <- calculate_avg(body_size_1, body_size_2)
      age_avg <- calculate_avg(age_1, age_2)

      df[i, 8] <- body_size_avg
      df[i, 9] <- age_avg
    }
    
    for (i in 1:nrow(df)) {
      if (any(!is.na(df[i, 1:7]))) {
        df[i, 1:7][is.na(df[i, 1:7])] <- 0 
      }
    }
    

    df$`Ant Pair` <- Ant_pair

    colnames(df) <- c("Q","N-N", "N-C", "N-F", "C-C", "C-F", "F-F", "avg_body_size", "avg_age", "Ant Pair")
    covariate_dfs[[day]] <- df
  }
  
  
  

}

{
 
  binary_matrix <- ifelse(!is.na(time_data[, -1]), 1, 0)
  
  ss1_common <- covariate_dfs[[1]]
  ss2_common <- covariate_dfs[[15]]
  ss3_common <- covariate_dfs[[25]]
  ss4_common <- covariate_dfs[[40]]
  normalize_z_score <- function(df) {
    
    df$avg_body_size <- (df$avg_body_size - mean(df$avg_body_size, na.rm = TRUE)) /
      sd(df$avg_body_size, na.rm = TRUE)
    
   
    df$avg_age <- (df$avg_age - mean(df$avg_age, na.rm = TRUE)) /
      sd(df$avg_age, na.rm = TRUE)
    
    return(df)
  }
  
 
  ss1_common <- normalize_z_score(ss1_common)
  ss2_common <- normalize_z_score(ss2_common)
  ss3_common <- normalize_z_score(ss3_common)
  ss4_common <- normalize_z_score(ss4_common)
}
  #####

  {
    num_data <- time_data[ , -1]
    
    ### 3. Keep only N-N
    {
      #keep_idx <-  which(ss1_common[ , "N-F"] == 1)
      keep_idx <- c(1:nrow(time_data))
    }
    
    # kk <- time_data[161,]
    # kk <- ss4_common[161,]
    # Data retained after filtering
    time_data <- time_data[keep_idx,]
    ss1_common <- ss1_common[keep_idx,]
    ss2_common <- ss2_common[keep_idx,]
    ss3_common <- ss3_common[keep_idx,]
    ss4_common <- ss4_common[keep_idx,]
    rownames(ss1_common) <- NULL
    rownames(ss2_common) <- NULL
    rownames(ss3_common) <- NULL
    rownames(ss4_common) <- NULL
  }



wide_to_event_long_with_death <- function(df_wide,
                                          pair_col   = "Ant Pair",
                                          day_prefix = "Day") {
  stopifnot(pair_col %in% names(df_wide))
  

  day_cols <- grep(paste0("^", day_prefix, "\\d+$"), names(df_wide), value = TRUE)
  if (!length(day_cols)) stop("No columns matching 'Day1, Day2, ...' were found.")
  day_idx  <- readr::parse_number(day_cols)
  ord      <- order(day_idx)
  day_cols <- day_cols[ord]
  day_idx  <- day_idx[ord]
  n_days   <- max(day_idx)
  

  get_first_na_day <- function(x) {
    pos <- which(is.na(x))
    if (length(pos)) pos[1] else (n_days + 1L)
  }
  death_day_vec <- apply(df_wide[, day_cols, drop = FALSE], 1, get_first_na_day)
  

  df_long <- df_wide %>%
    mutate(.row_id = row_number(),
           .death_day = death_day_vec) %>%
    pivot_longer(all_of(day_cols),
                 names_to = "day_label", values_to = "count") %>%
    mutate(day = readr::parse_number(day_label))
  

  df_interact_days <- df_long %>%
    filter(day < .death_day, !is.na(count), count > 0)
  

  df_interact_days <- df_interact_days %>%
    separate({{pair_col}}, into = c("ant1","ant2"), sep = "-", remove = FALSE)
  

  interactions <- df_interact_days %>%
    mutate(count_int = as.integer(round(count))) %>%
    tidyr::uncount(weights = count_int, .remove = FALSE, .id = "k") %>%
    group_by(.row_id, `Ant Pair`, ant1, ant2, day) %>%
    mutate(
      event_in_day = k,
      n_in_day     = n(),
      event_time   = (day - 1) + event_in_day / (n_in_day + 1),
      event_type   = "interaction",
      event        = 1L,
      death        = 0L
    ) %>%
    ungroup() %>%
    dplyr::select(`Ant Pair`, ant1, ant2, day, event_in_day,
                  event_time, event_type, event, death, .row_id)
  

  death_rows <- df_wide %>%
    mutate(.row_id = row_number(),
           .death_day = death_day_vec) %>%
    transmute(
      .row_id,
      `Ant Pair`,
      ant1 = sub("-.*$", "", .data[[pair_col]]),
      ant2 = sub("^.*?-", "", .data[[pair_col]]),
      day  = ifelse(.death_day <= n_days, .death_day, n_days),
      event_in_day = NA_integer_,
      event_time   = ifelse(.death_day <= n_days, .death_day - 1, n_days),
      event_type   = "death",
      event        = 1L,
      death        = 1L
    )
  
  bind_rows(interactions, death_rows) %>%
    arrange(.row_id, event_time, match(event_type, c("interaction","death"))) %>%
    dplyr::select(-.row_id)
}



augment_event_long_with_ag_columns <- function(event_long) {
  req <- c("Ant Pair", "event_time", "event_type")
  stopifnot(all(req %in% names(event_long)))
  

  key_df <- event_long %>%
    transmute(
      `Ant Pair`,
      a1 = sub("-.*$", "", `Ant Pair`),
      a2 = sub("^.*?-", "", `Ant Pair`),
      lo = pmin(a1, a2),
      hi = pmax(a1, a2),
      pair_key = paste0(lo, "-", hi)
    ) %>%
    distinct(pair_key) %>%
    arrange(pair_key) %>%
    mutate(id = as.integer(factor(pair_key, levels = pair_key))) %>%
    dplyr::select(pair_key, id)
  
  out <- event_long %>%
    mutate(
      a1 = sub("-.*$", "", `Ant Pair`),
      a2 = sub("^.*?-", "", `Ant Pair`),
      lo = pmin(a1, a2),
      hi = pmax(a1, a2),
      pair_key = paste0(lo, "-", hi)
    ) %>%
    left_join(key_df, by = "pair_key") %>%
    arrange(id, as.numeric(event_time),
            match(event_type, c("interaction", "death"))) %>%
    group_by(id) %>%
    mutate(
      stop  = as.numeric(event_time),
      start = dplyr::lag(stop, default = 0),
      event = ifelse(event_type == "interaction", 1L, 0L)
    ) %>%
    ungroup() %>%
    dplyr::select(-a1, -a2, -lo, -hi, -pair_key) %>%
    relocate(id, start, stop, event, .after = last_col())
  
  out
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

locf_covariates_within_edge <- function(event_long, extra_fill_cols = character()) {
  stopifnot(all(c("id","start","stop") %in% names(event_long)))
  

  non_fill_cols <- intersect(c(
    "Ant Pair","ant1","ant2",
    "id","day","event_time","event_type",
    "event_in_day","n_in_day",
    "start","stop","event","death"
  ), names(event_long))
  

  default_fill_patterns <- c("^Zedge_", "^body_size_", "^age_", "^group_",
                             "^g1$", "^g2$", "^pair_cat$", "^period$")
  

  default_fill_cols <- unique(unlist(lapply(default_fill_patterns, function(pat) {
    grep(pat, names(event_long), value = TRUE)
  })))
  

  fill_cols <- setdiff(unique(c(default_fill_cols, extra_fill_cols)), non_fill_cols)
  fill_cols <- intersect(fill_cols, names(event_long))
  

  if (length(fill_cols) == 0) {
    message("No covariate columns matched for filling; returning the original data.")
    out <- event_long %>% arrange(id, start, stop)
    return(list(data = out, na_cols = names(out)[sapply(out, function(x) any(is.na(x)))]))
  }
  

  out <- event_long %>%
    arrange(id, start, stop) %>%
    group_by(id) %>%
    tidyr::fill(all_of(fill_cols), .direction = "down") %>%
    ungroup()
  

  na_cols <- names(out)[sapply(out, function(x) any(is.na(x)))]
  




  
  list(data = out, na_cols = na_cols)
}



#   tag_id, body_size, age, group_period1, group_period2, group_period3, group_period4
add_covariates_to_events <- function(event_long, nodes_df) {
  
  library(dplyr); library(tidyr); library(stringr); library(readr)
  

  nodes_df <- nodes_df %>% mutate(tag_id = suppressWarnings(as.integer(tag_id)))
  

  day_to_period <- function(d) {
    d <- pmax(pmin(as.integer(d), 41L), 1L)
    cut(d, breaks = c(0,11,21,31,41),
        labels = c(1L,2L,3L,4L), right = TRUE, include.lowest = TRUE) |> as.integer()
  }
  

  nodes_groups_long <- nodes_df %>%
    rename(group1 = group_period1, group2 = group_period2,
           group3 = group_period3, group4 = group_period4) %>%
    pivot_longer(starts_with("group"), names_to = "period_lab", values_to = "group") %>%
    mutate(period = readr::parse_number(period_lab)) %>%
    dplyr::select(tag_id, body_size, age, period, group)
  
  nodes1 <- nodes_groups_long %>%
    rename(tag_id_1 = tag_id, body_size_1 = body_size, age_1 = age, group_1 = group)
  nodes2 <- nodes_groups_long %>%
    rename(tag_id_2 = tag_id, body_size_2 = body_size, age_2 = age, group_2 = group)
  

  ev <- event_long %>%
    mutate(
      id1 = suppressWarnings(as.integer(stringr::str_extract(ant1, "\\d+"))),
      id2 = suppressWarnings(as.integer(stringr::str_extract(ant2, "\\d+"))),
      period = day_to_period(day)
    ) %>%
    left_join(nodes1, by = c("id1" = "tag_id_1", "period" = "period")) %>%
    left_join(nodes2, by = c("id2" = "tag_id_2", "period" = "period"))
  

  normalize_group <- function(x) {
    x <- toupper(as.character(x))
    x[!x %in% c("N","C","F")] <- NA_character_
    x
  }
  ev <- ev %>%
    mutate(g1 = normalize_group(group_1),
           g2 = normalize_group(group_2))
  

  rank_map <- c(N = 1L, C = 2L, F = 3L)
  undirected_pair <- function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_character_)
    if (a == b) return(paste0(a, "-", b))
    if (rank_map[[a]] < rank_map[[b]]) paste0(a, "-", b) else paste0(b, "-", a)
  }
  ev$pair_cat <- mapply(undirected_pair, ev$g1, ev$g2, USE.NAMES = FALSE)
  ev$edge_group_na <- as.integer(is.na(ev$pair_cat))
  

  cats <- c("N-N","N-C","N-F","C-C","C-F","F-F")
  for (k in cats) {
    ev[[paste0("Zedge_", gsub("-", "", k))]] <- as.integer(ev$pair_cat == k)
  }
  

  ev %>%
    relocate(period, .after = day) %>%
    relocate(g1, g2, pair_cat, edge_group_na, .after = period) %>%
    relocate(body_size_1, age_1, body_size_2, age_2, .after = pair_cat)
}


rank_map <- c(N=1L, C=2L, F=3L)
undirected_pair <- function(a,b){
  if (is.na(a)||is.na(b)) return(NA_character_)
  if (a==b) return(paste0(a,"-",b))
  if (rank_map[[a]] < rank_map[[b]]) paste0(a,"-",b) else paste0(b,"-",a)
}


event_long <- wide_to_event_long_with_death(time_data)
event_long <- augment_event_long_with_ag_columns(event_long)
safe_num <- function(x) {

  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x[x %in% c("NA", "na", "Na", "")] <- NA_character_
  }
  suppressWarnings(as.numeric(x))
}

cov_data <- cov_data %>%
  mutate(
    body_size = safe_num(body_size),
    age       = safe_num(age)
  ) %>%
  mutate(
    body_size = as.numeric(scale(body_size)),
    age       = as.numeric(scale(age))
  )

event_with_cov <- add_covariates_to_events(event_long, cov_data)

res <- locf_covariates_within_edge(event_with_cov)
event_with_cov<- res$data

anyNA(event_with_cov)              # TRUE/FALSE
sum(is.na(event_with_cov))        

dplyr::glimpse(event_with_cov)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


event_with_cov <- event_with_cov %>%
  mutate(
    
    across(where(is.character), ~ replace(., . %in% c("NA","na","Na"), NA_character_)),

    across(where(is.factor), ~ {
      x <- as.character(.)
      x[x %in% c("NA","na","Na")] <- NA_character_
      factor(x)
    })
  )


event_with_cov <- event_with_cov %>%
  tidyr::fill(everything(), .direction = "down")
anyNA(event_with_cov)             
sum(is.na(event_with_cov))        






# save(event_with_cov,file="ant_data_longformat.Rdata")
# load("ant_data_longformat.Rdata")
event_with_cov <- event_with_cov %>% dplyr::select(-c("Zedge_NF"))   ######Set N-F type as reference 
event_with_cov <- event_with_cov %>% dplyr::select(-c("group_1","group_2"))


}


########## 2. Run the analysis -----------------------
###########======================
# check id 
{
  library(dplyr)
  thr <- 5000
  

  edge_nodes <- event_with_cov %>%
    group_by(id) %>%
    summarise(id1 = first(id1), id2 = first(id2), .groups = "drop")
  

  drop_ids <- edge_nodes %>%
    filter(id1 > thr & id2 > thr) %>%
    pull(id)
  

  event_sub <- event_with_cov %>%
    filter(!(id %in% drop_ids))
  

  id_map <- event_sub %>% distinct(id) %>% arrange(id) %>% mutate(id_new = row_number())
  event_sub <- event_sub %>% left_join(id_map, by = "id") %>%
    mutate(id = id_new) %>% dplyr::select(-id_new)
  

  c(before = n_distinct(event_with_cov$id), after = n_distinct(event_sub$id))
}

real_data <- list()
real_data$df <- event_sub

{
  library(dplyr)
  library(stringr)
  

  event_sub <- event_sub %>%
    mutate(id1 = as.integer(id1),
           id2 = as.integer(id2))
  

  node_levels <- sort(unique(c(event_sub$id1, event_sub$id2)))
  node_map <- tibble(old_node = node_levels,
                     node     = seq_along(node_levels))
  
  event_sub <- event_sub %>%

    left_join(node_map, by = c("id1" = "old_node")) %>%
    rename(id1_new = node) %>%

    left_join(node_map, by = c("id2" = "old_node")) %>%
    rename(id2_new = node) %>%

    mutate(id1 = id1_new, id2 = id2_new) %>%
    dplyr::select(-id1_new, -id2_new) %>%

    mutate(i = pmin(id1, id2),
           j = pmax(id1, id2),
           id1 = i, id2 = j) %>%
    dplyr::select(-i, -j)
  


  event_sub <- event_sub %>% rename(old_id = id)
  
  edge_map <- event_sub %>%
    distinct(old_id, id1, id2) %>%
    arrange(old_id) %>%
    mutate(id = row_number())
  
  event_sub <- event_sub %>%
    left_join(edge_map %>% dplyr::select(old_id, id), by = "old_id") %>%
    dplyr::select(-old_id) %>%
    arrange(id, start, stop)
  

  pairs_mat <- event_sub %>%
    distinct(id, id1, id2) %>%
    arrange(id) %>%
    dplyr::select(id1, id2) %>%
    as.matrix()
  storage.mode(pairs_mat) <- "integer"
  
  real_data$df    <- event_sub
  real_data$pairs <- pairs_mat
  

  N <- length(unique(c(event_sub$id1, event_sub$id2)))
  stopifnot(identical(sort(unique(c(event_sub$id1, event_sub$id2))), seq_len(N)))
  M <- nrow(real_data$pairs)
  stopifnot(identical(sort(unique(event_sub$id)), seq_len(M)))
  

  head(real_data$pairs)
  
}


df <- real_data$df
cols <- c("Zedge_NN","Zedge_NC","Zedge_CC","Zedge_CF","Zedge_FF")
cols <- c(cols, "Zedge_avg_age", "Zedge_avg_body_size")

library(dplyr)
real_data$df <- real_data$df %>%
  mutate(
    Zedge_avg_age = if_else(is.na(age_1) & is.na(age_2), 0,
                            rowMeans(cbind(age_1, age_2), na.rm = TRUE)),
    Zedge_avg_body_size = if_else(is.na(body_size_1) & is.na(body_size_2), 0,
                                  2*rowMeans(cbind(body_size_1, body_size_2), na.rm = TRUE)),
    across(any_of(cols), ~ . - mean(., na.rm = TRUE))
  )

library(Matrix)
library(survival)

source("real_data/cpp_run.R")
fit_real_c <- fit_edge_frailty_ppl_block_pcg_tv_fast_cpp(
  data_sim = real_data,           
  zi1_pattern = "_1$",    
  zi2_pattern = "_2$",     
  max_outer = 25,    
  beta_maxit =10,
  tol = 1e-3,
  damping = 0.9,
  ties_method = "efron"  
)

fit_real <- fit_real_c

#####check id
{

  df <- event_with_cov
  

  pairs_tbl <- df |> dplyr::distinct(id, id1, id2) |> dplyr::arrange(id)
  m_df <- nrow(pairs_tbl)
  cat("m_df =", m_df, "\n")
  

  Xe_cols <- grep("^Zedge_", names(df), value = TRUE)
  cat("Number of Zedge columns =", length(Xe_cols), "\n")
  

  if (exists("scan")) {
    m_scan <- if (!is.null(scan$m)) scan$m else NA_integer_
    K_scan <- if (!is.null(scan$K)) scan$K else NA_integer_
    cat("scan$m =", m_scan, " scan$K =", K_scan, "\n")
    if (!is.na(m_scan)) stopifnot(m_scan == m_df)
  } else {
    cat("scan is unavailable (or has not been updated)\n")
  }
  

  if (exists("ops") && !is.null(ops$S2)) {
    cat("dim(ops$S2) =", paste(dim(ops$S2), collapse="x"), "\n")
    stopifnot(all(dim(ops$S2) == c(m_df, m_df)))
  } else cat("ops$S2 is unavailable (or has not been updated)\n")
  

  if (exists("delta")) {
    len_delta <- length(delta)
    cat("length(delta) =", len_delta, "\n")
  }
  
}



df <-real_data$df
scan <- build_tv_scanline_index(real_data$df, zedge_pat="^Zedge_", 
                                zi1_pat="^body_size_", zi2_pat="^age_")
n    <- if (!is.null(real_data$Z_i)) nrow(real_data$Z_i) else max(real_data$pairs)
ops  <- make_linegraph_ops(n, real_data$pairs)
S2   <- make_S2(n, real_data$pairs)
ties_method  <- "efron"
kern  <- make_scan_kernel(scan, fit_real$beta_all,ties_method = ties_method)
lam <- compute_cumLambda(scan, fit_real$beta_all, fit_real$delta, method = ties_method)  # or "efron"

# Kdiag <- compute_K_diag_full(
#   scan, fit_real$beta_all, fit_real$delta,
#   cumLambda_by_time = lam$cumLambda_by_time
# )


Kdiag_real <- compute_K_diag_partial(scan,fit_real$beta_all,fit_real$delta,ties_method = "efron")
S2_real   <- make_S2(162, real_data$pairs)


ci_out_real  <- wald_ci_sigma2_rho(
  delta = fit_real$delta, n = n, S2 = S2_real,
  sigma2_hat =  fit_real$sigma^2, rho_hat = fit_real$rho,
  n_mc = 32L, n_subspace = 8L, 
  exact_trace = FALSE,
  K_diag = Kdiag_real,  chol_threshold = 5000
)



ci_out_real$ci_sigma2   
ci_out_real $ci_rho      



#### Regression parameter inference
{

  {
    Rcpp::sourceCpp(file="real_data/RCPP/infer_beta.cpp")
    Rcpp::sourceCpp(file="real_data/RCPP/schur_theta_block.cpp")
    
    
    

    build_theta_delta_ops_matmul <- function(scan, beta, delta, ties_method = c("breslow","efron")) {
      ties_method <- match.arg(ties_method)
      list(
        Htt_matrix = function() Htt_matrix_cpp(scan, beta, delta, ties_method),
        Htd_matmul = function(U)  Htd_matmul_cpp(scan, beta, delta, as.matrix(U), ties_method),
        Hdt_matmul = function(V)  Hdt_matmul_cpp(scan, beta, delta, as.matrix(V), ties_method)
      )
    }
    
    edge_joint_scores_fast <- function(scan, beta, delta, ties_method = c("breslow","efron")) {
      ties_method <- match.arg(ties_method)
      edge_joint_scores_cpp(scan, beta, delta, ties_method)
    }
    
    {
      cox_frailty_theta_vcov <- function(
    scan, pairs,
    beta_hat,
    delta_hat, sigma2, rho,
    K_diag ,
    ties_method = c("breslow","efron"),
    sand = FALSE,
    cluster = c("dyad","gray","edge","node"),
    beta_col_pattern  = "^Zedge_",
    gamma_col_pattern = "^Xi_plus_Xj_",
    ridge = 1e-8, chol_threshold = 5000L,
    pcg_tol = 1e-6, pcg_maxit = 2000
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
        if (!length(beta_idx))  stop("No beta columns matched: ", beta_col_pattern)
        if (!length(gamma_idx)) {
          gamma_idx <- setdiff(seq_len(p), beta_idx)
          # if (!length(gamma_idx)) stop("Gamma columns not found; set gamma_col_pattern.")
        }


        n_node <- max(pairs)
        ops    <- make_linegraph_ops(n_node, pairs)
        Msolve <- make_M_solver(ops, rho, chol_threshold, pcg_tol, pcg_maxit)
        OmegaInv_apply <- function(x) (1/sigma2) * Msolve(x)
        

        sk <- make_scan_kernel(scan, beta = theta_hat, ties_method = ties_method)
        K_apply <- function(x) sk$K_matvec(x, delta = delta_hat)
        

        Hdd_matvec <- function(x) K_apply(x) + OmegaInv_apply(x)
        
        Mdiag_pcg  <- K_diag + (1/sigma2)
        Hdd_solve  <- function(rhs) pcg_solve(Hdd_matvec, rhs, Mdiag = Mdiag_pcg,
                                              tol = pcg_tol, maxit = pcg_maxit)
        Hdd_solve_mat <- function(R){ R <- if (is.matrix(R)) R else cbind(R)
        out <- matrix(0.0, nrow(R), ncol(R))
        for (j in seq_len(ncol(R))) out[, j] <- Hdd_solve(R[, j])
        out
        }
        
        print(1)

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
        



        
        


          TDops <- build_theta_delta_ops_matmul(scan, beta_hat, delta_hat,
                                                ties_method = "efron")

          I_tt <- TDops$Htt_matrix() + diag(1e-8, ncol(scan$X_rows))

          return(I_tt)











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
            

            es <- edge_joint_scores(scan, theta_hat, delta_hat, ties_method)
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
        take <- function(M, r, c) M[r, c, drop = FALSE]
        Var_beta  <- take(Var_theta, beta_idx,  beta_idx)
        Var_gamma <- take(Var_theta, gamma_idx, gamma_idx)
        Cov_bg    <- take(Var_theta, beta_idx,  gamma_idx)
        
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
  
  
  inf_c <- cox_frailty_theta_vcov (
    scan, pairs = real_data$pairs,
    beta_hat =  fit_real$beta_all, 
    delta_hat = fit_real$delta,
    sand = FALSE, cluster = "node",
    K_diag = Kdiag_real,
    sigma2 = (fit_real$sigma)^2, rho = fit_real$rho,
    ties_method = "efron",          # or "efron"
    beta_col_pattern  = "^Zedge_",
    gamma_col_pattern = "^Xi_plus_Xj_"
  )
  beta_sd <- sqrt(diag(solve(inf_c)))
}

fit_real_c$beta_all


library(ggplot2)



########## Daily average number of events per ant pair
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  

  avg_by_day <- time_data %>%
    summarise(across(matches("^Day\\d+$"), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Day", values_to = "avg_per_pair") %>%
    mutate(day = as.integer(sub("Day", "", Day))) %>%
    arrange(day)
  
  ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
    geom_line() +
    geom_point() +
    labs(title = "Daily average number of events per ant pair (excluding NA)",
         x = "Day", y = "Average events per pair") +
    theme_minimal(base_size = 13)
  
}


########## Introduction figure: interaction counts on the left and time on the right
{

  {
    
    {
      library(ggplot2)
      library(dplyr)
      library(MASS)
      library(scales)
      
      x <- low_interactions$total_interaction
      x <- x[is.finite(x)]
      n <- length(x)
      
      mu <- mean(x)
      v  <- var(x)
      phi <- v / mu
      p_chi <- pchisq((n - 1) * v / mu, df = n - 1, lower.tail = FALSE)
      
      lab <- sprintf("Mean = %.4g\nVar = %.4g\nVar/Mean = %.3g\nPoisson λ̂ = %.4g\nChiSq p ≈ %.2g",
                     mu, v, phi, mu, p_chi)
      
      binwidth <- 10
      

      df <- data.frame(x = x)
      fit_nb <- MASS::glm.nb(x ~ 1, data = df)
      theta_hat <- fit_nb$theta
      mu_nb <- as.numeric(exp(coef(fit_nb))[1])
      

      x_cap <- ceiling(max(x) / binwidth) * binwidth
      bin_lefts  <- seq(0, x_cap - binwidth, by = binwidth)
      bin_rights <- bin_lefts + (binwidth - 1)
      bin_mids   <- bin_lefts + binwidth/2
      
      obs_df <- tibble(bin_left = floor(x / binwidth) * binwidth) %>%
        count(bin_left, name = "O") %>%
        right_join(tibble(bin_left = bin_lefts), by = "bin_left") %>%
        mutate(O = ifelse(is.na(O), 0, O)) %>%
        arrange(bin_left) %>%
        mutate(bin_mid = bin_left + binwidth/2)
      
      exp_df <- tibble(
        bin_left = bin_lefts,
        bin_mid  = bin_mids,
        E_pois = n * (ppois(bin_rights, mu) - ppois(bin_lefts - 1, mu)),
        E_nb   = n * (pnbinom(bin_rights, size = theta_hat, mu = mu_nb) -
                        pnbinom(bin_lefts - 1, size = theta_hat, mu = mu_nb))
      )
      
      plot_df <- left_join(obs_df, exp_df, by = c("bin_left","bin_mid")) %>%
        mutate(
          resid_pois = (O - E_pois) / sqrt(pmax(E_pois, 1e-8))
        )
      

      p1 <- ggplot(low_interactions, aes(x = total_interaction)) +
        geom_histogram(
          binwidth = binwidth, fill = "steelblue", color = "black",
          boundary = 0, closed = "left"
        ) +
        geom_line(data = exp_df, aes(x = bin_mid, y = E_pois),
                  inherit.aes = FALSE, linewidth = 1.1) +
        geom_line(data = exp_df, aes(x = bin_mid, y = E_nb),
                  inherit.aes = FALSE, linewidth = 1.1, linetype = 2) +
        annotate("text", x = Inf, y = Inf, label = lab,
                 hjust = 1.1, vjust = 1.1, size = 4.2) +
        labs(
          title = "Total Interactions: Observed vs Poisson/NB expected counts",
          subtitle = "Solid: Poisson; Dashed: Negative Binomial (captures overdispersion)",
          x = "Total Interactions",
          y = "Number of Edges (log1p scale)"
        ) +
        theme_minimal(base_size = 14) +
        scale_y_sqrt()
      scale_y_continuous(trans = "log1p")
      p1
      
      {
        library(ggplot2)
        library(dplyr)
        
        plot_total_interaction_like <- function(time_data,
                                                binwidth = 10,
                                                x_max = 900,
                                                y_max = 1050) {

          df <- time_data %>%
            mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
            filter(is.finite(total_interaction))
          

          mu <- mean(df$total_interaction)
          va <- var(df$total_interaction)
          
          lab <- sprintf("Mean = %.3f\n\nVariance = %.3f", mu, va)
          

          ggplot(df, aes(x = total_interaction)) +
            geom_histogram(
              binwidth = binwidth,
              boundary = 0, closed = "left",
              fill = "#BFD7EA", color = NA
            ) +
            annotate(
              "text", x = Inf, y = Inf, label = lab,
              hjust = 1.02, vjust = 1.05, size = 4.5
            ) +
            labs(
              title = "a)   Distribution of Total Interaction",
              x = "Interaction",
              y = "Count"
            ) +
            coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
            theme_classic(base_size = 14) +
            theme(

              panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
              panel.grid.minor = element_blank(),

              axis.line = element_line(color = "black", linewidth = 0.9),
              axis.ticks = element_line(color = "black"),
              axis.text = element_text(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.05)
            )
        }
        

        p <- plot_total_interaction_like(time_data, binwidth = 10, x_max = 900, y_max = 1050)
        print(p)
        
      }
      
      
    }
    
    {
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      

      avg_by_day <- time_data %>%
        summarise(across(matches("^Day\\d+$"), ~ mean(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Day", values_to = "avg_per_pair") %>%
        mutate(day = as.integer(sub("Day", "", Day))) %>%
        arrange(day)
      

      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      s <- summary(fit)
      
      slope <- coef(s)[ "day", "Estimate" ]
      pval  <- coef(s)[ "day", "Pr(>|t|)" ]
      
      lab <- sprintf("Slope = %.4g\np = %.4g", slope, pval)
      

      ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
        geom_line() +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.1, size = 4.5) +
        labs(title = "Daily average number of events per ant pair (excluding NA)",
             x = "Day", y = "Average events per pair") +
        theme_minimal(base_size = 13)
      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      

      daily_summary_df <- avg_by_day %>%
        transmute(
          day = day,
          avg_per_pair = avg_per_pair,
          fitted_avg_per_pair = as.numeric(predict(fit, newdata = avg_by_day))
        )
      
      daily_summary_df
      
    }
    library(writexl)
    
    write_xlsx(
      list(
        low_interactions = data.frame(low_interactions$total_interaction),
        avg_by_day       = avg_by_day
      ),
      path = "intro_data_plot.xlsx"
    )
    library(ggplot2)
    library(patchwork)
    

    p1 <- ggplot(low_interactions, aes(x = total_interaction)) +
      geom_histogram(binwidth = 10, fill = "steelblue", color = "black",
                     boundary = 0, closed = "left") +
      labs(title = "Distribution of Total Interactions",
           x = "Total Interactions", y = "Number of Edges") +
      theme_minimal(base_size = 14)
    

    p2 <- ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
      geom_line() + geom_point() +
      labs(title = "Average daily number of interactions per ant pair
 ",
           x = "Day", y = "Number of interactons") +
      theme_minimal(base_size = 13)
    

    (p1 | p2) +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(tag_levels = "A")
    
  }
  
  
  

  {
    corner_anno <- function(label, size = 4.6, h = 1.02, v = 1.05) {
      annotate("text", x = Inf, y = Inf, label = label,
               hjust = h, vjust = v, size = size)
    }
    
    
    plot_ant_figure_like_nobreak <- function(time_data,
                                             binwidth = 10,
                                             x_max_hist = 900,
                                             y_max_hist = 1050,
                                             y_pad_ratio = 0.06,
                                             widths = c(1.15, 1.00)) {
      
      suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(patchwork)
      })
      



      df_a <- time_data %>%
        mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
        filter(is.finite(total_interaction))
      
      mu <- mean(df_a$total_interaction)
      va <- var(df_a$total_interaction)
      lab_a <- sprintf("Mean = %.3f\n\nVariance = %.3f", mu, va)
      
      p1 <- ggplot(df_a, aes(x = total_interaction)) +
        geom_histogram(binwidth = binwidth, boundary = 0, closed = "left",
                       fill = "#BFD7EA", color = NA) +
        corner_anno(lab_a) +
        labs(title = "a)   Distribution of Total Interaction",
             x = "Interaction", y = "Count") +
        coord_cartesian(xlim = c(0, x_max_hist), ylim = c(0, y_max_hist), expand = FALSE) +
        theme_classic(base_size = 14) +
        theme(
          panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.9),
          axis.ticks = element_line(color = "black"),
          axis.text  = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.05),
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6)
        )
      
      



      mat <- as.matrix(time_data[, -1, drop = FALSE])
      avg_by_day <- tibble(
        day = seq_len(ncol(mat)),
        avg_per_pair = colMeans(mat, na.rm = TRUE)
      ) %>% filter(is.finite(avg_per_pair))
      
      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      slope <- unname(coef(fit)["day"])
      pval  <- summary(fit)$coefficients["day", "Pr(>|t|)"]
      
      p_text <- if (is.finite(pval) && pval < 0.001) "P value < 0.001" else sprintf("P value = %.3f", pval)
      lab_b  <- sprintf("Slope = %.3f\n\n%s", slope, p_text)
      
      blue_pt <- "#1F77B4"
      

      y_min <- min(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_max <- max(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_pad <- (y_max - y_min) * y_pad_ratio
      if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2
      

      x_anno <- max(avg_by_day$day) - 1
      y_anno <- y_max + y_pad * 0.2
      
      p2 <- ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
        geom_line(color = blue_pt, linewidth = 0.8) +
        geom_point(color = blue_pt, size = 2.6) +
        geom_smooth(method = "lm", se = FALSE,
                    color = "black", linetype = "dotted", linewidth = 0.9) +
        corner_anno(lab_b) +
        labs(title = "b)   Average Daily Number of Interactions per Ant Pair",
             x = "Day", y = "Number of Interactions") +
        scale_x_continuous(breaks = c(0, 10, 20, 30, 40),
                           limits = c(0, max(avg_by_day$day)),
                           expand = expansion(mult = c(0, 0.02))) +
        coord_cartesian(ylim = c(y_min - 0.05 * (y_max - y_min), y_max + y_pad), expand = FALSE) +
        theme_classic(base_size = 14) +
        theme(
          panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.9),
          axis.ticks = element_line(color = "black"),
          axis.text  = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.05),
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6)
        )
      



      fig <- p1 + p2 + plot_layout(ncol = 2, widths = widths)
      return(fig)
    }
    

    fig <- plot_ant_figure_like_nobreak(time_data)
    print(fig)

    
  }
  
  {
    plot_ant_figure_like_nobreak <- function(time_data,
                                             binwidth = 10,
                                             x_max_hist = 900,
                                             y_max_hist = NULL,
                                             y_pad_ratio = 0.06,
                                             font_family = "Times New Roman",
                                             font_scale = 2,
                                             show_curve_legend = FALSE,

                                             anno_x_frac_left  = 0.68,
                                             anno_x_frac_right = 0.78,
                                             anno_y_frac       = 0.995) {
      
      suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(patchwork)
      })
      
      base_size <- 14 * font_scale
      anno_size <- 4.6 * font_scale
      

      corner_anno_left <- function(label, x_pos, y_pos) {
        annotate("text",
                 x = x_pos, y = y_pos, label = label,
                 hjust = 0, vjust = 1,
                 size = anno_size, family = font_family)
      }
      



      df_a <- time_data %>%
        mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
        filter(is.finite(total_interaction))
      
      x_all <- as.numeric(df_a$total_interaction)
      x_all <- x_all[is.finite(x_all) & x_all >= 0]
      n <- length(x_all)
      
      mu <- mean(x_all)
      va <- var(x_all)
      lab_a <- sprintf("Mean = %.3f\n\nVariance = %.3f", mu, va)
      
      lambda_hat <- max(mu, 0)
      

      if (is.finite(va) && va > mu && mu > 0) {
        size_hat <- mu^2 / (va - mu)
      } else {
        size_hat <- 1e8
      }
      mu_nb_hat <- max(mu, 0)
      

      breaks_disp <- seq(0, x_max_hist + binwidth, by = binwidth)
      lefts  <- breaks_disp[-length(breaks_disp)]
      rights <- breaks_disp[-1]
      mids   <- (lefts + rights) / 2
      
      pbin_pois <- ppois(rights - 1, lambda_hat) - ppois(lefts - 1, lambda_hat)
      pbin_nb   <- pnbinom(rights - 1, size = size_hat, mu = mu_nb_hat) -
        pnbinom(lefts - 1, size = size_hat, mu = mu_nb_hat)
      
      curve_df <- rbind(
        data.frame(mid = mids, expected = n * pbin_pois, dist = "Poisson"),
        data.frame(mid = mids, expected = n * pbin_nb,   dist = "NegBin")
      )
      

      if (is.null(y_max_hist)) {
        x_clip <- x_all[x_all <= x_max_hist]
        h_counts <- hist(x_clip, breaks = breaks_disp, plot = FALSE, right = FALSE)$counts
        y_max_hist <- max(c(h_counts, curve_df$expected), na.rm = TRUE) * 1.06
      }
      

      x_anno_a <- x_max_hist * anno_x_frac_left
      y_anno_a <- y_max_hist * anno_y_frac
      
      p1 <- ggplot(df_a, aes(x = total_interaction)) +
        geom_histogram(
          binwidth = binwidth,
          boundary = 0, closed = "left",
          fill = "#BFD7EA", color = NA
        ) +
        geom_line(
          data = curve_df,
          aes(x = mid, y = expected, linetype = dist),
          inherit.aes = FALSE,
          color = "black", linewidth = 0.9
        ) +
        scale_linetype_manual(values = c(Poisson = "dashed", NegBin = "solid"), name = NULL) +
        corner_anno_left(lab_a, x_anno_a, y_anno_a) +
        labs(
          title = "a)   Distribution of Total Interaction",
          x = "Interaction",
          y = "Count"
        ) +
        coord_cartesian(xlim = c(0, x_max_hist), ylim = c(0, y_max_hist), expand = FALSE) +
        scale_y_sqrt() +
        theme_classic(base_size = base_size, base_family = font_family) +
        theme(
          panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.9),
          axis.ticks = element_line(color = "black"),
          axis.text  = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.05),
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6),
          legend.position = if (show_curve_legend) c(0.72, 0.35) else "none",
          legend.background = element_blank()
        )
      



      mat <- as.matrix(time_data[, -1, drop = FALSE])
      avg_by_day <- tibble(
        day = seq_len(ncol(mat)),
        avg_per_pair = colMeans(mat, na.rm = TRUE)
      ) %>% filter(is.finite(avg_per_pair))
      
      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      slope <- unname(coef(fit)["day"])
      pval  <- summary(fit)$coefficients["day", "Pr(>|t|)"]
      p_text <- if (is.finite(pval) && pval < 0.001) "P value < 0.001" else sprintf("P value = %.3f", pval)
      lab_b  <- sprintf("Slope = %.3f\n\n%s", slope, p_text)
      
      y_min <- min(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_max <- max(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_pad <- (y_max - y_min) * y_pad_ratio
      if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2
      

      x_max_b <- max(avg_by_day$day)
      y_upper_b <- y_max + y_pad
      y_lower_b <- y_min - 0.05 * (y_max - y_min)
      
      x_anno_b <- x_max_b * anno_x_frac_right
      y_anno_b <- y_upper_b * anno_y_frac
      
      blue_pt <- "#1F77B4"
      
      p2 <- ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
        geom_line(color = blue_pt, linewidth = 0.8) +
        geom_point(color = blue_pt, size = 2.6) +
        geom_smooth(method = "lm", se = FALSE,
                    color = "black", linetype = "dotted", linewidth = 0.9) +
        corner_anno_left(lab_b, x_anno_b, y_anno_b) +
        labs(
          title = "b)   Average Daily Number of Interactions per Ant Pair",
          x = "Day",
          y = "Number of Interactions"
        ) +
        scale_x_continuous(
          breaks = c(0, 10, 20, 30, 40),
          limits = c(0, x_max_b+1),
          expand = expansion(mult = c(0, 0.02))
        ) +
        coord_cartesian(ylim = c(y_lower_b, y_upper_b), expand = FALSE) +
        theme_classic(base_size = base_size, base_family = font_family) +
        theme(
          panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.9),
          axis.ticks = element_line(color = "black"),
          axis.text  = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.05),
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6)
        )
      



      p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
    }
    

    fig <- plot_ant_figure_like_nobreak(time_data, font_scale = 1.8,anno_y_frac       = 1)
    print(fig)
    
  }

  
  ggsave(
    filename = "Rplot.pdf",
    plot = fig,
    device = cairo_pdf,
    width = 24, height = 11, units = "in",
    limitsize = FALSE
  )
  {
    plot_ant_left_only <- function(time_data,
                                   binwidth = 10,
                                   x_max_hist = 900,
                                   y_max_hist = NULL,
                                   font_family = "Times New Roman",
                                   font_scale = 2) {
      
      suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
      })
      
      base_size <- 14 * font_scale
      

      df_a <- time_data %>%
        mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
        filter(is.finite(total_interaction))
      
      x_all <- as.numeric(df_a$total_interaction)
      x_all <- x_all[is.finite(x_all) & x_all >= 0]
      n <- length(x_all)
      

      mu <- mean(x_all)
      va <- var(x_all)
      
      lambda_hat <- max(mu, 0)
      

      if (is.finite(va) && va > mu && mu > 0) {
        size_hat <- mu^2 / (va - mu)
      } else {
        size_hat <- 1e8
      }
      mu_nb_hat <- max(mu, 0)
      

      breaks_disp <- seq(0, x_max_hist + binwidth, by = binwidth)
      lefts  <- breaks_disp[-length(breaks_disp)]
      rights <- breaks_disp[-1]
      mids   <- (lefts + rights) / 2
      
      pbin_pois <- ppois(rights - 1, lambda_hat) - ppois(lefts - 1, lambda_hat)
      pbin_nb   <- pnbinom(rights - 1, size = size_hat, mu = mu_nb_hat) -
        pnbinom(lefts - 1, size = size_hat, mu = mu_nb_hat)
      
      curve_df <- rbind(
        data.frame(mid = mids, expected = n * pbin_pois, dist = "Poisson"),
        data.frame(mid = mids, expected = n * pbin_nb,   dist = "NegBin")
      )
      

      if (is.null(y_max_hist)) {
        x_clip <- x_all[x_all <= x_max_hist]
        h_counts <- hist(x_clip, breaks = breaks_disp, plot = FALSE, right = FALSE)$counts
        y_max_hist <- max(c(h_counts, curve_df$expected), na.rm = TRUE) * 1.06
      }
      
      p_left <- ggplot(df_a, aes(x = total_interaction)) +
        geom_histogram(
          binwidth = binwidth,
          boundary = 0, closed = "left",
          fill = "#BFD7EA", color = NA
        ) +
        geom_line(
          data = subset(curve_df, dist == "Poisson"),
          aes(x = mid, y = expected),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black", linewidth = 0.9
        ) +
        geom_line(
          data = subset(curve_df, dist == "NegBin"),
          aes(x = mid, y = expected),
          inherit.aes = FALSE,
          linetype = "solid",
          color = "black", linewidth = 0.9
        ) +
        labs(x = "Interaction", y = "Count") +
        coord_cartesian(xlim = c(0, x_max_hist), ylim = c(0, y_max_hist), expand = FALSE) +
        scale_y_sqrt() +
        theme_classic(base_size = base_size, base_family = font_family) +
        theme(
          panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.9),
          axis.ticks = element_line(color = "black"),
          axis.text  = element_text(color = "black"),
          plot.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6)
        )
      
      return(p_left)
    }
    

    p_left <- plot_ant_left_only(time_data, binwidth = 10, x_max_hist = 900, font_scale = 2)
    print(p_left)
    ggsave(
      filename = "p_left.pdf",
      plot = p_left,
      device = cairo_pdf,
      width = 11, height = 10, units = "in",
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    ggsave(
      filename = "p_left.tiff",
      plot = p_left,
      device = ragg::agg_tiff,
      width = 11, height = 10, units = "in",
      dpi = 1000,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
    
  }
  ggsave("ant_fig.pdf", fig, width = 24, height = 11, units = "in", device = cairo_pdf)

  ggsave("ant_fig.png", fig, width = 14, height = 5.2, units = "in", dpi = 300)
  
}


#### Select 5 points and plot the distribution of delta_{ij}
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  



  delta <- as.numeric(fit_real$delta)
  pairs <- as.matrix(real_data$pairs)
  

  n_node <- max(pairs)
  




  Dmat <- matrix(NA_real_, n_node, n_node)
  
  for (e in seq_along(delta)) {
    i <- pairs[e, 1]
    j <- pairs[e, 2]
    Dmat[i, j] <- delta[e]
    Dmat[j, i] <- delta[e]
  }
  



  node_mean <- apply(Dmat, 1, function(x) mean(x, na.rm = TRUE))
  





  ord <- order(node_mean)
  k   <- length(ord)
  
  probs   <- c(0, 0.25, 0.5, 0.75, 1)
  pos     <- 1 + probs * (k - 1)
  idx_pos <- round(pos)
  idx_pos <- pmin(pmax(idx_pos, 1), k)
  
  id_very_low  <- ord[idx_pos[1]]
  id_low       <- ord[idx_pos[2]]
  id_mid       <- ord[idx_pos[3]]
  id_high      <- ord[idx_pos[4]]
  id_very_high <- ord[idx_pos[5]]
  
  c(id_very_low  = id_very_low,
    id_low       = id_low,
    id_mid       = id_mid,
    id_high      = id_high,
    id_very_high = id_very_high)

  




  vec_very_high <- Dmat[id_very_high, -id_very_high]
  vec_high      <- Dmat[id_high,      -id_high]
  vec_mid       <- Dmat[id_mid,       -id_mid]
  vec_low       <- Dmat[id_low,       -id_low]
  vec_very_low  <- Dmat[id_very_low,  -id_very_low]
  






  









  df_long <- tibble(
    very_high = vec_very_high,
    high      = vec_high,
    medium    = vec_mid,
    low       = vec_low,
    very_low  = vec_very_low
  ) %>%
    mutate(idx = row_number()) %>%
    pivot_longer(
      cols      = c(very_high, high, medium, low, very_low),
      names_to  = "group",
      values_to = "delta"
    ) %>%
    filter(!is.na(delta)) %>%
    mutate(group = recode(
      group,
      very_high = "Very high frailty",
      high      = "High frailty",
      medium    = "Medium frailty",
      low       = "Low frailty",
      very_low  = "Very low frailty"
    ))
  



  cols <- c(
    "Very high frailty" = "#B2182B",
    "High frailty"      = "#D55E00",
    "Medium frailty"    = "#CC79A7",
    "Low frailty"       = "#0072B2",
    "Very low frailty"  = "#009E73"
  )
  




  df_long$group <- factor(df_long$group, levels = rev(names(cols)))
  
  p_density <- ggplot(df_long, aes(x = delta, color = group, fill = group)) +
    geom_density(alpha = 0.25, adjust = 1) +
    labs(
      x = expression(delta[ij]),
      y = "Density",
      color = NULL,
      fill  = NULL
    ) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_bw(base_size = 18) +
    theme(
      text           = element_text(size = 18),
      legend.position = "top",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 16),
      axis.title      = element_text(size = 18),
      axis.text       = element_text(size = 16),
      panel.grid      = element_blank()
    )
  
  print(p_density)



  p_box <- ggplot(df_long, aes(x = group, y = delta, fill = group)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.alpha = 0.4) +
    scale_fill_manual(values = cols) +
    labs(
      x = NULL,
      y = expression(delta[ij])
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.text.x     = element_text(angle = 20, hjust = 1)
    )
  
  print(p_box)
  
  
}

df_D <- data.frame(
  very_high = as.numeric(vec_very_high),
  high      = as.numeric(vec_high),
  mid       = as.numeric(vec_mid),
  low       = as.numeric(vec_low),
  very_low  = as.numeric(vec_very_low)
)


