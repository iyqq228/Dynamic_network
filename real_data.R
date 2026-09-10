library(readxl)
{
  calculate_avg <- function(value1, value2) {
    # 提取实际值
    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    
    if (is.na(value1) & is.na(value2)) {
      return(mean(cov_data$body_size, na.rm = TRUE)) # 如果两项都缺失，返回总体平均值
    } else if (is.na(value1)) {
      return(value2) # 如果value1缺失，返回value2
    } else if (is.na(value2)) {
      return(value1) # 如果value2缺失，返回value1
    } else {
      return(mean(c(value1, value2))) # 如果都不缺失，返回两项的平均值
    }
  }
  
  
  
  time_data <- read_excel("C:/Users/wqy/Desktop/课题/network via survival/ant-data-forcode/time1.xlsx")
  cov_data <- read_excel("C:/Users/wqy/Desktop/课题/network via survival/ant-data-forcode/cov.xlsx")
  
  # 解析交互类型（N-N, N-C, N-F, C-C, C-F, F-F）
  interaction_types <- c("Q", "N-N", "N-C", "N-F", "C-C", "C-F", "F-F")
  
  # 生成41个时变协变量数据框
  covariate_dfs <- list()
  
  # 计算体型和年龄的平均值时处理缺失数据的函数
  calculate_avg <- function(value1, value2) {
    if (is.na(value1) & is.na(value2)) {
      return(mean(cov_data$body_size, na.rm = TRUE)) # 如果两项都缺失，返回总体平均值
    } else if (is.na(value1)) {
      return(value2) # 如果value1缺失，返回value2
    } else if (is.na(value2)) {
      return(value1) # 如果value2缺失，返回value1
    } else {
      return(mean(c(value1, value2))) # 如果都不缺失，返回两项的平均值
    }
  }
  
  # 迭代time_data中的每一天（每一列）
  for (day in c(1,15,25,40)) {
    # 初始化一个数据框
    df <- data.frame(matrix(NA, nrow = nrow(time_data), ncol = 8))
    
    # 确定该天属于哪个时间段（group_period1到group_period4）
    if (day <= 11) {
      period_column <- "group_period1"
    } else if (day <= 21) {
      period_column <- "group_period2"
    } else if (day <= 31) {
      period_column <- "group_period3"
    } else {
      period_column <- "group_period4"
    }
    
    # 生成前6列，代表交互类型的0-1变量
    Ant_pair <- time_data$`Ant Pair`
    for (i in 1:nrow(time_data)) {
      pair <- unlist(strsplit(as.character(Ant_pair[i]), "-"))
      ant1 <- gsub("Ant", "", pair[1])  # 去掉前面的"Ant"
      ant2 <- gsub("Ant", "", pair[2])  # 去掉前面的"Ant"
      
      # 获取交互的蚂蚁种类
      group_ant1 <- cov_data[which(cov_data$tag_id == ant1), period_column]
      group_ant2 <- cov_data[which(cov_data$tag_id == ant2), period_column]
      
      # 为每种交互类型分配0或1
      if ("Q" %in% c(group_ant1, group_ant2)) {
        interaction_column <- which(interaction_types == "Q")
      } else {
        # 为每种交互类型分配0或1
        # 定义顺序
        order <- c("N", "C", "F")
        
        # 将 group_ant1 和 group_ant2 转换为其在顺序中的位置
        interaction <- ifelse(match(group_ant1, order) < match(group_ant2, order), 
                              paste(group_ant1, group_ant2, sep = "-"), 
                              paste(group_ant2, group_ant1, sep = "-"))
        
        interaction_column <- which(interaction_types == interaction)
      }
      df[i, interaction_column] <- 1
      
      # 提取体型和年龄数据，并确保它们是数字
      body_size_1 <- cov_data[which(cov_data$tag_id == ant1), "body_size"][[1]]
      body_size_2 <- cov_data[which(cov_data$tag_id == ant2), "body_size"][[1]]
      
      age_1 <- cov_data[which(cov_data$tag_id == ant1), "age"][[1]]
      age_2 <- cov_data[which(cov_data$tag_id == ant2), "age"][[1]]
      
      # 计算体型和年龄的平均值
      body_size_avg <- calculate_avg(body_size_1, body_size_2)
      age_avg <- calculate_avg(age_1, age_2)
      
      # 填入体型和年龄的平均值
      df[i, 8] <- body_size_avg
      df[i, 9] <- age_avg
    }
    
    # 处理前六列中的NA
    for (i in 1:nrow(df)) {
      if (any(!is.na(df[i, 1:7]))) { # 如果前六列中有非NA项
        df[i, 1:7][is.na(df[i, 1:7])] <- 0 # 将NA替换为0
      }
    }
    
    # 添加"Ant Pair"列
    df$`Ant Pair` <- Ant_pair
    
    # 保存每一天的数据框
    colnames(df) <- c("Q","N-N", "N-C", "N-F", "C-C", "C-F", "F-F", "avg_body_size", "avg_age", "Ant Pair")
    covariate_dfs[[day]] <- df
  }
  
  
  
  # 结果是41个数据框，保存在covariate_dfs列表中
}


time_data <- read_excel("C:/Users/wqy/Desktop/课题/network via survival/ant-data-forcode/time1.xlsx")
time_data <- read_excel("C:/Users/23203/Desktop/共享课题/network via survival/ant-data-forcode/time1.xlsx")


##从这里开始
{
  
  ############1.处理数据-------------------
  # 创建矩阵，排除第一列
  binary_matrix <- ifelse(!is.na(time_data[, -1]), 1, 0)
  
  ss1_common <- covariate_dfs[[1]]
  ss2_common <- covariate_dfs[[15]]
  ss3_common <- covariate_dfs[[25]]
  ss4_common <- covariate_dfs[[40]]
  normalize_z_score <- function(df) {
    # 对 avg_body_size 列进行 Z-score 标准化
    df$avg_body_size <- (df$avg_body_size - mean(df$avg_body_size, na.rm = TRUE)) /
      sd(df$avg_body_size, na.rm = TRUE)
    
    # 对 avg_age 列进行 Z-score 标准化
    df$avg_age <- (df$avg_age - mean(df$avg_age, na.rm = TRUE)) /
      sd(df$avg_age, na.rm = TRUE)
    
    return(df)
  }
  
  # 对 ss1_common 到 ss4_common 进行标准化
  ss1_common <- normalize_z_score(ss1_common)
  ss2_common <- normalize_z_score(ss2_common)
  ss3_common <- normalize_z_score(ss3_common)
  ss4_common <- normalize_z_score(ss4_common)
}

#####为了保持数值稳定，去掉一些行
###记得把data和协变量一起去掉
{
  num_data <- time_data[ , -1]
  
  
  ##1.按照删失率去除
  {
    # # 计算每行缺失率和总和
    # na_prop <- rowMeans(is.na(num_data))
    # row_sums <- rowSums(num_data, na.rm = TRUE)
    # 
    # # 找出要保留的行索引
    # keep_idx <- which(na_prop <= 0.6 & row_sums >= 30 & row_sums <= 500)
    # #keep_idx <- which(na_prop <= 0)
    # 
    # cc_idx <- which(ss4_common[ , "C-C"] == 1)
    # #cc_idx <- which(ss4_common[ , "F-F"] == 1 & row_sums >300)
    # # 最终保留的索引
    # keep_idx <- setdiff(keep_idx, cc_idx)
  }
  
  
  ###2.去掉ant 380
  {
    # idx_380 <- which(str_detect(time_data$`Ant Pair`, "380"))
    # all_idx <- seq_len(nrow(time_data))
    # 
    # keep_idx <- setdiff(all_idx, idx_380)
  }
  
  ###3.只保留N-N
  {
    #keep_idx <-  which(ss1_common[ , "N-F"] == 1)
    keep_idx <- c(1:nrow(time_data))
  }
  
  # kk <- time_data[161,]
  # kk <- ss4_common[161,]
  # 用于筛选后的数据
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

##把数据转化成长数据
############=============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ========= ① 宽表 -> 逐事件长表（每对必有一条 death） =========
# df_wide 至少包含："Ant Pair", Day1..DayK
wide_to_event_long_with_death <- function(df_wide,
                                          pair_col   = "Ant Pair",
                                          day_prefix = "Day") {
  stopifnot(pair_col %in% names(df_wide))
  
  # 找 Day 列并排序
  day_cols <- grep(paste0("^", day_prefix, "\\d+$"), names(df_wide), value = TRUE)
  if (!length(day_cols)) stop("没有找到形如 'Day1, Day2, ...' 的列。")
  day_idx  <- readr::parse_number(day_cols)
  ord      <- order(day_idx)
  day_cols <- day_cols[ord]
  day_idx  <- day_idx[ord]
  n_days   <- max(day_idx)
  
  # 第一处 NA 的天；若无 NA，设为 n_days + 1
  get_first_na_day <- function(x) {
    pos <- which(is.na(x))
    if (length(pos)) pos[1] else (n_days + 1L)
  }
  death_day_vec <- apply(df_wide[, day_cols, drop = FALSE], 1, get_first_na_day)
  
  # 宽 -> 长（含 NA）
  df_long <- df_wide %>%
    mutate(.row_id = row_number(),
           .death_day = death_day_vec) %>%
    pivot_longer(all_of(day_cols),
                 names_to = "day_label", values_to = "count") %>%
    mutate(day = readr::parse_number(day_label))
  
  # death 之前的观测（严格小于 death_day），剔除 NA 与 0
  df_interact_days <- df_long %>%
    filter(day < .death_day, !is.na(count), count > 0)
  
  # 拆分蚂蚁对
  df_interact_days <- df_interact_days %>%
    separate({{pair_col}}, into = c("ant1","ant2"), sep = "-", remove = FALSE)
  
  # 展开为逐事件，等距放置到日内 (day-1, day]
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
    select(`Ant Pair`, ant1, ant2, day, event_in_day,
           event_time, event_type, event, death, .row_id)
  
  # 每对的唯一 death 行：若某天 NA，则 time = (NA 的那天)-1；否则 time = n_days
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
    select(-.row_id)
}


# 输入：包含以下列的 event_long
#   "Ant Pair", "event_time", "event_type"（取值 "interaction" 或 "death"）
# 输出：在原有列基础上，新增 id/start/stop/event 四列；行数不变
augment_event_long_with_ag_columns <- function(event_long) {
  req <- c("Ant Pair", "event_time", "event_type")
  stopifnot(all(req %in% names(event_long)))
  
  # 为无向边构造稳定 id
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
    select(pair_key, id)
  
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
    select(-a1, -a2, -lo, -hi, -pair_key) %>%
    relocate(id, start, stop, event, .after = last_col())  # 也可放在最前面：.before = 1
  
  out
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# 在每条边(id)内，按 start/stop 排序后，用上一条事件的协变量填充 NA
# 默认会填充：Zedge_*、body_size_*、age_*、group_*、g1/g2/pair_cat/period
# 你也可以通过 extra_fill_cols 手动再加需要填的列名
locf_covariates_within_edge <- function(event_long, extra_fill_cols = character()) {
  stopifnot(all(c("id","start","stop") %in% names(event_long)))
  
  # 不填充的“结构列”（其余列默认可以填充）
  non_fill_cols <- intersect(c(
    "Ant Pair","ant1","ant2",
    "id","day","event_time","event_type",
    "event_in_day","n_in_day",
    "start","stop","event","death"
  ), names(event_long))
  
  # 建议填充的模式
  default_fill_patterns <- c("^Zedge_", "^body_size_", "^age_", "^group_",
                             "^g1$", "^g2$", "^pair_cat$", "^period$")
  
  # 根据模式抓取可填充列
  default_fill_cols <- unique(unlist(lapply(default_fill_patterns, function(pat) {
    grep(pat, names(event_long), value = TRUE)
  })))
  
  # 真正要填充的列 = 默认 + 额外指定 − 非填充结构列
  fill_cols <- setdiff(unique(c(default_fill_cols, extra_fill_cols)), non_fill_cols)
  fill_cols <- intersect(fill_cols, names(event_long))  # 只保留存在的列
  
  # 若没有可填列，直接返回
  if (length(fill_cols) == 0) {
    message("没有匹配到可填充的协变量列；原数据返回。")
    out <- event_long %>% arrange(id, start, stop)
    return(list(data = out, na_cols = names(out)[sapply(out, function(x) any(is.na(x)))]))
  }
  
  # 按边内时间做 LOCF
  out <- event_long %>%
    arrange(id, start, stop) %>%
    group_by(id) %>%
    tidyr::fill(all_of(fill_cols), .direction = "down") %>%
    ungroup()
  
  # 检查是否仍有 NA
  na_cols <- names(out)[sapply(out, function(x) any(is.na(x)))]
  
  # 小贴士：如果还存在 NA，通常是该边首行本来就是 NA（没有“上一条”可填）
  # 你可以改用 downup：.direction="downup" 让首行再用“向上”补一次（即用将来一条）
  # out2 <- event_long %>% arrange(id,start,stop) %>% group_by(id) %>%
  #         tidyr::fill(all_of(fill_cols), .direction="downup") %>% ungroup()
  
  list(data = out, na_cols = na_cols)
}


# ========= ② 并入节点/边协变量 =========
# nodes_df 需包含：
#   tag_id, body_size, age, group_period1, group_period2, group_period3, group_period4
add_covariates_to_events <- function(event_long, nodes_df) {
  
  library(dplyr); library(tidyr); library(stringr); library(readr)
  
  # 统一 tag_id 为整数
  nodes_df <- nodes_df %>% mutate(tag_id = suppressWarnings(as.integer(tag_id)))
  
  # day -> period：1–11->1；12–21->2；22–31->3；32–41->4
  day_to_period <- function(d) {
    d <- pmax(pmin(as.integer(d), 41L), 1L)
    cut(d, breaks = c(0,11,21,31,41),
        labels = c(1L,2L,3L,4L), right = TRUE, include.lowest = TRUE) |> as.integer()
  }
  
  # 节点分组：宽 -> 长
  nodes_groups_long <- nodes_df %>%
    rename(group1 = group_period1, group2 = group_period2,
           group3 = group_period3, group4 = group_period4) %>%
    pivot_longer(starts_with("group"), names_to = "period_lab", values_to = "group") %>%
    mutate(period = readr::parse_number(period_lab)) %>%
    select(tag_id, body_size, age, period, group)
  
  nodes1 <- nodes_groups_long %>%
    rename(tag_id_1 = tag_id, body_size_1 = body_size, age_1 = age, group_1 = group)
  nodes2 <- nodes_groups_long %>%
    rename(tag_id_2 = tag_id, body_size_2 = body_size, age_2 = age, group_2 = group)
  
  # 从 "Ant57" 提取数值 ID，并映射 period
  ev <- event_long %>%
    mutate(
      id1 = suppressWarnings(as.integer(stringr::str_extract(ant1, "\\d+"))),
      id2 = suppressWarnings(as.integer(stringr::str_extract(ant2, "\\d+"))),
      period = day_to_period(day)
    ) %>%
    left_join(nodes1, by = c("id1" = "tag_id_1", "period" = "period")) %>%
    left_join(nodes2, by = c("id2" = "tag_id_2", "period" = "period"))
  
  # 标准化分组标签
  normalize_group <- function(x) {
    x <- toupper(as.character(x))
    x[!x %in% c("N","C","F")] <- NA_character_
    x
  }
  ev <- ev %>%
    mutate(g1 = normalize_group(group_1),
           g2 = normalize_group(group_2))
  
  # ✅ 无向归一化：使用自定义顺序 N < C < F（而非字母序）
  rank_map <- c(N = 1L, C = 2L, F = 3L)
  undirected_pair <- function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_character_)
    if (a == b) return(paste0(a, "-", b))
    if (rank_map[[a]] < rank_map[[b]]) paste0(a, "-", b) else paste0(b, "-", a)
  }
  ev$pair_cat <- mapply(undirected_pair, ev$g1, ev$g2, USE.NAMES = FALSE)
  ev$edge_group_na <- as.integer(is.na(ev$pair_cat))
  
  # 6 个边虚拟量（保持你的命名）
  cats <- c("N-N","N-C","N-F","C-C","C-F","F-F")
  for (k in cats) {
    ev[[paste0("Zedge_", gsub("-", "", k))]] <- as.integer(ev$pair_cat == k)
  }
  
  # 整理列顺序
  ev %>%
    relocate(period, .after = day) %>%
    relocate(g1, g2, pair_cat, edge_group_na, .after = period) %>%
    relocate(body_size_1, age_1, body_size_2, age_2, .after = pair_cat)
}



##########这段是否需要？（可能不需要了）
rank_map <- c(N=1L, C=2L, F=3L)
undirected_pair <- function(a,b){
  if (is.na(a)||is.na(b)) return(NA_character_)
  if (a==b) return(paste0(a,"-",b))
  if (rank_map[[a]] < rank_map[[b]]) paste0(a,"-",b) else paste0(b,"-",a)
}
event_with_cov <- event_with_cov %>%
  mutate(
    g1 = toupper(g1), g2 = toupper(g2),
    g1 = ifelse(g1 %in% c("N","C","F"), g1, NA),
    g2 = ifelse(g2 %in% c("N","C","F"), g2, NA),
    pair_cat = mapply(undirected_pair, g1, g2)
  )

for (nm in c("Zedge_NN","Zedge_NC","Zedge_NF","Zedge_CC","Zedge_CF","Zedge_FF")) {
  event_with_cov[[nm]] <- as.integer(event_with_cov$pair_cat ==
                                       sub("^Zedge_", "", nm) | FALSE)
}

# ================= 使用示例 =================
# 假设：
#   ant_df   <- 你的“宽表”（含 "Ant Pair", Day1..Day41）
#   nodes_df <- 你的节点表（含 tag_id, body_size, age, group_period1..4）
# 1) 宽 -> 长（逐事件 + 每对一条 death）
 event_long <- wide_to_event_long_with_death(time_data)
event_long <- augment_event_long_with_ag_columns(event_long)
safe_num <- function(x) {
  # 把 factor/character 统一转字符，清理 "NA"/"na"/"Na"/""，再安全转数值
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
# 2) 并入协变量（六类边虚拟量 + 两端体重/年龄）
 event_with_cov <- add_covariates_to_events(event_long, cov_data)
 # 常规：仅向下填充（death 行用上一行的协变量）
 res <- locf_covariates_within_edge( event_with_cov )
 event_with_cov<- res$data
 
 anyNA(event_with_cov)              # TRUE/FALSE
 sum(is.na(event_with_cov))         # NA 的总个数
 
 dplyr::glimpse(event_with_cov)
 
 ###填充完之后还是有协变量为NA
 #简单起见用上一行填上下一行
 suppressPackageStartupMessages({
   library(dplyr)
   library(tidyr)
 })
 
 # 1) 把字符串 "NA" 先转成真正的 NA
 event_with_cov <- event_with_cov %>%
   mutate(
     # 字符列："NA"/"na"/"Na" -> NA
     across(where(is.character), ~ replace(., . %in% c("NA","na","Na"), NA_character_)),
     # 因子列：转字符替换后再转回因子
     across(where(is.factor), ~ {
       x <- as.character(.)
       x[x %in% c("NA","na","Na")] <- NA_character_
       factor(x)
     })
   )
 
 # 2) 全表按当前顺序做“向下填充”（上一行顶上来）
 event_with_cov <- event_with_cov %>%
   tidyr::fill(everything(), .direction = "down")
 anyNA(event_with_cov)              # TRUE/FALSE
 sum(is.na(event_with_cov))         # NA 的总个数
 
save(event_with_cov,file="ant_data_longformat.Rdata")
load("ant_data_longformat.Rdata")
event_with_cov <- event_with_cov %>% select(-c("Zedge_NF"))
event_with_cov <- event_with_cov %>% select(-c("group_1","group_2"))

##########2.执行程序-----------------------
###########======================

#只取一块
{
  library(dplyr)
  thr <- 5000  # 自己设，例如 80
  
  # 每条边的端点编号（假设 event_with_cov 里已有 id1、id2）
  edge_nodes <- event_with_cov %>%
    group_by(id) %>%
    summarise(id1 = first(id1), id2 = first(id2), .groups = "drop")
  
  # 需要删除的边：id1>thr 且 id2>thr
  drop_ids <- edge_nodes %>%
    filter(id1 > thr & id2 > thr) %>%
    pull(id)
  
  # 过滤出保留的子集
  event_sub <- event_with_cov %>%
    filter(!(id %in% drop_ids))
  
  # （可选）把 id 重新连续编号 1..K
  id_map <- event_sub %>% distinct(id) %>% arrange(id) %>% mutate(id_new = row_number())
  event_sub <- event_sub %>% left_join(id_map, by = "id") %>%
    mutate(id = id_new) %>% select(-id_new)
  
  # 小检查
  c(before = n_distinct(event_with_cov$id), after = n_distinct(event_sub$id))
}

real_data <- list()
real_data$df <- event_sub
###重新编号id,因此最终结果的id是重编号过的，需要函数映射回去
{
  library(dplyr)
  library(stringr)
  
  # 如果 id1/id2 不是数值，先转
  event_sub <- event_sub %>%
    mutate(id1 = as.integer(id1),
           id2 = as.integer(id2))
  
  # 1) 节点重编号：old_node -> node (连续 1..N)
  node_levels <- sort(unique(c(event_sub$id1, event_sub$id2)))
  node_map <- tibble(old_node = node_levels,
                     node     = seq_along(node_levels))
  
  event_sub <- event_sub %>%
    # id1 映射
    left_join(node_map, by = c("id1" = "old_node")) %>%
    rename(id1_new = node) %>%
    # id2 映射
    left_join(node_map, by = c("id2" = "old_node")) %>%
    rename(id2_new = node) %>%
    # 覆盖为新节点编号
    mutate(id1 = id1_new, id2 = id2_new) %>%
    select(-id1_new, -id2_new) %>%
    # 无向边规范：小-大
    mutate(i = pmin(id1, id2),
           j = pmax(id1, id2),
           id1 = i, id2 = j) %>%
    select(-i, -j)
  
  # 2) 边重编号：old_id -> id (连续 1..m)
  # 保存旧 id 以免冲突
  event_sub <- event_sub %>% rename(old_id = id)
  
  edge_map <- event_sub %>%
    distinct(old_id, id1, id2) %>%
    arrange(old_id) %>%
    mutate(id = row_number())          # 连续 1..m
  
  event_sub <- event_sub %>%
    left_join(edge_map %>% select(old_id, id), by = "old_id") %>%
    select(-old_id) %>%
    arrange(id, start, stop)
  
  # 3) 生成 pairs 矩阵并写入 real_data
  pairs_mat <- event_sub %>%
    distinct(id, id1, id2) %>%
    arrange(id) %>%
    select(id1, id2) %>%
    as.matrix()
  storage.mode(pairs_mat) <- "integer"
  
  real_data$df    <- event_sub
  real_data$pairs <- pairs_mat
  
  # 4) 小检查：节点/边是否连续无间隔
  N <- length(unique(c(event_sub$id1, event_sub$id2)))
  stopifnot(identical(sort(unique(c(event_sub$id1, event_sub$id2))), seq_len(N)))
  M <- nrow(real_data$pairs)
  stopifnot(identical(sort(unique(event_sub$id)), seq_len(M)))
  
  # 看看结果
  head(real_data$pairs)
  
}

kk <- real_data$pairs

df <- real_data$df
fit_edge_frailty_ppl_block_pcg_tv_fast(
  data_sim = real_data,           # sim 需要包含：df, pairs, Z_i
  zi1_pattern = "e_1$",     # 端点1：匹配 body_size_1, age_1, 以及你以后加的其它 *_1
  zi2_pattern = "e_2$",     # 端点2：匹配 body_size_2, age_2, 以及你以后加的其它 *_2
  max_outer = 10,
  beta_maxit = 10,
  tol = 1e-3,
  damping = 0.8,
  chol_exact_threshold = 5000,
  ties_method = "breslow"  #breslow,efron
)


######用cpp的函数

Rcpp::sourceCpp("edge_frailty_core.cpp")
Rcpp::sourceCpp("approx_diag_K.cpp")

#后来改成了：
Rcpp::sourceCpp("RCPP/fit_edge_frailty.cpp")
# 仍然使用你现有的 build_tv_scanline_index() 与 update_sigma_rho_ratio_fix()
df <- real_data$df
cols <- c("Zedge_NN","Zedge_NC","Zedge_CC","Zedge_CF","Zedge_FF")

# dplyr 写法
library(dplyr)
real_data$df <- real_data$df %>%
  mutate(across(any_of(cols), ~ . - mean(., na.rm = TRUE)))

# 快速检查均值是否≈0
colMeans(real_data$df[intersect(cols, names(real_data$df))], na.rm = TRUE)

fit_real <- fit_edge_frailty_ppl_block_pcg_tv_fast_cpp(
  data_sim = real_data,           # sim 需要包含：df, pairs, Z_i
  zi1_pattern = "_1$",     # 端点1：匹配 body_size_1, age_1, 以及你以后加的其它 *_1
  zi2_pattern = "_2$",     # 端点2：匹配 body_size_2, age_2, 以及你以后加的其它 *_2
  max_outer = 25,
  beta_maxit =10,
  tol = 1e-3,
  damping = 0.9,
  ties_method = "efron"  #breslow,efron
)

real_data_10#
 save(fit_real,file="real_data_fit_11-13.Rdata")
 load("real_data_fit_11-13.Rdata")
#####check id
{
  # df：你准备用来拟合（或已经在拟合）的数据框（event_with_cov 或其清理版）
  df <- event_with_cov
  
  # 边数 m（按 id 统计）
  pairs_tbl <- df |> dplyr::distinct(id, id1, id2) |> dplyr::arrange(id)
  m_df <- nrow(pairs_tbl)
  cat("m_df =", m_df, "\n")
  
  # Zedge 列（仅用于确认 beta 维度是否也和 scan 一致；不是本次核心）
  Xe_cols <- grep("^Zedge_", names(df), value = TRUE)
  cat("Zedge 列数 =", length(Xe_cols), "\n")
  
  # 如果 scan 存在，读出它的 m,K；没有就先设占位
  if (exists("scan")) {
    m_scan <- if (!is.null(scan$m)) scan$m else NA_integer_
    K_scan <- if (!is.null(scan$K)) scan$K else NA_integer_
    cat("scan$m =", m_scan, " scan$K =", K_scan, "\n")
    if (!is.na(m_scan)) stopifnot(m_scan == m_df)
  } else {
    cat("scan 不存在（或未更新）\n")
  }
  
  # 如果 ops 存在，检查维度
  if (exists("ops") && !is.null(ops$S2)) {
    cat("dim(ops$S2) =", paste(dim(ops$S2), collapse="x"), "\n")
    stopifnot(all(dim(ops$S2) == c(m_df, m_df)))
  } else cat("ops$S2 不存在（或未更新）\n")
  
  # 如果 delta 存在，检查长度是否等于 m 或 m*K
  if (exists("delta")) {
    len_delta <- length(delta)
    cat("length(delta) =", len_delta, "\n")
  }
  
}
  


#########对他进行test
df <-real_data$df
scan <- build_tv_scanline_index(real_data$df, zedge_pat="^Zedge_", 
                                zi1_pat="^body_size_", zi2_pat="^age_")
n    <- if (!is.null(real_data$Z_i)) nrow(real_data$Z_i) else max(real_data$pairs)
ops  <- make_linegraph_ops(n, real_data$pairs)
S2   <- make_S2(n, real_data$pairs)
ties_method  <- "efron"
kern  <- make_scan_kernel(scan, fit_real$beta_all,ties_method = ties_method)
lam <- compute_cumLambda(scan, fit_real$beta_all, fit_real$delta, method = ties_method)  # 或 "efron"

Kdiag <- compute_K_diag_full(
  scan, fit_real$beta_all, fit_real$delta,
  cumLambda_by_time = lam$cumLambda_by_time
)
Kdiag_real <- compute_K_diag_partial(scan,fit_real$beta_all,fit_real$delta,ties_method = "efron")
S2_real   <- make_S2(162, real_data$pairs)

fit_real$beta_all
scan$X_rows


scan$times

#Rcpp::sourceCpp("RCPP/wald_ci_sigma2_rho_cpp.cpp")
S2_real <- methods::as(S2_real, "dgCMatrix")


ci_out_real  <- wald_ci_sigma2_rho(
  delta = fit_real$delta, n = n, S2 = S2_real,
  sigma2_hat =  fit_real$sigma^2, rho_hat = fit_real$rho,
  n_mc = 32L, n_subspace = 8L, 
  K_diag = Kdiag_real,  chol_threshold = 5000
)

# ci_out <- wald_ci_sigma2_rho_cpp(
#   delta = fit_real$delta, n = n, S2 = S2_real,
#   sigma2_hat =  fit_real$sigma^2, rho_hat = fit_real$rho,
#   K_diag = Kdiag_real,  chol_threshold = 5000
# )  这个有点问题

ci_out_real$ci_sigma2   # σ² 的 95% CI
ci_out_real $ci_rho      # ρ 的 95% CI（已投影到可行域）
c_sigma <- ci_out_real $ci_sigma2 
ci_out_real $se
ci_out_real $ci_rho


####回归参数的推断
{
  #####对回归参数推断
  {
    Rcpp::sourceCpp(file="RCPP/infer_beta.cpp")
    Rcpp::sourceCpp(file="RCPP/schur_theta_block.cpp")
    
   
    
    # 需要先 sourceCpp("src/frailty_ops_cpp.cpp") 或把它编进包里
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
      # =========================================================
      # Unified vcov for θ=(β,γ)
      #   - sand = FALSE: Var = [H^{-1}]_{θθ}
      #   - sand = TRUE : Var = [H^{-1}]_{θθ} * V_cluster * [H^{-1}]_{θθ}
      # cluster ∈ {"dyad"(default), "gray", "edge", "node"}
      # Fast: single scan for scores; batched Htd/Hdt; O(p) solves of H_dd.
      # =========================================================
      cox_frailty_theta_vcov <- function(
    scan, pairs,
    beta_hat,        # θ-hat: c(β, γ) in the order of scan$X_rows columns
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
        ## 定义一个本地帮助函数：稳定地计算 (B^T B) X
        apply_BTB <- function(X) {
          X <- if (is.matrix(X)) X else cbind(X)
          if (!is.null(ops$B)) {
            as.matrix(crossprod(ops$B, ops$B %*% X))   # == (S2 + 2I) X
          } else if (!is.null(ops$S2_matmat) && is.function(ops$S2_matmat)) {
            ops$S2_matmat(X) + 2 * X                   # 备用路径
          } else {
            stop("Neither ops$B nor ops$S2_matmat function is available.")
          }
        }
        
        X  <- scan$X_rows
        p  <- ncol(X)
        m  <- length(delta_hat)
        theta_hat <- beta_hat
        
        # ---- indices for β / γ ----
        beta_idx  <- grep(beta_col_pattern,  colnames(X), value = FALSE)
        gamma_idx <- grep(gamma_col_pattern, colnames(X), value = FALSE)
        if (!length(beta_idx))  stop("No beta columns matched: ", beta_col_pattern)
        if (!length(gamma_idx)) {
          gamma_idx <- setdiff(seq_len(p), beta_idx)
          if (!length(gamma_idx)) stop("Gamma columns not found; set gamma_col_pattern.")
        }
        
        # ---- line-graph ops & Ω^{-1} solver (Ω^{-1} = (1/σ^2) M^{-1}) ----
        n_node <- max(pairs)
        ops    <- make_linegraph_ops(n_node, pairs)
        Msolve <- make_M_solver(ops, rho, chol_threshold, pcg_tol, pcg_maxit)
        OmegaInv_apply <- function(x) (1/sigma2) * Msolve(x)
        
        # ---- scan kernel & K(·) ----
        sk <- make_scan_kernel(scan, beta = theta_hat, ties_method = ties_method)
        K_apply <- function(x) sk$K_matvec(x, delta = delta_hat)
        
        # H_dd = K + Ω^{-1}
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
        # ---- build batched θ–δ operators; fallback to vector version if needed ----
        Htt_fallback <- function() {
          H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
          H_ops$Htt_matrix()
        }
        Htd_matmul_fallback <- function(U){  # U: m×q
          H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
          U <- if (is.matrix(U)) U else cbind(U)
          do.call(cbind, lapply(seq_len(ncol(U)), function(j) H_ops$Htd_matvec(U[, j])))
        }
        Hdt_matmul_fallback <- function(V){  # V: p×q
          H_ops <- build_theta_delta_ops(scan, theta_hat, delta_hat, ties_method)
          V <- if (is.matrix(V)) V else cbind(V)
          do.call(cbind, lapply(seq_len(ncol(V)), function(j) H_ops$Hdt_matvec(V[, j])))
        }
        
        I_tt        <- Htt_fallback() + diag(ridge, p)
        Htd_matmul  <- Htd_matmul_fallback
        Hdt_matmul  <- Hdt_matmul_fallback
        #####C形式的函数只在实际数据中用到了
        # if (exists("build_theta_delta_ops_matmul") && is.function(build_theta_delta_ops_matmul)) {
        #   TDops <- build_theta_delta_ops_matmul(scan, theta_hat, delta_hat, ties_method)
        #   I_tt        <- TDops$Htt_matrix() + diag(ridge, p)
        #   Htd_matmul  <- TDops$Htd_matmul
        #   Hdt_matmul  <- TDops$Hdt_matmul
        # } else {
        #   I_tt        <- Htt_fallback() + diag(ridge, p)
        #   Htd_matmul  <- Htd_matmul_fallback
        #   Hdt_matmul  <- Hdt_matmul_fallback
        # }
        
        
        ###实际数据到这就行了
        # {
        #   TDops <- build_theta_delta_ops_matmul(scan, beta_hat, delta_hat,
        #                                         ties_method = "efron")
        #   
        #   I_tt <- TDops$Htt_matrix() + diag(1e-8, ncol(scan$X_rows))
        #   
        #   return(I_tt)
        #   res <- schur_theta_block_cpp_min(
        #     I_tt,
        #     TDops$Hdt_matmul,     # 注意：这里直接把 R 函数句柄传给 C++
        #     Hdd_solve_mat,        # 你已有的 RHS->解 的批量PCG函数
        #     TDops$Htd_matmul,
        #     ridge = 1e-8, ridge_tries = 5, block_size = 32, verbose = TRUE
        #   )
        #   S_inv <- res$S_inv
        #   
        # }
        # ---- bread: S = H_tt - H_td H_dd^{-1} H_dt  (θθ-Schur of H), S^{-1} = [H^{-1}]_{θθ}
        Ip <- diag(1, p)
        Q  <- Hdt_matmul(Ip)            # m×p
        Z  <- Hdd_solve_mat(Q)          # m×p
        S  <- I_tt - Htd_matmul(Z)      # p×p
        cf_S <- tryCatch(chol(S), error = function(e) NULL)
        if (is.null(cf_S)) { S <- S + diag(ridge, p); cf_S <- chol(S) }
        S_inv <- chol2inv(cf_S)         # [H^{-1}]_{θθ}
        
        # ---- sand = FALSE: return model-based and done ----
        if (!sand) {
          Var_theta <- S_inv
        } else {
          # --------- build θ-space meat Mid depending on 'cluster' ----------
          if (cluster == "gray") {
            # Gray: J = I_tt - H_td H_dd^{-1} I_dt - I_td H_dd^{-1} H_dt + H_td H_dd^{-1} I_dd H_dd^{-1} H_dt
            # where I_dt = Hdt_matmul, I_td = Htd_matmul, I_dd = K_apply
            Mid <- matrix(0.0, p, p)
            for (j in 1:p) {
              ej <- rep(0.0, p); ej[j] <- 1
              col1 <- I_tt[, j]
              z1   <- Hdd_solve(Hdt_matmul(ej))     # H_dd^{-1} I_dt e_j
              t12  <- Htd_matmul(z1)                # I_td H_dd^{-1} I_dt e_j
              z2   <- Hdd_solve(K_apply(z1))        # H_dd^{-1} I_dd H_dd^{-1} I_dt e_j
              t3   <- Htd_matmul(z2)
              Mid[, j] <- col1 - 2 * t12 + t3
            }
          } else {
            
            # edge-level joint scores U_e = (u_{θe}, u_{δe})  —— single scan
            es <- edge_joint_scores(scan, theta_hat, delta_hat, ties_method)
            Utheta_e <- as.matrix(es$U_theta_by_edge)  # m×p
            udelta   <- as.numeric(es$U_delta_by_edge) # m
            
            # Common pieces
            # H_dt(I) already computed: Q, Z (m×p)
            # θθ part accumulates into Mid; we'll add the cross terms and dd-term
            
            if (cluster == "edge") {
              M_tt <- crossprod(Utheta_e)                 # p×p
              M_dt <- Utheta_e * udelta                   # m×p
              Mdd_apply <- function(X) (udelta^2) * (if (is.matrix(X)) X else cbind(X))
            } else if (cluster == "node") {
              BTB_Uθ <- apply_BTB(Utheta_e)               # m×p
              M_tt   <- crossprod(Utheta_e, BTB_Uθ)       # p×p
              M_dt   <- BTB_Uθ * udelta                   # m×p
              Mdd_apply <- function(X) {
                X  <- if (is.matrix(X)) X else cbind(X)
                DX <- X * udelta
                udelta * apply_BTB(DX)                    # D (B^T B) (D X)
              }
            } else { # "dyad" = 2*nodes - edges
              c_v <- if (n_node > 1) n_node/(n_node - 1) else 1
              c_e <- if (m      > 1) m     /(m      - 1) else 1
              
              BTB_Uθ     <- apply_BTB(Utheta_e)           # m×p
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
                2 * (udelta * BTB_DX) - (udelta^2) * X    # 2 D (B^T B) D X - D^2 X
              }
            }
            
            
            # Map meat to θ-space: Mid = M_tt - H_td H_dd^{-1} M_dt - M_td H_dd^{-1} H_dt + H_td H_dd^{-1} M_dd H_dd^{-1} H_dt
            Mid <- M_tt
            Y   <- Hdd_solve_mat(M_dt)                 # m×p
            Mid <- Mid - Htd_matmul(Y)                 # - H_td H_dd^{-1} M_dt
            Mid <- Mid - crossprod(M_dt, Z)            # - M_td H_dd^{-1} H_dt
            T2  <- Hdd_solve_mat(Mdd_apply(Z))         # m×p
            Mid <- Mid + Htd_matmul(T2)                # + H_td H_dd^{-1} M_dd H_dd^{-1} H_dt
          }
          
          Var_theta <- S_inv %*% Mid %*% S_inv
        }
        Var_theta_unad <- S_inv
        take <- function(M, r, c) M[r, c, drop = FALSE]
        Var_beta  <- take(Var_theta, beta_idx,  beta_idx)
        Var_gamma <- take(Var_theta, gamma_idx, gamma_idx)
        Cov_bg    <- take(Var_theta, beta_idx,  gamma_idx)
        
        se <- sqrt(pmax(diag(Var_theta), 0))
        
        ####为什么这部分调整之后beta的方差会被估计小呢
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
    inf <- cox_frailty_theta_vcov (
      scan, pairs = real_data$pairs,
      beta_hat =  fit_real$beta_all, delta_hat = fit_real$delta,
      sand = FALSE, cluster = "node",
      K_diag = Kdiag_real,
      sigma2 = (fit_real$sigma)^2, rho = fit_real$rho,
      ties_method = "efron",          # 或 "efron"
      beta_col_pattern  = "^Zedge_",
      gamma_col_pattern = "^Xi_plus_Xj_"
    )
    save(inf,file="informationmatrix_real.RData")
    inf
    
    install.packages("RcppProgress")
    
    # 取 β/γ 的稳健 SE（dyadic two-way）
    ###这里的inf 是信息矩阵，接下来还要算按node聚合的方差
    Rcpp::sourceCpp("RCPP/edge_joint_scores_cpp.cpp")
    res <- edge_joint_scores_cpp(
      X_rows = scan$X_rows,
      id_row  = scan$id_row,
      start_sorted = scan$start_sorted,
      stop_sorted  = scan$stop_sorted,
      rows_start = scan$rows_start,
      rows_stop  = scan$rows_stop,
      times = scan$times,
      events_rowidx_by_time = scan$events_rowidx_by_time,
      beta  = fit_real$beta_all,
      delta = fit_real$delta,
      ties_method = "efron",
      S2_ = S2,          # 或者用 S2_ = ops$S2
      show_progress = TRUE,   # 显示进度条
      check_every   = 2048    # 大任务时可以调大，降低中断检查开销
    )
    
    Utheta_e <- res$U_theta_by_edge      # m×p
    udelta   <- res$U_delta_by_edge      # m
    BTB_Uθ   <- res$BTB_Utheta           # 若提供了 B 或 S2
    M_tt     <- res$M_tt                 # 若计算了 BTB_Uθ
 
    Var_theta_real <- solve(inf) %*%  M_tt %*% solve(inf)
    gamma_sd <- sqrt(diag(Var_theta_real))[c(6:7)]
    beta_sd <- sqrt(diag(solve(inf)))[c(1:5)]
    beta_real <- fit_real$beta_edge
    colnames(scan$X_rows)
    names(beta_real) <- c("Nurse-Nurse","Nurse-Cleaner","Cleaner-Cleaner","Cleaner-Forager","Forager-Forager")
    gamma_real <- fit_real$gamma
    names(gamma_real) <- c("bodysize","age")
    
    fit_real$sigma
    fit_real$rho
  }
}
####11-13的结果，回归系数的推断没改


####用rereg试一下

{
  library(reReg)
  df <- real_data$df
  df$body <- df$body_size_1+df$body_size_2
  #df$event[df$death == 1] <- 0   ###为了Rereg方法
  
  form_str <- "Recur(stop, id, event) ~ Zedge_NN+  Zedge_NC+ Zedge_CC+ Zedge_CF+ Zedge_FF +body "
  form <- as.formula(form_str)
  
  fit <- reReg(form, data = df, model = "cox.LWYY")
  fit$par1.se
   coef(fit)
   
   fit <- coxph(Surv(start, stop, event) ~ Zedge_NN+  Zedge_NC+ Zedge_CC+ Zedge_CF+ Zedge_FF +body ,
                data = df ,
                method = "efron",
                control = coxph.control(timefix = FALSE))   #ties = "efron"
   coef(fit)
}


######计算每条边 平均强度
###=================================从这里开始作图
################============================

#######==================制作这张表，参数的估计值和置信区间
{
  {
    ## --- inputs assumed available: beta_real, gamma_real, beta_sd, gamma_sd ---
    
    # 1) name alignment
    names(beta_sd)  <- names(beta_real)
    names(gamma_sd) <- names(gamma_real)
    
    # 2) CI
    conf_level <- 0.95
    z <- qnorm(1 - (1 - conf_level)/2)
    
    vars <- c(names(beta_real), names(gamma_real))
    ests <- c(as.numeric(beta_real), as.numeric(gamma_real))
    sds  <- c(as.numeric(beta_sd[names(beta_real)]),
              as.numeric(gamma_sd[names(gamma_real)]))
    
    lwr <- ests - z * sds
    upr <- ests + z * sds
    
    # 3) p-values (two-sided Wald); protect against 0 or NA SE
    zval <- ests / sds
    p_raw <- 2 * (1 - pnorm(abs(zval)))
    p_raw[!is.finite(p_raw)] <- NA
    
    fmt_p <- function(p) ifelse(is.na(p), NA,
                                ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
    
    # 4) English table
    tab_en <- data.frame(
      Variable = vars,
      Estimate = sprintf("%.4f", ests),
      `95% CI` = sprintf("[%.4f, %.4f]", lwr, upr),
      `p-value` = fmt_p(p_raw),
      check.names = FALSE
    )
    
    # 5) LaTeX code (booktabs)
    latex_code <- knitr::kable(
      tab_en, format = "latex", booktabs = TRUE,
      align = c("l","r","r","r"),
      caption = "Parameter estimates with 95\\% confidence intervals and p-values"
    )
    
    # 6) output
    cat(latex_code)
    dir.create("results", showWarnings = FALSE, recursive = TRUE)
    writeLines(latex_code, "results/table_parameters.tex")
  }
  
}



###可视化一下累积风险函数，发现近似线性的
{
  df_L <- data.frame(time = lam$times, value = lam$cumLambda_by_time)
  ggplot(df_L, aes(time, value)) +
    geom_step(direction = "hv", linewidth = 1) +
    labs(x = "time", y = expression(Lambda(t)), title = "累计危险强度 Λ(t)") +
    theme_minimal(base_size = 14)
  
}

####可视化fit$delta
{
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  plot_delta <- function(delta, bins = 5, trim_q = c(0.00001, 0.99999),
                         base_size = 16,
                         title = "Distribution of edge frailties δ") {
    df <- data.frame(delta = as.numeric(delta)) %>%
      filter(is.finite(delta))
    
    # 可选：去掉极端尾部，避免几条极端值把形状拉扁
    if (!is.null(trim_q)) {
      qs <- quantile(df$delta, probs = trim_q, na.rm = TRUE)
      df <- df %>% filter(delta >= qs[1], delta <= qs[2])
    }
    
    mu  <- mean(df$delta); med <- median(df$delta); sdv <- sd(df$delta)
    
    ggplot(df, aes(x = delta)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = bins, fill = "#F28E8C", color = "white",
                     alpha = 0.9, linewidth = 0.2) +
      geom_density(linewidth = 1.2, color = "#C0392B", adjust = 1.0) +
      geom_vline(xintercept = mu,  linetype = "dashed", color = "#2C7FB8", linewidth = 0.9) +
      geom_vline(xintercept = med, linetype = "dotted", color = "#2E8B57", linewidth = 0.9) +
      labs(x = expression(delta[ij]), y = "Density",
           title = title,
           subtitle = sprintf("n = %s, sd = %.3f", comma(nrow(df)), sdv)) +
      theme_minimal(base_size = base_size) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot")
  }
  
  p_delta <- plot_delta(fit_real$delta, bins = 60, trim_q = c(0.001, 0.9999))
  p_delta
  # ggsave("delta_dist.png", p_delta, width = 7, height = 4.5, dpi = 300)
  
}


# 计算每条边的“时间平均强度”：exp( beta^T * x_bar_ij + delta_ij )
#####分别计算1.带frailty的和不带的拟合值
{
edge_avg_intensity <- function(df, beta_all, delta_edge,
                               zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_") {
  idx <- build_tv_scanline_index(df, zedge_pat=zedge_pat, zi1_pat=zi1_pat, zi2_pat=zi2_pat)
  X <- idx$X_rows                         # 行级设计矩阵（已按 Zedge_* 与 Z1i_+Z2i_ 拼好）
  id <- idx$id_row                        # 每行对应的边 id（整数编码）
  dt <- df$stop - df$start               # 区间长度，作为时间权重
  if (any(dt <= 0)) stop("发现非正的区间长度（stop - start ≤ 0）。")
  
  # —— 逐边的时间权重和 —— #
  w_by_edge <- rowsum(dt, group = id)     # m×1
  # 逐列做加权和，然后除以时间权重 => 加权均值协变量
  Xw_sum <- rowsum(X * dt, group = id)    # m×p
  Xbar   <- Xw_sum / as.numeric(w_by_edge)
  
  # —— 确保 beta 与列名对齐（若 beta 有 names 则按列名重排）—— #
  if (!is.null(names(beta_all)) && !is.null(colnames(Xbar))) {
    miss <- setdiff(colnames(Xbar), names(beta_all))
    if (length(miss)) stop("beta_all 缺少这些系数：", paste(miss, collapse=", "))
    beta_use <- beta_all[colnames(Xbar)]
  } else {
    if (length(beta_all) != ncol(Xbar)) stop("beta_all 长度与设计矩阵列数不一致。")
    beta_use <- beta_all
  }
  
  # —— 将 delta 对齐到边 id 顺序 —— #
  # rowsum 的行名就是边 id（字符型）。如果 delta 有 names，就按行名重排；否则按 id 假定为 1..m 的顺序。
  rid <- rownames(Xbar)
  if (!is.null(names(delta_edge))) {
    if (!all(rid %in% names(delta_edge))) stop("delta_edge 的命名不含全部边 id。")
    delta_use <- as.numeric(delta_edge[rid])
  } else {
    if (length(delta_edge) != nrow(Xbar)) stop("delta_edge 长度与边数不一致。")
    delta_use <- as.numeric(delta_edge[as.integer(rid)])
  }
  
  # 线性预测 + 脆弱项，再取指数
  eta_bar <- as.vector(Xbar %*% beta_use) + delta_use
  lambda_bar <- exp(eta_bar)
  names(lambda_bar) <- rid
  lambda_bar
}
sd(fit_real$delta)
fit_real$sigma
df <- real_data$df
lambda_bar <- edge_avg_intensity(df, beta_all = fit_real$beta_all, delta_edge = fit_real$delta,
                                 zedge_pat="^Zedge_", zi1_pat="^body_size_", zi2_pat="^age_")


# 得到一个长度为“边数”的向量（带边 id 名称）
#save( lambda_bar,file="lambda_bar.RData")
# 计算 “平均强度 × 存活时长”
# method="scan" 表示使用 build_tv_scanline_index() 的并集时刻/分组来找每条边的死亡时刻
edge_integrated_intensity <- function(
    df, beta_all, delta_edge,
    zedge_pat="^Zedge_", zi1_pat="^Z1i_", zi2_pat="^Z2i_",
    method = c("scan","df","day")  # "day": 由 death==1 行的 day 列确定 T
) {
  method <- match.arg(method)
  
  idx <- build_tv_scanline_index(df, zedge_pat=zedge_pat, zi1_pat=zi1_pat, zi2_pat=zi2_pat)
  X   <- idx$X_rows
  id  <- idx$id_row
  
  # -------- 1) 每条边的死亡时刻 t_death_by_id --------
  rid_levels <- sort(unique(id))
  t_death_by_id <- rep(Inf, length(rid_levels))
  names(t_death_by_id) <- as.character(rid_levels)
  
  if (method == "scan") {
    times_all <- idx$times_all_events
    grp_all   <- idx$events_or_death_rowidx_by_time
    for (kk in seq_along(times_all)) {
      rows_k <- grp_all[[kk]]
      if (!length(rows_k)) next
      dead_rows <- rows_k[df$death[rows_k] == 1L]
      if (!length(dead_rows)) next
      dead_ids <- unique(id[dead_rows])
      need_set <- as.character(dead_ids[t_death_by_id[as.character(dead_ids)] == Inf])
      if (length(need_set)) t_death_by_id[need_set] <- times_all[kk]
    }
  } else if (method == "df") {
    death_stop <- ifelse(df$death == 1L, df$stop, NA_real_)
    tlist <- tapply(death_stop, id, function(v) {
      v <- v[is.finite(v)]
      if (length(v)) min(v) else Inf
    })
    t_death_by_id[names(tlist)] <- tlist
  } else { # method == "day"
    if (!("day" %in% names(df))) stop("method='day' 需要 df$day 列。")
    # 用 death==1 那行的 day 作为死亡时刻；没有死亡就用该边出现过的最大 day
    death_day <- ifelse(df$death == 1L, df$day, NA_real_)
    min_death_day_by_id <- tapply(death_day, id, function(v) {
      v <- v[is.finite(v)]
      if (length(v)) min(v) else NA_real_
    })
    max_day_by_id <- tapply(df$day, id, max, na.rm = TRUE)
    
    # 填入：有死亡用最小 death day；否则用该边的最大 day
    fill_names <- intersect(names(t_death_by_id), names(max_day_by_id))
    t_death_by_id[fill_names] <- ifelse(
      is.na(min_death_day_by_id[fill_names]),
      max_day_by_id[fill_names],
      min_death_day_by_id[fill_names]
    )
  }
  
  # 每行对应的死亡截断时刻（用于加权均值计算）
  t_death_row <- t_death_by_id[as.character(id)]
  
  # -------- 2) 存活期权重：截断到 t_death_row 之内 --------
  alive_len <- pmax(pmin(df$stop, t_death_row) - df$start, 0)
  Xw_sum <- rowsum(X * alive_len, group = id)       # m × p
  T_alive_sum <- rowsum(alive_len, group = id)      # m × 1（区间求和）
  
  # —— 用 sum(alive_len) 做均值协变量的分母 —— #
  T_sum_vec <- as.numeric(T_alive_sum)
  names(T_sum_vec) <- rownames(T_alive_sum)
  denom <- ifelse(T_sum_vec > 0, T_sum_vec, 1)
  Xbar  <- Xw_sum / denom
  if (any(T_sum_vec == 0)) Xbar[T_sum_vec == 0, ] <- 0
  
  # -------- 3) 系数/脆弱项对齐并计算 “平均强度” --------
  if (!is.null(names(beta_all)) && !is.null(colnames(Xbar))) {
    miss <- setdiff(colnames(Xbar), names(beta_all))
    if (length(miss)) stop("beta_all 缺少系数：", paste(miss, collapse=", "))
    beta_use <- beta_all[colnames(Xbar)]
  } else {
    if (length(beta_all) != ncol(Xbar)) stop("beta_all 长度与设计矩阵列数不一致")
    beta_use <- beta_all
  }
  
  rid <- rownames(Xbar)  # 边 id（字符）
  if (!is.null(names(delta_edge))) {
    if (!all(rid %in% names(delta_edge))) stop("delta_edge 命名不含全部边 id")
    delta_use <- as.numeric(delta_edge[rid])
  } else {
    if (length(delta_edge) != nrow(Xbar)) stop("delta_edge 长度与边数不一致")
    delta_use <- as.numeric(delta_edge[as.integer(rid)])
  }
  
  eta_bar <- as.vector(Xbar %*% beta_use) + delta_use
  lambda_bar <- exp(eta_bar)  # 存活期内的“平均强度”

  # -------- 4) 存活时长 T 的取法 --------
  # - method=="day"：T 直接取每条边的死亡 day（无死亡取其最大 day）
  # - 其他方法：T 取 sum(alive_len)（按截断后的实际覆盖时长）
  if (method == "day") {
    T_vec <- as.numeric(t_death_by_id[rid])
    names(T_vec) <- rid
  } else {
    T_vec <- T_sum_vec[rid]
  }
  
  # -------- 5) 返回 “平均强度 × 存活时长” 向量 --------
  lambda_int <- lambda_bar * T_vec
  names(lambda_int) <- rid
  
  # 附带组件，便于检查
  attr(lambda_int, "T_used") <- T_vec                # 本次用于相乘的 T
  attr(lambda_int, "T_sum_alive") <- T_sum_vec       # sum(alive_len)（供对比）
  attr(lambda_int, "lambda_bar") <- lambda_bar
  attr(lambda_int, "t_death_by_id") <- t_death_by_id
  attr(lambda_int, "covariate") <- exp(as.vector(Xbar %*% beta_use))
  
  lambda_int
}


df <- real_data$df
lambda_int <- edge_integrated_intensity(
  df,
  beta_all = fit_real$beta_all,
  delta_edge = fit_real$delta,
  method = "day" ,  # 或 "df"
  zedge_pat="^Zedge_", zi1_pat="^body_size_", zi2_pat="^age_"
)
# 得到长度为“边数”的向量（带边 id 名称）
# plot(lambda_int)
# lambda_bar <- lambda_int
# lambda_bar <- lambda_int
# save( lambda_bar,file="lambda_bar.RData")
# load("lambda_bar.RData")
# 右连续阶梯函数 Λ(t) 的评估器
eval_cumhaz_at <- function(lam) {
  stopifnot(!is.null(lam$times), !is.null(lam$cumLambda_by_time))
  tt  <- as.numeric(lam$times)
  cc  <- as.numeric(lam$cumLambda_by_time)
  o <- order(tt)
  tt <- tt[o]; cc <- cc[o]
  # 右连续：method="constant", f=1；域外取 0 和 max 值
  approxfun(tt, cc, method = "constant", f = 1,
            yleft = 0, yright = tail(cc, 1), ties = "ordered")
}

# --- 假设你已经从 edge_integrated_intensity 得到了 ---
lambda_bar <- attr(lambda_int, "lambda_bar")
lambda_cov <- attr(lambda_int, "covariate")
T_vec      <- attr(lambda_int, "T_used")   # 或你自己构造的 T 向量
T_vec <- attr(lambda_int, "T_used")
fΛ <- eval_cumhaz_at(lam)
Λ_T <- fΛ(T_vec)                 # 逐边的 Λ(T_i)
lambda_by_Λ <- lambda_bar * Λ_T  # 你要的结果
lambda_by_cov <- lambda_cov * Λ_T  ##没有frailty的结果
}

library(ggplot2)
library(dplyr)

hist(lambda_by_Λ,binwidth = 5) 
  
# ##########4.看看所有的边交互情况，有没有交互很少的
{
 
  {
    library(ggplot2)
    # 假设数据框名为 time_data
    time_data_long <- time_data  # 原始数据
    
    # 添加一列总交互次数（忽略 NA）
    time_data_long$total_interaction <- rowSums(time_data_long[ , -1], na.rm = TRUE)
    #time_data_long$total_interaction <- rowMeans(time_data_long[ , -1], na.rm = TRUE)
    #time_data_long$total_interaction <- rowSums(time_data_long[ , -c(1:30)], na.rm = TRUE)
    # 排序查看交互较少的边
    low_interactions <- time_data_long[order(time_data_long$total_interaction), ]
    head(low_interactions, 100)  # 查看交互最少的10条边
    ggplot(low_interactions, aes(x = total_interaction)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", boundary = 0, closed = "left") +
      labs(title = "Distribution of Total Interactions ",
           x = "Total Interactions",
           y = "Number of Edges") +
      theme_minimal(base_size = 14)
    
    
    
    ####只统计某个群体中的
    # {
    #   
    #   # Step 1：找出 ss1_common 中 N-N == 1 的蚂蚁对
    #   valid_pairs <- ss1_common %>%
    #     filter(`N-N` == 1) %>%
    #     pull(`Ant Pair`)
    #   
    #   # Step 2：从 time_data 中筛选出这些边
    #   filtered_edges <- time_data %>%
    #     filter(`Ant Pair` %in% valid_pairs)
    #   
    #   # Step 3：计算前 11 天的交互总和
    #   filtered_edges <- filtered_edges %>%
    #     mutate(total_interaction = rowSums(select(., Day1:Day11), na.rm = TRUE))
    #   }
    
  }
}



####三种情形的拟合图
{
  library(ggplot2)
  library(dplyr)
  library(grid)
  
  ## --- 数据 ---
  df_obs  <- data.frame(value = as.numeric(low_interactions$total_interaction)) %>%
    filter(is.finite(value))
  df_with <- data.frame(value = as.numeric(lambda_by_Λ)) %>%
    filter(is.finite(value))
  df_wo   <- data.frame(value = as.numeric(lambda_by_cov)) %>%
    filter(is.finite(value))
  
  df_all <- bind_rows(
    data.frame(x = df_obs$value,  panel = "Observed"),
    data.frame(x = df_with$value, panel = "Fitted with frailty"),
    data.frame(x = df_wo$value,   panel = "Fitted without frailty")
  ) %>%
    filter(is.finite(x), x <= 1000)
  
  ## --- 颜色 ---
  # col_obs   <- "#5DADE2"
  # col_with  <- "#E74C3C"
  # col_wo    <-  "#7DCEA0" 
  ## --- 颜色 ---
  col_obs      <- "#5DADE2"
  col_with     <- "#E74C3C"
  col_wo_fill  <- "#52BE80"   # 浅一点的填充绿（保持不变）
  col_wo_line  <- "#3DA766"   # 稍微深一点的线条绿
  
  cols_fill <- c(
    "Observed"               = col_obs,
    "Fitted with frailty"    = col_with,
    "Fitted without frailty" = col_wo_fill
  )
  
  cols_line <- c(
    "Observed"               = col_obs,
    "Fitted with frailty"    = col_with,
    "Fitted without frailty" = col_wo_line
  )
  
  #27AE60更深的绿色
  
  cols_named <- c(
    "Observed"               = col_obs,
    "Fitted with frailty"    = col_with,
    "Fitted without frailty" = col_wo
  )
  
  ## --- 参数 ---
  bw        <- 10          # ★ binwidth 改为 10
  boundary0 <- 0
  adjust_kde <- 0.5
  
  ## 样本量
  n_obs  <- nrow(df_obs)
  n_with <- nrow(df_with)
  n_wo   <- nrow(df_wo)
  
  ## 左轴上限（直方图最大计数）
  bin_df <- df_all %>%
    mutate(bin_id = floor((x - boundary0) / bw)) %>%
    group_by(panel, bin_id) %>%
    summarise(count = n(), .groups = "drop")
  y_hist_max <- max(bin_df$count, na.rm = TRUE)
  
  ## --- 边界反射的密度函数 ---
  reflect_density <- function(x, adjust, from = 0, to = 950) {
    x <- x[is.finite(x)]
    x_ext <- c(x, -x)                 # 0 处反射
    d <- density(x_ext, adjust = adjust, from = from, to = to)
    d$y <- 2 * d$y                    # 反射修正
    d
  }
  
  ## 计算三条密度（真密度 → 期望 count）
  dens_obs  <- reflect_density(df_obs$value,  adjust_kde)
  dens_with <- reflect_density(df_with$value, adjust_kde)
  dens_wo   <- reflect_density(df_wo$value,   adjust_kde)
  
  df_dens <- bind_rows(
    data.frame(x = dens_obs$x,  y = dens_obs$y  * n_obs  * bw,
               panel = "Observed"),
    data.frame(x = dens_with$x, y = dens_with$y * n_with * bw,
               panel = "Fitted with frailty"),
    data.frame(x = dens_wo$x,   y = dens_wo$y  * n_wo   * bw,
               panel = "Fitted without frailty")
  )
  
  # 密度→count 后的最大值
  max_dens_scaled <- max(df_dens$y, na.rm = TRUE)
  y_max <- max(y_hist_max, max_dens_scaled)
  
  ## 右轴换算用的 N（假定三组样本量差不多）
  N_ref <- n_obs
  
  ## --- 作图 ---
  p <- ggplot() +
    # 直方图（计数）
    geom_histogram(
      data = df_all,
      aes(x = x, y = after_stat(count), fill = panel),
      binwidth = bw, boundary = boundary0, closed = "left",
      position = "identity",
      color = NA, alpha = 0.5,
      show.legend = TRUE
    ) +
    # 反射修正后的密度曲线（期望 count）
    geom_line(
      data = df_dens,
      aes(x = x, y = y, colour = panel),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = cols_named,
      breaks = c("Observed", "Fitted with frailty", "Fitted without frailty"),
      name   = NULL
    ) +
    scale_colour_manual(
      values = cols_named,
      breaks = c("Observed", "Fitted with frailty", "Fitted without frailty"),
      name   = NULL
    ) +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis(~ . / (N_ref * bw), name = "Density")
    ) +
    labs(x = "Interaction") +
    theme_minimal(base_size = 20) +
    theme(
      # 背景网格虚线
      panel.grid.major = element_line(linetype = "dashed"),
      panel.grid.minor = element_blank(),
      
      axis.title.x         = element_text(size = 20),
      axis.title.y         = element_text(size = 20),
      axis.title.y.right   = element_text(size = 20, margin = margin(l = 8)),
      axis.text.x          = element_text(size = 18),
      axis.text.y          = element_text(size = 18),
      legend.text          = element_text(size = 18),
      legend.key.width     = unit(2.0, "lines"),
      legend.key.height    = unit(1.6, "lines"),
      legend.position      = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background    = element_rect(fill = "white", colour = NA)
    )+
    coord_cartesian(
      xlim = c(0, 950),
      ylim = c(0, y_max * 1.05)
    )
  
  p
  
}

wh <- dev.size("in")   # wh[1] = width, wh[2] = height

ggsave(
  filename = "interaction_hist_density.tiff",
  plot     = p,
  device   = "tiff",
  width    = wh[1],
  height   = wh[2],
  units    = "in",
  dpi      = 600,
  compression = "lzw"
)

library(writexl)  # 如果没装过：install.packages("writexl")

df_all <- data.frame(
  obs        = df_obs$value,
  with_frailty  = df_with$value,
  without_frailty = df_wo$value
)

# 2. 写入 Excel（当前工作目录下生成 lambda_compare.xlsx）
write_xlsx(df_all, "three_vectors.xlsx")

####做分箱和QQ图
{
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # ---------------- 1) 准备数据 ----------------
  y       <- as.numeric(df_obs$value)
  mu_fix  <- as.numeric(df_wo$value)     # 基线（无 frailty）的条件均值
  stopifnot(length(y) == length(mu_fix))
  keep <- is.finite(y) & is.finite(mu_fix) & mu_fix >= 0
  y <- y[keep]; mu_fix <- mu_fix[keep]
  
  # 总量配准：避免“整体偏移”影响校准
  rescale_to_total <- function(mu, y) { s <- sum(y)/max(sum(mu), .Machine$double.eps); mu * s }
  mu_fix <- rescale_to_total(mu_fix, y)
  
  # 估计/指定 frailty 方差（对数尺度的标准差 sigma）；你也可以用你文中的估计：sigma = 0.917
  sigma <- 0.917
  
  # ---------------- 2) 后验预测模拟 ----------------
  # A) 无 frailty：Poisson(mu_fix)
  # B) 有 frailty：delta ~ N(-sigma^2/2, sigma^2)，lambda = mu_fix * exp(delta)，Y ~ Poisson(lambda)
  nsim <- 500
  set.seed(1)
  
  sim_pois <- matrix(rpois(length(y)*nsim, lambda = rep(mu_fix, nsim)), nrow = length(y), ncol = nsim)
  
  delta_mat <- matrix(rnorm(length(y)*nsim, mean = -0.5*sigma^2, sd = sigma), nrow = length(y), ncol = nsim)
  lambda_fr <- mu_fix * exp(delta_mat)
  sim_frlty <- matrix(rpois(length(y)*nsim, lambda = lambda_fr), nrow = length(y), ncol = nsim)
  
  # ---------------- 3) 图A：均值–方差（按共同分箱） ----------------
  n_bins <- 20
  rank_var <- mu_fix                       # 共同排序变量
  bin_id <- as.integer(cut(rank_var, breaks = quantile(rank_var, probs = seq(0,1,length.out = n_bins+1),
                                                       na.rm = TRUE), include.lowest = TRUE, labels = FALSE))
  
  agg <- function(v) tapply(v, bin_id, function(z) c(mean = mean(z), var = var(z)))
  obs_mv <- do.call(rbind, agg(y))                       # 观测 mean/var
  fix_mv <- cbind(mean = tapply(mu_fix, bin_id, mean),   # Poisson 线：var = mean
                  var  = tapply(mu_fix, bin_id, mean))
  # frailty 预测方差：对每个 bin、对 nsim 次模拟求方差的均值（即 E[var] 的蒙特卡洛估计）
  fr_var <- tapply(seq_len(length(y)), bin_id, function(idx){
    # 每次模拟在该 bin 的方差
    apply(sim_frlty[idx, , drop = FALSE], 2, var) |> mean()
  })
  fr_mean <- tapply(mu_fix, bin_id, mean)
  fr_mv <- cbind(mean = fr_mean, var = fr_var)
  
  mv_df <- tibble(
    bin = seq_len(n_bins),
    obs_mean = obs_mv[,"mean"], obs_var = obs_mv[,"var"],
    poi_mean = fix_mv[,"mean"], poi_var = fix_mv[,"var"],
    fr_mean  = fr_mv[,"mean"],  fr_var  = fr_mv[,"var"]
  )
  
  pA <- ggplot(mv_df, aes(x = obs_mean, y = obs_var)) +
    geom_point(size = 2.2, color = "#374151") +
    geom_line(aes(x = poi_mean, y = poi_var), color = "#2563EB", linewidth = 1, linetype = 2) +
    geom_line(aes(x = fr_mean,  y = fr_var),  color = "#D1495B", linewidth = 1.1) +
    labs(x = "Bin mean of counts", y = "Bin variance of counts",
         title = "Mean–variance by predicted bins",
         subtitle = "Dots: observed; blue dashed: Poisson Var=Mean; red: frailty posterior-predictive Var") +
    theme_minimal(base_size = 13)
  
  # ---------------- 4) 图B：上尾覆盖（tail coverage） ----------------
  qs <- c(0.80, 0.90, 0.95, 0.99)
  # 每条边的预测分位（Poisson / Frailty），然后看观测是否超过
  q_pois  <- sapply(qs, function(q) apply(sim_pois, 1, quantile, probs = q))
  q_frlty <- sapply(qs, function(q) apply(sim_frlty, 1, quantile, probs = q))
  
  tail_cov <- tibble(
    q = rep(qs, times = 2),
    exceed = c(colMeans(y > q_pois), colMeans(y > q_frlty)),
    model  = rep(c("No frailty (Poisson)", "With frailty"), each = length(qs)),
    target = 1 - rep(qs, times = 2)
  )
  
  pB <- ggplot(tail_cov, aes(x = q, y = exceed, color = model)) +
    geom_hline(aes(yintercept = target), linetype = 3, color = "grey50") +
    geom_line(linewidth = 1.1) + geom_point(size = 2) +
    scale_x_continuous(breaks = qs, labels = paste0(qs*100, "%")) +
    labs(x = "Predictive quantile (q)", y = "Pr{ Y_obs > predictive q }",
         title = "Upper-tail coverage",
         subtitle = "Ideal line: 1 - q (grey). Frailty should be closer to target, esp. at 95–99%") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
  
  # ---------------- 5) 拼图 ----------------
  pA + pB + plot_layout(widths = c(1.1, 0.9))
  
}

##########============空间热力分布
{
  library(dplyr)
  library(ggplot2)
  
  # --- 参数 ---
  eps <- 1e-3
  bins <- 22
  spread <- 1.6
  pt_size <- 3.2
  
  # --- 预处理：N/C/F + 组成占比 ---
  df <- cov_data %>%
    transmute(
      group1 = case_when(
        group_period1 == "N" ~ "Nurse",
        group_period1 == "C" ~ "Cleaner",
        group_period1 == "F" ~ "Forage",
        TRUE ~ NA_character_
      ),
      b = as.numeric(visits_to_brood),
      e = as.numeric(visits_to_nest_entrance),
      r = as.numeric(visits_to_rubbishpile)
    ) %>%
    filter(!is.na(group1)) %>%
    mutate(total = b + e + r) %>%
    filter(total > 0) %>%
    mutate(
      b = (b + eps) / (total + 3*eps),
      e = (e + eps) / (total + 3*eps),
      r = (r + eps) / (total + 3*eps)
    ) %>%
    mutate(
      u1 = sqrt(1/2) * log(b / e),
      u2 = sqrt(1/6) * log((b * e) / (r^2)),
      x  = plogis(scale(u1) * spread),
      y  = plogis(scale(u2) * spread)
    )
  
  # --- 颜色（热力与点）---
  heat_cols <- c("#FFFDFE","#EFF5FF","#DCEBFF","#C6D9FF","#A5BFFF",
                 "#7EA0FF","#587EFF","#365AE3","#2537A6")
  
  p <- ggplot(df, aes(x, y)) +
    # 浅色背景
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "#FAFAFF", color = NA) +
    # 关键修复：使用 level（不是 ndensity）
    stat_density_2d_filled(
      aes(fill = after_stat(level)),
      bins = bins, contour_var = "ndensity", alpha = 0.88,
      show.legend = FALSE
    ) +
    scale_fill_gradientn(colours = heat_cols) +
    # 更大更深的点 + 仅保留 Period 1 图例
    geom_point(aes(color = group1), size = pt_size, alpha = 0.95) +
    scale_color_manual(
      values = c("Nurse"="#C2185B", "Cleaner"="#1D4ED8", "Forage"="#6D28D9"),
      name = "Period 1"
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range,
                    expand = expansion(mult = 0.02))+
    theme_void(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "#FAFAFF", color = NA),
      legend.position  = "right",
      legend.title     = element_text(face = "bold")
    )
  
  p
  
}

{
  # install.packages(c("igraph","ggplot2","ggrepel"))  # 如未安装
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  
  # pairs: 两列 (id1,id2)，与 lambda_vec 一一对应；节点编号从 1 开始
  # lambda_vec: 边强度（越大=越近）
  plot_nodes_2d_from_lambda <- function(pairs, lambda_vec, eps = 1e-9,
                                        alpha = 1,   # 距离变换：length = 1/(lambda^alpha + eps)，alpha<1可压缩动态范围
                                        label = TRUE, seed = 1L) {
    stopifnot(nrow(pairs) == length(lambda_vec))
    set.seed(seed)
    
    df <- data.frame(from = as.integer(pairs[,1]),
                     to   = as.integer(pairs[,2]),
                     lam  = as.numeric(lambda_vec))
    
    n <- max(df$from, df$to, na.rm = TRUE)
    g <- graph_from_data_frame(df, directed = FALSE,
                               vertices = data.frame(name = as.character(1:n)))
    # 把“越大越近”转成边长
    E(g)$length <- 1 / ( (E(g)$lam^alpha) + eps )
    
    memb <- components(g)$membership
    comps <- sort(unique(memb))
    
    coords <- matrix(NA_real_, n, 2)
    x_offset <- 0
    gap <- 3  # 分量之间的间隔
    
    for (cc in comps) {
      vids <- which(memb == cc)
      sg   <- induced_subgraph(g, vids = vids)
      
      if (vcount(sg) == 1L) {
        local <- matrix(c(0, 0), 1, 2)
      } else {
        D <- distances(sg, weights = E(sg)$length)
        # 经典 MDS（只要距离矩阵，不要边）
        local <- cmdscale(D, k = 2)
        # 标准化方便排版
        local <- scale(local, center = TRUE, scale = FALSE)
        sc <- sd(local)
        if (!is.na(sc) && sc > 0) local <- local / sc
      }
      # 并排摆放不同连通分量
      local[,1] <- local[,1] + x_offset
      coords[vids,] <- local
      width <- if (nrow(local) == 1) 1 else diff(range(local[,1]))
      x_offset <- x_offset + width + gap
    }
    
    # 画点（可按分量上色，或去掉 fill=comp 就是纯黑白点）
    plot_df <- data.frame(node = 1:n,
                          x = coords[,1], y = coords[,2],
                          comp = factor(memb))
    p <- ggplot(plot_df, aes(x, y, fill = comp)) +
      geom_point(shape = 21, size = 2.4, color = "black", stroke = 0.25) +
      { if (label) geom_text_repel(aes(label = node), size = 3, max.overlaps = 100) } +
      coord_equal() + theme_void() + theme(legend.position = "none") +
      labs(title = "2D embedding by λ (larger λ ⇒ closer)")
    print(p)
    
    invisible(plot_df)
  }
  
  # 用法：
  # plot_nodes_2d_from_lambda(real_data$pairs, lambda_by_L)
  # 如果强度范围跨度太大、点团挤在一起，可试：alpha = 0.5
   plot_nodes_2d_from_lambda(real_data$pairs, lambda_by_Λ, alpha = 1)
  
}
# 保存为 Excel，并把表当作 Excel 表对象（便于筛选/样式）
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
openxlsx::write.xlsx(cov_data, file = "cov_data.xlsx",
                     sheetName = "cov_data", asTable = TRUE, overwrite = TRUE)




#delta的分布
kk <- real_data$df
{
  # install.packages(c("igraph","ggplot2","ggrepel"))  # 如未安装
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  
  # pairs: 两列 (id1,id2)，与 lambda_vec 一一对应；节点编号从 1 开始
  # lambda_vec: 边强度（越大=越近）
  plot_nodes_2d_from_lambda <- function(pairs, lambda_vec, eps = 1e-9,
                                        alpha = 1,   # 距离变换：length = 1/(lambda^alpha + eps)，alpha<1可压缩动态范围
                                        label = TRUE, seed = 1L) {
    stopifnot(nrow(pairs) == length(lambda_vec))
    set.seed(seed)
    
    df <- data.frame(from = as.integer(pairs[,1]),
                     to   = as.integer(pairs[,2]),
                     lam  = as.numeric(lambda_vec))
    
    n <- max(df$from, df$to, na.rm = TRUE)
    g <- graph_from_data_frame(df, directed = FALSE,
                               vertices = data.frame(name = as.character(1:n)))
    # 把“越大越近”转成边长
    E(g)$length <- 1 / ( (E(g)$lam^alpha) + eps )
    
    memb <- components(g)$membership
    comps <- sort(unique(memb))
    
    coords <- matrix(NA_real_, n, 2)
    x_offset <- 0
    gap <- 3  # 分量之间的间隔
    
    for (cc in comps) {
      vids <- which(memb == cc)
      sg   <- induced_subgraph(g, vids = vids)
      
      if (vcount(sg) == 1L) {
        local <- matrix(c(0, 0), 1, 2)
      } else {
        D <- distances(sg, weights = E(sg)$length)
        # 经典 MDS（只要距离矩阵，不要边）
        local <- cmdscale(D, k = 2)
        # 标准化方便排版
        local <- scale(local, center = TRUE, scale = FALSE)
        sc <- sd(local)
        if (!is.na(sc) && sc > 0) local <- local / sc
      }
      # 并排摆放不同连通分量
      local[,1] <- local[,1] + x_offset
      coords[vids,] <- local
      width <- if (nrow(local) == 1) 1 else diff(range(local[,1]))
      x_offset <- x_offset + width + gap
    }
    
    # 画点（可按分量上色，或去掉 fill=comp 就是纯黑白点）
    plot_df <- data.frame(node = 1:n,
                          x = coords[,1], y = coords[,2],
                          comp = factor(memb))
    p <- ggplot(plot_df, aes(x, y, fill = comp)) +
      geom_point(shape = 21, size = 2.4, color = "black", stroke = 0.25) +
      { if (label) geom_text_repel(aes(label = node), size = 3, max.overlaps = 100) } +
      coord_equal() + theme_void() + theme(legend.position = "none") +
      labs(title = "2D embedding by λ (larger λ ⇒ closer)")
    print(p)
    
    invisible(plot_df)
  }
  
  # 用法：
  # plot_nodes_2d_from_lambda(real_data$pairs, lambda_by_L)
  # 如果强度范围跨度太大、点团挤在一起，可试：alpha = 0.5
   plot_nodes_2d_from_lambda(real_data$pairs, lambda_by_Λ, alpha = 1)
  
}

##########每天平均每对蚂蚁发生的次数
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  ## ---- 口径A：把 NA 视为“未观测”，计算当日所有已观测蚂蚁对的平均次数 ----
  avg_by_day <- time_data %>%
    summarise(across(matches("^Day\\d+$"), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Day", values_to = "avg_per_pair") %>%
    mutate(day = as.integer(sub("Day", "", Day))) %>%
    arrange(day)
  
  ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
    geom_line() +
    geom_point() +
    labs(title = "每日每对蚂蚁的平均事件次数（忽略 NA）",
         x = "Day", y = "Average events per pair") +
    theme_minimal(base_size = 13)
  
}


###########根据lambda_by_Λ画出网络的直观图
{
  
  #####先提取出身份
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # --- 1) 预处理：识别四个时期的身份列，并清洗为 N/C/F/NA ---
  grp_cols <- grep("^group_period\\d+$", names(cov_data), value = TRUE)
  stopifnot(length(grp_cols) >= 1)
  
  sanitize_group <- function(x) {
    x <- toupper(as.character(x))
    ifelse(x %in% c("N", "C", "F"), x, NA_character_)
  }
  
  # 长表：每个 tag_id × period 一行
  cov_long <- cov_data %>%
    mutate(tag_id_old = as.integer(tag_id)) %>%
    select(tag_id_old, all_of(grp_cols)) %>%
    pivot_longer(cols = all_of(grp_cols),
                 names_to = "period_col", values_to = "group_raw") %>%
    mutate(
      p = as.integer(str_extract(period_col, "\\d+")),
      group = sanitize_group(group_raw)
    ) %>%
    select(tag_id_old, p, group)
  
  # --- 2) 边表：把新编号 (node) 映射回原始编号 (old_node) ---
  edges <- as.data.frame(real_data$pairs) %>%
    as_tibble() %>%
    rename(node1 = id1, node2 = id2) %>%
    mutate(edge_index = row_number()) %>%
    # node_map: old_node -> node（新编号）
    left_join(node_map %>% transmute(node = as.integer(node),
                                     old_node1 = as.integer(old_node)),
              by = join_by(node1 == node)) %>%
    left_join(node_map %>% transmute(node = as.integer(node),
                                     old_node2 = as.integer(old_node)),
              by = join_by(node2 == node)) %>%
    # 附上每条边的 lambda 值（按 pairs 的顺序）
    mutate(lambda = `lambda_by_Λ`[edge_index]) %>%
    relocate(edge_index, node1, node2, old_node1, old_node2, lambda)
  
  # --- 3) 把双方在四个时期的身份并到每条边上 ---
  # 先把时期 p = 1..4 都铺开，再分别连到端点1和端点2
  ps <- sort(unique(cov_long$p))
  edge_long <- edges %>%
    tidyr::crossing(p = ps) %>%
    # 端点1身份
    left_join(cov_long %>% rename(group1 = group),
              by = c("old_node1" = "tag_id_old", "p")) %>%
    # 端点2身份
    left_join(cov_long %>% rename(group2 = group),
              by = c("old_node2" = "tag_id_old", "p"))
  
  # 宽表：group1_p1..p4, group2_p1..p4
  edge_groups <- edge_long %>%
    mutate(p = as.integer(p)) %>%
    pivot_wider(
      names_from = p,
      values_from = c(group1, group2),
      names_glue = "{.value}_p{p}"
    ) %>%
    arrange(edge_index)
  
  # 看看结果
  print(edge_groups, n = 10)
  
  ###做节点嵌入
  # ====== Packages ======
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(RSpectra)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  
  # ====== 0) 小工具：把身份标准化为 N/C/F/NA ======
  sanitize_group <- function(x) {
    x <- toupper(as.character(x))
    ifelse(x %in% c("N","C","F"), x, NA_character_)
  }
  
  # ====== 1) 边与权：由 pairs + lambda_by_Λ 生成加权无向图lambda_by_cov ======
  edges <- as_tibble(real_data$pairs) %>%
    rename(i = id1, j = id2) %>%
    mutate(w = `lambda_by_Λ`[row_number()])
  
  n_nodes <- max(node_map$node)               # 新编号的最大值就是节点总数
  W <- sparseMatrix(i = edges$i, j = edges$j, x = edges$w,
                    dims = c(n_nodes, n_nodes))
  W <- W + t(W)                               # 无向
  diag(W) <- 0
  
  deg <- Matrix::rowSums(W)                   # 度（权和）
  keep <- which(deg > 0)                      # 去掉孤立点，避免数值问题
  if (length(keep) < n_nodes) {
    message("有孤立节点（度=0），仅对有边的节点做嵌入；孤立点坐标将为 NA。")
  }
  
  # ====== 2) 图嵌入：拉普拉斯特征嵌入到二维 ======
  # 归一化拉普拉斯 L_sym = I - D^{-1/2} W D^{-1/2}
  Wk <- W[keep, keep, drop = FALSE]
  deg_k <- Matrix::rowSums(Wk)
  Dinv2 <- Diagonal(x = 1/sqrt(pmax(as.numeric(deg_k), 1e-12)))
  Lsym  <- Diagonal(n = length(keep)) - Dinv2 %*% Wk %*% Dinv2
  
  # 取最小的 3 个特征对（跳过最小的第1个常量特征，取第2/3作坐标）
  ev <- RSpectra::eigs_sym(Lsym, k = 3, which = "SM")
  Zk <- as.matrix(ev$vectors[, 2:3, drop = FALSE])
  colnames(Zk) <- c("x","y")
  
  # 把坐标放回完整节点集合（孤立点设 NA）
  coords <- tibble(node = keep, x = Zk[,1], y = Zk[,2]) %>%
    right_join(tibble(node = 1:n_nodes), by = "node") %>%
    arrange(node)
  
  # ====== 3) 节点 -> 原始编号（tag_id），拼上 4 期身份 ======
  coords2 <- coords %>%
    left_join(node_map %>% rename(tag_id = old_node), by = "node")
  
  grp_cols <- grep("^group_period\\d+$", names(cov_data), value = TRUE)
  stopifnot(length(grp_cols) >= 1)
  
  cov_groups <- cov_data %>%
    mutate(tag_id = as.integer(tag_id)) %>%
    select(tag_id, all_of(grp_cols)) %>%
    mutate(across(all_of(grp_cols), sanitize_group))
  
  # ====== 4) 生成四张图：每期只保留 group 非 NA 的个体 ======
  # ====== 4) 改版：排除 380，并统一四图坐标轴范围 ======
  exclude_ids <- 380L   # 如需排除多个：c(380, 999, ...)
  
  # 先把四期拼成长表
  cov_long <- cov_groups %>%
    tidyr::pivot_longer(all_of(grp_cols),
                        names_to = "period_col", values_to = "group") %>%
    mutate(p = as.integer(stringr::str_extract(period_col, "\\d+")))
  
  # 全时期联合数据，用于计算全局坐标范围
  df_all <- coords2 %>%
    dplyr::left_join(cov_long, by = "tag_id") %>%
    dplyr::filter(!is.na(group),                     # 该期有身份
                  !is.na(x), !is.na(y),              # 坐标存在（非孤立点）
                  !tag_id %in% exclude_ids)          # 排除离群编号
  
# 统一范围 + 留 2% 边距
xr <- range(df_all$x, na.rm = TRUE)
yr <- range(df_all$y, na.rm = TRUE)
pad <- 0.02
x_range_pad <- xr + c(-diff(xr)*pad,  diff(xr)*pad)
y_range_pad <- yr + c(-diff(yr)*pad,  diff(yr)*pad)

make_plot_for_period <- function(p_idx) {
  dfp <- df_all %>% dplyr::filter(p == p_idx)

  ggplot(dfp, aes(x = x, y = y, color = group)) +
    geom_point(size = 2.5, alpha = 0.95) +
    ggrepel::geom_text_repel(aes(label = tag_id),
                             size = 2.8, max.overlaps = 200,
                             show.legend = FALSE) +
    coord_cartesian(xlim = x_range_pad, ylim = y_range_pad, expand = FALSE) +
    labs(title = paste0("Period ", p_idx), x = NULL, y = NULL) +
    scale_color_manual(values = c(N = "#1f77b4", C = "#2ca02c", F = "#d62728")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

  plots <- lapply(seq_along(grp_cols), make_plot_for_period)
  final_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  final_plot
  
  
}

##########放在intro中的图，左边是统计交互次数，右边是时间
{
  # ##########4.看看所有的边交互情况，有没有交互很少的
  {
    
    {
      library(ggplot2)
      library(dplyr)
      library(MASS)   # glm.nb
      library(scales) # log1p_trans
      
      x <- low_interactions$total_interaction
      x <- x[is.finite(x)]
      n <- length(x)
      
      mu <- mean(x)
      v  <- var(x)
      phi <- v / mu
      p_chi <- pchisq((n - 1) * v / mu, df = n - 1, lower.tail = FALSE)  # 仅作描述性参考
      
      lab <- sprintf("Mean = %.4g\nVar = %.4g\nVar/Mean = %.3g\nPoisson λ̂ = %.4g\nChiSq p ≈ %.2g",
                     mu, v, phi, mu, p_chi)
      
      binwidth <- 10
      
      ## --- 用全数据拟合 NB（对比用；过离散时通常更贴合） ---
      df <- data.frame(x = x)
      fit_nb <- MASS::glm.nb(x ~ 1, data = df)
      theta_hat <- fit_nb$theta
      mu_nb <- as.numeric(exp(coef(fit_nb))[1])
      
      ## --- 构造 bin，并计算观测频数 O 与模型期望频数 E（Poisson / NB） ---
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
      
      ## ===================== 图 1：log1p y 轴下的“观测 vs 期望” =====================
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
        scale_y_continuous(trans = "log1p")  # ln(1+y)，底为 e
      p1
      
      {
        library(ggplot2)
        library(dplyr)
        
        plot_total_interaction_like <- function(time_data,
                                                binwidth = 10,
                                                x_max = 900,
                                                y_max = 1050) {
          # 1) 计算 total_interaction
          df <- time_data %>%
            mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
            filter(is.finite(total_interaction))
          
          # 2) 统计量
          mu <- mean(df$total_interaction)
          va <- var(df$total_interaction)
          
          lab <- sprintf("Mean = %.3f\n\nVariance = %.3f", mu, va)
          
          # 3) 作图（仿照示例样式）
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
              # 虚线网格（示例图的感觉）
              panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.5),
              panel.grid.minor = element_blank(),
              # 只保留左/下轴线（theme_classic 已经是这个效果）
              axis.line = element_line(color = "black", linewidth = 0.9),
              axis.ticks = element_line(color = "black"),
              axis.text = element_text(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.05)
            )
        }
        
        # 使用：
        p <- plot_total_interaction_like(time_data, binwidth = 10, x_max = 900, y_max = 1050)
        print(p)
        
      }
      
      
    }
    
    {
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      
      ## ---- 口径A：忽略 NA，算每一天的“已观测蚂蚁对”的平均事件次数 ----
      avg_by_day <- time_data %>%
        summarise(across(matches("^Day\\d+$"), ~ mean(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Day", values_to = "avg_per_pair") %>%
        mutate(day = as.integer(sub("Day", "", Day))) %>%
        arrange(day)
      
      # 1) 线性回归：y = avg_per_pair, x = day
      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      s <- summary(fit)
      
      slope <- coef(s)[ "day", "Estimate" ]
      pval  <- coef(s)[ "day", "Pr(>|t|)" ]
      
      lab <- sprintf("Slope = %.4g\np = %.4g", slope, pval)
      
      # 2) 作图：点+线+lm拟合线+右上角标注
      ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
        geom_line() +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.1, size = 4.5) +
        labs(title = "每日每对蚂蚁的平均事件次数（忽略 NA）",
             x = "Day", y = "Average events per pair") +
        theme_minimal(base_size = 13)
      fit <- lm(avg_per_pair ~ day, data = avg_by_day)
      
      ## 3) 汇总成三列数据框：day / observed / fitted
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
    
    # 图1：总交互次数直方图
    p1 <- ggplot(low_interactions, aes(x = total_interaction)) +
      geom_histogram(binwidth = 10, fill = "steelblue", color = "black",
                     boundary = 0, closed = "left") +
      labs(title = "Distribution of Total Interactions",
           x = "Total Interactions", y = "Number of Edges") +
      theme_minimal(base_size = 14)
    
    # 图2：每日每对蚂蚁的平均事件次数
    p2 <- ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
      geom_line() + geom_point() +
      labs(title = "Average daily number of interactions per ant pair
 ",
           x = "Day", y = "Number of interactons") +
      theme_minimal(base_size = 13)
    
    # 左右拼接
    (p1 | p2) +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(tag_levels = "A")  # 可选：给小图打 A/B 标签
    
  }
  
  
  
  #########模仿已经画好的
  {
    corner_anno <- function(label, size = 4.6, h = 1.02, v = 1.05) {
      annotate("text", x = Inf, y = Inf, label = label,
               hjust = h, vjust = v, size = size)
    }
    
    
    plot_ant_figure_like_nobreak <- function(time_data,
                                             binwidth = 10,
                                             x_max_hist = 900,
                                             y_max_hist = 1050,
                                             y_pad_ratio = 0.06,      # y轴上方留白比例
                                             widths = c(1.15, 1.00)) {
      
      suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(patchwork)
      })
      
      # -------------------------
      # Panel (a): Histogram
      # -------------------------
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
          plot.margin = margin(t = 8, r = 16, b = 6, l = 6)   # 关键：统一边距
        )
      
      
      # -------------------------
      # Panel (b): Avg per day (NO y break)
      # -------------------------
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
      
      # 自适应 y 轴范围 + 上方留白
      y_min <- min(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_max <- max(avg_by_day$avg_per_pair, na.rm = TRUE)
      y_pad <- (y_max - y_min) * y_pad_ratio
      if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2
      
      # 注释放在“实际坐标”避免任何分面/缩放带来的重复与漂移
      x_anno <- max(avg_by_day$day) - 1
      y_anno <- y_max + y_pad * 0.2
      
     p2 <- ggplot(avg_by_day, aes(x = day, y = avg_per_pair)) +
  geom_line(color = blue_pt, linewidth = 0.8) +
  geom_point(color = blue_pt, size = 2.6) +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linetype = "dotted", linewidth = 0.9) +
  corner_anno(lab_b) +   # 关键：与左图同一锚点
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
    plot.margin = margin(t = 8, r = 16, b = 6, l = 6)   # 关键：统一边距
  )

      # -------------------------
      # Combine
      # -------------------------
      fig <- p1 + p2 + plot_layout(ncol = 2, widths = widths)
      return(fig)
    }
    
    # 用法：
    fig <- plot_ant_figure_like_nobreak(time_data)
    print(fig)
    # ggsave("figure_like_nobreak.png", fig, width = 14, height = 4.8, dpi = 300)
    
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
                                             # 右上角注释放置位置（相对轴上限的比例，可微调）
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
      
      # 统一的“右上角附近 + 左对齐”注释
      corner_anno_left <- function(label, x_pos, y_pos) {
        annotate("text",
                 x = x_pos, y = y_pos, label = label,
                 hjust = 0, vjust = 1,                 # 左对齐 + 顶对齐
                 size = anno_size, family = font_family)
      }
      
      # =========================
      # (a) Histogram + Poisson/NB curves + sqrt y
      # =========================
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
      
      # NegBin 矩估计（稳健）
      if (is.finite(va) && va > mu && mu > 0) {
        size_hat <- mu^2 / (va - mu)
      } else {
        size_hat <- 1e8
      }
      mu_nb_hat <- max(mu, 0)
      
      # bins（只用于显示范围内的曲线）
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
      
      # y_max（只用 [0, x_max_hist] 内的直方图高度来估，避免 hist 报错）
      if (is.null(y_max_hist)) {
        x_clip <- x_all[x_all <= x_max_hist]
        h_counts <- hist(x_clip, breaks = breaks_disp, plot = FALSE, right = FALSE)$counts
        y_max_hist <- max(c(h_counts, curve_df$expected), na.rm = TRUE) * 1.06
      }
      
      # 注释位置：靠右上，且“左对齐”
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
      
      # =========================
      # (b) Average per day (no axis break)
      # =========================
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
      
      # 右图的 x/y 上限（用于注释定位）
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
      
      # =========================
      # Combine：左右各占一半
      # =========================
      p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
    }
    
    # 用法：
    fig <- plot_ant_figure_like_nobreak(time_data, font_scale = 1.8,anno_y_frac       = 1)
    print(fig)
    
  }
  #fig <- plot_ant_figure_like_nobreak(time_data)
  
  ggsave(
    filename = "Rplot.pdf",
    plot = fig,
    device = cairo_pdf,     # 可选，但更稳（尤其是字体）
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
      
      # ---- 数据：total interaction ----
      df_a <- time_data %>%
        mutate(total_interaction = rowSums(across(-1), na.rm = TRUE)) %>%
        filter(is.finite(total_interaction))
      
      x_all <- as.numeric(df_a$total_interaction)
      x_all <- x_all[is.finite(x_all) & x_all >= 0]
      n <- length(x_all)
      
      # ---- Poisson / NegBin 参数 ----
      mu <- mean(x_all)
      va <- var(x_all)
      
      lambda_hat <- max(mu, 0)
      
      # NegBin：矩估计（确保可用）
      if (is.finite(va) && va > mu && mu > 0) {
        size_hat <- mu^2 / (va - mu)
      } else {
        size_hat <- 1e8
      }
      mu_nb_hat <- max(mu, 0)
      
      # ---- 分箱（用于“每个bin的期望计数”曲线）----
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
      
      # ---- y 轴上限（只用显示范围内的数据估计，避免 hist 报错）----
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
        labs(x = "Interaction", y = "Count") +  # 不设 title
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
    
    # 用法示例：
    p_left <- plot_ant_left_only(time_data, binwidth = 10, x_max_hist = 900, font_scale = 2)
    print(p_left)
    ggsave(
      filename = "p_left.pdf",
      plot = p_left,
      device = cairo_pdf,     # 文字/线条更稳
      width = 11, height = 10, units = "in",
      dpi = 600,              # 对PDF影响不大，但写上无妨
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
  # 或者 PNG：
  ggsave("ant_fig.png", fig, width = 14, height = 5.2, units = "in", dpi = 300)
  
}


####挑选3个点，画出delta_{ij}的分布
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # ----------------------------------------
  # 1. 基础数据
  # ----------------------------------------
  delta <- as.numeric(fit_real$delta)
  pairs <- as.matrix(real_data$pairs)
  
  # 节点个数：默认用 pairs 里的最大编号
  n_node <- max(pairs)   # 如果你确定是 162，也可以直接写 n_node <- 162
  
  # ----------------------------------------
  # 2. 构造 δ 的对称邻接矩阵 Dmat
  #    行/列 = 节点，元素 = 对应边的 δ_ij
  # ----------------------------------------
  Dmat <- matrix(NA_real_, n_node, n_node)
  
  for (e in seq_along(delta)) {
    i <- pairs[e, 1]
    j <- pairs[e, 2]
    Dmat[i, j] <- delta[e]
    Dmat[j, i] <- delta[e]
  }
  
  # ----------------------------------------
  # 3. 每个端点的 incident δ_ij 平均值
  # ----------------------------------------
  node_mean <- apply(Dmat, 1, function(x) mean(x, na.rm = TRUE))
  
  # ----------------------------------------
  # 4. 选出 5 个分位点：0%, 25%, 50%, 75%, 100%
  #    注意：这里是按 node_mean 从小到大排序后，
  #          在排序序号上取 5 个等分位置
  # ----------------------------------------
  ord <- order(node_mean)   # 节点索引按均值从低到高排序
  k   <- length(ord)
  
  probs   <- c(0, 0.25, 0.5, 0.75, 1)
  pos     <- 1 + probs * (k - 1)        # 映射到 [1, k] 上的连续位置
  idx_pos <- round(pos)                 # 取整到最近的整数索引
  idx_pos <- pmin(pmax(idx_pos, 1), k)  # 边界保护
  
  id_very_low  <- ord[idx_pos[1]]  # 0% 分位
  id_low       <- ord[idx_pos[2]]  # 25% 分位
  id_mid       <- ord[idx_pos[3]]  # 50% 分位
  id_high      <- ord[idx_pos[4]]  # 75% 分位
  id_very_high <- ord[idx_pos[5]]  # 100% 分位
  
  c(id_very_low  = id_very_low,
    id_low       = id_low,
    id_mid       = id_mid,
    id_high      = id_high,
    id_very_high = id_very_high)
  # 可以先 print 看看是哪 5 个节点
  
  # ----------------------------------------
  # 5. 取出这五个端点对应的 (n_node - 1) 维向量（剔除自身）
  #    每个向量的第 j 个元素是该端点与节点 j 的 δ_ij
  # ----------------------------------------
  vec_very_high <- Dmat[id_very_high, -id_very_high]
  vec_high      <- Dmat[id_high,      -id_high]
  vec_mid       <- Dmat[id_mid,       -id_mid]
  vec_low       <- Dmat[id_low,       -id_low]
  vec_very_low  <- Dmat[id_very_low,  -id_very_low]
  
  # 如果你希望没有边的地方当作 0 而不是 NA，可以打开下面注释：
  # vec_very_high[is.na(vec_very_high)] <- 0
  # vec_high[is.na(vec_high)]           <- 0
  # vec_mid[is.na(vec_mid)]             <- 0
  # vec_low[is.na(vec_low)]             <- 0
  # vec_very_low[is.na(vec_very_low)]   <- 0
  
  # ----------------------------------------
  # 6. 整理为长格式，准备画图
  #    五组英文标签：
  #    "Very high vulnerability"
  #    "High vulnerability"
  #    "Medium vulnerability"
  #    "Low vulnerability"
  #    "Very low vulnerability"
  # ----------------------------------------
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
  
  # ----------------------------------------
  # 7. 配色：高脆弱性偏红，低脆弱性偏蓝/绿
  # ----------------------------------------
  cols <- c(
    "Very high frailty" = "#B2182B",  # 深红
    "High frailty"      = "#D55E00",  # 橙红
    "Medium frailty"    = "#CC79A7",  # 紫
    "Low frailty"       = "#0072B2",  # 蓝
    "Very low frailty"  = "#009E73"   # 绿
  )
  
  # ----------------------------------------
  # 8. 分布图（核密度）
  # ----------------------------------------
  # 在画图前加这一行
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
    theme_bw(base_size = 18) +  # 全局基准字号加大
    theme(
      text           = element_text(size = 18),  # 所有文字整体再放大一档
      legend.position = "top",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 16),
      axis.title      = element_text(size = 18),
      axis.text       = element_text(size = 16),
      panel.grid      = element_blank()
    )
  
  print(p_density)
  # ----------------------------------------
  # （可选）9. 箱线图版本：纵向比较 5 组
  # ----------------------------------------
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


# 2. 保存为 Excel（推荐 writexl 包）
# 安装一次：install.packages("writexl")
library(writexl)
write_xlsx(df_D, path = "D_vectors.xlsx")

head(real_data$pairs)
lambda_by_Λ

fit_real$beta_edge
colnames(scan$tim)
