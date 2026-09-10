// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

// -- 计算 a_row = exp(X %*% beta)，兼容 dgCMatrix / dense matrix
static arma::vec exp_X_beta(SEXP X_, const arma::vec& beta) {
  if (Rf_inherits(X_, "dgCMatrix")) {
    arma::sp_mat Xsp = Rcpp::as<arma::sp_mat>(X_);
    arma::vec eta = Xsp * beta;
    return arma::exp(eta);
  } else {
    arma::mat X = Rcpp::as<arma::mat>(X_);
    arma::vec eta = X * beta;
    return arma::exp(eta);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector approx_diag_K_cpp(Rcpp::List scan,
                                      Rcpp::NumericVector beta_,
                                      Rcpp::NumericVector delta_,
                                      std::string ties_method = "breslow") {
  // ---- 取出 scan 组件
  SEXP X_ = scan["X_rows"];                        // matrix or dgCMatrix
  IntegerVector id_row = scan["id_row"];           // length nR, 1..m
  NumericVector ss = scan["start_sorted"];         // length nR (sorted)
  NumericVector es = scan["stop_sorted"];          // length nR (sorted)
  IntegerVector rs = scan["rows_start"];           // length nR, 1..nR
  IntegerVector re = scan["rows_stop"];            // length nR, 1..nR
  NumericVector times = scan["times"];             // length K
  List grp = scan["events_rowidx_by_time"];        // length K; each IntegerVector of row ids (1..nR)
  
  // ---- 维度
  arma::vec beta = as<arma::vec>(beta_);
  arma::vec delta = as<arma::vec>(delta_);
  const int m  = delta.n_elem;
  
  // a_row = exp(X %*% beta)
  arma::vec a_row = exp_X_beta(X_, beta);
  
  // d_id = exp(delta)
  arma::vec d_id = arma::exp(delta);
  
  const int nR = a_row.n_elem;
  if (id_row.size() != nR || ss.size() != nR || es.size() != nR ||
      rs.size() != nR || re.size() != nR) {
    stop("scan fields have inconsistent lengths with nrow(X_rows).");
  }
  if (times.size() != grp.size()) {
    stop("length(times) must equal length(events_rowidx_by_time).");
  }
  
  // ---- 工作向量
  arma::vec Wid(m, arma::fill::zeros);   // ∑_{rows in risk set & id=i} a_row[r]*d_id[i]
  arma::vec diagK(m, arma::fill::zeros); // 输出
  double S0 = 0.0;                       // ∑_{rows in risk set} a_row[r]*d_id[id_row[r]]
  
  // 1-based -> 0-based 扫描指针
  int si = 0; // start index
  int ei = 0; // end index
  
  const bool use_breslow = (ties_method == "breslow");
  
  // ---- 主循环：按 time 逐点扫描
  const int K = times.size();
  for (int k = 0; k < K; ++k) {
    const double t = times[k];
    
    // 增加：所有 start_sorted <= t 的行进入风险集
    while (si < nR && ss[si] <= t) {
      int r = rs[si] - 1;                       // 1-based -> 0-based
      int eid = id_row[r] - 1;                  // 1-based -> 0-based
      double wr = a_row[r] * d_id[eid];
      S0 += wr;
      Wid[eid] += wr;
      ++si;
    }
    
    // 移除：所有 stop_sorted < t 的行离开风险集
    while (ei < nR && es[ei] < t) {
      int r = re[ei] - 1;                       // 1-based -> 0-based
      int eid = id_row[r] - 1;                  // 1-based -> 0-based
      double wr = a_row[r] * d_id[eid];
      S0 -= wr;
      Wid[eid] -= wr;
      ++ei;
    }
    
    if (!(S0 > 0.0) || !std::isfinite(S0)) continue;
    
    // 当前时点 k 的并列事件行号
    IntegerVector ev_rows = grp[k];
    const int mk = ev_rows.size();
    if (mk == 0) continue;
    
    if (use_breslow) {
      // w = Wid / S0; diagK += mk * (w * (1 - w))
      arma::vec w = Wid / S0;
      diagK += mk * (w % (1.0 - w));
    } else {
      // Efron：wevent[i] = 该时点属于 id=i 的事件行 wr_ev 之和
      arma::vec wevent(m, arma::fill::zeros);
      if (mk > 0) {
        for (int j = 0; j < mk; ++j) {
          int r = ev_rows[j] - 1;               // 1-based -> 0-based
          int eid = id_row[r] - 1;              // 1-based -> 0-based
          double wr_ev = a_row[r] * d_id[eid];
          wevent[eid] += wr_ev;
        }
      }
      const double Ek = arma::accu(wevent);
      
      if (!(Ek > 0.0) || !std::isfinite(Ek)) {
        arma::vec w = Wid / S0;
        diagK += mk * (w % (1.0 - w));
      } else {
        // for ell = 0 .. mk-1:
        for (int ell = 0; ell < mk; ++ell) {
          double alpha = static_cast<double>(ell) / static_cast<double>(mk);
          double denom = S0 - alpha * Ek;
          if (!(denom > 0.0) || !std::isfinite(denom)) continue;
          
          arma::vec w_ell = (Wid - alpha * wevent) / denom;
          diagK += (w_ell % (1.0 - w_ell));
        }
      }
    }
    
    // 可选：长循环里偶尔允许用户中断
    if ((k & 1023) == 0) Rcpp::checkUserInterrupt();
  }
  
  // pmax(diagK, 0)
  diagK.transform( [](double x) { return (x < 0.0) ? 0.0 : x; } );
  return wrap(diagK);
}
