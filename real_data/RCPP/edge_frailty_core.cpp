// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using arma::mat; using arma::vec;

inline double log_safe(double x){ return std::log(x); }

// ---- β block: value & grad wrt beta (Breslow/Efron) ----
// [[Rcpp::export]]
Rcpp::List cpp_value_grad_beta(
    const arma::vec& beta,
    const arma::vec& delta,
    const arma::mat& X,
    const IntegerVector& id,
    const NumericVector& start_sorted,
    const NumericVector& stop_sorted,
    const IntegerVector& rows_start,
    const IntegerVector& rows_stop,
    const NumericVector& times,
    const List& events_by_time,
    std::string ties_method = "efron"
){
  const int nR = X.n_rows, p = X.n_cols, K = times.size();
  const int m = delta.n_elem;
  
  arma::vec a_row = arma::exp(X * beta);
  arma::vec d_id  = arma::exp(delta);
  
  double S0 = 0.0, val = 0.0;
  arma::vec g_beta(p, arma::fill::zeros);
  arma::vec Xtw  (p, arma::fill::zeros);
  arma::vec Wid  (m, arma::fill::zeros);
  
  int si = 0, ei = 0;
  for (int k = 0; k < K; ++k){
    const double t = times[k];
    
    while (si < nR && start_sorted[si] <= t){
      int r = rows_start[si] - 1, eid = id[r] - 1;
      double wr = a_row[r] * d_id[eid];
      S0 += wr;  Xtw += wr * X.row(r).t();  Wid[eid] += wr;  ++si;
    }
    while (ei < nR && stop_sorted[ei] < t){
      int r = rows_stop[ei] - 1, eid = id[r] - 1;
      double wr = a_row[r] * d_id[eid];
      S0 -= wr;  Xtw -= wr * X.row(r).t();  Wid[eid] -= wr;  ++ei;
    }
    
    IntegerVector ev_rows = events_by_time[k];
    const int mk = ev_rows.size();
    if (mk == 0 || !(S0 > 0) || !R_finite(S0)) continue;
    
    // numerator
    for (int j = 0; j < mk; ++j){
      int r = ev_rows[j] - 1, eid = id[r] - 1;
      val += std::log(a_row[r]) + std::log(d_id[eid]);
    }
    arma::vec x_ev_sum(p, arma::fill::zeros);
    for (int j = 0; j < mk; ++j) x_ev_sum += X.row(ev_rows[j]-1).t();
    
    if (ties_method == "breslow"){
      val -= mk * log_safe(S0);
      g_beta += x_ev_sum - ( (double)mk * Xtw / S0 );
    } else { // efron
      arma::vec wr_ev(mk);
      for (int j = 0; j < mk; ++j){
        int r = ev_rows[j] - 1, eid = id[r] - 1;
        wr_ev[j] = a_row[r] * d_id[eid];
      }
      const double Ek = arma::accu(wr_ev);
      if (!(Ek > 0) || !R_finite(Ek)){
        val -= mk * log_safe(S0);
        g_beta += x_ev_sum - ( (double)mk * Xtw / S0 );
      } else {
        arma::vec Xev_w(p, arma::fill::zeros);
        for (int j = 0; j < mk; ++j){
          int r = ev_rows[j] - 1;
          Xev_w += wr_ev[j] * X.row(r).t();
        }
        for (int ell = 0; ell < mk; ++ell){
          double alpha = double(ell) / double(mk);
          double denom = S0 - alpha * Ek;
          if (!(denom > 0) || !R_finite(denom)) continue;
          val -= log_safe(denom);
          g_beta += -( (Xtw - alpha * Xev_w) / denom );
        }
        g_beta += x_ev_sum;
      }
    }
  }
  
  return List::create(_["value"]=val, _["grad"]=NumericVector(g_beta.begin(), g_beta.end()));
}

// ---- δ block: value & grad wrt delta (Breslow/Efron) ----
// [[Rcpp::export]]
Rcpp::List cpp_value_grad_delta(
    const arma::vec& delta,
    const arma::vec& beta,
    const arma::mat& X,
    const IntegerVector& id,
    const NumericVector& start_sorted,
    const NumericVector& stop_sorted,
    const IntegerVector& rows_start,
    const IntegerVector& rows_stop,
    const NumericVector& times,
    const List& events_by_time,
    std::string ties_method = "breslow"
){
  const int nR = X.n_rows, K = times.size();
  const int m  = delta.n_elem;
  
  arma::vec a_row = arma::exp(X * beta);
  arma::vec d_id  = arma::exp(delta);
  
  double S0 = 0.0, val = 0.0;
  arma::vec Wid(m, arma::fill::zeros);
  arma::vec g_delta(m, arma::fill::zeros);
  
  int si = 0, ei = 0;
  
  for (int k = 0; k < K; ++k){
    double t = times[k];
    
    while (si < nR && start_sorted[si] <= t){
      int r = rows_start[si] - 1, eid = id[r] - 1;
      double wr = a_row[r] * d_id[eid];
      S0 += wr;  Wid[eid] += wr;  ++si;
    }
    while (ei < nR && stop_sorted[ei] < t){
      int r = rows_stop[ei] - 1, eid = id[r] - 1;
      double wr = a_row[r] * d_id[eid];
      S0 -= wr;  Wid[eid] -= wr;  ++ei;
    }
    
    IntegerVector ev_rows = events_by_time[k];
    const int mk = ev_rows.size();
    if (mk == 0 || !(S0 > 0) || !R_finite(S0)) continue;
    
    // numerator
    for (int j = 0; j < mk; ++j){
      int r = ev_rows[j] - 1, eid = id[r] - 1;
      val += std::log(d_id[eid]);     // + sum xbeta handled in beta-block
      g_delta[eid] += 1.0;
    }
    
    if (ties_method == "breslow"){
      val -= mk * log_safe(S0);
      g_delta -= ( (double)mk * Wid / S0 );
    } else {
      // efron
      arma::vec wevent(m, arma::fill::zeros);
      double Ek = 0.0;
      for (int j = 0; j < mk; ++j){
        int r = ev_rows[j] - 1, eid = id[r] - 1;
        double wr = a_row[r] * d_id[eid];
        wevent[eid] += wr; Ek += wr;
      }
      if (!(Ek > 0) || !R_finite(Ek)){
        val -= mk * log_safe(S0);
        g_delta -= ( (double)mk * Wid / S0 );
      } else {
        for (int ell = 0; ell < mk; ++ell){
          double alpha = double(ell) / double(mk);
          double denom = S0 - alpha * Ek;
          if (!(denom > 0) || !R_finite(denom)) continue;
          val -= log_safe(denom);
          g_delta -= ( (Wid - alpha * wevent) / denom );
        }
      }
    }
  }
  
  return List::create(_["value"]=val, _["grad_delta"]=NumericVector(g_delta.begin(), g_delta.end()));
}

// ---- cumulative baseline hazard (Breslow/Efron) ----
// [[Rcpp::export]]
Rcpp::List cpp_compute_cumLambda(
    const arma::mat& X,
    const arma::vec& beta,
    const arma::vec& delta,
    const IntegerVector& id,
    const NumericVector& start_sorted,
    const NumericVector& stop_sorted,
    const IntegerVector& rows_start,
    const IntegerVector& rows_stop,
    const NumericVector& times,
    const List& events_by_time,
    std::string method = "breslow"
){
  const int nR = X.n_rows, K = times.size();
  const int m = delta.n_elem;
  
  arma::vec a_row = arma::exp(X * beta);
  arma::vec d_id  = arma::exp(delta);
  
  arma::vec Wid(m, arma::fill::zeros);
  arma::vec dLambda(K, arma::fill::zeros);
  arma::vec S0_by_time(K, arma::fill::zeros);
  arma::vec event_risk_by_time(K, arma::fill::zeros);
  
  double S0 = 0.0; int si = 0, ei = 0;
  
  for (int k = 0; k < K; ++k){
    double t = times[k];
    while (si < nR && start_sorted[si] <= t){
      int r = rows_start[si]-1, eid = id[r]-1;
      double wr = a_row[r] * d_id[eid];
      S0 += wr; Wid[eid] += wr; ++si;
    }
    while (ei < nR && stop_sorted[ei] < t){
      int r = rows_stop[ei]-1, eid = id[r]-1;
      double wr = a_row[r] * d_id[eid];
      S0 -= wr; Wid[eid] -= wr; ++ei;
    }
    
    S0_by_time[k] = S0;
    IntegerVector ev_rows = events_by_time[k];
    const int mk = ev_rows.size();
    if (!(S0 > 0) || !R_finite(S0) || mk == 0) continue;
    
    if (method == "efron"){
      double Ek = 0.0;
      for (int j = 0; j < mk; ++j){
        int r = ev_rows[j]-1, eid = id[r]-1;
        Ek += a_row[r] * d_id[eid];
      }
      event_risk_by_time[k] = Ek;
      for (int ell = 0; ell < mk; ++ell){
        double alpha = double(ell)/double(mk);
        dLambda[k] += 1.0/(S0 - alpha*Ek);
      }
    } else {
      dLambda[k] = double(mk)/S0;
      event_risk_by_time[k] = NA_REAL;
    }
  }
  
  arma::vec cumLambda = arma::cumsum(dLambda);
  return List::create(
    _["times"]=times,
    _["dLambda_by_time"]=NumericVector(dLambda.begin(), dLambda.end()),
    _["cumLambda_by_time"]=NumericVector(cumLambda.begin(), cumLambda.end()),
    _["S0_by_time"]=NumericVector(S0_by_time.begin(), S0_by_time.end()),
    _["event_risk_by_time"]=NumericVector(event_risk_by_time.begin(), event_risk_by_time.end())
  );
}

// ---- full diag(K) via Λ ----
// [[Rcpp::export]]
Rcpp::NumericVector cpp_K_diag_full(
    const arma::mat& X,
    const arma::vec& beta,
    const arma::vec& delta,
    const IntegerVector& id,
    const NumericVector& times,
    const IntegerVector& rows_start,
    const IntegerVector& rows_stop,
    const NumericVector& start_sorted,
    const NumericVector& stop_sorted,
    const NumericVector& cumLambda_by_time
){
  const int nR = X.n_rows, K = times.size();
  const int m = delta.n_elem;
  
  std::vector<double> start_by_row(nR), stop_by_row(nR);
  for (int i = 0; i < nR; ++i){
    start_by_row[ rows_start[i]-1 ] = start_sorted[i];
    stop_by_row [ rows_stop[i]-1 ]  = stop_sorted[i];
  }
  
  auto lower_count = [&](double s)->int{
    return std::lower_bound(times.begin(), times.end(), s) - times.begin();
  };
  auto upper_count = [&](double s)->int{
    return std::upper_bound(times.begin(), times.end(), s) - times.begin();
  };
  
  arma::vec xbeta = X * beta;
  arma::vec d_idv = arma::exp(delta);
  
  Rcpp::NumericVector diagK(m);
  for (int r = 0; r < nR; ++r){
    int k_start_left = lower_count(start_by_row[r]);
    int k_stop       = upper_count (stop_by_row[r]);
    
    double L_start = (k_start_left>0)? cumLambda_by_time[k_start_left-1] : 0.0;
    double L_stop  = (k_stop      >0)? cumLambda_by_time[k_stop-1]       : 0.0;
    double dL_row  = L_stop - L_start;
    if (dL_row <= 0) continue;
    
    int eid = id[r]-1;
    double exp_eta = std::exp(xbeta[r]) * d_idv[eid];
    diagK[eid] += dL_row * exp_eta;
  }
  for (int i = 0; i < m; ++i) if (diagK[i] < 0) diagK[i] = 0.0;
  return diagK;
}

// ---- Hessian-vector product K*x for δ ----
// [[Rcpp::export]]
Rcpp::NumericVector cpp_K_matvec(
    const arma::vec& x,
    const arma::vec& beta,
    const arma::vec& delta,
    const arma::mat& X,
    const IntegerVector& id,
    const NumericVector& start_sorted,
    const NumericVector& stop_sorted,
    const IntegerVector& rows_start,
    const IntegerVector& rows_stop,
    const NumericVector& times,
    const List& events_by_time,
    std::string ties_method = "breslow"
){
  const int nR = X.n_rows, K = times.size();
  const int m = delta.n_elem;
  
  arma::vec a_row = arma::exp(X * beta);
  arma::vec d_idv = arma::exp(delta);
  
  arma::vec out(m, arma::fill::zeros);
  arma::vec Wid(m, arma::fill::zeros);
  double S0 = 0.0; int si = 0, ei = 0;
  
  for (int k = 0; k < K; ++k){
    double t = times[k];
    while (si < nR && start_sorted[si] <= t){
      int r = rows_start[si]-1, eid = id[r]-1;
      double wr = a_row[r] * d_idv[eid];
      S0 += wr;  Wid[eid] += wr;  ++si;
    }
    while (ei < nR && stop_sorted[ei] < t){
      int r = rows_stop[ei]-1, eid = id[r]-1;
      double wr = a_row[r] * d_idv[eid];
      S0 -= wr;  Wid[eid] -= wr;  ++ei;
    }
    
    IntegerVector ev_rows = events_by_time[k];
    const int mk = ev_rows.size();
    if (!(S0 > 0) || !R_finite(S0) || mk == 0) continue;
    
    if (ties_method == "breslow"){
      arma::vec w = Wid / S0;
      double dot = arma::dot(w, x);
      out += (double)mk * w % (x - dot);
    } else {
      arma::vec wevent(m, arma::fill::zeros);
      double Ek = 0.0;
      for (int j = 0; j < mk; ++j){
        int r = ev_rows[j]-1, eid = id[r]-1;
        double wr = a_row[r] * d_idv[eid];
        wevent[eid] += wr;  Ek += wr;
      }
      if (!(Ek > 0) || !R_finite(Ek)){
        arma::vec w = Wid / S0;
        double dot = arma::dot(w, x);
        out += (double)mk * w % (x - dot);
      } else {
        for (int ell = 0; ell < mk; ++ell){
          double alpha = double(ell)/double(mk);
          double denom = S0 - alpha * Ek;
          if (!(denom > 0) || !R_finite(denom)) continue;
          arma::vec w_ell = (Wid - alpha * wevent) / denom;
          double dot = arma::dot(w_ell, x);
          out += w_ell % (x - dot);
        }
      }
    }
  }
  return Rcpp::NumericVector(out.begin(), out.end());
}

// [[Rcpp::export]]
NumericVector compute_K_diag_partial_cpp(const List& scan,
                                         const arma::vec& beta,
                                         const arma::vec& delta,
                                         const std::string& ties_method = "efron") {
  // 1) 取出基本对象 -------------------------------------------------
  arma::mat X = as<arma::mat>(scan["X_rows"]);        // nR × p
  IntegerVector id_r = scan["id_row"];                // 行 -> ID (1-based)
  NumericVector tt = scan["times"];                   // K
  NumericVector start_sorted = scan["start_sorted"];  // nR
  NumericVector stop_sorted  = scan["stop_sorted"];   // nR
  IntegerVector rows_start   = scan["rows_start"];    // nR, 1-based row index
  IntegerVector rows_stop    = scan["rows_stop"];     // nR, 1-based row index
  List events_by_time        = scan["events_rowidx_by_time"]; // length K, each is Int vec (1-based)
  
  const int nR = X.n_rows;
  const int m  = delta.n_elem;                        // number of IDs (edges)
  const int K  = tt.size();
  
  // 2) 行级 exp(X beta) --------------------------------------------
  arma::vec a_row = exp(X * beta);                    // length nR
  
  // 3) ID 级 exp(delta) --------------------------------------------
  arma::vec d_id = exp(delta);                        // length m
  
  // 4) 扫描状态量 --------------------------------------------------
  std::vector<double> Wid(m, 0.0);    // 风险集内按ID聚合的权重
  double S0 = 0.0;                    // 风险集总和
  NumericVector diagK(m);             // 输出
  
  int si = 0; // index into start_sorted / rows_start (0-based in C++)
  int ei = 0; // index into stop_sorted / rows_stop
  
  // 5) 主时间循环 --------------------------------------------------
  for (int k = 0; k < K; ++k) {
    double t = tt[k];
    
    // ---- 有行进入风险集：start < t ----
    while (si < nR && start_sorted[si] < t) {
      int r   = rows_start[si] - 1;   // rows_start 是 1-based
      int eid = id_r[r] - 1;          // id_r 是 1-based
      double wr = a_row[r] * d_id[eid];
      
      S0        += wr;
      Wid[eid]  += wr;
      ++si;
    }
    
    // ---- 有行离开风险集：stop < t ----
    while (ei < nR && stop_sorted[ei] < t) {
      int r   = rows_stop[ei] - 1;
      int eid = id_r[r] - 1;
      double wr = a_row[r] * d_id[eid];
      
      S0        -= wr;
      Wid[eid]  -= wr;
      ++ei;
    }
    
    // 该时刻发生的事件行
    IntegerVector ev_rows = events_by_time[k];
    int mk = ev_rows.size();          // 这个时刻的事件数
    
    if (mk == 0 || !(S0 > 0.0) || !R_finite(S0)) {
      continue;
    }
    
    // =============== Breslow ===============
    if (ties_method == "breslow") {
      for (int j = 0; j < m; ++j) {
        double w = Wid[j] / S0;
        diagK[j] += mk * (w * (1.0 - w));
      }
    } else {
      // =============== Efron ===============
      // 计算事件部分的 wevent (按ID聚合)
      std::vector<double> wevent(m, 0.0);
      double Ek = 0.0;
      for (int ii = 0; ii < mk; ++ii) {
        int r  = ev_rows[ii] - 1;     // 事件行 1-based
        int id = id_r[r] - 1;
        double wri = a_row[r] * d_id[id];
        wevent[id] += wri;
        Ek += wri;
      }
      
      if (!(Ek > 0.0) || !R_finite(Ek)) {
        // 退化为 Breslow
        for (int j = 0; j < m; ++j) {
          double w = Wid[j] / S0;
          diagK[j] += mk * (w * (1.0 - w));
        }
      } else {
        // 正常的 Efron：ell = 0..mk-1
        for (int ell = 0; ell < mk; ++ell) {
          double alpha = static_cast<double>(ell) / static_cast<double>(mk);
          double denom = S0 - alpha * Ek;
          if (!(denom > 0.0) || !R_finite(denom)) continue;
          
          for (int j = 0; j < m; ++j) {
            double w_ell = (Wid[j] - alpha * wevent[j]) / denom;
            diagK[j] += (w_ell * (1.0 - w_ell));
          }
        }
      }
    }
  }
  
  // 6) 数值防护：负数截断为 0 --------------------------------------
  for (int j = 0; j < m; ++j) {
    if (diagK[j] < 0.0) diagK[j] = 0.0;
  }
  
  return diagK;
}
