// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
#include <iomanip>   // 进度显示需要

static inline int compute_next_report(int K, int progress_every){
  if (progress_every > 0) return std::max(1, progress_every);
  return std::max(1, K / 100); // 默认1%刷新一次
}

static inline void maybe_progress(const char* label, int k, int K, int next_report){
  
  if ((k % next_report) == 0 || (k + 1) == K){
    double pct = 100.0 * (k + 1) / (double)K;
    Rcpp::Rcout << "\r[" << label << "] " << std::fixed << std::setprecision(1) << pct << "%";
    Rcpp::Rcout.flush();
  }
  // 也让用户可中断
  if ((k & 255) == 0) Rcpp::checkUserInterrupt();
}
// ---- 小工具：安全取 list 元素 ----
template <typename T>
T get_field(const List& L, const char* name) {
  if (!L.containsElementNamed(name))
    stop(std::string("scan list missing field: ") + name);
  return as<T>(L[name]);
}

// 预先声明内部核函数
static void scan_sweep_prepare(
    const List& scan,
    const NumericMatrix& X,           // nR x p
    const NumericVector& beta,        // p
    const NumericVector& delta,       // m
    NumericVector& a_row,             // nR
    NumericVector& d_id               // m
) {
  const int nR = X.nrow();
  const int p  = X.ncol();
  if (beta.size() != p) stop("beta length mismatch with ncol(X)");
  // a_row = exp(X %*% beta)
  a_row = NumericVector(nR);
  for (int r = 0; r < nR; ++r) {
    double s = 0.0;
    for (int j = 0; j < p; ++j) s += X(r, j) * beta[j];
    a_row[r] = std::exp(s);
  }
  // d_id = exp(delta)
  const int m = delta.size();
  d_id = NumericVector(m);
  for (int e = 0; e < m; ++e) d_id[e] = std::exp(delta[e]);
}

// ---- 共用：一次“进入/离开”事件行的增量更新 Wid/S0/Wxe/WX_colsum ----
static inline void add_row(
    int r0,                          // 0-based row id
    const IntegerVector& id0,        // nR (1-based ids in R)
    const NumericVector& a_row,
    const NumericVector& d_id,
    const NumericMatrix& X,          // nR x p
    double sign,                     // +1 (enter) / -1 (leave)
    double& S0,
    NumericVector& Wid,              // m
    NumericMatrix& Wxe,              // m x p
    NumericVector& WX_colsum         // p
) {
  const int e = id0[r0] - 1;         // 0-based edge id
  if (e < 0) return; // 容错
  const double wr = a_row[r0] * d_id[e] * sign;
  S0 += wr;
  Wid[e] += wr;
  const int p = X.ncol();
  for (int j = 0; j < p; ++j) {
    const double inc = wr * X(r0, j);
    Wxe(e, j) += inc;
    WX_colsum[j] += inc;
  }
}

// ---- Htt: Breslow的 S2nd 累加 ----
static void Htt_add_breslow(
    const double mk,
    const double S0,
    const NumericVector& Wid,        // m
    const NumericMatrix& Wxe,        // m x p
    const NumericVector& WX_colsum,  // p
    NumericMatrix& Htt               // p x p, in-place +=
) {
  const int m = Wid.size();
  const int p = Htt.nrow();
  // barx = colSums(Wxe)/S0  (我们维护了 WX_colsum)
  NumericVector barx(p);
  for (int j = 0; j < p; ++j) barx[j] = WX_colsum[j] / S0;
  
  // S2nd = sum_e ( Wxe_e * Wxe_e^T / (S0 * Wid_e) )
  // 累加到 Htt:  Htt += mk*(S2nd - barx%*%t(barx))
  // 先做 Htt += mk*S2nd
  for (int e = 0; e < m; ++e) {
    const double we = Wid[e];
    if (!(we > 0.0) || !R_finite(we)) continue;
    const double coeff = mk / (S0 * we);
    for (int j = 0; j < p; ++j) {
      const double xej = Wxe(e, j);
      for (int l = 0; l < p; ++l) {
        Htt(j, l) += coeff * xej * Wxe(e, l);
      }
    }
  }
  // 再做 Htt -= mk*(barx outer barx)
  for (int j = 0; j < p; ++j) {
    const double bj = barx[j] * mk;
    for (int l = 0; l < p; ++l) {
      Htt(j, l) -= bj * barx[l];
    }
  }
}

// ---- Htt: Efron 单层 ell 的累加 ----
static void Htt_add_efron_layer(
    const double denom,              // S0 - alpha*Ek
    const double mk_or_1,            // Breslow是mk; Efron每层用 1
    const NumericVector& Wid_ell,    // m
    const NumericMatrix& Wxe_ell,    // m x p
    const NumericVector& WXcol_ell,  // p
    NumericMatrix& Htt               // p x p
) {
  const int m = Wid_ell.size();
  const int p = Htt.nrow();
  NumericVector barx(p);
  for (int j = 0; j < p; ++j) barx[j] = WXcol_ell[j] / denom;
  
  // Htt += (S2nd - barx barx^T)
  for (int e = 0; e < m; ++e) {
    const double we = Wid_ell[e];
    if (!(we > 0.0) || !R_finite(we)) continue;
    const double coeff = mk_or_1 / (denom * we);
    for (int j = 0; j < p; ++j) {
      const double xej = Wxe_ell(e, j);
      for (int l = 0; l < p; ++l) {
        Htt(j, l) += coeff * xej * Wxe_ell(e, l);
      }
    }
  }
  for (int j = 0; j < p; ++j) {
    const double bj = barx[j] * mk_or_1;
    for (int l = 0; l < p; ++l) {
      Htt(j, l) -= bj * barx[l];
    }
  }
}

// [[Rcpp::export]]
NumericMatrix Htt_matrix_cpp(
    List scan,
    NumericVector beta,
    NumericVector delta,
    std::string ties_method = "breslow"
) {
  NumericMatrix X   = get_field<NumericMatrix>(scan, "X_rows");      // nR x p
  IntegerVector id0 = get_field<IntegerVector>(scan, "id_row");      // nR, 1-based
  NumericVector ss  = get_field<NumericVector>(scan, "start_sorted");
  NumericVector es  = get_field<NumericVector>(scan, "stop_sorted");
  IntegerVector rs  = get_field<IntegerVector>(scan, "rows_start");  // 1-based
  IntegerVector re  = get_field<IntegerVector>(scan, "rows_stop");   // 1-based
  NumericVector times = get_field<NumericVector>(scan, "times");
  List groups = get_field<List>(scan, "events_rowidx_by_time");      // list of int[]
  
  const int nR = X.nrow();
  const int p  = X.ncol();
  const int m  = delta.size();
  const int K  = times.size();
  
  NumericVector a_row, d_id;
  scan_sweep_prepare(scan, X, beta, delta, a_row, d_id);
  
  NumericVector Wid(m);                 // m
  NumericMatrix Wxe(m, p);              // m x p
  NumericVector WX_colsum(p);           // p
  double S0 = 0.0;
  int si = 0, ei = 0;                   // 0-based pointers over ss/es/rs/re
  
  NumericMatrix Htt(p, p);              // output
  
  const bool is_breslow = (ties_method == "breslow");
  
  for (int k = 0; k < K; ++k) {
    maybe_progress("Htt", k, K, 100);
    const double t = times[k];
    while (si < nR && ss[si] < t) {
      const int r0 = rs[si] - 1;
      add_row(r0, id0, a_row, d_id, X, +1.0, S0, Wid, Wxe, WX_colsum);
      ++si;
    }
    while (ei < nR && es[ei] < t) {
      const int r0 = re[ei] - 1;
      add_row(r0, id0, a_row, d_id, X, -1.0, S0, Wid, Wxe, WX_colsum);
      ++ei;
    }
    // mk, events @ time k
    IntegerVector ev = groups[k];
    const int mk = ev.size();
    if (!(S0 > 0.0) || !R_finite(S0) || mk == 0) continue;
    
    if (is_breslow) {
      Htt_add_breslow(static_cast<double>(mk), S0, Wid, Wxe, WX_colsum, Htt);
    } else {
      // EFRON
      NumericVector wevent(m);              // Σ wr_ev per edge
      NumericMatrix WXevent(m, p);          // Σ wr_ev * x per edge
      NumericVector WXevent_colsum(p);      // col sums of WXevent
      
      if (mk > 0) {
        for (int ii = 0; ii < mk; ++ii) {
          const int r0 = ev[ii] - 1;
          const int e = id0[r0] - 1;
          const double wrev = a_row[r0] * d_id[e];
          wevent[e] += wrev;
          for (int j = 0; j < p; ++j) {
            const double inc = wrev * X(r0, j);
            WXevent(e, j) += inc;
            WXevent_colsum[j] += inc;
          }
        }
      }
      double Ek = 0.0;
      for (int e = 0; e < m; ++e) Ek += wevent[e];
      
      // 每一层 ell
      for (int ell = 0; ell < mk; ++ell) {
        const double alpha = static_cast<double>(ell) / static_cast<double>(mk);
        const double denom = S0 - alpha * Ek;
        if (!(denom > 0.0) || !R_finite(denom)) continue;
        
        // Wid_ell, Wxe_ell, WXcol_ell(=colSums Wxe_ell)
        NumericVector Wid_ell(m);
        NumericMatrix Wxe_ell(m, p);
        NumericVector WXcol_ell(p);
        
        for (int e = 0; e < m; ++e) {
          Wid_ell[e] = Wid[e] - alpha * wevent[e];
        }
        for (int e = 0; e < m; ++e) {
          for (int j = 0; j < p; ++j) {
            const double val = Wxe(e, j) - alpha * WXevent(e, j);
            Wxe_ell(e, j) = val;
            WXcol_ell[j] += val;
          }
        }
        
        Htt_add_efron_layer(denom, 1.0, Wid_ell, Wxe_ell, WXcol_ell, Htt);
      }
    }
  }
  return Htt;
}

// [[Rcpp::export]]
NumericMatrix Htd_matmul_cpp( // 输入 U: m×q, 输出: p×q
    List scan,
    NumericVector beta,
    NumericVector delta,
    NumericMatrix U,
    std::string ties_method = "breslow"
) {
  NumericMatrix X   = get_field<NumericMatrix>(scan, "X_rows");      // nR x p
  IntegerVector id0 = get_field<IntegerVector>(scan, "id_row");
  NumericVector ss  = get_field<NumericVector>(scan, "start_sorted");
  NumericVector es  = get_field<NumericVector>(scan, "stop_sorted");
  IntegerVector rs  = get_field<IntegerVector>(scan, "rows_start");
  IntegerVector re  = get_field<IntegerVector>(scan, "rows_stop");
  NumericVector times = get_field<NumericVector>(scan, "times");
  List groups = get_field<List>(scan, "events_rowidx_by_time");
  
  const int nR = X.nrow();
  const int p  = X.ncol();
  const int m  = delta.size();
  const int K  = times.size();
  const int q  = U.ncol();
  if (U.nrow() != m) stop("U (m x q): row mismatch with m");
  
  NumericVector a_row, d_id;
  scan_sweep_prepare(scan, X, beta, delta, a_row, d_id);
  
  NumericVector Wid(m);
  NumericMatrix Wxe(m, p);
  NumericVector WX_colsum(p);
  double S0 = 0.0;
  int si = 0, ei = 0;
  
  NumericMatrix Out(p, q);             // 累加输出
  
  const bool is_breslow = (ties_method == "breslow");
  
  for (int k = 0; k < K; ++k) {
    maybe_progress("Htd", k, K, 100);
    const double t = times[k];
    while (si < nR && ss[si] < t) { add_row(rs[si]-1, id0, a_row, d_id, X, +1.0, S0, Wid, Wxe, WX_colsum); ++si; }
    while (ei < nR && es[ei] < t) { add_row(re[ei]-1, id0, a_row, d_id, X, -1.0, S0, Wid, Wxe, WX_colsum); ++ei; }
    
    IntegerVector ev = groups[k];
    const int mk = ev.size();
    if (!(S0 > 0.0) || !R_finite(S0) || mk == 0) continue;
    
    if (is_breslow) {
      // barx
      NumericVector barx(p);
      for (int j = 0; j < p; ++j) barx[j] = WX_colsum[j] / S0;
      
      // s_vec = t(Wxe) %*% U / S0    (p×q)
      // s_scl = t(Wid)  %*% U / S0    (1×q)
      std::vector<double> s_scl(q, 0.0);
      NumericMatrix s_vec(p, q);
      for (int e = 0; e < m; ++e) {
        for (int jj = 0; jj < q; ++jj) s_scl[jj] += Wid[e] * U(e, jj);
        for (int j = 0; j < p; ++j) {
          const double xej = Wxe(e, j);
          for (int jj = 0; jj < q; ++jj) s_vec(j, jj) += xej * U(e, jj);
        }
      }
      for (int jj = 0; jj < q; ++jj) s_scl[jj] /= S0;
      for (int j = 0; j < p; ++j) for (int jj = 0; jj < q; ++jj) s_vec(j, jj) /= S0;
      
      // Out += mk * (s_vec - barx %*% s_scl)
      for (int j = 0; j < p; ++j) {
        const double bj = barx[j] * mk;
        for (int jj = 0; jj < q; ++jj) {
          Out(j, jj) += mk * s_vec(j, jj) - bj * s_scl[jj];
        }
      }
    } else {
      // EFRON
      NumericVector wevent(m);
      NumericMatrix WXevent(m, p);
      NumericVector WXevent_colsum(p);
      if (mk > 0) {
        for (int ii = 0; ii < mk; ++ii) {
          const int r0 = ev[ii] - 1;
          const int e  = id0[r0] - 1;
          const double wrev = a_row[r0] * d_id[e];
          wevent[e] += wrev;
          for (int j = 0; j < p; ++j) {
            const double inc = wrev * X(r0, j);
            WXevent(e, j) += inc;
            WXevent_colsum[j] += inc;
          }
        }
      }
      double Ek = 0.0;
      for (int e = 0; e < m; ++e) Ek += wevent[e];
      
      for (int ell = 0; ell < mk; ++ell) {
        const double alpha = static_cast<double>(ell) / static_cast<double>(mk);
        const double denom = S0 - alpha * Ek;
        if (!(denom > 0.0) || !R_finite(denom)) continue;
        
        // Wid_ell / Wxe_ell / WXcol_ell
        NumericVector Wid_ell(m);
        NumericVector WXcol_ell(p);
        for (int e = 0; e < m; ++e) Wid_ell[e] = Wid[e] - alpha * wevent[e];
        
        // s_vec = t(Wxe_ell) %*% U / denom,  s_scl = t(Wid_ell) %*% U / denom
        std::vector<double> s_scl(q, 0.0);
        NumericMatrix s_vec(p, q);
        
        for (int e = 0; e < m; ++e) {
          const double we = Wid[e] - alpha * wevent[e];
          for (int jj = 0; jj < q; ++jj) s_scl[jj] += we * U(e, jj);
        }
        for (int e = 0; e < m; ++e) {
          for (int j = 0; j < p; ++j) {
            const double val = Wxe(e, j) - alpha * WXevent(e, j);
            WXcol_ell[j] += val;
            for (int jj = 0; jj < q; ++jj) s_vec(j, jj) += val * U(e, jj);
          }
        }
        for (int jj = 0; jj < q; ++jj) s_scl[jj] /= denom;
        for (int j = 0; j < p; ++j) for (int jj = 0; jj < q; ++jj) s_vec(j, jj) /= denom;
        
        // Out += (s_vec - barx %*% s_scl)
        for (int j = 0; j < p; ++j) {
          const double bj = (WXcol_ell[j] / denom);
          for (int jj = 0; jj < q; ++jj) {
            Out(j, jj) += s_vec(j, jj) - bj * s_scl[jj];
          }
        }
      }
    }
  }
  return Out;
}

// [[Rcpp::export]]
NumericMatrix Hdt_matmul_cpp( // 输入 V: p×q, 输出: m×q
    List scan,
    NumericVector beta,
    NumericVector delta,
    NumericMatrix V,
    std::string ties_method = "breslow"
) {
  NumericMatrix X   = get_field<NumericMatrix>(scan, "X_rows");
  IntegerVector id0 = get_field<IntegerVector>(scan, "id_row");
  NumericVector ss  = get_field<NumericVector>(scan, "start_sorted");
  NumericVector es  = get_field<NumericVector>(scan, "stop_sorted");
  IntegerVector rs  = get_field<IntegerVector>(scan, "rows_start");
  IntegerVector re  = get_field<IntegerVector>(scan, "rows_stop");
  NumericVector times = get_field<NumericVector>(scan, "times");
  List groups = get_field<List>(scan, "events_rowidx_by_time");
  
  const int nR = X.nrow();
  const int p  = X.ncol();
  const int m  = delta.size();
  const int K  = times.size();
  const int q  = V.ncol();
  if (V.nrow() != p) stop("V (p x q): row mismatch with p");
  
  NumericVector a_row, d_id;
  scan_sweep_prepare(scan, X, beta, delta, a_row, d_id);
  
  NumericVector Wid(m);
  NumericMatrix Wxe(m, p);
  NumericVector WX_colsum(p);
  double S0 = 0.0;
  int si = 0, ei = 0;
  
  NumericMatrix Out(m, q);
  
  const bool is_breslow = (ties_method == "breslow");
  
  for (int k = 0; k < K; ++k) {
    maybe_progress("Hdt", k, K, 100);
    const double t = times[k];
    while (si < nR && ss[si] < t) { add_row(rs[si]-1, id0, a_row, d_id, X, +1.0, S0, Wid, Wxe, WX_colsum); ++si; }
    while (ei < nR && es[ei] < t) { add_row(re[ei]-1, id0, a_row, d_id, X, -1.0, S0, Wid, Wxe, WX_colsum); ++ei; }
    
    IntegerVector ev = groups[k];
    const int mk = ev.size();
    if (!(S0 > 0.0) || !R_finite(S0) || mk == 0) continue;
    
    if (is_breslow) {
      // barx, s1 = Wxe%*%V / S0, s2 = t(barx)%*%V
      NumericVector barx(p);
      for (int j = 0; j < p; ++j) barx[j] = WX_colsum[j] / S0;
      
      NumericMatrix s1(m, q);
      std::vector<double> s2(q, 0.0);
      
      // s1 = Wxe * V
      for (int e = 0; e < m; ++e) {
        for (int jj = 0; jj < q; ++jj) {
          double acc = 0.0;
          for (int j = 0; j < p; ++j) acc += Wxe(e, j) * V(j, jj);
          s1(e, jj) = acc / S0;
        }
      }
      for (int jj = 0; jj < q; ++jj) {
        double acc = 0.0;
        for (int j = 0; j < p; ++j) acc += barx[j] * V(j, jj);
        s2[jj] = acc;
      }
      for (int e = 0; e < m; ++e) {
        const double fac = static_cast<double>(mk) * (Wid[e] / S0);
        for (int jj = 0; jj < q; ++jj) {
          Out(e, jj) += static_cast<double>(mk) * s1(e, jj) - fac * s2[jj];
        }
      }
    } else {
      // EFRON
      NumericVector wevent(m);
      NumericMatrix WXevent(m, p);
      NumericVector WXevent_colsum(p);
      if (mk > 0) {
        for (int ii = 0; ii < mk; ++ii) {
          const int r0 = ev[ii] - 1;
          const int e  = id0[r0] - 1;
          const double wrev = a_row[r0] * d_id[e];
          wevent[e] += wrev;
          for (int j = 0; j < p; ++j) {
            const double inc = wrev * X(r0, j);
            WXevent(e, j) += inc;
            WXevent_colsum[j] += inc;
          }
        }
      }
      double Ek = 0.0;
      for (int e = 0; e < m; ++e) Ek += wevent[e];
      
      for (int ell = 0; ell < mk; ++ell) {
        const double alpha = static_cast<double>(ell) / static_cast<double>(mk);
        const double denom = S0 - alpha * Ek;
        if (!(denom > 0.0) || !R_finite(denom)) continue;
        
        // s1 = Wxe_ell %*% V / denom; s2 = t(barx_ell) %*% V
        NumericMatrix s1(m, q);
        std::vector<double> s2(q, 0.0);
        
        // barx_ell colsum:
        NumericVector WXcol_ell(p);
        NumericVector Wid_ell(m);
        for (int e = 0; e < m; ++e) Wid_ell[e] = Wid[e] - alpha * wevent[e];
        for (int e = 0; e < m; ++e) {
          for (int j = 0; j < p; ++j) {
            const double val = Wxe(e, j) - alpha * WXevent(e, j);
            WXcol_ell[j] += val;
          }
        }
        NumericVector barx_ell(p);
        for (int j = 0; j < p; ++j) barx_ell[j] = WXcol_ell[j] / denom;
        
        // s1
        for (int e = 0; e < m; ++e) {
          for (int jj = 0; jj < q; ++jj) {
            double acc = 0.0;
            for (int j = 0; j < p; ++j) {
              const double val = Wxe(e, j) - alpha * WXevent(e, j);
              acc += val * V(j, jj);
            }
            s1(e, jj) = acc / denom;
          }
        }
        for (int jj = 0; jj < q; ++jj) {
          double acc = 0.0;
          for (int j = 0; j < p; ++j) acc += barx_ell[j] * V(j, jj);
          s2[jj] = acc;
        }
        for (int e = 0; e < m; ++e) {
          const double fac = Wid_ell[e] / denom;
          for (int jj = 0; jj < q; ++jj) {
            Out(e, jj) += s1(e, jj) - fac * s2[jj];
          }
        }
      }
    }
  }
  return Out;
}

// ---------------------------------------------------------
// edge_joint_scores: 返回 U_theta_by_edge (m×p) 与 U_delta_by_edge (m)
// ---------------------------------------------------------
// [[Rcpp::export]]
List edge_joint_scores_cpp(
    List scan,
    NumericVector beta,
    NumericVector delta,
    std::string ties_method = "breslow"
) {
  NumericMatrix X   = get_field<NumericMatrix>(scan, "X_rows");
  IntegerVector id0 = get_field<IntegerVector>(scan, "id_row");
  NumericVector ss  = get_field<NumericVector>(scan, "start_sorted");
  NumericVector es  = get_field<NumericVector>(scan, "stop_sorted");
  IntegerVector rs  = get_field<IntegerVector>(scan, "rows_start");
  IntegerVector re  = get_field<IntegerVector>(scan, "rows_stop");
  NumericVector times = get_field<NumericVector>(scan, "times");
  List groups = get_field<List>(scan, "events_rowidx_by_time");
  
  const int nR = X.nrow();
  const int p  = X.ncol();
  const int m  = delta.size();
  const int K  = times.size();
  
  NumericVector a_row, d_id;
  scan_sweep_prepare(scan, X, beta, delta, a_row, d_id);
  
  NumericMatrix Utheta_by_edge(m, p);
  NumericVector Udelta_by_edge(m);
  
  NumericVector Wid(m);
  NumericMatrix Wxe(m, p);
  double S0 = 0.0;
  NumericVector WX_colsum(p);
  int si = 0, ei = 0;
  
  const bool is_breslow = (ties_method == "breslow");
  
  for (int k = 0; k < K; ++k) {
    const double t = times[k];
    while (si < nR && ss[si] < t) { add_row(rs[si]-1, id0, a_row, d_id, X, +1.0, S0, Wid, Wxe, WX_colsum); ++si; }
    while (ei < nR && es[ei] < t) { add_row(re[ei]-1, id0, a_row, d_id, X, -1.0, S0, Wid, Wxe, WX_colsum); ++ei; }
    
    IntegerVector ev = groups[k];
    const int mk = ev.size();
    if (!(S0 > 0.0) || !R_finite(S0) || mk == 0) continue;
    
    // 事件端：+x(ev), +1(ev)
    if (mk > 0) {
      for (int ii = 0; ii < mk; ++ii) {
        const int r0 = ev[ii] - 1;
        const int e  = id0[r0] - 1;
        for (int j = 0; j < p; ++j) Utheta_by_edge(e, j) += X(r0, j);
        Udelta_by_edge[e] += 1.0;
      }
    }
    
    if (is_breslow) {
      // 风险集端：- (mk/S0)*[Wxe , Wid]
      const double fac = static_cast<double>(mk) / S0;
      for (int e = 0; e < m; ++e) {
        Udelta_by_edge[e] -= fac * Wid[e];
        for (int j = 0; j < p; ++j) Utheta_by_edge(e, j) -= fac * Wxe(e, j);
      }
    } else {
      // Efron
      NumericVector wevent(m);
      NumericMatrix WXevent(m, p);
      if (mk > 0) {
        for (int ii = 0; ii < mk; ++ii) {
          const int r0 = ev[ii] - 1;
          const int e  = id0[r0] - 1;
          const double wrev = a_row[r0] * d_id[e];
          wevent[e] += wrev;
          for (int j = 0; j < p; ++j) WXevent(e, j) += wrev * X(r0, j);
        }
      }
      double Ek = 0.0; for (int e = 0; e < m; ++e) Ek += wevent[e];
      for (int ell = 0; ell < mk; ++ell) {
        const double alpha = static_cast<double>(ell) / static_cast<double>(mk);
        const double denom = S0 - alpha * Ek;
        if (!(denom > 0.0) || !R_finite(denom)) continue;
        for (int e = 0; e < m; ++e) {
          Udelta_by_edge[e] -= (Wid[e] - alpha * wevent[e]) / denom;
          for (int j = 0; j < p; ++j) {
            Utheta_by_edge(e, j) -= (Wxe(e, j) - alpha * WXevent(e, j)) / denom;
          }
        }
      }
    }
  }
  
  return List::create(
    _["U_theta_by_edge"] = Utheta_by_edge,
    _["U_delta_by_edge"] = Udelta_by_edge
  );
}
