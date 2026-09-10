// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace Rcpp;

// --- sanitize possible bad macros from other headers (Windows etc.) ---
#ifdef parameter
#undef parameter
#endif
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#ifdef ERROR
#undef ERROR
#endif

// -------- small helpers --------
static inline void symmetrize_inplace(NumericMatrix& A) {
  const int n = A.nrow();
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double v = 0.5 * (A(i, j) + A(j, i));
      A(i, j) = v;
      A(j, i) = v;
    }
  }
}

static bool chol_lower(const NumericMatrix& A, NumericMatrix& L) {
  const int n = A.nrow();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) L(i, j) = 0.0;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = A(i, j);
      for (int k = 0; k < j; ++k) sum -= L(i, k) * L(j, k);
      if (i == j) {
        if (!(sum > 0.0) || !R_finite(sum)) return false;
        L(i, j) = std::sqrt(sum);
      } else {
        if (L(j, j) == 0.0) return false;
        L(i, j) = sum / L(j, j);
      }
    }
  }
  return true;
}

static NumericMatrix invert_via_cholL(const NumericMatrix& L) {
  const int n = L.nrow();
  // Linv: solve L * Linv = I  (forward substitution for each column)
  NumericMatrix Linv(n, n);
  for (int col = 0; col < n; ++col) {
    for (int r = 0; r < n; ++r) {
      double rhs = (r == col) ? 1.0 : 0.0;
      for (int k = 0; k < r; ++k) rhs -= L(r, k) * Linv(k, col);
      Linv(r, col) = rhs / L(r, r);
    }
  }
  // S^{-1} = Linv^T * Linv
  NumericMatrix Sinv(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double acc = 0.0;
      for (int k = 0; k < n; ++k) acc += Linv(k, i) * Linv(k, j);
      Sinv(i, j) = acc;
      Sinv(j, i) = acc;
    }
  }
  return Sinv;
}

// -------- R callbacks we will call from C++ --------
// NOTE: we pass them in from R as Rcpp::Function, so we don't include their headers.
static inline NumericMatrix call_Hdt(const Function& Hdt, const NumericMatrix& V) {
  return as<NumericMatrix>(Hdt(V));
}
static inline NumericMatrix call_Hdd_solve(const Function& HddSolve, const NumericMatrix& RHS) {
  return as<NumericMatrix>(HddSolve(RHS));
}
static inline NumericMatrix call_Htd(const Function& Htd, const NumericMatrix& U) {
  return as<NumericMatrix>(Htd(U));
}

// --------------------------------------------------
// S = I_tt - Htd * Hdd^{-1} * Hdt  (blocked assembly) + chol inverse
// --------------------------------------------------
// [[Rcpp::export]]
Rcpp::List schur_theta_block_cpp_min(
    Rcpp::NumericMatrix I_tt,     // p x p
    Rcpp::Function Hdt_matmul,    // V(p×q) -> m×q
    Rcpp::Function Hdd_solve_mat, // RHS(m×q) -> sol(m×q)
    Rcpp::Function Htd_matmul,    // U(m×q) -> p×q
    double ridge = 1e-8,
    int ridge_tries = 5,
    int block_size = 32,
    bool verbose = true
) {
  const int p = I_tt.nrow();
  if (I_tt.ncol() != p) stop("I_tt must be square (p x p).");
  if (block_size < 1) block_size = 1;
  
  // Start from I_tt (already may include small ridge outside)
  NumericMatrix S = clone(I_tt);
  
  // Blocked identity to save memory
  const int nb = (p + block_size - 1) / block_size;
  for (int b = 0; b < nb; ++b) {
    const int c0 = b * block_size;
    const int k  = std::min(block_size, p - c0);
    
    // V = I(:, cols)
    NumericMatrix V(p, k);
    for (int j = 0; j < k; ++j) {
      V(c0 + j, j) = 1.0;
    }
    
    // Q = Hdt * V   (m x k)
    NumericMatrix Q = call_Hdt(Hdt_matmul, V);
    // Z = Hdd^{-1} Q
    NumericMatrix Z = call_Hdd_solve(Hdd_solve_mat, Q);
    // T = Htd * Z   (p x k)
    NumericMatrix T = call_Htd(Htd_matmul, Z);
    
    // S(:, cols) -= T
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < p; ++i) {
        S(i, c0 + j) -= T(i, j);
      }
    }
    
    if (verbose) {
      double pct = 100.0 * static_cast<double>(b + 1) / static_cast<double>(nb);
      Rcpp::Rcout << "\r[Schur] " << std::fixed << std::setprecision(1) << pct << "%";
      Rcpp::Rcout.flush();
    }
    if ((b & 3) == 0) Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcpp::Rcout << "\n";
  
  // Symmetrize & chol with ridge fallback
  symmetrize_inplace(S);
  NumericMatrix L(p, p);
  bool ok = chol_lower(S, L);
  double cur_ridge = ridge;
  int tried = 0;
  while (!ok && tried < ridge_tries) {
    if (verbose) Rcpp::Rcout << "[schur] chol failed; add ridge=" << cur_ridge << "\n";
    for (int i = 0; i < p; ++i) S(i, i) += cur_ridge;
    ok = chol_lower(S, L);
    cur_ridge *= 10.0;
    ++tried;
  }
  if (!ok) stop("Cholesky failed even after ridge attempts.");
  
  NumericMatrix S_inv = invert_via_cholL(L);
  return List::create(_["S"] = S, _["S_inv"] = S_inv);
}
