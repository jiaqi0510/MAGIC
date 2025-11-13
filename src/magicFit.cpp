// MAGIC: Methylation Analysis with Genomic Inferred Contexts
// 
// Module: Beta-binomial mixture model optimization
// Description: Maximum likelihood estimation of mixture model parameters using
//              L-BFGS-B optimization with OpenMP parallelization. Implements
//              efficient likelihood and gradient computations for methylation data.
//
// Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
// Contact: mjthompson69@gmail.com
// License: MIT
// Date: 2025

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

const double EPSILON = 1e-15;     
const double MIN_PROB = 1e-12;    
const double MAX_PROB = 1.0 - 1e-12; 
const int MAX_FACTORIAL_CACHE = 500;
const int MAX_N_LOOKUP = 500;

const double MIN_ALPHA_BETA = 0.01;
const double MAX_ALPHA_BETA = 1000.0;
const double MIN_PI = 0.001;
const double MAX_PI = 0.999;
const double MAX_LOG_PARAM = 700.0;

static std::vector<double> log_factorial_table;
static std::vector<std::vector<double>> log_choose_table;
static bool tables_initialized = false;

void init_lookup_tables() {
  if (tables_initialized) return;
  
  log_factorial_table.resize(MAX_FACTORIAL_CACHE + 1);
  log_factorial_table[0] = 0.0;
  for (int i = 1; i <= MAX_FACTORIAL_CACHE; i++) {
    log_factorial_table[i] = log_factorial_table[i-1] + std::log(static_cast<double>(i));
  }
  
  log_choose_table.resize(MAX_N_LOOKUP + 1);
  for (int n = 0; n <= MAX_N_LOOKUP; n++) {
    log_choose_table[n].resize(n + 1);
    for (int m = 0; m <= n; m++) {
      log_choose_table[n][m] = log_factorial_table[n] - 
                               log_factorial_table[m] - 
                               log_factorial_table[n - m];
    }
  }
  
  tables_initialized = true;
}

inline double fast_log_factorial(int n) {
  if (n < 0) return -INFINITY;
  if (n <= MAX_FACTORIAL_CACHE) return log_factorial_table[n];
  return std::lgamma(static_cast<double>(n + 1));
}

inline double fast_log_choose(int n, int m) {
  if (n < 0 || m < 0 || m > n) return -INFINITY;
  if (n <= MAX_N_LOOKUP) return log_choose_table[n][m];
  return fast_log_factorial(n) - fast_log_factorial(m) - fast_log_factorial(n - m);
}

inline double fast_digamma(double x) {
  if (x <= 0) return -INFINITY;
  
  if (x > 8.0) {
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    return std::log(x) - 0.5 * inv_x - inv_x2 * (1.0/12.0 - inv_x2 * (1.0/120.0 - inv_x2/252.0));
  }
  
  double result = 0.0;
  while (x < 8.0) {
    result -= 1.0 / x;
    x += 1.0;
  }
  
  double inv_x = 1.0 / x;
  double inv_x2 = inv_x * inv_x;
  result += std::log(x) - 0.5 * inv_x - inv_x2 * (1.0/12.0 - inv_x2 * (1.0/120.0 - inv_x2/252.0));
  
  return result;
}

inline double fast_lgamma(double x) {
  if (x <= 0) return INFINITY;
  
  static const double log_factorials[] = {
    0.0, 0.0, 0.6931471805599453, 1.791759469228055, 3.1780538303479458,
    4.787491742782046, 6.579251212010101, 8.525161361065414, 10.604602902745251,
    12.801827480081469, 15.104412573075516, 17.502307845873887, 19.987214495661885,
    22.552163853123421, 25.191221182738683, 27.899271383840894, 30.671860106080672,
    33.505073450136891, 36.395445208033053, 39.339884187199495, 42.335616460753485
  };
  
  if (x == std::floor(x) && x >= 1.0 && x <= 20.0) {
    return log_factorials[static_cast<int>(x) - 1];
  }
  
  if (x > 8.0) {
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    return (x - 0.5) * std::log(x) - x + 0.5 * std::log(2.0 * 3.14159265358979323846) + 
           inv_x * (1.0/12.0 - inv_x2 * (1.0/360.0 - inv_x2/1260.0));
  }
  
  double result = 0.0;
  while (x < 8.0) {
    result -= std::log(x);
    x += 1.0;
  }
  
  double inv_x = 1.0 / x;
  double inv_x2 = inv_x * inv_x;
  result += (x - 0.5) * std::log(x) - x + 0.5 * std::log(2.0 * 3.14159265358979323846) + 
            inv_x * (1.0/12.0 - inv_x2 * (1.0/360.0 - inv_x2/1260.0));
  
  return result;
}

inline double fast_lbeta(double a, double b) {
  return fast_lgamma(a) + fast_lgamma(b) - fast_lgamma(a + b);
}

inline void log_to_natural(const NumericVector& log_params, int k_comp,
                          NumericVector& alpha, NumericVector& beta, NumericVector& pi,
                          bool& boundary_hit) {
  boundary_hit = false;
  
  for (int k = 0; k < k_comp; k++) {
    if (log_params[k] > MAX_LOG_PARAM) {
      stop("Parameter overflow: log(alpha[%d]) = %f exceeds safe range (max %.0f)", k, log_params[k], MAX_LOG_PARAM);
    }
    if (log_params[k + k_comp] > MAX_LOG_PARAM) {
      stop("Parameter overflow: log(beta[%d]) = %f exceeds safe range (max %.0f)", k, log_params[k + k_comp], MAX_LOG_PARAM);
    }
    
    double alpha_raw = std::exp(log_params[k]);
    double beta_raw = std::exp(log_params[k + k_comp]);
    
    alpha[k] = std::max(MIN_ALPHA_BETA, std::min(MAX_ALPHA_BETA, alpha_raw));
    beta[k] = std::max(MIN_ALPHA_BETA, std::min(MAX_ALPHA_BETA, beta_raw));
    
    if (alpha_raw < MIN_ALPHA_BETA || alpha_raw > MAX_ALPHA_BETA ||
        beta_raw < MIN_ALPHA_BETA || beta_raw > MAX_ALPHA_BETA) {
      boundary_hit = true;
    }
  }
  
  if (k_comp > 1) {
    int offset = 2 * k_comp;
    double sum_exp = 1.0;
    
    for (int k = 0; k < k_comp - 1; k++) {
      if (log_params[offset + k] > MAX_LOG_PARAM) {
        stop("Parameter overflow: log(pi_transform[%d]) = %f exceeds safe range (max %.0f)", k, log_params[offset + k], MAX_LOG_PARAM);
      }
      sum_exp += std::exp(log_params[offset + k]);
    }
    
    pi[k_comp - 1] = 1.0 / sum_exp;
    for (int k = 0; k < k_comp - 1; k++) {
      pi[k] = std::exp(log_params[offset + k]) / sum_exp;
    }
    
    for (int k = 0; k < k_comp; k++) {
      double pi_raw = pi[k];
      pi[k] = std::max(MIN_PI, std::min(MAX_PI, pi[k]));
      if (pi_raw < MIN_PI || pi_raw > MAX_PI) {
        boundary_hit = true;
      }
    }
    
    double pi_sum = 0.0;
    for (int k = 0; k < k_comp; k++) pi_sum += pi[k];
    for (int k = 0; k < k_comp; k++) pi[k] /= pi_sum;
    
  } else {
    pi[0] = 1.0;
  }
}

inline void check_params(const NumericVector& log_params, int k_comp) {
  int exp_size = (k_comp == 1) ? 2 : (3 * k_comp - 1);
  
  if (log_params.size() != exp_size) {
    stop("Parameter vector size mismatch: expected %d, got %d", exp_size, log_params.size());
  }
}

inline double log_sum_exp(const double* ll, int k_comp) {
  double max_ll = ll[0];
  for (int k = 1; k < k_comp; k++) {
    if (ll[k] > max_ll) max_ll = ll[k];
  }
  
  if (!std::isfinite(max_ll)) return -INFINITY;
  
  double sum_exp = 0.0;
  for (int k = 0; k < k_comp; k++) {
    if (std::isfinite(ll[k])) {
      double exp_diff = std::exp(ll[k] - max_ll);
      if (std::isfinite(exp_diff)) sum_exp += exp_diff;
    }
  }
  
  if (sum_exp <= EPSILON) return -INFINITY;
  return max_ll + std::log(sum_exp);
}

struct CompContext {
  std::vector<double> log_beta_denom;
  std::vector<double> digamma_alpha;
  std::vector<double> digamma_beta;
  std::vector<double> digamma_alpha_beta;
  std::vector<double> log_pi_vals;
  
  CompContext(int k_comp) : 
    log_beta_denom(k_comp), digamma_alpha(k_comp), digamma_beta(k_comp),
    digamma_alpha_beta(k_comp), log_pi_vals(k_comp) {}
    
  void update_cache(const NumericVector& alpha, const NumericVector& beta, 
                   const NumericVector& pi, int k_comp) {
    for (int k = 0; k < k_comp; k++) {
      log_beta_denom[k] = fast_lbeta(alpha[k], beta[k]);
      digamma_alpha[k] = fast_digamma(alpha[k]);
      digamma_beta[k] = fast_digamma(beta[k]);
      digamma_alpha_beta[k] = fast_digamma(alpha[k] + beta[k]);
      log_pi_vals[k] = std::log(std::max(MIN_PROB, pi[k]));
    }
  }
};

// [[Rcpp::export]]
NumericVector compute_ll(NumericMatrix meth_data,
                        NumericMatrix total_data,
                        LogicalMatrix valid_mask,
                        NumericVector log_params) {
  init_lookup_tables();
  
  const int n_features = meth_data.nrow();
  const int n_samples = meth_data.ncol();
  
  if (n_features == 0 || n_samples == 0) {
    stop("Invalid dimensions: n_features=%d, n_samples=%d", n_features, n_samples);
  }
  
  if (total_data.nrow() != n_features || total_data.ncol() != n_samples) {
    stop("Dimension mismatch: meth_data is %dx%d but total_data is %dx%d", 
         n_features, n_samples, total_data.nrow(), total_data.ncol());
  }
  
  if (valid_mask.nrow() != n_features || valid_mask.ncol() != n_samples) {
    stop("Dimension mismatch: meth_data is %dx%d but valid_mask is %dx%d",
         n_features, n_samples, valid_mask.nrow(), valid_mask.ncol());
  }
  
  int k_comp;
  int param_size = log_params.size();
  
  if (param_size == 2) {
    k_comp = 1;
  } else if ((param_size + 1) % 3 == 0 && param_size > 2) {
    k_comp = (param_size + 1) / 3;
  } else {
    stop("Cannot determine k_comp from parameter vector size: %d", param_size);
  }
  
  check_params(log_params, k_comp);
  
  NumericVector alpha(k_comp), beta(k_comp), pi(k_comp);
  bool boundary_hit = false;
  log_to_natural(log_params, k_comp, alpha, beta, pi, boundary_hit);
  
  CompContext ctx(k_comp);
  ctx.update_cache(alpha, beta, pi, k_comp);
  
  NumericVector feature_ll(n_features);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int f = 0; f < n_features; f++) {
    std::vector<double> comp_ll(k_comp);
    
    for (int k = 0; k < k_comp; k++) {
      double site_ll = 0.0;
      bool has_valid_data = false;
      
      for (int s = 0; s < n_samples; s++) {
        const double m = meth_data(f, s);
        const double n = total_data(f, s);
        
        if (valid_mask(f, s) && n > 0 && m >= 0 && m <= n && std::isfinite(m) && std::isfinite(n)) {
          const double log_choose = fast_log_choose(static_cast<int>(n), static_cast<int>(m));
          const double log_beta_num = fast_lbeta(m + alpha[k], n - m + beta[k]);
          site_ll += log_choose + log_beta_num - ctx.log_beta_denom[k];
          has_valid_data = true;
        }
      }
      
      if (has_valid_data) {
        comp_ll[k] = ctx.log_pi_vals[k] + site_ll;
      } else {
        comp_ll[k] = -INFINITY;
      }
      
      if (!std::isfinite(comp_ll[k])) comp_ll[k] = -INFINITY;
    }
    
    feature_ll[f] = log_sum_exp(comp_ll.data(), k_comp);
  }
  
  return feature_ll;
}

// [[Rcpp::export]]
NumericVector compute_gradients(NumericMatrix meth_data,
                               NumericMatrix total_data,
                               LogicalMatrix valid_mask,
                               NumericVector log_params) {
  init_lookup_tables();
  
  const int n_features = meth_data.nrow();
  const int n_samples = meth_data.ncol();
  
  if (n_features == 0 || n_samples == 0) {
    stop("Invalid dimensions: n_features=%d, n_samples=%d", n_features, n_samples);
  }
  
  if (total_data.nrow() != n_features || total_data.ncol() != n_samples) {
    stop("Dimension mismatch: meth_data is %dx%d but total_data is %dx%d",
         n_features, n_samples, total_data.nrow(), total_data.ncol());
  }
  
  if (valid_mask.nrow() != n_features || valid_mask.ncol() != n_samples) {
    stop("Dimension mismatch: meth_data is %dx%d but valid_mask is %dx%d",
         n_features, n_samples, valid_mask.nrow(), valid_mask.ncol());
  }
  
  int k_comp;
  int param_size = log_params.size();
  
  if (param_size == 2) {
    k_comp = 1;
  } else if ((param_size + 1) % 3 == 0 && param_size > 2) {
    k_comp = (param_size + 1) / 3;
  } else {
    stop("Cannot determine k_comp from parameter vector size: %d", param_size);
  }
  
  check_params(log_params, k_comp);
  
  NumericVector alpha(k_comp), beta(k_comp), pi(k_comp);
  bool boundary_hit = false;
  log_to_natural(log_params, k_comp, alpha, beta, pi, boundary_hit);
  
  CompContext ctx(k_comp);
  ctx.update_cache(alpha, beta, pi, k_comp);
  
  NumericVector grad_log_params(log_params.size(), 0.0);
  
  #ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<double> grad_alpha_local(k_comp, 0.0);
    std::vector<double> grad_beta_local(k_comp, 0.0);
    std::vector<double> grad_pi_local(std::max(0, k_comp - 1), 0.0);
    
    #pragma omp for schedule(dynamic)
    for (int f = 0; f < n_features; f++) {
      
      std::vector<double> comp_ll(k_comp);
      for (int k = 0; k < k_comp; k++) {
        double site_ll = 0.0;
        bool has_valid_data = false;
        
        for (int s = 0; s < n_samples; s++) {
          const double m = meth_data(f, s);
          const double n = total_data(f, s);
          
          if (valid_mask(f, s) && n > 0 && m >= 0 && m <= n && std::isfinite(m) && std::isfinite(n)) {
            const double log_choose = fast_log_choose(static_cast<int>(n), static_cast<int>(m));
            const double log_beta_num = fast_lbeta(m + alpha[k], n - m + beta[k]);
            site_ll += log_choose + log_beta_num - ctx.log_beta_denom[k];
            has_valid_data = true;
          }
        }
        
        if (has_valid_data) {
          comp_ll[k] = ctx.log_pi_vals[k] + site_ll;
        } else {
          comp_ll[k] = -INFINITY;
        }
        
        if (!std::isfinite(comp_ll[k])) comp_ll[k] = -INFINITY;
      }
      
      const double log_sum = log_sum_exp(comp_ll.data(), k_comp);
      
      if (!std::isfinite(log_sum)) continue;
      
      std::vector<double> resp(k_comp);
      for (int k = 0; k < k_comp; k++) {
        if (std::isfinite(comp_ll[k])) {
          resp[k] = std::exp(comp_ll[k] - log_sum);
          if (!std::isfinite(resp[k])) resp[k] = 0.0;
        } else {
          resp[k] = 0.0;
        }
      }
      
      for (int k = 0; k < k_comp; k++) {
        const double resp_k = resp[k];
        if (resp_k < EPSILON) continue;
        
        for (int s = 0; s < n_samples; s++) {
          if (valid_mask(f, s)) {
            const double m = meth_data(f, s);
            const double n = total_data(f, s);
            
            if (m >= 0 && n > 0 && m <= n && std::isfinite(m) && std::isfinite(n)) {
              const double digamma_m_alpha = fast_digamma(m + alpha[k]);
              const double digamma_nm_beta = fast_digamma(n - m + beta[k]);
              const double digamma_n_alpha_beta = fast_digamma(n + alpha[k] + beta[k]);
              
              if (std::isfinite(digamma_m_alpha) && std::isfinite(digamma_nm_beta) && 
                  std::isfinite(digamma_n_alpha_beta)) {
                
                double grad_alpha_natural = digamma_m_alpha - ctx.digamma_alpha[k] + 
                                          ctx.digamma_alpha_beta[k] - digamma_n_alpha_beta;
                grad_alpha_local[k] += resp_k * alpha[k] * grad_alpha_natural;
                
                double grad_beta_natural = digamma_nm_beta - ctx.digamma_beta[k] + 
                                         ctx.digamma_alpha_beta[k] - digamma_n_alpha_beta;
                grad_beta_local[k] += resp_k * beta[k] * grad_beta_natural;
              }
            }
          }
        }
      }
      
      if (k_comp > 1) {
        for (int k = 0; k < k_comp - 1; k++) {
          grad_pi_local[k] += resp[k] - pi[k];
        }
      }
    }
    
    #pragma omp critical
    {
      for (int k = 0; k < k_comp; k++) {
        grad_log_params[k] += grad_alpha_local[k];
        grad_log_params[k_comp + k] += grad_beta_local[k];
      }
      
      if (k_comp > 1) {
        int offset = 2 * k_comp;
        for (int k = 0; k < k_comp - 1; k++) {
          grad_log_params[offset + k] += grad_pi_local[k];
        }
      }
    }
  }
  #else
  stop("This code requires OpenMP for performance on large datasets. Please compile with -fopenmp");
  #endif
  
  return grad_log_params;
}

// [[Rcpp::export]]
List check_parameter_boundaries(NumericVector log_params) {
  int param_size = log_params.size();
  int k_comp;
  
  if (param_size == 2) {
    k_comp = 1;
  } else if ((param_size + 1) % 3 == 0 && param_size > 2) {
    k_comp = (param_size + 1) / 3;
  } else {
    stop("Cannot determine k_comp from parameter vector size: %d", param_size);
  }
  
  NumericVector alpha(k_comp), beta(k_comp), pi(k_comp);
  bool boundary_hit = false;
  log_to_natural(log_params, k_comp, alpha, beta, pi, boundary_hit);
  
  const double BOUNDARY_TOL = 0.02;
  
  LogicalVector alpha_at_min(k_comp);
  LogicalVector alpha_at_max(k_comp);
  LogicalVector beta_at_min(k_comp);
  LogicalVector beta_at_max(k_comp);
  LogicalVector pi_at_min(k_comp);
  LogicalVector pi_at_max(k_comp);
  
  for (int k = 0; k < k_comp; k++) {
    alpha_at_min[k] = (alpha[k] < MIN_ALPHA_BETA + BOUNDARY_TOL);
    alpha_at_max[k] = (alpha[k] > MAX_ALPHA_BETA - BOUNDARY_TOL);
    beta_at_min[k] = (beta[k] < MIN_ALPHA_BETA + BOUNDARY_TOL);
    beta_at_max[k] = (beta[k] > MAX_ALPHA_BETA - BOUNDARY_TOL);
    pi_at_min[k] = (pi[k] < MIN_PI + BOUNDARY_TOL);
    pi_at_max[k] = (pi[k] > MAX_PI - BOUNDARY_TOL);
  }
  
  return List::create(
    Named("alpha_at_min") = alpha_at_min,
    Named("alpha_at_max") = alpha_at_max,
    Named("beta_at_min") = beta_at_min,
    Named("beta_at_max") = beta_at_max,
    Named("pi_at_min") = pi_at_min,
    Named("pi_at_max") = pi_at_max,
    Named("alpha") = alpha,
    Named("beta") = beta,
    Named("pi") = pi
  );
}
