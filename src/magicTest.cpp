// MAGIC: Methylation Analysis with Genomic Inferred Contexts
//
// Module: Bayesian differential methylation testing
// Description: Implements Bayesian hypothesis testing for differential methylation
//              between conditions using mixture-weighted dispersion estimation.
//              Includes correlation analysis and component-wise testing.
//
// Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
// Contact: mjthompson69@gmail.com
// License: MIT
// Date: 2025

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

const double LOG_ZERO = -std::numeric_limits<double>::infinity();
const double TOLERANCE = 1e-10;

struct TestResult {
    double mean1, mean2, difference, p_value, wald_stat;
    double bf_10, prob_diff;
    int dom_comp1, dom_comp2;
    double entropy1, entropy2;
    double pearson_r, pearson_p, spearman_rho, spearman_p, magic_r, magic_p;
    bool valid;
    
    TestResult() : mean1(0), mean2(0), difference(0), p_value(1), wald_stat(0),
                   bf_10(1), prob_diff(0.5), dom_comp1(1), dom_comp2(1),
                   entropy1(0), entropy2(0), pearson_r(0), pearson_p(1),
                   spearman_rho(0), spearman_p(1), magic_r(0), magic_p(1),
                   valid(false) {}
};

// Beta-binomial density (DSS exact implementation)
inline double dbb_exact(double size, double x, double mu, double phi, bool log_result = true) {
    double tmp = 1.0/phi - 1.0;
    double alpha = mu * tmp;
    double beta = tmp - alpha;
    
    double lchoose = std::lgamma(size + 1.0) - std::lgamma(x + 1.0) - std::lgamma(size - x + 1.0);
    double lbeta1 = std::lgamma(beta) + std::lgamma(alpha) - std::lgamma(beta + alpha);
    double lbeta2 = std::lgamma(size - x + beta) + std::lgamma(x + alpha) - std::lgamma(size + beta + alpha);
    
    double result = lchoose - lbeta1 + lbeta2;
    return log_result ? result : std::exp(result);
}

// Penalized likelihood for dispersion estimation
inline double plik_logN(double phi_log, const std::vector<double>& total, 
                       const std::vector<double>& meth, const std::vector<double>& mu, 
                       double m0, double tau) {
    double phi = std::exp(phi_log);
    double dbb_sum = 0.0;
    
    for (size_t i = 0; i < total.size(); i++) {
        dbb_sum += dbb_exact(total[i], meth[i], mu[i], phi, true);
    }
    
    double dnorm = -0.5 * std::log(2.0 * M_PI * tau * tau) - 0.5 * std::pow((phi_log - m0) / tau, 2.0);
    return -(dbb_sum + dnorm);
}

// Golden section optimization for dispersion
inline double optimize_dispersion(const std::vector<double>& total, const std::vector<double>& meth,
                                 const std::vector<double>& mu, double m0, double tau) {
    double a = -5.0, b = std::log(0.99), tol = 1e-3;
    const double golden = (std::sqrt(5.0) - 1.0) / 2.0;
    
    double c = b - golden * (b - a);
    double d = a + golden * (b - a);
    
    while (std::abs(b - a) > tol) {
        double fc = plik_logN(c, total, meth, mu, m0, tau);
        double fd = plik_logN(d, total, meth, mu, m0, tau);
        
        if (fc < fd) { b = d; d = c; c = b - golden * (b - a); }
        else { a = c; c = d; d = a + golden * (b - a); }
    }
    
    return std::exp((a + b) / 2.0);
}

// Global dispersion prior estimation - MATRIX VERSION
// [[Rcpp::export]]
NumericVector estimate_global_prior_from_matrix(const NumericMatrix& meth_matrix,
                                                const NumericMatrix& total_matrix) {
    int n_features = meth_matrix.nrow();
    int n_samples = meth_matrix.ncol();
    
    if (n_features == 0 || n_samples <= 1) {
        warning("Insufficient data for prior estimation (n_features=%d, n_samples=%d), using defaults", n_features, n_samples);
        return NumericVector::create(-3.0, 1.0);
    }
    
    if (total_matrix.nrow() != n_features || total_matrix.ncol() != n_samples) {
        stop("Dimension mismatch: meth_matrix is %dx%d but total_matrix is %dx%d",
             n_features, n_samples, total_matrix.nrow(), total_matrix.ncol());
    }
    
    std::vector<double> log_phi_vals;
    log_phi_vals.reserve(n_features / 10);
    
    for (int site = 0; site < n_features; site++) {
        bool valid = true;
        int valid_samples = 0;
        for (int j = 0; j < n_samples; j++) {
            if (total_matrix(site, j) >= 10) {
                valid_samples++;
            }
        }
        
        if (valid_samples < std::max(3, n_samples / 4)) continue;
        
        double site_mean = 0.0, site_var = 0.0;
        int count = 0;
        
        for (int j = 0; j < n_samples; j++) {
            if (total_matrix(site, j) >= 10) {
                double prop = meth_matrix(site, j) / total_matrix(site, j);
                site_mean += prop;
                count++;
            }
        }
        
        if (count < 3) continue;
        site_mean /= count;
        
        if (site_mean <= 0.0) site_mean = 1e-5;
        if (site_mean >= 1.0) site_mean = 1.0 - 1e-5;
        
        for (int j = 0; j < n_samples; j++) {
            if (total_matrix(site, j) >= 10) {
                double prop = meth_matrix(site, j) / total_matrix(site, j);
                site_var += (prop - site_mean) * (prop - site_mean);
            }
        }
        site_var /= (count - 1);
        
        if (site_var > 0.0) {
            double phi = site_var / (site_mean * (1.0 - site_mean));
            if (phi > 0.0 && std::isfinite(phi)) {
                log_phi_vals.push_back(std::log(phi));
            }
        }
    }
    
    if (log_phi_vals.size() < 50) return NumericVector::create(-3.0, 1.0);
    
    std::sort(log_phi_vals.begin(), log_phi_vals.end());
    size_t n = log_phi_vals.size();
    double median = (n % 2 == 0) ? (log_phi_vals[n/2-1] + log_phi_vals[n/2]) / 2.0 : log_phi_vals[n/2];
    double q1 = log_phi_vals[n/4], q3 = log_phi_vals[3*n/4];
    
    return NumericVector::create(median, (q3 - q1) / 1.39);
}

// Posterior computation for mixture analysis
std::vector<double> compute_posteriors(const std::vector<double>& meth, const std::vector<double>& total,
                                      const std::vector<double>& pi, const std::vector<double>& alpha,
                                      const std::vector<double>& beta) {
    const int k_comp = pi.size();
    std::vector<double> log_probs(k_comp, LOG_ZERO);
    
    for (int k = 0; k < k_comp; k++) {
        double log_lik = 0.0;
        bool valid = true;
        
        for (size_t i = 0; i < meth.size(); i++) {
            if (total[i] > 0) {
                double y = meth[i], n = total[i], a = alpha[k], b = beta[k];
                if (n <= 0 || y < 0 || y > n || a <= 0 || b <= 0) { valid = false; break; }
                
                double meth_lik = std::lgamma(n + 1.0) - std::lgamma(y + 1.0) - std::lgamma(n - y + 1.0) +
                                 std::lgamma(y + a) + std::lgamma(n - y + b) - std::lgamma(n + a + b) +
                                 std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
                
                if (std::isfinite(meth_lik)) log_lik += meth_lik;
                else { valid = false; break; }
            }
        }
        
        if (valid) log_probs[k] = std::log(pi[k]) + log_lik;
    }
    
    double max_log = *std::max_element(log_probs.begin(), log_probs.end());
    if (!std::isfinite(max_log)) {
        return std::vector<double>(k_comp, 1.0 / k_comp);
    }
    
    double sum_exp = 0.0;
    for (double& lp : log_probs) {
        if (lp != LOG_ZERO) {
            lp = std::exp(lp - max_log);
            sum_exp += lp;
        } else {
            lp = 0.0;
        }
    }
    
    if (sum_exp > 0.0) {
        for (double& p : log_probs) p /= sum_exp;
    } else {
        std::fill(log_probs.begin(), log_probs.end(), 1.0 / k_comp);
    }
    
    return log_probs;
}

// MIXTURE-WEIGHTED DISPERSION TEST (formerly w3/multicomp)
TestResult compute_mixture_test(const std::vector<double>& meth1, const std::vector<double>& total1,
                                const std::vector<double>& meth2, const std::vector<double>& total2,
                                const std::vector<double>& pi, const std::vector<double>& alpha,
                                const std::vector<double>& beta, double global_m0, double global_tau) {
    TestResult result;
    
    // Calculate means
    double mean1 = 0.0, mean2 = 0.0, sum_total1 = 0.0, sum_total2 = 0.0;
    for (size_t i = 0; i < meth1.size(); i++) {
        mean1 += meth1[i];
        sum_total1 += total1[i];
    }
    for (size_t i = 0; i < meth2.size(); i++) {
        mean2 += meth2[i];
        sum_total2 += total2[i];
    }
    mean1 /= sum_total1;
    mean2 /= sum_total2;
    
    std::vector<double> mu1(meth1.size(), mean1), mu2(meth2.size(), mean2);
    
    // Compute posteriors for mixture-weighted dispersion
    auto posts1 = compute_posteriors(meth1, total1, pi, alpha, beta);
    auto posts2 = compute_posteriors(meth2, total2, pi, alpha, beta);
    
    double mix_phi1 = 0.0, mix_phi2 = 0.0;
    for (size_t k = 0; k < pi.size(); k++) {
        double comp_phi = 1.0 / (alpha[k] + beta[k] + 1.0);
        mix_phi1 += posts1[k] * comp_phi;
        mix_phi2 += posts2[k] * comp_phi;
    }
    
    double phi1 = optimize_dispersion(total1, meth1, mu1, std::log(mix_phi1), 0.5);
    double phi2 = optimize_dispersion(total2, meth2, mu2, std::log(mix_phi2), 0.5);
    
    // Weight calculation
    std::vector<double> wt1(meth1.size()), wt2(meth2.size());
    for (size_t i = 0; i < total1.size(); i++) wt1[i] = 1.0 / (1.0 + (total1[i] - 1.0) * phi1);
    for (size_t i = 0; i < total2.size(); i++) wt2[i] = 1.0 / (1.0 + (total2[i] - 1.0) * phi2);
    
    double mean_wt1 = 0.0, mean_wt2 = 0.0;
    for (double w : wt1) mean_wt1 += w; mean_wt1 /= wt1.size();
    for (double w : wt2) mean_wt2 += w; mean_wt2 /= wt2.size();
    for (size_t i = 0; i < wt1.size(); i++) wt1[i] /= mean_wt1;
    for (size_t i = 0; i < wt2.size(); i++) wt2[i] /= mean_wt2;
    
    // Apply weights and re-estimate means
    double sum_x1_wt = 0.0, sum_n1_wt = 0.0, sum_x2_wt = 0.0, sum_n2_wt = 0.0;
    for (size_t i = 0; i < meth1.size(); i++) {
        sum_x1_wt += meth1[i] * wt1[i];
        sum_n1_wt += total1[i] * wt1[i];
    }
    for (size_t i = 0; i < meth2.size(); i++) {
        sum_x2_wt += meth2[i] * wt2[i];
        sum_n2_wt += total2[i] * wt2[i];
    }
    
    result.mean1 = sum_x1_wt / sum_n1_wt;
    result.mean2 = sum_x2_wt / sum_n2_wt;
    result.difference = result.mean2 - result.mean1;
    
    // Variance calculation
    double total_cov1 = 0.0, total_cov2 = 0.0;
    for (double t : total1) total_cov1 += t;
    for (double t : total2) total_cov2 += t;
    
    double var1_sum = 0.0, var2_sum = 0.0;
    for (size_t i = 0; i < total1.size(); i++) {
        var1_sum += total1[i] * result.mean1 * (1.0 - result.mean1) * (1.0 + (total1[i] - 1.0) * phi1);
    }
    for (size_t i = 0; i < total2.size(); i++) {
        var2_sum += total2[i] * result.mean2 * (1.0 - result.mean2) * (1.0 + (total2[i] - 1.0) * phi2);
    }
    
    double var1 = var1_sum / (total_cov1 * total_cov1);
    double var2 = var2_sum / (total_cov2 * total_cov2);
    double total_var = std::max(var1 + var2, 1e-5);
    
    result.wald_stat = result.difference / std::sqrt(total_var);
    result.p_value = 2.0 * 0.5 * std::erfc(std::abs(result.wald_stat) / std::sqrt(2.0));
    result.valid = true;
    
    return result;
}

// MAGIC Bayes Factor test
TestResult compute_magic_bayes(const std::vector<double>& meth_counts, const std::vector<double>& total_counts,
                              const std::vector<int>& group_labels, const std::vector<double>& pi,
                              const std::vector<double>& alpha, const std::vector<double>& beta,
                              const std::vector<double>& alpha_0) {
    TestResult result;
    const int k_comp = pi.size();
    
    std::vector<double> n_tilde_1(k_comp, 0.0), n_tilde_2(k_comp, 0.0), n_tilde_pooled(k_comp, 0.0);
    
    // Separate data by groups for posterior calculations
    std::vector<double> meth1, total1, meth2, total2;
    for (size_t i = 0; i < meth_counts.size(); i++) {
        if (total_counts[i] > 0) {
            auto posts = compute_posteriors({meth_counts[i]}, {total_counts[i]}, pi, alpha, beta);
            
            for (int k = 0; k < k_comp; k++) {
                if (group_labels[i] == 0) n_tilde_1[k] += posts[k];
                else if (group_labels[i] == 1) n_tilde_2[k] += posts[k];
                n_tilde_pooled[k] += posts[k];
            }
            
            if (group_labels[i] == 0) {
                meth1.push_back(meth_counts[i]);
                total1.push_back(total_counts[i]);
            } else if (group_labels[i] == 1) {
                meth2.push_back(meth_counts[i]);
                total2.push_back(total_counts[i]);
            }
        }
    }
    
    // Calculate marginal likelihoods
    auto log_mvbeta = [](const std::vector<double>& alpha) {
        double sum_alpha = 0.0, sum_lgamma = 0.0;
        for (double a : alpha) { sum_alpha += a; sum_lgamma += std::lgamma(a); }
        return sum_lgamma - std::lgamma(sum_alpha);
    };
    
    std::vector<double> n_pooled_plus_alpha0(k_comp), n_group1_plus_alpha0(k_comp), n_group2_plus_alpha0(k_comp);
    for (int k = 0; k < k_comp; k++) {
        n_pooled_plus_alpha0[k] = n_tilde_pooled[k] + alpha_0[k];
        n_group1_plus_alpha0[k] = n_tilde_1[k] + alpha_0[k];
        n_group2_plus_alpha0[k] = n_tilde_2[k] + alpha_0[k];
    }
    
    double log_marg_lik_m0 = log_mvbeta(n_pooled_plus_alpha0) - log_mvbeta(alpha_0);
    double log_marg_lik_m1 = log_mvbeta(n_group1_plus_alpha0) + log_mvbeta(n_group2_plus_alpha0) - 2.0 * log_mvbeta(alpha_0);
    
    double log_bf_10 = log_marg_lik_m1 - log_marg_lik_m0;
    
    if (std::isfinite(log_bf_10)) {
        result.bf_10 = std::exp(log_bf_10);
        result.prob_diff = result.bf_10 / (1.0 + result.bf_10);
    } else {
        result.bf_10 = (log_bf_10 > 0) ? std::numeric_limits<double>::infinity() : 0.0;
        result.prob_diff = (log_bf_10 > 0) ? 1.0 : 0.0;
    }
    
    // Calculate group posteriors and dominant components/entropy
    if (!meth1.empty() && !meth2.empty()) {
        auto posts1 = compute_posteriors(meth1, total1, pi, alpha, beta);
        auto posts2 = compute_posteriors(meth2, total2, pi, alpha, beta);
        
        // Find dominant components
        result.dom_comp1 = std::distance(posts1.begin(), std::max_element(posts1.begin(), posts1.end())) + 1;
        result.dom_comp2 = std::distance(posts2.begin(), std::max_element(posts2.begin(), posts2.end())) + 1;
        
        // Calculate entropies
        result.entropy1 = 0.0;
        result.entropy2 = 0.0;
        for (int k = 0; k < k_comp; k++) {
            if (posts1[k] > 0.0) result.entropy1 -= posts1[k] * std::log(posts1[k]);
            if (posts2[k] > 0.0) result.entropy2 -= posts2[k] * std::log(posts2[k]);
        }
        
        // Calculate group means
        result.mean1 = 0.0;
        result.mean2 = 0.0;
        for (int k = 0; k < k_comp; k++) {
            result.mean1 += posts1[k] * (alpha[k] / (alpha[k] + beta[k]));
            result.mean2 += posts2[k] * (alpha[k] / (alpha[k] + beta[k]));
        }
    }
    
    result.valid = true;
    return result;
}

// CORRELATION FUNCTIONS
TestResult compute_correlations(const std::vector<double>& meth_counts, const std::vector<double>& total_counts,
                               const std::vector<double>& trait_values, const std::vector<double>& pi = {},
                               const std::vector<double>& alpha = {}, const std::vector<double>& beta = {}) {
    TestResult result;
    const size_t n = meth_counts.size();
    
    std::vector<double> proportions(n);
    for (size_t i = 0; i < n; i++) proportions[i] = meth_counts[i] / total_counts[i];
    
    // Rank transform for Spearman
    auto rank_transform = [](const std::vector<double>& vals) {
        std::vector<std::pair<double, size_t>> indexed(vals.size());
        for (size_t i = 0; i < vals.size(); i++) indexed[i] = {vals[i], i};
        std::sort(indexed.begin(), indexed.end());
        
        std::vector<double> ranks(vals.size());
        for (size_t i = 0; i < indexed.size(); i++) ranks[indexed[i].second] = i + 1.0;
        return ranks;
    };
    
    // Pearson correlation
    double mean_prop = 0.0, mean_trait = 0.0;
    for (size_t i = 0; i < n; i++) {
        mean_prop += proportions[i];
        mean_trait += trait_values[i];
    }
    mean_prop /= n;
    mean_trait /= n;
    
    double cov = 0.0, var_prop = 0.0, var_trait = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dp = proportions[i] - mean_prop;
        double dt = trait_values[i] - mean_trait;
        cov += dp * dt;
        var_prop += dp * dp;
        var_trait += dt * dt;
    }
    
    if (var_prop > 0 && var_trait > 0) {
        result.pearson_r = cov / std::sqrt(var_prop * var_trait);
        double t_stat = result.pearson_r * std::sqrt((n - 2) / (1.0 - result.pearson_r * result.pearson_r));
        result.pearson_p = 2.0 * 0.5 * std::erfc(std::abs(t_stat) / std::sqrt(2.0)) * std::sqrt(n - 2.0);
    }
    
    // Spearman correlation
    auto ranks_prop = rank_transform(proportions);
    auto ranks_trait = rank_transform(trait_values);
    
    double mean_rank_prop = 0.0, mean_rank_trait = 0.0;
    for (size_t i = 0; i < n; i++) {
        mean_rank_prop += ranks_prop[i];
        mean_rank_trait += ranks_trait[i];
    }
    mean_rank_prop /= n;
    mean_rank_trait /= n;
    
    double cov_ranks = 0.0, var_rank_prop = 0.0, var_rank_trait = 0.0;
    for (size_t i = 0; i < n; i++) {
        double drp = ranks_prop[i] - mean_rank_prop;
        double drt = ranks_trait[i] - mean_rank_trait;
        cov_ranks += drp * drt;
        var_rank_prop += drp * drp;
        var_rank_trait += drt * drt;
    }
    
    if (var_rank_prop > 0 && var_rank_trait > 0) {
        result.spearman_rho = cov_ranks / std::sqrt(var_rank_prop * var_rank_trait);
        double t_stat = result.spearman_rho * std::sqrt((n - 2) / (1.0 - result.spearman_rho * result.spearman_rho));
        result.spearman_p = 2.0 * 0.5 * std::erfc(std::abs(t_stat) / std::sqrt(2.0)) * std::sqrt(n - 2.0);
    }
    
    // MAGIC correlation (mixture-based)
    if (!pi.empty() && !alpha.empty() && !beta.empty()) {
        std::vector<double> magic_props(n);
        for (size_t i = 0; i < n; i++) {
            auto posts = compute_posteriors({meth_counts[i]}, {total_counts[i]}, pi, alpha, beta);
            magic_props[i] = 0.0;
            for (size_t k = 0; k < pi.size(); k++) {
                magic_props[i] += posts[k] * (alpha[k] / (alpha[k] + beta[k]));
            }
        }
        
        double mean_magic = 0.0;
        for (double mp : magic_props) mean_magic += mp;
        mean_magic /= n;
        
        double cov_magic = 0.0, var_magic = 0.0;
        for (size_t i = 0; i < n; i++) {
            double dm = magic_props[i] - mean_magic;
            double dt = trait_values[i] - mean_trait;
            cov_magic += dm * dt;
            var_magic += dm * dm;
        }
        
        if (var_magic > 0 && var_trait > 0) {
            result.magic_r = cov_magic / std::sqrt(var_magic * var_trait);
            double t_stat = result.magic_r * std::sqrt((n - 2) / (1.0 - result.magic_r * result.magic_r));
            result.magic_p = 2.0 * 0.5 * std::erfc(std::abs(t_stat) / std::sqrt(2.0)) * std::sqrt(n - 2.0);
        }
    }
    
    result.valid = true;
    return result;
}

// [[Rcpp::export]]
DataFrame magic_tests_chunked(const NumericMatrix& meth_matrix,
                             const NumericMatrix& total_matrix,
                             const LogicalMatrix& val_matrix,
                             const IntegerMatrix& group_matrix,
                             const NumericMatrix& trait_matrix,
                             const NumericVector& pi_universal,
                             const NumericVector& alpha,
                             const NumericVector& beta,
                             const NumericVector& alpha_0,
                             double global_m0,
                             double global_tau,
                             std::string test_type,
                             int n_threads,
                             bool is_correlation) {
    
    int n_features = meth_matrix.nrow();
    int n_samples = meth_matrix.ncol();
    int k_comp = pi_universal.size();
    
    if (n_features == 0 || n_samples == 0) {
        stop("Invalid dimensions: n_features=%d, n_samples=%d", n_features, n_samples);
    }
    
    if (total_matrix.nrow() != n_features || total_matrix.ncol() != n_samples) {
        stop("Dimension mismatch: meth_matrix is %dx%d but total_matrix is %dx%d",
             n_features, n_samples, total_matrix.nrow(), total_matrix.ncol());
    }
    
    if (val_matrix.nrow() != n_features || val_matrix.ncol() != n_samples) {
        stop("Dimension mismatch: meth_matrix is %dx%d but val_matrix is %dx%d",
             n_features, n_samples, val_matrix.nrow(), val_matrix.ncol());
    }
    
    if (alpha.size() != k_comp || beta.size() != k_comp) {
        stop("Parameter size mismatch: k_comp=%d but alpha.size()=%d, beta.size()=%d",
             k_comp, alpha.size(), beta.size());
    }
    
    for (int k = 0; k < k_comp; k++) {
        if (alpha[k] <= 0 || !std::isfinite(alpha[k])) {
            stop("Invalid alpha[%d] = %f (must be positive and finite)", k, alpha[k]);
        }
        if (beta[k] <= 0 || !std::isfinite(beta[k])) {
            stop("Invalid beta[%d] = %f (must be positive and finite)", k, beta[k]);
        }
        if (pi_universal[k] < 0 || pi_universal[k] > 1 || !std::isfinite(pi_universal[k])) {
            stop("Invalid pi[%d] = %f (must be in [0,1])", k, pi_universal[k]);
        }
    }
    
    std::vector<double> pi_vec(pi_universal.begin(), pi_universal.end());
    std::vector<double> alpha_vec(alpha.begin(), alpha.end());
    std::vector<double> beta_vec(beta.begin(), beta.end());
    std::vector<double> alpha_0_vec(alpha_0.begin(), alpha_0.end());
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #endif
    
    bool run_mixture = (test_type == "mixture" || test_type == "all");
    bool run_magic = (test_type == "magic" || test_type == "all");
    
    std::unique_ptr<std::vector<TestResult>> mixture_results;
    std::unique_ptr<std::vector<TestResult>> magic_results;
    std::unique_ptr<std::vector<TestResult>> corr_results;
    
    if (is_correlation) {
        corr_results = std::make_unique<std::vector<TestResult>>(n_features);
    } else {
        if (run_mixture) mixture_results = std::make_unique<std::vector<TestResult>>(n_features);
        if (run_magic) magic_results = std::make_unique<std::vector<TestResult>>(n_features);
    }
    
    #pragma omp parallel for schedule(dynamic) if(n_features > 100)
    for (int f = 0; f < n_features; f++) {
        std::vector<double> meth_vec, total_vec, trait_vec;
        std::vector<int> group_vec;
        
        for (int s = 0; s < n_samples; s++) {
            if (val_matrix(f, s)) {
                meth_vec.push_back(meth_matrix(f, s));
                total_vec.push_back(total_matrix(f, s));
                if (is_correlation) {
                    trait_vec.push_back(trait_matrix(f, s));
                } else {
                    group_vec.push_back(group_matrix(f, s));
                }
            }
        }
        
        if (is_correlation && meth_vec.size() >= 3) {
            if (corr_results) {
                (*corr_results)[f] = compute_correlations(meth_vec, total_vec, trait_vec, pi_vec, alpha_vec, beta_vec);
            }
        } else if (!is_correlation) {
            std::vector<double> meth1, total1, meth2, total2;
            for (size_t i = 0; i < group_vec.size(); i++) {
                if (group_vec[i] == 0) {
                    meth1.push_back(meth_vec[i]);
                    total1.push_back(total_vec[i]);
                } else if (group_vec[i] == 1) {
                    meth2.push_back(meth_vec[i]);
                    total2.push_back(total_vec[i]);
                }
            }
            
            if (!meth1.empty() && !meth2.empty()) {
                if (run_mixture) {
                    (*mixture_results)[f] = compute_mixture_test(meth1, total1, meth2, total2, 
                                                                 pi_vec, alpha_vec, beta_vec, global_m0, global_tau);
                }
                if (run_magic) {
                    (*magic_results)[f] = compute_magic_bayes(meth_vec, total_vec, group_vec, 
                                                              pi_vec, alpha_vec, beta_vec, alpha_0_vec);
                }
            }
        }
    }

    // Build output DataFrame
    if (is_correlation) {
        std::vector<double> pears(n_features), pear_ps(n_features), spears(n_features), spear_ps(n_features);
        std::vector<double> magics(n_features), magic_ps(n_features);
        
        for (int f = 0; f < n_features; f++) {
            if (corr_results && (*corr_results)[f].valid) {
                pears[f] = (*corr_results)[f].pearson_r;
                pear_ps[f] = (*corr_results)[f].pearson_p;
                spears[f] = (*corr_results)[f].spearman_rho;
                spear_ps[f] = (*corr_results)[f].spearman_p;
                magics[f] = (*corr_results)[f].magic_r;
                magic_ps[f] = (*corr_results)[f].magic_p;
            } else {
                pears[f] = pear_ps[f] = spears[f] = spear_ps[f] = magics[f] = magic_ps[f] = NA_REAL;
            }
        }
        
        return DataFrame::create(
            Named("pearR") = NumericVector(pears.begin(), pears.end()),
            Named("pearP") = NumericVector(pear_ps.begin(), pear_ps.end()),
            Named("spearRho") = NumericVector(spears.begin(), spears.end()),
            Named("spearP") = NumericVector(spear_ps.begin(), spear_ps.end()),
            Named("magicR") = NumericVector(magics.begin(), magics.end()),
            Named("magicP") = NumericVector(magic_ps.begin(), magic_ps.end())
        );
    }
    
    // Multi-method group comparison output
    List df_columns;
    
    // Magic results
    if (run_magic && magic_results) {
        std::vector<double> magic_bfs(n_features), magic_probs(n_features), magic_means1(n_features), magic_means2(n_features);
        std::vector<int> magic_comps1(n_features), magic_comps2(n_features);
        std::vector<double> magic_ents1(n_features), magic_ents2(n_features);
        
        for (int f = 0; f < n_features; f++) {
            if ((*magic_results)[f].valid) {
                magic_bfs[f] = (*magic_results)[f].bf_10;
                magic_probs[f] = (*magic_results)[f].prob_diff;
                magic_means1[f] = (*magic_results)[f].mean1;
                magic_means2[f] = (*magic_results)[f].mean2;
                magic_comps1[f] = (*magic_results)[f].dom_comp1;
                magic_comps2[f] = (*magic_results)[f].dom_comp2;
                magic_ents1[f] = (*magic_results)[f].entropy1;
                magic_ents2[f] = (*magic_results)[f].entropy2;
            } else {
                magic_bfs[f] = magic_probs[f] = magic_means1[f] = magic_means2[f] = magic_ents1[f] = magic_ents2[f] = NA_REAL;
                magic_comps1[f] = magic_comps2[f] = NA_INTEGER;
            }
        }
        
        df_columns["magic_bf_10"] = NumericVector(magic_bfs.begin(), magic_bfs.end());
        df_columns["magic_prob"] = NumericVector(magic_probs.begin(), magic_probs.end());
        df_columns["magic_mean_group1"] = NumericVector(magic_means1.begin(), magic_means1.end());
        df_columns["magic_mean_group2"] = NumericVector(magic_means2.begin(), magic_means2.end());
        df_columns["magic_dom_comp1"] = IntegerVector(magic_comps1.begin(), magic_comps1.end());
        df_columns["magic_dom_comp2"] = IntegerVector(magic_comps2.begin(), magic_comps2.end());
        df_columns["magic_entropy1"] = NumericVector(magic_ents1.begin(), magic_ents1.end());
        df_columns["magic_entropy2"] = NumericVector(magic_ents2.begin(), magic_ents2.end());
    }
    
    // Mixture results
    if (run_mixture && mixture_results) {
        std::vector<double> mixture_means1(n_features), mixture_means2(n_features), mixture_diffs(n_features);
        std::vector<double> mixture_walds(n_features), mixture_ps(n_features);
        
        for (int f = 0; f < n_features; f++) {
            if ((*mixture_results)[f].valid) {
                mixture_means1[f] = (*mixture_results)[f].mean1;
                mixture_means2[f] = (*mixture_results)[f].mean2;
                mixture_diffs[f] = (*mixture_results)[f].difference;
                mixture_walds[f] = (*mixture_results)[f].wald_stat;
                mixture_ps[f] = (*mixture_results)[f].p_value;
            } else {
                mixture_means1[f] = mixture_means2[f] = mixture_diffs[f] = mixture_walds[f] = mixture_ps[f] = NA_REAL;
            }
        }
        
        df_columns["mixture_mean_group1"] = NumericVector(mixture_means1.begin(), mixture_means1.end());
        df_columns["mixture_mean_group2"] = NumericVector(mixture_means2.begin(), mixture_means2.end());
        df_columns["mixture_difference"] = NumericVector(mixture_diffs.begin(), mixture_diffs.end());
        df_columns["mixture_wald_stat"] = NumericVector(mixture_walds.begin(), mixture_walds.end());
        df_columns["mixture_p_value"] = NumericVector(mixture_ps.begin(), mixture_ps.end());
    }
    
    return DataFrame(df_columns);
}
