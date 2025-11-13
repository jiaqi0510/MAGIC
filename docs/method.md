# Mathematical Method

MAGIC models DNA methylation using beta-binomial mixture models with Bayesian differential testing.

## Beta-Binomial Mixture Model

### Model Specification

For CpG site *i* with methylation counts *m<sub>i</sub>* out of *n<sub>i</sub>* total reads:

```
P(m_i | n_i, Θ) = Σ_{k=1}^K π_k · BetaBinomial(m_i | n_i, α_k, β_k)
```

Where:
- *K* is the number of mixture components
- *π<sub>k</sub>* are mixing proportions with Σ*π<sub>k</sub>* = 1
- *α<sub>k</sub>*, *β<sub>k</sub>* > 0 are beta-binomial shape parameters for component *k*

### Beta-Binomial Distribution

```
BetaBinomial(m | n, α, β) = C(n,m) · B(m + α, n - m + β) / B(α, β)
```

Where:
- *C(n,m)* is the binomial coefficient
- *B(·,·)* is the beta function

### Parameterization

**Natural Space:**
- Shape parameters: *α<sub>k</sub>*, *β<sub>k</sub>*
- Mean methylation: *μ<sub>k</sub>* = *α<sub>k</sub>* / (*α<sub>k</sub>* + *β<sub>k</sub>*)
- Dispersion: *φ<sub>k</sub>* = 1 / (*α<sub>k</sub>* + *β<sub>k</sub>* + 1)

**Log Space (for optimization):**
- log(*α<sub>k</sub>*), log(*β<sub>k</sub>*), log(*π<sub>k</sub>* / *π<sub>K</sub>*)
- Ensures positivity constraints
- Improves numerical stability

### Properties

**Mean:**
```
E[m_i / n_i] = Σ_k π_k · μ_k
```

**Variance:**
```
Var[m_i / n_i] = Σ_k π_k · (μ_k(1-μ_k) / (n_i + 1/φ_k)) + Σ_k π_k(μ_k - μ)²
```

The mixture captures:
- Within-component variance (beta-binomial overdispersion)
- Between-component variance (biological heterogeneity)

## Parameter Estimation

### Maximum Likelihood

Estimate parameters *Θ* = {*α*, *β*, *π*} by maximizing log-likelihood:

```
ℓ(Θ) = Σ_{i=1}^N log[Σ_{k=1}^K π_k · BetaBinomial(m_i | n_i, α_k, β_k)]
```

Over *N* CpG sites.

### Optimization

**Algorithm:** L-BFGS-B quasi-Newton method

**Parameter Space:** Log-transformed for unconstrained optimization

**Initialization:** 
- Multiple random starts for global optimization
- Mixture proportions: Dirichlet(1,...,1)
- Shape parameters: Sample from prior distribution based on data

**Convergence:**
- Relative change in objective function < tolerance
- Gradient norm check
- Maximum iterations

### Gradient Computation

Analytic gradients computed via chain rule:

```
∂ℓ/∂θ = Σ_i Σ_k [π_k · BB(m_i|n_i,α_k,β_k) / Σ_j π_j · BB(m_i|n_i,α_j,β_j)] · ∂log BB/∂θ
```

Using:
- Digamma functions: ψ(*x*) = d/d*x* log Γ(*x*)
- Beta-binomial log-density derivatives

## Model Selection

### Bayesian Information Criterion

```
BIC = -2·ℓ(Θ̂) + p·log(N)
```

Where:
- *ℓ(Θ̂)* is maximized log-likelihood
- *p* = 3*K* - 1 is number of free parameters
- *N* is number of CpG sites

Lower BIC indicates better model fit with penalty for complexity.

### Akaike Information Criterion

```
AIC = -2·ℓ(Θ̂) + 2p
```

Alternative to BIC with different penalty term.

## Differential Methylation Testing

### Bayesian Framework

For each CpG site, test difference in mean methylation between groups.

### Posterior Component Assignment

For group *g* and site *i*, posterior probability of component *k*:

```
γ_ik^(g) = π_k · BB(m_i^(g) | n_i^(g), α_k, β_k) / Σ_j π_j · BB(m_i^(g) | n_i^(g), α_j, β_j)
```

### Bayesian Dispersion Estimation

**Prior:** Mixture-informed log-normal prior on dispersion

For each site, estimate dispersion *φ* via penalized likelihood:

```
φ̂ = argmax_φ {Σ_s log BB(m_s|n_s,μ,φ) + log N(log φ | m₀, τ²)}
```

Where:
- *s* indexes samples within site
- *μ* is empirical mean methylation
- (*m₀*, *τ²*) are global prior parameters estimated from data

**Global Prior Estimation:**

From genome-wide data:
1. Estimate site-specific dispersions
2. Fit log-normal distribution to dispersion estimates
3. Use as prior for individual sites

### Test Statistics

**Wald Test:**

For mean difference *Δμ* = *μ*<sub>1</sub> - *μ*<sub>2</sub>:

```
W = Δμ̂ / SE(Δμ̂)
```

Where standard error accounts for:
- Beta-binomial variance within groups
- Sample sizes
- Estimated dispersions

**Bayes Factor:**

```
BF₁₀ = P(D | H₁) / P(D | H₀)
```

Where:
- *H₀*: no difference between groups
- *H₁*: difference exists
- Computed via Laplace approximation

Convert to probability:

```
P(H₁ | D) = BF₁₀ / (1 + BF₁₀)
```

Assuming equal prior odds.

## Numerical Implementation

### Log-Space Calculations

All likelihood calculations in log space to prevent:
- Underflow with small probabilities
- Overflow with large factorials

### Log-Sum-Exp Trick

For numerical stability:

```
log(Σ_k exp(x_k)) = x_max + log(Σ_k exp(x_k - x_max))
```

Where *x_max* = max{*x<sub>k</sub>*}

### Special Functions

**Log-Gamma:** Stirling's approximation for large arguments

**Digamma:** Asymptotic expansion for *x* > 8, recursion for smaller *x*

**Log-Beta:** Computed as log Γ(*a*) + log Γ(*b*) - log Γ(*a*+*b*)

### Parallelization

**OpenMP:** Parallel processing of sites

Sites are independent conditional on parameters, enabling:
- Parallel likelihood evaluation
- Parallel gradient computation
- Parallel differential testing

## Convergence and Stability

### Multiple Initializations

Run optimization from multiple starting points:
- Different random seeds
- Select best fit (lowest BIC)
- Assess convergence stability via coefficient of variation

### Component Alignment

For multi-run comparison, align components by mean methylation:
1. Compute mean methylation per component
2. Sort components in ascending order
3. Compare aligned parameters across runs

### Identifiability

Mixture models have label-switching non-identifiability:
- Component labels arbitrary
- Alignment required for comparison
- Resolved by ordering constraint

## Statistical Properties

### Consistency

Under regularity conditions:
- Parameter estimates converge to true values as *N* → ∞
- Rate: *O(N<sup>-1/2</sup>)*

### Asymptotic Normality

For large *N*:
```
√N(Θ̂ - Θ) →ᵈ N(0, I⁻¹(Θ))
```

Where *I(Θ)* is Fisher information matrix.

### Model Misspecification

If true model has *K'* > *K* components:
- Fitted model approximates true mixture
- BIC favors true *K'* for sufficient *N*

If *K'* < *K*:
- Extra components have small mixing proportions
- BIC penalizes unnecessary complexity

## Computational Complexity

**Optimization:**
- Per iteration: *O(NKS)* where *S* is samples
- Typical convergence: 50-200 iterations
- Multiple runs: multiply by number of initializations

**Testing:**
- Per site: *O(KS)*
- Total: *O(NKS)*
- Parallelizes linearly with available cores

**Memory:**
- Storage: *O(NS)*
- Working memory: *O(K)* per thread

## References

**Beta-Binomial Distribution:**
- Johnson, Kotz & Kemp (1992). *Univariate Discrete Distributions*

**EM Algorithm for Mixtures:**
- Dempster, Laird & Rubin (1977). *JRSS-B*

**Bayesian Dispersion Estimation:**
- DSS method: Feng, Conneely & Wu (2014). *Biostatistics*

**Numerical Optimization:**
- Nocedal & Wright (2006). *Numerical Optimization*
