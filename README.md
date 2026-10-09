# undeRfactor

What does a single factor measure when the items come from two correlated
traits, and when do fit indices notice? `undeRfactor` implements the results of
the paper *When one factor represents two traits: symmetry breaking and
population comparative fit* (manuscript in preparation).

| Function | What it returns |
|---|---|
| `phi_star()` | The trait correlation below which the one-factor solution stops being an equal composite and follows one trait, for ML, ULS and WLSMV |
| `pop_fit()` | Exact population discrepancy, CFI, TLI, RMSEA, SRMR, loadings and trait shares of the one-factor model, with every distinct minimum |
| `blind_spot()` | The trait correlations above which conventional cutoffs no longer reject the one-factor model |
| `factor_share()` | The share of each trait in a fitted factor |
| `order_check()` | Whether item order or starting values change a one-factor solution in a sample |

`two_trait_cor()` builds the population correlation matrix of a two-trait design.

## Installation

```r
# install.packages("remotes")
remotes::install_github("atomashevic/underfactor")
```

## Example

Two traits, four items each, all loadings .7.

**Where the factor stops being a composite.** Under maximum likelihood the
threshold is φ\* = .409. Unweighted least squares keeps the composite down to
1/p = .125, and WLSMV lower still.

```r
library(undeRfactor)
phi_star(m = 4, lambda = 0.7)
#> [1] 0.4089404
phi_star(m = 4, lambda = 0.7, estimator = "ULS")
#> [1] 0.125
phi_star(m = 4, lambda = 0.5, estimator = "WLSMV", thresholds = c(-0.8, 0, 0.8))
#> [1] 0.08064138
```

**Population fit.** At a trait correlation of .3 the one-factor ML solution
follows one trait (78.6% of its weight), and its population fit is poor.

```r
pop_fit(two_trait_cor(phi = 0.3, lambda_A = 0.7, m_A = 4))
#> Population one-factor fit (ML, closed form), 8 items
#>   F = 0.9408   CFI = 0.589   TLI = 0.425   RMSEA = 0.217   SRMR = 0.161
#>   The factor follows one trait: share 0.786. The two minima are mirror images,
#>   one per trait, with the same fit.
```

**Blind spots of the usual cutoffs.** Above these trait correlations each
cutoff accepts one factor even with unlimited data (Table 1 of the paper).

```r
blind_spot(lambda_A = 0.7, m_A = 4)
#>       CFI       TLI     RMSEA      SRMR
#> 0.8160628 0.8478085 0.8753186 0.6259122
```

**Item order in a sample.** With slightly unequal loadings (.68 and .72) and
weakly correlated traits, the one-factor discrepancy has two minima. lavaan's
default start reaches the worse one, which follows the weaker trait, whenever
that block is listed first.

```r
set.seed(1)
R <- two_trait_cor(phi = 0.1, lambda_A = 0.68, lambda_B = 0.72)
X <- as.data.frame(matrix(rnorm(2000 * 8), 2000) %*% chol(unclass(R)))
names(X) <- rownames(R)
order_check(X, blocks = list(A = paste0("A", 1:4), B = paste0("B", 1:4)), nstart = 10)
#> One-factor fits: 12 of 12 converged, 2 distinct solutions.
#> Two-factor model: mean loadings 0.68 (A) and 0.72 (B), trait correlation 0.15.
#>  solution   fmin    chisq    cfi  rmsea share_A follows
#>         1 0.4888 1955.086 0.5776 0.2199  0.0677       B
#>         2 0.6175 2470.027 0.4652 0.2475  0.8770       A
#>                     found_by
#>  default (B first), 6 random
#>  default (A first), 4 random
#> Item order changes the solution: with A listed first the default start reaches solution 2, with B first solution 1.
#> A default fit stopped at a local minimum; solution 1 has the lowest discrepancy.
```

## Reproducing the paper

The tests in `tests/testthat` check the package against the paper's numbers:
the thresholds of Theorems 1 and 2 and Corollary 1, Table 1, the closed forms
against numerical minimization, the spurious minimum of Proposition 1, and the
PHQ-ADS shares of Section 8.

## License

MIT
