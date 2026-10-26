# undeRfactor 0.2.0

This is the version cited in version 1 of the preprint *Fitting One Factor to
Two Traits: Population Solutions and Fit Indices*. Theorem, section and
appendix numbers refer to that version.

* `factor_share()` keeps absolute-loading preprocessing by default and adds
  `mode = "signed"` for signed covariance matching. Negative loading inputs
  without an explicit mode now warn; explicit `mode = "absolute"` opts into
  the legacy interpretation.
* Signed coefficient ratios are distinguished from new bounded magnitude
  shares. Zero sums return undefined signed shares, and zero composites return
  undefined signed and magnitude shares, with an explicit status.
* Sum/contrast evaluation avoids unnecessary cancellation near singular trait
  correlations. Continuous condition diagnostics expose remaining sensitivity;
  no near-boundary rejection threshold is added.
* Invalid/nonfinite inputs, zero reference loadings, singular correlations, and
  mismatched vector lengths now fail with informative errors. Scalars still
  recycle against a common vector length. Five original result columns retain
  their order; diagnostic columns are appended. Existing well-defined values
  are unchanged up to floating-point roundoff. Callers that relied on malformed
  inputs, implicit partial recycling, or a fixed column count must adapt.
* `pop_fit()`, its minima summaries, and `order_check()` now pass signed block
  means and consistently oriented reference parameters to the signed mode.
  Positive-loading manuscript results are preserved. Unavailable shares have
  status/reason fields and do not discard usable fits.
* Documentation of `phi_star()` now states that the WLSMV weights use
  numerically differentiated cell probabilities and that uniqueness of the
  WLSMV root is established only in the evaluated designs.
* `pop_fit()` no longer aborts on some inputs, for example the PHQ-ADS layout
  with loading .70 and `phi = .8`. L-BFGS-B could evaluate the ML profile at
  negative uniquenesses from the start `(.1, .9)`; uniquenesses are now
  optimized on the log scale, and a failing start is dropped instead of
  ending the search. Results agree with the previous optimizer to about 1e-8
  wherever it ran.
* For designs with one loading per block and `phi != 0`, `pop_fit()` now takes
  every minimum from the exact quartic of Appendix D (`method = "exact (block
  quartic)"`) rather than from multistart; `exact = FALSE` keeps the numerical
  search. The exchangeable closed forms now also cover `phi < 0`.
* Negative `phi`: minima that differ only in their overall sign are merged, so
  a bipolar composite is no longer reported as two mirror minima. When the
  signed coefficients have opposite signs, the printout reports the magnitude
  share instead of a signed share that can be arbitrarily large; `minima` gains
  `magnitude_share_A`. Signed shares keep the generating orientation, as before.
* `pop_fit()` returns `scores`: the trait composition of the regression or
  Bartlett factor score and of the sum score.
* New `block_solutions()` and `fold_phi()`: all block-uniform stationary points
  of a two-block design (minima and saddles) from the quartic, and the trait
  correlations at which the number of minima changes.
* New `score_share()`: trait shares and trait correlations of factor scores
  and sum scores, defined item by item for heterogeneous loadings.
* New `predict_one_factor()`: from a lavaan two-factor fit or a loading
  matrix, the population one-factor solution from a pooled loading (T1), from
  block summaries (T2, exact) and from the full implied matrix (T3), with the
  Appendix E bound.
* `order_check()` with `estimator = "ULS"` now standardizes continuous items,
  so it fits ULS on correlations (Theorem 2); before, rescaling one item could
  change the solution. It also accepts a covariance matrix with `nobs`, and
  reports `min_uniqueness` for each solution.
* `phi_star(estimator = "WLSMV")` refuses item-specific thresholds with a
  pointer to `pop_fit()`, instead of failing with an uninformative error.

# undeRfactor 0.1.0

* Initial version.
