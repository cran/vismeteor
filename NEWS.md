# vismeteor 2.0.2

## Changes

- Clarified documentation and vignettes for `vmgeom` and `vmideal` models, including improved mathematical notation and explanations.
- Hidden internal parameters from exported function documentation to reduce user confusion.

# vismeteor 2.0.1

## Changes

- Updated README guidance, vignettes, and roxygen references to clarify data sources, models, and to correct documentation errors.
- Improved the performance of `vmgeomVstFromMagn()` and `vmidealVstFromMagn()` by tightening the interpolation steps used during the variance-stabilising transforms.
- Streamlined the `vmtable()` rounding routine to reduce allocations.
- Added a GitHub Actions workflow to run package checks automatically.
- Documented the derivation scripts in `inst/derivation/` to improve reproducibility.
- Added `load_vmdb()` regression tests.

# vismeteor 2.0.0

## Highlights
This release introduces variance-stabilizing transformations for the ideal distribution of visual meteor magnitudes (`vmidealVstFromMagn()`) as well as for visual meteor magnitudes under a geometric distribution (`vmgeomVstFromMagn()`).

## Other changes
The function `vmperception()` now better matches the perception probabilities of *Koschack & Rendtel (1990b)*.  
The argument `deriv.degree` has been removed, as it was only intended for internal testing and had no practical relevance for regular use.

Laplace-transformed perception probabilities have been replaced by variance-stabilizing transformations, which also means that the function `vmperception.l()` has been removed.

**Note:** This release includes breaking changes and is **not fully backward compatible** due to the removal of parameters and functions.

# vismeteor 1.8.5
- Initial CRAN release
