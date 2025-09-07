# vismeteor 2.0.0

## Highlights
This release introduces variance-stabilizing transformations for ideally distributed visual meteor magnitudes 
(`vmidealVstFromMagn()`) as well as for visual meteor magnitudes under a geometric distribution (`vmgeomVstFromMagn()`).

## Highlights
This release introduces variance-stabilizing transformations for the ideal distribution of visual meteor magnitudes (`vmidealVstFromMagn()`) as well as for visual meteor magnitudes under a geometric distribution (`vmgeomVstFromMagn()`).

## Other changes
The function `vmperception()` now better matches the perception probabilities of *Koschack & Rendtel (1990b)*.  
The argument `deriv.degree` has been removed, as it was only intended for internal testing and had no practical relevance for regular use.

Laplace-transformed perception probabilities have been replaced by variance-stabilizing transformations, which also means that the function `vmperception.l()` has been removed.

**Note:** This release includes breaking changes and is **not fully backward compatible** due to the removal of parameters and functions.

# vismeteor 1.8.5
- Initial CRAN release
