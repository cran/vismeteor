# vismeteor 3.0.0

## New features

- `select_knots()` performs forward/backward stepwise selection of spline
  knots from a candidate set, scoring each fit with a user-supplied function
  (e.g. AIC/BIC). Backward selection supports a "bulk removal" mode, and
  scoring can run in parallel via the `parallel` package. See
  `?select_knots` and `vignette("select_knots")`.

## Breaking changes

### Naming overhaul (snake_case)

All public identifiers were migrated from dotted / camelCase to `snake_case` to
match the tidyverse / r-lib convention and to remove visual confusion with S3
dispatch (`print.foo`). This affects exported functions, function parameters,
and `data.frame` columns returned by `load_vmdb_*()` and the example datasets
`PER_2015_rates` / `PER_2015_magn` (which were regenerated).

The d/p/q/r distribution prefixes (`dvmgeom`, `pvmideal`, ...) are unchanged,
and `lower.tail` is preserved to match base R conventions.

Renamed exported functions:

| 2.1.0                | 3.0.0                    |
|----------------------|--------------------------|
| `freq.quantile`      | `freq_quantile`          |
| `vmgeomVstFromMagn`  | `vmgeom_vst_from_magn`   |
| `vmgeomVstToR`       | `vmgeom_vst_to_r`        |
| `vmidealVstFromMagn` | `vmideal_vst_from_magn`  |
| `vmidealVstToPsi`    | `vmideal_vst_to_psi`     |

Renamed parameters: `lim.magn`, `magn.id`, `rate.id`, `session.id`,
`perception.fun`, `sun.alt.max`, `moon.alt.max`, `deriv.degree`, `withSessions`,
`withMagnitudes` → their snake_case forms.

### Closer alignment with the imo-vmdb API

The imo-vmdb API already uses snake_case, so column remapping is reduced to a
few semantic renames: `id` → `rate_id` / `magn_id` / `session_id` (so foreign
keys in the same data frame stay unambiguous) and `mean` → `magn_mean` (to
avoid shadowing base R `mean()`). Three columns that previously carried
R-specific names now pass through unchanged from the API:

| 2.1.0           | 3.0.0       |
|-----------------|-------------|
| `shower.code`   | `shower`    |
| `radiant.alt`   | `rad_alt`   |
| `radiant.az`    | `rad_az`    |

### Other breaking changes

- `perception_fun` now defaults to `vmperception` instead of `NULL` with an
  internal fallback. Callers that explicitly passed `perception.fun = NULL`
  must drop the argument or supply a function.
- `vmperception(m)` is now `vmperception(dm)` — the parameter is the
  difference between the limiting magnitude and the meteor magnitude, and
  the new name reflects that. Positional calls remain compatible.

### Migration

A `sed` sweep of your scripts is usually sufficient:

```sh
sed -i '' -E '
  s/\bfreq\.quantile\b/freq_quantile/g;
  s/\bvmgeomVstFromMagn\b/vmgeom_vst_from_magn/g;
  s/\bvmgeomVstToR\b/vmgeom_vst_to_r/g;
  s/\bvmidealVstFromMagn\b/vmideal_vst_from_magn/g;
  s/\bvmidealVstToPsi\b/vmideal_vst_to_psi/g;
  s/\blim\.magn\b/lim_magn/g;
  s/\bmagn\.id\b/magn_id/g;
  s/\brate\.id\b/rate_id/g;
  s/\bsession\.id\b/session_id/g;
  s/\bperception\.fun\b/perception_fun/g;
  s/\bderiv\.degree\b/deriv_degree/g;
  s/\bwithSessions\b/with_sessions/g;
  s/\bwithMagnitudes\b/with_magnitudes/g;
  s/\bshower\.code\b/shower/g;
  s/\bperiod\.start\b/period_start/g;
  s/\bperiod\.end\b/period_end/g;
  s/\bsl\.start\b/sl_start/g;
  s/\bsl\.end\b/sl_end/g;
  s/\bt\.eff\b/t_eff/g;
  s/\btime\.sidereal\b/sidereal_time/g;
  s/\bsun\.alt\b/sun_alt/g;
  s/\bsun\.az\b/sun_az/g;
  s/\bmoon\.alt\b/moon_alt/g;
  s/\bmoon\.az\b/moon_az/g;
  s/\bmoon\.illum\b/moon_illum/g;
  s/\bfield\.alt\b/field_alt/g;
  s/\bfield\.az\b/field_az/g;
  s/\bradiant\.alt\b/rad_alt/g;
  s/\bradiant\.az\b/rad_az/g;
  s/\bmagn\.mean\b/magn_mean/g;
  s/\blocation\.name\b/location_name/g;
  s/\bobserver\.id\b/observer_id/g;
  s/\bobserver\.name\b/observer_name/g;
' your-script.R
```

# vismeteor 2.1.0

## Breaking changes

- `load_vmdb_rates()` and `load_vmdb_magnitudes()` now connect to an
  [imo-vmdb](https://pypi.org/project/imo-vmdb/) REST API instead of a
  direct database connection. 
- Multi-range filters for `period`, `sl`, and `lim_magn` (previously a
  matrix with multiple rows that were OR-joined) are no longer supported.
  Each filter is collapsed to a single bounding min/max. Use multiple calls
  combined with `rbind()` if disjoint ranges are needed.
- Minimum R version raised to **4.1.0** (was 3.5.0) due to use of the native pipe operator `|>` and `\()`.

# vismeteor 2.0.2

## Changes

- Clarified documentation and vignettes for `vmgeom` and `vmideal` models, including improved mathematical notation and explanations.
- Hidden internal parameters from exported function documentation to reduce user confusion.

# vismeteor 2.0.1

## Changes

- Updated README guidance, vignettes, and roxygen references to clarify data sources, models, and to correct documentation errors.
- Improved the performance of `vmgeom_vst_from_magn()` and `vmideal_vst_from_magn()` by tightening the interpolation steps used during the variance-stabilising transforms.
- Streamlined the `vmtable()` rounding routine to reduce allocations.
- Added a GitHub Actions workflow to run package checks automatically.
- Documented the derivation scripts in `inst/derivation/` to improve reproducibility.
- Added `load_vmdb()` regression tests.

# vismeteor 2.0.0

## Highlights
This release introduces variance-stabilizing transformations for the ideal distribution of visual meteor magnitudes (`vmideal_vst_from_magn()`) as well as for visual meteor magnitudes under a geometric distribution (`vmgeom_vst_from_magn()`).

## Other changes
The function `vmperception()` now better matches the perception probabilities of *Koschack & Rendtel (1990b)*.  
The argument `deriv_degree` has been removed, as it was only intended for internal testing and had no practical relevance for regular use.

Laplace-transformed perception probabilities have been replaced by variance-stabilizing transformations, which also means that the function `vmperception.l()` has been removed.

**Note:** This release includes breaking changes and is **not fully backward compatible** due to the removal of parameters and functions.

# vismeteor 1.8.5
- Initial CRAN release
