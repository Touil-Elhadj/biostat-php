# Changelog

All notable changes to this library are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-05-18

Initial public release.

### Added

- Pure-PHP biostatistics library (`TouilElhadj\BiostatPhp\BiostatAnalysis`)
  implementing twenty-two methods across five statistical families:
  - **Descriptive**: `mean`, `std`, `median`, `quantile`.
  - **Categorical**: `chi2Test2x2` (with optional Yates' correction),
    `oddsRatio` (with Haldane–Anscombe correction for empty cells).
  - **Means**: `tTest` (Welch with Welch–Satterthwaite df), `anova`.
  - **Correlation**: `pearson`, `spearman`.
  - **Binomial**: `binomialTest`.
  - **Logistic**: `logisticRegression`, `logisticRegressionMulti`,
    `hosmerLemeshow`.
  - **Multiple testing**: `benjaminiHochberg` (FDR).
  - **Diagnostics**: `vif`, `boxTidwell`.
  - **Clustered**: `glmmLogistic` (PQL), `geeLogistic` (Liang–Zeger
    sandwich variance).
  - **Missing data**: `mice` (chained Gibbs with Predictive Mean Matching),
    `rubinPool` (Rubin's rules with Barnard–Rubin df and FMI).
- PSR-4 namespaced packaging (`TouilElhadj\BiostatPhp\`) and a fallback
  `autoload.php` for users who prefer not to use Composer.
- Two re-usable traits: `LinearAlgebra` (matMul, matTranspose, matVec,
  matrixInverse, olsRegression) and `Distributions` (normal, χ²,
  Student-*t* and Fisher-*F* CDFs via the regularised incomplete beta
  function).
- `ConvergenceException` for iterative solvers that fail to meet the
  tolerance.
- PHPUnit test suite (40+ assertions across 9 test classes) cross-checked
  against R 4.3 reference values.
- GitHub Actions CI matrix across PHP 8.0, 8.1, 8.2 and 8.3.
- Complete mathematical specification in `docs/statistical-methods.md`
  and a per-method numerical validation table in
  `docs/validation-tables.md`.
- Stand-alone executable examples in `examples/`.
- Project metadata: `LICENSE` (MIT), `README.md`, `CITATION.cff`,
  `paper.md`, `paper.bib`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  `SECURITY.md`.

### Provenance

This library is extracted from the analytical engine of the
[`chlef-touilelhadj`](https://github.com/Touil-Elhadj/chlef-touilelhadj)
research platform (master-thesis project, Hassiba Benbouali University
of Chlef, 2025–2026 academic year). The extraction added namespacing,
PSR-4 autoloading, strict type declarations, PHPUnit tests, the
`LinearAlgebra` and `Distributions` traits, and a fully documented
mathematical specification.

[1.0.0]: https://github.com/Touil-Elhadj/biostat-php/releases/tag/v1.0.0
