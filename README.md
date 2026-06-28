# biostat-php

> A pure-PHP biostatistics library implementing descriptive, bivariate and
> multivariate methods for survey-based epidemiological studies — including
> logistic regression, VIF, Box–Tidwell, GLMM, GEE, MICE and Rubin's pooling.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PHP Version](https://img.shields.io/badge/PHP-%E2%89%A58.0-777BB4?logo=php&logoColor=white)](https://www.php.net/)
[![Packagist](https://img.shields.io/badge/packagist-touilelhadj%2Fbiostat--php-orange)](https://packagist.org/packages/touilelhadj/biostat-php)
[![CI](https://github.com/Touil-Elhadj/biostat-php/actions/workflows/ci.yml/badge.svg)](https://github.com/Touil-Elhadj/biostat-php/actions)
[![Tests](https://img.shields.io/badge/tests-PHPUnit-blue)](tests/)
[![Status](https://img.shields.io/badge/status-active-success)](https://github.com/Touil-Elhadj/biostat-php)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21013661.svg)](https://doi.org/10.5281/zenodo.21013661)

---

## 📖 Overview

`biostat-php` brings the analytic side of a cross-sectional epidemiological
study to PHP. It implements, in pure PHP and without external
dependencies, the same family of statistical methods normally available
only in R (`stats`, `car`, `lme4`, `geepack`, `mice`) or SPSS.

The library was extracted from a real research instrument — a trilingual
data-collection platform used to study adolescent overweight in the
Wilaya of Chlef, Algeria — where it ran the 48 pre-registered hypotheses
of the underlying master thesis on a low-cost shared-hosting environment
that did not allow R or Python.

Every method is cross-checked against R 4.x and IBM SPSS 25; the
quantitative comparison is documented in
[`docs/validation-tables.md`](docs/validation-tables.md) and exercised by
the PHPUnit suite in [`tests/`](tests/).

---

## 🚀 Installation

```bash
composer require touilelhadj/biostat-php
```

Requirements:

- PHP **≥ 8.0**
- No external PHP extension beyond the defaults (`mbstring`, `json` only
  for the test fixtures)
- No system library, no Composer dependency at runtime

---

## ⚡ Quick start

```php
<?php
require_once __DIR__ . '/vendor/autoload.php';

use TouilElhadj\BiostatPhp\BiostatAnalysis;

$stats = new BiostatAnalysis();

// Two-by-two table: 45 exposed cases, 18 exposed non-cases,
//                   30 unexposed cases, 22 unexposed non-cases.
$chi = $stats->chi2Test2x2(45, 18, 30, 22);
echo "chi² = {$chi['chi2']}, p = {$chi['p']}\n";
//      chi² = 1.8027, p = 0.1794

$or = $stats->oddsRatio(45, 18, 30, 22);
echo "OR = {$or['or']} [{$or['ci_low']}, {$or['ci_high']}]\n";
//      OR = 1.83 [0.84, 3.98]
```

```php
// Welch's t-test
$t = $stats->tTest(
    [10, 12, 11, 13, 14, 12, 11],
    [15, 17, 16, 18, 16, 17, 15]
);
echo "t = {$t['t']}, df = {$t['df']}, p = {$t['p']}\n";

// Logistic regression (one continuous covariate)
$y = [0,0,0,1,1,1,1,1,0,1,1,0,1,1,1,0,0,1,1,1];
$x = [1,1,2,2,3,3,4,4,1,2,3,4,5,5,4,2,1,3,4,5];
$lr = $stats->logisticRegression($y, $x);
echo "OR = {$lr['or']}, p = {$lr['p']}\n";

// Benjamini–Hochberg FDR adjustment
$adj = BiostatAnalysis::benjaminiHochberg([0.001, 0.01, 0.04, 0.2, 0.5]);
//   [0.005, 0.025, 0.067, 0.25, 0.5]
```

More examples — including GLMM and MICE — are in [`examples/`](examples/).
A complete per-method walkthrough with R-equivalent commands is in
[`docs/usage-examples.md`](docs/usage-examples.md).

---

## 🧪 Method catalogue

| Family | Methods |
|---|---|
| Descriptive | `mean`, `std`, `median`, `quantile` |
| 2 × 2 tables | `chi2Test2x2` (with Yates), `oddsRatio` (with Haldane–Anscombe) |
| Means | `tTest` (Welch), `anova` |
| Correlation | `pearson`, `spearman` |
| Binomial | `binomialTest` |
| Logistic regression | `logisticRegression`, `logisticRegressionMulti`, `hosmerLemeshow` |
| Multiple testing | `benjaminiHochberg` (FDR) |
| Multicollinearity | `vif` |
| Logit linearity | `boxTidwell` |
| Mixed / clustered | `glmmLogistic` (PQL), `geeLogistic` (Liang–Zeger sandwich) |
| Missing data | `mice` (PMM + chained Gibbs), `rubinPool` (Rubin's rules) |

Mathematical formulations and bibliographic references for every method
are in [`docs/statistical-methods.md`](docs/statistical-methods.md).

---

## ✅ Verification

Reference values for every method were computed independently in R 4.3.0
and IBM SPSS 25. The tolerances used in the test suite are:

| Family | Tolerance |
|---|---|
| *p*-values | ± 0.01 |
| Odds ratios | ± 0.01 |
| Correlation coefficients | ± 0.001 |
| Regression coefficients | ± 0.05 |
| Variance components (GLMM) | ± 0.05 |

A complete numerical comparison table — R command, R value, PHP value,
|Δ| — lives in [`docs/validation-tables.md`](docs/validation-tables.md).

---

## 🧰 Running the tests

```bash
composer install
composer test                      # run PHPUnit
composer test:coverage             # with HTML coverage report
composer analyse                   # PHPStan level 6
composer cs                        # PSR-12 check
```

---

## 📁 Repository layout

```
biostat-php/
├── src/
│   ├── BiostatAnalysis.php       main class (≈ 1 460 LOC)
│   ├── LinearAlgebra.php         matMul / transpose / OLS  (trait)
│   ├── Distributions.php         normal / χ² / t / F CDFs  (trait)
│   └── Exceptions/
│       └── ConvergenceException.php
├── tests/                        PHPUnit tests against R reference values
├── docs/
│   ├── statistical-methods.md    formal mathematical specification
│   ├── validation-tables.md      R vs PHP numerical comparison
│   └── usage-examples.md         worked example per public method
├── examples/                     stand-alone runnable scripts
├── composer.json                 PSR-4 autoload + dev tools
├── phpunit.xml                   PHPUnit configuration
├── paper.md                      JOSS paper (≈ 1 000 words)
├── paper.bib                     BibTeX references
├── CITATION.cff                  machine-readable citation
├── LICENSE                       MIT
└── README.md
```

---

## 📚 Citation

If you use this library in a research publication, please cite it. A
machine-readable citation file is provided in
[`CITATION.cff`](CITATION.cff).

### Recommended citation

> TOUIL, E. (2026). *biostat-php: a pure-PHP biostatistics library for
> survey-based epidemiological studies.* Version 1.0.0.
> https://github.com/Touil-Elhadj/biostat-php

### BibTeX

```bibtex
@software{touil_biostat_php_2026,
  author       = {TOUIL, Elhadj},
  title        = {biostat-php: a pure-PHP biostatistics library for
                  survey-based epidemiological studies},
  year         = {2026},
  version      = {1.0.0},
  url          = {https://github.com/Touil-Elhadj/biostat-php}
}
```

---

## 🌍 Real-world use

`biostat-php` is the analytic engine of the
[**chlef-touilelhadj**](https://github.com/Touil-Elhadj/chlef-touilelhadj)
platform, a trilingual web instrument used to enrol 1 220 adolescents
(14–19 years) in the Wilaya of Chlef, Algeria, during the 2025–2026
academic year. All statistical results of the underlying master thesis
were produced by this library and independently verified against R 4.x
and SPSS 25.

---

## 🤝 Contributing

Contributions are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for
the development workflow and the **statistical-contribution checklist**
(every new method must ship with a closed-form reference value from R
or SPSS plus a PHPUnit assertion).

---

## 🛡 Security

To report a vulnerability please follow the procedure in
[`SECURITY.md`](SECURITY.md).

---

## 📄 License

Released under the [MIT License](LICENSE) — © 2026 TOUIL Elhadj.
