---
title: 'biostat-php: a pure-PHP biostatistics library for survey-based
  epidemiological studies'
tags:
  - PHP
  - biostatistics
  - epidemiology
  - logistic-regression
  - GLMM
  - GEE
  - multiple-imputation
  - Rubin-pooling
authors:
  - name: Elhadj TOUIL
    orcid: 0009-0000-2400-459X
    affiliation: 1
affiliations:
  - name: Faculty of Medicine, Hassiba Benbouali University of Chlef (UHBC), Algeria
    index: 1
date: 28 June 2026
bibliography: paper.bib
---

# Summary

`biostat-php` is a self-contained library that implements, in pure PHP,
the family of statistical methods commonly used in cross-sectional and
clustered survey-based epidemiological studies. The library covers
descriptive statistics, χ² and Fisher's exact tests with the odds ratio
and 95 % confidence interval, Welch's *t*-test, one-way ANOVA, Pearson
and Spearman correlations, logistic regression by Newton–Raphson with
the Hosmer–Lemeshow goodness-of-fit test, Benjamini–Hochberg FDR control
[@benjamini1995controlling], Variance Inflation Factor [@allison2012logistic],
the Box–Tidwell test for linearity of the logit [@box1962transformation],
a generalised linear mixed model with random intercept fitted by Penalised
Quasi-Likelihood [@breslow1993approximate], Generalised Estimating
Equations with Liang–Zeger sandwich variance [@liang1986longitudinal], and
Multiple Imputation by Chained Equations [@vanbuuren2011mice] with
inference pooling under Rubin's rules [@rubin1987multiple]. Every routine
is documented with its closed-form mathematical specification and is
cross-checked against R 4.x and IBM SPSS Statistics 25; the reference
values are part of the test suite that ships with the library.

# Statement of need

The R and Python ecosystems contain mature implementations of every
statistical method covered by this library, but they are unavailable in
many real-world deployment environments encountered by epidemiologists
working in low- or middle-income settings. A common pattern is a
PHP/MySQL web platform — installed on a low-cost shared-hosting plan
that disables shell access, Python and R — used to collect questionnaire
data from schools, primary-care clinics or community surveys. The
analytic stage of such projects is therefore conventionally handled by
exporting the data and processing it on a separate workstation running
R or SPSS. The export/re-import workflow is slow, error-prone and
difficult to reproduce, and forces the investigator to maintain two
disjoint environments.

`biostat-php` removes this gap by providing, in a single language and
without any external dependency, the analytic operations needed to take
a survey from raw responses to publication-ready estimates. The library
is intentionally minimal: it does not provide a domain-specific
front-end and makes no assumption about the structure of the data
beyond the row-oriented format returned by `PDOStatement::fetchAll()`.
It is designed to be embedded both in interactive PHP web pages (live
analytical dashboards) and in command-line scripts (batch analyses,
nightly reports). To the best of our knowledge, no comparable library
exists in the PHP ecosystem; the major statistical packages bundled
with REDCap [@harris2009redcap] and LimeSurvey, for instance, only
provide descriptive and frequency statistics and delegate inferential
analysis to external software.

# Functionality

The single public class `BiostatAnalysis` exposes twenty-two methods
grouped in five families: descriptive (`mean`, `std`, `median`,
`quantile`); two-by-two tables (`chi2Test2x2` with Yates' correction,
`oddsRatio` with Haldane–Anscombe correction for empty cells);
comparison of means (`tTest` with Welch–Satterthwaite degrees of
freedom, `anova` with one-way *F* test); correlation (`pearson`,
`spearman`); logistic regression (`logisticRegression`,
`logisticRegressionMulti`, `hosmerLemeshow`); multiple-testing control
(`benjaminiHochberg`); regression diagnostics (`vif`, `boxTidwell`);
clustered designs (`glmmLogistic`, `geeLogistic`); and missing data
(`mice`, `rubinPool`).

Linear-algebra primitives (matrix multiplication, transposition,
Gauss–Jordan inverse, ordinary least squares by normal equations) and
probability-distribution CDFs (standard normal by Zelen–Severo, χ² by
Wilson–Hilferty, Student-*t* and Fisher-*F* by the regularised
incomplete beta function) are factored out as two reusable traits
`LinearAlgebra` and `Distributions`. Both traits are protected by
contract and never appear in the public API, but they are documented
in `docs/architecture.md` for the benefit of downstream maintainers.

The advanced methods deserve a brief comment on their numerical
strategy. `glmmLogistic` fits the random-intercept logistic mixed model
by Penalised Quasi-Likelihood, solving Henderson's mixed-model
equations [@henderson1953estimation] at each iteration and updating
the variance component σ²_u by restricted ML on the residuals of the
random effects; the intra-cluster correlation on the latent scale is
reported as $\sigma^2_u / (\sigma^2_u + \pi^2/3)$. PQL is known to
under-estimate σ²_u relative to the Laplace approximation used by
`lme4::glmer` when cluster sizes are small and the response is rare;
the library reports the convergence flag and the within-cluster
variance estimate so that the user can detect this regime and switch
to a Laplace-based alternative if necessary. `geeLogistic` follows the
Liang–Zeger framework with an exchangeable working correlation matrix
and reports both the model-based and the robust sandwich standard
errors so that the user can assess the adequacy of the working
correlation. `mice` performs $m$ imputations of mixed-type data by
chained equations with Predictive Mean Matching for continuous
variables and proportional-odds / multinomial logit for categorical
variables; `rubinPool` aggregates the $m$ estimates using Rubin's
total variance $T = U + (1 + 1/m) B$, the Barnard–Rubin adjusted
degrees of freedom and the fraction of missing information.

# Quality assurance

A PHPUnit test suite (37 assertions across nine test classes) checks
every public method against a reference value computed in R 4.3 or
IBM SPSS 25. The numerical comparison — R command, R value, PHP value,
absolute difference, tolerance — is documented in
`docs/validation-tables.md` so that a reviewer can independently
reproduce each assertion via the R script `tests/fixtures/reference-values.R`
that is shipped with the library. The tolerance budget is ± 0.001 on
χ² statistics and unadjusted *p*-values, ± 0.01 on odds ratios, ±
0.001 on correlation coefficients and on Benjamini-Hochberg adjusted
*p*-values, ± 0.01 on logistic-regression coefficients, and ± 0.05 on
GLMM variance components. As of release 1.0.1 every assertion passes
within these tolerances against R 4.3.0. A GitHub Actions workflow
runs the suite on PHP 8.0, 8.1, 8.2 and 8.3 on every push.

# Limitations and intended use

`biostat-php` is targeted at the typical sample sizes of a
cross-sectional school or community survey (n in the hundreds to low
thousands). The linear-algebra implementation stores matrices as
arrays of arrays and is not competitive with R, Python/NumPy or C++
on much larger problems. For datasets exceeding 10⁵ observations or
for designs with deeply nested random effects, the user is advised to
fall back to the corresponding R packages. PQL convergence is also
known to be slow and possibly biased for weakly-clustered data with
rare outcomes; this caveat is reported by the routine itself.

# Real-world use

The library was extracted from the analytic engine of the
`chlef-data-analysis` platform, a trilingual web instrument used to
enrol 1 220 adolescents (14–19 years) in the Wilaya of Chlef, Algeria,
during the 2025–2026 academic year. All statistical results of the
underlying master thesis at the Hassiba Benbouali University of Chlef
were produced by `biostat-php` and independently verified against R
4.3 and IBM SPSS 25 before publication.

# Acknowledgements

The author thanks Dr Ali Haimoud S. (UHBC, Faculty of Medicine) for
the supervision of the master thesis from which this library was
extracted.

# References
