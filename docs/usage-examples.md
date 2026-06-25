# Usage examples

This document gives a **complete, runnable** example for every public
method of `BiostatAnalysis`. Each section shows:

- the **call signature**,
- the **input** with a concrete example,
- the **output** (actual values produced),
- the **R equivalent** so you can verify independently,
- a **note** on when to use it and any caveats.

> All examples assume:
> ```php
> require_once 'vendor/autoload.php';
> use TouilElhadj\BiostatPhp\BiostatAnalysis;
> $stats = new BiostatAnalysis();
> ```

---

## Table of contents

1. [Descriptive statistics](#1-descriptive-statistics)
2. [2 × 2 contingency tables](#2-22-contingency-tables)
3. [Means comparison](#3-means-comparison)
4. [Correlation](#4-correlation)
5. [Binomial test](#5-binomial-test)
6. [Logistic regression](#6-logistic-regression)
7. [Goodness of fit](#7-hosmerlemeshow-goodness-of-fit)
8. [Multiple testing](#8-benjaminihochberg-fdr-adjustment)
9. [Multicollinearity](#9-variance-inflation-factor-vif)
10. [Logit linearity](#10-boxtidwell-test)
11. [Mixed model](#11-glmmlogistic-random-intercept-logistic-glmm)
12. [GEE](#12-geelogistic-generalised-estimating-equations)
13. [Multiple imputation](#13-mice-multiple-imputation-by-chained-equations)
14. [Rubin pooling](#14-rubinpool-rubins-rules)

---

## 1. Descriptive statistics

### `mean($array)`, `std($array)`, `median($array)`, `quantile($array, $q)`

```php
$x = [2.3, 4.1, 5.7, 6.2, 7.8, 8.4, 9.1, 10.5, 11.2, 12.9];

$stats->mean($x);                  // 7.82
$stats->std($x);                   // 3.296732  (n-1 denominator, like R sd())
$stats->median($x);                // 8.10
$stats->quantile($x, 0.25);        // 5.825     (R type=7 default)
$stats->quantile($x, 0.75);        // 10.150
```

**R equivalent**:
```r
mean(x); sd(x); median(x); quantile(x, c(0.25, 0.75))
```

**Notes**: `std()` returns the sample standard deviation
(n − 1 denominator) to match R's `sd()`. `quantile()` uses linear
interpolation matching R's default `type = 7`. Non-numeric values are
silently filtered out; an empty array returns 0.

---

## 2. 2 × 2 contingency tables

### `chi2Test2x2($a, $b, $c, $d, $yates = true)`

```
            Exposed   Not exposed
  Cases       a            c
  Controls    b            d
```

```php
// 45 exposed cases, 18 exposed controls,
// 30 unexposed cases, 22 unexposed controls
$r = $stats->chi2Test2x2(45, 18, 30, 22, false);   // no Yates
// → ['chi2' => 2.3695, 'p' => 0.1237, 'df' => 1, 'significant' => false]

$r = $stats->chi2Test2x2(45, 18, 30, 22, true);    // with Yates
// → ['chi2' => 1.8027, 'p' => 0.1794, 'df' => 1, 'significant' => false]
```

**R equivalent**:
```r
m <- matrix(c(45, 18, 30, 22), nrow = 2, byrow = TRUE)
chisq.test(m, correct = FALSE)$statistic   # 2.3695
chisq.test(m, correct = TRUE)$statistic    # 1.8027
```

**Notes**: Yates' correction (`$yates = true`) is the default and
matches R's `chisq.test(..., correct = TRUE)`. Set `$yates = false`
for the uncorrected statistic.

---

### `oddsRatio($a, $b, $c, $d)`

```php
$r = $stats->oddsRatio(45, 18, 30, 22);
// → ['or' => 1.83, 'ci_low' => 0.84, 'ci_high' => 3.98,
//    'correction_applied' => false]

// With a zero cell — Haldane-Anscombe correction is automatic:
$r = $stats->oddsRatio(10, 0, 5, 10);
// → ['or' => 40.09, 'ci_low' => 1.96, 'ci_high' => 820.56,
//    'correction_applied' => true]
```

**R equivalent**:
```r
epitools::oddsratio.wald(m)$measure
#                     estimate     lower    upper
#  Exposed = 1            1.000        NA       NA
#  Exposed = 2            1.833     0.844    3.984
```

**Notes**: when any of the four cells is zero, +0.5 is added to every
cell (Haldane–Anscombe). The flag `correction_applied` tells you
whether the correction was triggered.

---

## 3. Means comparison

### `tTest($group1, $group2)`  — Welch's two-sample *t*-test

```php
$a = [10, 12, 11, 13, 14, 12, 11];
$b = [15, 17, 16, 18, 16, 17, 15];

$r = $stats->tTest($a, $b);
// → ['t' => -6.7117, 'df' => 11.6, 'p' => 0.0001,
//    'm1' => 11.86, 'm2' => 16.29, ...]
```

**R equivalent**:
```r
t.test(a, b)   # t = -6.7117, df = 11.6, p = 2.71e-05
```

**Notes**: Welch's correction relaxes the equal-variance assumption,
giving non-integer degrees of freedom via the Welch–Satterthwaite
formula. Use this whenever variances may differ — it reduces to
Student's *t* when variances are equal.

---

### `anova($groups)` — one-way ANOVA

```php
$r = $stats->anova([
    [10, 12, 11],
    [15, 17, 16],
    [20, 22, 21],
]);
// → ['F' => 75, 'p' => 0.0001, 'dfB' => 2, 'dfW' => 6, ...]
```

**R equivalent**:
```r
g1 <- c(10, 12, 11); g2 <- c(15, 17, 16); g3 <- c(20, 22, 21)
summary(aov(values ~ ind, data = stack(list(g1=g1, g2=g2, g3=g3))))
# F(2, 6) = 75, p = 6.34e-05
```

**Notes**: input is an array of groups, not a long-format design
matrix. The library does not run post-hoc pairwise comparisons; the
caller should do these separately (e.g. pairwise Welch t-tests with
Benjamini–Hochberg adjustment).

---

## 4. Correlation

### `pearson($x, $y)`

```php
$x = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
$y = [15, 18, 17, 22, 24, 26, 25, 29, 32, 30];

$r = $stats->pearson($x, $y);
// → ['r' => 0.967, 'p' => 0.0000, 'n' => 10, 't' => 10.69, ...]
```

**R equivalent**:
```r
cor.test(x, y, method = "pearson")   # r = 0.9668, p = 5.16e-06
```

---

### `spearman($x, $y)`

```php
$r = $stats->spearman($x, $y);
// → ['r' => 0.976, 'p' => 0.0000, 'n' => 10, ...]
```

**R equivalent**:
```r
cor.test(x, y, method = "spearman")   # rho ≈ 0.976
```

**Notes**: Spearman's ρ is Pearson's *r* of the mid-ranks. It is
robust to non-linear monotonic relationships and to outliers. The
implementation handles ties via mid-rank averaging.

---

## 5. Binomial test

### `binomialTest($successes, $n, $p0 = 0.5)`

```php
// Out of 250 surveyed teens, 65 reported daily soda consumption.
// Test if the prevalence is significantly above 20 %.
$r = $stats->binomialTest(65, 250, 0.20);
// → ['p_obs' => 0.26, 'p0' => 0.20, 'z' => 2.37, 'p_value' => 0.018, ...]
```

**R equivalent**:
```r
binom.test(65, 250, p = 0.20)
prop.test(65, 250, p = 0.20, correct = FALSE)
```

---

## 6. Logistic regression

### `logisticRegression($y, $x)` — single continuous predictor

```php
$y = [0,0,0,0,0, 0,0,1,0,1, 0,1,0,1,1, 1,1,0,1,1,
      1,0,1,1,1, 1,1,1,1,1];
$x = [1,2,2,3,3, 4,4,4,5,5, 5,6,6,6,7, 7,7,8,8,8,
      9,9,9,10,10, 10,11,11,12,12];

$r = $stats->logisticRegression($y, $x);
// → ['intercept' => -3.9015,
//    'coef'      => 0.6814,
//    'or'        => 1.98,
//    'ci_low'    => 1.69, 'ci_high' => 2.31,
//    'p'         => 0.0000,
//    'significant' => true]
```

**R equivalent**:
```r
m <- glm(y ~ x, family = binomial(link = "logit"))
summary(m)
# (Intercept) -3.9015
# x            0.6814
exp(coef(m))   # 1.9767
```

---

### `logisticRegressionMulti($y, $X, $names = [])` — multiple predictors

```php
// y = obesity (0/1)
// X = n × k matrix; intercept will be added automatically
$X = [];
foreach ($subjects as $s) {
    $X[] = [
        $s['age'],
        $s['screen_hours'],
        $s['daily_soda'] ? 1 : 0,
    ];
    $y[] = $s['obese'] ? 1 : 0;
}

$r = $stats->logisticRegressionMulti($y, $X, ['age', 'screen', 'soda']);
// → [
//      'coef'      => [β0, β_age, β_screen, β_soda],
//      'se'        => [...],
//      'or'        => [...], 'ci_low' => [...], 'ci_high' => [...],
//      'p'         => [...],
//      'aic'       => 412.3,
//      'auc'       => 0.748,
//      'hl'        => ['chi2' => 4.21, 'df' => 8, 'p' => 0.838],
//      'converged' => true,
//      'iter'      => 6,
//      'predicted' => [...],
//    ]
```

**R equivalent**:
```r
m <- glm(obese ~ age + screen + soda,
         data = d, family = binomial)
summary(m)
ResourceSelection::hoslem.test(d$obese, fitted(m), g = 10)
pROC::auc(d$obese, fitted(m))
```

**Notes**:
- Intercept is added automatically — do **not** include a constant
  column in `$X`.
- AUC is the Mann–Whitney *U* statistic on the fitted probabilities.
- AIC = −2 · log-lik + 2(k + 1).

---

## 7. Hosmer–Lemeshow goodness-of-fit

### `hosmerLemeshow($y, $p, $g = 10)`

```php
// $p must be a vector of predicted probabilities,
// typically from logisticRegressionMulti()['predicted'].
$hl = $stats->hosmerLemeshow($y, $p, 10);
// → ['chi2' => 4.21, 'df' => 8, 'p' => 0.838,
//    'groups' => [...details per decile...]]
```

**R equivalent**:
```r
ResourceSelection::hoslem.test(y, p_hat, g = 10)
```

**Interpretation**: a *large* p-value means **no evidence of misfit**
(good calibration). A small p-value (< 0.05) suggests the model
predictions deviate from observed frequencies within at least one
decile.

---

## 8. Benjamini–Hochberg FDR adjustment

### `BiostatAnalysis::benjaminiHochberg($pvalues)`  (static)

```php
$p = [
    'sex'    => 0.001,
    'age'    => 0.008,
    'income' => 0.039,
    'region' => 0.042,
    'screen' => 0.205,
];

$adj = BiostatAnalysis::benjaminiHochberg($p);
// → ['sex' => 0.005, 'age' => 0.02, 'income' => 0.0525,
//    'region' => 0.0525, 'screen' => 0.205]
```

**R equivalent**:
```r
p.adjust(p, method = "BH")
```

**Notes**: keys are preserved, so you can keep your variable names
attached to their adjusted *p*-values throughout the analysis.
Monotonicity is enforced (the cumulative minimum step from Benjamini
& Hochberg 1995).

---

## 9. Variance Inflation Factor (VIF)

### `vif($X, $names = [])`

```php
$X = [
    [25, 70, 1.75, 22.9],
    [13, 45, 1.50, 20.0],
    // ... one row per subject
];

$r = $stats->vif($X, ['age', 'weight', 'height', 'bmi']);
// → ['age' =>    ['vif' => 1.12, 'r2' => 0.10, 'tolerance' => 0.90, 'flag' => 'OK'],
//    'weight' => ['vif' => 8.34, 'r2' => 0.88, 'tolerance' => 0.12, 'flag' => 'PROBLEMATIC'],
//    'height' => ['vif' => 3.21, 'r2' => 0.69, 'tolerance' => 0.31, 'flag' => 'NOTICEABLE'],
//    'bmi'    => ['vif' => 9.10, 'r2' => 0.89, 'tolerance' => 0.11, 'flag' => 'PROBLEMATIC']]
```

**R equivalent**:
```r
car::vif(lm(y ~ age + weight + height + bmi, data = d))
```

**Interpretation**:
- VIF < 2.5: no concern.
- 2.5 ≤ VIF < 5: noticeable.
- 5 ≤ VIF < 10: problematic.
- VIF ≥ 10: severe multicollinearity — consider dropping or combining
  predictors.

---

## 10. Box–Tidwell test

### `boxTidwell($y, $X, $continuous_idx, $names = [])`

```php
// continuous_idx tells which columns of X are continuous variables
// that should be tested for log-linearity of the logit.
$r = $stats->boxTidwell(
    $obese,
    $X,                          // [age, sex, bmi]
    [0, 2],                      // age and bmi are continuous
    ['age', 'sex', 'bmi']
);
// → ['age' => ['p' => 0.421, 'linear' => true],
//    'bmi' => ['p' => 0.018, 'linear' => false]]
```

**R equivalent**:
```r
car::boxTidwell(y ~ age + bmi, other.x = ~ sex, data = d)
```

**Interpretation**: a small p-value on the X · log(X) interaction
suggests the logit is **not linear** in that continuous predictor —
consider a transformation (log, square, splines).

---

## 11. `glmmLogistic` — random-intercept logistic GLMM

### Signature

```php
glmmLogistic($y, $X, $cluster, $names = [], $max_iter = 50)
```

```php
// y          : n binary outcomes
// X          : n × k predictor matrix (no intercept column)
// cluster    : n integers identifying the cluster of each obs
// names      : labels of the k predictors (optional)

$r = $stats->glmmLogistic($obese, $X, $school_id, ['age', 'screen'], 50);
// → [
//      'coef'      => [β_intercept, β_age, β_screen],
//      'se'        => [...],
//      'or'        => [...], 'ci_low' => [...], 'ci_high' => [...],
//      'p'         => [...],
//      'sigma2_u'  => 0.342,    // between-school variance
//      'icc'       => 0.094,    // σ²_u / (σ²_u + π²/3)
//      'converged' => true,
//      'iter'      => 8,
//    ]
```

**R equivalent**:
```r
m <- lme4::glmer(obese ~ age + screen + (1 | school),
                 family = binomial, data = d)
summary(m)
performance::icc(m)
```

**Caveats**:
- Estimation is by **PQL** (Breslow & Clayton 1993). This is known to
  underestimate σ²_u relative to `lme4::glmer` (which uses Laplace)
  when cluster sizes are small or the outcome is rare. For studies
  with > 50 clusters and prevalence between 10 % and 90 %, the bias
  is usually < 10 %.
- The routine returns `converged = false` if it does not meet the
  tolerance within `max_iter`; the partial estimate is still returned.

---

## 12. `geeLogistic` — Generalised Estimating Equations

### Signature

```php
geeLogistic($y, $X, $cluster, $names = [], $max_iter = 30)
```

```php
$r = $stats->geeLogistic($obese, $X, $school_id, ['age', 'screen']);
// → [
//      'coef'      => [β_intercept, β_age, β_screen],
//      'se_model'  => [...],   // naive (sandwich denominator only)
//      'se_robust' => [...],   // Liang–Zeger sandwich
//      'or'        => [...], 'ci_low' => [...], 'ci_high' => [...],
//      'p'         => [...],
//      'alpha'     => 0.082,   // exchangeable working correlation
//      'converged' => true,
//      'n_clusters' => 38,
//    ]
```

**R equivalent**:
```r
m <- geepack::geeglm(obese ~ age + screen,
                     id = school, family = binomial,
                     corstr = "exchangeable", data = d)
summary(m)   # reports the robust SE by default
```

**When to use vs. GLMM**:
- **GEE** is appropriate when interest is in **population-averaged**
  effects (marginal model). The robust SE remains valid even if the
  working correlation is mis-specified.
- **GLMM** gives **cluster-specific** effects on the latent scale,
  and decomposes total variance into within- and between-cluster
  components (ICC).

---

## 13. `mice` — Multiple Imputation by Chained Equations

### Signature

```php
mice($data, $var_types, $m = 20, $max_iter = 20, $donors = 5)
```

```php
$data = [
    ['age' => 14, 'bmi' => 22.1, 'screen' => '>4h', 'sex' => 'F'],
    ['age' => 15, 'bmi' => null, 'screen' => '<2h', 'sex' => 'M'],
    // ... 1218 more rows, some with missing values
];

$types = [
    'age'    => 'continuous',
    'bmi'    => 'continuous',
    'screen' => 'ordered',     // ordered factor
    'sex'    => 'factor',      // unordered factor
];

$imputed = $stats->mice($data, $types, m: 20, max_iter: 20, donors: 5);
// → array of 20 complete datasets, each the same shape as $data
```

**R equivalent**:
```r
imp <- mice::mice(d, m = 20, maxit = 20, method = c(
    age    = "pmm",
    bmi    = "pmm",
    screen = "polr",
    sex    = "polyreg"
))
mice::complete(imp, action = "long")
```

**Notes**:
- Continuous variables use **Predictive Mean Matching** with `donors`
  nearest donors (default 5).
- Ordered factors use proportional-odds logistic regression.
- Unordered factors use multinomial logistic regression.
- Default `m = 20` follows the recommendation of Bodner (2008) for
  studies with up to 50 % missingness.

---

## 14. `rubinPool` — Rubin's rules for pooling

### Signature

```php
rubinPool($estimates, $standard_err)
```

After fitting the same model on each of the `m` imputed datasets,
collect the point estimates and standard errors and pool them:

```php
// 20 imputations × 4 parameters (intercept, age, screen, soda)
$estimates = [
    [-3.21, 0.142, 1.045, 0.872],
    [-3.18, 0.139, 1.052, 0.881],
    // ... 18 more rows
];
$standard_err = [
    [0.420, 0.012, 0.092, 0.084],
    [0.418, 0.012, 0.091, 0.083],
    // ... 18 more rows
];

$pooled = $stats->rubinPool($estimates, $standard_err);
// → [
//      'm'         => 20,
//      'beta'      => [-3.19, 0.140, 1.049, 0.876],
//      'se'        => [...],         // sqrt(T) where T = U + (1+1/m)B
//      'or'        => [...],
//      'ci_low'    => [...], 'ci_high' => [...],
//      'p'         => [...],
//      'df'        => [...],         // Barnard-Rubin adjusted df
//      'within_U'  => [...],
//      'between_B' => [...],
//      'total_T'   => [...],
//      'fmi'       => [...],         // fraction of missing information
//    ]
```

**R equivalent**:
```r
fits <- with(imp, glm(obese ~ age + screen + soda, family = binomial))
pool(fits)
```

---

## Putting it together — a realistic workflow

Below is the analytic pipeline of a cross-sectional study (n = 1 220
adolescents, clustered by school, with some missing data):

```php
require 'vendor/autoload.php';
use TouilElhadj\BiostatPhp\BiostatAnalysis;

$stats = new BiostatAnalysis();

// 1. Multiple imputation of the missing data
$imputed = $stats->mice($data, $var_types, m: 20);

// 2. For each imputed dataset, fit the adjusted logistic GLMM
$est = []; $se = [];
foreach ($imputed as $dataset) {
    $y       = array_column($dataset, 'is_overweight');
    $X       = build_design_matrix($dataset, ['age', 'screen', 'soda', 'sex']);
    $cluster = array_column($dataset, 'school_id');

    $fit = $stats->glmmLogistic($y, $X, $cluster, ['age','screen','soda','sex']);

    $est[] = $fit['coef'];
    $se[]  = $fit['se'];
}

// 3. Pool the m fits via Rubin's rules
$pooled = $stats->rubinPool($est, $se);

// 4. Adjust the four predictor p-values for multiple testing
$p_raw = ['age' => $pooled['p'][1],
          'screen' => $pooled['p'][2],
          'soda' => $pooled['p'][3],
          'sex' => $pooled['p'][4]];

$p_adj = BiostatAnalysis::benjaminiHochberg($p_raw);

// 5. Report
foreach (['age', 'screen', 'soda', 'sex'] as $i => $name) {
    printf("%-7s  OR = %.2f [%.2f, %.2f]  raw p = %.4f  BH p = %.4f\n",
        $name,
        $pooled['or'][$i + 1],
        $pooled['ci_low'][$i + 1],
        $pooled['ci_high'][$i + 1],
        $p_raw[$name],
        $p_adj[$name]
    );
}
```

This is essentially the analytic pipeline of the underlying master
thesis at UHBC (Chlef 2025–2026), reproduced here as a worked example.
