# Statistical methods — mathematical specification

This document gives the formal specification of every public method of
`BiostatAnalysis`. Each section contains:

1. the **model / statistic**,
2. the **algorithm** (closed-form or iterative),
3. the **peer-reviewed reference**,
4. cross-references to the unit test(s) and to the validation table
   row(s) that verify the implementation.

> Notation. $y$ is the outcome vector of length $n$; $X$ is the
> $n \times p$ design matrix; $\beta$ is the coefficient vector;
> $\sigma(z) = 1/(1 + e^{-z})$.

---

## 1. Descriptive statistics

| Method | Formula | Verification |
|---|---|---|
| `mean(x)` | $\bar x = \tfrac{1}{n} \sum x_i$ | `DescriptiveTest::testMean` |
| `std(x)` | $s = \sqrt{\tfrac{1}{n-1} \sum (x_i - \bar x)^2}$ (R-compatible *n*−1) | `DescriptiveTest::testStd` |
| `median(x)` | order statistic $x_{((n+1)/2)}$ (or mean of two middle) | `DescriptiveTest::testMedian*` |
| `quantile(x, p)` | linear interpolation between consecutive order statistics, matching R `type = 7` | `DescriptiveTest::testQuantile*` |

---

## 2. χ² test for 2 × 2 tables — `chi2Test2x2`

For observed counts $O_{ij}$ and expected counts $E_{ij}$:

$$
\chi^2 = \sum_{i,j} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}, \quad df = (r-1)(c-1) = 1.
$$

With Yates' continuity correction enabled:

$$
\chi^2_{\text{Yates}} = \sum_{i,j} \frac{(|O_{ij} - E_{ij}| - 0.5)^2}{E_{ij}}.
$$

The *p*-value is computed from the upper-tail χ² survival function
implemented in the `Distributions` trait.

**Verification**: `CategoricalTest` rows 1–2; validation table 2.

---

## 3. Odds ratio with 95 % CI — `oddsRatio`

For the 2 × 2 table $(a, b, c, d)$:

$$
\widehat{OR} = \frac{a d}{b c}, \qquad
\widehat{SE}(\ln \widehat{OR}) = \sqrt{\tfrac{1}{a}+\tfrac{1}{b}+\tfrac{1}{c}+\tfrac{1}{d}},
$$

$$
CI_{95\%} = \exp\bigl(\ln \widehat{OR} \pm 1.96 \cdot \widehat{SE}\bigr).
$$

When any cell is zero, the Haldane–Anscombe continuity correction
(`+0.5` added to every cell) is applied automatically.

**Verification**: `CategoricalTest::testOddsRatio*`; validation table 2.

---

## 4. Welch's *t*-test — `tTest`

For two independent samples with possibly unequal variances:

$$
t = \frac{\bar x_1 - \bar x_2}{\sqrt{s_1^2/n_1 + s_2^2/n_2}}, \quad
\nu = \frac{(s_1^2/n_1 + s_2^2/n_2)^2}{\frac{(s_1^2/n_1)^2}{n_1-1} + \frac{(s_2^2/n_2)^2}{n_2-1}}.
$$

The two-sided *p*-value uses the Student-*t* survival with $\nu$ df.
**Note**: $\nu$ is in general non-integer (Welch–Satterthwaite); the
implementation accepts fractional df.

**Verification**: `MeansTest::testWelch*`; validation table 3.

---

## 5. One-way ANOVA — `anova`

$$
F = \frac{MS_{\text{between}}}{MS_{\text{within}}} =
\frac{\sum_g n_g (\bar y_g - \bar y)^2 / (k-1)}
     {\sum_g (n_g - 1) s_g^2 / (n - k)}, \quad df_1 = k-1,\ df_2 = n-k.
$$

**Verification**: `MeansTest::testAnova*`; validation table 4.

---

## 6. Correlation

`pearson(x, y)`:

$$
r = \frac{\sum (x_i - \bar x)(y_i - \bar y)}
        {\sqrt{\sum (x_i - \bar x)^2 \sum (y_i - \bar y)^2}}, \quad
t = r\sqrt{\frac{n-2}{1 - r^2}}.
$$

`spearman(x, y)`: Pearson correlation of the **mid-ranks** of $x$ and
$y$, with the same *t*-based significance test.

**Verification**: `CorrelationTest`; validation table 5.

---

## 7. Logistic regression — `logisticRegression`, `logisticRegressionMulti`

Model: $\Pr(Y_i = 1 | x_i) = \sigma(x_i^\top \beta)$. Estimation by
**Newton–Raphson on the score equation**:

$$
\beta^{(k+1)} = \beta^{(k)} + (X^\top W^{(k)} X)^{-1} X^\top (y - \mu^{(k)}),
$$

with $W^{(k)} = \mathrm{diag}\bigl(\mu_i^{(k)}(1 - \mu_i^{(k)})\bigr)$.
Iteration stops when $\max_j |\beta^{(k+1)}_j - \beta^{(k)}_j| < 10^{-6}$
or after 50 iterations. The covariance matrix is $(X^\top W X)^{-1}$;
SEs are its diagonal square roots.

**AUC** is the Mann–Whitney $U$ statistic on the predicted
probabilities. **Hosmer–Lemeshow** divides the data into $g = 10$
deciles of predicted probability and refers
$\sum_j (O_j - E_j)^2 / (E_j (1 - E_j/n_j))$ to $\chi^2_{g-2}$.

**Verification**: `LogisticRegressionTest`; validation table 7.

---

## 8. Benjamini–Hochberg FDR — `benjaminiHochberg`

Given $m$ ordered *p*-values $p_{(1)} \leq \cdots \leq p_{(m)}$, the BH
adjusted *p*-values are

$$
\tilde p_{(i)} = \min_{j \geq i} \frac{m \, p_{(j)}}{j},
$$

with monotonicity enforced by the cumulative minimum.

**Reference**: @benjamini1995controlling.
**Verification**: `BenjaminiHochbergTest`; validation table 6.

---

## 9. Variance Inflation Factor — `vif`

For each predictor $X_j$, the auxiliary OLS regression of $X_j$ on the
other predictors yields $R_j^2$, and

$$
\mathrm{VIF}_j = \frac{1}{1 - R_j^2}.
$$

Thresholds reported (Allison 2012): VIF > 2.5 = noticeable, VIF > 5 =
problematic, VIF > 10 = severe multicollinearity.

**Verification**: `VifTest`; validation table 8.

---

## 10. Box–Tidwell test — `boxTidwell`

For continuous predictors $X_j$, the term $X_j \ln X_j$ is added to a
logistic regression that includes $X_j$ and the other covariates. The
Wald *p*-value on the coefficient of $X_j \ln X_j$ assesses departure
from linearity of the logit in $X_j$.

**Reference**: @box1962transformation.

---

## 11. GLMM (logistic, random intercept) — `glmmLogistic`

$$
\Pr(Y_{ij} = 1 | x_{ij}, u_i) = \sigma(x_{ij}^\top \beta + u_i), \quad
u_i \sim \mathcal N(0, \sigma_u^2).
$$

Estimation by **Penalised Quasi-Likelihood** (Breslow & Clayton, 1993)
using Henderson's mixed-model equations at each iteration:

$$
\begin{pmatrix} X^\top W X & X^\top W Z \\ Z^\top W X & Z^\top W Z + \sigma_u^{-2} I \end{pmatrix}
\begin{pmatrix} \beta \\ u \end{pmatrix} =
\begin{pmatrix} X^\top W y^* \\ Z^\top W y^* \end{pmatrix},
$$

where $y^* = \eta^{(k)} + (y - \mu^{(k)})/W^{(k)}$. The variance
component $\sigma_u^2$ is updated by ML on the residuals of $u$; the
intra-cluster correlation on the latent scale is $\sigma_u^2 / (\sigma_u^2 + \pi^2/3)$.

**Caveat**: PQL is known to underestimate $\sigma_u^2$ when cluster
sizes are small and the outcome is rare; the routine returns
`converged` and the within-cluster variance estimate so the caller
can detect this regime.

**Reference**: @breslow1993approximate, @henderson1953estimation.

---

## 12. GEE (logistic, exchangeable) — `geeLogistic`

Marginal model $\Pr(Y_{ij} = 1) = \sigma(x_{ij}^\top \beta)$ with the
exchangeable working correlation

$$
R(\alpha)_{jk} = \begin{cases} 1 & j = k \\ \alpha & j \neq k \end{cases}.
$$

The estimating equation is

$$
\sum_i D_i^\top V_i^{-1} (y_i - \mu_i(\beta)) = 0,
\quad V_i = A_i^{1/2} R(\alpha) A_i^{1/2},
$$

with $A_i = \mathrm{diag}(\mu_{ij}(1-\mu_{ij}))$. The **Liang–Zeger
sandwich variance** is

$$
\widehat{\mathrm{Var}}(\hat \beta) = M_0^{-1} M_1 M_0^{-1},
$$

with $M_0 = \sum_i D_i^\top V_i^{-1} D_i$ and
$M_1 = \sum_i D_i^\top V_i^{-1} (y_i - \mu_i)(y_i - \mu_i)^\top V_i^{-1} D_i$.

Both the model-based ($M_0^{-1}$) and the robust standard errors are
reported.

**Reference**: @liang1986longitudinal.

---

## 13. Multiple imputation by chained equations — `mice`

For each variable $X_j$ with missing values:

1. Fit a model $X_j \sim X_{-j}$ on the observed rows:
   - **continuous**: linear regression then Predictive Mean Matching
     with $d = 5$ donors,
   - **ordered factor**: proportional-odds logistic regression,
   - **unordered factor**: multinomial logistic regression.
2. Draw an imputation from the predictive distribution.
3. Iterate over all variables → one Gibbs sweep. Default 20 sweeps to
   reach approximate stationarity.
4. Repeat the whole procedure $m = 20$ times to produce $m$ imputed
   datasets.

**Reference**: @vanbuuren2011mice.

---

## 14. Rubin's pooling — `rubinPool`

Given $m$ point estimates $\hat\beta^{(k)}$ and within-imputation
variances $U^{(k)}$:

$$
\bar\beta = \tfrac{1}{m} \sum_k \hat\beta^{(k)}, \quad
\bar U = \tfrac{1}{m} \sum_k U^{(k)}, \quad
B = \tfrac{1}{m-1} \sum_k (\hat\beta^{(k)} - \bar\beta)^2,
$$

$$
T = \bar U + (1 + 1/m) B, \quad
\nu = (m-1) \left(1 + \frac{\bar U}{(1+1/m) B}\right)^2.
$$

The Barnard–Rubin adjusted df is reported. The fraction of missing
information is

$$
\mathrm{FMI} = \frac{(1 + 1/m) B + 2/(\nu+3)}{T}.
$$

**Reference**: @rubin1987multiple.

---

## Linear algebra helpers (protected)

`matMul`, `matVec`, `matTranspose`, `matrixInverse`, `olsRegression`.
All operate on arrays of arrays. `matrixInverse` uses Gauss–Jordan
elimination with partial pivoting; `olsRegression` solves the normal
equations $(X^\top X) \beta = X^\top y$ directly. Tolerance threshold:
pivot magnitude < $10^{-12}$ is treated as singular.

These are protected methods of the `LinearAlgebra` trait used
internally; they are not part of the public API but are documented
for the benefit of downstream maintainers.
