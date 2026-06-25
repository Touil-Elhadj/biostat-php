<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for the Newton–Raphson logistic regression and the
 * Hosmer–Lemeshow goodness-of-fit test.
 *
 * Reference values were obtained in R 4.3.0 (and reproduced by manual
 * IRLS / Newton-Raphson in Python — see tests/fixtures/reference-values.R):
 *
 *   y <- c(0,0,0,0,0, 0,0,1,0,1, 0,1,0,1,1, 1,1,0,1,1, 1,0,1,1,1, 1,1,1,1,1)
 *   x <- c(1,2,2,3,3, 4,4,4,5,5, 5,6,6,6,7, 7,7,8,8,8, 9,9,9,10,10,
 *          10,11,11,12,12)
 *   m <- glm(y ~ x, family = binomial(link = "logit"))
 *   summary(m)
 *   # (Intercept) -3.9015   beta_x 0.6814
 *   # exp(0.6814) = OR_x = 1.9767
 *   # p(x) < 0.001
 */
final class LogisticRegressionTest extends BiostatTestCase
{
    /**
     * @var array<int, int>
     */
    private array $y = [
        0, 0, 0, 0, 0,  0, 0, 1, 0, 1,
        0, 1, 0, 1, 1,  1, 1, 0, 1, 1,
        1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
    ];

    /**
     * @var array<int, int|float>
     */
    private array $x = [
        1, 2, 2, 3, 3,   4, 4, 4, 5, 5,
        5, 6, 6, 6, 7,   7, 7, 8, 8, 8,
        9, 9, 9, 10, 10, 10, 11, 11, 12, 12,
    ];

    public function testSimpleLogisticReturnsArray(): void
    {
        $r = $this->stats->logisticRegression($this->y, $this->x);
        $this->assertIsArray($r);
        $this->assertArrayHasKey('coef', $r);
        $this->assertArrayHasKey('or', $r);
        $this->assertArrayHasKey('p', $r);
    }

    public function testSimpleLogisticCoefficient(): void
    {
        $r = $this->stats->logisticRegression($this->y, $this->x);
        $this->assertNear(0.6814, (float) $r['coef'], 0.01, 'beta_x');
    }

    public function testSimpleLogisticOddsRatio(): void
    {
        $r = $this->stats->logisticRegression($this->y, $this->x);
        $this->assertNear(1.9767, (float) $r['or'], 0.02, 'OR_x = exp(beta_x)');
    }

    public function testSimpleLogisticIsSignificant(): void
    {
        $r = $this->stats->logisticRegression($this->y, $this->x);
        $this->assertLessThan(0.05, (float) $r['p'], 'p-value should be < 0.05');
    }

    public function testHosmerLemeshow(): void
    {
        // Build calibrated predictions from the fitted simple logistic
        $simple = $this->stats->logisticRegression($this->y, $this->x);
        $b0     = $simple['intercept'] ?? 0.0;
        $b1     = $simple['coef'];
        $p      = array_map(
            fn($xi) => 1 / (1 + exp(-($b0 + $b1 * $xi))),
            $this->x
        );
        $hl = $this->stats->hosmerLemeshow($this->y, $p, 5);
        $this->assertIsArray($hl);
        $this->assertArrayHasKey('chi2', $hl);
        $this->assertArrayHasKey('p', $hl);
        // Well-fitted model → HL p should be > 0.05 (no evidence of misfit)
        $this->assertGreaterThan(0.05, (float) $hl['p'], 'HL p-value');
    }

    public function testMultivariateLogisticRunsOnTwoCovariates(): void
    {
        // We deliberately give a 2nd covariate that is correlated with the
        // first but does add some information — strong separation otherwise
        // makes the multivariate logistic diverge on this small sample.
        $X = [];
        foreach ($this->x as $i => $xi) {
            $X[] = [$xi, $i % 2];   // x1 = predictor, x2 = parity flag
        }
        $r = $this->stats->logisticRegressionMulti($this->y, $X, ['x1', 'x2']);
        $this->assertIsArray($r);
        $this->assertArrayHasKey('coef', $r);
        $this->assertArrayHasKey('auc', $r);
        $this->assertGreaterThan(0.5, (float) $r['auc'], 'AUC should be > 0.5');
    }
}
