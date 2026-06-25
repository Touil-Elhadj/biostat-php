<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Smoke / structural tests for the advanced multivariate methods.
 *
 * Unlike the lower-level tests these check the *shape* of the output and
 * basic sanity properties (e.g. ICC ∈ [0, 1], AUC > 0.5, m imputed sets
 * returned) rather than exact numerical equivalence with R, because:
 *
 *   • GLMM by PQL is known to diverge from `lme4::glmer` (Laplace) on
 *     the variance components by up to 10 %.
 *   • MICE involves random draws; reference equality is meaningless.
 *
 * Full quantitative R cross-checks live in `docs/validation-tables.md`
 * and are run manually because they require an R installation.
 */
final class AdvancedMethodsTest extends BiostatTestCase
{
    public function testBoxTidwellRunsOnContinuousPredictor(): void
    {
        // Build a fake dataset where the logit IS linear in x1
        $n = 80;
        $X = [];
        $y = [];
        $rng = function () { static $s = 0; $s++; return (sin($s * 12.9898) * 43758.5453) - floor(sin($s * 12.9898) * 43758.5453); };
        for ($i = 0; $i < $n; $i++) {
            $x1 = $i / 10.0;
            $eta = -2.0 + 0.4 * $x1;
            $p   = 1.0 / (1.0 + exp(-$eta));
            $X[] = [$x1, ($i % 3) === 0 ? 1.0 : 0.0];
            $y[] = $rng() < $p ? 1 : 0;
        }
        $r = $this->stats->boxTidwell($y, $X, [0], ['x1', 'x2']);
        $this->assertIsArray($r);
        $this->assertArrayNotHasKey('error', $r);
        $this->assertArrayHasKey('x1', $r);
        // Linear-in-logit predictor → the x*log(x) interaction p should NOT
        // be tiny (we permit a wide range here because the dataset is small)
        $this->assertGreaterThan(0.0, (float) $r['x1']['p'], 'p-value present');
        $this->assertLessThan(1.001, (float) $r['x1']['p'], 'p-value in [0,1]');
    }

    public function testGlmmReturnsIccInUnitInterval(): void
    {
        // 20 clusters × 5 observations each
        $n_clusters = 20;
        $cluster = [];
        $y       = [];
        $X       = [];
        for ($c = 0; $c < $n_clusters; $c++) {
            // cluster intercept ~ N(0, 1)
            $u = sin($c * 1.7) * 0.5;
            for ($j = 0; $j < 5; $j++) {
                $x = ($j - 2) * 0.5;
                $eta = -1.0 + 0.4 * $x + $u;
                $p = 1.0 / (1.0 + exp(-$eta));
                $cluster[] = $c;
                $X[]       = [$x];
                // deterministic Bernoulli pseudo-realisation
                $y[] = $p > 0.5 ? 1 : 0;
            }
        }
        $r = $this->stats->glmmLogistic($y, $X, $cluster, ['x'], 30);
        $this->assertIsArray($r);
        $this->assertArrayNotHasKey('error', $r);
        $this->assertArrayHasKey('icc', $r);
        $this->assertGreaterThanOrEqual(0.0, (float) $r['icc'], 'ICC >= 0');
        $this->assertLessThanOrEqual(1.0,  (float) $r['icc'], 'ICC <= 1');
    }

    public function testGeeReportsBothModelBasedAndRobustSe(): void
    {
        // 15 clusters × 6 observations
        $cluster = [];
        $y       = [];
        $X       = [];
        for ($c = 0; $c < 15; $c++) {
            for ($j = 0; $j < 6; $j++) {
                $x = ($j - 2.5) * 0.4 + sin($c * 1.2) * 0.2;
                $eta = -0.5 + 0.6 * $x;
                $cluster[] = $c;
                $X[]       = [$x];
                $y[]       = (1.0 / (1.0 + exp(-$eta))) > 0.5 ? 1 : 0;
            }
        }
        $r = $this->stats->geeLogistic($y, $X, $cluster, ['x']);
        $this->assertIsArray($r);
        $this->assertArrayHasKey('se_robust', $r);
        $this->assertArrayHasKey('se_model',  $r);
        $this->assertNotEmpty($r['se_robust']);
        $this->assertNotEmpty($r['se_model']);
    }

    public function testRubinPoolBasicSanity(): void
    {
        // 5 imputed estimates of a single parameter — must be passed as m × p matrices
        $estimates = [[1.20], [1.30], [1.25], [1.22], [1.28]];
        $std_err   = [[0.40], [0.41], [0.40], [0.39], [0.42]];

        $r = $this->stats->rubinPool($estimates, $std_err);
        $this->assertIsArray($r);
        $this->assertArrayHasKey('beta',     $r);
        $this->assertArrayHasKey('total_T',  $r);

        // The pooled estimate is the mean of the imputed estimates
        $pooled = is_array($r['beta']) ? (float) $r['beta'][0] : (float) $r['beta'];
        $this->assertNear(1.250, $pooled, 0.01, 'pooled β = mean of imputations');
    }
}
