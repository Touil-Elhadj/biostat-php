<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp;

/**
 * Probability-distribution CDF helpers used by the inferential methods.
 *
 * All routines return upper-tail probabilities (i.e. p-values) except the
 * standard-normal CDF, which is two-sided as conventionally defined.
 *
 * Numerical sources:
 *   • Standard normal — Zelen & Severo (1964) approximation, |error| < 7.5e-8.
 *   • Chi-squared — Wilson–Hilferty cube-root transformation for df > 1.
 *   • Student-t and Fisher-F — regularised incomplete beta function
 *     (Press et al., "Numerical Recipes in C", §6.4).
 */
trait Distributions
{
    /**
     * Two-sided cumulative distribution of the standard normal: Φ(z).
     *
     * Uses the Zelen–Severo polynomial approximation (1964).
     *
     * @param float $z standardised value
     * @return float P(Z ≤ z), in [0, 1]
     */
    protected function normalCDF(float $z): float
    {
        $t = 1.0 / (1.0 + 0.2316419 * abs($z));
        $d = 0.3989423 * exp(-$z * $z / 2.0);
        $p = $d * $t * (
            0.3193815 + $t * (
                -0.3565638 + $t * (
                    1.781478 + $t * (
                        -1.821256 + $t * 1.330274
                    )
                )
            )
        );

        return $z > 0 ? 1.0 - $p : $p;
    }

    /**
     * Upper-tail χ² survival function: P(χ²_df > x).
     *
     * Special-cases df = 1 (where the survival reduces to the normal tail).
     * For df > 1, uses the Wilson–Hilferty cube-root transformation that
     * maps χ² to an approximately standard-normal variable.
     *
     * @param float     $x   chi-squared statistic (must be ≥ 0)
     * @param int|float $df  degrees of freedom (must be ≥ 1)
     * @return float P(X > x), in [0, 1]
     */
    protected function chi2CDF(float $x, $df): float
    {
        if ($x <= 0) {
            return 1.0;
        }

        if ((int)$df === 1) {
            // Survival of χ²₁ equals 2·(1 − Φ(√x))
            return 2.0 * (1.0 - $this->normalCDF(sqrt($x)));
        }

        // Wilson–Hilferty: ((x/df)^(1/3) − (1 − 2/(9df))) / sqrt(2/(9df)) ~ N(0,1)
        $z = pow($x / $df, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * $df));
        $z /= sqrt(2.0 / (9.0 * $df));

        return 1.0 - $this->normalCDF($z);
    }

    /**
     * Two-sided Student-t survival function: P(|T_df| > |t|).
     *
     * Uses the identity
     *   P(T > t) = ½ · I_{df/(df+t²)}(df/2, 1/2)
     * where I_x is the regularised incomplete beta function.
     *
     * @param float     $t   t-statistic
     * @param int|float $df  degrees of freedom (must be ≥ 1; fractional df
     *                       arising from Welch–Satterthwaite are accepted)
     * @return float two-sided p-value, in [0, 1]
     */
    protected function studentTCDF(float $t, $df): float
    {
        if ($df < 1) {
            return 0.5;
        }
        $t = abs($t);
        $x = $df / ($df + $t * $t);
        return $this->betaIncReg($x, $df / 2.0, 0.5);
    }

    /**
     * Fisher-F upper-tail survival function: P(F_{df1,df2} > f).
     *
     * Uses the identity
     *   P(F > f) = I_{df2/(df2+df1·f)}(df2/2, df1/2).
     *
     * @param float     $F    F-statistic
     * @param int|float $df1  numerator degrees of freedom
     * @param int|float $df2  denominator degrees of freedom
     * @return float p-value, in [0, 1]
     */
    protected function fCDF(float $F, $df1, $df2): float
    {
        if ($F <= 0) {
            return 1.0;
        }
        $x = $df2 / ($df2 + $df1 * $F);
        return $this->betaIncReg($x, $df2 / 2.0, $df1 / 2.0);
    }

    /**
     * Regularised incomplete beta function I_x(a, b).
     *
     * Continued-fraction expansion of Press et al. (Numerical Recipes 3e,
     * §6.4). Absolute accuracy ≈ 3 × 10⁻⁷ across the full (0,1) range.
     *
     * @param float $x  point of evaluation, in [0, 1]
     * @param float $a  first shape parameter (must be > 0)
     * @param float $b  second shape parameter (must be > 0)
     * @return float I_x(a, b), in [0, 1]
     */
    protected function betaIncReg(float $x, float $a, float $b): float
    {
        if ($x <= 0.0) {
            return 0.0;
        }
        if ($x >= 1.0) {
            return 1.0;
        }

        $bt = exp(
            $this->lnGamma($a + $b) - $this->lnGamma($a) - $this->lnGamma($b)
            + $a * log($x) + $b * log(1.0 - $x)
        );

        if ($x < ($a + 1.0) / ($a + $b + 2.0)) {
            return $bt * $this->betacf($a, $b, $x) / $a;
        }
        return 1.0 - $bt * $this->betacf($b, $a, 1.0 - $x) / $b;
    }

    /**
     * Lentz's algorithm for the continued-fraction component of the
     * incomplete beta function (Numerical Recipes, §6.4).
     *
     * @param float $a
     * @param float $b
     * @param float $x
     * @return float continued fraction at (a, b, x)
     */
    private function betacf(float $a, float $b, float $x): float
    {
        $maxIt = 200;
        $eps   = 3e-7;
        $fpmin = 1e-30;

        $qab = $a + $b;
        $qap = $a + 1.0;
        $qam = $a - 1.0;
        $c   = 1.0;
        $d   = 1.0 - $qab * $x / $qap;
        if (abs($d) < $fpmin) {
            $d = $fpmin;
        }
        $d = 1.0 / $d;
        $h = $d;

        for ($m = 1; $m <= $maxIt; $m++) {
            $m2 = 2 * $m;
            $aa = $m * ($b - $m) * $x / (($qam + $m2) * ($a + $m2));
            $d  = 1.0 + $aa * $d;
            if (abs($d) < $fpmin) {
                $d = $fpmin;
            }
            $c = 1.0 + $aa / $c;
            if (abs($c) < $fpmin) {
                $c = $fpmin;
            }
            $d  = 1.0 / $d;
            $h *= $d * $c;
            $aa = -($a + $m) * ($qab + $m) * $x / (($a + $m2) * ($qap + $m2));
            $d  = 1.0 + $aa * $d;
            if (abs($d) < $fpmin) {
                $d = $fpmin;
            }
            $c = 1.0 + $aa / $c;
            if (abs($c) < $fpmin) {
                $c = $fpmin;
            }
            $d   = 1.0 / $d;
            $del = $d * $c;
            $h  *= $del;
            if (abs($del - 1.0) < $eps) {
                break;
            }
        }
        return $h;
    }

    /**
     * Lanczos approximation of ln Γ(x), good for x > 0.
     *
     * Source: Numerical Recipes 3e, §6.1.
     *
     * @param float $x must be > 0
     * @return float natural logarithm of the gamma function
     */
    private function lnGamma(float $x): float
    {
        static $cof = [
            76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5,
        ];

        $y    = $x;
        $tmp  = $x + 5.5;
        $tmp -= ($x + 0.5) * log($tmp);
        $ser  = 1.000000000190015;
        foreach ($cof as $c) {
            $y++;
            $ser += $c / $y;
        }
        return -$tmp + log(2.5066282746310005 * $ser / $x);
    }
}
