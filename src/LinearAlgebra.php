<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp;

use InvalidArgumentException;

/**
 * Linear-algebra primitives used by the higher-level statistical methods.
 *
 * All matrices are stored as PHP arrays-of-arrays in row-major order:
 *   $A[i][j]  →  i-th row, j-th column (0-indexed).
 *
 * The trait is consumed by {@see BiostatAnalysis}; client code never
 * accesses these methods directly.
 */
trait LinearAlgebra
{
    /**
     * Matrix multiplication C = A · B.
     *
     * @param array<int, array<int, float>> $A n × m
     * @param array<int, array<int, float>> $B m × p
     * @return array<int, array<int, float>>      n × p
     *
     * @throws InvalidArgumentException when inner dimensions disagree.
     */
    protected function matMul(array $A, array $B): array
    {
        $nA = count($A);
        $mA = count($A[0]);
        $nB = count($B);
        $mB = count($B[0]);

        if ($mA !== $nB) {
            throw new InvalidArgumentException(
                "matMul: incompatible dimensions ({$nA}×{$mA} · {$nB}×{$mB})"
            );
        }

        $C = [];
        for ($i = 0; $i < $nA; $i++) {
            $row = array_fill(0, $mB, 0.0);
            for ($k = 0; $k < $mA; $k++) {
                $aik = $A[$i][$k];
                if ($aik == 0.0) {
                    continue;
                }
                for ($j = 0; $j < $mB; $j++) {
                    $row[$j] += $aik * $B[$k][$j];
                }
            }
            $C[$i] = $row;
        }
        return $C;
    }

    /**
     * Matrix-vector product: A · v.
     *
     * @param array<int, array<int, float>> $A n × m
     * @param array<int, float>             $v m
     * @return array<int, float>                n
     *
     * @throws InvalidArgumentException when dimensions disagree.
     */
    protected function matVec(array $A, array $v): array
    {
        $n = count($A);
        $m = count($v);

        if (count($A[0]) !== $m) {
            throw new InvalidArgumentException('matVec: incompatible dimensions');
        }

        $out = array_fill(0, $n, 0.0);
        for ($i = 0; $i < $n; $i++) {
            $s = 0.0;
            for ($j = 0; $j < $m; $j++) {
                $s += $A[$i][$j] * $v[$j];
            }
            $out[$i] = $s;
        }
        return $out;
    }

    /**
     * Transpose: returns A^T.
     *
     * @param array<int, array<int, float>> $A n × m
     * @return array<int, array<int, float>>      m × n
     */
    protected function matTranspose(array $A): array
    {
        $n = count($A);
        $m = count($A[0]);

        $T = [];
        for ($j = 0; $j < $m; $j++) {
            $row = array_fill(0, $n, 0.0);
            for ($i = 0; $i < $n; $i++) {
                $row[$i] = $A[$i][$j];
            }
            $T[$j] = $row;
        }
        return $T;
    }

    /**
     * Matrix inverse via Gauss–Jordan elimination with partial pivoting.
     *
     * Returns null if the matrix is singular (largest pivot below 1e-12).
     *
     * @param array<int, array<int, float>> $A square matrix
     * @return array<int, array<int, float>>|null
     */
    protected function matrixInverse(array $A): ?array
    {
        $n = count($A);
        $M = [];

        // Build the augmented matrix [A | I]
        for ($i = 0; $i < $n; $i++) {
            $M[$i] = array_merge($A[$i], array_fill(0, $n, 0.0));
            $M[$i][$n + $i] = 1.0;
        }

        for ($i = 0; $i < $n; $i++) {
            // Partial pivot
            $maxRow = $i;
            for ($k = $i + 1; $k < $n; $k++) {
                if (abs($M[$k][$i]) > abs($M[$maxRow][$i])) {
                    $maxRow = $k;
                }
            }
            if (abs($M[$maxRow][$i]) < 1e-12) {
                return null;
            }
            if ($maxRow !== $i) {
                [$M[$i], $M[$maxRow]] = [$M[$maxRow], $M[$i]];
            }

            $piv = $M[$i][$i];
            for ($j = 0; $j < 2 * $n; $j++) {
                $M[$i][$j] /= $piv;
            }
            for ($k = 0; $k < $n; $k++) {
                if ($k === $i) {
                    continue;
                }
                $f = $M[$k][$i];
                for ($j = 0; $j < 2 * $n; $j++) {
                    $M[$k][$j] -= $f * $M[$i][$j];
                }
            }
        }

        $inv = [];
        for ($i = 0; $i < $n; $i++) {
            $inv[$i] = array_slice($M[$i], $n);
        }
        return $inv;
    }

    /**
     * Ordinary least-squares regression β = (X'X)⁻¹ X'y, with full diagnostics.
     *
     * The intercept must be provided as the first column of $X if desired —
     * the routine does not add one automatically.
     *
     * @param array<int, float>             $y outcome, length n
     * @param array<int, array<int, float>> $X design matrix, n × k
     *
     * @return array{
     *     beta: array<int, float>,
     *     se: array<int, float|null>,
     *     r2: float,
     *     ss_tot: float,
     *     ss_res: float,
     *     sigma2: float,
     *     pred: array<int, float>,
     *     invXtX: array<int, array<int, float>>,
     *     n: int,
     *     k: int
     * }|null   null if singular or under-determined
     */
    protected function olsRegression(array $y, array $X): ?array
    {
        $n = count($y);
        $k = count($X[0] ?? []);
        if ($n < $k + 1) {
            return null;
        }

        $Xt    = $this->matTranspose($X);
        $XtX   = $this->matMul($Xt, $X);
        $invXtX = $this->matrixInverse($XtX);
        if ($invXtX === null) {
            return null;
        }

        // β = (X'X)⁻¹ X'y
        $Xty = array_fill(0, $k, 0.0);
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $k; $j++) {
                $Xty[$j] += $X[$i][$j] * $y[$i];
            }
        }

        $beta = array_fill(0, $k, 0.0);
        for ($j = 0; $j < $k; $j++) {
            for ($l = 0; $l < $k; $l++) {
                $beta[$j] += $invXtX[$j][$l] * $Xty[$l];
            }
        }

        // Predictions, residuals, R²
        $ybar = array_sum($y) / $n;
        $ss_tot = 0.0;
        $ss_res = 0.0;
        $pred = array_fill(0, $n, 0.0);
        for ($i = 0; $i < $n; $i++) {
            $yi_hat = 0.0;
            for ($j = 0; $j < $k; $j++) {
                $yi_hat += $X[$i][$j] * $beta[$j];
            }
            $pred[$i] = $yi_hat;
            $ss_tot += ($y[$i] - $ybar) ** 2;
            $ss_res += ($y[$i] - $yi_hat) ** 2;
        }

        $r2     = $ss_tot > 0 ? max(0.0, 1.0 - $ss_res / $ss_tot) : 0.0;
        $df_res = $n - $k;
        $sigma2 = $df_res > 0 ? $ss_res / $df_res : 0.0;

        // Standard errors from diag((X'X)⁻¹) · σ²
        $se = array_fill(0, $k, null);
        for ($j = 0; $j < $k; $j++) {
            $v = $invXtX[$j][$j] * $sigma2;
            $se[$j] = $v >= 0 ? sqrt($v) : null;
        }

        return [
            'beta'   => $beta,
            'se'     => $se,
            'r2'     => $r2,
            'ss_tot' => $ss_tot,
            'ss_res' => $ss_res,
            'sigma2' => $sigma2,
            'pred'   => $pred,
            'invXtX' => $invXtX,
            'n'      => $n,
            'k'      => $k,
        ];
    }
}
