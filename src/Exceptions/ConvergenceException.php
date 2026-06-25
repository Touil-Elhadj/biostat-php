<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Exceptions;

use RuntimeException;

/**
 * Thrown by iterative routines (logistic regression, GLMM by PQL, GEE)
 * when the maximum number of iterations is reached before the
 * convergence tolerance is met.
 *
 * Callers can still recover the last iterate from the exception's
 * properties when partial results are acceptable.
 */
class ConvergenceException extends RuntimeException
{
    /** @var array<int, float>|null Final coefficient vector at last iteration. */
    public ?array $lastEstimate;

    /** @var int Number of iterations actually executed. */
    public int $iterations;

    /** @var float Final residual change between the last two iterates. */
    public float $finalDelta;

    /**
     * @param string                   $message       human-readable description
     * @param int                      $iterations    how many iterations ran
     * @param float                    $finalDelta    last between-iter change
     * @param array<int, float>|null   $lastEstimate  the un-converged coefficients
     */
    public function __construct(
        string $message,
        int $iterations,
        float $finalDelta,
        ?array $lastEstimate = null
    ) {
        parent::__construct($message);
        $this->iterations   = $iterations;
        $this->finalDelta   = $finalDelta;
        $this->lastEstimate = $lastEstimate;
    }
}
