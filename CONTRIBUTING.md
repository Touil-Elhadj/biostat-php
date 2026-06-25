# Contributing to biostat-php

Thank you for considering a contribution! This document describes the
development workflow and the special requirements for **statistical**
contributions.

## Reporting bugs

Open an issue using the bug-report template. Include:

- the version of PHP (`php -v`),
- the exact steps to reproduce the problem,
- the expected vs. actual behaviour,
- the inputs you used (vectors / matrices) so we can replay them.

If you have a reference value computed in R or SPSS that disagrees with
this library's output, please also include the R / SPSS command used
and the value you obtained.

## Proposing a feature

Open an issue using the feature-request template. Briefly describe:

- the use case,
- the relevant peer-reviewed reference(s),
- the equivalent implementation in R, Python or SPSS, if any.

## Pull-request workflow

1. Fork the repository.
2. Create a topic branch from `main`:
   ```bash
   git checkout -b feature/short-description
   ```
3. Make your changes (see *Coding style* and *Statistical
   contributions* below).
4. Run the test suite:
   ```bash
   composer install
   composer test
   composer analyse
   composer cs
   ```
   Everything must pass before submitting the PR.
5. Update `docs/` and `CHANGELOG.md` if your change is user-visible.
6. Push and open a PR against `main`.

## Coding style

- **PSR-12**, enforced by `composer cs`.
- 4-space indentation, LF line endings, UTF-8 without BOM.
- `declare(strict_types=1);` on every source file.
- **Type hints everywhere**: parameters, return types and properties.
- PHPDoc with `@param`, `@return` and the relevant mathematical formula.
- Early returns instead of deeply nested `if / else`.

## Statistical contributions ★ (read carefully)

The library's credibility rests on **every public method being
numerically verifiable** against an independent reference implementation.
Therefore every new statistical method must include:

1. **A formal specification.**
   Add an entry to `docs/statistical-methods.md` with:
   - the mathematical formulation in LaTeX,
   - the algorithm (PQL iteration, Newton–Raphson step, etc.),
   - the canonical bibliographic reference (peer-reviewed; not a blog
     post).
2. **A reference value.**
   Compute the expected result in R (preferred) or SPSS for a small,
   fully-specified input, and document it in `docs/validation-tables.md`
   in the format:
   | Test | Input | PHP | R | \|Δ\| | Tolerance | Status |
3. **A PHPUnit assertion.**
   Add a test in `tests/` that calls `assertNear($expected, $actual,
   $tolerance, '<short label>')`. Tolerances are documented in the
   `README.md` table; do not silently relax them — discuss in the PR.
4. **A docblock in the source.**
   The method's PHPDoc must include the reference formula and at least
   one citation.

If your contribution adds a new caveat (e.g. "this method is biased for
rare outcomes"), please document it in the docblock **and** in
`paper.md` so that the JOSS paper stays accurate.

## What is *not* in scope

- **Domain-specific scores** (BMI Z-scores, KIDMED, IOTF cut-offs, etc.):
  these belong in downstream applications, not in the core library.
- **Plotting or visualisation**: PHP is not a plotting language; the
  library returns plain arrays that the caller can feed to Chart.js,
  ApexCharts, or whatever else they prefer.
- **I/O routines**: the library never reads from disk or a database;
  it expects PHP arrays.

## Code of conduct

By participating in this project you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Security

To report a vulnerability, please follow the procedure in
[SECURITY.md](SECURITY.md). Do not open a public issue.

## Questions

For anything that is not a bug or a feature request, please open a
[Discussion](https://github.com/Touil-Elhadj/biostat-php/discussions)
rather than an issue.
