## Summary

<!-- One-paragraph description of what this PR does and why. -->

## Type of change

- [ ] Bug fix (numerical correctness)
- [ ] New statistical method
- [ ] Improvement to an existing method (performance, numerical stability)
- [ ] Documentation update
- [ ] Refactor / cleanup
- [ ] Other

## Checklist

- [ ] I ran `composer test` and every assertion passes.
- [ ] I ran `composer analyse` (PHPStan level 6) and no new issues.
- [ ] I ran `composer cs` (PSR-12) and the code is compliant.
- [ ] I updated `CHANGELOG.md` under the relevant version section.

### If this PR adds or modifies a statistical method

- [ ] I added the reference value to `docs/validation-tables.md`,
      with the R / SPSS command used to compute it.
- [ ] I added a corresponding `assertNear()` call in `tests/`.
- [ ] I added the mathematical specification to
      `docs/statistical-methods.md`.
- [ ] I cited the canonical peer-reviewed reference in the PHPDoc
      and in `paper.bib`.

## Related issues

<!-- Closes #..., refs #... -->
