# Security policy

## Supported versions

Only the latest minor release of the library is actively maintained for
security patches.

| Version  | Supported |
|----------|-----------|
| `1.x`    | ✅        |
| `< 1.0`  | ❌        |

## Reporting a vulnerability

If you discover a security vulnerability in `biostat-php` please **do
not** open a public issue. Instead, send the details by e-mail to the
address listed in `CITATION.cff`, with the subject line
`[security] biostat-php`.

Please include:

1. A description of the vulnerability and its potential impact.
2. Step-by-step instructions to reproduce it.
3. Affected versions / commit SHAs.
4. Any suggested mitigation, if you have one.

You should receive an acknowledgement within **7 days**. We will keep
you updated on the progress of the fix and credit you in the release
notes unless you prefer to remain anonymous.

## Scope

Because this is a numerical library that performs no I/O, the
plausible attack surface is small:

In scope:

- arithmetic-related vulnerabilities, including denial-of-service via
  pathological input (e.g. an input vector that causes an iterative
  solver to loop indefinitely),
- code-execution vulnerabilities (e.g. unsafe `eval`, `unserialize` —
  there should be none, but please report any),
- memory-exhaustion issues triggered by inputs of unbounded size.

Out of scope:

- correctness bugs in statistical methods — these are open issues, not
  vulnerabilities. Please file them through the normal issue tracker.
- vulnerabilities in third-party software (PHP, Composer, PHPUnit).
- issues that require the attacker to already have full PHP code
  execution on the host.
