# Security Policy

## Supported versions

Security fixes are released for the current major version. Users should run the
latest available patch release before reporting an issue that may already have
been corrected.

| Version | Security support |
| ------- | ---------------- |
| 1.x     | Yes              |
| < 1.0   | No               |

CHIMERA supports the Python versions declared in `pyproject.toml`. A Python
interpreter that has reached upstream end of life is not a supported deployment
target even if an older CHIMERA installation still starts on it.

## Reporting a vulnerability

Do not open a public issue for a suspected vulnerability. Submit a private
report through [GitHub Security Advisories](https://github.com/Alexander-Mitrofanov/CHIMERA/security/advisories/new).
If that form is unavailable, contact a repository maintainer privately and ask
for a secure reporting channel without including exploit details in the first
message.

Include, where possible:

- the affected CHIMERA version, Python version, operating system, and install
  method;
- a minimal reproduction and the security impact;
- whether untrusted FASTA, metadata, configuration, archive, or output paths are
  required to trigger the issue;
- relevant logs with credentials, access tokens, private genome data, and human
  subject information removed; and
- any known mitigations or evidence that the issue is already being exploited.

Do not attach confidential or controlled-access sequence data. A synthetic
reproducer is strongly preferred.

## Response and disclosure

Maintainers aim to acknowledge a complete report within three business days,
provide an initial assessment within ten business days, and send status updates
at least every seven days while remediation is active. These are response goals,
not guarantees. Timing depends on severity, reproducibility, and coordination
with affected upstream projects.

Please allow maintainers a reasonable embargo period to reproduce the issue,
prepare a fix, test supported releases, and coordinate disclosure. A security
advisory and patched release will credit reporters who request attribution.

## Scope

Security reports include, but are not limited to:

- path traversal or unintended file overwrite from crafted inputs;
- command or argument injection into optional external tools;
- unsafe archive handling, deserialization, or temporary-file behavior;
- dependency or container vulnerabilities with a practical CHIMERA impact;
- provenance or integrity defects that allow benchmark inputs or outputs to be
  silently substituted; and
- accidental disclosure of local paths, credentials, or controlled data.

A scientific-validity defect, unexpected biological result, or documentation
error without a confidentiality, integrity, or availability impact should be
reported through the public issue tracker instead. If the distinction is
unclear, use the private security channel.

## Operational guidance

Treat downloaded references and user-supplied metadata as untrusted input. Run
CHIMERA with the least filesystem privileges required, keep reference manifests
and checksums with published datasets, avoid embedding secrets in configuration
files, and use the non-root container image for isolated workflows.

Last updated: 2026-08-08.
