# Changelog

## v1.1.0 - September 2026

- Pin R package versions in `installrequiredRpackages.R` for reproducible
  container builds (Bioconductor 3.21, dated CRAN snapshot 2025-08-07)
- Write `sessionInfo.txt` into the container at build time
- Move run settings from the helper scripts into `config/params.conf`
- Add an example dataset (public LL-100 cell lines)
- Expand the documentation (output table, clinical table)
- Bundle the BLAT v35 binaries in `vendor/blat/` and fix the container
  build (CA certificates, moved download URLs)

## v1.0.0 - August 2025

- Initial release
