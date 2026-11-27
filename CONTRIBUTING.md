# Contributing to ethogramR

Thank you for your interest in ethogramR. Contributions of all sizes are welcome.

## Reporting problems or asking for help

Open an issue at <https://github.com/YOUR-GITHUB-USERNAME/ethogramR/issues>. Please include:

- what you ran (a minimal R example, ideally using the bundled example data or a
  few rows of your own data),
- what you expected and what happened instead,
- the output of `sessionInfo()`.

Questions about how to use a function or interpret a result are also welcome as issues.

## Suggesting features

Open an issue describing the analysis you need, the data you have (scan, focal or
video coding) and, if possible, a reference for the method.

## Contributing code

1. Fork the repository and create a branch for your change.
2. Add or update tests in `tests/testthat/` for any change in behaviour. Where a
   function computes a published index, include a test against a value worked
   out by hand.
3. Document new arguments with roxygen2 comments and run `devtools::document()`.
4. Run `devtools::check()` and make sure it finishes with no errors or warnings.
5. Add a line to `NEWS.md` and open a pull request describing the change.

## Code of conduct

Please be respectful and constructive. Harassment or discrimination of any kind
is not tolerated.
