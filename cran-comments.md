## Resubmission
This is a resubmission. In this version I have followed the recommandations:
* Updated description
* Removed the unnecessary license file
* changed the package name to MoMPCA (old name was MMPCA) because of conflicting name
  with a CRAN package submitted in-between.
* Changed \dontrun by \donttest for long examples
* 

__Note__ : I had to remove some cross-references due to a warning on Debian : "Non-file package-anchored link(s) in documentation object". It seems related to this issue : 
https://github.com/r-lib/roxygen2/issues/707


## Test environments
* local linux ubuntu 16.04 xenial, R 3.4.4
* ubuntu 14.04 (on travis-ci), devel and release
* OS X xcode10.2 (on travis-ci), devel and release
* Windows Server 2008 R2 SP1 (R-hub), R-devel, 32/64 bit 
* Fedora Linux (R-hub), R-devel, clang, gfortran
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Windows (win-builder), R Under development (unstable) (2019-07-05 r76784)

## R CMD check results
── R CMD check results ─────────────────────────────────────── MoMPCA 1.0.0 ────
Duration: 1m 49.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
