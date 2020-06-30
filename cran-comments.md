## Resubmission
This is a resubmission. In this version I have followed the recommandations:

* Updated description
* Removed the unnecessary license file
* changed the package name to MoMPCA (old name was MMPCA) because of conflicting name
  with a CRAN package submitted in-between.
* Changed \dontrun by \donttest for long examples


The first submission came back with the following review :

__Thanks, Please always explain all acronyms (e.g. PCA) in the description field. Please only capitalize sentence beginnings and names in the description text.__

  * All acronyms were explained in DESCRIPTION and capitalization checked.


__"Currently only the branch and bound Classification Variational Expectation  Maximisation is
implemented." Do you realize that CRAN is a platform for released software products, not a development platform? The Description field is intended to be a (one paragraph) description of what the package does and why it may be useful. Please elaborate.__

 * Sorry for the unfortunate choice of words. This sentence was removed, the descriptions field improved and the link with a published research paper made explicit.

__The LICENSE file is only needed if you have additional restrictions to the GPL-3 which you have not? In that case omit the file and its reference in the DESCRIPTION file.__

 * Thank you for this precision. I omitted the LICENSE file but the kept the field License: GPL-3 in DESCRIPTION because of a R CMD check error : 
 ❯ checking for file ‘MoMPCA/DESCRIPTION’ ... ERROR
  Champ requis mais manquant ou vide :
    ‘License’

__You write information messages to the console that cannot be easily suppressed. Instead of print()/cat() rather use message()/warning() if you really have to write text to the console.__

 * Thank you for this advice : the package now uses message() instead of print() and cat(). 

__`\dontrun{}` should be only used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing, (or replace `\dontrun{}` with `\donttest{}`).__

 * Thank you for this advice, `\dontrun{}` examples are now replaced by `\donttest{}` 


__Note__ : I had to remove some cross-references due to a warning on r-devel: "Non-file package-anchored link(s) in documentation object". It seems related to this issue : 
https://github.com/r-lib/roxygen2/issues/707


## Test environments
All the following test have succeded without warnings:

  * Local linux ubuntu 16.04 xenial, R 3.6.1
  * On Travis CI:
    * Ubuntu 14.04 (on travis-ci), devel and release
    * OS X xcode10.2 (on travis-ci), devel and release
 * Using `devtools::check_win_*()`:
    * **Oldrelease**: using R version 3.6.3 (2020-02-29) using platform: x86_64-w64-mingw32 (64-bit) (2f0aYhVLuWNV)
    * **Devel **:  using R Under development (unstable) (2020-06-29 r78751) using platform: x86_64-w64-mingw32 (64-bit)
    * **Release**: using R version 4.0.2 (2020-06-22) using platform: x86_64-w64-mingw32 (64-bit)
    
Possibly mis-spelled words in DESCRIPTION:
  Classication (8:328)
  ICL (8:353)
  Jouvin (8:505)
  VEM (8:298)
  Variational (8:258)
  al (8:516)
  et (8:512)
  
  * Using `rhub::check`
    * **Fedora**: NOTE   Note: found 16 marked UTF-8 strings
    *
    *
    
## Local R CMD check results
── R CMD check results ─────────────────────────────────────── MoMPCA 1.0.0 ────
Duration: 1m 49.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
