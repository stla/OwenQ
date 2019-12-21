# Release version 1.0.3 (2019-12-21)

## Release summary 

Fixed the "additional issues" of CRAN checks.

## Test environments

   * Windows 10 64bit, R 3.6.2
   * Ubuntu 18.04.3 LTS, R 3.6.1
   * online win-builder.r-project.org 

## R CMD check results

Status: OK except on win-devel, because there is no package testthat


____

# Release version 1.0.2 (2019-08-30)

## Release summary 

Attempted to fix the C++ problems raised by the CRAN checks on Solaris, and 
a test which failed.

## Test environments

   * windows 7 64bit, R 3.6.1
   * online win-builder.r-project.org 

## R CMD check results

Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'St�phane Laurent <laurent_step@outlook.fr>'

Days since last update: 0

____

# Release version 1.0.1 (2019-08-29)

## Release summary 

Attempted to fix the problems raised by the CRAN checks on Solaris; 
I replaced log(2) with log(2.0).

## Test environments

   * windows 7 64bit, R 3.6.1
   * online win-builder.r-project.org 

## R CMD check results

Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'St�phane Laurent <laurent_step@outlook.fr>'

Days since last update: 1

___

# Release version 1.0.0 (2019-08-24)

## Release summary

- This is the first submission.

## Test environments

   * ubuntu 14.04, R 3.4.0
   * windows 7 64bit, R 3.3.3
   * windows 7 64bit, R 3.5.3
   * windows 7 64bit, R 3.6.1
   * online win-builder.r-project.org 

## R CMD check results

   * NOTE

Maintainer: 'Stéphane Laurent <laurent_step@outlook.fr>'

New submission

## Resubmission comments

  * I made clearer that the commented lines in the examples are not commented examples.
  * In the vignettes I replaces the `par(.......)` with `oldpar <- par(......)` followed by `par(oldpar)`.
