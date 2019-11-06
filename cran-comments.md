## Resubmission
This is a resubmission. In this version I have:
* changed all T, F  to TRUE, FALSE, respectively

## Resubmission
This is a resubmission. In this version I have:
* wrapped the (working) example of the function TC_CAR1_sc inside \donttest
  to avoid NOTE about elapsed time (> 10s)

## Resubmission
This is a resubmission. In this version I have:
* added the publication year to reference in the Description field of DESCRIPTION file
* included toy runable examples for the main function TC_CAR1
* set name_dir_sim = NULL in the example of the function sc_TC_CAR1 
  so that the function does not modify the user's home directory. 
  This example is set to \donttest since its elapsed time > 5s
* replaced all cat() by message()


## Resubmission
This is a resubmission. In this version I have:
* fixed typo CRAN URL link in the README.md file

## Resubmission
This is a resubmission. In this version I have:

* provided reference describing the method in the description field of DESCRIPTION file
* added all authors (no copyright holders except me) in the Authors@R field with the appropriate roles
* updated license to GPL (>=2)

## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
