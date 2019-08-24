# Release version 1.0.0 (2017-06-26)

## Release summary

- This is the first submission.

- The two NOTES from R CMD check are expected.

- The WARNING from R CMD check is related to the boost library in the BH package. More precisely, it is related to the file `bernoulli_details.hpp` which is not used by my package (it uses only `owens-t.hpp` which does not include `bernoulli_detail.hpp`). I do not know how I could fix this warning.

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

   * WARNING
   
Found the following significant warnings:

  d:/RCompile/CRANpkg/lib/3.7/BH/include/boost/math/special_functions/detail/bernoulli_details.hpp:72:36: warning: ISO C++ 1998 does not support 'long long' [-Wlong-long]
  
  d:/RCompile/CRANpkg/lib/3.7/BH/include/boost/math/special_functions/detail/bernoulli_details.hpp:96:9: warning: ISO C++ 1998 does not support 'long long' [-Wlong-long]
  
See 'd:/RCompile/CRANguest/R-devel/OwenQ.Rcheck/00install.out' for details.
