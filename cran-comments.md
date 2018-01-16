## Test environments
* local ubuntu 16.04 install, R 3.4.3
* win-builder (devel and release)

## Fixes to CRAN Package Check Results

* WARNINGS from the flag -Wreorder have been fixed.
* WARNINGS about unused package (gamlss.dist) import have been fixed by removing 
  the imported package. Note: gamlss.dist is used (with namespace qualifier) in 
  a dontrun example.
* The NOTE about package size is explained below (and was present when the 
  previous version of the package was submitted).
* The ERROR on Flavor: r-oldrel-osx-x86_64 states "this R is version 3.3.2, 
  package 'scanstatistics' requires R >=  3.4". Section 1.1.3 of Writing R 
  Extensions does not mention that a minimum backward compatability in terms of
  R versions is needed; this ERROR thus seems to be an acceptable one.

## R CMD check results

### Using Ubuntu 16.04
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* installed size is 7.0Mb, sub-directories of 1Mb or more: libs 6.4Mb

Explanation: The size is due to use of templated classes and functions, and 
virtual functions. The cost, in terms of code duplication and inability to add 
new functionality, that would result from not using these features of C++ is
simply too high. Thus, I think the larger size of the installed package is 
warranted.

### Using win-builder
There were no ERRORs, WARNINGs or NOTEs. 

# Downstream dependencies
There are currently no downstream dependencies for this package.

