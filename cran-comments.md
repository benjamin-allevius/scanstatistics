## Test environments
* local ubuntu 16.04 install, R 3.4.1
* win-builder (devel and release)

## R CMD check results

### Using Ubuntu 16.04
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* installed size is 6.8Mb, sub-directories of 1Mb or more: libs 6.3Mb

Explanation: The size is due to use of templated classes and functions, and 
virtual functions. The cost, in terms of code duplication and inability to add 
new functionality, that would result from not using these features of C++ is
simply too high. Thus, I think the larger size of the installed package is 
warranted.

### Using win-builder
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* New maintainer: Benjamin Allévius <benjak@math.su.se>,
  Old maintainer(s): Benjamin Kjellson <benjak@math.su.se>

Explanation: Since the last release of the package, I have changed my last name 
from Kjellson to Allévius.

# Downstream dependencies
There are currently no downstream dependencies for this package.

