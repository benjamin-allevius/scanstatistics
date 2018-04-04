## Release summary

This release removes a number of internal functions that caused the package to
not load in some circumstances.

## Test environments
* local Kubuntu 16.04 install, R 3.4.4
* ubuntu 14.04.5 (on travis-ci), R 3.4.4
* win-builder (devel)


## R CMD check results

### Using Ubuntu 16.04
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* installed size is 6.9Mb, sub-directories of 1Mb or more: libs 6.4Mb

Explanation (same as given on previous releases): The size is due to use of 
templated classes and functions, and virtual functions. The cost, in terms of 
code duplication and inability to add new functionality, that would result from 
not using these features of C++ is simply too high. Thus, I think the larger 
size of the installed package is warranted.

### Using Ubuntu 14.04.5
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* installed size is 6.2Mb, sub-directories of 1Mb or more: libs 5.7Mb

Explanation: see above.

### Using win-builder
There were no ERRORs, WARNINGs or NOTEs. 

# Downstream dependencies
There are currently no downstream dependencies for this package.

