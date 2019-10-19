# PCMBaseCpp 0.1.6

* Fixed an crash situation on Windows due to exceptions bypassing the 
`try(..., silent=TRUE)` construct. Added a TraversalTaskWrapper class for 
handling exceptions on the C++ level and returning a status string to the R-level. 
See https://github.com/venelin/PCMBaseCpp/issues/1#issue-507813590 and 
https://github.com/RcppCore/Rcpp/issues/1001#issue-509047134 for details.
* Added a `NEWS.md` file to track changes to the package.
