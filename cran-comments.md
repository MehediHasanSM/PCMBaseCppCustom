Hello, 

Thanks for reviewing the package. According to your comments, the following 
updates have been done:

- Shortened package title to get below 65 characters;
- Added references to the description field in DESCRIPTION;
- Removed \\dontrun statement from examples for MiniBenchmarkRvsCpp;
- To reduce the execution time, examples for this function are performed only on
the first row of the benchmarkData table. This is the minimal possible variant 
and runs within 5 seconds on my laptop.
- Added examples for the exported functions: PCMParamGetFullVector, PCMInfoCpp 
and PCMTreePreorder;
- Added returned value for the function PCMTreePreorder;
- Updated .Rd files;