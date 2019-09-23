Hello, 

Thanks for reviewing the package. According to your comments, the following 
updates have been done:

- Shortened package title to get below 65 characters;
- Added references to the description field in DESCRIPTION;
- Removed \\dontrun statement from examples for MiniBenchmarkRvsCpp;
- To reduce the execution time, examples for this function are performed only on
the first row of the benchmarkData table. This is the minimal possible variant.
- Added examples for the exported functions: PCMParamGetFullVector, PCMInfoCpp 
and PCMTreePreorder;
- Added returned value for the function PCMTreePreorder;
- Updated .Rd files;

Now, I expect that the check on CRAN would raise spell-checking notes as follows:
  Mitov (11:61, 13:16)
  Stadler (13:26)
  al (11:70)
  et (11:67)
  
These should be ignored, since these are simply my and my co-author's family name
and "et al." is commonly used in references.

Best regards, 
Venelin