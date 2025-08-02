
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCMBaseCppCustom

This is a customized version of the fast C++ backend for the PCMBase
R-package. This version adds a C++ implementation for the Early Burst
(EB) phylogenetic comparative model.

# Installation

The package needs a C++ 11 or later compiler and Rcpp to be installed in
your R-environment. Once this is done, you can install the most recent
version of the package from github:

``` r
devtools::install_github("MehediHasanSM/PCMBaseCppCustom", force = TRUE)
```

# Citing PCMBase

To give credit to the PCMBase package in a publication, please cite the
following articles:

Mitov, V., & Stadler, T. (2018). Parallel likelihood calculation for
phylogenetic comparative models: The SPLITT C++ library. Methods in
Ecology and Evolution, 2041–210X.13136.
<http://doi.org/10.1111/2041-210X.13136>

Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2019). Fast
likelihood calculation for multivariate Gaussian phylogenetic models
with shifts. Theor. Popul. Biol.
<https://doi.org/10.1016/j.tpb.2019.11.005>

# Used 3rd party libraries

The PCMBaseCppCustom R-package uses the following R-packages and C++
libraries:

- For tree processing in C++: The [SPLITT
  library](https://venelin.github.io/SPLITT/) (Mitov and Stadler 2018);
- For data processing in R: data.table v1.17.8 (Dowle and Srinivasan
  2019);
- For algebraic manipulation: The [Armadillo C++ template
  library](https://arma.sourceforge.net/) (Sanderson and Curtin 2016)
  and its port to R RcppArmadillo v14.6.0.1 (Eddelbuettel et al. 2019);
- For unit-testing: testthat v3.2.3 (Wickham 2019), covr v3.6.4 (Hester
  2018);
- For documentation and web-site generation: roxygen2 v7.3.2 (Wickham,
  Danenberg, and Eugster 2018), pkgdown v2.1.3 (Wickham and Hesselberth
  2018);

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-R-data.table" class="csl-entry">

Dowle, Matt, and Arun Srinivasan. 2019. *Data.table: Extension of
‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

</div>

<div id="ref-R-RcppArmadillo" class="csl-entry">

Eddelbuettel, Dirk, Romain Francois, Doug Bates, and Binxiang Ni. 2019.
*RcppArmadillo: ’Rcpp’ Integration for the ’Armadillo’ Templated Linear
Algebra Library*. <https://CRAN.R-project.org/package=RcppArmadillo>.

</div>

<div id="ref-R-covr" class="csl-entry">

Hester, Jim. 2018. *Covr: Test Coverage for Packages*.
<https://CRAN.R-project.org/package=covr>.

</div>

<div id="ref-Mitov:2018dqa" class="csl-entry">

Mitov, Venelin, and Tanja Stadler. 2018. “<span class="nocase">Parallel
likelihood calculation for phylogenetic comparative models: The SPLITT
C++ library</span>.” *Methods in Ecology and Evolution*, December,
2041–210X.13136.

</div>

<div id="ref-Sanderson:2016cs" class="csl-entry">

Sanderson, Conrad, and Ryan Curtin. 2016.
“<span class="nocase">Armadillo: a template-based C++ library for linear
algebra</span>.” *Journal of Open Source Software* 1 (2).

</div>

<div id="ref-R-testthat" class="csl-entry">

Wickham, Hadley. 2019. *Testthat: Unit Testing for r*.
<https://CRAN.R-project.org/package=testthat>.

</div>

<div id="ref-R-roxygen2" class="csl-entry">

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2:
In-Line Documentation for r*.
<https://CRAN.R-project.org/package=roxygen2>.

</div>

<div id="ref-R-pkgdown" class="csl-entry">

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static HTML
Documentation for a Package*.
<https://CRAN.R-project.org/package=pkgdown>.

</div>

</div>
