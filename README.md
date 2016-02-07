### General Remarks

Grand-total number of CRAN downloads since initial release to CRAN (2015-07-28), 
as logged by [RStudio CRAN mirror](http://cran-logs.rstudio.com/):

[![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

Number of CRAN downloads in the last month:

[![](http://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

Number of CRAN downloads in the last week:

[![](http://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

=======
Travis CI build result:

[![Build Status](https://travis-ci.org/jedazard/PRIMsrc.png?branch=master)](https://travis-ci.org/jedazard/PRIMsrc)

CRAN checks:

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PRIMsrc)](https://cran.r-project.org/web/checks/check_results_PRIMsrc.html).

===============
### Description

PRIMsrc performs a unified treatment of Bump 
Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and 
Classification settings (SRC). The method generates decision rules 
delineating a region in the predictor space, where the response is larger 
than its average over the entire space. The region is shaped as a 
hyperdimensional box or hyperrectangle that is not necessarily contiguous. 


Assumptions are that the multivariate input covariates can be discrete or 
continuous and the univariate response variable can be discrete 
(Classification), continuous (Regression) or a time-to event, possibly 
censored (Survival). It is intended to handle low and high-dimensional 
multivariate datasets, including the situation where the number of covariates 
exceeds or dominates that of samples (p > n or p >> n paradigm). 

============
### Branches

The default branch (master) hosts the current development release (version 
0.6.3) of the survival bump hunting procedure that implements the case of a 
survival response. At this point, this version is also restricted to a 
directed peeling search of the first box covered by the recursive coverage 
(outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm. New 
features will be added soon as they are available. 

The main function relies on an optional variable pre-selection procedure that is run before the PRSP algorithm. 
At this point, this is done by a cross-validated penalization of the partial likelihood using the R package [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html).

In this version, the bump hunting procedure and the cross-validation procedures that control the model size and model peeling length are carried out by two separate procedures within a single main function `sbh()` that generates a unique S3-class object 'PRSP'.  


- The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. 
Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. 
Specifically, it does a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion of choice: 
Concordance Error Rate (CER), Log-Rank Test (LRT), or Log-Hazard Ratio (LHR) (see companion paper below for details).

- The second branch (unified) will host the future complete version of the code (version 1.0.0), including undirected peeling search by Patient Rule Induction Method (PRIM) that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).

===========
### License

PRIMsrc is Open Source / Free Software, available under the GNU General Public License, version 3. 
See details [here](https://github.com/jedazard/PRIMsrc/blob/master/LICENSE).

==============
### Wiki

See Wiki page [here](https://github.com/jedazard/PRIMsrc/wiki) for Publications, Roadmap, Documentation and Manual, Usage and Examples, Installation, Requirements and Support.

