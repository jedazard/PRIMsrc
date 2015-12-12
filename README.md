#### General Remarks

Grand-total number of downloads since initial release to CRAN (2015-07-28), as logged by [RStudio CRAN mirror](http://cran-logs.rstudio.com/):

![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc) 

Number of downloads last month:

![](http://cranlogs.r-pkg.org/badges/PRIMsrc)

Travis CI build result: 

![](https://travis-ci.org/jedazard/PRIMsrc.svg)

============
Description:
============
Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. Assumptions are that the multivariate input covariates can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event, possibly censored (Survival). It is intended to handle low and high-dimensional multivariate datasets, including the situation where the number of covariates exceeds or dominates that of samples (p > n or p >> n paradigm).

=========
Branches:
=========
- The default branch (master) hosts the current development release (version 0.6.3) of the survival bump hunting procedure that implements the case of a survival response. At this point, this version is also restricted to a directed peeling search of the first box covered by the recursive coverage (outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm. New features will be added soon as they are available.

	The main function relies on an optional variable pre-selection procedure that is run before the PRSP algorithm. At this point, this is done by a cross-validated penalization of the partial likelihood using the R package 'glmnet'.

	In this version, the bump hunting procedure and the cross-validation procedures that control the model size and model peeling length are carried out by two separate procedures within a single main function 'sbh()' that generates a unique S3-class object 'PRSP'.  


- The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. Specifically, it does a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion of choice: CER, LRT, or LHR (see companion paper below for details).

	In addition, this cross-validation procedure is carried out separately from the main function 'sbh()' in a cross-validation function called 'cv.sbh()'. Altogether, this allows a more rigorous treatment of model validation, a better control on the user-end and an improvement of the maintenance on the back-end. In the process, two S3-class objects are created instead of one: an additional S3-class object 'CV' is output by the cross-validation function cv.sbh() and used as input in the main function 'sbh()'. 


- The second branch (unified) will host the future complete version of the code (version 1.0.0), including undirected peeling search by Patient Rule Induction Method (PRIM) that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).

========
License:
========
PRIMsrc is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PRIMsrc.pdf) details the end-user functions. At this stage and for simplicity, there is a unique end-user main function for fitting a cross-validated Survival Bump Hunting model (sbh(...)). There are 4 end-user specific plotting functions (`plot_****(...)`) along with 4 S3-generic functions: `summary(...)`, `predict(...)`, `plot(...)` and `print(...)`. Available are also 5 synthetic datasets and 2 real datasets including altogether low and high-dimensional situations (for p < n, p > n and p >> n cases). See the "PRIMsrc-package" introduction section of the manual for more details and examples.

===========
References:
===========
CRAN release:
https://cran.r-project.org/web/packages/PRIMsrc/index.html


The companion papers (accepted and submitted for publication) can be accessed here:

- ASA-IMS JSM Proceedings (2014): 

https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf

- Archives arXiv:

http://arxiv.org/abs/1501.03856.

- Statistical Analysis and Data Mining. The ASA Data Science Journal (to appear):

http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1932-1872

=============
Requirements:
=============
PRIMsrc 0.6.3 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-11-04 r69597) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. See for instance the 'CRAN Package Check Results' here:

https://cran.r-project.org/web/checks/check_results_PRIMsrc.html

=============
Installation: 
=============
- To install PRIMsrc from CRAN, simply download and install the current version (0.6.3) from the CRAN repository:

install.packages("PRIMsrc")

- Alternatively, you can install the most up-to-date version (0.6.3) from GitHub, using devtools:

install.packages("devtools")

library("devtools")

devtools::install_github("jedazard/PRIMsrc")

======
Usage: 
======
- To load the PRIMsrc library in an R session and start using it:

library("PRIMsrc")

- Check the package news with the R command:

PRIMsrc.news()

- Check on how to cite the package with the R command:

citation("PRIMsrc")

etc...
