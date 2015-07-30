============
Description:
============
Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. Assumptions are that the multivariate input covariates can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event, possibly censored (Survival). It is intended to handle low and high-dimensional multivariate datasets, including the situation where the number of covariates exceeds or dominates that of samples (p > n or p >> n paradigm).

=========
Branches:
=========
- The default branch (master) hosts the current development release of the survival bump hunting procedure that implements the case of a survival response. At this point, this version is also restricted to a directed peeling search of the first box covered by the recursive coverage (outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm. New features will be added soon as they are available.

	The main function relies on an internal variable pre-selection procedure before the PRSP algorithm is run. At this point, this is done either by regular Cox-regression (from the R package 'survival') or cross-validated Elasticnet Regularized Cox-Regression (from the R package 'glmnet'), depending on whether the number of covariates is less (p <= n) or greater (p > n) than the number of samples, respectively.
	
	In this version, the bump hunting procedure and the cross-validation procedures that control the model size and model peeling length are carried out by two separate procedures within a single main function 'sbh()' that generates an S3-class object 'PRSP'.  


- The first branch (devel) hosts an alternative development version that offers a more rigorous and modular version of the code. Here, an alternative single internal cross-validation procedure is carried out in a single cross-validation function called 'cv.sbh()' to simultaneously control the model size and model peeling length before the PRSP algorithm is run. Specifically, it includes a univariate bump hunting variable selection procedure, where model size and model peeling length are simultaneously optimized by cross-validation to moptimize the cross-validation criterion of choice: CER, LRT, or LHR (see companion paper below for details).

	In addition, this cross-validation procedure is carried out separately of the main function 'sbh()'. Altogether, this allows a more rigorous treatment of model validation, a better control on the user-end and an improvement of the maintenance on the back-end. In the process, two S3-class objects are created instead of one: an additional S3-class object 'CV' is output by the cross-validation function cv.sbh() and used as input in the main function 'sbh()'. 


- The second branch (unified) will host the undirected peeling search version by Patient Rule Induction Method (PRIM) that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).

========
License:
========
PRIMsrc is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PRIMsrc.pdf) details the end-user functions. At this stage and for simplicity, there is a unique end-user main function for fitting a cross-validated Survival Bump Hunting model (sbh(...)). There are 5 end-user plotting functions (plot_****(...)) along with two S3 generic functions: summary(...) and predict(...). Available are also 5 synthetic datasets and 1 real dataset including altogether low and high-dimensional situations (for p < n, p > n and p >> n cases). See the "PRIMsrc-package" introduction section of the manual for more details and examples.

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
PRIMsrc 0.5.7 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-07-20 r68705) and Travis CI. 

Installation has been tested on Windows, Linux and OSX platforms. See for instance the 'CRAN Package Check Results' here:

https://cran.r-project.org/web/checks/check_results_PRIMsrc.html

=============
Installation: 
=============
- To install PRIMsrc from CRAN, simply download and install the current version (0.5.7) from the CRAN repository:

install.packages("PRIMsrc")

- Alternatively, you can install the most up-to-date development version (0.5.7) from GitHub, using devtools:

library(devtools)
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
