# PRIMsrc
Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification 


=============
### License

PRIMsrc is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


=============
### Downloads

CRAN downloads since initial release to CRAN (2015-07-28):
<<<<<<< HEAD
[![](https://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
as tracked by [RStudio CRAN mirror](http://cran-logs.rstudio.com/)

CRAN downloads in the last month:
[![](https://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)

CRAN downloads in the last week:
[![](https://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
=======
[![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
as tracked by [RStudio CRAN mirror](http://cran-logs.rstudio.com/)

CRAN downloads in the last month:
[![](http://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)

CRAN downloads in the last week:
[![](http://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
>>>>>>> origin/master


============
### Branches

- The default branch (master) hosts the current development release (version 0.7.0) of the survival bump hunting procedure that implements the case of a survival response. At this point, this version is also restricted to a directed peeling search of the first box covered by the recursive coverage (outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm (Dazard et al., 2014, 2015, 2016). New features will be added soon as they are available. 

<<<<<<< HEAD
The package relies on an optional variable screening (pre-selection) procedure that is run before the PRSP algorithm and final variable usage (selection) procedure is done. This is done by four possible cross-validated variable screening (pre-selection) procedures offered to the user from the main end-user survival Bump Hunting function `sbh()`. At this point, the user can choose between:
=======
The main function relies on an optional variable pre-selection procedure that is run before the PRSP algorithm. 
At this point, this is done by a cross-validated penalization of the partial likelihood using the R package [`glmnet`](https://CRAN.R-project.org/package=glmnet).
>>>>>>> origin/master

   + Univariate Patient Recursive Survival Peeling algorithm (default of package `PRIMsrc`)
   + Penalized Censored Quantile Regression (by Semismooth Newton Coordinate Descent algorithm adapted from package `hqreg`)
   + Penalized Partial Likelihood (adapted from package `glmnet`)
   + Supervised Principal Component Analysis (adapted from package `superpc`)
   
In this version, the Cross-Validation (CV) procedure and Bump Hunting procedures that control model size (#covariates) and model complexity (#peeling steps), respectively, to fit the Survival Bump Hunting model, are carried out internally by two consecutive tasks within a single main end-user survival Bump Hunting function `sbh()`. The returned S3-class `sbh` object contains cross-validated estimates of all the decision-rules of used covariates and all other statistical quantities of interest at each iteration of the peeling sequence (inner loop of the PRSP algorithm). This enables the graphical display of results of profiling curves for model selection/tuning, peeling trajectories, covariate traces and survival distributions (see companion papers Dazard et al., 2014, 2015, 2016 for details). 

The package `PRIMsrc` offers a number of options for the number of replications of the fitting procedure to be perfomed: \eqn{B}; the type of \eqn{K}-fold cross-validation desired: (replicated)-averaged or-combined; as well as the peeling and cross-validation critera for model selection/tuning, and a few more parameters for the PRSP algorithm. The package takes advantage of the R packages `parallel` and `snow`, which allows users to create a parallel backend within an R session, enabling access to a cluster of compute cores and/or nodes on a local and/or remote machine(s) with either. The package supports two types of communication mechanisms between master and worker processes: 'Socket' or  'Message-Passing Interface' ('MPI').

- The second branch (unified) will host the future complete version of the code (version 1.0.0), including undirected peeling search derived from the Patient Rule Induction Method (PRIM), that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).


================
### Requirements

PRIMsrc 0.7.0 requires R-3.0.2 (2013-09-25). It was built and tested under R version 3.4.0 (2017-04-21) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:
[![Build Status](https://travis-ci.org/jedazard/PRIMsrc.png?branch=master)](https://travis-ci.org/jedazard/PRIMsrc)

See CRAN checks:
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PRIMsrc)](https://cran.r-project.org/web/checks/check_results_PRIMsrc.html).


================
### Installation

* To install [`PRIMsrc` from CRAN repository](https://CRAN.R-project.org/package=PRIMsrc), simply download and install the current version (0.7.0) from the CRAN repository:

```{r}
install.packages("PRIMsrc")
```

* Alternatively, you can install the most up-to-date development version (>= 0.7.0) of [`PRIMsrc` from GitHub repository](https://github.com/jedazard/PRIMsrc) using devtools, simply run:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jedazard/PRIMsrc")
```

=========
### Usage

* To load the PRIMsrc library in an R session and start using it:

```{r}
library("PRIMsrc")
```

* Version 0.7.0 is a major release with significant user-visible changes.
Check details of new features, changes, and bug fixes with the following R command:

```{r}
PRIMsrc.news()
```

* Check on how to cite the package with the R command:

```{r}
citation("PRIMsrc")
```

etc...


==================
### Website - Wiki

- See Project [Website](http://www.primsrc.com) for General Remarks, Goal and Why Use PRIMsrc.
- See Project [Wiki](https://github.com/jedazard/PRIMsrc/wiki) for Roadmap, Documentation and Examples, Publications, Case Studies, Support and How to contribute (code and documentation).


===================
### Acknowledgments

Authors: 
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)
   + Michael Choe, M.D. [(mjc206@case.edu)](mjc206@case.edu)
   + Michael LeBlanc, Ph.D. [(mleblanc@fhcrc.org)](mleblanc@fhcrc.org)
   + Alberto Santana, MBA. [(ahs4@case.edu)](ahs4@case.edu)

Maintainers: 
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)

Funding/Provision/Help:   
   + This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
   + This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.


==============
### References

   + Dazard J-E. and Rao J.S. (2017). 
      *Variable Selection Strategies for High-Dimensional Survival Bump Hunting using Recursive Peeling Methods*. 
      (in prep).
   + Yi C. and Huang J. (2016).
      *Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression*. 
      J. Comp Graph. Statistics, DOI: 10.1080/10618600.2016.1256816.
   + Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2016). 
      *Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods*. 
      Statistical Analysis and Data Mining, 9(1):12-42. 
   + Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015). 
      *R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification*. 
      In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA. 
      American Statistical Association IMS - JSM, p. 650-664. 
   + Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2014).
      *Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods*. 
      In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA. 
      American Statistical Association IMS - JSM, p. 3366-3380. 
   + Dazard J-E. and J.S. Rao (2010).
      *Local Sparse Bump Hunting*. 
      J. Comp Graph. Statistics, 19(4):900-92.
