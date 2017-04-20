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
[![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
as tracked by [RStudio CRAN mirror](http://cran-logs.rstudio.com/)

CRAN downloads in the last month:
[![](http://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)

CRAN downloads in the last week:
[![](http://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)


============
### Branches

The default branch (master) hosts the current development release (version 
0.6.3) of the survival bump hunting procedure that implements the case of a 
survival response. At this point, this version is also restricted to a 
directed peeling search of the first box covered by the recursive coverage 
(outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm. New 
features will be added soon as they are available. 

The main function relies on an optional variable pre-selection procedure that is run before the PRSP algorithm. 
At this point, this is done by a cross-validated penalization of the partial likelihood using the R package [`glmnet`](https://CRAN.R-project.org/package=glmnet).

In this version, the bump hunting procedure and the cross-validation procedures that control the model size and model peeling length are carried out by two separate procedures within a single main function `sbh()` that generates a unique S3-class object `PRSP`.  

- The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. 
Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. 
Specifically, it does a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion of choice: 
Concordance Error Rate (CER), Log-Rank Test (LRT), or Log-Hazard Ratio (LHR) (see companion paper below for details).

- The second branch (unified) will host the future complete version of the code (version 1.0.0), including undirected peeling search by Patient Rule Induction Method (PRIM) that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).


================
### Requirements

PRIMsrc 0.6.3 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-11-04 r69597) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:
[![Build Status](https://travis-ci.org/jedazard/PRIMsrc.png?branch=master)](https://travis-ci.org/jedazard/PRIMsrc)

See CRAN checks:
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PRIMsrc)](https://cran.r-project.org/web/checks/check_results_PRIMsrc.html).


================
### Installation

* To install PRIMsrc from CRAN, simply download and install the current version (0.6.3) from the CRAN repository:

```{r}
install.packages("PRIMsrc")
```

* Alternatively, you can install the most up-to-date development version (0.6.3) from GitHub, using devtools:

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

* Check the package news with the R command:

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
