# PRIMsrc
Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification in a multivariate setting and in high-dimensional data.


===============
### Description

The general problem in "Bump Hunting" (BH) is to identify, characterize and predict hidden structures in the data that are informative and 
significant. In practice, "Bump Hunting" refers to the task of mapping out a local region of the multi-dimensional input space where a target 
function of interest, usually unknown, assumes smaller or larger values than its average over the entire space. In general, the region could 
be any smooth shape (e.g. a convex hull) possibly disjoint.
    
`PRIMsrc` implements a unified treatment of "Bump Hunting" (BH) by algorithms derived from the Patient Rule Induction Method (PRIM) 
(Friedman and Fisher, 1999) for Survival, Regression and Classification outcomes (SRC). To estimate the region, PRIM generates decision rules 
delineating hyperdimensional boxes (hyperrectangles) of the input space, not necessarily contiguous, where the outcome is smaller or larger 
than its average over the entire space.

Assumptions are that the multivariate input variables can be discrete or continuous, and the univariate outcome variable can be discrete 
(Classification), continuous (Regression), or a time-to-event, possibly censored (Survival). It is intended to handle low and high-dimensional 
multivariate datasets, including the paradigm where the number of covariates (_p_) exceeds or dominates that of samples (_n_): _p_ > _n_ or 
_p_ >> _n_.    

Please note that the current version is a development release, that only implements the case of a survival outcome. At this point, 
this version of `PRIMsrc` is also restricted to a directed peeling search of the first box covered by the recursive coverage (outer) 
loop of our PRSP or PRGSP algorithm (see details below). New features will be added as soon as available. 


===============
### Details

In a direct application, "Bump Hunting" (BH) can identify subgroups of observations for which their outcome is as extreme as possible. 
Similarly to this traditional goal of subgroup finding, `PRIMsrc` also implements the goal of mapping out a region (possibly disjointed) 
of the input space where the outcome _difference_ between existing (fixed) groups of observations is as extreme as possible. We refer to the 
later goal as "Group Bump Hunting" (GBH).

In the case of a time-to event outcome, possibly censored (as in survival or risk analysis), "Survival Bump Hunting" (SBH) is done by our 
Patient Recursive Survival Peeling (PRSP) algorithm. See Dazard and Rao (2014, 2015, 2016, 2021a) for details, as well as Dazard et al. (2021c)
for an application in Patient Survival Subtyping. Alternatively, "Group Survival Bump Hunting" (GSBH) is done by using a derivation of PRSP with 
specifc peeling and cross-validation criterion, called Patient Recursive Group Survival Peeling (PRGSP). See Dazard and Rao (2021b) for details, 
as well as Rao et al. (2020) for an application in Survival Disparity Subtyping.

The package relies on an optional variable screening (pre-selection) procedure that is run before the PRSP algorithm and final variable usage 
(selection) procedure is done. This is done by four possible cross-validated variable screening (pre-selection) procedures offered to the user 
from the main end-user survival Bump Hunting function `sbh()`. At this point, the user can choose between:

   + Univariate Patient Recursive Survival Peeling algorithm (default of package `PRIMsrc`)
   + Penalized Censored Quantile Regression (by Semismooth Newton Coordinate Descent algorithm adapted from package `hqreg`)
   + Penalized Partial Likelihood (adapted from package `glmnet`)
   + Supervised Principal Component Analysis (adapted from package `superpc`)
   
In this version, the Cross-Validation (CV) procedure and Bump Hunting procedures that control model size (#covariates) 
and model complexity (#peeling steps), respectively, to fit the Survival Bump Hunting model, are carried out internally by two consecutive 
tasks within a single main end-user survival Bump Hunting function `sbh()`. The returned S3-class `sbh` object contains cross-validated 
estimates of all the decision-rules of used covariates and all other statistical quantities of interest at each iteration of the peeling 
sequence (inner loop of the PRSP algorithm). This enables the graphical display of results of profiling curves for model selection/tuning, 
peeling trajectories, covariate traces and survival distributions.

The package `PRIMsrc` offers a number of options for the number of replications of the fitting procedure to be perfomed: _B_; 
the type of _K_-fold cross-validation desired: (replicated)-averaged or-combined; as well as the peeling and cross-validation critera 
for model selection/tuning, and a few more parameters for the PRSP algorithm. The package takes advantage of the 
R packages `parallel` and `snow`, which allows users to create a parallel backend within an R session, enabling access to a cluster 
of compute cores and/or nodes on a local and/or remote machine(s) with either. The package supports two types of communication mechanisms 
between master and worker processes: 'Socket' or  'Message-Passing Interface' ('MPI').


============
### Branches

This branch (master) is the  default one, that hosts the current development release (version 0.9.0) of the  "Survival Bump Hunting" (SBH)
(or "Group Survival Bump Hunting" (GSBH)) procedure. Note that `PRIMsrc` is still a non-production release and that version 0.9.0 implements 
significant user-visible changes. Check details of new features, changes, and bug fixes in the "Usage" section below.

The future branch (unified) will host the complete version of the code (version 1.0.0), including undirected peeling search derived 
from the  Patient Rule Induction Method (PRIM), and unified treatment of bump hunting for every type of common outcome: Survival, Regression, 
and Classification (SRC).


===========
### License

PRIMsrc is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](https://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](https://www.gnu.org/licenses/gpl-3.0.html).

=============
### Downloads

CRAN downloads since initial release to CRAN (2015-07-28):
[![](https://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)
as tracked by [RStudio CRAN mirror](http://cran-logs.rstudio.com/)

CRAN downloads in the last month:
[![](https://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)

CRAN downloads in the last week:
[![](https://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](https://CRAN.R-project.org/package=PRIMsrc)


================
### Requirements

PRIMsrc (>= 0.9.0) requires R-3.5.0 (2018-04-23). It was built and tested under R version 4.1.0 (2021-05-18) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:
[![Build Status](https://app.travis-ci.com/jedazard/PRIMsrc.svg?branch=master)](https://app.travis-ci.com/github/jedazard/PRIMsrc)

See CRAN checks:
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PRIMsrc)](https://cran.r-project.org/web/checks/check_results_PRIMsrc.html)


================
### Installation

* To install the stable version of `PRIMsrc`, simply download and install the current version (0.8.2) from the [CRAN](https://CRAN.R-project.org/package=PRIMsrc) 
repository:

```{r}
install.packages("PRIMsrc")
```

* Alternatively, you can install the most up-to-date development version (>= 0.9.0) of `PRIMsrc` from the [GitHub](https://github.com/jedazard/PRIMsrc) repository:

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

* Check details of new features, changes, and bug fixes with the following R command:

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

- See Project [Website](https://www.primsrc.com) for General Remarks, Goal and Why Use PRIMsrc.
- See Project [Wiki](https://github.com/jedazard/PRIMsrc/wiki) for Roadmap, Documentation and Examples, Publications, Case Studies, Support and How to contribute (code and documentation).


===================
### Acknowledgments

Authors: 
   + Jean-Eudes Dazard, Ph.D. <jean-eudes.dazard@case.edu>
   + Michael Choe, M.D. <mjc206@case.edu>
   + Michael LeBlanc, Ph.D. <mleblanc@fhcrc.org>
   + J. Sunil Rao, Ph.D. <JRao@biostat.med.miami.edu>
   + Alberto Santana, MBA. <ahs4@case.edu>

Maintainers: 
   + Jean-Eudes Dazard, Ph.D.<jean-eudes.dazard@case.edu>

Funding/Provision/Help:   
   + This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
   + This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.


==============
### References

   + Dazard J-E. and Rao J.S. 
      *Variable Selection Strategies for High-Dimensional Recursive Peeling-Based Survival Bump Hunting Models*. 
      [2021a (in prep)].
      
   + Dazard J-E. and Rao J.S. 
      *Group Bump Hunting by Recursive Peeling-Based Methods: Application to Survival/Risk Predictive Models*. 
      [2021b (in prep)].
      
   + Dazard J-E., Choe M., Pawitan Y., and Rao J.S. 
      *Identification and Characterization of Informative Prognostic Subgroups by Survival Bump Hunting*. 
      [2021c (in prep)].
      
   + Rao J.S., Huilin Y., and Dazard J-E. 
      *Disparity Subtyping: Bringing Precision Medicine Closer to Disparity Science*. 
      [Cancer Epidemiology Biomarkers & Prevention (2020), 29(6 Suppl):C018](https://cebp.aacrjournals.org/content/29/6_Supplement_1/C018).
      
   + Dazard J-E., Choe M., LeBlanc M., and Rao J.S. 
      *Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods*. 
      [Statistical Analysis and Data Mining (2016), 9(1):12-42](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11301). 
      (The American Statistical Association Data Science Journal)

   + Dazard J-E., Choe M., LeBlanc M., and Rao J.S. 
      *R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification*. 
      In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA. 
      American Statistical Association IMS - JSM.
      [JSM (2015), 650-664](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4718587/).

   + Dazard J-E., Choe M., LeBlanc M., and Rao J.S.
      *Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods*. 
      In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA. 
      American Statistical Association IMS - JSM. 
      [JSM (2014), 3366-3380](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4795911/).
      
   + Diaz-Pachon D.A., Saenz J.P., Rao J.S., and Dazard J-E. 
     *Mode Hunting Through Active Information*. 
     [Applied Stochastic Models in Business and Industry (2019), 35(2):376-393](https://onlinelibrary.wiley.com/doi/abs/10.1002/asmb.2430).

   + Diaz-Pachon D.A., Dazard J-E., and Rao J.S. 
     *Unsupervised Bump Hunting Using Principal Components*. 
     In: Ahmed SE, editor. Big and Complex Data Analysis: Methodologies and Applications. 
     Contributions to Statistics, vol. Edited Refereed Volume. 
     [Springer International Publishing, Cham Switzerland (2017), 325-345](https://link.springer.com/chapter/10.1007/978-3-319-41573-4_16).
            
   + Dazard J-E. and Rao J.S.
      *Local Sparse Bump Hunting*. 
      [J. Comp. Graph. Statistics (2010), 19(4):900-92](https://amstat.tandfonline.com/doi/abs/10.1198/jcgs.2010.09029).
   
   + Friedman J. and Fisher N.
      *Bump Hunting in High-Dimensional Data*. 
      [Stat. Computing (1999), 9:123-143](https://link.springer.com/article/10.1023/A:1008894516817).

   + Yi C. and Huang J.
      *Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression*. 
      [J. Comp. Graph. Statistics (2017), 26(3):547-557](https://amstat.tandfonline.com/doi/full/10.1080/10618600.2016.1256816#.W2ybVOhKiHs). 

