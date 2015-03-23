=======
PRIMsrc
=======
Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. Assumptions are that the multivariate input variables can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event, possibly censored (Survival). It is intended to handle low and high-dimensional multivariate datasets (p < n and p > n cases), including the situation where the number of variables far exceeds that of samples (p >> n paradigm).

The current version is a developmental release that only includes the search of the first box of the recursive coverage loop and the case of a survival response. Ultimately, it will include all the features described above. New features will be added soon as they are available. At this time, the main function for fitting a Survival Bump Hunting Model in high-dimensional data relies on an internal cross-validated fitting procedure by Elasticnet Regularized Cox-Regression from the R package 'glmnet'. This variable selection procedure is made optional and recommended only for situation where the number of variables exceeds that of samples (p > n and p >> n cases).  

See also below the package news with the R command: PRIMsrc.news().

PRIMsrc is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

==========
References
==========
The companion papers (accepted and submitted for publication) can be accessed here:

ASA-IMS JSM Proceedings (2014): 
https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf

Archives arXiv:
http://arxiv.org/abs/1501.03856.

See also below on how to cite the package with the R command: citation("PRIMsrc").

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PRIMsrc.pdf) details the end-user functions. At this stage and for simplicity, there is a unique end-user main function for fitting a cross-validated Survival Bump Hunting model (sbh(...)). There are 5 end-user S3 generic plotting functions (plot.****(...)) along with two S3 generic functions: summary(...) and predict(...). Available are also 5 synthetic datasets and 2 real datasets including altogether low and high-dimensional situations (for p < n, p > n and p >> n cases). See the "PRIMsrc-package" introduction section of the manual for more details and examples.

=============
Installation: 
=============
PRIMsrc 0.5.5 was built under R version 3.0.2 (2013-09-25).
Installation has been tested on Windows, Linux and Mac platforms.
To install PRIMsrc, simply type:

library(devtools)

devtools::install_github("jedazard/PRIMsrc")

=============
Usage: 
=============
To load the PRIMsrc library in an R session and start using it:

library("PRIMsrc")

PRIMsrc.news()

citation("PRIMsrc")

etc...
