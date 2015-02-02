=======
PrimSRC
=======
Performs a unified treatment of Bump Hunting by Patient Rule Induction Method in Survival, Regression and Classification settings (SRC). The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. The region is shaped as a hyperdimensional box that is not necessarily contiguous.
Assumptions are that the multivariate input variables can be discrete or continuous and the
univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event,
possibly censored (Survival). It is intended to handle high-dimensional multivariate datasets,
where the number of variables far exceeds that of the samples (p >> n paradigm).

The current version is the development release that only includes the case of low dimensional
situations and a survival response. Ultimately, it will include all the features described above. New
features will be added soon as they are available.

See also below the package news with the R command: PrimSRC.news().

PrimSRC is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

==========
References
==========
The companion papers (accepted and submitted for publication) can be accessed here:

ASA-IMS JSM Proceedings (2014): 
https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf

Archives arXiv:
http://arxiv.org/abs/1501.03856.

See also below on how to cite the package with the R command: citation("PrimSRC").

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PrimSRC.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only one end-user function (sbh), 5 end-user plotting functions (plot_****) and 2 end-user datasets (synthetic and real). See the "PrimSRC-package" introduction section of the manual for more details and examples of use.

=============
Installation: 
=============
PrimSRC 0.3.0 was built under R version 3.1.2 (2014-10-31).
Installation has been tested on Windows, Linux and Mac platforms.
To install this development version 0.3.0, simply type:

library(devtools)

devtools::install_github("jedazard/PrimSRC")

=============
Usage: 
=============
To load the PrimSRC library in an R session and start using it:

library("PrimSRC")

PrimSRC.news()

citation("PrimSRC")

etc...
