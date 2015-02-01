=======
PrimSRC
=======
This is the development release of the future PrimSRC CRAN package for Bump Hunting by Patient Rule Induction Method in Survival, Regression and Classification settings. 

==========
References
==========
The companion papers (accepted and submitted for publication) can be accessed here (see also PrimSRC package citation below):

JSM Proceedings 2014 (ASA-IMS): https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf

Archives arXiv: http://arxiv.org/abs/1501.03856.

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PrimSRC.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only one end-user function (sbh), 5 end-user plotting functions (plot_****) and 2 end-user datasets (synthetic and real). See the "PrimSRC-package" introduction section of the manual for more details and examples of use.

=============
Installation: 
=============
PrimSRC 0.3.0 was built under the latest R version 3.1.2 (2014-10-31).
Installation has been tested on Windows, Linux and Mac platforms.
To install the software and load the PrimSRC library in an R session, simply type:

library(devtools)

devtools::install_github("jedazard/PrimSRC")

library("PrimSRC")

#PrimSRC package news:

PrimSRC.news()

#PrimSRC package citation:

citation("PrimSRC")
