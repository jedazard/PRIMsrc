=======
PrimSRC
=======
This is the development release of the future PrimSRC CRAN package for Bump Hunting by Patient Rule Induction Method in Survival, Regression and Classification settings. The companion paper (submitted for publication) can be accessed here: http://arxiv.org/abs/1501.03856.

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (PrimSRC.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only one end-user function (sbh) and 5 end-user plotting functions (plot_****). See the "PrimSRC-package" introduction section of the manual for more details and examples of use.

=============
Installation: 
=============
The latest R version 3.1.2 (2014-10-31) is recommended.
Installation has been tested on Windows, Linux and Mac platforms.
To install the software and load the PrimSRC library in an R session, simply type:

library(devtools)

devtools::install_github("jedazard/PrimSRC")

library("PrimSRC")

