\name{plot_profile}

\alias{plot_profile}

\title{Visualization for Model Selection/Tuning}

\description{
    Function for plotting the cross-validated model selection/tuning profiles of a \code{sbh} object. 
    It uses the user's choice of cross-validation criterion statistics among the Log Hazard Ratio (LHR), 
    Log-Rank Test (LRT) or Concordance Error Rate (CER). The function plots (as it applies) both profiles 
    of cross-validation criterion as a function of variables screening size (cardinal subset of top-screened 
    variables in the PRSP variable screening procedure), and peeling length (number of peeling steps
    of the peeling sequence in the inner loop of the PRSP or PRGSP algorithm).
}

\usage{
  plot_profile(object,
               main = "Profile Plots",
               xlim = NULL,
               ylim = NULL,
               add.sd = TRUE, 
               add.profiles = TRUE,
               add.caption = TRUE, 
               text.caption = c("Mean","Std. Error"),
               pch = 20, 
               col = 1, 
               lty = 1, 
               lwd = 0.5, 
               cex = 0.5,
               device = NULL, 
               file = "Profile Plots", 
               path = getwd(), 
               horizontal = FALSE, 
               width = 8.5, 
               height = 5.0, ...)
}

\arguments{
  \item{object}{
    Object of class \code{sbh} as generated by the main function \code{\link[PRIMsrc]{sbh}}.
  }
  \item{main}{
    \code{Character} \code{vector}. Main Title. 
    Defaults to "Profile Plots".
  }
  \item{xlim}{
    \code{Numeric} \code{vector} of length 2. The x limits [x1, x2] of the plot. 
    Defaults to \code{NULL}. 
  }
  \item{ylim}{
    \code{Numeric} \code{vector} of length 2. The y limits [y1, y2] of the plot. 
    Defaults to \code{NULL}. 
  }
  \item{add.sd}{
    \code{Logical} scalar. Shall the standard error bars be plotted? 
    Defaults to \code{TRUE}.
  }
  \item{add.profiles}{
    \code{Logical} scalar. Shall the individual profiles (for all replicates) be plotted?
    Defaults to \code{TRUE}.
  }
  \item{add.caption}{
    \code{Logical} scalar. Should the caption be plotted? 
    Defaults to \code{TRUE}.
  }
  \item{text.caption}{
    \code{Character} \code{vector} of caption content. 
    Defaults to \{"Mean","Std. Error"\}.
  }
  \item{pch}{
    \code{Integer} scalar of symbol number for all the profiles. 
    Defaults to 20.
  }
  \item{col}{
    \code{Integer} scalar of line color of the mean profile. 
    Defaults to 1.
  }
  \item{lty}{
    \code{Integer} scalar of line type of the mean profile. 
    Defaults to 1.
  }
  \item{lwd}{
    \code{Numeric} scalar of line width of the mean profile. 
    Defaults to 0.5.
  }
  \item{cex}{
    \code{Numeric} scalar of symbol expansion for all the profiles. 
    Defaults to 0.5.
  }
  \item{device}{
    Graphic display device in \{\code{NULL}, "PS", "PDF"\}. 
    Defaults to \code{NULL} (standard output screen).
    Currently implemented graphic display devices are "PS" (Postscript) or "PDF" (Portable Document Format).
  }
  \item{file}{
    File name for output graphic. Defaults to "Profile Plot".
  }
  \item{path}{
    Absolute path (without final (back)slash separator). 
    Defaults to working directory path.
  }
  \item{horizontal}{
    \code{Logical} scalar. Orientation of the printed image. 
    Defaults to \code{FALSE}, that is potrait orientation.
  }
  \item{width}{
    \code{Numeric} scalar. Width of the graphics region in inches. 
    Defaults to 8.5.
  }
  \item{height}{
    \code{Numeric} scalar. Height of the graphics region in inches. 
    Defaults to 5.0.
  }
  \item{\dots}{
    Generic arguments passed to other plotting functions.
  }
}

\details{
  Model tuning is done by applying the cross-validation criterion defined by the user's choice of specific statistic.
  The goal is to find the optimal value of model parameters by maximization of LHR or LRT, or minimization of CER.
  The parameters to optimize are (i) the cardinal of top-ranked variables subsets (if the "prsp" variable screening is chosen), 
  and (ii) the number of peeling steps of the peeling sequence (inner loop of our PRSP algorithm) 
  in any case of variable screening method.
  
  Currently, this is done internally for visualization purposes, but it will ultimately offer the option to be done
  interactively with the end-user as well for parameter choosing/model selection.
}

\value{
  Invisible. None. Displays the plot(s) on the specified \code{device}.
}

\author{
    \itemize{
        \item "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
        \item "Michael Choe, M.D." \email{mjc206@case.edu}
        \item "Michael LeBlanc, Ph.D." \email{mleblanc@fhcrc.org}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
        \item "J. Sunil Rao, Ph.D." \email{Rao@biostat.med.miami.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
}

\references{
    \itemize{
        \item Dazard J-E. and Rao J.S. (2021a).
              "\emph{Variable Selection Strategies for High-Dimensional Recursive Peeling-Based Survival Bump Hunting Models.}"
              (in prep).
        \item Dazard J-E. and Rao J.S. (2021b).
              "\emph{Group Bump Hunting by Recursive Peeling-Based Methods: Application to Survival/Risk Predictive Models.}"
              (in prep).
        \item Dazard J-E., Choe M., Pawitan Y., and Rao J.S. (2021c).
              "\emph{Identification and Characterization of Informative Prognostic Subgroups by Survival Bump Hunting.}"
              (in prep).
        \item Rao J.S., Huilin Y., and Dazard J-E. (2020).
              "\emph{Disparity Subtyping: Bringing Precision Medicine Closer to Disparity Science.}"
              Cancer Epidemiology Biomarkers & Prevention, 29(6 Suppl):C018.
        \item Yi C. and Huang J. (2017).
              "\emph{Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression}."
              J. Comp Graph. Statistics, 26(3):547-557.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2016).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              Statistical Analysis and Data Mining, 9(1):12-42.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2015).
              "\emph{R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification.}"
              In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA.
              American Statistical Association IMS - JSM, p. 650-664.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
    }
}

\section{Acknowledgments}{
    This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
    This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\note{
  End-user plotting function.
}

\keyword{Exploratory Survival/Risk Analysis}
\keyword{Survival/Risk Estimation & Prediction}
\keyword{Non-Parametric Method}
\keyword{Cross-Validation}
\keyword{Bump Hunting}
\keyword{Rule-Induction Method}
