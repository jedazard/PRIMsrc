\name{Real.2}

\docType{data}

\alias{Real.2}

\title{Real Dataset #2: Genomic Dataset (\eqn{p >> n} case)}

\description{
  Publicly available lung cancer genomic data from the Chemores Cohort Study. This data is part of an integrated study 
  of mRNA, miRNA and clinical variables to characterize the molecular distinctions between squamous cell carcinoma (SCC) 
  and adenocarcinoma (AC) in Non Small Cell Lung Cancer (NSCLC) aside large cell lung carcinoma (LCC). Tissue samples 
  were analysed from a cohort of 123 patients, who underwent complete surgical resection at the Institut Mutualiste 
  Montsouris (Paris, France) between 30 January 2002 and 26 June 2006. The studied outcome was the "Disease-Free Survival Time". 
  Patients were followed until the first relapse occurred or administrative censoring. In this genomic dataset, 
  the expression levels of Agilent miRNA probes (\eqn{p=939}) were included from the \eqn{n=123} cohort samples. 
  The miRNA data contains normalized expression levels. See below the paper by Lazar et al. (2013) and Array Express 
  data repository for complete description of the samples, tissue preparation, Agilent array technology, and data normalization. 
  In addition to the genomic data, five clinical variables, also evaluated on the cohort samples, are included as 
  continuous variable ('Age') and nominal variables ('Type','KRAS.status','EGFR.status','P53.status'). 
  This dataset represents a situation where the number of covariates dominates the number of complete observations, 
  or \eqn{p >> n} case.  See Lazar et al. (2013) and the \href{http://www.chemores.ki.se/index.html}{CHEMORES Consortium} 
  website for more details.
}

\usage{Real.2}

\format{
  Dataset consists of a \code{numeric} \code{data.frame} containing \eqn{n=123} complete observations (samples) 
  by rows and \eqn{p=939} genomic covariates by columns, not including the censoring indicator and (censored) time-to-event variables. 
  It comes as a compressed Rda data file.
}

\source{See real data application in Dazard et al., 2015.}

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
        \item Dazard J-E. and Rao J.S. (2018).
              "\emph{Variable Selection Strategies for High-Dimensional Survival Bump Hunting using Recursive Peeling Methods.}"
              (in prep).
        \item Rao J.S., Huilin Y. and Dazard J-E. (2018).
              "\emph{Disparity Subtyping: Bringing Precision Medicine Closer to Disparity Science.}"
              (in prep).
        \item Diaz-Pachon D.A., Saenz J.P., Dazard J-E. and Rao J.S. (2018).
              "\emph{Mode Hunting through Active Information.}"
              (in press).
        \item Diaz-Pachon D.A., Dazard J-E. and Rao J.S. (2017).
              "\emph{Unsupervised Bump Hunting Using Principal Components.}"           
              In: Ahmed SE, editor. Big and Complex Data Analysis: Methodologies and Applications. 
              Contributions to Statistics, vol. Edited Refereed Volume. 
              Springer International Publishing, Cham Switzerland, p. 325-345.
        \item Yi C. and Huang J. (2017).
              "\emph{Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression}."
              J. Comp Graph. Statistics, 26(3):547-557.
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2016).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              Statistical Analysis and Data Mining, 9(1):12-42.
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification.}"
              In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA.
              American Statistical Association IMS - JSM, p. 650-664.
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
        \item Dazard J-E. and J.S. Rao (2010).
              "\emph{Local Sparse Bump Hunting.}"
              J. Comp Graph. Statistics, 19(4):900-92.
        \item Lazar V. et al. (2013).
               "\emph{Integrated molecular portrait of non-small cell lung cancers.}"
              BMC Medical Genomics 6:53-65.
    }
}

\section{Acknowledgments}{
    This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
    This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\keyword{Real Dataset}
\keyword{Tumor sample comparisons}

\seealso{
    \itemize{
       Array Express data repository at the European Bioinformatics Institute.
       Accession number: #E-MTAB-1134 (MIR): \url{www.ebi.ac.uk/arrayexpress/}
    }
    \itemize{
       CHEMORES Consortium website: \url{http://www.chemores.ki.se/index.html}
    }
}
