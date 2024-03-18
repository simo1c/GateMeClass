\name{GateMeClass_train}
\alias{GateMeClass_train}
\title{Extracting the marker table from a reference annotated dataset}
\description{The training module of GateMeClass inputs an annotated cytometry dataset in the form of an expression matrix and outputs a marker table with a pseudo gating strategy for each cell type}
\usage{
GateMeClass_train(
    reference,
    labels,
    GMM_parameterization = "E",
    verbose = TRUE, 
    seed = 1)
}

\arguments{
  \item{reference}{
  The expression matrix of the reference annotated dataset.
  }
  \item{labels}{
  A character vector with the labels of the reference dataset.
  }
  \item{GMM_parameterization}{
  A character vector with the GMM (Gaussian-Mixture-Model) parameter: "V" (Variable) or "E" (Equal).
  }
  \item{verbose}{
  If TRUE show output information.
  }
  \item{seed}{
  Seed for randomization.
  }
}
\value{
    The function returns a data.frame containing the marker table extracted.
}
\examples{
# For examples refer to https://github.com/simo1c/GateMeClass
}
