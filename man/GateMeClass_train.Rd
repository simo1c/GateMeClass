\name{GateMeClass_train}
\alias{GateMeClass_train}
\title{Extracting the marker table from a reference annotated dataset}
\description{The training module of GateMeClass inputs an annotated cytometry dataset in the form of an expression matrix and outputs a marker table with a pseudo gating strategy for each cell type}
\usage{
GateMeClass_train(reference,
    labels,
    GMM_parameterization = "E",
    verbose = T, 
    seed = 1)
}

\arguments{
  \item{reference}{
  Expression matrix of the reference annotated dataset, class = matrix
  }
  \item{labels}{
  Labels of the reference dataset, class = character
  }
  \item{GMM_parameterization}{
  GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character
  }
  \item{verbose}{
  Show output information
  }
  \item{seed}{
  Seed for randomization
  }
}
\value{
    The function returns a data.frame containing the marker table extracted.
}
\examples{
# For examples refer to https://github.com/simo1c/GateMeClass
}
