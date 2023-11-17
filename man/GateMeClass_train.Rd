\name{GateMeClass_train}
\alias{GateMeClass_train}
\title{Extracting the marker table from a reference annotated dataset}
\usage{

GateMeClass_train(reference,
    labels,
    RSS = T,
    GMM_parameterization = "V",
    verbose = T, 
    seed = 1)
}

\arguments{
  \item{reference}{
  Expression matrix of the reference annotated dataset
  }
  \item{labels}{
  Labels of the reference dataset, class = character
  }
  \item{RSS}{
  RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution
  }
  \item{GMM_parameterization}{
  GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal)
  }
  \item{verbose}{
  Show output information
  }
  \item{seed}{
  Seed for randomization
  }
}
\examples{
For usage examples refer to https://github.com/simo1c/GateMeClass
}



