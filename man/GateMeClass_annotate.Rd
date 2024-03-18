\name{GateMeClass_annotate}
\alias{GateMeClass_annotate}
\title{Annotation of a cytometry dataset using a marker table}
\description{The annotation module of GateMeClass inputs a dataset in the form of a matrix and a marker table and returns a vector of cell annotations}
\usage{
GateMeClass_annotate(
    exp_matrix,
    marker_table,
    train_parameters = = list(reference = NULL, labels = NULL),
    GMM_parameterization = "E",
    reject_option = FALSE,
    sampling = 0.1,
    k = 20,				
    verbose = TRUE,
    narrow_marker_table = TRUE,
    seed = 1)
}
\arguments{
  \item{exp_matrix}{
  An expression matrix.
  }
  \item{marker_table}{
  A data.frame with a manually curated or extracted marker table.
  }
  \item{train_parameters}{
  A list with the parameters for the training function.
  }
  \item{GMM_parameterization}{
  A character vector with the GMM (Gaussian-Mixture-Model) parameter: "V" (Variable) or "E" (Equal).
  }
  \item{reject_option}{
  If TRUE this parameter tries to detect cell types not defined in the marker table using MNN algorithm.
  }
  \item{sampling}{
  Percentage of the cells used for the annotation.
  }
  \item{k}{
  k parameter of k-NN (k-Nearest-Neighbour) used to refine uncertain labels to the most similar already annotated.
  }
  \item{narrow_marker_table}{
  Format of marker table. TRUE for using a marker table with two columns, "Cell" with the name of cell type and "Gate" with the gating strategy (e.g., CD3+CD4+). FALSE for using a column for each marker.
  }
  \item{verbose}{
  If TRUE show output information.
  }
  \item{seed}{
  Seed for randomization.
  }
}

\value{
    The function returns a list with the following elements:
    \item{labels}{
        A character vector with the labels returned by GateMeClass.
    }    
    \item{marker_table}{
        A data.frame with the marker table used by GateMeClass.
    }    
    \item{cell_signatures}{
        A data.frame with the marker signatures of each cell attributed by GateMeClass, NA for cells unsampled.
    }
}
\examples{
# For examples refer to https://github.com/simo1c/GateMeClass
}
