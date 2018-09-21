\name{Input_preprocess}
\alias{Input_preprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The preprocess of input data.
}
\description{
Integrating the expression profile and corresponding sample information,
}
\usage{
Input_preprocess(expression_profile, perturb_information, batch_information = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the original expression profile whose Row is gene and column is sample.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample. The vector's name is the sample name.
}
  \item{batch_information}{
  NULL by default. A character vector showing the batch of each sample. The vector's name is the sample name.
}
}
\value{
A list containing "expression profile", "perturb_information" and "batch_information" if it is not NULL.
}
\author{
Bin Duan
}
