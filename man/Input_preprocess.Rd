\name{Input_preprocess}
\alias{Input_preprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking the input data.
}
\description{
Integrating the expression profile and corresponding sample information, including each sample's knockout gene, detected sgRNA and batch information.
}
\usage{
Input_preprocess(expression_profile, sample_info_gene, sample_info_sgRNA, sample_info_batch = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{expression_profile}{
  A dataframe showing the original expression profile with sample names and gene names.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample. The vector's column name is the sample name.
}
  \item{sample_info_sgRNA}{
  A character vector showing the sgRNA detected in each sample. The vector's column name is the sample name.
}
  \item{sample_info_batch}{
  NULL by default. A character vector showing the batch of each sample. The vector's column name is the sample name.
}
}
\value{
A list containing "expression profile", "sample_info_gene", "sample_info_sgRNA" and "sample_info_batch" if it is not NULL.
}
\author{
Bin Duan
}
