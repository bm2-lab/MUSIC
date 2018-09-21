\name{Data_imputation}
\alias{Data_imputation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data imputation by SAVER.
}
\description{
To improve the quality of single cell CRISPR screening data, MASCOT performed data imputation by SAVER.
}
\usage{
Data_imputation(expression_profile, cpu_num = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after performing cell quality control,that is, the function of "Cell_qc()" in this package.
}
  \item{cpu_num}{
  The cpu number for parallel computation.
}
\value{
A dataframe showing the expression profile after performing this step.
}
\references{
Huang, M. et al. SAVER: gene expression recovery for single-cell RNA sequencing. Nature methods 15, 539-542 (2018).
}
\author{
Bin Duan
}
