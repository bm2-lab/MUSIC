\name{Data_imputation}
\alias{Data_imputation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data imputation by SAVER.
}
\description{
To improve the quality of single cell CRISPR screening data, MUSIC performed data imputation by SAVER.
}
\usage{
Data_imputation(expression_profile, perturb_information, cpu_num = 4, split = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after performing cell quality control,that is, the function of "Cell_qc()" in this package.
}
  \item{perturb_information}{
  A character vector showing the information of sample after performing the function of "Cell_qc()" in this package.
}
  \item{cpu_num}{
  The cpu number for parallel computation. The default is 4. Parallel computation is strongly recommeneded to use because this step may take long time without parallel computation.
}
  \item{split}{
  The number of data block. The default is 2. When the data is too big for your computer, you can divide the data into multiple blocks and then calculate. The more blocks, the more space is saved and the longer time will take.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{expression_profile }{The expression profile after performing this step.}
  \item{perturb_information }{The information (perturbation names and sample names) of cells retained after performing this step.}
}
\references{
Huang, M. et al. SAVER: gene expression recovery for single-cell RNA sequencing. Nature methods 15, 539-542 (2018).
}
\author{
Bin Duan
}
