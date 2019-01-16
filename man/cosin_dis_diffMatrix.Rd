\name{cosin_dis_diffMatrix}
\alias{cosin_dis_diffMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parallel calculation for cosin similarity between two matrices.
}
\description{
To calculate the cosine simularity between any row between two matrices with parallel calculation strategy.
}
\usage{
cosin_dis_diffMatrix(matrix1, matrix2, cpu_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix1}{
  A matrix. Taking gene expression profile as example, the column is sample and the row is gene. 
}
  \item{matrix2}{
  A matrix. Taking gene expression profile as example, the column is sample and the row is gene. 
}
  \item{cpu_num}{
  The cpu number for parallel computation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix to show the cosin similarity between any points between two matrices.
}

\author{
Bin Duan
}
