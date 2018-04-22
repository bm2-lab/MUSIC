\name{Get_distribution_diff}
\alias{Get_distribution_diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcualting topic distribution difference between case and control.
}
\description{
With the result of topic model, we MASCOT can calcualte topic distribution difference between case and control.
}
\usage{
Get_distribution_diff(model, sample_info_gene, KO_efficiency, index = -4.35e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  Object of class "LDA" with the optimal topic number.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample after all the filterings.
}
  \item{index}{
   A constant fitted by simulation in our analysis. The result is -4.35e-05.
}
  \item{KO_efficiency}{
  A numeric vector showing the overall knockout efficiency of each knockout obtained from "Cell_filtering()".
}
}
\value{
A dataframe showing some information, including ratio difference, background difference and so on for each knockout.
}
\author{
Bin Duan
}
