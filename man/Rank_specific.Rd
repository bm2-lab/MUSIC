\name{Rank_specific}
\alias{Rank_specific}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the topic-specific perturbation effect ranking list.
}
\description{
MASCOT prioritizes the gene perturbation effect either as an overall perturbation effect on the cells, or in a functional topic-specific way. Here, the function calcultes and obtains the topic-specific perturbation effect ranking list.
}
\usage{
Rank_specific(distri_diff, output = FALSE, file_path = "./rank_specific.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distri_diff}{
  A dataframe showing the topic distribution difference between case and control which can obtain from "Diff_topic_distri()" in this package.
}
  \item{output}{
  FALSE by default. If it is TRUE, then the result will be outputted to a file.
}
  \item{file_path}{
  The save path of the result. It works only when the parameter "output" is TRUE.
}
}
\value{
A dataframe, the topic-specific perturbation effect ranking list 
}
\author{
Bin Duan
}
