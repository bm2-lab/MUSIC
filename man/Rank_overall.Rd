\name{Rank_overall}
\alias{Rank_overall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the overall perturbation effect ranking list.
}
\description{
MASCOT prioritizes the gene perturbation effect either as an overall perturbation effect on the cells, or in a functional topic-specific way. Here, the function calcultes and obtains the overall perturbation effect ranking list.}
\usage{
Rank_overall(distri_diff, offTarget_hash = hash(), output = FALSE, file_path = "./rank_overall.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distri_diff}{
  A dataframe showing the topic distribution difference between case and control which can obtain from "Diff_topic_distri()" in this package.
}
  \item{offTarget_hash}{
  A hash table showing offtarget information which can obtain from "Get_Offtarget()". The default is a null hash.
}
  \item{output}{
  FALSE by default. If it is TRUE, then the result will be outputted to a file.
}
  \item{file_path}{
  The save path of the result. It works only when the parameter "output" is TRUE.
}
}
\value{
A dataframe, the overall perturbation effect ranking list.
}
\author{
Bin Duan
}
