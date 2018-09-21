\name{Diff_perturb_effect}
\alias{Diff_perturb_effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating perturbation effect difference under different experimental conditions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Diff_perturb_effect(rank_overall_condition_1, rank_overall_condition_2, fold_change = 2, plot = FALSE, plot_path = "./perturbation_effect_difference.pdf", output = FALSE,file_path = "./perturbation_effect_difference.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank_overall_condition_1}{
  A dataframe. The overall perturbation effect ranking list, under condition 1. It can be obtained bu function "Rank_overall()" in this package.
}
  \item{rank_overall_condition_2}{
  A dataframe. The overall perturbation effect ranking list, under condition 2. It can be obtained bu function "Rank_overall()" in this package.
}
  \item{fold_change}{
  The cutoff to show the difference. 2 by default.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
  The save path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
  \item{output}{
  FALSE by default. If it is TRUE, then the result will be outputted to a file.
}
  \item{file_path}{
  The save path of the result. It works only when the parameter "output" is TRUE.
}
}
\value{
A dataframe.
}
\author{
Bin Duan
}
