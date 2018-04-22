\name{Rank_diff}
\alias{Rank_diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining important knockouts between two different conditions.
}
\description{
Obtaining knockouts with highly different ranking between two different conditions for the same type cellls.
}
\usage{
Rank_diff(rankOverall_result_summary_condition1, rankOverall_result_summary_condition2, difference_threshold = 0.5, plot = TRUE, plot_path = "~/rank_diff.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rankOverall_result_summary_condition1}{
  The overall perturbation impact ranking list when knockout cells are under condition 1.
}
  \item{rankOverall_result_summary_condition2}{
  The overall perturbation impact ranking list when knockout cells are under condition 2.
}
  \item{difference_threshold}{
  The threshold ratio of the ranking difference in all knockouts between two conditions. The default is 0.5 .
}
  \item{plot}{
  TURE by default. If FAALSE, don't plot the graph.
}
  \item{plot_path}{
  The path of graph you plot. It works only when "plot" is TURE.
}
}
\value{
A dataframe showing knockouts with the highly different ranking and their corresponding ranking difference.
}
\author{
Bin Duan
}
