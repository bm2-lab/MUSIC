\name{Diff_topic_distri}
\alias{Diff_topic_distri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcualting topic distribution difference between case and control.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Diff_topic_distri(model, perturb_information, plot = FALSE, plot_path = "./distribution_of_topic.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  Object of class "LDA" with the optimal topic number selected.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample after all the filterings.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
  The save path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\value{
A dataframe.
}
\author{
Bin Duan
}
