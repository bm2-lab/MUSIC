\name{Correlation_perturbation}
\alias{Correlation_perturbation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the relationships between different perturbations.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Correlation_perturbation(distri_diff, cutoff = 0.9, gene = "all", plot = FALSE, plot_path = "./correlation_network.pdf", output = FALSE, file_path = "./correlation_perturbation.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distri_diff}{
  A dataframe . This dataframe can be obtained by the function "Get_distri_diff()" in this package.
}
  \item{cutoff}{
  The cutoff of correlation which will be used to plot the network graph. The default is 0.9.
}
  \item{gene}{
  A character or vector. The perturbation choosed to plot the network graph. The default is "all" which means plot all perturbations.
}
  \item{plot}{
  FALSE by default. If it is TRUE, then plot the network graph based on the perturbation correlation.
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
A dataframe showing the perturbation correlation.
}
\author{
Bin Duan
}
