\name{Plot_cell_topic}
\alias{Plot_cell_topic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting a heatmap showing the distribution of cells in each topic.
}
\description{
After obtaining topic distribution in each cell with certain topic number, a heatmap can directory show the distribution of cells in each topic.
}
\usage{
Plot_cell_topic(model, plot_path = "~/distribution_cell_in_topics.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  Object of class "LDA" with the optimal topic number.
}
  \item{plot_path}{
  The path of the graph.
}
}
\author{
Bin Duan
}
