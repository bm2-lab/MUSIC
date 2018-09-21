\name{Select_topic_number}
\alias{Select_topic_number}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To select the optimal topic number.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Select_topic_number(topic_model_list, plot = FALSE, plot_path = "./select_topic_number.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{topic_model_list}{
  A list of "LDA" class with different topic number. It can be obtained by the function "Get_topics()" in this package.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
  The save path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\value{
An object, "LDA" class with the optimal topic number selected. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Bin Duan
}
