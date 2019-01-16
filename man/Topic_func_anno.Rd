\name{Topic_func_anno}
\alias{Topic_func_anno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotating each topic's function.
}
\description{
After calculating topics, MASCOT obtains the occurring probabilities of genes available in each topic. For each topic, genes with large probabilities were used to annotate the topic’s function. Specifically, MASCOT firstly selects the n\% top-ranked (default 20\%) genes of each topic based on their occurring probabilities. Then these genes were used to perform the functional enrichment annotation with R package clusterProfiler together with a FDR correction. Finally the top-ranked n (default 3) GO terms were selected to represent each topic’s key functions.
}
\usage{
Topic_func_anno(model, species = "Hs", topNum = 5, plot = TRUE, plot_path = "./topic_annotation_GO.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  Object of class "LDA" with the optimal topic number.
}
  \item{species}{
  The species the cells belong to. Currently species="Hs" for homo sapiens or species="Mm" for mus musculus are available.
}
  \item{topNum}{
  The top-ranked number (default 5) GO terms were selected to represent each topic’s key functions.
}
  \item{plot}{
  TRUE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
  The save path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\value{
A dataframe showing the annotation result.

}
\references{
Yu G, Wang LG, Han Y, He Q. Y. clusterProfiler: an R package for comparing biological themes among gene clusters. Omics. 2012;16:284–287. doi: 10.1089/omi.2011.0118.
}
\author{
Bin Duan
}
