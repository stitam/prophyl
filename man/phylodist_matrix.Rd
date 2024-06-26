% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylodist_matrix.R
\name{phylodist_matrix}
\alias{phylodist_matrix}
\title{Temporal (cophenetic) distance matrix}
\usage{
phylodist_matrix(tree, df, id_var, focus_on = NULL, focus_by = NULL)
}
\arguments{
\item{tree}{treedater; a dated phylogenetic tree}

\item{df}{data.frame; a data frame containing sample metadata}

\item{id_var}{character; a variable of the df we want to use as row and
column names. Must be unique for each row of the data frame. Must match
tree tip labels.}

\item{focus_on}{character; one or more values of focus_by to focus on}

\item{focus_by}{character; a variable of the df to focus on}
}
\value{
a matrix
}
\description{
This function calculates a matrix which tells the cophenetic distance between
pairs of samples. This matrix is one of the "distance" matrices used in risk
calculations.
}
\examples{
df <- data.frame(
  assembly = LETTERS[1:5],
  continent = c("europe", "europe", "asia", "europe", "asia"),
  country = c("hungary", "hungary", "china", "serbia", "china"),
  mlst = c("ST1", "ST2", "ST1", "ST2", "ST3")
)
set.seed(0)
tr <- ape::rtree(6, tip.label = df$assembly)
phylodist_matrix(tree = tr, df = df,id_var = "assembly")
}
