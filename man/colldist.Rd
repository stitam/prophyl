% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colldist.R
\name{colldist}
\alias{colldist}
\title{Distances between sample collection dates}
\usage{
colldist(tree_tbl, date_var = "collection_day")
}
\arguments{
\item{tree_tbl}{tibble; the phylogenetic tree in tibble format.}

\item{date_var}{character; the name of the variable which stores sample
collection dates.}
}
\value{
a tibble with three columns.
}
\description{
This functions calculates distances between sample collection dates for all
pairs of tips within a phylogenetic tree.
}
\examples{
\dontrun{
# TODO put a tree_tbl into data so it can be imported.
colldist(tree_tbl)
}
}
