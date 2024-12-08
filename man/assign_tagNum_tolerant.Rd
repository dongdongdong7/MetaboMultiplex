% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak_group.R
\name{assign_tagNum_tolerant}
\alias{assign_tagNum_tolerant}
\title{Assign tagNum with a tolerant method}
\usage{
assign_tagNum_tolerant(data, isoMzDiff = 0.01, isoRtDiff = 1, thread = 3)
}
\arguments{
\item{data}{data list.}

\item{isoMzDiff}{The mz difference between a compound and its isotope.}

\item{isoRtDiff}{The rt difference between a compound and its isotope.}

\item{thread}{thread.}
}
\value{
A data list.
}
\description{
This function should be execute optionally after assign_tagNum function.
Assign tagNum in a tolerant way to peaks where tagNum is NA.
}
\examples{
assign_tagNum_tolerant(data, plexPara = plexPara)
}
