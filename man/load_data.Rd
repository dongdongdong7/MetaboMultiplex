% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialize_data.R
\name{load_data}
\alias{load_data}
\title{Load data}
\usage{
load_data(dataPara)
}
\arguments{
\item{dataPara}{dataPara.}
}
\value{
Data list.
}
\description{
Load data.
}
\examples{
dataPara <- set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/AP/mix1/",
                         res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
                         sampleData = data.frame(sample_id = c("mix1_1", "mix1_2", "mix1_3"), injection_index = 1:3))
data <- load_data(dataPara = dataPara)
}
