# Initialize data from loading data to data annotation
# 241120
# Barry Song

#' @title Load data
#' @description
#' Load data.
#'
#' @param dataPara dataPara.
#'
#' @return Data list.
#' @export
#'
#' @examples
#' dataPara <- set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/AP/mix1/",
#'                          res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
#'                          sampleData = data.frame(sample_id = c("mix1_1", "mix1_2", "mix1_3"), injection_index = 1:3))
#' data <- load_data(dataPara = dataPara)
load_data <- function(dataPara){
  rawData <- MsExperiment::readMsExperiment(spectraFiles = dataPara$sampleInfo$sample_path,
                                            sampleData = dataPara$sampleInfo)
  return(list(rawData = rawData))
}

#' @title Peak picking
#' @description
#' Execute the peak picking algorithm.
#'
#' @param data data.
#' @param xcmsPara xcmsPara.
#' @param chunkSize chunkSize.
#' @param BPPARAM BPPARAM.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' data <- peakPicking(data, xcmsPara = set_xcmsPara())
peakPicking <- function(data, xcmsPara, chunkSize = 3L, BPPARAM = BiocParallel::SnowParam(workers = 3)){
  rawData_new <- xcms::findChromPeaks(object = data$rawData, param = xcmsPara$centwavePara, msLevel = 1L, chunkSize = chunkSize, BPPARAM = BPPARAM)
  rawData_new <- xcms::refineChromPeaks(object = rawData_new, param = xcmsPara$mergepeakPara, msLevel = 1L, chunkSize = chunkSize, BPPARAM = BPPARAM)
  peaksInfo <- dplyr::as_tibble(cbind(xcms::chromPeaks(rawData_new), xcms::chromPeakData(rawData_new)), rownames = "cpid")
  return(list(rawData = rawData_new,
              peaksInfo = peaksInfo))
}
#' @title Peak Annotation
#' @description
#' Peak Annotation using cliqueMS.
#'
#' @param data data list with rawData after peak picking.
#' @param polarity positive or negative.
#' @param adinfo adduct information.
#' @param ppm Relative error in ppm to consider that two features have the mass difference of an isotope.
#' @param chunkSize Sample size per thread.
#' @param thread Parallel thread.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' dataPara <- set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/AP/mix1/",
#' res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
#' sampleData = data.frame(sample_id = c("mix1_1", "mix1_2", "mix1_3"), injection_index = 1:3))
#' data <- load_data(dataPara = dataPara)
#' data <- peakPicking(data, xcmsPara = set_xcmsPara())
#' data(positive.adinfo, package = "cliqueMS")
#' positive.adinfo <- positive.adinfo[positive.adinfo$adduct %in%
#'                                      c("[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+H-NH3]+", "[M+K]+", "[M+NH4]+"),]
#' data <- peakAnnotation(data, polarity = "positive", adinfo = positive.adinfo, chunkSize = 1, thread = 1)
peakAnnotation <- function(data, polarity = "positive", adinfo, ppm = 10, chunkSize = 1, thread = 1){
  chunks <- split(data$rawData, cut(seq_along(data$rawData), length(data$rawData), labels = FALSE))
  chunks <- split(chunks, rep(1:(length(chunks) %/% chunkSize + 1), each = chunkSize, length.out = length(chunks)))
  for(i in 1:length(chunks)){
    for(j in 1:length(chunks[[i]])){
      attributes(chunks[[i]][[j]])$order <- i
    }
  }
  loop <- function(data_n){
    data_n <- xcms:::.XCMSnExp2xcmsSet(data_n)
    set.seed(2)
    ex.cliqueGroups <- cliqueMS::getCliques(data_n, filter = TRUE, silent = FALSE)
    ex.Isotopes <- cliqueMS::getIsotopes(ex.cliqueGroups, ppm = ppm, maxCharge = 4)
    ex.Adducts <- cliqueMS::getAnnotation(ex.Isotopes, ppm = ppm,
                                          adinfo = adinfo, polarity = polarity,
                                          normalizeScore = TRUE)
    resTable <- ex.Adducts@peaklist[, c("mz", "rt", "maxo", "sample", "cliqueGroup", "isotope", "mass1", "an1")]
    return(resTable)
  }
  # pb <- progress::progress_bar$new(
  #   format = paste0("[:bar] :current/:",
  #                   "total (:percent) in ",
  #                   ":elapsed"),
  #   total = length(data$rawData),
  #   clear = FALSE
  # )
  # opts <- list(progress = function(n) {pb$tick()})
  pb <- utils::txtProgressBar(max = length(chunks), style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl) #doSNOW:::doSNOW()
  resList <- foreach::`%dopar%`(foreach::foreach(chunk = chunks, n = 1:length(chunks),
                                                 .packages = c("cliqueMS", "xcms"),
                                                 #.verbose = TRUE,
                                                 .options.snow = opts),
                                {
                                  resTable <- purrr::list_rbind(lapply(chunk, loop))
                                  return(resTable)
                                })
  snow::stopCluster(cl)
  gc()
  resTable <- purrr::list_rbind(resList)
  resTable <- dplyr::as_tibble(resTable, rownames = "cpid") %>%
    dplyr::select(cpid, cliqueGroup, isotope, mass1, an1)
  resTable[which(resTable$an1 == ""),]$an1 <- NA
  data$peaksInfo <- dplyr::left_join(data$peaksInfo, resTable, by = c("cpid" = "cpid"))
  return(data)
}

# Visualize ChromPeaks using  xcms.
# It will return a XChromatograms or MChromatograms object.
# xchr_test <- xcms::chromPeakChromatograms(data$rawData, peaks = data$peaksInfo$cpid[7000:9000])
# xcms::plot(xchr_test[1])
# xcms::plot(xchr_test[1000])
