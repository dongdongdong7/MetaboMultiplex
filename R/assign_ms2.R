# This script is used to store function assigning ms2.
# 241120
# Barry Song

# .get_spectra2(sps = sps, cpid = data$peaksInfo[7000, ]$cpid,
#               mzmin = data$peaksInfo[7000, ]$mzmin, mzmax = data$peaksInfo[7000, ]$mzmax,
#               rtmin = data$peaksInfo[7000, ]$rtmin, rtmax = data$peaksInfo[7000, ]$rtmax,
#               spectraOrigin = "D:\\fudan\\Projects\\2024\\MultichannelR\\Progress\\build_package\\test_data\\AP\\mix1\\CXYC136A240712-dx-DEANS-AP_mix1_0712-1.mzML")
.get_spectra2 <- function(sps, cpid, mzmin, mzmax, rtmin, rtmax, spectraOrigin, mzdiff = 0.001, rtdiff = 0){
  sp2 <- sps %>%
    Spectra::filterDataOrigin(dataOrigin = spectraOrigin) %>%
    Spectra::filterMsLevel(2L) %>%
    Spectra::filterPrecursorMzRange(c(mzmin - mzdiff, mzmin + mzdiff)) %>%
    Spectra::filterRt(c(rtmin - rtdiff, rtmax + rtdiff))

  if(length(sp2) == 0){
    sp2 <- NULL
  }else{
    sp2$cpid <- cpid
  }

  return(sp2)
}
.get_spectra2_closeRT <- function(sps, cpid, rt, mzmin, mzmax, rtmin, rtmax, spectraOrigin,mzdiff = 0.001, rtdiff = 0){
  sp2 <- .get_spectra2(sps = sps, cpid = cpid, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin, rtmax = rtmax, spectraOrigin = spectraOrigin, mzdiff = mzdiff, rtdiff = rtdiff)
  if(is.null(sp2)){
    return(NULL)
  }else{
    idx <- which.min(abs(Spectra::rtime(sp2) - rt))
    sp2 <- sp2[idx]
  }
  return(sp2)
}
# sp2 is a Spectra object whose length equals 1.
.standardizeSpectra <- function(sp2, th = 0.01){
  if(is.null(sp2)){
    return(sp2)
  }

  norm_fun <- function(z, ...) {
    z[, "intensity"] <- z[, "intensity"] /
      max(z[, "intensity"], na.rm = TRUE) * 100
    z
  }

  sp2 <- Spectra::filterIntensity(sp2, intensity = max(Spectra::intensity(sp2)) * th)
  sp2 <- Spectra::addProcessing(sp2, FUN = norm_fun)

  return(sp2)
}

#' @title Assign MS2
#' @description
#' Each peak is assigned a ms2.
#'
#' @param data data list with rawData after peak annotation.
#' @param thread thread.
#' @param mzdiff mzdiff.
#' @param rtdiff rtdiff.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' getSpectra2(data = data)
getSpectra2 <- function(data, thread = 1, mzdiff = 0.001, rtdiff = 0){
  sps <- MsExperiment::spectra(data$rawData)
  spectraOrigin_vec <- MsExperiment::sampleData(data$rawData)$spectraOrigin
  peaksInfo <- data$peaksInfo
  loop <- function(i){
    peakInfo <- peaksInfo[i, ]
    sp2_tmp <- .get_spectra2_closeRT(sps = sps, cpid = peakInfo$cpid, rt = peakInfo$rt,
                                     mzmin = peakInfo$mzmin, mzmax = peakInfo$mzmax,
                                     rtmin = peakInfo$rtmin, rtmax = peakInfo$rtmax,
                                     spectraOrigin = spectraOrigin_vec[peakInfo$sample],
                                     mzdiff = mzdiff, rtdiff = rtdiff)
    .standardizeSpectra(sp2_tmp)
  }
  pb <- utils::txtProgressBar(max = nrow(peaksInfo), style = 3)
  if(thread == 1){
    sp2List <- lapply(1:nrow(peaksInfo), function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    sp2List <- foreach::`%dopar%`(foreach::foreach(i = 1:nrow(peaksInfo),
                                                   .packages = c("Spectra"),
                                                   .options.snow = opts),
                                  {
                                    loop(i)
                                  })
    snow::stopCluster(cl)
    gc()
  }else stop("Thread wrong!")
  if(length(sp2List) == nrow(peaksInfo)) peaksInfo$spectra <- sp2List
  data$peaksInfo <- peaksInfo
  return(data)
}
