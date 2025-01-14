# This is the script to create parameter list for MetaboMultiplex
# 241118
# Barry Song

#' @title Set plexPara for program.
#' @description
#' Derivated metabolites produce specific secondary mass spectra which can be used to extract and identify.
#'
#' @param targetGroup Target compound group type.
#' @param specialFrag Special fragments from dervited metabolites.
#' @param deltaMz Mass difference from one peak group.
#' @param deltaRt Retention time difference range from one peak group.
#' @param tolMz1 Mass difference tolerance from same fragment on ms1.
#' @param tolMz2 Mass difference tolerance from same fragment on ms2.
#' @param pps Peak to peak similarity threshold. Only peaks greater than pps are retained during the missing value imputation.
#'
#' @return A parameter list.
#' @export
#'
#' @examples
#' plexPara <- set_plexPara(targetGroup = "Amine")
set_plexPara <- function(targetGroup = c("Amine", "Phenol", "Alcohol", "Custom")[1],
                         specialFrag = NULL,
                         deltaMz = 2.0125,
                         deltaRt = 30,
                         tolMz1 = 0.01,
                         tolMz2 = 0.01,
                         pps = 0.8){
  if(targetGroup == "Amine"){
    # specialFrag <- list(
    #   c(171.1038, 172.1101, 173.1165, 174.1229, 175.1291, 176.1351),
    #   c(198.1276, 200.1398, 202.1527, 204.1652, 206.1776, 208.1899),
    #   c(381.1267, 383.1393, 385.1581, 387.1644, 389.1770, 391.1895)
    # )
    specialFrag <- list(
      c(198.1276, 200.1398, 202.1527, 204.1652, 206.1776, 208.1899)
    )
  }else if(targetGroup == "Phenol"){
    specialFrag <- list(
      c(199.1351, 201.1482, 203.1604, 205.1734, 207.1857, 209.1986)
    )
  }else if(targetGroup == "Alcohol"){
    specialFrag <- list(
      c(280.1000, 282.1136, 284.1247, 286.1369, 288.1494, 290.1631)
    )
  }else if(targetGroup == "Custom"){
    if(is.null(specialFrag)) stop("Please input specialFrag!")
    if(!is.list(specialFrag)) stop("Please input a list contains numeric vector!")
  }
  if(!all(diff(sapply(specialFrag, length)) == 0)){
    stop("Plex of specialFrag should be fixed!")
  }
  plexNumber <- sapply(specialFrag, length)[1]
  plexPara <- list(
    specialFrag = specialFrag,
    plexNumber = plexNumber,
    deltaMz = deltaMz,
    deltaRt = deltaRt,
    tolMz1 = tolMz1,
    tolMz2 = tolMz2,
    pps = pps
  )
  return(plexPara)
}

#' @title Set xcmsPara
#' @description
#' See xcms::CentWaveParam function help.
#' See xcms::MergeNeighboringPeakParam function help.
#'
#' @param mergePeak Whether to merge neighbouring peaks.
#'
#' @return A xcmsPara list.
#' @export
#'
#' @examples
#' xcmsPara <- set_xcmsPara()
set_xcmsPara <- function(ppm = 30, peakwidth = c(4, 30), snthresh = 1,
                             noise = 100, prefilter = c(3, 100), firstBaselineCheck = FALSE,
                         mergePeak = TRUE,
                         expandRt = 2, expandMz = 0.01, minProp = 0.75){
  centwavePara <- xcms::CentWaveParam(ppm = ppm, peakwidth = peakwidth, snthresh = snthresh,
                                      noise = noise, prefilter = prefilter, firstBaselineCheck = firstBaselineCheck)
  mergepeakPara <- xcms::MergeNeighboringPeaksParam(expandRt = expandRt, expandMz = expandMz, ppm = ppm, minProp = minProp)
  return(list(centwavePara = centwavePara, mergePeak = mergePeak, mergepeakPara = mergepeakPara))
}

#' @title Set data parameters
#' @description
#' Set data parameters.
#'
#' @param data_dir Sample directory.
#' @param res_dir Results directory.
#' @param sampleData Sample information data frame.
#'
#' @return A parameter list.
#' @export
#'
#' @examples
#' set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/AP/mix1/",
#'              res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
#'              sampleData = data.frame(sample_id = c("mix1_1", "mix1_2", "mix1_3"), injection_index = 1:3))
set_dataPara <- function(data_dir, res_dir, sampleData){
  files_path <- list.files(data_dir, pattern = ".mzML")
  files_path <- paste0(data_dir, files_path)
  sampleData$sample_path <- files_path
  sampleData <- sampleData[order(sampleData$injection_index), ]

  dataPara <- list(sampleInfo = sampleData,
                   sampleNumber = nrow(sampleData),
                   resultDir = res_dir)
  return(dataPara)
}
