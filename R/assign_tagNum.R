#' @title Assign tagNum
#' @description
#' Assign tagNum using cliqueMS isotope result.
#'
#' @param data A data list.
#' @param isoMzDiff The mz difference between a compound and its isotope.
#' @param thread thread.
#'
#' @return A data list.
#'
#' @examples
#' data <- .assign_tagNum_cliqueMS(data, isoMzDiff = 0.01, thread = 3)
.assign_tagNum_cliqueMS <- function(data, isoMzDiff = 0.01, thread = 1){
  peaksInfo <- data$peaksInfo
  if(is.null(data$peaksInfo$tagNum)) idx_no_tagNum <- 1:nrow(peaksInfo)
  else idx_no_tagNum <- which(is.na(data$peaksInfo$tagNum))

  deltaIsotope1 <- 1.0033
  deltaIsotope2 <- deltaIsotope1 / 2
  deltaIsotope3 <- deltaIsotope1 / 3
  deltaIsotope4 <- deltaIsotope1 / 4
  peaksInfo$isoGroup <- stringr::str_extract(peaksInfo$isotope, "(?<=\\[)\\d+(?=\\])")

  loop <- function(i){
    peakInfo <- peaksInfo[i, ]
    isoIdx <- peakInfo$isoGroup
    if(is.na(isoIdx)) return(NA)
    isoGroup <- peaksInfo %>%
      dplyr::filter(sample == peakInfo$sample) %>%
      dplyr::filter(isoGroup == isoIdx)
    if(nrow(isoGroup) >= 2){
      tagNum <- which(dplyr::near(mean(abs(diff(sort(isoGroup$mz)))),
                                  c(deltaIsotope1, deltaIsotope2, deltaIsotope3, deltaIsotope4),
                                  tol = isoMzDiff))
      if(length(tagNum) == 0) tagNum <- NA
    }else tagNum <- NA
    return(tagNum)
  }

  pb <- utils::txtProgressBar(max = length(idx_no_tagNum), style = 3)
  if(thread == 1){
    tagNum_vec <- sapply(1:length(idx_no_tagNum), function(x) {
      utils::setTxtProgressBar(pb, x)
      i <- idx_no_tagNum[x]
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    tagNum_vec <- foreach::`%dopar%`(foreach::foreach(i = 1:nrow(peaksInfo),
                                                      .packages = c("dplyr"),
                                                      .options.snow = opts,
                                                      .combine = "c"),
                                     {
                                       i <- idx_no_tagNum[x]
                                       loop(i)
                                     })
    snow::stopCluster(cl)
    gc()
  }
  if(is.null(data$peaksInfo$tagNum)) peaksInfo$tagNum <- tagNum_vec
  else peaksInfo$tagNum[idx_no_tagNum] <- tagNum_vec
  data$peaksInfo <- peaksInfo
  return(data)
}

#' @title Assign tagNum
#' @description
#' Assign tagNum using MetaboMultiplex method.
#'
#' @param data A data list.
#' @param isoMzDiff The mz difference between a compound and its isotope.
#' @param isoRtDiff The rt difference between a compound and its isotope.
#' @param thread thread.
#'
#' @return A data list.
#'
#' @examples
#' data <- .assign_tagNum_MetaboMultiplex(data, thread = 3)
.assign_tagNum_MetaboMultiplex <- function(data, isoMzDiff = 0.01, isoRtDiff = 1,thread = 3){
  # TODO: need chunk parameter
  peaksInfo <- data$peaksInfo
  if(is.null(data$peaksInfo$tagNum)) idx_no_tagNum <- 1:nrow(peaksInfo)
  else idx_no_tagNum <- which(is.na(data$peaksInfo$tagNum))

  deltaIsotope1 <- 1.0033
  deltaIsotope2 <- deltaIsotope1 / 2
  deltaIsotope3 <- deltaIsotope1 / 3
  deltaIsotope4 <- deltaIsotope1 / 4

  loop <- function(i){
    sampleIdx <- peaksInfo[i, ]$sample
    mz_isotope1 <- peaksInfo[i, ]$mz + deltaIsotope1
    mz_isotope2 <- peaksInfo[i, ]$mz + deltaIsotope2
    mz_isotope3 <- peaksInfo[i, ]$mz + deltaIsotope3
    mz_isotope4 <- peaksInfo[i, ]$mz + deltaIsotope4

    idx <- sapply(c(mz_isotope1, mz_isotope2, mz_isotope3, mz_isotope4), function(mz_isotope) {
      idx_tmp <- which(peaksInfo$sample == sampleIdx &
                         dplyr::near(peaksInfo$mz, mz_isotope, tol = isoMzDiff) &
                         dplyr::near(peaksInfo$rt, peaksInfo[i, ]$rt, tol = isoRtDiff) &
                         peaksInfo$maxo < peaksInfo[i, ]$maxo)
      if(length(idx_tmp) == 0) idx_tmp <- NA
      else if(length(idx_tmp) > 1){
        ref_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[i, ]$cpid)[1]
        ppc_vec <- sapply(idx_tmp, function(j) {
          if(is.na(j)) return(NA)
          tmp_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[j, ]$cpid)[1]
          #.compare_peaks(tmp_chr, ref_chr)
          # Retention time shift of natural isotope 13C is weak
          # So apexAlign is FALSE, and method is pearson.
          .comparePeaks(chr1 = tmp_chr, chr2 = ref_chr, method = "pearson", apexAlign = FALSE, plot = FALSE)
        })
        idx_tmp <- idx_tmp[which.max(ppc_vec)]
      }
      return(idx_tmp)
    })
    if(length(idx[!is.na(idx)]) > 1){
      ref_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[i, ]$cpid)[1]
      ppc_vec <- sapply(idx, function(j) {
        if(is.na(j)) return(NA)
        tmp_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[j, ]$cpid)[1]
        #.compare_peaks(tmp_chr, ref_chr)
        # Retention time shift of natural isotope 13C is weak
        # So apexAlign is FALSE, and method is pearson.
        .comparePeaks(chr1 = tmp_chr, chr2 = ref_chr, method = "pearson", apexAlign = FALSE, plot = FALSE)
      })
      if(all(is.na(ppc_vec))) return(NA) # ppc_vec all NA
      else{
        tmp <- 1:4
        idx[tmp[-which.max(ppc_vec) ]] <- NA
      }
    }
    else if(all(is.na(idx))) return(NA)
    return(which(!is.na(idx)))
  }

  pb <- utils::txtProgressBar(max = length(idx_no_tagNum), style = 3)
  if(thread == 1){
    tagNum_vec <- sapply(1:length(idx_no_tagNum), function(x) {
      utils::setTxtProgressBar(pb, x)
      #print(paste0(x, "/", length(idx_no_tagNum)))
      i <- idx_no_tagNum[x]
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    tagNum_vec <- foreach::`%dopar%`(foreach::foreach(x = 1:length(idx_no_tagNum),
                                                                .packages = c("xcms", "dplyr"),
                                                                .combine = "c",
                                                                .export = c(".comparePeaks"),
                                                                .options.snow = opts),
                                               {
                                                 i <- idx_no_tagNum[x]
                                                 loop(i)
                                               })
    snow::stopCluster(cl)
    gc()
  }else stop("Thread wrong!")
  if(is.null(data$peaksInfo$tagNum)) peaksInfo$tagNum <- tagNum_vec
  else peaksInfo$tagNum[idx_no_tagNum] <- tagNum_vec
  data$peaksInfo <- peaksInfo
  return(data)
}

#' @title Assign tagNum for peaks
#' @description
#' Assign tagNum using MetaboMultiplex method or cliqueMS method.
#'
#' @param data A data list.
#' @param isoMzDiff The mz difference between a compound and its isotope.
#' @param isoRtDiff The rt difference between a compound and its isotope.
#' @param method MetaboMultiplex or cliqueMS.
#' @param thread thread.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' data <- tagNumAssigning(data)
tagNumAssigning <- function(data, isoMzDiff = 0.01, isoRtDiff = 1,
                            method = c("MetaboMultiplex", "cliqueMS", "cliqueMS+MetaboMultiplex", "MetaboMultiplex+cliqueMS")[1],
                            thread = 3){
  if(method == "MetaboMultiplex"){
    message("You are using MetaboMultiplex to assign tagNum...")
    data <- .assign_tagNum_MetaboMultiplex(data, isoMzDiff = isoMzDiff, isoRtDiff = isoRtDiff, thread = thread)
  }else if(method == "cliqueMS"){
    message("You are using cliqueMS to assign tagNum...")
    data <- .assign_tagNum_cliqueMS(data, isoMzDiff = isoMzDiff, thread = thread)
  }else if(method == "cliqueMS+MetaboMultiplex"){
    message("You are using cliqueMS to assign tagNum...")
    data <- .assign_tagNum_cliqueMS(data, isoMzDiff = isoMzDiff, thread = thread)
    message("You are using MetaboMultiplex to assign tagNum...")
    data <- .assign_tagNum_MetaboMultiplex(data, isoMzDiff = isoMzDiff, isoRtDiff = isoRtDiff, thread = thread)
  }else if(method == "MetaboMultiplex+cliqueMS"){
    message("You are using MetaboMultiplex to assign tagNum...")
    data <- .assign_tagNum_MetaboMultiplex(data, isoMzDiff = isoMzDiff, isoRtDiff = isoRtDiff, thread = thread)
    message("You are using cliqueMS to assign tagNum...")
    data <- .assign_tagNum_cliqueMS(data, isoMzDiff = isoMzDiff, thread = thread)
  }else stop("Method is wrong!")
  return(data)
}
