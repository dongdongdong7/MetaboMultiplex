#' @title Assign plex information
#' @description
#' Assign plex information for each peak containing precursorMz, specialFragment and plexIdx.
#'
#' @param data data list.
#' @param plexPara plex parameter list.
#' @param thread thread.
#'
#' @return data list.
#' @export
#'
#' @examples
#' data <- assign_plexInfo(data, plexPara, thread = 3)
assign_plexInfo <- function(data, plexPara, thread = 1){
  sp2List <- data$peaksInfo$spectra
  loop <- function(i){
    x <- sp2List[[i]]
    if(!is.null(x)){
      precursorMz <- Spectra::precursorMz(x)
      mz <- unlist(Spectra::mz(x))
      int <- unlist(Spectra::intensity(x))
      sf_df <- purrr::list_rbind(lapply(plexPara$specialFrag, function(specialFrag) {
        sf <- specialFrag[sapply(specialFrag, function(Frag) {
          any(dplyr::near(Frag, mz, tol = plexPara$tolMz2))
        })]
        if(length(sf) == 0) return(NULL)
        max_int <- sapply(sf, function(a) {
          idx <- which(dplyr::near(mz, a, tol = plexPara$tolMz2))
          max(int[idx])
        })
        sf <- sf[which.max(max_int)]
        sf_int <- max(max_int)
        plex_idx <- which(dplyr::near(specialFrag, sf, tol = plexPara$tolMz2))
        return(data.frame(sf = sf, sf_int = sf_int, plex_idx = plex_idx))
      }))
      if(nrow(sf_df) == 0) return(data.frame(precursorMz = NA, specialFragment = NA, plexIdx = NA))
      else{
        sf_df <- sf_df[which.max(sf_df$sf_int), ]
        return(data.frame(precursorMz = precursorMz, specialFragment = sf_df$sf, plexIdx = sf_df$plex_idx))
      }
    }else return(data.frame(precursorMz = NA, specialFragment = NA, plexIdx = NA))
  }
  pb <- utils::txtProgressBar(max = length(sp2List), style = 3)
  if(thread == 1){
    plexList <- lapply(1:length(sp2List), function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    plexList <- foreach::`%dopar%`(foreach::foreach(i = 1:length(sp2List),
                                                   .packages = c("Spectra", "dplyr"),
                                                   .options.snow = opts),
                                  {
                                    loop(i)
                                  })
    snow::stopCluster(cl)
    gc()
  }else stop("Thread wrong!")
  plexDf <- dplyr::as_tibble(purrr::list_rbind(plexList))
  if(nrow(data$peaksInfo) == nrow(plexDf)) data$peaksInfo <- dplyr::as_tibble(cbind(data$peaksInfo, plexDf))
  else stop("length(sp2List) != nrow(peaksInfo)")
  return(data)
}
#' @title Assign tagNum
#' @description
#' Assign tagNum using cliqueMS's results.
#'
#' @param data data list.
#' @param plexPara plexPara list.
#' @param thread thread.
#'
#' @return data list.
#' @export
#'
#' @examples
#' data <- assign_tagNum(data, plexPara = plexPara, thread = 3)
assign_tagNum <- function(data, plexPara, thread = 1){
  peaksInfo <- data$peaksInfo
  deltaIsotope1 <- 1.0033
  deltaIsotope2 <- deltaIsotope1 / 2
  deltaIsotope3 <- deltaIsotope1 / 3
  deltaIsotope4 <- deltaIsotope1 / 4

  loop <- function(i){
    peakInfo <- peaksInfo[i, ]
    cliqueGroup <- peaksInfo %>%
      dplyr::filter(sample == peakInfo$sample) %>%
      dplyr::filter(cliqueGroup == peakInfo$cliqueGroup)
    isoIdx <- stringr::str_extract(peakInfo$isotope, "(?<=\\[)\\d+(?=\\])")
    if(!is.na(isoIdx)){
      cliqueGroup <- cliqueGroup[which(stringr::str_extract(cliqueGroup$isotope, "(?<=\\[)\\d+(?=\\])") == isoIdx), ]
      tagNum <- which(dplyr::near(mean(abs(diff(cliqueGroup$mz))),
                                  c(deltaIsotope1, deltaIsotope2, deltaIsotope3, deltaIsotope4),
                                  tol = plexPara$tolMz1))
      if(length(tagNum) == 0) tagNum <- NA
    }else tagNum <- NA
    return(tagNum)
  }

  pb <- utils::txtProgressBar(max = nrow(peaksInfo), style = 3)
  if(thread == 1){
    tagNum_vec <- sapply(1:nrow(peaksInfo), function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    tagNum_vec <- foreach::`%dopar%`(foreach::foreach(i = 1:nrow(peaksInfo),
                                                    .packages = c("stringr", "dplyr"),
                                                    .options.snow = opts,
                                                    .combine = "c"),
                                   {
                                     loop(i)
                                   })
    snow::stopCluster(cl)
    gc()
  }
  peaksInfo$tagNum <- tagNum_vec
  data$peaksInfo <- peaksInfo
  return(data)
}
# assign_tagNum(data, plexPara = plexPara)
# assign_tagNum <- function(data, plexPara, isoRt = 1){
#   browser()
#
#   peaksInfo <- data$peaksInfo
#
#   deltaIsotope1 <- 1.0033
#   deltaIsotope2 <- deltaIsotope1 / 2
#   deltaIsotope3 <- deltaIsotope1 / 3
#   deltaIsotope4 <- deltaIsotope1 / 4
#
#   isoTb <- lapply(1:nrow(peaksInfo), function(i) {
#     print(paste0(i, " / ", nrow(peaksInfo)))
#     sampleIdx <- peaksInfo[i, ]$sample
#     mz_isotope1 <- peaksInfo[i, ]$mz + deltaIsotope1
#     mz_isotope2 <- peaksInfo[i, ]$mz + deltaIsotope2
#     mz_isotope3 <- peaksInfo[i, ]$mz + deltaIsotope3
#     mz_isotope4 <- peaksInfo[i, ]$mz + deltaIsotope4
#     iso_tmp1 <- peaksInfo %>%
#       dplyr::filter(sample == sampleIdx) %>%
#       dplyr::filter(dplyr::near(mz, mz_isotope1, tol = plexPara$tolMz1) & dplyr::near(rt, peaksInfo[i, ]$rt, tol = isoRt)) %>%
#       dplyr::filter(maxo < peaksInfo[i, ]$maxo)
#     iso_tmp2 <- peaksInfo %>%
#       dplyr::filter(sample == sampleIdx) %>%
#       dplyr::filter(dplyr::near(mz, mz_isotope2, tol = plexPara$tolMz1) & dplyr::near(rt, peaksInfo[i, ]$rt, tol = isoRt)) %>%
#       dplyr::filter(maxo < peaksInfo[i, ]$maxo)
#     iso_tmp3 <- peaksInfo %>%
#       dplyr::filter(sample == sampleIdx) %>%
#       dplyr::filter(dplyr::near(mz, mz_isotope3, tol = plexPara$tolMz1) & dplyr::near(rt, peaksInfo[i, ]$rt, tol = isoRt)) %>%
#       dplyr::filter(maxo < peaksInfo[i, ]$maxo)
#     iso_tmp4 <- peaksInfo %>%
#       dplyr::filter(sample == sampleIdx) %>%
#       dplyr::filter(dplyr::near(mz, mz_isotope4, tol = plexPara$tolMz1) & dplyr::near(rt, peaksInfo[i, ]$rt, tol = isoRt)) %>%
#       dplyr::filter(maxo < peaksInfo[i, ]$maxo)
#     iso_tmp <- rbind(iso_tmp1, iso_tmp2, iso_tmp3, iso_tmp4)
#     if(nrow(iso_tmp) != 0){
#       iso_tmp <- rbind(peaksInfo[i, ], iso_tmp)
#     }
#     return(iso_tmp)
#   })
#   nrow_isoTb <- sapply(isoTb, nrow)
# }

#' @title Peak Grouping
#' @description
#' Peak grouping step.
#'
#' @param data data list.
#' @param plexPara plexPara list.
#' @param thread thread.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' peakGrouping(data = data, plexPara = plexPara, thread = 3)
peakGrouping <- function(data, plexPara, thread = 1){
  peaksInfo <- data$peaksInfo
  peakGroupList <- lapply(unique(peaksInfo$sample), function(n) {
    message(paste0(n, "/", length(unique(peaksInfo$sample))))
    peaksInfo_n <- peaksInfo %>%
      dplyr::filter(sample == n)
    peaksInfo_n_plex <- peaksInfo_n %>%
      dplyr::filter(!is.na(plexIdx)) %>%
      dplyr::arrange(mz)
    peaksInfo_n_unplex <- peaksInfo_n %>%
      dplyr::filter(is.na(plexIdx)) %>%
      dplyr::arrange(mz)
    delete_idx <<- c()
    print(length(delete_idx))
    peakGroupListAll <- lapply(1:(plexPara$plexNumber - 1), function(start_loop) {
      peakGroupList <- lapply(1:nrow(peaksInfo_n_plex), function(i) {
        if(i %in% delete_idx) return(NULL)
        if(peaksInfo_n_plex[i, ]$plexIdx == start_loop){
          idx <- sapply(1:(plexPara$plexNumber - start_loop), function(j) {
            idx <- which(dplyr::near(peaksInfo_n_plex$mz - peaksInfo_n_plex[i, ]$mz, plexPara$deltaMz * j, tol = plexPara$tolMz1) &
                           dplyr::near(peaksInfo_n_plex$rt, peaksInfo_n_plex[i, ]$rt, tol = plexPara$deltaRt) &
                           peaksInfo_n_plex$plexIdx - peaksInfo_n_plex[i, ]$plexIdx == j)
            # if length(idx) > 1, ppc will be used to filter peaks.
            if(length(idx) > 1){
              target_cpid <- peaksInfo_n_plex[i, ]$cpid
              other_cpid <- peaksInfo_n_plex[idx, ]$cpid
              target_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = target_cpid)[1]
              other_chrs <- xcms::chromPeakChromatograms(data$rawData, peaks = other_cpid)
              ppc_vec <- sapply(1:nrow(other_chrs), function(k) {
                other_chr <- other_chrs[k]
                .compare_peaks(chr1 = other_chr, target_chr)
              })
              # idx <- idx[which.min(peaksInfo_n_plex$mz[idx] - peaksInfo_n_plex[i, ]$mz)]
              # idx <- idx[which.min(abs(peaksInfo_n_plex$rt[idx] - peaksInfo_n_plex[i, ]$rt))]
              idx <- idx[which.max(ppc_vec)]
            }
            if(length(idx) == 0) return(NA)
            else return(idx)
          })
          idx[idx %in% delete_idx] <- NA
          c(rep(NA, start_loop - 1), i, idx)
        }else return(NULL)
      })
      peakGroupList <- peakGroupList[!sapply(peakGroupList, is.null)]
      delete_idx <<- c(delete_idx, unique(purrr::list_c(peakGroupList)))
      delete_idx <<- delete_idx[!is.na(delete_idx)]
      return(peakGroupList)
    })
    peakGroupListAll <- purrr::list_flatten(peakGroupListAll)
    peakGroupListAll <- lapply(peakGroupListAll, function(x) {peaksInfo_n_plex[x, ]})
    peakGroup_num_vec <- sapply(peakGroupListAll, function(x) {length(which(!is.na(x$cpid)))})
    peakGroupListAll <- peakGroupListAll[which(peakGroup_num_vec >= 3)]
    # Use peaksInfo_n_unplex to fill peakGroup
    peakGroupListAll <- lapply(peakGroupListAll, function(x) {
      missing_plexIdx <- which(is.na(x$plexIdx))
      if(length(missing_plexIdx) == 0) return(x)
      ref_plexIdx <- which(!is.na(x$plexIdx))[1]
      ref_cpid <- x$cpid[ref_plexIdx]
      ref_mz <- x$mz[ref_plexIdx]
      ref_rt <- x$rt[ref_plexIdx]
      fill_idx <- sapply(missing_plexIdx, function(i) {
        if(i > ref_plexIdx) fill_mz <- ref_mz + plexPara$deltaMz * (i - ref_plexIdx)
        else fill_mz <- ref_mz - plexPara$deltaMz * (ref_plexIdx - i)
        idx <- which(dplyr::near(peaksInfo_n_unplex$mz, fill_mz, tol = plexPara$tolMz1) &
                       dplyr::near(peaksInfo_n_unplex$rt, ref_rt, tol = plexPara$deltaRt))
        if(length(idx) > 0){
          fill_cpid <- peaksInfo_n_unplex[idx, ]$cpid
          ref_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = ref_cpid)[1]
          fill_chrs <- xcms::chromPeakChromatograms(data$rawData, peaks = fill_cpid)
          ppc_vec <- sapply(1:nrow(fill_chrs), function(k) {
            fill_chr <- fill_chrs[k]
            .compare_peaks(chr1 = fill_chr, ref_chr)
          })
          idx <- idx[which.max(ppc_vec)]
          ppc <- ppc_vec[which.max(ppc_vec)]
          if(ppc < 0.8) idx <- NA
        }
        else idx <- NA
        return(idx)
      })
      x[missing_plexIdx, ] <- peaksInfo_n_unplex[c(1, fill_idx), ][-1, ]
      return(x)
    })
    #rm(delete_idx)
    print(length(delete_idx))
    return(peakGroupListAll)
  })
  data$peakGroupList <- peakGroupList
  return(data)
}

# Fill blank plex.
# test <- xcms::manualChromPeaks(data$rawData,
#                                chromPeaks = cbind(
#                                  mzmin = c(243.13, 251.18), mzmax = c(243.14, 251.19),
#                                  rtmin = c(740, 740), rtmax = c(750, 750)
#                                ),
#                                sample = 1)
