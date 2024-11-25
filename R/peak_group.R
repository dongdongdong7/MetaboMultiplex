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
  peaksInfo$isoGroup <- stringr::str_extract(peaksInfo$isotope, "(?<=\\[)\\d+(?=\\])")

  loop <- function(i){
    peakInfo <- peaksInfo[i, ]
    isoIdx <- peakInfo$isoGroup
    if(is.na(isoIdx)) return(NA)
    isoGroup <- peaksInfo %>%
      dplyr::filter(sample == peakInfo$sample) %>%
      dplyr::filter(isoGroup == isoIdx)
    if(nrow(isoGroup) >= 2){
      tagNum <- which(dplyr::near(mean(abs(diff(isoGroup$mz))),
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
                                                    .packages = c("dplyr"),
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
#' @title Assign tagNum with a tolerant method
#' @description
#' This function should be execute optionally after assign_tagNum function.
#' Assign tagNum in a tolerant way to peaks where tagNum is NA.
#'
#' @param data data list.
#' @param plexPara plexPara list.
#' @param isoRt Isotope retention time search range.
#' @param thread thread.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' assign_tagNum_tolerant(data, plexPara = plexPara)
assign_tagNum_tolerant <- function(data, plexPara, isoRt = 1, thread = 3){
  # TODO: need chunk parameter
  peaksInfo <- data$peaksInfo
  idx_no_tagNum <- which(is.na(data$peaksInfo$tagNum))

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
                      dplyr::near(peaksInfo$mz, mz_isotope, tol = plexPara$tolMz1) &
                      dplyr::near(peaksInfo$rt, peaksInfo[i, ]$rt, tol = isoRt) &
                      peaksInfo$maxo < peaksInfo[i, ]$maxo)
      if(length(idx_tmp) == 0) idx_tmp <- NA
      return(idx_tmp)
    })
    if(length(idx[!is.na(idx)]) > 1){
      ref_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[i, ]$cpid)[1]
      ppc_vec <- sapply(idx, function(j) {
        if(is.na(j)) return(NA)
        tmp_chr <- xcms::chromPeakChromatograms(data$rawData, peaks = peaksInfo[j, ]$cpid)[1]
        .compare_peaks(tmp_chr, ref_chr)
      })
      tmp <- 1:4
      idx[tmp[-which.max(ppc_vec) ]] <- NA
    }
    else if(all(is.na(idx))) return(NA)
    return(which(!is.na(idx)))
  }

  pb <- utils::txtProgressBar(max = length(idx_no_tagNum), style = 3)
  if(thread == 1){
    tagNum_vec_no_tagNum <- sapply(1:length(idx_no_tagNum), function(x) {
      utils::setTxtProgressBar(pb, x)
      i <- idx_no_tagNum[x]
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    tagNum_vec_no_tagNum <- foreach::`%dopar%`(foreach::foreach(x = 1:length(idx_no_tagNum),
                                                    .packages = c("xcms", "dplyr"),
                                                    .combine = "c",
                                                    .export = c(".compare_peaks"),
                                                    .options.snow = opts),
                                   {
                                     i <- idx_no_tagNum[x]
                                     loop(i)
                                   })
    snow::stopCluster(cl)
    gc()
  }else stop("Thread wrong!")
  peaksInfo$tagNum[idx_no_tagNum] <- tagNum_vec_no_tagNum
  data$peaksInfo <- peaksInfo
  return(data)
}

#' @title Peak Grouping
#' @description
#' Peak grouping step.
#'
#' @param data data list.
#' @param plexPara plexPara list.
#' @param thread thread.
#' @param extra_formula Additional molecular formula from derivatisation.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' peakGrouping(data = data, plexPara = plexPara, thread = 3)
peakGrouping <- function(data, plexPara, thread = 1, extra_formula = "C14H15NO2S"){
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
  peakGroupList <- purrr::list_flatten(peakGroupList)

  reagent <- MetaboCoreUtils::standardizeFormula(extra_formula)
  reagent_mz <- MetaboCoreUtils::formula2mz(reagent, adduct = "[M+H]+")[[1]]
  reagent_mass <- MetaboCoreUtils::mz2mass(reagent_mz)[[1]]

  peakGroup <- purrr::list_rbind(lapply(1:length(peakGroupList), function(i) {
    peakGroup <- peakGroupList[[i]]
    table_tmp <- table(peakGroup$tagNum)
    if(length(table_tmp) == 0){
      tagNum <- 1;tagNum_real <- NA
    }else{
      tagNum <- as.integer(names(table_tmp[which.max(table_tmp)]))
      tagNum_real <- tagNum
    }
    ref_plexIdx <- which(!is.na(peakGroup$mz))[1]
    ref_mz <- peakGroup[ref_plexIdx, ]$mz
    ref_rt <- peakGroup[ref_plexIdx, ]$rt
    peaksNum <- length(which(!is.na(peakGroup$mz)))
    sample <- unique(peakGroup$sample[!is.na(peakGroup$sample)])
    reagent_mz_i <- reagent_mz + (ref_plexIdx - 1) * plexPara$deltaMz
    mass <- (ref_mz - reagent_mz_i) * tagNum
    adduct <- paste0(unique(peakGroup$an1[!is.na(peakGroup$an1)]), collapse = ";")
    if(mass < 0) mass <- NA
    peakGroup_new <- dplyr::tibble(mass = mass, rt = ref_rt, peaksNum = peaksNum, sample = sample, tagNum = tagNum_real, adduct = adduct, peaks = list(peakGroup))
    return(peakGroup_new)
  }))
  peakGroup$pgid <- paste0("pg", formatC(1:length(peakGroupList), flag = "0", width = ceiling(log10(length(peakGroupList)))))
  peakGroup <- peakGroup %>%
    dplyr::select(pgid, mass, rt, peaksNum, sample, tagNum, adduct, peaks)
  data$peakGroup <- peakGroup
  return(data)
}

# Fill blank plex.
# test <- xcms::manualChromPeaks(data$rawData,
#                                chromPeaks = cbind(
#                                  mzmin = c(243.13, 251.18), mzmax = c(243.14, 251.19),
#                                  rtmin = c(740, 740), rtmax = c(750, 750)
#                                ),
#                                sample = 1)
