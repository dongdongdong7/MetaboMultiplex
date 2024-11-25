#' @title Peak Aligning
#'
#' @param data data list.
#' @param plexPara plexPara list.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' data <- peakAligning(data = data, plexPara = plexPara)
peakAligning <- function(data, plexPara){
  browser()
  peakGroup <- data$peakGroup %>%
    dplyr::filter(!is.na(mass))
    #dplyr::filter(stringr::str_detect(adduct, pattern = "\\[M\\+H\\]\\+") | is.na(adduct))
  delete_idx <<- c()
  sampleIdx_all <- unique(peakGroup$sample)
  featureGroupIdxList <- lapply(1:nrow(peakGroup), function(i) {
    if(i %in% delete_idx) return(NULL)
    sampleIdx_i <- peakGroup[i, ]$sample
    sampleIdx_o <- sampleIdx_all[sampleIdx_all != sampleIdx_i]
    idx <- sapply(sampleIdx_o, function(j) {
      idx_j <- which(dplyr::near(peakGroup$mass, peakGroup[i, ]$mass, tol = plexPara$tolMz1) &
                       dplyr::near(peakGroup$rt, peakGroup[i, ]$rt, tol = plexPara$deltaRt) &
                       peakGroup$sample == j)
      idx_j <- idx_j[which.min(abs(peakGroup[i, ]$rt - peakGroup[idx_j, ]$rt))]
      if(length(idx_j) == 0) return(NA)
      if(idx_j %in% delete_idx) return(NA)
      return(idx_j)
    })
    idx_ini <- rep(NA, length(sampleIdx_all))
    idx_ini[sampleIdx_i] <- i;idx_ini[sampleIdx_o] <- idx
    delete_idx <<- c(delete_idx, idx_ini[!is.na(idx_ini)])
    return(idx_ini)
  })
  featureGroupIdxList <- featureGroupIdxList[!sapply(featureGroupIdxList, is.null)]
  featureGroup <- purrr::list_rbind(lapply(featureGroupIdxList, function(i){
    pgid_vec <- peakGroup[i, ]$pgid
    mass_i <- mean(peakGroup[i, ]$mass, na.rm = TRUE)
    rt_i <- mean(peakGroup[i, ]$rt, na.rm = TRUE)
    peaksNum_i <- sum(peakGroup[i, ]$peaksNum, na.rm = TRUE)
    table_tagNum <- table(peakGroup[i, ]$tagNum)
    if(length(table_tagNum) == 0) tagNum_i <- NA
    else tagNum_i <- as.integer(names(table_tagNum[which.max(table_tagNum)]))
    adduct_i <- peakGroup[i, ]$adduct
    adduct_i <- list(unique(purrr::list_c(adduct_i)))
    featureGroup <- dplyr::tibble(mass = mass_i, rt = rt_i, pgid = list(pgid_vec),peaksNum = peaksNum_i, tagNum = tagNum_i, adduct = adduct_i)
    return(featureGroup)
  }))
  data$featureGroup <- featureGroup
  return(data)
}
