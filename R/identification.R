#' @title Metabolite Identification
#' @description
#' Metabolite identification using retention time and ms2.
#'
#' @param data data list.
#' @param cmpLibrary cmpLibrary.
#' @param ms2Library ms2Library.
#' @param rt_weight rt_weight.
#' @param ms2_weight ms2_weight.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' data <- fgIdentification(data = data, cmpLibrary = cmpLibrary, ms2Library = ms2Library)
fgIdentification <- function(data, cmpLibrary, ms2Library, rt_weight = 0.7, ms2_weight = 0.3){
  if(!"featureGroup" %in% names(data)) stop("You should get featureGroup first.")
  featureGrop <- data$featureGroup
  message("Metabolite identification using retention time and ms2...")
  pb <- utils::txtProgressBar(max = nrow(featureGrop), style = 3)
  idenResList <- lapply(1:nrow(featureGrop), function(i) {
    utils::setTxtProgressBar(pb, i)
    featureGroup_i <- featureGrop[i, ]
    candidatesIdx <- which(dplyr::near(cmpLibrary$monisotopic_molecular_weight, featureGroup_i$mass, tol = 0.03))
    if(length(candidatesIdx) != 0){
      candidates <- cmpLibrary[candidatesIdx, ]
      # Retention Time Score
      featureGroup_i_rt <- featureGroup_i$rt
      candidates_rt <- candidates$RTP * 60
      candidates_rt_score <-  sapply(candidates_rt, function(x) {
        max(0, 1 - abs(x - featureGroup_i_rt) / 600)
      })

      # Spectra Score
      featureGroup_i_sp <- featureGroup_i$spectra[[1]]
      if(is.null(featureGroup_i_sp)) candidates_ms2_score <- rep(0, length(candidatesIdx))
      else{
        featureGroup_i_sp <- Spectra::peaksData(featureGroup_i_sp)[[1]]
        featureGroup_i_sp[, "mz"] <- round(featureGroup_i_sp[, "mz"])
        candidates_accession <- candidates$accession
        ms2Library_idx <- sapply(candidates_accession, function(acc) {
          idx <- which(ms2Library$accession == acc)
          if(length(idx) == 0) return(NA)
          idx
        })
        candidates_ms2 <- ms2Library[c(ms2Library_idx), ]
        candidates_ms2_score <- sapply(1:nrow(candidates_ms2), function(j) {
          if(is.null(candidates_ms2[j, ]$mz) | is.null(candidates_ms2[j, ]$intensity)) return(0)
          sp_library <- Spectra::Spectra(object = candidates_ms2[j, ])
          sp_library <- .standardizeSpectra(sp_library)
          sp_library <- Spectra::peaksData(sp_library)[[1]]
          tmpList <- Spectra::joinPeaks(x = featureGroup_i_sp, y = sp_library, type = "outer", tolerance = 0.02)
          MsCoreUtils::ndotproduct(x = tmpList$x, y = tmpList$y)
        })
      }
      score <- sapply(1:length(candidates_rt_score), function(j) {
        rt_weight * candidates_rt_score[j] + ms2_weight * candidates_ms2_score[j]
      })
      candidates <- candidates[order(score, decreasing = TRUE), ]
      candidates_rt_score <- round(candidates_rt_score[order(score, decreasing = TRUE)], digits = 4)
      candidates_ms2_score <- round(candidates_ms2_score[order(score, decreasing = TRUE)], digits = 4)
      score <- round(score[order(score, decreasing = TRUE)], digits = 4)
      data.frame(accession = candidates$accession, score = score, score_rt = candidates_rt_score, score_ms2 = candidates_ms2_score)
    }
    else return(NULL)
  })
  featureGroup_new <- purrr::list_rbind(lapply(1:nrow(featureGrop), function(i) {
    utils::setTxtProgressBar(pb, i)
    featureGrop_i <- featureGrop[i, ]
    idenRes <- idenResList[[i]]
    if(!is.null(idenRes)){
      accession_vec <- paste0(idenRes$accession, collapse = ";")
      score_vec <- paste0(idenRes$score, collapse = ";")
      score_rt_vec <- paste0(idenRes$score_rt, collapse = ";")
      score_ms2_vec <- paste0(idenRes$score_ms2, collapse = ";")
    }else{
      accession_vec <- ""
      score_vec <- ""
      score_rt_vec <- ""
      score_ms2_vec <- ""
    }
    featureGrop_i$accession <- accession_vec
    featureGrop_i$score <- score_vec
    featureGrop_i$score_rt <- score_rt_vec
    featureGrop_i$score_ms2 <- score_ms2_vec
    return(featureGrop_i)
  }))
  data$featureGroup <- featureGroup_new
  return(data)
}
