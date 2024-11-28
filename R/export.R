#' @title Get Quant Result
#'
#' @param data data list.
#' @param plexPara plexPara.
#' @param quant Quant type.
#'
#' @return A data list.
#' @export
#'
#' @examples
#' data <- getQuantRes(data, plexPara = plexPara)
getQuantRes <- function(data, plexPara, quant = c("into", "intb", "maxo")[1]){
  featureGroup <- data$featureGroup
  peakGroup <- data$peakGroup
  featureGroup_new <- purrr::list_rbind(lapply(1:nrow(featureGroup), function(i) {
    featureGroup_i <- featureGroup[i, ]
    pgid_i <- featureGroup_i$pgid[[1]]
    sampleIdx_exists <- which(!is.na(pgid_i))
    sampleIdx_missing <- which(is.na(pgid_i))
    peakGroup_i <- peakGroup[sapply(pgid_i[sampleIdx_exists], function(pgid) {which(peakGroup$pgid == pgid)}), ]
    quant_vec <- rep(NA, length(pgid_i) * plexPara$plexNumber)
    idx_exists <- purrr::list_c(lapply(sampleIdx_exists, function(j) {
      (((j - 1) * plexPara$plexNumber) + 1):((((j - 1) * plexPara$plexNumber) + 1) + plexPara$plexNumber - 1)
    }))
    quant_vec[idx_exists] <- unlist(purrr::list_rbind(peakGroup_i$peaks)[, quant])
    sampleIdx_plexIdx <- unlist(purrr::map(1:length(pgid_i), function(x){
      paste0(x, "_", 1:plexPara$plexNumber)
    }))
    report <- as.data.frame(matrix(quant_vec, nrow = 1))
    colnames(report) <- sampleIdx_plexIdx
    featureGroup_i <- dplyr::as_tibble(cbind(featureGroup_i, report))
  }))
  data$featureGroup <- featureGroup_new
  return(data)
}

#' @title Export Restable
#'
#' @param data data list.
#' @param file file path.
#'
#' @export
#'
#' @examples
#' exportRes(data = data, file = "D:/fudan/resTable.xlsx")
exportRes <- function(data, file = "resTable.xlsx"){
  featureGroup <- data$featureGroup
  featureGroup$adduct <- sapply(featureGroup$adduct, function(x) {paste0(x, collapse = ";")})
  featureGroup <- featureGroup[, -which(colnames(featureGroup) %in% c("pgid", "spectra"))]
  openxlsx::write.xlsx(featureGroup, file = file)
}
