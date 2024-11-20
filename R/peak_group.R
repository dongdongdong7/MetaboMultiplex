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
