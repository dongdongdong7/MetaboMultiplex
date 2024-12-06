#' @title Visualise a peak group
#' @description
#' Visualise a peak group with four types, chromatograms, spectra, cross-correlation and pearson.
#'
#' @param data data list.
#' @param pgid pgid.
#' @param type Visualisation type.
#'
#' @return No return.
#' @export
#'
#' @examples
#' plot_peakGroup(data = data_AP_mix1_Amine, pgid = "PG049", type = "chromatograms")
plot_peakGroup <- function(data, pgid, type = c("chromatograms", "spectra", "cross-correlation", "pearson")[1]){
  peaksGroup <- data$peakGroup
  i <- which(peaksGroup$pgid == pgid)
  if(length(i) == 0) stop("Can not find pgid!")
  peaks <- peaksGroup[which(peaksGroup$pgid == pgid), ]$peaks[[1]]
  chrList <- lapply(1:nrow(peaks), function(i) {
    cpid <- peaks[i, ]$cpid
    if(is.na(cpid)) return(NULL)
    else{
      chr_tmp <- xcms::chromPeakChromatograms(object = data$rawData, peaks = cpid)[1]
      return(chr_tmp)
    }
  })
  if(type == "chromatograms"){
    rtime_vec <- lapply(chrList, function(x){
      if(is.null(x)) return(mean(peaks$rt, na.rm = TRUE))
      return(x@rtime)
    })
    intensity_vec <- lapply(chrList, function(x){
      if(is.null(x)) return(0)
      return(x@intensity)
    })
    plex_chr <- sapply(1:length(chrList), function(x) paste0("PLEX", x))
    plex_vec <- lapply(1:length(chrList), function(x) {
      rep(plex_chr[x], length(rtime_vec[[x]]))
    })
    rtime_vec <- purrr::list_c(rtime_vec)
    intensity_vec <- purrr::list_c(intensity_vec)
    plex_vec <- purrr::list_c(plex_vec)
    plex_vec <- factor(plex_vec, levels = sort(plex_chr, decreasing = TRUE))
    df <- data.frame(rtime = rtime_vec, intensity = intensity_vec, plex = plex_vec) %>%
      dplyr::filter(!is.na(intensity))
    ggplot2::ggplot(df, mapping = ggplot2::aes(x = rtime, y = plex, height = intensity, group = plex, fill = plex, color = plex)) +
      ggridges::geom_density_ridges(stat = "identity", scale = 1) +
      ggplot2::xlab(label = "Retention Time") +
      ggplot2::ylab(label = "Plex Index") +
      ggplot2::theme_bw() +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
      ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
      ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 15),plot.margin = ggplot2::unit(rep(0.2, 4), "cm"))
  }
  else if(type == "spectra"){
    specialFrag_vec <- peaks$specialFragment
    sp_list <- peaks$spectra
    plex_vec <- sapply(1:length(chrList), function(x) paste0("PLEX", x))
    par(mfrow = c(ceiling(length(plex_vec) / 2),2), mar = c(2,2,2,2))
    for(i in 1:6){
      if(is.null(sp_list[[i]])) plot(0,0, type = "n", main = plex_vec[i])
      else{
        test_label <- rep(NA, length(Spectra::mz(sp_list[[i]])[[1]]))
        test_label[which(dplyr::near(Spectra::mz(sp_list[[i]])[[1]], specialFrag_vec[i], tol = 0.01))] <- specialFrag_vec[i]
        Spectra::plotSpectra(sp_list[[i]], labels = test_label, main = plex_vec[i], labelCol = "red",cex.axis = 1, labelCex = 1.5, xlab = NULL, ylab = NULL)
      }
    }
  }
  else if(type == "cross-correlation"){
    plex_vec <- sapply(1:length(chrList), function(x) paste0("PLEX", x))
    cc_matrix <- matrix(NA, nrow = length(chrList), ncol = length(chrList))
    rownames(cc_matrix) <- plex_vec
    colnames(cc_matrix) <- plex_vec
    for(i in 1:length(chrList)){
      for(j in 1:length(chrList)){
        if(is.null(chrList[[i]]) | is.null(chrList[[j]])){
          cc_matrix[i, j] <- 0
        }else{
          cc_score <- .comparePeaks(chr1 = chrList[[i]], chr2 = chrList[[j]], method = "cross-correlation", apexAlign = TRUE, plot = FALSE)
          cc_matrix[i, j] <- cc_score
        }
      }
    }
    get_lower_tri <- function(mat){
      mat[upper.tri(mat)] <- NA
      return(mat)
    }
    cc_matrix_lower_tri <- get_lower_tri(cc_matrix)
    cc_df <- reshape2::melt(cc_matrix_lower_tri, na.rm = TRUE)
    ggplot2::ggplot(data = cc_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "gray", high = "red",
                                    midpoint = 0, limit = c(0, 1), space = "Lab",
                                    name = "Cross-Correlation") +
      ggplot2::theme_minimal() +
      ggplot2::geom_text(ggplot2::aes(Var1, Var2, label = round(value, 2)), color = "black", size = 4) +
      ggplot2::coord_fixed() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.45, 0.75),
        legend.direction = "horizontal"
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", barwidth = 8, barheight = 1, title.hjust = 0.5))
  }
  else if(type == "pearson"){
    plex_vec <- sapply(1:length(chrList), function(x) paste0("PLEX", x))
    ps_matrix <- matrix(NA, nrow = length(chrList), ncol = length(chrList))
    rownames(ps_matrix) <- plex_vec
    colnames(ps_matrix) <- plex_vec
    for(i in 1:length(chrList)){
      for(j in 1:length(chrList)){
        if(is.null(chrList[[i]]) | is.null(chrList[[j]])){
          ps_matrix[i, j] <- 0
        }else{
          cc_score <- .comparePeaks(chr1 = chrList[[i]], chr2 = chrList[[j]], method = "pearson", apexAlign = TRUE, plot = FALSE)
          ps_matrix[i, j] <- cc_score
        }
      }
    }
    get_lower_tri <- function(mat){
      mat[upper.tri(mat)] <- NA
      return(mat)
    }
    ps_matrix_lower_tri <- get_lower_tri(ps_matrix)
    ps_df <- reshape2::melt(ps_matrix_lower_tri, na.rm = TRUE)
    ggplot2::ggplot(data = ps_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "gray", high = "red",
                                    midpoint = 0, limit = c(0, 1), space = "Lab",
                                    name = "Pearson") +
      ggplot2::theme_minimal() +
      ggplot2::geom_text(ggplot2::aes(Var1, Var2, label = round(value, 2)), color = "black", size = 4) +
      ggplot2::coord_fixed() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.45, 0.75),
        legend.direction = "horizontal"
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", barwidth = 8, barheight = 1, title.hjust = 0.5))
  }
  else stop("Type is wrong!")

}
