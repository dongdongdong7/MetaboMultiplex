#' @title Load AmidoLibrary
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' cmpLibrary <- load_AmidoLibrary()
load_AmidoLibrary <- function(){
  AmidoLibraryPath <- system.file("database/AmidoLibrary.xlsx", package = "MetaboMultiplex")
  AmidoLibrary <- openxlsx::read.xlsx(AmidoLibraryPath)
  return(AmidoLibrary)
}

#' @title Load HydroxylLibrary
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' cmpLibrary <- load_HydroxylLibrary()
load_HydroxylLibrary <- function(){
  HydroxyLibraryPath <- system.file("database/HydroxylLibrary.xlsx", package = "MetaboMultiplex")
  HydroxylLibrary <- openxlsx::read.xlsx(HydroxyLibraryPath)
  return(HydroxylLibrary)
}

#' @title Load PheHydroLibrary
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' cmpLibrary <- load_PheHydroLibrary()
load_PheHydroLibrary <- function(){
  PheHydroLibraryPath <- system.file("database/PheHydroLibrary.xlsx", package = "MetaboMultiplex")
  PheHydroLibrary <- openxlsx::read.xlsx(PheHydroLibraryPath)
  return(PheHydroLibrary)
}

#' @title Load DEANS-BANK
#'
#' @param thread thread.
#'
#' @return A DataFrame object.
#' @export
#'
#' @examples
#' ms2Library <- load_DEANSBANK(thread = 3)
load_DEANSBANK <- function(thread = 1){
  DEANSBANKPath <- system.file("database/spectraLibrary.msp", package = "MetaboMultiplex")
  DEANSBANK <- MsBackendMsp::readMsp(DEANSBANKPath,
                                     mapping = c(name = "Name", accession = "DB#", mslevel = "MsLevel", num_pekas = "Num Peaks"),
                                     BPPARAM = BiocParallel::SnowParam(workers = thread))
  return(DEANSBANK)
}

