# AP Test
{
  # Prepare data
  {
    plexPara_Amine <- set_plexPara(targetGroup = "Amine", deltaRt = 40)
    plexPara_Phenol <- set_plexPara(targetGroup = "Phenol", deltaRt = 40)
    xcmsPara = set_xcmsPara(ppm = 30, peakwidth = c(4, 30), snthresh = 1,
                            noise = 100, prefilter = c(3, 100), firstBaselineCheck = FALSE,
                            expandRt = 2, expandMz = 0.01, minProp = 0.75)
    dataPara <- set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/AP/mix1/",
                             res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
                             sampleData = data.frame(sample_id = c("mix1_1", "mix1_2", "mix1_3"), injection_index = 1:3))
    data <- load_data(dataPara = dataPara)
    data <- peakPicking(data, xcmsPara = xcmsPara, chunkSize = 3, BPPARAM = BiocParallel::SnowParam(workers = 3))
    data(positive.adinfo, package = "cliqueMS")
    positive.adinfo <- positive.adinfo[positive.adinfo$adduct %in%
                                         c("[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+H-NH3]+", "[M+K]+", "[M+NH4]+"),]
    data <- peakAnnotation(data, polarity = "positive", adinfo = positive.adinfo, chunkSize = 2, thread = 2)
    data <- getSpectra2(data = data, thread = 5)
    data <- assign_tagNum(data, thread = 3)
    data <- assign_tagNum_tolerant(data, thread = 3)
    ms2Library <- load_DEANSBANK(thread = 3)
  }
  # Amine
  {
    data_Amine <- assign_plexInfo(data, plexPara = plexPara_Amine, thread = 3)
    data_Amine <- peakGrouping(data = data_Amine, plexPara = plexPara_Amine, thread = 1)
    data_Amine <- peakAligning(data = data_Amine, plexPara = plexPara_Amine)
    cmpLibrary_Amine <- load_AmidoLibrary()
    data_Amine <- fgIdentification(data = data_Amine, cmpLibrary = cmpLibrary_Amine, ms2Library = ms2Library)
    data_Amine <- getQuantRes(data_Amine, plexPara = plexPara_Amine, quant = "into")
    exportRes(data = data_Amine, file = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp/resTable_Amine.xlsx")
  }
  # Phenol
  {
    data_Phenol <- assign_plexInfo(data, plexPara = plexPara_Phenol, thread = 3)
    data_Phenol <- peakGrouping(data = data_Phenol, plexPara = plexPara_Phenol, thread = 1)
    data_Phenol <- peakAligning(data = data_Phenol, plexPara = plexPara_Phenol)
    cmpLibrary_Phenol <- load_PheHydroLibrary()
    data_Phenol <- fgIdentification(data = data_Phenol, cmpLibrary = cmpLibrary_Phenol, ms2Library = ms2Library)
    data_Phenol <- getQuantRes(data_Phenol, plexPara = plexPara_Phenol, quant = "into")
    exportRes(data = data_Phenol, file = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp/resTable_Phenol.xlsx")
  }
}

# Hy Test
{
  # Prepare data
  {
    plexPara_Alcohol <- set_plexPara(targetGroup = "Alcohol", deltaRt = 40)
    xcmsPara = set_xcmsPara(ppm = 30, peakwidth = c(4, 30), snthresh = 1,
                            noise = 100, prefilter = c(3, 100), firstBaselineCheck = FALSE,
                            expandRt = 2, expandMz = 0.01, minProp = 0.75)
    dataPara <- set_dataPara(data_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/test_data/Hy/mix1/",
                             res_dir = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp",
                             sampleData = data.frame(sample_id = c("mix1_1"), injection_index = 1))
    data <- load_data(dataPara = dataPara)
    data <- peakPicking(data, xcmsPara = xcmsPara, chunkSize = 1, BPPARAM = BiocParallel::SerialParam())
    data(positive.adinfo, package = "cliqueMS")
    positive.adinfo <- positive.adinfo[positive.adinfo$adduct %in%
                                         c("[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+H-NH3]+", "[M+K]+", "[M+NH4]+"),]
    data <- peakAnnotation(data, polarity = "positive", adinfo = positive.adinfo, chunkSize = 1, thread = 1)
    data <- getSpectra2(data = data, thread = 5)
    data <- assign_tagNum(data, thread = 3)
    data <- assign_tagNum_tolerant(data, thread = 3)
    ms2Library <- load_DEANSBANK(thread = 3)
  }
  # Alcohol
  {
    data_Alcohol <- assign_plexInfo(data, plexPara = plexPara_Alcohol, thread = 3)
    data_Alcohol <- peakGrouping(data = data_Alcohol, plexPara = plexPara_Alcohol, thread = 1)
    data_Alcohol <- peakAligning(data = data_Alcohol, plexPara = plexPara_Alcohol)
    cmpLibrary_Alcohol <- load_HydroxylLibrary()
    data_Alcohol <- fgIdentification(data = data_Alcohol, cmpLibrary = cmpLibrary_Alcohol, ms2Library = ms2Library)
    data_Alcohol <- getQuantRes(data_Alcohol, plexPara = plexPara_Alcohol, quant = "into")
    exportRes(data = data_Alcohol, file = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/tmp/resTable_Alcohol.xlsx")
  }
}
