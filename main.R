plexPara <- set_plexPara(targetGroup = "Amine", deltaRt = 40)
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
data <- peakAnnotation(data, polarity = "positive", adinfo = positive.adinfo, thread = 3)
data <- getSpectra2(data = data, thread = 5)
data <- assign_plexInfo(data, plexPara, thread = 3)
data <- assign_tagNum(data, plexPara = plexPara, thread = 3)
data <- assign_tagNum_tolerant(data, plexPara = plexPara, isoRt = 1, thread = 3)
data <- peakGrouping(data = data, plexPara = plexPara, thread = 1)
data <- peakAligning(data = data, plexPara = plexPara)
