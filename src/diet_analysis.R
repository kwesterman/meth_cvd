foodLib <- c(satfat="NUT_SATFAT", monfat="NUT_MONFAT", polyfat="NUT_POLY", beans="FFD60",
             alcohol="NUT_ALCO", sucrose="NUT_SUCR", sodium="NUT_SODIUM", n3="NUT_OMEGA",
             folate="NUT_FOLEQ", sugars="NUT_SUGTOT", transfat="NUT_TRN02", ssb="FFD145",
             proantho="NUT_PROMON", ecg="NUT_UECG", kale="FFD66", prcsdmeat="FFD78")

library(readxl)
dietDat_all <- read_excel("../data/diet/phs000007.v28.pht002350.v4.p10.c1.vr_ffreq_ex08_1_0615s.HMB-IRB-MDS_ex8_diet.xlsx", sheet=2)
dietDat <- dplyr::select(dietDat_all, shareid, foodLib)
names(dietDat) <- c("shareid", names(foodLib))

load("../int/mrsDat2.RData")
dietDat <- left_join(dplyr::select(mrsDat, shareid, mrs), dietDat, by="shareid")
dietMRScors <- cor(dietDat[-1], use="pairwise.complete.obs")
heatmap.2(dietMRScors, Rowv=F, Colv=F)
