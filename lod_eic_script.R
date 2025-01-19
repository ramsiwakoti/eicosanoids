lod_eic <- read.csv("lod_eic.csv")

# Define the mapping
mapping <- c(
  `13,14-dihydro-15-keto Prostaglandin E2` = "dipge2",
  `13,14-dihydro-15-keto Prostaglandin D2` = "dipgd2",
  `13,14-dihydro-15-keto Prostaglandin F2α` = "dipgf2",
  `13,14-dihydro-15-keto Prostaglandin J2` = "dipgj2",
  `15-deoxy-Δ12,14-Prostaglandin J2` = "deopgj2",
  `Prostaglandin A2` = "pga2",
  `Prostaglandin B2` = "pgb2",
  `Prostaglandin E3` = "pge3",
  `Prostaglandin D3` = "pgd3",
  `Prostaglandin J2` = "pgj2",
  `Bicyclo Prostaglandin E2` = "bc_pge2",
  `Bicyclo Prostaglandin E1` = "bc_pge1",
  `Prostaglandin E1 (power)` = "pge1",
  `Prostaglandin D2` = "pgd2",
  `Prostaglandin E2` = "pge2",
  `Thromboxane B2` = "txb2",
  `9-OxoODE` = "x9_oxoode",
  `11(S)-HETE` = "x11hete",
  `16(S)-HETE` = "x16hete",
  `17(S)-HETE` = "x17hete",
  `(±)18-HETE` = "x18hete",
  `20-carboxy Arachidonic Acid` = "x20_caa",
  `20(S)-HETE` = "x20_hete",
  `(±)12,13-DiHOME` = "x12_13_dihome",
  `(±)9,10-DiHOME` = "x9_10_dihome",
  `11(12)-EET` = "x11_12_eet",
  `14(15)-EET` = "x14_15_eet",
  `5(6)-EET` = "x5_6_eet",
  `8(9)-EET` = "x8_9_eet",
  `(±)8,9-DHET` = "x8_9_dhet",
  `(±)11,12-DHET` = "x11_12_dhet",
  `9(10)-EpOME` = "x9epome",
  `12(13)-EpOME` = "x12epome",
  `(±)5,6-DHET` = "x5_6_dhet",
  `9s-HODE` = "x9s_hode",
  `8(S)-HETE` = "x8hete",
  `5(S)-HETE` = "x5_hete",
  `12(S)-HETE` = "x12_hete",
  `15(S)-HETE` = "x15_hete",
  `Leukotriene B4` = "ltb4",
  `Leukotriene C4 methyl ester` = "ltc4.me",
  `Leukotriene D4` = "ltd4",
  `Leukotriene E4` = "lte4",
  `Resolvin D1` = "rvd1",
  `Resolvin D2` = "rvd2",
  `13-OxoODE` = "x13_oxoode",
  `13S-HODE` = "x13s_hode",
  `5-OxoETE` = "x5oxoete",
  `12-OxoETE` = "x12oxoete",
  `15-OxoETE` = "x15oxoete",
  `Eicosapentaenoic Acid` = "epa",
  `Docosahexaenoic Acid` = "dha",
  `Arachidonic Acid` = "aa",
  `(α)-Linolenic Acid` = "ala",
  `Linoleic Acid` = "la"
)

# Add a new variable to `lod_eic` based on the mapping
lod_eic <- lod_eic %>%
  mutate(
    analyte_format = names(mapping)[match(analyte, mapping)]
  ) %>% select(-X)

#Since lte4 and Leukotriene C4 methyl ester are not in the lod file, add them

extra_lod <- data.frame(analyte = c("lte4", "ltc4.me"),
                        LOD = c(0.2, 2),
                        analyte_format = c("Leukotriene E4", "Leukotriene C4 methyl ester"))

lod_eic_updated <- rbind.data.frame(lod_eic, extra_lod)

write.csv(lod_eic_updated, "lod_eic_updated.csv")
