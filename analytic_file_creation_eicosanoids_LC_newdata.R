#In this script, I will create analytic files necessary to run PFAS-eicosanoids analysis project in LIFECODES only

#Load some libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(naniar)
library(corrplot)


#Import data

#Lifecodes - eicosanoids 
eic_lc_old <- readxl::read_xlsx(path = "C:\\Users\\rsiwa\\Dropbox (University of Michigan)\\Past Analysis and Data\\data archive\\Phase 1 data\\Eicosanoids\\174 Samples for Eicosanoids (1).xlsx",
                            sheet = "Data")
eic_lc_old2 <- eic_lc_old %>% filter(!is.na(ID)) %>% 
  select(-c(`Sample No....1`, Batch, Visit, `Box position`, `Sample No....6`))

eic_lc_old2_batch <- eic_lc_old %>% select(ID, Batch)
names(eic_lc_old2_batch) <- eic_lc_old


eic_lc_old2 <- eic_lc_old2 %>% mutate_if(is.character, as.numeric)
eic_lc_old2 <- as.data.frame(eic_lc_old2)
names(eic_lc_old2) <- tolower(names(eic_lc_old2))

eic_lc_old3 <- eic_lc_old2 %>%
  rename(
    dipgf2       = `12_13,14dhkpgf2a`,
    dipgd2       = `8_13,14dpgd2`,
    pgb2         = `53_pgb2`,
    x15_hete     = `150_15hete`,
    x12_hete     = `148_12hete`,
    x5_hete      = `146_5hete`,
    x20_hete     = `152_20hete`,
    x11hete      = `68_11hete *`,
    x16hete      = `69_16hete *`,
    x17hete      = `70_17hete *`,
    x18hete      = `71_18hete*`,
    x8hete       = `66_8hete`,
    aa           = `144_aa`,
    epa          = `115epa`,
    dha          = `118dha`,
    x20_caa      = `121caa*`,
    txb2         = `134_ txb2`,
    pge2         = `130_ pge2`,
    pgd2         = `128_ pgd2`,
    pga2         = `51_pga2`,
    deopgj2      = `19_15do12.14pgj2`,
    pge1         = `60_pge1`,
    bc_pge1      = `59_bcpge1`,
    bc_pge2      = `58_bcpge2`,
    pgd3         = `55_pgd3`,
    pge3         = `56_pge3`,
    pgj2         = `57_pgj2`,
    rvd2         = `43 rvd2`,
    rvd1         = `41 rvd1`,
    x12_13_dihome = `1_12,13-dihome`,
    x9_10_dihome  = `3_9,10-dihome`,
    x11_12_dhet   = `94_11,12_dhet`,
    x8_9_dhet     = `93_8,9_dhet`,
    x13_oxoode    = `98_13-oxoode`,
    x9_oxoode     = `96_9-oxoode`,
    x5_6_dhet     = `110_5_6_dhet`,
    x8_9_eet      = `25-8(9)-eet`,
    x14_15_eet    = `15_14(15)-eet`,
    x11_12_eet    = `5_11(12)-eet`,
    x5_6_eet      = `21_5(6)-eet`,
    x12epome      = `105_12_13_epome`,
    x9epome       = `103_9_10_epome`,
    ltd4          = `36_ltd4`,
    lte4          = `38_lte4`,
    ltb4          = `30_ltb4`,
    ltc4.me       = `33_ltc4 me`,
    x13s_hode     = `140_13s-hode`,
    x9s_hode      = `138_9s-hode`,
    x5oxoete      = `155_5oxoete`,
    x15oxoete     = `157_15oxoete`,
    x12oxoete     = `156_12oxoete`,
    ala           = `119_ala`,
    la            = `142_la`
  )

eic_lc_new <- read.csv(file= "C:\\Users\\rsiwa\\Dropbox (University of Michigan)\\Past Analysis and Data\\data archive\\Phase 1 data\\Eicosanoids phase 1 2024\\eicosanoids_cleanlong_010725.csv")
eic_lc_new$analyte <- tolower(eic_lc_new$analyte)
eic_lc_new <- eic_lc_new %>% mutate(analyte = ifelse(analyte == "x18_epa", "epa", analyte)) %>%
  rename("id" = "participant_id")

#Look at LOD

det_rate_eic <- eic_lc_new %>% group_by(analyte) %>% summarize(detrate = 1 - mean(BLOD, na.rm = TRUE))
new_eic <- det_rate_eic %>% select(analyte)

lod_eic <- eic_lc_new %>% select(analyte, analyte_full, LOD) %>% distinct()



#Compare the old vs new data

old_eic <- data.frame(analyte = names(eic_lc_old3))

new_old_comp <- old_eic %>% inner_join(new_eic, by = "analyte") #50 common eicosanoids
inold_notnew <- old_eic %>% anti_join(new_eic, by = "analyte")
innew_notold <- new_eic %>% anti_join(old_eic, by = "analyte")

#New file has DiPGE2 and DiPGJ2, which are also present in PROTECT data but not older LIFECODES data
#Old file has lte4 and ltc4.me, which are not present in newer LIFECODES data

#See how many have duplicate PFAS data

id_eic_new <- eic_lc_new %>% select(id, visit_id) %>% 
  group_by(id, visit_id) %>% tally() # 382 ids with visitid = 1 and 156 ids with visitid = 3
table(id_eic_new$visit_id)


#Check if the same ID has two records from 1 visit 

id_eic_new_chk0 <- eic_lc_new %>% group_by(id, visit_id) %>% tally()

#The answer is Yes (IDs 1408 and 668)
#
##Check values for these IDs

id_1408_668 <- eic_lc_new %>% filter(id %in% c(1408, 668))


#Check how many have one vs two records
id_eic_new_chk <- id_eic_new %>% group_by(id) %>% tally()
table(id_eic_new_chk$n)

#Check how many are duplicated with the old records 
id_eic_old <- eic_lc_old3 %>% select(id)
id_overlap <- id_eic_new_chk0 %>% inner_join(id_eic_old, by = c("id"))

#IDs in old records but not new
id_only_old <- id_eic_old %>% anti_join(id_eic_new_chk0, by = c("id"))
id_only_old %>% inner_join(eic_lc_new, by = "id")

#128 IDs overlap. But none of these IDs have duplicate records in the new file. All the overlaps have Visit 1 data

# Out of 264 participants with non-duplicates, check the distributions of visit 1 vs 3

eic_new_n_chk <- id_eic_new_chk %>% filter(n == 1) %>% anti_join(id_eic_old, by = "id")
eic_new_n_chk2 <- eic_new_n_chk %>% inner_join(id_eic_new, by = "id")
table(eic_new_n_chk2$visit_id)

#Check how many records have notes on low reliability
notes_eic_new <- eic_lc_new %>% select(id, visit_id, analyte, notes) %>% filter(!is.na(notes))
notes_eic_new_analyte <- notes_eic_new %>% group_by(analyte) %>% tally()


#Number of participants with this warning
notes_wanrning_n <- notes_eic_new %>% group_by(id) %>% tally()
#These warnings apply for all samples or participants. 

#Check the detection rate by visit
#
det_rate_eic_byvisit <- eic_lc_new %>% group_by(analyte, visit_id) %>% 
  summarize(detrate = 1 - mean(BLOD, na.rm = TRUE)) %>%
  arrange(analyte, visit_id, desc(detrate))

#IDs without detection rate because they do not have LODs

det_rate_nan <- det_rate_eic_byvisit %>% filter(is.nan(detrate))

#Det rate wide 

det_rate_eic_wide <- det_rate_eic_byvisit %>% 
  pivot_wider(names_from = visit_id, values_from = detrate, names_prefix = "visit") %>%
  mutate(det_comp = ifelse(visit1 >= 0.7 & visit3 >= 0.7, 1, 0))

#See if the analytes detected below 70% are also the ones with unreliable peaks

det_peak_chk <- det_rate_eic_wide %>% inner_join(notes_eic_new_analyte, by = "analyte")

## Now look at the correlation between repeated measurements

## In new file, we will look at correlation between repeated measurements 

eic_to_include <- det_rate_eic_wide %>% filter((det_comp == 1|is.na(det_comp)))
eic_to_exclude <- det_rate_eic_wide %>% filter((det_comp == 0))

eic_lc_new_short <- eic_lc_new %>% select(id, visit_id, analyte, conc) %>%
  pivot_wider(names_from = c("visit_id", "analyte"), values_from = conc, names_prefix = "v") %>%
  select(id, starts_with("v1"), starts_with("v3"))

eic_lc_new_short2 <-  eic_lc_new %>% select(id, visit_id, analyte, conc) %>% filter(!(id %in% c(1408, 668))) %>%
  pivot_wider(names_from = c("visit_id", "analyte"), values_from = conc, names_prefix = "v") %>%
  select(id, starts_with("v1"), starts_with("v3"))

# Visualize missing values

visdat::vis_miss(eic_lc_new_short2) #looks little weird because for v3, we dont have data for all participants 

#Look at the correlation now

cor_all_eic <- cor((eic_lc_new_short2 %>% select(starts_with("v1"))), (eic_lc_new_short2 %>% 
                                                                         select(starts_with("v3"))), method = "spearman",
    use = "pairwise.complete.obs")

cor_all_eic_diag <- data.frame(analyte = rownames(cor_all_eic), corr = diag(cor_all_eic))

cor(eic_lc_new_short2$v1_aa, eic_lc_new_short2$v3_aa, method = "spearman", use = "pairwise.complete.obs")

## Now calculate correlation after combining old and new data


##Convert old data to the long-form

eic_lc_old3_long <- eic_lc_old3 %>%
  pivot_longer(cols = dipgf2:la, names_to = "analyte", values_to = "conc") %>%
  mutate(visit_id = 





#Lifecodes - PFAS
pfas_lc <- read.csv("C:/Users/rsiwa/OneDrive/Michigan_course/Lifecodes/pfas_phase2/PFASultrasound/pfas_covars_phase1.csv") 
names(pfas_lc) <- tolower(names(pfas_lc))


pfas_lc2 <- pfas_lc  %>% 
  select(id, starts_with("conc_"), starts_with("blod_")) %>%  
  mutate(cohort = "LC") 
names(pfas_lc2)