#In this script, I will create analytic files necessary to run PFAS-eicosanoids analysis project

#Load some libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(naniar)
library(corrplot)


#Import data

#Lifecodes - eicosanoids 
eic_lc <- readxl::read_xlsx(path = "C:\\Users\\rsiwa\\Dropbox (University of Michigan)\\Past Analysis and Data\\Eicosanoids\\174 Samples for Eicosanoids (1).xlsx",
                            sheet = "Data")
eic_lc2 <- eic_lc %>% filter(!is.na(ID)) %>% 
  select(-c(`Sample No....1`, Batch, Visit, `Box position`, `Sample No....6`))

eic_lc2 <- eic_lc2 %>% mutate_if(is.character, as.numeric)
eic_lc2 <- as.data.frame(eic_lc2)
names(eic_lc2) <- tolower(names(eic_lc2))

#Lifecodes - PFAS
pfas_lc <- read.csv("C:/Users/rsiwa/OneDrive/Michigan_course/Lifecodes/pfas_phase2/PFASultrasound/pfas_covars_phase1.csv") 
names(pfas_lc) <- tolower(names(pfas_lc))


pfas_lc2 <- pfas_lc  %>% 
  select(id, starts_with("conc_"), starts_with("blod_")) %>%  
  mutate(cohort = "LC") 
names(pfas_lc2)

#Protect-eicosanoids

eic_pr <- read.csv("C:\\Users\\rsiwa\\OneDrive\\Michigan_course\\eicosanoids\\eic.data.csv") %>%
  select(-X, -visitid) %>% rename("id" = "studyid")
names(eic_pr) <- tolower(names(eic_pr))
names(eic_pr) 

#63 eicosanoids including the parent compound

#Protect-PFAS
# pfas_pr <- read.csv("C:\\Users\\rsiwa\\OneDrive\\Michigan_course\\PROTECTpfas\\os_inf_pfas_wide_short3.csv") %>%
#   select(-X) %>% 
#   select(id, starts_with("conc_"), starts_with("blod_")) %>% 
#   distinct()

pfas_pr <- read.csv("C:\\Users\\rsiwa\\OneDrive\\Michigan_course\\PROTECTpfas\\pfas_data.csv") %>%
  select(-X) %>% 
  select(id, starts_with("conc_"), starts_with("blod_")) %>% 
  rename("blod_me.pfosa.acoh" = "blod_mpah",
         "blod_pfdea" = "blod_pfde",
         "blod_pfhpa" = "blod_pfhp",
         "blod_pfhxs" = "blod_pfhs") %>%
  distinct()

names(pfas_pr) <- tolower(names(pfas_pr))

pfas_pr2 <- pfas_pr %>% 
  mutate(cohort = "PR")
names(pfas_pr2) 

pfas_lc2_to_combine <- pfas_lc2 %>% select(all_of(names(pfas_pr2)))

#Look at the names

names(pfas_lc2)
names(pfas_lc2_to_combine)

pfas_lc_pr <- rbind.data.frame(pfas_lc2_to_combine, pfas_pr2)


#Rename variables in the LC eicosanoids file 

eic_lc3 <- eic_lc2 %>%
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

dim(eic_lc)
dim(eic_lc2)
dim(eic_lc3)
dim(eic_pr)

# Variables in eic_pr, but not in eic_lc

eic_pr_vars <- data.frame(var = names(eic_pr))
eic_lc_vars <- data.frame(var = names(eic_lc3))

eic_pr_lc <- eic_pr_vars %>% inner_join(eic_lc_vars, by = "var")
eic_pr_nolc <- eic_pr_vars %>% anti_join(eic_lc_vars, by = "var")

eic_pr_nolc_dat  <- eic_pr %>% select(id, all_of(eic_pr_nolc$var))

#Look at missing values
visdat::vis_miss(eic_pr_nolc_dat)

#85% missing values for all variables there except dipge2 and dipgj2 

eic_pr2 <- eic_pr %>% select(all_of(names(eic_lc3)), dipge2, dipgj2) %>%
  mutate(cohort = "PR")


# Look at missing values in other datasets
visdat::vis_miss(eic_lc3)
visdat::vis_miss(eic_pr2)

#ltc4.me has 15% missing data, but all others do not

#Since LC file does not have two variables: pge2 and pgj2, 

#I will add variables pge2 and pgj2 to the LC file

eic_lc4 <- eic_lc3 %>% mutate(dipge2 = NA, dipgj2 = NA, cohort = "LC")

eic_lc_pr <- rbind.data.frame(eic_pr2, eic_lc4)

dim(eic_lc4)
dim(eic_pr2)
dim(eic_lc_pr)

#we need to also harmonize covariate files

#LC

covars_lc <- pfas_lc %>%
  mutate(ga_samp = ifelse(visitid == 1, ga_t1, ga_t2)) %>%
  select(id, visitid, year_t1, age, bmi_prepreg, mat_edu_cat, race, ins_cat, alc, alc_cat,
         smpreg, gender1, nullip, parity, starts_with("ptb"), ga_samp, ga_t1, ga_t2) 
covars_lc2 <- covars_lc %>% select(-ga_t1, -ga_t2, -visitid) %>%
  mutate(race4 = ifelse(race %in% c("Asian", "Other", "Mixed"), "Asian_other", race))

covars_pr <- read.csv("C:\\Users\\rsiwa\\OneDrive\\Michigan_course\\PROTECTpfas\\covars_new2b.csv") %>% select(-X)
names(covars_pr)

#To match some of the covariates name

covars_pr2 <- covars_pr %>%
  rename("mat_edu_cat" = "ed_stat_cat",
         "gender1b" = "ppsex",
         "ptb" = "ptb_status",
         "year_t1" = "year",
         "bmi_prepreg" = "prebmi",
         "age" = "isage",
         "ins_cat" = "ins_stat2",
         "smpregb" = "smksmellcat",
         "ga_samp" = "ga_samp_pfas")


covars_pr2_tokeep <- covars_pr2 %>%
  select(id, mat_edu_cat, gender1b, ptb, year_t1, bmi_prepreg, age,  ins_cat, pregnum_recode, smpregb, ga_samp)

#Harmonization
#Education: Not needed - consistent recoding in both files
#gender1: 1 = female, 2 = male
#PTB: Not needed
#year_t1: Not needed
#bmi_prepre: Not needed
#age: Not needed
#ins_cat: 1 = private; 2 = public (Not needed)
#pregnum_recode: 0 = not nulliparous, 1 = nulliparous [Needs recoding - right now, pregnum = 1 means first time pregnant, and 
#pregnum = 2 means not first time pregnant]

#smpreg = 1 = yes smoking during pregnnacy; 2 = No, not smoking during pregnancy. Needs reocding - right now 1 means never 2,
#means ever, and 3 means current smoker
#Just keep recode smpreg = 3 to smpreg = 1, else 0 since they are not current smoker

#ga_samp: Not needed

covars_pr2_tokeep2 <- covars_pr2_tokeep %>%
  mutate(gender1 = case_when(gender1b == 1 ~ 0,
                             gender1b == 2 ~ 1,
                             .default = gender1b),
         nullip = case_when(pregnum_recode == 1 ~ 1,
                            pregnum_recode == 2 ~ 0,
                            .default = pregnum_recode),
         smpreg = case_when(smpregb == 3 ~ 1,
                            smpregb == 1 ~ 0,
                            smpregb == 2 ~ 0,
                            .default = smpregb),
         race4 = "Hispanic-PR") %>%
  select(-gender1b, -pregnum_recode, -smpregb)

##Now lets combine PFAS, eicosanoids, and covariates file for each cohort

covars_lc2_tokeep <- covars_lc2 %>% select(id, age, bmi_prepreg, mat_edu_cat, ins_cat,  race4, 
                                    ga_samp, nullip, gender1, year_t1, smpreg, ptb) %>%
  mutate(cohort = "LC")


covars_pr2_tokeep2 <- covars_pr2_tokeep2 %>% select(id, age, bmi_prepreg, mat_edu_cat, ins_cat,  race4, 
                                    ga_samp, nullip, gender1, year_t1, smpreg, ptb) %>%
  mutate(cohort = "PR")

covars_lc_pr <- rbind.data.frame(covars_lc2_tokeep, covars_pr2_tokeep2)


## Combine PFAS and eicosanoids file

pfas_eic_combined <- pfas_lc_pr %>% 
  inner_join(eic_lc_pr, by = c("id", "cohort"))

#export analytic files

write.csv(pfas_eic_combined, "pfas_eic_pr_lc_comb.csv")
write.csv(covars_lc_pr, "covars_pr_lc_comb.csv")
write.csv(covars_lc, "covars_lc.csv")
write.csv(covars_pr2, "covars_pr.csv")

















      

