##### Bioenergetic capacity - Script 5 #####
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Dec 18, 2024              ##
#                                         ##
############################################

#### Front matter #-------------------------------------------------------------
library(ggpubr)
library(lmerTest)
library(dplyr)
library(effects)
library(magrittr)
library(grid)
library(gridExtra)
library(emmeans)
library(metap)

# clear environment
rm(list = ls())

# load helper function for multi-SNP association tests
source("helper_scripts/CalculateGenoCombinations_Function.R")

#### load required data #-------------------------------------------------------
load("data/Script-2.ADNI.genetics.RData") # Contents described in Script 2
load("data/store/ADNI.bcage.RData") # ADNI predicted bioenergetic age
load("data/store/AGES.dt.bcage.RData") # AGES data with predicted bioenergetic age
load("data/Script-5.ADNI.ROSMAP.longitudinal.RData") # contains three data frames:
# ----------
# adni.long: Longitudinal phenotypic data from ADNI from baseline to 5 years
# ----------
# df with dimensions 6496 (longitudinal data, identified by RID and VISCODE) x 15
# (IDs, basic covariates, and outcomes)
#
# --------
# AGES.dt: preprocessed acylcarnitine (AC) profiles, phenotypes and covariates in AGES
# --------
# df with dimensions 575 (samples, identified by sample.id) x 61
# (IDs, covariates, phenotypes, + ACs)

#### Derive data/formats required by this script #------------------------------
## ADNI ------
# Add alleles of 2-SNP combination
adni.data <- merge(adni.bcage, 
                   data.frame(RID = rownames(ADNI.dt) %>% as.numeric(),
                              GenoClust = paste0(ADNI.dt$rs924135, 
                                                 ADNI.dt$rs17806888)), 
                   by = "RID", all = T) %>%
  # set rare SNP combinations (less frequent than 11 observations) to NA 
  group_by(GenoClust) %>% 
  mutate(GenoClust = ifelse(n() > 11, GenoClust, NA)) %>% 
  ungroup() %>% as.data.frame()

# Merge with longitudinal data
adni.long <- merge(adni.data, adni.long, by = "RID", all = T) %>% 
  select(
    RID,
    VISCODE,
    ADAS13,
    ADNI_EF,
    ADNI_MEM,
    MMSE,
    CDRSB,
    FAQ,
    GenoClust,
    Age,
    BCAge,
    Sex,
    ApoE4,
    Education,
    Cohort,
    Diagnosis,
    DX
  )

# recode visit codes to continuous variable
adni.long$VISCODE <- case_when(adni.long$VISCODE == "bl" ~ 0, 
                               adni.long$VISCODE == "m12" ~ 1, 
                               adni.long$VISCODE == "m24" ~ 2, 
                               adni.long$VISCODE == "m36" ~ 3, 
                               adni.long$VISCODE == "m48" ~ 4, 
                               adni.long$VISCODE == "m60" ~ 5)

# add median split of bioenergetic age
adni.long$percentile <- ifelse(adni.long$BCAge < summary(adni.long %>% 
                                                           filter(VISCODE==0) %>% 
                                                           pull(BCAge))[3], 
                               "0-50%", "50-100%")

# introduce baseline values for outcome variables
adni.long <- adni.long %>% group_by(RID) %>% 
  mutate(DX.base = DX[VISCODE==0],
         ADAS13.base = ADAS13[VISCODE==0],
         CDRSB.base = CDRSB[VISCODE==0],
         MMSE.base = MMSE[VISCODE==0],
         FAQ.base = FAQ[VISCODE==0],
         ADNI_EF.base = ADNI_EF[VISCODE==0],
         ADNI_MEM.base = ADNI_MEM[VISCODE==0]) %>% ungroup()

# clean up variables not needed anymore
rm(ADNI.dt, adni.data, adni.bcage, ac.haplotypes, snpAnno)

################################################################################
#### Get set of alleles most significantly associated with cognition -----------
################################################################################

# calculate all possible allele combinations
tests <- calculateGenoCombinations(adni.long) %>% arrange(metaP)

# check effect direction for the most significant allele combination
effdir <- if ({
  tests %>% head(1) %>% pull(Estimate.ADNI_MEM) %>% sign()
} > 0) {
  c("faster", "slower")
} else{
  c("slower", "faster")
}

# store alleles of the most significant combination
alleles <- tests %>% arrange(metaP) %>% head(1) %>% pull(alleles) %>% 
  as.character %>% strsplit(., ",") %>% unlist

# group ADNI participants into faster and slower group
adni.long$GxG <- ifelse(is.na(adni.long$GenoClust), NA, 
                        ifelse(adni.long$GenoClust %in% alleles, 
                               effdir[1], effdir[2]))

# same for ROS/MAP
ROSMAP.long$GxG <- ifelse(is.na(ROSMAP.long$GenoClust), NA, 
                        ifelse(ROSMAP.long$GenoClust %in% alleles, 
                               effdir[1], effdir[2]))

################################################################################
#### Calculate longitudinal associations ---------------------------------------
################################################################################

## For predicted bioenergetic age [BCAge]
# ADNI
statsBC <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    # return results for continuous trait
    form <-
      paste0(y,
             " ~ BCAge * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)")
    summary(lmer(as.formula(form), adni.long))$coefficients %>% tail(1) %>% 
      cbind(data.frame(trait = niceTrait), .) %>% as.data.frame
  }) %>% bind_rows()

# Replication in AGES
stats.ages <-
  rbind(
    summary(lm(
      BCAge ~ cogdx.fu + Age + Sex + Education + ApoE4, AGES.dt[AGES.dt$cogdx.fu %in%
                                                                  c("CN", "MCI"), ]
    ))$coefficients[2, ],
    summary(lm(
      BCAge ~ cogdx.fu + Age + Sex + Education + ApoE4, AGES.dt[AGES.dt$cogdx.fu %in%
                                                                  c("MCI", "AD"), ]
    ))$coefficients[2, ],
    summary(lm(
      BCAge ~ cogdx.fu + Age + Sex + Education + ApoE4, AGES.dt[AGES.dt$cogdx.fu %in%
                                                                  c("CN", "AD"), ]
    ))$coefficients[2, ]
  ) %>%
  as.data.frame %>% arrange(`Pr(>|t|)`) %>% pull(`Pr(>|t|)`) %>% 
  as.numeric() %>% sprintf("%.2g", .)

# Prepare residuals for plotting below
AGES.dt$BCAge.res <- lm(BCAge ~ Age + Sex + Education + ApoE4, AGES.dt)$residuals

## For allele groupings [GXG]
# ADNI
statsGxG <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    # return results for continuous trait
    form <-
      paste0(y,
             " ~ GxG * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)")
    summary(lmer(as.formula(form), adni.long))$coefficients %>% tail(1) %>% 
      cbind(data.frame(trait = niceTrait), .) %>% as.data.frame
  }) %>% bind_rows()

# Replication in ROS/MAP
statsGxG_r <- summary(
  lmer(
    cogn_global ~ GxG * fu_year + Sex + ApoE4 + Education + age_bl + dx_bl + (1 | projid),
    ROSMAP.long
  ))$coefficients %>% tail(1) %>% 
  cbind(data.frame(trait = "Global cognition"), .) %>% as.data.frame

## For the interaction of allele groupings and bioenergetic age [GXGxE]
# Use median split to extract the stats for the contrast of interest:
# those with the protective genotype but high vs. low bioenergetic age
statsGxGxE <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    x <-
      lmer(as.formula(
        paste0(
          y,
          " ~ percentile * GxG * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)"
        )
      ), adni.long)
    EMMs <-
      pairs(
        emmeans(
          x,
          ~ VISCODE * percentile * GxG,
          at = list(VISCODE = c(0, 5)),
          pbkrtest.limit = 5100
        ),
        by = NULL,
        adjust = NULL
      )
    idx <- which(summary(EMMs)$contrast == "(VISCODE0 0-50% slower) - (VISCODE5 0-50% slower)")
    Diff <-
      contrast(EMMs,
               by = NULL,
               "trt.vs.ctrl",
               ref = idx,
               adjust = NULL)
    idx <- which(summary(Diff)$contrast == "((VISCODE0 50-100% slower) - (VISCODE5 50-100% slower)) - ((VISCODE0 0-50% slower) - (VISCODE5 0-50% slower))")
    cbind(data.frame(trait = niceTrait), summary(Diff)[idx, -1])
  }) %>% bind_rows()


################################################################################
#### Simulated clinical trial in ADNI ------------------------------------------
################################################################################

## Create dataset
# filter out CN group (but keep SMCs) and only include subjects 
# with all 3 visits
adni.trial <- adni.long %>% filter(Diagnosis != "CN", VISCODE<3)
adni.trial <-
  adni.trial %>% filter(RID %in% {
    adni.trial %>% group_by(RID) %>% summarise(N = n()) %>% filter(N == 3) %>% pull(RID)
  }) 

# Select lowest 25% and highest 25% of bioenergetic age as treatment arms
adni.trial$quartile <-
  ifelse(
    adni.trial$BCAge < summary(adni.trial %>% filter(VISCODE == 0) %>% pull(BCAge))[2],
    "0-25%",
    ifelse(
      adni.trial$BCAge > summary(adni.trial %>% filter(VISCODE == 0) %>% pull(BCAge))[5],
      "75-100%",
      NA
    )
  )

# Add information on anti-dementia drug intake, CDR, and racial/ethnic group
adni.trial <- merge(adni.trial, adni.meds.cdr.re, by = "RID")
rm(adni.meds.cdr.re) # not needed anymore

# Trial characteristics
trial.characteristics <-
  adni.trial %>% filter(VISCODE == 0,!is.na(quartile)) %>% group_by(quartile) %>%
  summarize(
    Total.n = n(),
    Age = paste0(
      paste0(sprintf("%.2f", mean(Age, na.rm = T)), " (", sprintf("%.2f", Agesd = sd(Age, na.rm = T)), "; "),
      paste(Agerange = range(Age, na.rm = T), collapse = " - "), ")"),
    Female = paste0(sum(Sex == 2), " (", sprintf("%.2f", sum(Sex == 2) / n() * 100), " %)"),
    Non.Hispanic.Whites = paste0(
      sum(raceth == "NH-White"),
      " (", sprintf("%.2f", sum(raceth == "NH-White") / n() * 100), " %)"),
    Education = paste0(
      paste0(
        sprintf("%.2f", mean(Education, na.rm = T)), " (",
        sprintf("%.2f", Educationsd = sd(Education, na.rm = T)), "; "),
      paste(
        Educationrange = range(Education, na.rm = T),
        collapse = " - "), ")"),
    APOE4.carrier = paste0(sum(ApoE4 > 0), " (", sprintf("%.2f", sum(ApoE4 > 0) / n() * 100), " %)"),
    AD.Medication = paste0(sum(antiDementiaDrugs), " (",
      sprintf("%.2f", sum(antiDementiaDrugs) / n() * 100), " %)"),
    Genotype.slower = paste0(
      sum(GxG == "slower", na.rm = T), " (",
      sprintf("%.2f", sum(GxG == "slower", na.rm = T) / n() * 100), " %)"),
    CDR0 = paste0(sum(CDGLOBAL == 0), " (", sprintf("%.2f", sum(CDGLOBAL == 0) / n() * 100), " %)"),
    CDR0.5 = paste0(sum(CDGLOBAL == 0.5), " (", sprintf("%.2f", sum(CDGLOBAL == 0.5) / n() * 100), " %)"),
    CDR1 = paste0(sum(CDGLOBAL == 1), " (", sprintf("%.2f", sum(CDGLOBAL == 1) / n() * 100), " %)"),
    SMC = paste0(sum(DX == "CN"), " (", sprintf("%.2f", sum(DX == "CN") / n() * 100), " %)"),
    MCI = paste0(sum(DX == "MCI"), " (", sprintf("%.2f", sum(DX == "MCI") / n() * 100), " %)"),
    AD = paste0(sum(DX == "Dementia"), " (", sprintf("%.2f", sum(DX == "Dementia") / n() * 100), " %)"),
    CDRSB = paste0(
      paste0(sprintf("%.2f", mean(CDRSB, na.rm = T)), 
             " (", sprintf("%.2f", CDRSBsd = sd(CDRSB, na.rm = T)), "; "),
      paste(CDRSBrange = range(CDRSB, na.rm = T), collapse = " - "), ")"),
    ADAS13 = paste0(
      paste0(sprintf("%.2f", mean(ADAS13, na.rm = T)), 
             " (", sprintf("%.2f", ADASsd = sd(ADAS13, na.rm = T)), "; "),
      paste(ADASrange = range(ADAS13, na.rm = T), collapse = " - "), ")"),
    MMSE = paste0(
      paste0(sprintf("%.2f", mean(MMSE, na.rm = T)),
             " (", sprintf("%.2f", MMSEsd = sd(MMSE, na.rm = T)), "; "),
      paste(MMSErange = range(MMSE, na.rm = T), collapse = " - "), ")"),
    FAQ = paste0(
      paste0(sprintf("%.2f", mean(FAQ, na.rm = T)),
             " (", sprintf("%.2f", FAQsd = sd(FAQ, na.rm = T)), "; "),
      paste(FAQrange = range(FAQ, na.rm = T), collapse = " - "), ")")
  ) %>% t() %>% as.data.frame()  %>%
  rename_at(c("V1", "V2"), ~ c("group 1", "group 2")) %>% 
  tibble::rownames_to_column("Characteristic") %>% 
  relocate("Characteristic", .before = "group 1") %>%
  mutate(
    P.value = c(NA, NA,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(Age, quartile) %>% 
          wilcox.test(Age ~ quartile, .)
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(Sex, quartile) %>% 
          table %>% as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(raceth, quartile) %>% 
          mutate(raceth = factor(raceth == "NH-White")) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(Education, quartile) %>% 
          wilcox.test(Education ~ quartile, .)
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(ApoE4, quartile) %>% 
          mutate(ApoE4 = as.numeric(ApoE4 > 0)) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% 
          select(antiDementiaDrugs, quartile) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(GxG, quartile) %>% 
          table %>% as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(CDGLOBAL, quartile) %>% 
          mutate(CDGLOBAL = factor(CDGLOBAL == 0)) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(CDGLOBAL, quartile) %>% 
          mutate(CDGLOBAL = factor(CDGLOBAL == 0.5)) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(CDGLOBAL, quartile) %>% 
          mutate(CDGLOBAL = factor(CDGLOBAL == 1)) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(DX, quartile) %>% 
          mutate(DX = factor(DX == "CN")) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(DX, quartile) %>% 
          mutate(DX = factor(DX == "MCI")) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(DX, quartile) %>% 
          mutate(DX = factor(DX == "Dementia")) %>% table %>% 
          as.data.frame.matrix %>% fisher.test()
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(CDRSB, quartile) %>% 
          wilcox.test(CDRSB ~ quartile, .)
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(ADAS13, quartile) %>% 
          wilcox.test(ADAS13 ~ quartile, .)
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(MMSE, quartile) %>% 
          wilcox.test(MMSE ~ quartile, .)
      }$p.value,
      {
        adni.trial %>% filter(VISCODE == 0) %>% select(FAQ, quartile) %>% 
          wilcox.test(FAQ ~ quartile, .)
      }$p.value
    )
  ) %>%
  mutate(P.value = ifelse(
    P.value < 0.01,
    sprintf("%.2e", P.value),
    sprintf("%.3f", P.value)
  ))

# Test outcomes
outcomes <- c("CDRSB", "ADAS13", "MMSE", "FAQ")
names(outcomes) <- outcomes

trial.outcomes <- outcomes %>% lapply(., function(y) {
  out <-
    adni.trial %>% filter(VISCODE == 0,!is.na(quartile)) %>% 
    group_by(quartile) %>% summarize(n = sum(!is.na(!!ensym(y))))
  ## Effect of bioenergetic age alone
  # Here, selection of relevant entries from EMM models is trivial: 
  # reference is first, outcome of interest is last entry
  x <-
    lmer(as.formula(
      paste0(
        y,
        " ~ ",
        y,
        ".base * VISCODE + DX.base + antiDementiaDrugs + Age * VISCODE + Education + Sex + ApoE4 + quartile * VISCODE + raceth + (1 | RID)"
      )
    ), adni.trial)
  EMMs <-
    pairs(emmeans(x, ~ VISCODE * quartile, at = list(VISCODE = c(0, 1.5))),
          by = NULL,
          adjust = NULL)
  out <-
    out %>% mutate(Adjusted.Mean.Change = sprintf("%.3f", summary(EMMs)$estimate[c(1, nrow(summary(EMMs)))] * (-1)))
  Diff <-
    contrast(EMMs,
             by = NULL,
             "trt.vs.ctrl",
             ref = 1,
             adjust = NULL)
  out <-
    out %>% mutate(Adjusted.Mean.Difference = c(
      sprintf("%.3f", confint(Diff)[nrow(summary(EMMs)) - 1, 2]),
      paste0(
        "(",
        sprintf("%.3f", confint(Diff)[nrow(summary(EMMs)) - 1, 5]),
        " to ",
        sprintf("%.3f", confint(Diff)[nrow(summary(EMMs)) - 1, 6]),
        ")"
      )
    )) %>%
    mutate(P.value = c(summary(Diff)[nrow(summary(EMMs)) - 1, 6], NA)) %>%
    mutate(P.value = ifelse(
      P.value < 0.01,
      sprintf("%.2e", P.value),
      sprintf("%.3f", P.value)
    )) %>%
    ## Effect of bioenergetic age in combination with GxG effect
    # Here, selection of relevant results from EMM models is not as trivial, 
    # so matching reference by contrast string, outcome of interest is still last
    mutate(n.GxG = {
      adni.trial %>% filter(VISCODE == 0,!is.na(quartile), GxG == "slower") %>% 
        group_by(quartile) %>% summarize(n.GxG = sum(!is.na(!!ensym(y)))) %>% 
        pull(n.GxG)
    })
  x2 <-
    lmer(as.formula(
      paste0(
        y,
        " ~ ",
        y,
        ".base * VISCODE + DX.base + antiDementiaDrugs + Age * VISCODE + Education + Sex + ApoE4 + quartile * VISCODE * GxG + raceth + (1 | RID)"
      )
    ), adni.trial)
  EMMs2 <-
    pairs(emmeans(x2, ~ VISCODE * quartile * GxG, at = list(VISCODE = c(0, 1.5))),
          by = NULL,
          adjust = NULL)
  idx <- which(summary(EMMs2)$contrast == "(VISCODE0 0-25% slower) - (VISCODE1.5 0-25% slower)")
  out <-
    out %>% mutate(Adjusted.Mean.Change.GxG = summary(EMMs2)$estimate[c(idx, nrow(summary(EMMs2)))] * (-1))
  Diff2 <-
    contrast(EMMs2,
             by = NULL,
             "trt.vs.ctrl",
             ref = idx,
             adjust = NULL)
  out <-
    out %>% mutate(Adjusted.Mean.Difference.GxG = c(
      sprintf("%.3f", confint(Diff2)[nrow(summary(EMMs2)) - 1, 2]),
      paste0(
        "(",
        sprintf("%.3f", confint(Diff2)[nrow(summary(EMMs2)) - 1, 5]),
        " to ",
        sprintf("%.3f", confint(Diff2)[nrow(summary(EMMs2)) - 1, 6]),
        ")"
      )
    )) %>%
    mutate(P.value.GxG = c(summary(Diff2)[nrow(summary(EMMs2)) - 1, 6], NA)) %>%
    mutate(P.value.GxG = ifelse(
      P.value.GxG < 0.01,
      sprintf("%.2e", P.value.GxG),
      sprintf("%.3f", P.value.GxG)
    ))
} %>% t)


################################################################################
#### Format and save outcomes #-------------------------------------------------
################################################################################

#### Table 1 -------------------------------------------------------------------
out <- trial.outcomes %>% names %>%
  lapply(., function(x) {
    trial.outcomes[[x]][3:5,] %>% 
      as.data.frame() %>% `colnames<-`(c("0-25%_All", "75-100%_All")) %>% 
      cbind(trial.outcomes[[x]][7:9,] %>% 
              as.data.frame() %>% 
              `colnames<-`(c("0-25%_slower", "75-100%_slower"))) %>%
      tibble::rownames_to_column("Outcome") %>% rbind(c(x, NA, NA, NA, NA), .)
  }) %>% 
  bind_rows()
openxlsx::write.xlsx(out, file = "results/Table 1.xlsx")

#### Supplementary Table 16 ----------------------------------------------------
out <-
  data.frame("Bioenergetic age (ADNI)", NA, NA, NA, NA, NA) %>% 
  rename_all( ~ colnames(statsBC)) %>% rbind(statsBC) %>%
  rbind(data.frame("Genotype (ADNI)", NA, NA, NA, NA, NA) %>% 
          rename_all( ~ colnames(statsBC))) %>%
  rbind(statsGxG) %>% rbind(
    data.frame("Genotype (ROS/MAP)", NA, NA, NA, NA, NA) %>% 
      rename_all( ~ colnames(statsBC))) %>%
  rbind(statsGxG_r) %>% 
  rbind(data.frame("Genotype/bioenergetic age interaction (ADNI)", NA, NA, NA, NA, NA) %>% 
          rename_all( ~ colnames(statsBC))
        ) %>% rbind(rbind(data.frame(
          t(colnames(statsGxGxE %>% rename_at(c("SE"),  ~ c("Std. Error")))
  )) %>% rename_all( ~ colnames(statsGxGxE)), statsGxGxE) %>% 
    rename_all( ~ c(colnames(statsBC))))
openxlsx::write.xlsx(out, file = "results/Supplementary Table 16.xlsx")

#### Supplementary Table 17 ----------------------------------------------------
colnames(tests) = colnames(tests) %>% gsub(
  x = .,
  pattern = "r...t...",
  replacement = "(t) ",
  fixed = T
) %>% gsub(
  x = .,
  pattern = "..",
  replacement = ".",
  fixed = T
)
openxlsx::write.xlsx(tests, file = "results/Supplementary Table 17.xlsx")

#### Supplementary Table 18 ----------------------------------------------------
openxlsx::write.xlsx(trial.characteristics, 
                     file = "results/Supplementary Table 18.xlsx")

#### Supplementary Figure 8 ----------------------------------------------------
p1 <- adni.long %>% filter(VISCODE == 0 & !is.na(GenoClust)) %>%
  group_by(GxG, GenoClust) %>% summarise(N = n()) %>% 
  mutate(GenoClust = as.character(GenoClust)) %>%
  as.data.frame %>% rbind(.,
                          data.frame(
                            GxG = "total",
                            GenoClust = c("faster", "slower"),
                            N = c(sum(.$N[.$GxG == "faster"]), sum(.$N[.$GxG == "slower"])),
                            stringsAsFactors = F
                          )) %>%
  ggplot(aes(x = GenoClust, y = N, fill = GxG)) +
  geom_bar(stat = 'identity', color = 'black') +
  facet_grid( ~ GxG, scales = "free", space = "free_x") +
  scale_y_continuous(labels = scales::comma_format(accuracy = 2)) +
  geom_text(aes(label = N), vjust = -0.25, size = rel(4.09)) +
  theme_bw() +
  labs(tag = "a") +
  xlab("alleles") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(color = 'black', size = 12),
    strip.text = element_text(size = 12),
    plot.tag = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) + ggsci::scale_fill_jama()

p2 <- ROSMAP.long %>% filter(fu_year == 0 & !is.na(GenoClust)) %>%
  group_by(GxG, GenoClust) %>% summarise(N = n()) %>% 
  mutate(GenoClust = as.character(GenoClust)) %>%
  as.data.frame %>% rbind(.,
                          data.frame(
                            GxG = "total",
                            GenoClust = c("faster", "slower"),
                            N = c(sum(.$N[.$GxG == "faster"]), sum(.$N[.$GxG == "slower"])),
                            stringsAsFactors = F
                          )) %>%
  ggplot(aes(x = GenoClust, y = N, fill = GxG)) +
  geom_bar(stat = 'identity', color = 'black') +
  facet_grid( ~ GxG, scales = "free", space = "free_x") +
  scale_y_continuous(labels = scales::comma_format(accuracy = 2)) +
  geom_text(aes(label = N), vjust = -0.25, size = rel(4.09)) +
  theme_bw() +
  labs(tag = "b") +
  xlab("alleles") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(color = 'black', size = 12),
    strip.text = element_text(size = 12),
    plot.tag = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) + ggsci::scale_color_jama() + ggsci::scale_fill_jama()

pdf("results/Supplementary Figure 8.pdf", width = 14, height = 7)
ggpubr::ggarrange(p1,p2,ncol = 2)
dev.off()

# clean up
rm(p1, p2)

#### Figure 4 ------------------------------------------------------------------
## Bioenergetic age - AGES
p1 <-
  ggboxplot(
    data = AGES.dt,
    x = "cogdx.fu",
    y = "BCAge.res",
    notch = T,
    fill = "cogdx.fu",
    add.params = list(
      alpha = .5,
      fill = "cogdx.fu",
      shape = 21
    ),
    add = "jitter",
    xlab = "Diagnosis at 5-year follow up",
    ylab = "Predicted bioenergetic Age\n(baseline)"
  ) +
  stat_compare_means(
    comparisons = list(c("CN", "MCI"), c("MCI", "AD"), c("CN", "AD")),
    method = "t.test",
    method.args = list(var.equal = T),
    symnum.args = list(
      cutpoints = c(0.00001, 0.001, 0.1, Inf),
      symbols = stats.ages
    )
  )  +
  labs(tag = "d") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(color = 'black', size = 12),
    axis.title.x = element_text(hjust = 1),
    plot.tag = element_text(size = 13, face = "bold")
  ) + ggsci::scale_fill_jama()

## Bioenergetic age - ADNI
plotsBC <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    # prepare plot for 50% split
    form <-
      paste0(y,
             " ~ percentile * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)")
    x <- lmer(as.formula(form), adni.long)
    tmpe <- effect(c("percentile*VISCODE"), x) %>% as.data.frame
    ggplot(tmpe, aes(x = VISCODE, y = fit, group = percentile)) + 
      geom_line(aes(color = percentile)) +
      geom_ribbon(aes(
        ymin = fit - se,
        ymax = fit + se,
        fill = percentile
      ), alpha = .2) +
      ylab(niceTrait) + xlab("Years from baseline") +
      annotate(
        geom = "text",
        size = rel(4.09),
        x = 2.5,
        y = max(tmpe$fit),
        label = paste0("P = ", sprintf("%.2g", statsBC[which(statsBC$trait == niceTrait), ncol(statsBC)]))
      ) +
      ggsci::scale_fill_jama(name = "Bioenergetic age [percentile]") +
      ggsci::scale_color_jama(name = "Bioenergetic age [percentile]") +
      theme_bw()
  })
# Polish plots
plotsBC <- 1:length(plotsBC) %>% lapply(., function(y){
  plotsBC[[y]] + labs(tag = letters[y]) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(color='black', size = 12),
          plot.tag = element_text(size = 13, face = "bold"))
})

## Allele combination effects - ADNI
plotsGxG <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    # prepare plot for 50% split
    form <-
      paste0(y,
             " ~ GxG * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)")
    x <- lmer(as.formula(form), adni.long)
    tmpe <- effect(c("GxG*VISCODE"), x) %>% as.data.frame
    ggplot(tmpe, aes(x = VISCODE, y = fit, group = GxG)) + 
    geom_line(aes(color = GxG)) + geom_ribbon(aes(ymin = fit - se,
                                                  ymax = fit + se,
                                                  fill = GxG), alpha = .2) +
      ylab(niceTrait) + xlab("Years from baseline") + annotate(
        geom = "text",
        size = rel(4.09),
        x = 2.5,
        y = max(tmpe$fit),
        label = paste0("P = ", sprintf("%.2g", statsGxG[which(statsGxG$trait == niceTrait), ncol(statsGxG)]))
      ) +
      #ggsci::scale_fill_jama() + ggsci::scale_color_jama()
      scale_fill_manual(values = ggsci::pal_jama()(7)[c(6, 7)], name = "Genotype") +
      scale_color_manual(values = ggsci::pal_jama()(7)[c(6, 7)], name = "Genotype") +
      theme_bw()
  })
# Polish plots
plotsGxG <- 1:length(plotsGxG) %>% lapply(., function(y) {
  plotsGxG[[y]] + labs(tag = letters[y]) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(color = 'black', size = 12),
      plot.tag = element_text(size = 13, face = "bold")
    )
})

## Allele combination effects - ROS/MAP
p2 <-
  lmer(
    cogn_global ~ GxG * fu_year + Sex + ApoE4 + Education + age_bl + dx_bl + (1 | projid),
    ROSMAP.long
  ) %>%
  effect(c("GxG*fu_year"), .) %>% as.data.frame %>%
  ggplot(., aes(x = fu_year, y = fit, group = GxG)) +
  geom_line(aes(color = GxG)) +
  geom_ribbon(aes(
    ymin = fit - se,
    ymax = fit + se,
    fill = GxG
  ), alpha = .2) +
  annotate(
    geom = "text",
    x = 6.5,
    y = 0.134144,
    size = rel(4.09),
    label = paste0("P = ", sprintf("%.2g", statsGxG_r[ncol(statsGxG_r)]))
  ) +
  ylab("Global cognition") + xlab("Years from baseline") + labs(tag = "d") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(color = 'black', size = 12),
    plot.tag = element_text(size = 13, face = "bold")
  ) +
  scale_fill_manual(values = ggsci::pal_jama()(7)[c(6, 7)], name = "Genotype") +
  scale_color_manual(values = ggsci::pal_jama()(7)[c(6, 7)], name = "Genotype") +
  theme_bw()

## Interactions of allele combination with bioenergetic age - ADNI
plotsGxGxE <-
  c("ADAS13", "ADNI_MEM", "ADNI_EF") %>% lapply(., function(y) {
    niceTrait <-
      ifelse(y == "ADAS13",
             "ADAS-Cog. 13",
             ifelse(y == "ADNI_MEM", "Memory", "Executive Function"))
    # prepare plot for 50% split
    form <-
      paste0(y,
             " ~ percentile * GxG * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)")
    x <- lmer(as.formula(form), adni.long)
    tmpe <- effect(c("percentile*GxG*VISCODE"), x) %>% as.data.frame
    ggplot(tmpe, aes(
      x = VISCODE,
      y = fit,
      group = interaction(GxG, percentile)
    )) + geom_line(aes(color = interaction(GxG, percentile))) + geom_ribbon(aes(
      ymin = fit - se,
      ymax = fit + se,
      fill = interaction(GxG, percentile)
    ), alpha = .2) +
      ylab(niceTrait) + xlab("Years from baseline") + annotate(
        geom = "text",
        size = rel(4.09),
        x = 2.5,
        y = max(tmpe$fit),
        label = paste0("P = ", sprintf("%.2g", statsGxGxE[which(statsGxGxE$trait == niceTrait), ncol(statsGxGxE)]))
      ) +
      scale_fill_manual(
        values = ggsci::pal_jama()(5)[c(1, 5, 3, 4)],
        name = "Genotype/bioenergetic age interaction (genotype : bioenergetic age percentile)",
        labels = c(
          "faster : 0-50%",
          "slower : 0-50%",
          "faster : 50-100%",
          "slower : 50-100%"
        )
      ) +
      scale_color_manual(
        values = ggsci::pal_jama()(5)[c(1, 5, 3, 4)],
        name = "Genotype/bioenergetic age interaction (genotype : bioenergetic age percentile)",
        labels = c(
          "faster : 0-50%",
          "slower : 0-50%",
          "faster : 50-100%",
          "slower : 50-100%"
        )
      ) +
      theme_bw()
  })
# Polish plots
plotsGxGxE <- 1:length(plotsGxGxE) %>% lapply(., function(y){
  plotsGxGxE[[y]] + labs(tag = letters[y]) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(color='black', size = 12),
          plot.tag = element_text(size = 13, face = "bold"))
})

x1 <- ggpubr::ggarrange(
  plotsBC[[1]],
  plotsBC[[2]],
  plotsBC[[3]],
  p1 + coord_cartesian(clip = "off"),
  common.legend = T,
  ncol = 4,
  legend = "bottom"
)
x1 <- annotate_figure(x1, top = text_grob("Bioenergetic age",
                                          color = "black", size = 14))
x2 <- ggpubr::ggarrange(
  plotsGxG[[1]] + labs(tag = "e"),
  plotsGxG[[2]] + labs(tag = "f"),
  plotsGxG[[3]] + labs(tag = "g"),
  p2 + labs(tag = "h") + theme(plot.tag = element_text(face = "bold")),
  common.legend = T,
  ncol = 4,
  legend = "bottom"
)
x2 <- annotate_figure(x2, top = text_grob("Genotype",
                                          color = "black", size = 14))
x3 <- ggpubr::ggarrange(
  plotsGxGxE[[1]] + labs(tag = "i"),
  plotsGxGxE[[2]] + labs(tag = "j"),
  plotsGxGxE[[3]] + labs(tag = "k"),
  common.legend = T,
  ncol = 3,
  legend = "bottom"
)
x3 <-
  annotate_figure(x3,
                  top = text_grob(
                    "Genotype/bioenergetic age interaction",
                    color = "black",
                    size = 14
                  ))

pdf("results/Figure4.pdf", width = 11, height = 11)
grid.arrange(x1, x2, x3, layout_matrix = matrix(
  c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3),
  ncol = 1,
  byrow = T
))
grid.rect(
  x = 0.5,
  y = 8 / 44,
  height = 8 / 22,
  width = 1,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  )
)
grid.rect(
  x = 0.5,
  y = 23 / 44,
  height = 7 / 22,
  width = 1,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  )
)
grid.rect(
  x = 0.5,
  y = 37 / 44,
  height = 7 / 22,
  width = 1,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  )
)
dev.off()
