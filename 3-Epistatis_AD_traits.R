##### Bioenergetic capacity - Script 3 #####
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Nov 29, 2024              ##
#                                         ##
############################################

#### Front matter #-------------------------------------------------------------
library(magrittr)
library(dplyr)
library(openxlsx)
library(survival)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(ggpubr)
library(grid)
library(gridExtra)

# clear environment
rm(list = ls())

# load helper functions for multi-SNP association tests and plotting
source("helper_scripts/Epistasis_Helper_Functions.R")
source("helper_scripts/Epistasis_Plotting_Functions.R")

#### load required data #-------------------------------------------------------
load("data/Script-2.ADNI.genetics.RData") # Contents described in Script 2
load("data/Script-3.ADNI.Mayo.ROSMAP.geno-pheno.RData") # contains three data frames:
# -----------
# ADNI.gp.dt: ADNI phenotype data for all individuals in ADNI with genetics data
# -----------
# df with dimensions 1548 (samples, identified by column RID) x 30 
# (covariates & ATN(C)/diagnosis readouts) -> row order indentical with ADNI.dt
#
# ----------
# ROSMAP.dt: Genotype and phenotype data for ROS/MAP 
# ----------
# df with dimensions 2059 (samples, identified by column projid) x 47 
# (SNPs, covariates & AT(C)/diagnosis readouts)
#
# --------
# Mayo.dt: Genotype and phenotype data for MayoLOAD
# --------
# df with dimensions 2067 (samples, identified by column PTID) x 43 
# (SNPs, covariates & diagnosis readouts)

#### Derive data/formats required by this script #------------------------------
## Global ----
# Extract SNPs to be tested for epistasis
snps <- snpAnno %>% filter(interaction.model == "yes") %>% pull(SNP)
# Create list of unique single SNPs/combinations for haplotype analysis
haps <- ac.haplotypes %>% pull(metabolite) %>% as.list() %>% 
  `names<-`(c(.)) %>% 
  lapply(function(x){
    snps <- ac.haplotypes %>% filter(metabolite == x) %>% pull(snps) %>% 
      strsplit(split = ", ", fixed = T) %>% unlist 
    do.call("c", lapply(seq_along(snps), function(i) combn(snps, i, FUN = list))) %>% 
      lapply(paste, collapse = "_") %>% `names<-`(unlist(.)) %>% 
      lapply(function(y){
        strsplit(y, split = "_", fixed = T) %>% unlist
      })
  })

## ADNI ------
# Merge ADNI genetics data with phenotypic data
ADNI.dt <- cbind(ADNI.gp.dt, ADNI.dt)
# List of traits to be tested
ADNI.tnames <-
  c("ABETA142", "PTAU", "FDG.bl", "ADAS.Cog13", "Diagnosis")
# List of transformations to be applied
ttrans <-
  c(
    "scale(ABETA142)",
    "scale(log2(PTAU))",
    "scale(FDG.bl)",
    "scale(sqrt(ADAS.Cog13))",
    "Diagnosis"
  )
# List of covariates
ADNI.covar <-
  c("Education", "Cohort", "ApoE4", "Age", "Sex", "bmi_in_kg_p_m2")
# select relevant data
ADNI.dt %<>% select(all_of(c(ADNI.tnames, ADNI.covar, snps)))
ADNI.covar <- ADNI.covar[-c(1:2)]
ADNI.covar[1] <- "factor(ApoE4)"
# build trait map
ADNI.traits <- data.frame(
  trait = ttrans,
  type = "linear_regression",
  covars = c(NA, "Education", "Education + Cohort", rep("Education", 2)),
  stringsAsFactors = F
)

# ROS/MAP ----
ROSMAP.tnames <-
  c("Global_Cognition",
    "Amyloid",
    "Tangles",
    "Gpath",
    "Clinical_Diagnosis")
ROSMAP.covar <- c("PMI", "ApoE4", "Age", "Sex", "Education")
ROSMAP.dt %<>% select(all_of(c(
  ROSMAP.tnames, ROSMAP.covar, "Education", snps
)))
ROSMAP.covar <- ROSMAP.covar[-1]
ROSMAP.covar[1] <- "factor(ApoE4)"

# build trait map
ROSMAP.traits <- data.frame(
  trait = ROSMAP.tnames,
  type = "linear_regression",
  covars = NA,
  stringsAsFactors = F
)
ROSMAP.traits$covars[ROSMAP.traits$trait %in% c("Amyloid", "Tangles", "Gpath")] <-
  "PMI"

# MayoLOAD ---
# List of traits to be tested
Mayo.tnames <- c("numDiagnosis", "Dementia")
# List of covariates
Mayo.covar <- c("ApoE4", "Age", "Sex")
# select relevant data
Mayo.dt %<>% select(all_of(c(Mayo.tnames, Mayo.covar, snps)))
Mayo.covar[1] <- "factor(ApoE4)"

# build trait map
Mayo.traits <- data.frame(
  trait = Mayo.tnames,
  type = c("linear_regression", "logistic_regression"),
  covars = NA,
  stringsAsFactors = F
)

## Cleanup
rm(list = c("ttrans", "snps", "snpAnno", "ADNI.gp.dt"))

################################################################################
#### Perform epistatic modeling #-----------------------------------------------
################################################################################

#### ADNI ####
## lm-glm
# run the analysis for lm-glm 
re <-
  epistasisGLM(
    dt = ADNI.dt,
    traits = ADNI.traits,
    covar = ADNI.covar,
    tnames = ADNI.tnames,
    haps = haps
  )
ADNI.lm <- epistasisANOVA(re, logtransform = T)

## ranking based model
# fit ranking model for all of the variables except the binary ones for which 
# glm-binomial is fitted
re0 <-
  epistasisRanking(
    dt = ADNI.dt,
    traits = ADNI.traits,
    covar = ADNI.covar,
    tnames = ADNI.tnames,
    haps = haps,
    re = re
  )
ADNI.ranking <- epistasisANOVA(re0, logtransform =T)
# Store for plotting
ADNI.re0 <- re0

## Cleanup
rm(re, re0) 
gc()

#### ROS/MAP ####
## lm-glm
# run the analysis for lm-glm 
re <-
  epistasisGLM(
    dt = ROSMAP.dt,
    traits = ROSMAP.traits,
    covar = ROSMAP.covar,
    tnames = ROSMAP.tnames,
    haps = haps
  )
ROSMAP.lm <- epistasisANOVA(re, logtransform = T)


## ranking based model
# fit ranking model for all of the variables except the binary ones for which 
# glm-binomial is fitted
re0 <-
  epistasisRanking(
    dt = ROSMAP.dt,
    traits = ROSMAP.traits,
    covar = ROSMAP.covar,
    tnames = ROSMAP.tnames,
    haps = haps,
    re = re
  )
ROSMAP.ranking <- epistasisANOVA(re0, logtransform =T)

# Cleanup
rm(re, re0) 
gc()

#### Mayo ####
## lm-glm
# run the analysis for lm-glm 
re <-
  epistasisGLM(
    dt = Mayo.dt,
    traits = Mayo.traits,
    covar = Mayo.covar,
    tnames = Mayo.tnames,
    haps = haps
  )
Mayo.lm <- epistasisANOVA(re, logtransform = T)


## ranking based model
# fit ranking model for all of the variables except the binary ones for which 
# glm-binomial is fitted
re0 <-
  epistasisRanking(
    dt = Mayo.dt,
    traits = Mayo.traits,
    covar = Mayo.covar,
    tnames = Mayo.tnames,
    haps = haps,
    re = re
  )
Mayo.ranking <- epistasisANOVA(re0, logtransform =T)

# Cleanup
rm(re, re0) 
gc()


################################################################################
#### Format and save outcomes #-------------------------------------------------
################################################################################

#### Supplementary table 10 ----------------------------------------------------
out <- ADNI.ranking %>% names %>% lapply(., function(x) {
  ADNI.ranking[[x]] %>% as.data.frame %>% 
    tibble::rownames_to_column(var = "SNPs") %>% mutate(SNPs = gsub(
      x = SNPs,
      pattern = "_",
      replacement = " x "
      )) %>%  mutate("Genetic Model" = x)
  }) %>% bind_rows() %>% 
  mutate(`Genetic Model` = factor(.$`Genetic Model`, 
                                  levels = ac.haplotypes$metabolite[1:nrow(ac.haplotypes)])) %>%
  cbind(
    .,
    ROSMAP.ranking %>% lapply(., as.data.frame) %>% bind_rows,
    Mayo.ranking %>% lapply(., as.data.frame) %>% bind_rows
  ) %>%
  group_by(`Genetic Model`) %>% arrange(`Genetic Model`) %>% as.data.frame %>% 
  relocate("Genetic Model", .before = "SNPs") %>%
  rename_at(
    colnames(.)[-c(1:2)],
    ~ c(
      "CSF Abeta1-42",
      "CSF p-tau",
      "FDG-PET",
      "ADAS-Cog. 13",
      "Diagnosis [ADNI]",
      "Global cognition",
      "Amyloid load",
      "PHF tangle load",
      "Global pathology",
      "Diagnosis [ROS/MAP]",
      "Diagnosis [staged]",
      "Diagnosis [binary]"
    )
  ) %>%
  group_by(`Genetic Model`) %>% 
  mutate("Genetic Model" = c(`Genetic Model` %>% 
                               as.character %>% 
                               head(1), rep(NA, length(`Genetic Model`) - 1))) %>% 
  as.data.frame(row.names = NULL) %>% 
  mutate_at(colnames(.)[-c(1:2)], function(x) { 10 ^ x })

write.xlsx(out, file = "results/Supplementary Table 10.xlsx")

#### Supplementary Figure 3 ----------------------------------------------------
# results in a nutshell as heatmaps 
# create plot matrix for ADNI
plot.adni <- ADNI.ranking %>%
  {
    lapply(structure(names(.), names = names(.)), function(i) {
      m1 = apply(.[[i]], 2, min) < log10(0.05 / ncol(.[[i]]))
      m2 = apply(.[[i]], 2, function(x)
        grep(names(which.min(x)), pattern = "_") %>% length)
      m3 = -1 * (apply(.[[i]], 2, min))
      m4 <- t(m1 * m2 * m3)
      rownames(m4) <- i
      as.data.frame(m4)
    })
  } %>% bind_rows()

# create plot matrix for ROS/MAP
plot.rosmap <- ROSMAP.ranking %>%
  {
    lapply(structure(names(.), names = names(.)), function(i) {
      m1 = apply(.[[i]], 2, min) < log10(0.05 / ncol(.[[i]]))
      m2 = apply(.[[i]], 2, function(x)
        grep(names(which.min(x)), pattern = "_") %>% length)
      m3 = -1 * (apply(.[[i]], 2, min))
      m4 <- t(m1 * m2 * m3)
      rownames(m4) <- i
      as.data.frame(m4)
    })
  } %>% bind_rows()

# create plot matrix for Mayo
plot.mayo <- Mayo.ranking %>%
  {
    lapply(structure(names(.), names = names(.)), function(i) {
      m1 = apply(.[[i]], 2, min) < log10(0.05 / ncol(.[[i]]))
      m2 = apply(.[[i]], 2, function(x)
        grep(names(which.min(x)), pattern = "_") %>% length)
      m3 = -1 * (apply(.[[i]], 2, min))
      m4 <- t(m1 * m2 * m3)
      rownames(m4) <- i
      as.data.frame(m4)
    })
  } %>% bind_rows()

# combined df
plot.dat <- t(cbind(plot.adni, plot.rosmap, plot.mayo))

### Plot Heatmap
sfig3a <- epiHeatmap(plot.dat)

### Phenotype scatter plots
man.plot <- ADNI.ranking %>% names %>%
  lapply(., function(x) {
    ADNI.ranking[[x]] %>%
      abs %>% as.data.frame %>%
      tibble::rownames_to_column(var = "rwnames") %>% rowwise %>%
      mutate(combo = (strsplit(rwnames, split = "_") %>% unlist %>%
                        length > 1)) %>% as.data.frame %>%
      mutate(ac = x) %>% mutate(label = paste0(ac, "\n",
                                               gsub(
                                                 rwnames,
                                                 pattern = "_",
                                                 replacement = " x "
                                               ), "")) %>%
      select(-rwnames) %>% group_by(ac) %>%
      mutate(ab.max = ABETA142 == max(.$ABETA142)) %>%
      mutate(tau.max = PTAU == max(.$PTAU)) %>% 
      mutate(fdg.max = FDG.bl == max(.$FDG.bl)) %>%
      mutate(adas.max = ADAS.Cog13 == max(.$ADAS.Cog13)) %>%
      mutate(dx.max = Diagnosis == max(.$Diagnosis)) %>%
      as.data.frame}) %>% bind_rows() %>% 
  mutate(ac = factor(.$ac, levels = ac.haplotypes$metabolite[1:nrow(ac.haplotypes)])) %>%
  group_by(ac) %>% arrange(ac) %>% as.data.frame %>% 
  mutate(xcord = 1:nrow(.)) %>% group_by(ac) %>% 
  mutate(xcord = xcord + (5 * (which(levels(ac) == ac[1]) -1 ))) %>%
  mutate(xticks = mean(xcord)) %>% as.data.frame

sfig3b <-
  starManhattanGeno(data = man.plot,
                trait = "ABETA142",
                traitmax = "ab.max") +
  labs(tag = "b") + theme(plot.tag = element_text(face = "bold"))
sfig3c <-
  starManhattanGeno(data = man.plot,
                trait = "PTAU",
                traitmax = "tau.max") +
  labs(tag = "c") + theme(plot.tag = element_text(face = "bold"))
sfig3d <-
  starManhattanGeno(data = man.plot,
                trait = "FDG.bl",
                traitmax = "fdg.max") +
  labs(tag = "d") + theme(plot.tag = element_text(face = "bold"))
sfig3e <-
  starManhattanGeno(data = man.plot,
                trait = "ADAS.Cog13",
                traitmax = "adas.max") +
  labs(tag = "e") + theme(plot.tag = element_text(face = "bold"))
sfig3f <-
  starManhattanGeno(data = man.plot,
                trait = "Diagnosis",
                traitmax = "dx.max") + 
  labs(tag = "f") + theme(plot.tag = element_text(face = "bold"))


pdf("results/Supplementary Figure 3.pdf", width = 11, height = 14)
grid.arrange(sfig3a, 
             sfig3b, 
             sfig3c, 
             sfig3d, 
             sfig3e, 
             sfig3f, 
             layout_matrix = matrix(c(1,2,3,4,5,6), ncol=2, byrow = T))
dev.off()
