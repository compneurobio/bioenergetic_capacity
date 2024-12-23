##### Bioenergetic capacity - Script 4 #####
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Dec 18, 2024              ##
#                                         ##
############################################

#### Front matter #-------------------------------------------------------------
library(openxlsx)
library(ggplot2)
library(caret)
library(ggpubr)
library(grid)
library(gridExtra)
library(dplyr)
library(magrittr)
library(lmerTest)
library(effects)
library(ggfortify)
library(qqplotr)

# clear environment
rm(list=ls())

#### load required data #-------------------------------------------------------
load("data/Script-1.ADNI.p180.RData") # Contents described in Script 1
load("data/store/SGI_AC_clustering.RData") # SGI clustering results from Script 1
load("data/Script-4.KORA.AGES.AC.RData") # contains two data frames:
# --------
# KORA.dt: preprocessed acylcarnitine (AC) profiles and covariates in KORA
# --------
# df with dimensions 3029 (samples, identified by sample.id) x 27
# (IDs, basic covariates, + 22 ACs)
#
# --------
# AGES.dt: preprocessed acylcarnitine (AC) profiles, phenotypes and covariates in AGES
# --------
# df with dimensions 575 (samples, identified by sample.id) x 61
# (IDs, covariates, phenotypes, + ACs)

#### Derive data/formats required by this script #------------------------------
## ADNI ------
# Pull together AC levels and phenotypes
adni.data <- data.frame(RID = rownames(p180) %>% as.numeric(), 
                        cbind(p180, p180.pheno), check.names = F)

# Regress out cohort effects from acylcarnitine profiles
p180.adj <- p180[, -1] %>% apply(2, function(x){ 
  lm(x ~ p180$Cohort)$residuals %>% scale 
}) %>% as.data.frame %>% `rownames<-`(rownames(p180))

# Limit list of ACs to those available in all cohorts (KORA is the limiting factor here)
metnames <- colnames(p180.adj)[colnames(p180.adj) %in% colnames(KORA.dt)]

# Cleanup - keep only data that we need here
rm(list = setdiff(
  ls(),
  c(
    "adni.data",
    "AGES.dt",
    "KORA.dt",
    "metpath",
    "SGI_clusters",
    "metnames",
    "p180.adj"
  )
))


################################################################################
#### Select reference individuals ----------------------------------------------
################################################################################

### Retrieve reference individuals while masking cohorst-specific IDs
## Thresholds based on representative ranges in healthy females for  
## age (ADNI + AGES) and BMI (ADNI only) while approximately balancing 
## sample size across the two cohorts

# KORA
ref.kora <- KORA.dt %>% filter(Sex == 2 &
                                 Age > 59.9 & Age < 72 &
                                 BMI <= 29.72 &
                                 BMI >= 23.31) %>%
  rename_at("sample.id", ~ "RID") %>% select(all_of(c("RID", metnames))) %>%
  mutate(cohort = "KORA") %>% mutate(uID = paste0(cohort, ".", 1:nrow(.)))
ref.kora[, metnames] <- apply(ref.kora[, metnames], 2, scale)

# ADNI
ref.adni <-
  # Age in ADNI is quite accurate, so be a bit more restrictive "<"
  adni.data %>% filter(Diagnosis.raw == "CN" &
                         Sex == 2 & Age < 72) %>%
  select(all_of(c("RID", metnames))) %>%
  mutate(cohort = "ADNI") %>% mutate(uID = paste0(cohort, ".", 1:nrow(.)))
ref.adni[, metnames] <- apply(ref.adni[, metnames], 2, scale)

# AGES
ref.ages <-
  # Age in AGES is rounded to years, so be a bit less restrictive "<="
  AGES.dt %>% filter(COGADNEW == 0, Sex == 2 &
                       Age <= 72) %>% rename_at("id", ~ "RID") %>%
  select(all_of(c("RID", metnames)))  %>% mutate(cohort = "AGES") %>%
  mutate(uID = paste0(cohort, ".", 1:nrow(.)))
ref.ages[, metnames] <- apply(ref.ages[, metnames], 2, scale)

reference <- ref.kora %>% rbind(., ref.adni) %>% rbind(., ref.ages)

# clean up
rm(list = c("ref.kora", "ref.adni", "ref.ages"))

### Identify multivariable outliers in reference individuals
# Using Mahalanobis distance
mahal.d <-
  mahalanobis(reference[, metnames], center = F,
              cov = cov(reference[, metnames]))
nrid <- nrow(reference)
nmets <- length(metnames)
o <- which(mahal.d > qchisq(1 - 0.05 / nrid, nmets))

# Generate QQ plot
plot.dat <- data.frame(
  MD = mahal.d,
  outlier = factor(ifelse(
    mahal.d > qchisq(1 - 0.05 / nrid, nmets), "yes", "no"
  ), levels = c("yes", "no")),
  label = reference$uID
)
plot.dat$label[-o] <- ""
plot.dat <- plot.dat[order(plot.dat$MD),]
plot.dat$qt = c(qchisq((1:(nrid - 3)) / (nrid - 2), nmets),
                qchisq(((nrid - 2):nrid) / (nrid + 1), nmets))

p1 <- ggplot(plot.dat, aes(sample = MD)) +
  geom_qq_band(
    distribution = "chisq",
    dparams = list(df = nmets),
    bandType = "pointwise",
    conf = 1 - 0.05 / (nrid - 3),
    alpha = 0.25
  ) +
  stat_qq_line(distribution = "chisq", dparams = list(df = nmets)) +
  geom_point(aes(x = qt, y = MD, color = outlier)) +
  geom_text(aes(x = qt, y = MD, label = label), hjust = 1.1) +
  xlab(bquote({
    chi ^ 2
  }[p - quantiles])) + ylab(bquote({
    italic(D) ^ 2
  }[Mahalanobis])) + ggsci::scale_color_jama()

# Remove outliers
reference <- reference[, -o]

# Generate PCA plot (dim 1 vs. 2)
pca.ref <- prcomp(reference[, metnames], scale. = T, center = T)
p2 <- autoplot(pca.ref, data = reference, colour = "cohort") +
  ggsci::scale_color_jama()

# clean up
rm(list = c("plot.dat", "nmets", "nrid", "mahal.d", "o", "pca.ref"))


################################################################################
#### Learn age predictor in KORA -----------------------------------------------
################################################################################
### Do reference-based rescale
ref.subj <- which(KORA.dt$sample.id %in% reference$RID[reference$cohort == "KORA"])
KORA.dt[,metnames] <- apply(KORA.dt[,metnames], 2, function(x){
  tr <- scale(x[ref.subj])
  x * attr(tr, 'scaled:scale') + attr(tr, 'scaled:center')
})

# calculate global model (lm is cross-validated, so we will use this later)
KORA.dt$Age.scaled <- scale(KORA.dt$Age)
form <- paste0("Age.scaled ~  ", paste(metnames, collapse = " + "))
global.model <- lm(as.formula(form), KORA.dt) # use biglm to share

## Perform repeated 10-fold CV on original age scale to report on model stability
# setting seed to generate a reproducible random sampling
set.seed(1234)

# defining training control as repeated cross-validation with
# value of K = 10 and repeats = 3 times
train_control <- trainControl(method = "repeatedcv",
                              number = 10, repeats = 3)

form <- paste0("Age ~ ", paste(metnames, collapse = " + "))
CVmodel <- train(as.formula(form), data = KORA.dt,
                 method = "lm",
                 trControl = train_control)

KORA.dt$BCAge <- predict(CVmodel$finalModel, newdata = KORA.dt)

cvresults <-
  data.frame(
    CVmodel$resample,
    N_train = CVmodel$control$index %>% lapply(length) %>% as.numeric(),
    N_control = (
      dim(KORA.dt)[1] - CVmodel$control$index %>% lapply(length) %>% as.numeric()
    ),
    stringsAsFactors = F
  )
cvresults <-
  rbind(cvresults[, c(4, 1:3, 5:6)], c("Average", NA, NA, NA, NA, NA))
cvresults <-
  rbind(cvresults, names(CVmodel$results[, -1])[c(1, 4, 2, 5, 3, 6)])
cvresults <-
  rbind(cvresults, as.character(CVmodel$results[, -1])[c(1, 4, 2, 5, 3, 6)])

# clean up
rm(list = c("train_control", "CVmodel", "ref.subj", "form"))


################################################################################
### APPLY TO ADNI --------------------------------------------------------------
################################################################################

adni.data[, metnames] <- apply(adni.data[, metnames], 2, scale)
### Do reference-based rescale
ref.subj.adni <-
  which(adni.data$RID %in% reference$RID[reference$cohort == "ADNI"])
adni.data[, metnames] <- apply(adni.data[, metnames], 2, function(x) {
  tr <- scale(x[ref.subj.adni])
  x * attr(tr, 'scaled:scale') + attr(tr, 'scaled:center')
})
tr <- scale(adni.data$Age)

BCAge <- predict(global.model, newdata = adni.data)
BCAge <-
  BCAge * attr(tr, 'scaled:scale') + attr(tr, 'scaled:center')
adni.data$BCAge <- BCAge

adni.data <- merge(adni.data, SGI_clusters, by = "RID")
adni.data %<>% mutate(AgeDiff = as.numeric(scale(scale(BCAge) - scale(Age))))

# clean up
rm(list = c("tr", "ref.subj.adni", "BCAge"))


################################################################################
### APPLY TO AGES --------------------------------------------------------------
################################################################################

AGES.dt[,metnames] <- apply(AGES.dt[,metnames],2,scale)

#### Do the same demographics-based rescale ####
ref.subj.ages <- which(AGES.dt$id %in% reference$RID[reference$cohort=="AGES"])
AGES.dt[,metnames] <- apply(AGES.dt[,metnames], 2, function(x){
  tr <- scale(x[ref.subj.ages])
  x * attr(tr, 'scaled:scale') + attr(tr, 'scaled:center')
})
tr <- scale(AGES.dt$Age)

BCAge <- predict(global.model, newdata = AGES.dt)
BCAge <- BCAge * attr(tr, 'scaled:scale') + attr(tr, 'scaled:center')
AGES.dt$BCAge <- BCAge

AGES.dt %<>% mutate(AgeDiff = as.numeric(scale(scale(BCAge) - scale(Age))))

# clean up
rm(list = c("tr", "ref.subj.ages", "BCAge"))


################################################################################
#### Investigate bioenergetic age (cross-sectional) #---------------------------
################################################################################

#### ADNI ----------------------------------------------------------------------
## Define diagnostic groups
DXgroups <- list(All = paste0("Diagnosis.raw %in% c(\"",
                              paste(levels(adni.data$Diagnosis.raw), 
                                    collapse = "\", \""),"\")"),
                 CN = "Diagnosis.raw %in% c(\"CN\", \"SMC\")",
                 MCI = "Diagnosis.raw %in% c(\"EMCI\", \"LMCI\")",
                 AD = "Diagnosis.raw %in% c(\"AD\")",
                 "cluster 7" = "l4 == 7")

## Correlation analysis
adni.age.cor <- DXgroups %>% names %>% head(-1) %>% lapply(. , function(x){
  summary(lm(scale(Age) ~ scale(BCAge), 
             data = adni.data %>% 
               filter(!! rlang::parse_expr(DXgroups[[x]]))
  ))$coefficients[2, ] %>% 
    t %>% cbind(data.frame(Cohort = "ADNI", Group = x), .) %>%
    rename_at("Estimate", ~"Rho")
}) %>% bind_rows

## Trait associations
traits <- c("CSF Ab1-42", "CSF p-tau", "FDG-PET", "ADAS-Cog. 13")

out.adni <- DXgroups %>% names %>% lapply(. , function(x) {
  traits %>% lapply(function(trait) {
    data.frame(sample = x, trait = trait) %>%
      cbind(.,
            summary(lm(
              as.formula(
                paste0(
                  "scale(`",
                  trait,
                  "`) ~ BCAge + Age + Sex + factor(`APOE e4`) + Education + BMI"
                )
              ),
              data = adni.data %>% filter(!!rlang::parse_expr(DXgroups[[x]]))
            ))$coefficients[2, ] %>% t)
  })
}) %>% bind_rows

## Prepare plot for differences in age measures across diagnostic groups
# get color code
fill.c <- ggsci::pal_jama()(3)

plot.dat <- adni.data
plotlist <- list()
traits <- list(Age = "Chronologcial age", 
               BCAge = "Predicted bioenergetic age", 
               AgeDiff = expression(Delta('Bioenergetic - chronological age')))
for(bp in c("l2", "l4", "l5")) {
  for(trait in names(traits)){
    plot.dat[is.na(plot.dat[, bp]), bp] <- "other"
    plot.dat[, bp] <- factor(plot.dat[, bp])
    
    lvls <- levels(plot.dat[, bp])
    lvlcombs <- combn(lvls, 2)
    lvlslist <- list()
    for (i in 1:ncol(lvlcombs)) {
      lvlslist[[length(lvlslist) + 1]] <- lvlcombs[, i]
    }
    taglist <- list(l2 = "c", l4 = "d", l5 = "e")
    
    pvals <- lvlslist %>% lapply(function(x){ 
      t.test(as.formula(paste0(trait, " ~ ", bp)), 
             data = plot.dat[plot.dat[,bp] %in% x,],
             var.equal = T)$p.value
    }) %>% unlist() %>% sort()
    
    bpt <-
      ggboxplot(
        data = plot.dat,
        x = bp,
        y = trait,
        add = "jitter",
        notch = T,
        fill = bp,
        add.params = list(
          alpha = .5,
          shape = 21,
          color = "black",
          fill = bp
        )
      ) +
      scale_fill_manual(values = fill.c[1:length(lvls)]) +
      stat_compare_means(
        comparisons = lvlslist,
        label = "p.format",
        method = "t.test",
        method.args = list(var.equal = T),
        label.x = 1.35,
        size = rel(4.09),
        symnum.args = list(cutpoints = c(0, pvals + 1e-11),
                           symbols = sapply(pvals, function(p) {
                             ifelse(p < 1e-4, 
                                    sprintf("%.1e", p), 
                                    sprintf("%.2g", p))
                           }))
      ) +
      xlab("Cluster") + ylab(traits[[trait]]) +
      theme(
        legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.margin = unit(c(.05, .03, .03, .03), units = "npc"),
        plot.tag.position = c(.01, 1)
      ) + coord_cartesian(clip = "off")
    
    if(trait == "Age"){
      bpt <- bpt  + labs(tag = taglist[[bp]])
    }
    
    if(trait == "AgeDiff"){
      bpt <- bpt +
        geom_hline(
          yintercept = 0,
          color = "black",
          linetype = "dashed",
          alpha = .5
        )
    }
    
    plotlist[[length(plotlist) + 1]] <- bpt
  }
}

## Calculate differences in age measures across diagnostic groups
out.clusters <- plot.dat %>% select(l2, l4, l5) %>% lapply(function(x) {
  lvls <- levels(x)
  lvlcombs <- combn(lvls, 2)
  lvlslist <- list()
  for (i in 1:ncol(lvlcombs)) {
    lvlslist[[length(lvlslist) + 1]] <- lvlcombs[, i]
  }
  lvlslist %>% lapply(function(y) {
    z1 <-
      t.test(plot.dat$Age[x == y[[1]]], plot.dat$Age[x == y[[2]]], var.equal = T)
    z2 <-
      t.test(plot.dat$BCAge[x == y[[1]]], plot.dat$BCAge[x == y[[2]]], var.equal = T)
    z3 <-
      t.test(plot.dat$AgeDiff[x == y[[1]]], plot.dat$AgeDiff[x == y[[2]]], var.equal = T)
    comp <- paste(y, collapse = " vs. ")
    rbind(
      data.frame(
        co = comp,
        measure = z1$data.name,
        m1 = z1$estimate[1],
        m2 = z1$estimate[2],
        t = z1$statistic,
        p = z1$p.value,
        n1 <- plot.dat$Age[x == y[[1]]] %>% length,
        n2 <- plot.dat$Age[x == y[[2]]] %>% length
      ),
      data.frame(
        co = comp,
        measure = z2$data.name,
        m1 = z2$estimate[1],
        m2 = z2$estimate[2],
        t = z2$statistic,
        p = z2$p.value,
        n1 <- plot.dat$Age[x == y[[1]]] %>% length,
        n2 <- plot.dat$Age[x == y[[2]]] %>% length
      ),
      data.frame(
        co = comp,
        measure = z3$data.name,
        m1 = z3$estimate[1],
        m2 = z3$estimate[2],
        t = z3$statistic,
        p = z3$p.value,
        n1 <- plot.dat$Age[x == y[[1]]] %>% length,
        n2 <- plot.dat$Age[x == y[[2]]] %>% length
      )
    ) %>% mutate(measure = c("chronological age", "bioenergetic age", "age delta")) %>%
      rename_at(
        colnames(.),
        ~ c(
          "clusters in comparison",
          "measure",
          "mean in group 1",
          "mean in group 2",
          "t-statistic",
          "P-value",
          "N in group 1",
          "N in group 2"
        )
      )
  })
}) %>% bind_rows %>% relocate(measure, .before = `clusters in comparison`) %>% 
  group_by(measure) %>% arrange(measure) %>% as.data.frame

## Prepare statistics for PCA plot
ac_expl_var <- {prcomp(p180.adj, scale. = T) %>% 
    summary()}$importance[2, 1:2]

ac_pca <- prcomp(p180.adj, scale. = T)$x %>%
  cbind(., data.frame("Pred. Age" = adni.data$BCAge, check.names = F)) %>% 
  cbind(., data.frame("Cluster ID" = factor(rowSums(adni.data[, c("l4", "l3")], 
                                                    na.rm = T)), check.names = F)) 

# clean up
rm(fill.c, traits, DXgroups, bp, trait, i, pvals, taglist, bpt, 
   lvls, lvlcombs, lvlslist)

#### AGES ----------------------------------------------------------------------
## Define diagnostic groups
DXgroups <- list(All = "c(0,1,3)",
                 CN = "c(0)",
                 MCI = "c(1)")

## Correlation analysis
ages.age.cor <- DXgroups %>% names %>% lapply(. , function(x){
  summary(lm(scale(Age) ~ scale(BCAge), AGES.dt %>% 
               filter(!! rlang::parse_expr(paste0("COGADNEW %in%",DXgroups[[x]])))
             ))$coefficients[2, ] %>% 
    t %>% cbind(data.frame(Cohort = "AGES", Group = x), .) %>%
    rename_at("Estimate", ~"Rho")
}) %>% bind_rows

# Recode diagnosis - limit to CN, MCI & AD
AGES.dt$cogdx.bl <- AGES.dt$COGADNEW
AGES.dt$cogdx.bl <-
  case_when(AGES.dt$cogdx.bl == 0 ~ "CN", AGES.dt$cogdx.bl == 1 ~ "MCI")
AGES.dt$cogdx.bl <- factor(AGES.dt$cogdx.bl)

AGES.dt$cogdx.fu <- AGES.dt$A2COGADNEW
AGES.dt$cogdx.fu <-
  case_when(AGES.dt$cogdx.fu == 0 ~ "CN",
            AGES.dt$cogdx.fu == 1 ~ "MCI",
            AGES.dt$cogdx.fu == 2 ~ "AD")
AGES.dt$cogdx.fu <-
  factor(AGES.dt$cogdx.fu, levels = c("CN", "MCI", "AD"))

# Remove one individual with missing Education for correct case numbers
AGES.dt <- AGES.dt %>% filter(!is.na(Education))

## Trait associations
out.ages <-
  data.frame(
    "time point" = c("baseline", NA, NA, NA, NA, "5-year follow up", NA, NA),
    trait = c(
      "Cognition [speed]",
      "Cognition [working]",
      "Cognition [memory]",
      "Grey matter volume",
      "Diagnosis [CN vs. MCI]",
      "Diagnosis [CN vs. MCI]",
      "Diagnosis [MCI vs. AD]",
      "Diagnosis [CN vs. AD]"
    ),
    check.names = F
  ) %>%
  cbind(
    rbind(
      summary(lm(
        SPEED2 ~ BCAge + Age + Sex + Education + factor(ApoE4), AGES.dt
      ))$coefficients[2, ],
      summary(lm(
        WORKING2 ~ BCAge + Age + Sex + Education + factor(ApoE4), AGES.dt
      ))$coefficients[2, ],
      summary(lm(
        MEMORY2 ~ BCAge + Age + Sex + Education + factor(ApoE4), AGES.dt
      ))$coefficients[2, ],
      summary(lm(
        GM ~ BCAge + Age + Sex + Education + factor(ApoE4) + ICV, AGES.dt
      ))$coefficients[2, ],
      summary(lm(
        BCAge ~ cogdx.bl + Age + Sex + Education + ApoE4, AGES.dt
      ))$coefficients[2, ],
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
    )
  )%>%
  cbind(data.frame("N in group 1" = c(NA, NA, NA, NA, 
                                      sum(AGES.dt$cogdx.bl == "CN", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "CN", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "MCI", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "CN", na.rm = T)),
                   "N in group 2" = c(NA, NA, NA, NA, 
                                      sum(AGES.dt$cogdx.bl == "MCI", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "MCI", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "AD", na.rm = T),
                                      sum(AGES.dt$cogdx.fu == "AD", na.rm = T))
                   , check.names = F))

# clean up
rm(DXgroups)

#### KORA ----------------------------------------------------------------------
## Correlation analysis
# combine with ADNI and AGES
age.cor.all <-
  summary(lm(scale(Age) ~ scale(BCAge), KORA.dt))$coefficients[2,] %>%
  t %>% cbind(data.frame(Cohort = "KORA", Group = "All"), .) %>%
  rename_at("Estimate", ~ "Rho") %>% rbind(., adni.age.cor) %>% 
  rbind(., ages.age.cor)

# clean up
rm(adni.age.cor, ages.age.cor)


################################################################################
#### Format and save outcomes #-------------------------------------------------
################################################################################

#### Supplementary Table 11 ----------------------------------------------------
write.xlsx(cvresults, file = "results/Supplementary Table 11.xlsx")

#### Supplementary Table 12 ----------------------------------------------------
write.xlsx(age.cor.all, file = "results/Supplementary Table 12.xlsx")

#### Supplementary Table 13 ----------------------------------------------------
write.xlsx(out.clusters, file = "results/Supplementary Table 13.xlsx")

#### Supplementary Table 14 ----------------------------------------------------
write.xlsx(out.adni, file = "results/Supplementary Table 14.xlsx")

#### Supplementary Table 15 ----------------------------------------------------
write.xlsx(out.ages, file = "results/Supplementary Table 15.xlsx")

#### Supplementary Table 19 ----------------------------------------------------
# This is added upon acceptance of the paper, hence the numbering is a bit off
write.xlsx(
  global.model$coefficients %>% as.data.frame() %>% 
  tibble::rownames_to_column("variable") %>% 
  rename_at(c("."), ~c("coefficient")) %>% 
  merge(., metpath %>% 
          select(- Biocrates.classification), 
        by.x = "variable", by.y = "ID", all.x = T) %>% 
  mutate(NiceName = case_when(is.na(NiceName) ~ variable, TRUE ~ NiceName)) %>% 
  select(- variable) %>% rename_at(c("NiceName"), ~c("variable")) %>% 
  relocate(variable, .before = "coefficient"),
  file = "results/Supplementary Table 19.xlsx")

#### Supplementary Figure 4 ----------------------------------------------------
pdf("results/Supplementary Figure 4.pdf",
    width = 11,
    height = 14 / 3)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#### Supplementary Figure 5 ----------------------------------------------------
pdf("results/Supplementary Figure 5.pdf")
ggstatsplot::ggscatterstats(KORA.dt,
  x = Age,
  y = BCAge,
  type = "pearson",
  conf.level = 0.95,
  ylab = "Predicted Bioenergetic Age",
  title = "Chronological Age vs. Predicted Bioenergetic Age [KORA]",
  point.args = list(
    color = as.numeric(KORA.dt$sample.id %in% reference$RID[reference$cohort ==
                                                              "KORA"]) + 1,
    size = 2,
    alpha = 0.4
  )
)
dev.off()

#### Supplementary Figure 6 ----------------------------------------------------
pdf("results/Supplementary Figure 6.pdf")
ggstatsplot::ggscatterstats(ac_pca, 
                            x=PC1, 
                            y=`Pred. Age`, 
                            type = "spearman", 
                            conf.level = 0.95, 
                            ylab = "Predicted Bioenergetic Age", 
                            title = "Predicted Bioenergetic Age vs. PC1 [ADNI]", 
                            point.args = list(
                              size = 2, 
                              alpha = 0.4)) + ggsci::scale_color_jama() + ggsci::scale_fill_jama()
dev.off()

#### Supplementary Figure 7 ----------------------------------------------------
adni.data$Dx2 <-
  ifelse(
    adni.data$Diagnosis.raw %in% c("CN", "SMC"),
    "CN",
    ifelse(adni.data$Diagnosis.raw %in% c("EMCI", "LMCI"), "MCI", "AD")
  )
adni.data$Dx2 <-
  factor(adni.data$Dx2, levels = c("CN", "MCI", "AD"))

p0 <-
  ggboxplot(
    adni.data,
    x = "Dx2",
    y = "Age",
    notch = T,
    fill = "Dx2",
    add.params = list(
      alpha = .5,
      fill = "Dx2",
      shape = 21
    ),
    add = "jitter",
    xlab = "Diagnosis",
    ylab = "Chronological Age"
  ) + 
  stat_compare_means(comparisons = list(c("CN", "MCI"), 
                                        c("CN", "AD"), c("MCI", "AD")),
                     method.args = list(var.equal = T))  + 
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold")) + labs(tag = "a") + 
  ggsci::scale_fill_jama()

p01 <-
  ggboxplot(
    adni.data,
    x = "Dx2",
    y = "BCAge",
    notch = T,
    fill = "Dx2",
    add.params = list(
      alpha = .5,
      fill = "Dx2",
      shape = 21
    ),
    add = "jitter",
    xlab = "Diagnosis",
    ylab = "Predicted bioenergetic age"
  ) + 
  stat_compare_means(comparisons = list(c("CN", "MCI"), c("CN", "AD"),
                                        c("MCI", "AD")),
                     method.args = list(var.equal = T)) +
  theme(legend.position = "none") +
  ggsci::scale_fill_jama()

p02 <-
  ggboxplot(
    adni.data,
    x = "Dx2",
    y = "AgeDiff",
    notch = T,
    fill = "Dx2",
    add.params = list(
      alpha = .5,
      fill = "Dx2",
      shape = 21
    ),
    add = "jitter",
    xlab = "Diagnosis",
    ylab = ""
  ) + ylab(expression(Delta('Bioenergetic - chronological age'))) + 
  stat_compare_means(comparisons = list(c("CN", "MCI"), c("CN", "AD"), 
                                        c("MCI", "AD")), 
                     method.args = list(var.equal = T)) + 
  theme(legend.position = "none") + 
  ggsci::scale_fill_jama()

adni.data$resot <- ifelse(is.na(adni.data$l4),
                          "other",
                          ifelse(adni.data$l4 == 6, "other", "cluster 7"))

p1 <-
  ggboxplot(
    adni.data,
    x = "resot",
    y = "AgeDiff",
    notch = T,
    facet.by = "Dx2",
    fill = "resot",
    add.params = list(
      alpha = .5,
      fill = "resot",
      shape = 21
    ),
    add = "jitter",
    xlab = "Cluster",
    ylab = ""
  ) + ylab(expression(Delta('Bioenergetic - chronological age'))) + 
  stat_compare_means(comparisons = list(unique(adni.data$resot)), 
                     method.args = list(var.equal = T)) + 
  coord_cartesian(clip = F) + 
  theme(legend.position = "none", plot.tag = element_text(face = "bold")) + 
  labs(tag = "d") + 
  ggsci::scale_fill_jama()

# P-value in p2 for MCI too small, to show the exact value needs override
pvals <- c("CN", "MCI", "AD") %>% 
  lapply(function(x){ 
    wilcox.test(as.formula(paste0("BCAge ~ resot")), 
           data = adni.data[adni.data$Dx2 == x,],
           var.equal = T)$p.value
  }) %>% unlist() %>% sort()

p2 <-
  ggboxplot(
    adni.data,
    x = "resot",
    y = "BCAge",
    notch = T,
    facet.by = "Dx2",
    fill = "resot",
    add.params = list(
      alpha = .5,
      fill = "resot",
      shape = 21
    ),
    add = "jitter",
    xlab = "Cluster",
    ylab = "Predicted bioenergetic age"
  ) + 
  stat_compare_means(comparisons = list(unique(adni.data$resot)),
                     method.args = list(var.equal = T),
                     symnum.args = list(cutpoints = c(0, pvals + 1e-11),
                                        symbols = sapply(pvals, function(p) {
                                          ifelse(p < 1e-4, 
                                                 sprintf("%.1e", p), 
                                                 sprintf("%.2g", p))
                                        }))) + 
  coord_cartesian(clip = F) + theme(legend.position = "none", 
                                    plot.tag = element_text(face = "bold")) + 
  labs(tag = "c") + 
  ggsci::scale_fill_jama()

p3 <-
  ggboxplot(
    adni.data,
    x = "resot",
    y = "Age",
    notch = T,
    facet.by = "Dx2",
    fill = "resot",
    add.params = list(
      alpha = .5,
      fill = "resot",
      shape = 21
    ),
    add = "jitter",
    xlab = "Cluster",
    ylab = "Chronological age"
  ) + 
  stat_compare_means(comparisons = list(unique(adni.data$resot)),
                     method.args = list(var.equal = T)) + 
  coord_cartesian(clip = F) + 
  theme(legend.position = "none", plot.tag = element_text(face = "bold")) + 
  labs(tag = "b") + 
  ggsci::scale_fill_jama()

pdf("results/Supplementary Figure 7.pdf", height = 20,width=11)
grid.arrange(p0, p01, p02, p3, p2, p1, layout_matrix = matrix(
  c(1, 2, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
  ncol = 3,
  byrow = T
))
dev.off()

#### Figure 3 ------------------------------------------------------------------
p3 <-
  ggplot(ac_pca, aes(x = PC1, y = PC2, color = `Cluster ID`)) + geom_point() +
  ggsci::scale_color_jama() +
  theme_bw() + 
  xlab(paste0("PC1 (", sprintf("%.2f", ac_expl_var[1] * 100), "%)")) +
  ylab(paste0("PC2 (", sprintf("%.2f", ac_expl_var[2] * 100), "%)")) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  ) +
  labs(tag = "a")

p4 <-
  ggplot(ac_pca, aes(x = PC1, y = PC2, color = `Pred. Age`)) + geom_point() +
  scale_color_gradient2(
    low = ggsci::pal_jama()(1),
    mid = "grey95",
    high = ggsci::pal_jama()(4)[4],
    midpoint = mean(ac_pca$`Pred. Age`, na.rm = T)
  ) +
  theme_bw() + 
  xlab(paste0("PC1 (", sprintf("%.2f", ac_expl_var[1] * 100), "%)")) +
  ylab(paste0("PC2 (", sprintf("%.2f", ac_expl_var[2] * 100), "%)")) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  ) +
  labs(tag = "b")

## Print figure
pdf("results/Figure3.pdf", width=11,height=14)
grid.arrange(p3, p4,
             plotlist[[1]], plotlist[[2]], plotlist[[3]], 
             plotlist[[4]], plotlist[[5]], plotlist[[6]],
             plotlist[[7]], plotlist[[8]], plotlist[[9]], 
             layout_matrix = matrix(c(1,1,1,2,2,2,
                                      1,1,1,2,2,2,
                                      1,1,1,2,2,2,
                                      3,3,4,4,5,5,
                                      3,3,4,4,5,5,
                                      6,6,7,7,8,8,
                                      6,6,7,7,8,8,
                                      9,9,10,10,11,11,
                                      9,9,10,10,11,11), ncol=6, byrow=T)
             )
grid.rect(
  width = 1,
  height = 3 / 9,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  ),
  x = 0.5,
  y = 5 / 6
)
grid.rect(
  width = 1,
  height = 2 / 9,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  ),
  x = 0.5,
  y = 5 / 9
)
grid.rect(
  width = 1,
  height = 2 / 9,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  ),
  x = 0.5,
  y = 3 / 9
)
grid.rect(
  width = 1,
  height = 2 / 9,
  gp = gpar(
    lwd = 0.5,
    col = "black",
    fill = NA
  ),
  x = 0.5,
  y = 1 / 9
)
dev.off()

#### Store BCAge variables for use in Script 5 ---------------------------------
# store only minimal data with follow-up diagnosis information available for AGES
AGES.dt <- AGES.dt %>% filter(!is.na(cogdx.fu)) %>% 
  select(BCAge, cogdx.fu, Age, Sex, Education, ApoE4)
save(AGES.dt, file = "data/store/AGES.dt.bcage.RData")
# store only minimal data for ADNI
adni.bcage <- adni.data %>% select(RID, BCAge)
save(adni.bcage, file = "data/store/ADNI.bcage.RData")
