##### Bioenergetic capacity - Script 1 #####
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Dec 17, 2024              ##
#                                         ##
############################################

#### Front matter #-------------------------------------------------------------
library(sgi)
library(MASS)
library(dplyr)
library(ggpubr)
library(ggstatsplot)
library(gridExtra)
library(grid)
library(performance)
library(openxlsx)
library(pheatmap)
library(ggplotify)

# clear environment
rm(list = ls())

# source association tests and plotting functions for SGI
source("helper_scripts/SGI_Association_Functions.R")
source("helper_scripts/SGI_Plotting_Functions.R")

#### load required data #-------------------------------------------------------
load("data/Script-1.ADNI.p180.RData") # contains three data frames:
# -----
# p180: preprocessed acylcarnitine (AC) profiles
# -----
# df with dimensions 1531 (samples, identified by RIDs in rownames) x 24 
# (cohort information + 23 ACs)
#
# --------
# metpath: basic annotations of metabolites
# --------
# df with dimensions 23 (ACs) x 3 (ID, formatted name, pathway)
#
# -----------
# p180.pheno: cross-sectional demographic/clinical data
# -----------
# df with dimensions 1531 (samples, identified by RIDs in rownames) x 11
# (5 phenotypes: CSF Ab[1-42] and p-tau, FDG-PET, ADAS-Cog. 13, and 2x Diagnosis;
#  5 covariates: Age, sex, BMI, years of education, copies of APOE e4)

#### Derive data/formats required by this script #------------------------------
## ADNI ------
# Regress out cohort effects from acylcarnitine profiles
p180.adj <- p180[, -1] %>% apply(2, function(x){ 
  lm(x ~ p180$Cohort)$residuals %>% scale 
  }) %>% as.data.frame %>% `rownames<-`(rownames(p180))

rm(p180) # not needed anymore

################################################################################
#### SGI Analysis and phenotype associations #----------------------------------
################################################################################

# hierarchical clustering
hc = hclust(dist(p180.adj), method = "ward.D2")

# initialize SGI structure; minsize is set to 5% of sample size
sg = sgi_init(
  hc,
  minsize = ceiling(0.05 * nrow(p180.adj)),
  outcomes = p180.pheno,
  user_defined_tests = list(
    numeric = my_lm,
    factor = my_ord_logit,
    character = my_logit
  )
)
# run SGI
as = sgi_run(sg)


################################################################################
#### Sensitivity analysis with Covariates  #------------------------------------
################################################################################

# Covariates: sex, age, BMI, copies of APOE e4, and years of education
p180.pheno2 <- data.frame(row.names = row.names(p180.pheno))

# Define traits and their corresponding classes
traits <- list(
  CSFAb = list(column = "CSF Ab1-42", class = "ADtrait_and_covars"),
  CSFPtau = list(column = "CSF p-tau", class = "ADtrait_and_covars"),
  FDG = list(column = "FDG-PET", class = "ADtrait_and_covars"),
  ADAS = list(column = "ADAS-Cog. 13", class = "ADtrait_and_covars"),
  DX = list(column = "Diagnosis", class = "DX_and_covars")
)

# Create data frames with trait-specific classes to select correct SGI test
trait_dfs <- lapply(names(traits), function(trait) {
  info <- traits[[trait]]
  df <- data.frame(
    Trait = p180.pheno[[info$column]],
    Sex = p180.pheno$Sex,
    Age = p180.pheno$Age,
    BMI = p180.pheno$BMI,
    APOE4 = factor(p180.pheno$`APOE e4`),
    Education = p180.pheno$Education
  )
  class(df) <- c(info$class, class(df))
  df
})

# Assign the created data frames to p180.pheno2
names(trait_dfs) <- paste0(names(traits), "_cov")
for(trait in names(trait_dfs)) {
  p180.pheno2[[trait]] <- trait_dfs[[trait]]
}

# initialize SGI structure; minsize is set to 5% of sample size
# --> use covariate-adjusted version of association tests
sg_cov = sgi_init(
  hc,
  minsize = ceiling(0.05 * nrow(p180.adj)),
  outcomes = p180.pheno2,
  user_defined_tests =
    list(ADtrait_and_covars = my_cov_test,
         DX_and_covars = my_ord_cov_test)
)
as_cov = sgi_run(sg_cov)

# cleanup
rm(list = c("p180.pheno2", "traits", "trait_dfs"))

################################################################################
#### Format and save outcomes #-------------------------------------------------
################################################################################

# Create results & data/store folders if they do not exist already
# NOTE: All main scripts expect a "data" folder to exist in the working dir
if (!dir.exists("results")) {
  dir.create("results")
}
if (!dir.exists("data/store")) {
  dir.create("data/store")
}

# Store cluster assignment for follow-up analyses ------------------------------
SGI_clusters <- get_vcps(sg) %>% as.data.frame %>%
  tibble::rownames_to_column(var = "RID") %>% mutate(RID = as.numeric(RID))
save(sg, as, SGI_clusters, file = "data/store/SGI_AC_clustering.RData")

## Supplementary table 3 -------------------------------------------------------
# Results of all association tests between cluster pairs
as_cN <- apply(get_vcps(sg), 2, function(x) {
  table(x)
}) %>% t %>%
  as.data.frame %>% rename_at(c(1, 2), 
                              ~ c("Total N cluster 1", "Total N cluster 2"))

Pheno_assoc <- as$results %>% names %>%
  lapply(function(x) {
    y <- as$results[[x]]
    y$trait <- x
    y <- cbind(y, as_cN)
  }) %>%
  bind_rows() %>% 
  relocate(level,
           cid1,
           cid2,
           trait,
           stat,
           pval,
           padj,
           test,
           `Total N cluster 1`,
           `Total N cluster 2`) %>%
  rename_at(
    c("cid1", "cid2", "level", "padj", "pval", "stat"),
    ~ c(
      "ID cluster 1",
      "ID cluster 2",
      "branching point",
      "P-value (adjusted)",
      "P-value (raw)",
      "statistic"
    )
  ) %>%
  mutate_at("branching point", function(x) {
    (x - 1)
  })

# Annotate exact numbers of samples available for each phenotype/cluster
tmp <- merge(SGI_clusters, p180.pheno, by.x = "RID", by.y = 0)
tmp <- names(p180.pheno) %>% lapply(function(x) {
  names(tmp)[grep(x = names(tmp), pattern = "l\\d")] %>% lapply(function(y) {
    z <-
      tmp %>% select_at(c(x, y)) %>% na.omit() %>% group_by(get(y)) %>% 
      summarize(nonmiss = sum(!is.na(get(x))))
    data.frame(
      cid1 = as.numeric(z[1, 1]),
      cid2 = as.numeric(z[2, 1]),
      trait = x,
      nonmiss1 = z$nonmiss[1],
      nonmiss2 = z$nonmiss[2]
    )
  }) %>% bind_rows()
}) %>% bind_rows() %>% rename_at(
  c("cid1", "cid2", "nonmiss1", "nonmiss2"),
  ~ c(
    "ID cluster 1",
    "ID cluster 2",
    "N with data cluster 1",
    "N with data cluster 2"
  )
)
Pheno_assoc <-
  merge(Pheno_assoc, tmp, by = c("ID cluster 1", "ID cluster 2", "trait")) %>% 
  relocate(`branching point`, .before = `ID cluster 1`)

rm(tmp)

openxlsx::write.xlsx(Pheno_assoc, file = "results/Supplementary Table 3.xlsx")

## Supplementary table 4 -------------------------------------------------------
# Results of sensitivity analysis including covariates
Pheno_assoc_cov <- as_cov$results %>% names %>%
  lapply(function(x) {
    y <- as_cov$results[[x]]
    y$trait <- {as$results %>% names}[which({as_cov$results %>% names} == x)]
    y <- cbind(y, as_cN)
  }) %>%
  bind_rows() %>% 
  # filter down to associations significant in the original analysis
  merge(as$results %>% names %>% lapply(function(x) {
            y <- as$results[[x]]
            y$trait <- x
            y <- cbind(y, as_cN)
          }) %>%
          bind_rows() %>% filter(padj <= 0.05) %>% select(cid1, cid2, trait), 
        ., by = c("cid1", "cid2", "trait")) %>%
  relocate(level,
           cid1,
           cid2,
           trait,
           stat,
           pval,
           padj,
           test,
           `Total N cluster 1`,
           `Total N cluster 2`) %>% select(-padj) %>%
  rename_at(
    c("cid1", "cid2", "level", "pval", "stat"),
    ~ c(
      "ID cluster 1",
      "ID cluster 2",
      "branching point",
      "P-value",
      "statistic"
    )
  ) %>%
  mutate_at("branching point", function(x) {
    (x - 1)
  })

# Annotate exact numbers of samples available for each phenotype/cluster
tmp <- merge(SGI_clusters, p180.pheno, by.x = "RID", by.y = 0)
tmp <- names(p180.pheno) %>% lapply(function(x) {
  names(tmp)[grep(x = names(tmp), pattern = "l\\d")] %>% lapply(function(y) {
    z <- tmp %>% 
      select_at(c(x, y, "Sex", "Age", "BMI", "Education", "APOE e4")) %>% 
      na.omit() %>% group_by(get(y)) %>% 
      # Covariates fixed to sex, age, BMI, education, and APOE4 genotype
      summarize(nonmiss = sum(
        !is.na(get(x)) &
        !is.na(get("Sex")) &
        !is.na(get("Age")) &
        !is.na(get("BMI")) &
        !is.na(get("Education")) & 
        !is.na(get("APOE e4"))
      ))
    data.frame(
      cid1 = as.numeric(z[1, 1]),
      cid2 = as.numeric(z[2, 1]),
      trait = x,
      nonmiss1 = z$nonmiss[1],
      nonmiss2 = z$nonmiss[2]
    )
  }) %>% bind_rows()
}) %>% bind_rows() %>% rename_at(
  c("cid1", "cid2", "nonmiss1", "nonmiss2"),
  ~ c(
    "ID cluster 1",
    "ID cluster 2",
    "N with data cluster 1",
    "N with data cluster 2"
  )
)
Pheno_assoc_cov <-
  merge(Pheno_assoc_cov, tmp, by = c("ID cluster 1", "ID cluster 2", "trait")) %>% 
  relocate(`branching point`, .before = `ID cluster 1`)

rm(tmp)

openxlsx::write.xlsx(Pheno_assoc_cov, file = "results/Supplementary Table 4.xlsx")

## Figure 2 --------------------------------------------------------------------
pheno.cluster <- p180.pheno %>% cbind(get_vcps(sg) %>% 
                                        data.frame %>% 
                                        rename_at(colnames(.),function(x){ 
                                          paste0("level.", c(1:length(x))) 
                                          })) %>% 
  mutate_at("BMI", function(x){ scale(x) }) %>% 
  mutate_at("Education", function(x){ scale(x) })

# Prepare single plots
p1 <-
  SgiBoxScatter(
    pheno.cluster,
    xvar = "level.1",
    yvar = "`CSF p-tau`",
    tag = "b",
    ylbl = "CSF p-tau*"
  ) + ggsci::scale_fill_jama()
p2 <-
  SgiBoxScatter(
    pheno.cluster,
    xvar = "level.1",
    yvar = "`FDG-PET`",
    tag = "c",
    ylbl = "FDG-PET*"
  ) + ggsci::scale_fill_jama()
p3 <-
  SgiBoxScatter(
    pheno.cluster,
    xvar = "level.1",
    yvar = "`ADAS-Cog. 13`",
    tag = "d",
    ylbl = "ADAS-Cog. 13*"
  ) + ggsci::scale_fill_jama()
p4 <-
  SgiBoxScatter(pheno.cluster, xvar = "level.3", yvar = "Age") + 
  labs(tag = "e") + theme(plot.tag = element_text(face = "bold")) + 
  ggsci::scale_fill_jama()
p5 <-
  SgiBoxScatter(
    pheno.cluster,
    xvar = "level.3",
    yvar = "`CSF Ab1-42`",
    tag = "f",
    ylbl = bquote(paste("CSF A", beta[1 - 42], "*"))
  ) + ggsci::scale_fill_jama()
p6 <-
  SgiBoxScatter(pheno.cluster,
                xvar = "level.4",
                yvar = "BMI",
                ylbl = "BMI*") + labs(tag = "g") + 
  theme(plot.tag = element_text(face = "bold")) + ggsci::scale_fill_jama()
p7 <-
  SgiStackedBar(
    pheno.cluster,
    yvar = "level.4",
    xvar = "Sex",
    level = 4,
    lbls = c("Female", "Male"),
    as = as
  ) + labs(tag = "h") + theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(.01, 0.95)
  )  + ggsci::scale_fill_jama(labels = c("Female", "Male"))
p8 <-
  SgiStackedBar(
    pheno.cluster,
    yvar = "level.4",
    xvar = "Diagnosis",
    level = 4,
    as = as
  ) + labs(tag = "i") + theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(.01, 0.95)
  )  + ggsci::scale_fill_jama()

# Prepare plot of tree-phenotype associations with AC heatmap
# generate tree plot, show results for adjusted p-values <0.05
gg_tree = plot(
  as,
  padj_th = 0.05,
  cluster_ids_as_label = T,
  tree_opt = "up_rect_down_tri",
  branching_point.size = 8,
  branching_point.color = "black",
  outcomes_at_middle = F
) + ggsci::scale_fill_jama() + ggsci::scale_color_jama(na.value = "darkgrey")

# Format layers in the tree plot
gg_tree$layers[[4]]$aes_params$alpha <- 0.6
gg_tree$layers[[4]]$aes_params$size <- 4.333333333
gg_tree$layers[[5]]$aes_params$fill <- "white"
gg_tree$layers[[5]]$aes_params$stroke <- 0.25
gg_tree$layers[[5]]$aes_params$shape <- 21
gg_tree$layers[[6]]$aes_params$fill <- rgb(1, 1, 1, .6)
gg_tree$layers[[6]]$aes_params$size <- 4.333333333
gg_tree$layers[[6]]$geom_params$parse <- TRUE
gg_tree$layers[[6]]$data$outcome[1] <- "paste(\"CSF p-tau\")"
gg_tree$layers[[6]]$data$outcome[2] <- "paste(\"FDG-PET\")"
gg_tree$layers[[6]]$data$outcome[3] <- "paste(\"ADAS-Cog.\",~13)"
gg_tree$layers[[6]]$data$outcome[5] <-
  "paste(\"CSF A\", beta[1-42])"

fig2a <- plot_overview(
  gg_tree = gg_tree,
  as = as,
  #outcomes = p180.pheno,
  xdata    = p180.adj,
  data.color = colorRampPalette(
    c(
      ggsci::pal_jama()(1),
      ggsci::pal_jama()(1),
      "white",
      ggsci::pal_jama()(4)[4],
      ggsci::pal_jama()(4)[4]
    )
  )(100),
  data_title = "acylcarnitine profiles",
  draw_legends = T
) +
  labs(tag = "a") +
  theme(
    plot.margin = unit(c(-0.09,-0.01,-0.04, -0.04), units = "npc"),
    plot.tag = element_text(face = "bold", size = rel(1.309091)),
    plot.tag.position = c(.05, .90)
  )

# Format layers in the overview plot
fig2a$layers[[3]]$aes_params$size <- 4.333333333
fig2a$layers[[6]]$aes_params$size <- 4.333333333
fig2a$layers[[5]] <- NULL
fig2a$layers[[5]]$geom_params$xmin <- 104
fig2a$layers[[5]]$geom_params$xmax <- 112


#### Save Figure 2 panels for completion in next script ------------------------
figure2 <- list(fig2a,
                p1,
                p2,
                p4,
                p5,
                p6,
                p3,
                p7,
                p8)
save(figure2, file = "data/store/figure2_part1.RData")