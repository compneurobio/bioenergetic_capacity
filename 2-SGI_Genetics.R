##### Bioenergetic capacity - Script 2 #####
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
library(magrittr)
library(ggh4x)
library(variancePartition)

# clear environment
rm(list = ls())

# load association tests and plotting functions for SGI
source("helper_scripts/SGI_Association_Functions.R")
source("helper_scripts/SGI_Plotting_Functions.R")

# load helper functions for multi-SNP association tests and plotting
source("helper_scripts/Epistasis_Helper_Functions.R")
source("helper_scripts/Epistasis_Plotting_Functions.R")

#### load required data #-------------------------------------------------------
load("data/Script-1.ADNI.p180.RData")      # Contents described in Script 1
load("data/store/SGI_AC_clustering.RData") # SGI clustering results from Script 1
load("data/store/figure2_part1.RData")     # Figure 2 panels from Script 1
load("data/Script-2.ADNI.genetics.RData")  # contains three data frames:
# --------
# ADNI.dt: ADNI genotype data for SNPs identified by Shin et al., 2014, to be linked to AC levels
# --------
# df with dimensions 1548 (samples, identified by RIDs in rownames) x 32 (SNPs)
#
# --------------
# ac.haplotypes: List of SNPs per acylcarnitine for haplotype analysis 
# --------------
# df with dimensions 12 (ACs) x 3 (reported AC, corresponding p180 AC, SNP list)
#
# --------
# snpAnno: Annotation for SNPs reported by Shin et al., 2014
# --------
# df with dimensions 32 (SNPs) x 6 (rsID, chr, position, top-associated AC(s),
# interaction model to be tested (yes/no), annotated gene)

#### Derive data/formats required by this script #------------------------------
## ADNI ------
# Regress out cohort effects from acylcarnitine profiles
p180.adj <- p180[, -1] %>% apply(2, function(x){ 
  lm(x ~ p180$Cohort)$residuals %>% scale 
}) %>% as.data.frame %>% `rownames<-`(rownames(p180))

rm(p180) # not needed anymore

# Filter down genetic data to individuals with acylcarnitine data
# Once numeric for additive genetic effects on acylcarnitines
ADNI.geno.alleles <- ADNI.dt %>% 
  merge(p180.pheno %>% select(Age), ., by = 0, all.x = T) %>% 
  select(-Age) %>% arrange(as.numeric(Row.names)) %>% 
  tibble::column_to_rownames("Row.names")
# Once with allele codes for haplotype (factor/genotype effects) analysis
ADNI.geno <- ADNI.geno.alleles %>% 
  mutate(across(everything(), function(x){ x %>% factor %>% as.numeric }))

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

################################################################################
#### SGI Analysis using single SNPs #-------------------------------------------
################################################################################
# extract stored hc
hc <- sg$hc

# do SGI association for single SNPs (all SNPs per locus)
sg2 <- sgi_init(
  hc,
  minsize = ceiling(0.05 * nrow(p180.adj)),
  outcomes = ADNI.geno,
  user_defined_tests = list(numeric = my_logit)
)
# run SGI
as2 = sgi_run(sg2)

## do SGI association for single loci (best hit per locus)
sig.results <-
  as2$results %>% names %>% lapply(function(x) {
    y <- as2$results[[x]]
    y$SNP <- x
    y
  }) %>%
  bind_rows %>% filter(padj <= 0.05) %>% 
  merge(., (snpAnno %>% select(SNP, annotated.gene)), by = "SNP")

best.hits <-
  sig.results %>% group_by(annotated.gene) %>% slice_min(padj) %>% 
  as.data.frame %>% select(SNP, annotated.gene)

ADNI.loci <- ADNI.geno[, best.hits$SNP]
colnames(ADNI.loci) <- best.hits$annotated.gene

sg3 <- sgi_init(
  hc,
  minsize = ceiling(0.05 * nrow(p180.adj)),
  outcomes = ADNI.loci,
  user_defined_tests = list(numeric = my_logit)
)
# run SGI
as3 = sgi_run(sg3)

################################################################################
#### SNP-acylcarnitine associations #-------------------------------------------
################################################################################

SNP_metabolite_assocs <-
  p180.adj %>% colnames %>% lapply(function(x) {
    ADNI.geno %>% colnames %>% lapply(function(y) {
      met <- p180.adj[, x]
      snp <- as.numeric(factor(ADNI.geno[, y])) - 1
      res <- summary(lm(met ~ snp + p180.pheno$Age + p180.pheno$Sex))
      return(
        data.frame(
          SNP = y,
          acylcarnitine = x,
          beta = res$coefficients[2, 1],
          "SE(beta)" = res$coefficients[2, 2],
          "t-statistic" = res$coefficients[2, 3],
          pvalue = res$coefficients[2, 4],
          stringsAsFactors = F
        )
      )
    }) %>% bind_rows
  }) %>% bind_rows %>% 
  mutate("P-value (BH-adjusted)" = p.adjust(pvalue, method = "BH")) %>%
  as.data.frame %>% rename_at(c("pvalue"),  ~ c("P-value (raw)")) %>%
  mutate(acylcarnitine = (inner_join(., metpath, 
                                     by = c("acylcarnitine" = "ID")) %>% 
                            pull(NiceName)))

################################################################################
#### Analysis of explained variance #-------------------------------------------
################################################################################

## Explained variance of clustering by genetics
geno.cluster <- ADNI.geno %>% cbind(get_vcps(sg2) %>%
                                      data.frame %>%
                                      rename_at(colnames(.), function(x) {
                                        paste0("level.",
                                               x %>% gsub(pattern = "l", 
                                                          replacement = "") %>% 
                                                 as.numeric() - 1)
                                      }))

SNP_explained_var <-
  sig.results %>% mutate(R2_McKelvey = (apply(., 1, function(x) {
    lvl <- paste0("level.", (as.numeric(x[6]) - 1))
    form <- paste0("factor(", lvl, ") ~ ", x[1])
    r2_mckelvey(glm(as.formula(form), data = geno.cluster, family = "binomial"))
  }))) %>% as.data.frame %>% group_by(level) %>% 
  mutate(total_R2_McKelvey = (r2_mckelvey(
    glm(as.formula(
      paste0(
        "factor(level.",
        unique(as.numeric(level) - 1),
        ") ~ ",
        paste(as.character(SNP), collapse = " + ")
      )
    ), data = geno.cluster, family = "binomial")
  ))) %>% as.data.frame %>% select(-test) %>%
  mutate_at(c("level"), ~ (. - 1)) %>%
  arrange(level) %>% 
  select(SNP, level, annotated.gene, R2_McKelvey, total_R2_McKelvey) %>%
  rename_at(c("level", "annotated.gene"),
            ~ c("branching point", "annotated gene"))

#### Variance partition / Explained variance of clustering by AC levels ####
int_rs17806888_rs924135 <-
  data.frame(
    int_rs17806888_rs924135 = interaction(ADNI.geno.alleles$rs17806888, 
                                          ADNI.geno.alleles$rs924135) %>%
      as.character(),
    stringsAsFactors = F
  )
# filter rare genotype combinations
int_rs17806888_rs924135$int_rs17806888_rs924135[which(int_rs17806888_rs924135$int_rs17806888_rs924135 %in% {
  table(int_rs17806888_rs924135$int_rs17806888_rs924135) %>% names
}[which(table(int_rs17806888_rs924135$int_rs17806888_rs924135) <= 10)])] <-
  NA

exprObj <-
  t(p180.adj[rownames(p180.adj) %in% rownames(na.omit(ADNI.geno)), ])
incl.levels <- as$results %>% bind_rows %>% filter(padj <= 0.05) %>%
  rbind(., (as2$results %>% bind_rows %>% filter(padj <= 0.05))) %>%
  arrange(level) %>%   pull(level) %>% unique %>% paste0("l", .)
covars <- c("Age", "Sex", "BMI", "Education")
info <-
  cbind((get_vcps(sg) %>% as.data.frame %>% select(any_of(incl.levels))),
        (ADNI.geno %>% select(all_of(
          sig.results %>% pull(SNP) %>% unique
        ))),
        int_rs17806888_rs924135,
        (p180.pheno %>% select(all_of(
          c("Age", "Sex", "BMI", "Education")
        )))
  ) %>% mutate_at("Sex", ~ as.numeric(as.character(Sex))) %>% 
  filter(rownames(.) %in% rownames(na.omit(ADNI.geno)))
form <-
  as.formula(paste0(
    "~",
    paste((sig.results %>% pull(SNP) %>% unique), collapse = "+"),
    " + ",
    paste(covars, collapse = "+"),
    " + ",
    "(1 | int_rs17806888_rs924135) + ",
    paste(paste0("(1 | ", incl.levels, ")"), collapse = " + ")
  ))
tmp2 <- info[,-1]
tmp2[is.na(tmp2)] <- 0
info[,-1] <- tmp2
info[, incl.levels] <- lapply(info[, incl.levels], factor)
rm(tmp2)

varPart <- fitExtractVarPartModel(exprObj, form, info)
vardata <-
  data.frame(metabolite = rownames(varPart), varPart[, incl.levels])
vardata$total <- apply(vardata[,-1], 1, sum)
varPart$Residuals <- vardata$total

varPart <- varPart[, c(incl.levels, "Residuals")]

varPart2 <- varPart
colnames(varPart2) <-
  c(as.character(sapply(incl.levels, function(x) {
    paste0("cluster ", names(table(
      get_vcps(sg) %>% as.data.frame %>% select(all_of(x))
    ))[1], " vs. ", names(table(
      get_vcps(sg) %>% as.data.frame %>% select(all_of(x))
    ))[2])
  })), "total")
colnames(varPart) <- colnames(varPart2)
rownames(varPart2) <-
  as.character(unlist(sapply(rownames(varPart2), function(x)
    metpath$NiceName[metpath$ID == x])))

## Cleanup
rm(
  list = c(
    "exprObj",
    "info",
    "int_rs17806888_rs924135",
    "vardata",
    "form",
    "incl.levels"
  )
)

################################################################################
#### SGI Analysis using multiple SNPs #-----------------------------------------
################################################################################

# prepare data
snps <- colnames(ADNI.dt)
SGI_clusters <- apply(SGI_clusters, 2, as.numeric)
ADNI.dt <-
  merge(
    ADNI.dt,
    SGI_clusters,
    by.x = 0,
    by.y = "RID",
    all.x = T
  )

# List of traits to be tested
ADNI.tnames <- colnames(SGI_clusters)[-1][1:4]
ADNI.dt %<>% select(all_of(c(ADNI.tnames, snps)))
ADNI.covar <- NA

# build trait map
ADNI.traits <- data.frame(
  trait = ADNI.tnames,
  type = "logistic_regression",
  covars = NA,
  stringsAsFactors = F
)

#### ADNI ####
# lm-glm -----------------------------------------------------------------------

# run the analysis for lm-glm
re <-
  epistasisGLM(
    dt = ADNI.dt,
    traits = ADNI.traits,
    covar = ADNI.covar,
    tnames = ADNI.tnames,
    haps = haps,
    transform = T
  )
rep <- epistasisANOVA(re)
ADNI.lm <- epistasisANOVA(re, logtransform = T)


# ranking based model ----------------------------------------------------------

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
rep0 <- epistasisANOVA(re0)
ADNI.ranking <- epistasisANOVA(re0, logtransform = T)
# Store for plotting
ADNI.re0 <- re0

# Cleanup
rm(rep0, re0)
rm(rep, re)
gc()

# create output table
multi_SNP_results <-
  ADNI.ranking %>% names %>% lapply(., function(x) {
    ADNI.ranking[[x]] %>% as.data.frame %>% 
      tibble::rownames_to_column(var = "SNPs") %>% mutate(SNPs = gsub(
        x = SNPs,
        pattern = "_",
        replacement = " x "
    )) %>%
      mutate("Genetic Model" = x)
  }) %>% bind_rows() %>% 
  mutate(`Genetic Model` = factor(.$`Genetic Model`, 
                                  levels = ac.haplotypes$metabolite[1:nrow(ac.haplotypes)])) %>%
  group_by(`Genetic Model`) %>% arrange(`Genetic Model`) %>% as.data.frame %>% 
  relocate("Genetic Model", .before = "SNPs") %>%
  group_by(`Genetic Model`) %>% 
  mutate("Genetic Model" = c(`Genetic Model` %>% 
           as.character %>% head(1), rep(NA, length(`Genetic Model`) - 1))) %>% 
  as.data.frame(row.names = NULL) %>% mutate_at(colnames(.)[-c(1:2)], function(x) {
    10 ^ x
    })

# calculate R² 
ADNI.r2 <- ADNI.re0 %>% lapply(function(x) {
  x$ms %>% lapply(function(y) {
    y %>% lapply(function(z) {
      r2_mckelvey(z)
    }) %>% bind_rows()
  }) %>% bind_rows()
})


################################################################################
#### Format and save outcomes #-------------------------------------------------
################################################################################

plot.r2 <-
  data.frame(x = unlist(ADNI.r2),
             Clusters = c(
               rep("2 vs. 3", times = sapply(ADNI.r2$l2, length)),
               rep("4 vs. 5", times = sapply(ADNI.r2$l2, length)),
               rep("6 vs. 7", times = sapply(ADNI.r2$l2, length)),
               rep("8 vs. 9", times = sapply(ADNI.r2$l2, length))
             ))

r2.filter <-
  c(
    which(multi_SNP_results$l2 <= 0.05 / ncol(get_vcps(sg))),
    2 * which(multi_SNP_results$l3 <= 0.05 / ncol(get_vcps(sg))),
    3 * which(multi_SNP_results$l4 <= 0.05 / ncol(get_vcps(sg))),
    4 * which(multi_SNP_results$l5 <= 0.05 / ncol(get_vcps(sg)))
  )
plot.r2$x[-r2.filter] <- NA

multi_SNP_results %<>%
  mutate(l2_R2_McKelvey = plot.r2$x[c(1:sapply(ADNI.r2$l2, length))]) %>%
  relocate(l2_R2_McKelvey, .after = "l2") %>%
  mutate(l3_R2_McKelvey = plot.r2$x[c(1:sapply(ADNI.r2$l2, length)) + 
                                      nrow(multi_SNP_results)]) %>%
  relocate(l3_R2_McKelvey, .after = "l3") %>%
  mutate(l4_R2_McKelvey = plot.r2$x[c(1:sapply(ADNI.r2$l2, length)) + 
                                      nrow(multi_SNP_results) * 2]) %>%
  relocate(l4_R2_McKelvey, .after = "l4") %>%
  mutate(l5_R2_McKelvey = plot.r2$x[c(1:sapply(ADNI.r2$l2, length)) + 
                                      nrow(multi_SNP_results) * 3]) %>%
  relocate(l5_R2_McKelvey, .after = "l5")

#### Supplementary table 5 -----------------------------------------------------
openxlsx::write.xlsx(snpAnno, file = "results/Supplementary Table 5.xlsx")

#### Supplementary table 6 -----------------------------------------------------
openxlsx::write.xlsx(SNP_metabolite_assocs, 
                     file = "results/Supplementary Table 6.xlsx")

#### Supplementary table 7 -----------------------------------------------------
as2_cN <-
  (get_vcps(sg) * ifelse(is.na(ADNI.geno[, 1]), NA, 1)) %>% 
  apply(2, function(x) {
    table(x)
  }) %>% t %>% as.data.frame %>% 
  rename_at(c(1, 2), ~ c("N cluster 1", "N cluster 2"))

SNP_assoc <-
  as2$results %>% names %>% lapply(function(x) {
    y <- as2$results[[x]]
    y$trait <- x
    y <- cbind(y, as2_cN)
  }) %>%
  bind_rows() %>% select_if(colnames(.) != "test") %>% relocate(level,
                                                                cid1,
                                                                cid2,
                                                                `N cluster 1`,
                                                                `N cluster 2`,
                                                                trait,
                                                                stat,
                                                                pval,
                                                                padj) %>%
  rename_at(
    c("cid1", "cid2", "level", "padj", "pval", "stat", "trait"),
    ~ c(
      "ID cluster 1",
      "ID cluster 2",
      "branching point",
      "P-value (adjusted)",
      "P-value (raw)",
      "statistic",
      "SNP"
    )
  ) %>%
  mutate_at("branching point", function(x) {
    (x - 1)
  }) %>%
  mutate("annotated gene" = (
    inner_join(.[, "SNP"], snpAnno, by = "SNP") %>%
      as.data.frame %>% pull(annotated.gene)
  ))

openxlsx::write.xlsx(SNP_assoc, file = "results/Supplementary Table 7.xlsx")

#### Supplementary table 8 -----------------------------------------------------
openxlsx::write.xlsx(SNP_explained_var, 
                     file = "results/Supplementary Table 8.xlsx")

#### Supplementary table 9 -----------------------------------------------------
colnames(multi_SNP_results)[-c(1:2)] <-
  c(
    "bp_1",
    "bp_1_R2_McKelvey",
    "bp_2",
    "bp_2_R2_McKelvey",
    "bp_3",
    "bp_3_R2_McKelvey",
    "bp_4",
    "bp_4_R2_McKelvey"
  )
openxlsx::write.xlsx(multi_SNP_results, 
                     file = "results/Supplementary Table 9.xlsx")

#### Supplementary Figure 1 ----------------------------------------------------
pV <-
  plotVarPart(varPart, label.angle = 45) + 
  ggtitle("Variance decomposition per branching point") + 
  ylab("Variance explained (%)") + 
  theme(
    plot.title = element_text(size = rel(1.09)),
    axis.title = element_text(size = rel(1.12)),
    axis.text.x = element_text(size = rel(1.308)),
    axis.text.y = element_text(size = rel(1.308)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(
      t = 5.5,
      l = 5.5,
      r = 5.5,
      b = -10,
      unit = "pt"
    ),
    plot.tag = element_text(face = "bold", size = rel(1.272727)),
    plot.tag.position = c(.01, 0.995)
  ) + labs(tag = "a") + ggsci::scale_fill_jama()
pV2 <-
  as.ggplot(
    pheatmap(
      silent = T,
      border_color = "white",
      varPart2[order(varPart2$total, decreasing = T), ],
      fontsize_col = rel(12),
      fontsize_row = rel(12),
      cluster_rows = F,
      cluster_cols = F,
      display_numbers = T,
      fontsize_number = rel(10),
      legend_breaks = seq(from = 0, to = 1, by = 0.1),
      angle_col = 45
    )
  ) + ggtitle("Explained variance by acylcarnitine") + theme(
    plot.title = element_text(
      face = "plain",
      hjust = 0.225,
      size = rel(1.09)
    ),
    plot.tag = element_text(face = "bold", size = rel(1.272727)),
    plot.tag.position = c(-0.02, .983)
  ) + labs(tag = "b")

pdf("results/Supplementary Figure 1.pdf",
    height = 7,
    width = 11)
grid.arrange(pV, pV2, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 1, 
                                      1, 1, 1, 
                                      3, 2, 2, 
                                      2, 2, 2, 
                                      2, 2, 2), 
                                    nrow = 1))
dev.off()


#### Supplementary Figure 2 ----------------------------------------------------
geno.cluster <- ADNI.geno.alleles %>% 
                  cbind(get_vcps(sg2) %>%
                          data.frame %>%
                          rename_at(colnames(.), function(x) {
                            paste0("level.",
                                   x %>% 
                                     gsub(pattern = "l", replacement = "") %>% 
                                     as.numeric() - 1)
                            }))

# Plot Tree per Locus - SNPs in Outcomes
gg_tree2 = plot(
  as3,
  padj_th = 0.05,
  cluster_ids_as_label = T,
  tree_opt = "up_rect_down_tri",
  branching_point.size = 8,
  branching_point.color = "black",
  outcomes_at_middle = F
) +
  ggsci::scale_color_jama() + ggsci::scale_fill_jama()

# Format layers in the tree plot
gg_tree2$layers[[4]]$aes_params$alpha <- 0.6
gg_tree2$layers[[4]]$aes_params$size <- 4.333333333
gg_tree2$layers[[5]]$aes_params$fill <- "white"
gg_tree2$layers[[5]]$aes_params$stroke <- 0.25
gg_tree2$layers[[5]]$aes_params$shape <- 21
gg_tree2$layers[[6]]$aes_params$fill <- rgb(1, 1, 1, .6)
gg_tree2$layers[[6]]$aes_params$size <- 4.333333333
gg_tree2$layers[[6]]$aes_params$fontface <- "italic"

# plot overview
snps <-
  lapply(as2$results, function(x)
    x[x$padj <= 0.05]) %>% unlist %>% names %>% gsub(
      x = (.),
      pattern = "\\.\\w+",
      replacement = "",
      perl = T
    ) %>% unique()
snps <- snps[c(13, 4:6, 3, 2, 7, 8, 11, 10, 9, 12, 1)]

levelmap <-
  lapply(as2$results, function(x)
    x[x$padj <= 0.05]) %>% unlist
levelmap <- levelmap[grep(names(levelmap), pattern = "level")]
names(levelmap) <-
  gsub(pattern = ".level",
       replacement = "",
       x = names(levelmap))

for (ssnp in names(levelmap)) {
  slev <- paste0("l", as.character(levelmap[ssnp]))
  ADNI.geno[is.na(data.frame(sgi::get_vcps(sg2))[, slev]), ssnp] <- NA
}

# plot overview, including locus associations and SNP profiles
sfig2a <- plot_overview(
  gg_tree = gg_tree2,
  as = as2,
  outcomes = ADNI.geno[, snps] %>% apply(., 2, function(x)
    as.numeric(factor(x))) %>% as.data.frame,
  outcomes_nmax = 13,
  outcomes_at_middle = F
) +
  labs(tag = "a") +
  theme(
    plot.margin = unit(c(-0.09, 0.08,-0.07, -0.12), units = "npc"),
    plot.tag = element_text(face = "bold", size = rel(1.309091)),
    plot.tag.position = c(.1255, .90)
  ) +
  coord_cartesian(clip = "off")
for (y in 1:length(sfig2a$layers)) {
  if (exists(where = sfig2a$layers[[y]]$aes_params, x = "label")) {
    sfig2a$layers[[y]]$aes_params$fontface <- "italic"
    sfig2a$layers[[y]]$aes_params$size <- 4.3333333
    sfig2a$layers[[y]]$aes_params$hjust <- 0
    sfig2a$layers[[y]]$aes_params$x <- 100.5
  }
}

# Multi-SNP SGI scatter plots
man.plot <- ADNI.ranking %>% names %>%
  lapply(., function(x) {
    ADNI.ranking[[x]] %>%
      abs %>% as.data.frame %>%
      tibble::rownames_to_column(var = "rwnames") %>% rowwise %>%
      mutate(combo = (strsplit(rwnames, split = "_") %>% unlist %>% 
                        length > 1)) %>% as.data.frame %>%
      mutate(ac = x) %>% mutate(label = paste0(ac, "\n", gsub(
        rwnames, pattern = "_", replacement = " x "
      ), "")) %>% select(-rwnames) %>%
      group_by(ac) %>%
      mutate(l2.max = l2 == max(.$l2)) %>%
      mutate(l3.max = l3 == max(.$l3)) %>%
      mutate(l4.max = l4 == max(.$l4)) %>%
      mutate(l5.max = l5 == max(.$l5)) %>%
      as.data.frame
  }) %>% bind_rows() %>% 
  mutate(ac = factor(.$ac, levels = ac.haplotypes$metabolite[1:nrow(ac.haplotypes)])) %>%
  group_by(ac) %>% arrange(ac) %>% as.data.frame %>% mutate(xcord = 1:nrow(.)) %>% 
  group_by(ac) %>% mutate(xcord = xcord + (5 * (which(levels(ac) == ac[1]) - 1))) %>%
  mutate(xticks = mean(xcord)) %>% as.data.frame

# create star plots for multi-SNP effects
sfig2b <-
  starManhattan(data = man.plot,
                trait = "l2",
                traitmax = "l2.max") + labs(tag = "b") + 
  theme(plot.tag = element_text(face = "bold"))
sfig2c <-
  starManhattan(data = man.plot,
                trait = "l5",
                traitmax = "l5.max") + labs(tag = "c") + 
  theme(plot.tag = element_text(face = "bold"))

pdf("results/Supplementary Figure 2.pdf",
    width = 11,
    height = 35 / 3)
grid.arrange(sfig2a,
             sfig2b,
             sfig2c,
             layout_matrix = matrix(
               c(1, 1,
                 1, 1,
                 1, 1,
                 1, 1,
                 1, 1,
                 1, 1,
                 2, 3,
                 2, 3,
                 2, 3,
                 2, 3),
               nrow = 10,
               byrow = T
             ))
dev.off()

#### Figure 2 ------------------------------------------------------------------
# Generate Figure 2, panel b; other panels have been generated in Script 1
plotdata <-
  rbind(
    varPart %>% select(-total) %>% colnames %>% lapply(., function(x) {
      data.frame(
        x = varPart[, x],
        Clusters = gsub(
          x = x,
          pattern = "cluster ",
          replacement = "",
          fixed = T
        ),
        analysis = "Acylcarnitines"
      )
    }) %>% bind_rows,
    data.frame(
      x = SNP_explained_var$R2_McKelvey,
      Clusters = (
        SNP_explained_var %>% pull(`branching point`) %>% sapply(., function(x) {
          paste0(names(table(
            get_vcps(sg) %>% as.data.frame %>% select(any_of(paste0("l", x + 1)))
          ))[1], " vs. ", names(table(
            get_vcps(sg) %>% as.data.frame %>% select(any_of(paste0("l", x + 1)))
          ))[2])
        })
      ),
      analysis = "Single SNPs"
    ),
    data.frame(plot.r2[r2.filter, ], analysis = "Epistatic effects")
  ) %>% filter(!Clusters %in% c("40 vs. 41", "16 vs. 17"))



fig2b <- ggplot(plotdata, aes(x = Clusters, y = x)) +
  geom_violin(aes(group = Clusters, fill = Clusters), scale = "width") +
  geom_boxplot(aes(group = Clusters), width = 0.05) + ggsci::scale_fill_jama() +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(analysis ~ ., switch = "y", scales = c("free_y")) +
  ylab("Variance explained") + theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = rel(1.272727)),
    plot.tag.position = c(.02, 1),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    strip.text = element_text(size = 12),
    plot.margin = unit(c(0.02, .03, .01, .01), units = "npc"),
    axis.line = element_line(colour = "black")
  ) + labs(tag = "j") + force_panelsizes(rows = c(1, 0.6, 0.35))

pdf("results/Figure2.pdf",
    width = 11,
    height = 14)
grid.arrange(
  figure2[[1]],
  fig2b,
  figure2[[2]],
  figure2[[3]],
  figure2[[7]],
  figure2[[4]],
  figure2[[5]],
  figure2[[6]],
  figure2[[8]],
  figure2[[9]],
  layout_matrix = matrix(
    c(1, 1, 1, 2,
      1, 1, 1, 2,
      3, 4, 5, 6,
      7, 8, 9, 10),
    nrow = 4,
    byrow = T
  )
)
dev.off()
