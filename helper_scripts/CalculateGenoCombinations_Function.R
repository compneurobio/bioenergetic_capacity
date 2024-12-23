#### Function for Calculating Allelic Effects ####
#                                               ##
#      Author: Matthias Arnold, PhD             ##
#                                               ##
#              Dec 17, 2024                     ##
#                                               ##
##################################################

# function for testing allele combinations -------------------------------------
calculateGenoCombinations <-
  function(data,
           upper = 4,
           traits = c("ADNI_MEM", "ADNI_EF", "ADAS13")) {
    c(1:upper) %>% lapply(., function(i) {
      data %>% pull(GenoClust) %>% as.character %>%
        na.omit %>% unique %>% combn(i) %>%
        apply(., 2, list) %>% lapply(., unlist) %>%
        lapply(., function(j) {
          traits %>% lapply(., function(trait) {
            tmp <- data
            tmp$test <-
              ifelse(!is.na(tmp$GenoClust),
                     ifelse(tmp$GenoClust %in% j, "s1", "s2"),
                     NA)
            form <-
              as.formula(
                paste0(
                  trait,
                  " ~ test * VISCODE + Age + Sex + ApoE4 + Education + Cohort + Diagnosis + (1 | RID)"
                )
              )
            x <-
              summary(lmer(form, tmp %>% select_at(
                .,
                c(
                  trait,
                  "RID",
                  "test",
                  "VISCODE",
                  "Age",
                  "Sex",
                  "ApoE4",
                  "Education",
                  "Cohort",
                  "Diagnosis"
                )
              ) %>% na.omit()))$coefficients
            x <- x[nrow(x), ] %>% t
            colnames(x) <- paste0(colnames(x), ".", trait)
            data.frame(x)
          }) %>% bind_cols %>% 
            cbind(data.frame(alleles = paste(j, collapse = ",")), .)
        }) %>% bind_rows
    }) %>% bind_rows %>%
      mutate(metaP = (select(., paste0(
        "Pr...t...", traits
      )) %>% apply(1, function(x)
        sumlog(x)$p)))
  }