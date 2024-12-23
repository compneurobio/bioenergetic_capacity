#### Functions for Calculating Epistatic Effects ###
#                                                 ##
#      Authors: Matthias Arnold, PhD,             ##
#               Mustafa Buyukozkan, PhD           ##
#                                                 ##
#               Dec 17, 2024                      ##
#                                                 ##
####################################################

# function for running lm-glm --------------------------------------------------
epistasisGLM <-
  function(dt, traits, haps, covar, tnames, transform = F) {
    re <- seq(nrow(traits)) %>% {
      names(.) = tnames
      .
    } %>%
      lapply(function(trow) {
        itrait = traits$trait[trow]
        if (traits$type[trow] == "linear_regression") {
          func <- lm
        } else if (traits$type[trow] == "logistic_regression" &
                   itrait != "Progression") {
          func <- function(...)
            glm(..., family = "binomial")
          if (transform) {
            dt[, itrait] <- factor(dt[, itrait] - min(dt[, itrait], na.rm = T))
          }
        } else{
          return(NULL)
        }
        # fit without SNPs =  base model =  null model
        covars = paste(na.omit(c(covar, traits$covars[trow])), collapse = " + ")
        covars = ifelse(covars == "", 1, covars)
        formula_null = as.formula(paste0(itrait, " ~ ", covars))
        m_null <- func(formula_null, dt)
        
        # fits for each combination of SNPs for each haplotypes
        fits <- lapply(haps, lapply, function(snps) {
          # if a dummy variables has <10 1s then remove that one
          # problematic factor levels with too few occurance
          # for better fit and anova test
          hap = apply(dt[, snps, drop = F], 1, paste, collapse = "")
          hap[apply(dt[, snps, drop = F], 1, function(x) {
            any(is.na(x))
          })] <- "A0missing"
          hap[!hap %in% names(table(hap[!is.na(dt[, tnames[trow]])]))[table(hap[!is.na(dt[, tnames[trow]])]) >
                                                                        9]]  <- "A0missing"
          hap = model.matrix( ~ hap)
          hap = hap[, colSums(hap) > 9]
          hap = hap[, colnames(hap) != "hapA0missing"]
          hap = hap[, colnames(hap) != "(Intercept)"]
          
          formula1 = as.formula(paste0(itrait, " ~ ", "hap", " + ", covars))
          m1 <- func(formula1, data.frame(dt))
          m1
        })
        list(m_null = m_null, ms = fits)
      })
    return(re)
  }

# function for ANOVA test ------------------------------------------------------
epistasisANOVA <- function(re, logtransform = F) {
  rep <-
    lapply(re, function(nod)
      lapply(nod$ms, sapply, function(m1)
        try(anova(m1, nod$m_null, test = "LRT") %>% {
          .[[length(.)]][2]
        })
        %>% {
          if (inherits(., "try-error"))
            NA
          else
            .
        }))
  
  if (logtransform) {
    # log10 pvalues
    rep <-
      haps %>% names %>% {
        names(.) = .
        .
      } %>% lapply(function(x)
        try(sapply(rep, function(y)
          y[[x]] %>% log10)))
  }
  return(rep)
}

# function for ranking-based model ---------------------------------------------
epistasisRanking <- function(dt, traits, haps, covar, tnames, re) {
  # for the trick of obtaining transformed variable in arbitrary environment
  env =  as.environment(as.list(dt))
  parent.env(env) = environment()
  
  # fit ranking model for all of the variables except the binary ones for which
  # glm-binomial is fitted
  re0 <- seq(nrow(traits)) %>% {
    names(.) = tnames
    .
  } %>%
    lapply(function(trow) {
      itrait = traits$trait[trow]
      if (traits$type[trow] == "linear_regression") {
        func <- coxph
      } else if (traits$type[trow] == "logistic_regression" &
                 itrait != "Progression") {
        # glms-binary fits already obtained in previous analysis
        return(NULL)
      } else{
        return(NULL)
      }
      # null model
      covars = paste(na.omit(c(covar, traits$covars[trow])), collapse = " + ")
      formula_null = as.formula(paste0("S ~ ", covars))
      S = eval(parse(text = itrait), envir = env) %>% {
        . - min(., na.rm = T) + 1
      }
      S = Surv(S, S > 0)
      m_null <- func(formula_null, data.frame(S = S, dt))
      
      # fits with SNPs
      fits <- lapply(haps, lapply, function(snps) {
        # if a dummy variables has <10 1s then remove that one
        # problematic factor levels with too few occurance
        # for better fit and anova test
        hap = apply(dt[, snps, drop = F], 1, paste, collapse = "")
        hap[apply(dt[, snps, drop = F], 1, function(x) {
          any(is.na(x))
        })] <- "A0missing"
        hap[!hap %in% names(table(hap[!is.na(dt[, tnames[trow]])]))[table(hap[!is.na(dt[, tnames[trow]])]) >
                                                                      9]]  <- "A0missing"
        hap = model.matrix( ~ hap)
        hap = hap[, colSums(hap) > 9]
        hap = hap[, colnames(hap) != "hapA0missing"]
        hap = hap[, colnames(hap) != "(Intercept)"]
        
        # to reduce risk of overfitting if a haplotypes has too many( i.e. >15)
        # degrees of freedom, reduce it with penalization (df -> 15)
        formula1 <-
          if (ncol(hap) > 15)
            as.formula(paste0("S ~ ", "ridge( hap, df =15)", " + ", covars))
        else
          as.formula(paste0("S ~ ", "hap", " + ", covars))
        
        m1 <-
          suppressWarnings(try(func(formula1, data.frame(S = S, dt)))
          )
        m1
      })
      
      list(m_null = m_null, ms = fits)
    })
  
  # exports glm-binomial fits from previous results
  re0 = re0[!sapply(re0, is.null)]
  re0[setdiff(names(re), names(re0))] <-
    re[setdiff(names(re), names(re0))]
  
  return(re0)
}