#### Association Functions in SGI format ###
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Nov 19, 2024              ##
#                                         ##
############################################

## Note: requires package 'rms' for ordinal regression

# lm function for SGI
my_lm <- function(v, y) {
  # try to run logit
  res <- tryCatch({
    wfit <- summary(lm(v ~ y))$coefficients
    # return results if test
    list(pval = wfit[2, 4], stat = wfit[2, 1])
    
  }, error = function(e) {
    # something didn't work with the test (probably no samples in one of the groups); return NA vector
    list(pval = NA, stat = NA)
  })
  
  # return results
  return(res)
  
}

# ordinal rms logit function for SGI
my_ord_logit <- function(x, y) {
  tt = try({
    y = as.numeric(y)
    df =  data.frame(y = x, x = y)
    covars = "1"
    fit  = rms::orm(y ~ x, df, family = "probit")
    pp <-
      expm1(pchisq(
        -diff(fit$deviance),
        df = fit$stats["d.f."],
        lower.tail = F,
        log.p = T
      )) + 1.00
    c(p = unname(pp), B = unname(fit$coefficients["x"]))
  })
  if (inherits(tt, "try-error")) {
    pval = NA
    stat = NA
  } else{
    pval = tt["p"]
    stat = tt["B"]
  }
  list(pval = unname(pval), stat = unname(stat))
}

# lm function with covariate-adjustment for SGI
# Note that in this case, the argument v will be a data.frame.
# Significance is determined by an ANOVA likelihood ratio test
my_cov_test <- function(v, y) {
  # ----- input -----------------------
  # v: outcome
  # y: valid cluster pairs as a vector of factors
  # -----------------------------------
  
  v$cluster_pair = y
  # null model with only sex
  form <-
    paste0(names(v)[1], "~", paste(names(v)[2:(ncol(v) - 1)], collapse = "+"))
  mnull = lm(as.formula(form), data = v)
  # model testing the cluster pair
  form <-
    paste0(names(v)[1], "~", paste(names(v)[2:ncol(v)], collapse = "+"))
  m1 = lm(as.formula(form), data = v)
  anv = anova(mnull, m1, test = "LRT")
  pval = anv$`Pr(>Chi)`[2]
  stat = anv$`Sum of Sq`[2]
  
  # return
  return(list(pval = pval, stat = stat))
}

# ordinal logit function with covariate-adjustment for SGI
# Note that in this case, the argument v will be a data.frame.
# Significance is determined by an ANOVA Chi-square test
my_ord_cov_test <- function(v, y) {
  # ----- input -----------------------
  # v: outcome
  # y: valid cluster pairs as a vector of factors
  # -----------------------------------
  
  v$cluster_pair = as.numeric(y)
  # null model with only sex
  form <-
    paste0(names(v)[1], "~", paste(names(v)[2:ncol(v)], collapse = "+"))
  fit  = rms::orm(as.formula(form), v, family = "probit")
  anv = anova(fit, cluster_pair, test = "Chisq")
  pval = anv["cluster_pair", 3]
  stat = anv["cluster_pair", 1]
  
  # return
  return(list(pval = pval, stat = stat))
}

# logit function for SGI
my_logit <- function(v, y) {
  # try to run logit
  res <- tryCatch({
    #wfit <- summary(glm((y-min(y))~v, family="binomial"))$coefficients
    wfit <- summary(glm(y ~ v, family = "binomial"))$coefficients
    # return results if test
    list(pval = wfit[2, 4], stat = wfit[2, 1])
    
  }, error = function(e) {
    # something didn't work with the test (probably no samples in one of the groups); return NA vector
    list(pval = NA, stat = NA)
  })
  
  # return results
  return(res)
  
}

## Set attribute "type" to be printed in the SGI output
attr(my_lm, "type") <- "linear regression"
attr(my_logit, "type") <- "logistic regression"
attr(my_ord_logit, "type") <- "ordinal probit regression"
attr(my_cov_test, "type") <- "linear regression with covariate adjustment"
attr(my_ord_cov_test, "type") <- "ordinal probit regression with covariate adjustment"
