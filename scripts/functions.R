# Box-cox transformation for a given lambda:
# Input a vector and a non-zero lambda
bc_transform <- function(y, lambda) (y^lambda - 1)/lambda
bc_reverse <- function(y_bc, lambda) (lambda*y_bc + 1)^(1/lambda)


# Input training dataset (non-parameter columns removed) to perform
# stepwise selection using BIC 
run_stepwise <- function(dataset){
  min_mod <- lm(error~1, data = dataset)
  full_mod <- lm(error~., data = dataset)
  selected <- step(object = min_mod,
                   scope = list(lower = min_mod, upper = full_mod),
                   direction = "both",
                   k = log(nrow(dataset)))
  return(selected)
}

# Input training dataset (non-parameter columns removed) to run regsubsets
# (from leaps package), find the model with smallest BIC, and return it
# as an lm object
run_subset <- function(dataset){
  subsets <- regsubsets(error~., data = dataset, nvmax = 15, really.big = T)
  terms <- names(coef(subsets, which.min(summary(subsets)$bic)))[-1]
  terms <- gsub("str[A-Z]*", "str", terms)
  terms <- paste(terms, collapse = " + ")
  f <- as.formula(paste("error ~ ", terms))
  lm_mod <- lm(formula = f, data = dataset)
  return(lm_mod)
}


# LOOCV with a modification:
# Instead of leaving one data out, it leaves a group of data points belonging to a system.
# The idea is that this will provide more realistic CV result for the nature of error model.

# This function takes as input:
#    1. lm formula to test
#    2. data frame to test (any one of the organized data tables like "bind", "fold", etc)
# and outputs a list containing:
#    1. $formula: lm formula that was tested
#    2. $means: mean coverage, and mean interval width
#    3. $df: data frame that contains predicted error, prediction interval, and 
#            test result of whether or not the actual error falls under upr pred interval

loocv <- function(f, dataset){ #inpu: lm formula, dataframe (e.g. fold or bind)
  df <- dataset[c("systems", "mut", "ddg_exp", "total", "error")]
  df$err_pred <- NA
  df$lwr <- NA
  df$upr <- NA
  df$test <- NA
  sys <- unique(dataset$systems)
  str_values <- unique(dataset$str) #vector of all existing str as produced by dssp
  for (i in 1:length(sys)){
    tset <- dataset[dataset$systems != sys[i], ]
    vset <- dataset[dataset$systems == sys[i], ]
    lm_i <- lm(formula = f, data = tset)
    lm_i$xlevels[["str"]] <- str_values #to prevent errors when some levels don't exist
    pred_i <- predict(lm_i, newdata = vset, interval = "prediction")
    for (r in 1:nrow(pred_i)){
      row_i <- rownames(pred_i)[r]
      df[rownames(df) == row_i, "err_pred"] <- pred_i[rownames(pred_i) == row_i, "fit"]
      df[rownames(df) == row_i, "lwr"] <- pred_i[rownames(pred_i) == row_i, "lwr"]
      df[rownames(df) == row_i, "upr"] <- pred_i[rownames(pred_i) == row_i, "upr"]
    }
  }
  df$test <- df$error < df$upr #not testing df$error > df$lwr as in original code
  out <- list()
  out$formula <- f
  out$df <- df
  out$coverage <- mean(df$test)
  out$med_upr_bd <- median(df$upr) 
  return(out)
}  


# modified version (response variable is named "error_t" instead of "error") to use 
# on box-cox transformed dataset (used in boxcox_transform.Rmd)

loocv_t <- function(f, dataset){ #inpu: lm formula, dataframe (e.g. fold or bind)
  df <- dataset[c("systems", "mut", "ddg_exp", "total", "error_t")]
  df$err_pred <- NA
  df$lwr <- NA
  df$upr <- NA
  df$test <- NA
  sys <- unique(dataset$systems)
  str_values <- unique(dataset$str) #vector of all existing str as produced by dssp
  for (i in 1:length(sys)){
    tset <- dataset[dataset$systems != sys[i], ]
    vset <- dataset[dataset$systems == sys[i], ]
    lm_i <- lm(formula = f, data = tset)
    lm_i$xlevels[["str"]] <- str_values #to prevent errors when some levels don't exist
    pred_i <- predict(lm_i, newdata = vset, interval = "prediction")
    for (r in 1:nrow(pred_i)){
      row_i <- rownames(pred_i)[r]
      df[rownames(df) == row_i, "err_pred"] <- pred_i[rownames(pred_i) == row_i, "fit"]
      df[rownames(df) == row_i, "lwr"] <- pred_i[rownames(pred_i) == row_i, "lwr"]
      df[rownames(df) == row_i, "upr"] <- pred_i[rownames(pred_i) == row_i, "upr"]
    }
  }
  df$test <- df$error_t < df$upr #not testing df$error > df$lwr as in original code
  out <- list()
  out$formula <- f
  out$df <- df
  out$coverage <- mean(df$test)
  out$med_upr_bd <- median(df$upr) 
  return(out)
}  


# Bias correction function: scatterplots of ddg_exp and ddg_predicted (FoldX) show a bias. 
# The slope of the linear fitted line results in systematic bias in error (ddg_pred - ddg_exp) 
# as well. This function calculates a "corrected" ddg_pred.
#     Input: x = char string for ddg_exp, 
#            y = char string for ddg_predicted ("total" in my datasets), 
#            dataset = dataframe to use
#     Output: corrected ddg_pred to reduce systemic bias.
correct_bias <- function(x, y, dataset){
  mod <- lm(formula(paste(y, "~", x)), data = dataset)
  intercept <- coefficients(mod)[1]
  slope <- coefficients(mod)[2]
  y_prime <- (dataset[,y] - intercept)/slope
  return(y_prime)
}





