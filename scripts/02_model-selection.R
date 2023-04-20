### Summary ###
# input: csv files (organized data tables)
# Model selection using stepwise selection and best subset selection
# output: models.rds (list of 12 models selected from above)


library(leaps)
source("scripts/functions.R")


###   Read in organized data tables    ###

# Full tables:
fold <- read.csv("outputs/table_full_f.csv")
bind <- read.csv("outputs/table_full_b.csv")
fold_tr <- read.csv("outputs/tr_table_full_f.csv")
bind_tr <- read.csv("outputs/tr_table_full_b.csv")

# FoldX-only tables:
foldx <- read.csv("outputs/table_fxonly_f.csv")
bindx <- read.csv("outputs/table_fxonly_b.csv")
foldx_tr <- read.csv("outputs/tr_table_fxonly_f.csv")
bindx_tr <- read.csv("outputs/tr_table_fxonly_b.csv")

# remove rows that had no match with DSSP
foldx <- foldx[-c(449, 642, 643),]
foldx_tr <- foldx_tr[-c(449, 642, 643),]


###   Stepwise model selection ###

mod_step_f <- run_stepwise(fold_tr)
mod_step_b <- run_stepwise(bind_tr)

mod_step_bic_fx <- run_stepwise(foldx_tr)
mod_step_bic_bx <- run_stepwise(bindx_tr)



###    Best Subset Selection    ###


best_subset_f <- regsubsets(error~., data = fold_tr, nvmax = 15, really.big = T)
best_subset_b <- regsubsets(error~., data = bind_tr, nvmax = 15, really.big = T)

best_subset_fx <- regsubsets(error~., data = foldx_tr, nvmax = 15)
best_subset_bx <- regsubsets(error~., data = bindx_tr, nvmax = 15)

# # Coefficients of the smallest bic
# subset_terms_f <- coef(mod_subset_f, which.min(summary(mod_subset_f)$bic))
# subset_terms_b <- coef(mod_subset_b, which.min(summary(mod_subset_b)$bic))
# subset_terms_fx <- coef(mod_subset_fx, which.min(summary(mod_subset_fx)$bic))
# subset_terms_bx <- coef(mod_subset_bx, which.min(summary(mod_subset_bx)$bic))

# Select the smallest BIC and get lm objects
subsets <- list(best_subset_f, best_subset_b, best_subset_fx, best_subset_bx)
datasets <- list(f = fold, b = bind, fx = foldx, bx = bindx)
subset_lms <- list()
for (i in 1:length(subsets)){
  subset_i <- subsets[[i]]
  terms <- names(coef(subset_i, which.min(summary(subset_i)$bic)))[-1] # extract names of coef, without "(intercept)" term
  terms <- gsub("str[A-Z]*", "str", terms) # coef contains str elements if selected, we need them to be just "str" as param
  terms <- paste(terms, collapse = " + ")
  f <- as.formula(paste("error ~ ", terms))
  lm_i <- lm(formula = f, data = datasets[[i]])
  subset_lms[[i]] <- lm_i
}


### Save models as RDS ###

# In stepwise selection, we go with BIC for following reasons:
#   1. Be consistent with best subset method
#   2. Since there are so many parameters, we want a method with a larger penalty to prevent overfitting

mods <- list(mod_step_f = mod_step_f,
             mod_step_b = mod_step_b,
             mod_step_fx = mod_step_bic_fx,
             mod_step_bx = mod_step_bic_bx,
             
             mod_subset_f = subset_lms[[1]],
             mod_subset_b = subset_lms[[2]],
             mod_subset_fx = subset_lms[[3]],
             mod_subset_bx = subset_lms[[4]])


saveRDS(mods, "outputs/models.rds")









