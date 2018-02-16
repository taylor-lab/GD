#!/usr/bin/env Rscript

library(data.table)
library(car)

##########################################################################################

# Load data: genomic features w/ N >=30 events
wgd.dat <- fread("FeatureMatrix.txt")
dim(wgd.dat) # 9181 x 270

# Run logistic regression
model.full <-
  glm(WGD ~ ., family = binomial(link = 'logit'),
      data = wgd.dat)

# Check variance inflation factors
model.vif = as.data.table(car::vif(model.full))
model.vif[, Test := (`GVIF^(1/(2*Df))`) ^ 2]
any(model.vif$Test > 10) # FALSE

# Extract model variables
model.vars <-
  as.data.table(summary(model.full)$coefficients, 
                keep.rownames = T)

# Compute 95% confidence intervals for odds ratios
model.vars[, OR := exp(Estimate)]
model.vars[, OR_lower := exp(Estimate - 1.96 * `Std. Error`)]
model.vars[, OR_upper := exp(Estimate + 1.96 * `Std. Error`)]
model.vars = model.vars[order(`Pr(>|z|)`)]
