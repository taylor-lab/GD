#!/usr/bin/env Rscript

library(data.table)
library(readxl)
library(survival)

##########################################################################################

# Load per-patient WGD and outcome data
os.dat <-
  as.data.table(readxl::read_excel(
    "Extended_Data_Table_2.xlsx",
    skip = 2,
    n_max = 9692
  ))

# Pan-cancer analysis
model.pancan <-
  coxph(formula = Surv(OS, as.numeric(factor(VitalStatus))) ~ WGD,
        data = os.dat[Use == T])
summary(model.pancan)


# Load curated breast cancer data: HR+/HER2- TP53wt tumors
breast.dat <- fread("BreastData.txt")

# Run Cox model
model.breast <- coxph(
  formula = Surv(as.numeric(os_days) / 30,
                 as.numeric(factor(VitalStatus))) ~ gd + DetailedTumorType2 + age_at_diagnosis + menopause_status + stage + grade + esr1_mutant,
  data = breast.dat
)
summary(model.breast)


# Load curated colorectal cancer data: KRAS-mutant tumors
colo.dat <- fread("ColorectalData.txt")

# Run Cox model
model.colo <- coxph(
  formula = Surv(as.numeric(os_days) / 30,
                 as.numeric(factor(VitalStatus))) ~ gd + age_at_diagnosis + loc + msi,
  data = colo.dat
)


# Load curated pancreas cancer data: primary pancreatic adenocarcinomas
panc.dat <- fread("PancreaticData.txt")

# Run Cox model
model.panc <- coxph(
  formula = Surv(as.numeric(os_days) / 30,
                 as.numeric(factor(VitalStatus))) ~ gd + age + resected,
  data = panc.dat
)
