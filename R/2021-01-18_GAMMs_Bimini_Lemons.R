# what script does
# who wrote it
# their email addresses
# date started


# load packages ####
# RB: still not complete,bad internet limiting downloads
library(nlme)
library(mgcv)
library(plotfunctions)
library(MuMIn)
library(ggplot2)
# RB: current internet prohibiting download of below packages
library(lubridate)
# install.packages("ggpubr")
library(ggpubr)
# install.packages("itsadug")
library(itsadug)
# install.packages("tidyverse")
library(tidyverse)

# latest edit date
today() # "2021-01-17"

# set WD ####
root <- "C:/Users/sosf-/Documents/Rob PhD paper"
root <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/"
setwd(root)
getwd()

# read data either 1h or 5min ####
GAMMdata <- read.csv("./data/sharksgamm1h.csv") # 1h
GAMMtimeRes <- "1h"
GAMMdata <- read.csv("./data/sharksgamm5min.csv") # 5 min
GAMMtimeRes <- "5m"
GAMMdata

# set categorical predictors as factors and view data ####
GAMMdata$Tidal_Phase <- factor(GAMMdata$Tidal_Phase, ordered = TRUE)
GAMMdata$ID <- factor(GAMMdata$ID, ordered = TRUE)
GAMMdata$Season <- factor(GAMMdata$Season, ordered = TRUE)
GAMMdata$Diel_Phase <- factor(GAMMdata$Diel_Phase, ordered = TRUE)
str(GAMMdata)

# chi square tests for independence of predictors - all are independent ####
chisq.test(GAMMdata$Diel_Phase, GAMMdata$Tidal_Phase, correct = FALSE)
# X-squared = 54.58, df = 6, p-value = 5.635e-10
chisq.test(GAMMdata$Diel_Phase, GAMMdata$Season, correct = FALSE)
# X-squared = 61.122, df = 3, p-value = 3.384e-13
chisq.test(GAMMdata$Diel_Phase, GAMMdata$Size, correct = FALSE)
# X-squared = 227.69, df = 48, p-value < 2.2e-16
chisq.test(GAMMdata$Tidal_Phase, GAMMdata$Season, correct = FALSE)
# X-squared = 0.090629, df = 2, p-value = 0.9557
chisq.test(GAMMdata$Size, GAMMdata$Tidal_Phase, correct = FALSE)
# X-squared = 28.35, df = 32, p-value = 0.652
chisq.test(GAMMdata$Season, GAMMdata$Size, correct = FALSE)
# X-squared = 22023, df = 16, p-value < 2.2e-16


# KW tests, factor expvars vs all resvars, ODBA response ####
attach(GAMMdata)
kruskal.test(ODBA ~ Tidal_Phase)
# 1h: significant
# 5m: Kruskal-Wallis chi-squared = 354.58, df = 2, p-value < 2.2e-16
pairwise.wilcox.test(ODBA_hour, Tidal_Phase, p.adjust.method = "BH")
# bug ####
# object 'ODBA_hour' not found
kruskal.test(ODBA ~ Diel_Phase)
# 1h: significant but only just
# 5m: Kruskal-Wallis chi-squared = 150.87, df = 3, p-value < 2.2e-16
pairwise.wilcox.test(ODBA, Diel_Phase, p.adjust.method = "BH")
# 1h:
#
# 5m:
# 1       2       3
# 2 0.11    -       -
# 3 5.4e-08 < 2e-16 -
# 4 1.9e-06 2.2e-09 < 2e-16
# P value adjustment method: BH
kruskal.test(ODBA ~ Season)
# 1h: no significance
# Kruskal-Wallis chi-squared = 2.2751, df = 1, p-value = 0.1315


# KW tests, factor expvars vs all resvars, resting response ####
kruskal.test(resting ~ Tidal_Phase)
# 1h: significant
# 5m: Kruskal-Wallis chi-squared = 202.54, df = 2, p-value < 2.2e-16
pairwise.wilcox.test(resting, Tidal_Phase, p.adjust.method = "BH")
# 1h:
#
# 5m:
#   1       2
# 2 6.8e-05 -
# 3 < 2e-16 < 2e-16
kruskal.test(resting ~ Diel_Phase)
# 1h: no significance
# 5m: Kruskal-Wallis chi-squared = 40.078, df = 3, p-value = 1.026e-08
kruskal.test(resting ~ Season)
# 1h: no significance
# 5m: Kruskal-Wallis chi-squared = 73.554, df = 1, p-value < 2.2e-16


# KW tests, factor expvars vs all resvars, bursting response ####
kruskal.test(bursting ~ Tidal_Phase)
# 1h: significant
# 5m: Kruskal-Wallis chi-squared = 37.274, df = 2, p-value = 8.055e-09
pairwise.wilcox.test(bursting, Tidal_Phase, p.adjust.method = "BH")
# 1h:
#
# 5m:
#   1       2
# 2 0.24    -
# 3 1.9e-08 5.0e-06
kruskal.test(bursting ~ Diel_Phase)
# 1h: significant
# 5m: Kruskal-Wallis chi-squared = 92.943, df = 3, p-value < 2.2e-16
pairwise.wilcox.test(bursting, Diel_Phase, p.adjust.method = "BH")
# 1h:
#
# 5m:
#   1       2       3
# 2 0.24479 -       -
# 3 1.4e-10 6.4e-13 -
# 4 0.34759 0.00019 < 2e-16
kruskal.test(bursting ~ Season)
# 1h: no significance
# 5m: Kruskal-Wallis chi-squared = 14.055, df = 1, p-value = 0.0001776






# GAMMs: ODBA response: GAMM1 all predictors ####
# Tidal Phase, Diel Phase, Season, Size and ID as random effect

# ACF for GAMM 1 ODBA
ACFODBAgammfit1 <- gam((ODBA) ~ s(Size) + s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Season, bs = "fs"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_ODBA_GAMM1.png"))
acfODBAgamm1 <- acf(residuals(ACFODBAgammfit1), plot = T)
dev.off()
# @RB ####
# auto-export plots to image files. you get the results in a folder and don't have to re-run everything to check plot results when you come back to this months/years later.
acfODBAgammout1 <- acfODBAgamm1$acf[2]

# GAMM1 ODBA
start_time <- Sys.time() # 2021-01-15 started ~ 11:30 ended before 20:34 (maybe up to 1 hour more)
gammdataODBA1 <- gamm(ODBA ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr") + s(Season, bs = "fs"),
  correlation = corCAR1(value = acfODBAgammout1),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_ODBA_1.txt"))
save(gammdataODBA1, file = "gammdataODBA1.RData", compress = "xz")

summary(gammdataODBA1$gam)
# 5m:
# Family: gaussian
# Link function: identity
#
# Formula:
#   ODBA ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Size, bs = "cr") + s(Season, bs = "fs")
#
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.074876   0.002993   25.02   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value
# s(Tidal_Phase) 1.958e+00      2 56.147  < 2e-16 ***
#   s(Diel_Phase)  2.663e+00      3  7.321 2.59e-05 ***
#   s(Size)        1.000e+00      1 33.573  < 2e-16 ***
#   s(Season)      5.792e-06      1  0.000     0.78
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# R-sq.(adj) =  0.111
# Scale est. = 0.00099289  n = 25075
gam.check(gammdataODBA1$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.00e+00 1.96e+00      NA      NA
# s(Diel_Phase)  4.00e+00 2.66e+00      NA      NA
# s(Size)        9.00e+00 1.00e+00    0.93  <2e-16 ***
#   s(Season)      2.00e+00 5.79e-06      NA      NA
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# GAMMs: ODBA response: GAMM2 no season ####
# Tidal Phase, Diel Phase and Size predictors (ommitting Season)

# ACF for GAMM2 ODBA
ACFODBAgammfit2 <- gam((ODBA) ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_ODBA_GAMM2.png"))
acfODBAgamm2 <- acf(residuals(ACFODBAgammfit2), plot = T)
dev.off()
acfODBAgammout2 <- acfODBAgamm2$acf[2]

# GAMM2 ODBA
start_time <- Sys.time()
gammdataODBA2 <- gamm(ODBA ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr"),
  correlation = corCAR1(value = acfODBAgammout2),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_ODBA_2.txt"))
save(gammdataODBA2, file = "gammdataODBA2.RData", compress = "xz")
summary(gammdataODBA2$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   ODBA ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Size, bs = "cr")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.074876   0.002993   25.02   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value
# s(Tidal_Phase) 1.958      2 56.147  < 2e-16 ***
#   s(Diel_Phase)  2.663      3  7.321 2.59e-05 ***
#   s(Size)        1.000      1 33.572  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.111
# Scale est. = 0.00099289  n = 25075
gam.check(gammdataODBA2$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                  k'  edf k-index p-value
# s(Tidal_Phase) 3.00 1.96      NA      NA
# s(Diel_Phase)  4.00 2.66      NA      NA
# s(Size)        9.00 1.00    0.95   0.005 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# GAMMs: ODBA response: GAMM3 no season or diel ####
# Tidal Phase and Size predictors (ommitting season and Diel phase)

# ACF3 ODBA
ACFODBAgammfit3 <- gam((ODBA) ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_ODBA_GAMM3.png"))
acfODBAgamm3 <- acf(residuals(ACFODBAgammfit3), plot = T)
dev.off()
acfODBAgammout3 <- acfODBAgamm3$acf[2]

# GAMM3 ODBA
start_time <- Sys.time()
gammdataODBA3 <- gamm(ODBA ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr"),
  correlation = corCAR1(value = acfODBAgammout3),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_ODBA_3.txt"))
save(gammdataODBA3, file = "gammdataODBA3.RData", compress = "xz")
summary(gammdataODBA3$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   ODBA ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.07408    0.00275   26.93   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(Tidal_Phase) 1.957      2 56.17  <2e-16 ***
#   s(Size)        1.000      1 33.61  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# R-sq.(adj) =  0.109
# Scale est. = 0.00099499  n = 25075
gam.check(gammdataODBA3$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                  k'  edf k-index p-value
# s(Tidal_Phase) 3.00 1.96      NA      NA
# s(Size)        9.00 1.00    0.95  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

png(paste0(root, today(), "_plot.gam_ACFODBAgammfit3_", GAMMtimeRes, ".png"))
plot.gam(ACFODBAgammfit3)
dev.off()
png(paste0(root, today(), "_ODBA-Size_", GAMMtimeRes, ".png"))
plot(ODBA ~ Size, data = GAMMdata)
dev.off()
png(paste0(root, today(), "_ODBA-Tidal_", GAMMtimeRes, ".png"))
plot(ODBA ~ Tidal_Phase, data = GAMMdata) # not useful, will run some interaction.plots of means to represent factor interactions.
dev.off()

# model selection and fit via AIC
AIC(gammdataODBA1, gammdataODBA2, gammdataODBA3)
# 5m:
#               df       AIC
# gammdataODBA1  9 -107403.6
# gammdataODBA2  8 -107405.6
# gammdataODBA3  7 -107393.3
# plotting smooths for size against ODBA







# GAMMs: Resting response: GAMM1 all predictors ####
# Tidal Phase, Diel Phase, Season, Size and ID as random effect

# ACF for GAMM 1 resting
ACFrestinggammfit1 <- gam((resting) ~ s(Size, bs = "cr") + s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Season, bs = "fs"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Resting_GAMM1.png"))
acfrestinggamm1 <- acf(residuals(ACFrestinggammfit1), plot = T)
dev.off()
acfrestinggammout1 <- acfrestinggamm1$acf[2]

# GAMM1 resting
start_time <- Sys.time()
gammdataresting1 <- gamm(resting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr") + s(Season, bs = "fs"),
  correlation = corCAR1(value = acfrestinggammout1),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Resting_1.txt")) # as.character strips the unit, don't know if number is seconds minutes or hours
save(gammdataresting1, file = "gammdataresting1.RData", compress = "xz")
# load("gammdataresting1.RData")
summary(gammdataresting1$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   resting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Size, bs = "cr") + s(Season, bs = "fs")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    8.076      1.747   4.624 3.79e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value
# s(Tidal_Phase) 1.930e+00      2 30.821 < 2e-16 ***
#   s(Diel_Phase)  1.396e+00      3  1.033 0.08781 .
# s(Size)        1.000e+00      1 10.770 0.00103 **
#   s(Season)      1.458e-05      1  0.000 0.95981
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0413
# Scale est. = 456.08    n = 25075
gam.check(gammdataresting1$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.00e+00 1.93e+00      NA      NA
# s(Diel_Phase)  4.00e+00 1.40e+00      NA      NA
# s(Size)        9.00e+00 1.00e+00    0.99    0.34
# s(Season)      2.00e+00 1.46e-05      NA      NA



# GAMMs: Resting response: GAMM2 no season ####
# Tidal Phase, Diel Phase and Size predictors (ommitting Season)

# ACF for GAMM2 resting
ACFrestinggammfit2 <- gam((resting) ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Resting_GAMM2.png"))
acfrestinggamm2 <- acf(residuals(ACFrestinggammfit2), plot = T)
dev.off()
acfrestinggammout2 <- acfrestinggamm2$acf[2]

# GAMM2 resting
start_time <- Sys.time()
gammdataresting2 <- gamm(resting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr"),
  correlation = corCAR1(value = acfrestinggammout2),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Resting_2.txt"))
save(gammdataresting2, file = "gammdataresting2.RData", compress = "xz") # Try to use saveRDS instead of save.
# load("gammdataresting2.RData")
summary(gammdataresting2$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   resting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Size, bs = "cr")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    8.077      1.746   4.626 3.74e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value
# s(Tidal_Phase) 1.929      2 30.819 < 2e-16 ***
#   s(Diel_Phase)  1.393      3  1.031 0.08777 .
# s(Size)        1.000      1 10.769 0.00103 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0413
# Scale est. = 456.08    n = 25075
gam.check(gammdataresting2$gam)
#5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                  k'  edf k-index p-value
# s(Tidal_Phase) 3.00 1.93      NA      NA
# s(Diel_Phase)  4.00 1.39      NA      NA
# s(Size)        9.00 1.00    0.95  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# GAMMs: Resting response: GAMM3 no season or diel ####
# Tidal Phase and Size predictors (ommitting season and Diel phase)
# ACF3 resting
ACFrestinggammfit3 <- gam((resting) ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Resting_GAMM3.png"))
acfrestinggamm3 <- acf(residuals(ACFrestinggammfit3), plot = T)
dev.off()
acfrestinggammout3 <- acfrestinggamm3$acf[2]

# GAMM3 resting
start_time <- Sys.time()
gammdataresting3 <- gamm(resting ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr"),
  correlation = corCAR1(value = acfrestinggammout3),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Resting_3.txt"))
save(gammdataresting3, file = "gammdataresting3.RData", compress = "xz")
# load("gammdataresting3.RData")
summary(gammdataresting3$gam)
#5m:
# Family: gaussian
# Link function: identity
# Formula:
#   resting ~ s(Tidal_Phase, bs = "fs") + s(Size, bs = "cr")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    8.235      1.713   4.808 1.53e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(Tidal_Phase) 1.929      2 30.84 < 2e-16 ***
#   s(Size)        1.000      1 10.74 0.00105 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0404
# Scale est. = 456.44    n = 25075
gam.check(gammdataresting3$gam)
#5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                  k'  edf k-index p-value
# s(Tidal_Phase) 3.00 1.93      NA      NA
# s(Size)        9.00 1.00    0.97   0.005 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##### plotting smooths for size against resting
png(paste0(root, today(), "_plot.gam_ACFrestinggammfit3_", GAMMtimeRes, ".png"))
plot.gam(ACFrestinggammfit3)
dev.off()
png(paste0(root, today(), "_resting-Size_", GAMMtimeRes, ".png"))
plot(resting ~ Size, data = GAMMdata)
dev.off()

##### model selection and fit via AIC
AIC(gammdataresting1, gammdataresting2, gammdataresting3)
#5m:
#                  df      AIC
# gammdataresting1  9 209591.3
# gammdataresting2  8 209589.3
# gammdataresting3  7 209588.4









# GAMMs: Bursting response: GAMM1 all predictors ####
# Tidal Phase, Diel Phase, Season, Size and ID as random effect

# ACF for GAMM 1 bursting
ACFburstinggammfit1 <- gam((bursting) ~ s(Size, bs = "cr") + s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Season, bs = "fs"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Bursting_GAMM1.png"))
acfburstinggamm1 <- acf(residuals(ACFburstinggammfit1), plot = T)
dev.off()
acfburstinggammout1 <- acfburstinggamm1$acf[2]

# GAMM1 bursting
start_time <- Sys.time()
gammdatabursting1 <- gamm(bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Size, bs = "cr") + s(Season, bs = "fs"),
  correlation = corCAR1(value = acfburstinggammout1),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Bursting_1.txt"))
save(gammdatabursting1, file = "gammdatabursting1.RData", compress = "xz")
summary(gammdatabursting1$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Size, bs = "cr") + s(Season, bs = "fs")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   4.7725     0.5616   8.498   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(Tidal_Phase) 4.306e-07      2 0.000  0.4842
# s(Diel_Phase)  2.103e+00      3 2.397  0.0176 *
# s(Size)        1.000e+00      1 3.573  0.0588 .
# s(Season)      2.862e-05      1 0.000  0.4829
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.00536
# Scale est. = 194.51    n = 25075
gam.check(gammdatabursting1$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.00e+00 4.31e-07      NA      NA
# s(Diel_Phase)  4.00e+00 2.10e+00      NA      NA
# s(Size)        9.00e+00 1.00e+00    0.97   0.015 *
# s(Season)      2.00e+00 2.86e-05      NA      NA
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# GAMMs: Bursting response: GAMM2 no season ####
# Tidal Phase, Diel Phase and Size predictors (ommitting Season)

# ACF for GAMM2 bursting
ACFburstinggammfit2 <- gam((bursting) ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Season, bs = "fs"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Bursting_GAMM2.png"))
acfburstinggamm2 <- acf(residuals(ACFburstinggammfit2), plot = T) ## could you check this, doesnt plot residuals here and gives an error
dev.off()
acfburstinggammout2 <- acfburstinggamm2$acf[2]

# GAMM2 bursting
start_time <- Sys.time()
gammdatabursting2 <- gamm(bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") + s(Season, bs = "fs"),
  correlation = corCAR1(value = acfburstinggammout2),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Bursting_2.txt"))
save(gammdatabursting2, file = "gammdatabursting2.RData", compress = "xz")
summary(gammdatabursting2$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs") +
#   s(Season, bs = "fs")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   4.7502     0.5988   7.933 2.23e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value
# s(Tidal_Phase) 6.319e-07      2 0.00  0.4889
# s(Diel_Phase)  2.117e+00      3 2.41  0.0176 *
# s(Season)      2.473e-05      1 0.00  0.5506
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.00134
# Scale est. = 194.51    n = 25075
gam.check(gammdatabursting2$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.00e+00 6.32e-07      NA      NA
# s(Diel_Phase)  4.00e+00 2.12e+00      NA      NA
# s(Season)      2.00e+00 2.47e-05      NA      NA

# GAMMs: Bursting response: GAMM3 no season or diel ####
# Tidal Phase and Size predictors (ommitting season and Diel phase)

# ACF3 bursting
ACFburstinggammfit3 <- gam((bursting) ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs"), random = list(ID = ~1), data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals_", GAMMtimeRes, "_Bursting_GAMM3.png"))
acfburstinggamm3 <- acf(residuals(ACFburstinggammfit3), plot = T) ## could you check this, doesnt plot residuals here and gives an error
dev.off()
acfburstinggammout3 <- acfburstinggamm3$acf[2]

# GAMM3 bursting
start_time <- Sys.time()
gammdatabursting3 <- gamm(bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs"),
  correlation = corCAR1(value = acfburstinggammout3),
  random = list(ID = ~1),
  data = GAMMdata,
  na.action = na.omit
)
end_time <- Sys.time()
run_time <- end_time - start_time
writeLines(paste0(as.character(run_time), " secs"), con = paste0("gammruntime_Bursting_3.txt"))
save(gammdatabursting3, file = "gammdatabursting3.RData", compress = "xz")
summary(gammdatabursting3$gam)
# 5m:
# Family: gaussian
# Link function: identity
# Formula:
#   bursting ~ s(Tidal_Phase, bs = "fs") + s(Diel_Phase, bs = "fs")
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   4.7502     0.5988   7.933 2.22e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value
# s(Tidal_Phase) 0.0001793      2 0.00  0.4889
# s(Diel_Phase)  2.1166872      3 2.41  0.0176 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.00134
# Scale est. = 194.51    n = 25075
gam.check(gammdatabursting3$gam)
# 5m:
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.000000 0.000179      NA      NA
# s(Diel_Phase)  4.000000 2.116687      NA      NA

##### plotting factors against bursting
png(paste0(root, today(), "_plot.gam_ACFburstinggammfit1_", GAMMtimeRes, ".png"))
plot.gam(ACFburstinggammfit1)
dev.off()
png(paste0(root, today(), "_bursting-Size_", GAMMtimeRes, ".png"))
plot(bursting ~ Size, data = GAMMdata) #### just for fun
dev.off()

##### model selection and fit via AIC
AIC(gammdatabursting1, gammdatabursting2, gammdatabursting3)
# 5m:
# df      AIC
# gammdatabursting1  9 193708.4
# gammdatabursting2  7 193707.7
# gammdatabursting3  6 193705.7
