root <- "C:/Users/sosf-/Documents/Rob PhD paper"
root <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/"
setwd(root)
getwd()

GAMMdata <- read.csv("sharksgamm1h.csv")
GAMMdata

# set predictors as factors
GAMMdata$Tidal_Phase <- factor(GAMMdata$Tidal_Phase, ordered = TRUE)
GAMMdata$ID <- factor(GAMMdata$ID, ordered = TRUE)
GAMMdata$Season <- factor(GAMMdata$Season, ordered = TRUE)
GAMMdata$Diel_Phase <- factor(GAMMdata$Diel_Phase, ordered = TRUE)
str(GAMMdata)


# get libraries - still not complete,bad internet limiting downloads
library(nlme)
library(mgcv)
# install.packages("ggplot2")
library(ggplot2) # ggsave
library(lubridate) # today
# install.packages("ggpubr")
library(ggpubr)
# install.packages(MuMIn)
library(MuMIn)
# install.packages("itsadug") # not available for R3.6.3, needs >=4
library(itsadug)

### GAMMs for ODBA response######
#### skip to below and see "code for GAMMs with 'factor smooths' on row 75"

# ACF
ACFsharksgammfit <- gam(ODBA_hour ~ s(Tidal_Phase, # formula: response ~ parameters inc smoothers
                                      # SHOULD THIS BE GAM OR GAMM???
                                      bs = "cr",
                                      k = 2
), # smoothing basis. cyclic cubic regression splines (see cyclic.cubic.spline). i.e. a penalized cubic regression splines whose ends match, up to second derivative.
# Error in place.knots(x, nk): more knots than unique data values is not allowed (same for any value of k inc 1, when in s()), because:
# In smooth.construct.cc.smooth.spec(object, dk$data, dk$knots): basis dimension, k, increased to minimum possible
# Works with CR
# length(unique(GAMMdata$Tidal_Phase)) # 3
# length(unique(GAMMdata$ODBA_hour)) # 990
# https://stackoverflow.com/questions/40056566/mgcv-how-to-set-number-and-or-locations-of-knots-for-splines
random = list(ID ~ Tidal_Phase), # not a gam parameter not gam.fit...? Where is this from?
data = GAMMdata
)
# Warning message:In smooth.construct.cr.smooth.spec(object, dk$data, dk$knots): basis dimension, k, increased to minimum possible

png(paste0(root, today(), "_GAMResiduals.png"))
acfsharkgamm <- acf(residuals(ACFsharksgammfit),
                    plot = T
)
dev.off()

acfsharkgammout <- acfsharkgamm$acf[2]

# R requires more memory to run the GAMM. The following should provide it.
memory.limit(size = 500000) # SD not run

# run the GAMM
gammdatagamm1 <- gamm(ODBA_hour ~ s(Tidal_Phase,
                                    bs = "cr",
                                    k = 2
),
correlation = corCAR1(value = acfsharkgammout),
random = list(ID = ~ODBA_hour),
data = GAMMdata,
na.action = na.omit
)
# Error in solve.default(pdMatrix(a, factor = TRUE)): system is computationally singular: reciprocal condition number = 2.5695e-224
# In addition: Warning messages:
# 1: In smooth.construct.cr.smooth.spec(object, dk$data, dk$knots) : basis dimension, k, increased to minimum possible
# 2: In lme.formula(y ~ X - 1, random = rand, data = strip.offset(mf),: nlminb problem, convergence error code = 1. message = false convergence (8)

# It means your design matrix is not invertible and therefore can't be used to develop a regression model.
# This results from linearly dependent columns, i.e. strongly correlated variables.
# Examine the pairwise covariance (or correlation) of your variables to investigate if there are any variables that can potentially be removed.
# You're looking for covariances (or correlations) >> 0.
# Alternatively, you can probably automate this variable selection by using a forward stepwise regression.

###### code for GAMMs with 'factor smooths' run on 16/12/20##########

#### run the GAMMs#####

# gamm1 with all predictors (Tidal Phase, Diel Phase, Season, Size and ID as random effect)
# ACF1
ACFsharksgammfit1 <- gam(ODBA_hour ~
                           s(Size) +
                           s(Tidal_Phase, bs = "fs") +
                           s(Diel_Phase, bs = "fs") +
                           s(Season, bs = "fs"),
                         random = list(ID = ~1),
                         data = GAMMdata)

png(paste0(root, today(), "_GAMResiduals.png"))
acfsharkgamm1 <- acf(residuals(ACFsharksgammfit1),
                     plot = T)
dev.off()

acfsharkgammout1 <- acfsharkgamm1$acf[2] # 0.1833622


# GAMM1
gammdatagamm1 <- gamm(ODBA_hour ~
                        s(Tidal_Phase, bs = "fs") +
                        s(Diel_Phase, bs = "fs") +
                        s(Size, bs = "cr") +
                        s(Season, bs = "fs"),
                      correlation = corCAR1(value = acfsharkgammout1),
                      random = list(ID = ~1),
                      data = GAMMdata,
                      na.action = na.omit)
summary(gammdatagamm1$gam)
# Family: gaussian # Link function: identity
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.07385    0.00268   27.56   <2e-16 ***
#
# Approximate significance of smooth terms:
#                      edf Ref.df      F  p-value
# s(Tidal_Phase) 1.803e+00      2 10.836 6.56e-06 ***
# s(Diel_Phase)  6.309e-01      3  0.306    0.221
# s(Size)        1.000e+00      1 33.128 9.86e-09 ***
# s(Season)      1.731e-07      1  0.000    0.680
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# R-sq.(adj) =  0.166
# Scale est. = 0.00063539  n = 2089
png(paste0(root, today(), "_GAM_ResponseVsFittedVariables1.png"))
gam.check(gammdatagamm1$gam)
dev.off()
# 'gamm' based fit - care required with interpretation.
# Checks based on working residuals may be misleading.
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
#
#                      k'      edf k-index p-value
# s(Tidal_Phase) 3.00e+00 1.80e+00      NA      NA
# s(Diel_Phase)  4.00e+00 6.31e-01      NA      NA
# s(Size)        9.00e+00 1.00e+00    0.75  <2e-16 ***
# s(Season)      2.00e+00 1.73e-07      NA      NA


# gamm2 with Tidal Phase, Diel Phase and Size
# ACF2
ACFsharksgammfit2 <- gam(ODBA_hour ~
                           s(Tidal_Phase, bs = "fs") +
                           s(Diel_Phase, bs = "fs") +
                           s(Size),
                         random = list(ID = ~1),
                         data = GAMMdata)
png(paste0(root, today(), "_GAMResiduals2.png"))
acfsharkgamm2 <- acf(residuals(ACFsharksgammfit2),
                     plot = T) ## could you check this, doesnt plot residuals here and gives an error #SD: looks fine to me
dev.off()
acfsharkgammout2 <- acfsharkgamm2$acf[2]
# GAMM2
gammdatagamm2 <- gamm(ODBA_hour ~
                        s(Tidal_Phase, bs = "fs") +
                        s(Diel_Phase, bs = "fs") +
                        s(Size),
                      correlation = corCAR1(value = acfsharkgammout2),
                      random = list(ID = ~1),
                      data = GAMMdata,
                      na.action = na.omit
)
summary(gammdatagamm2$gam)
# Family: gaussian # Link function: identity
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.07385    0.00268   27.56   <2e-16 ***
# Approximate significance of smooth terms:
#                  edf Ref.df      F  p-value
# s(Tidal_Phase) 1.803      2 10.836 6.56e-06 ***
# s(Diel_Phase)  0.631      3  0.306    0.221
# s(Size)        1.000      1 33.129 9.86e-09 ***
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.166
# Scale est. = 0.00063539  n = 2089
png(paste0(root, today(), "_GAM_ResponseVsFittedVariables2.png"))
gam.check(gammdatagamm2$gam)
dev.off()
#                   k'   edf k-index p-value
# s(Tidal_Phase) 3.000 1.803      NA      NA
# s(Diel_Phase)  4.000 0.631      NA      NA
# s(Size)        9.000 1.000    0.75  <2e-16 ***

# gamm3 with Tidal Phase and Size
# ACF3
ACFsharksgammfit3 <- gam(ODBA_hour ~
                           s(Tidal_Phase, bs = "fs") +
                           s(Size),
                         random = list(ID = ~1),
                         data = GAMMdata)
acfsharkgamm3 <- acf(residuals(ACFsharksgammfit3),
                     plot = T) ## could you check this, doesn't plot residuals here and gives an error
acfsharkgammout3 <- acfsharkgamm3$acf[2]
# GAMM3
gammdatagamm3 <- gamm(ODBA_hour ~
                        s(Tidal_Phase, bs = "fs") +
                        s(Size),
                      correlation = corCAR1(value = acfsharkgammout3),
                      random = list(ID = ~1),
                      data = GAMMdata,
                      na.action = na.omit)
summary(gammdatagamm3$gam)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.073782   0.002643   27.92   <2e-16 ***
#                  edf Ref.df     F  p-value
# s(Tidal_Phase) 1.802      2 10.83 6.56e-06 ***
# s(Size)        1.000      1 33.19 9.57e-09 ***
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.165
# Scale est. = 0.00063582  n = 2089
png(paste0(root, today(), "_GAM_ResponseVsFittedVariables3.png"))
gam.check(gammdatagamm3$gam)
dev.off()





# how to plot GAMM objects?####
# https://cran.r-project.org/web/packages/itsadug/vignettes/inspect.html




# model selection and fit via AIC####
AIC(gammdatagamm1, gammdatagamm2, gammdatagamm3)
#               df       AIC
# gammdatagamm1  9 -9458.013
# gammdatagamm2  8 -9460.013
# gammdatagamm3  7 -9461.856
logLik(gammdatagamm1) # 4738.007 (df=9)
logLik(gammdatagamm2) # 4738.007 (df=8)
logLik(gammdatagamm3) # 4737.928 (df=7)
