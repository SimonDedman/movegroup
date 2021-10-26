root <- "C:/Users/sosf-/Documents/Rob PhD paper"
root <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/"
setwd(root)
getwd()

GAMMdata <- read.csv("sharksgamm1h.csv")
GAMMdata
library(nlme)
library(mgcv)
library(ggplot2) # ggsave
library(lubridate) #today

#GAMM for ODBA response
#ACF
ACFsharksgammfit <- gam(ODBA_hour ~ s(Tidal_Phase, # formula: response ~ parameters inc smoothers
                                      # SHOULD THIS BE GAM OR GAMM???
                                      bs = "cr",
                                      k = 2), # smoothing basis. cyclic cubic regression splines (see cyclic.cubic.spline). i.e. a penalized cubic regression splines whose ends match, up to second derivative.
                        # Error in place.knots(x, nk): more knots than unique data values is not allowed (same for any value of k inc 1, when in s()), because:
                        # In smooth.construct.cc.smooth.spec(object, dk$data, dk$knots): basis dimension, k, increased to minimum possible
                        # Works with CR
                        # length(unique(GAMMdata$Tidal_Phase)) # 3
                        # length(unique(GAMMdata$ODBA_hour)) # 990
                        # https://stackoverflow.com/questions/40056566/mgcv-how-to-set-number-and-or-locations-of-knots-for-splines
                        random = list(ID ~ Tidal_Phase), # not a gam parameter not gam.fit...? Where is this from?
                        data = GAMMdata)
# Warning message:In smooth.construct.cr.smooth.spec(object, dk$data, dk$knots): basis dimension, k, increased to minimum possible

png(paste0(root, today(), "_GAMResiduals.png"))
acfsharkgamm <- acf(residuals(ACFsharksgammfit),
                    plot = T)
dev.off()

acfsharkgammout <- acfsharkgamm$acf[2]

#R requires more memory to run the GAMM. The following should provide it.
memory.limit(size = 500000) # SD not run

#run the GAMM
gammdatagamm1 <- gamm(ODBA_hour ~ s(Tidal_Phase,
                                    bs = "cr",
                                    k = 2),
                  correlation = corCAR1(value = acfsharkgammout),
                  random = list(ID = ~ODBA_hour),
                  data = GAMMdata,
                  na.action = na.omit)
# Error in solve.default(pdMatrix(a, factor = TRUE)): system is computationally singular: reciprocal condition number = 2.5695e-224
# In addition: Warning messages:
# 1: In smooth.construct.cr.smooth.spec(object, dk$data, dk$knots) : basis dimension, k, increased to minimum possible
# 2: In lme.formula(y ~ X - 1, random = rand, data = strip.offset(mf),: nlminb problem, convergence error code = 1. message = false convergence (8)

# It means your design matrix is not invertible and therefore can't be used to develop a regression model.
# This results from linearly dependent columns, i.e. strongly correlated variables.
# Examine the pairwise covariance (or correlation) of your variables to investigate if there are any variables that can potentially be removed.
# You're looking for covariances (or correlations) >> 0.
# Alternatively, you can probably automate this variable selection by using a forward stepwise regression.



summary(gammdatagamm1$gam)
gam.check(gammdatagamm1$gam)
