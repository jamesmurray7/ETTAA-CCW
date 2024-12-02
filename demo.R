# /////////////////////////////////////////////////////////////////////////////
# demo.R                                                                      /
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -/
# Demonstration use of functions for "Trial emulation to assess the effect of /
# surgery on aneurysm-related survival when there are competing risks in      /
# patients with thoracic aortic aneurysm"                                     /
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -/
# NB: The results don't align with the main paper, due to the simulated data  /
# set not being created in a particularly clever way. Serves as illustration  /
# only.                                                                       /
# /////////////////////////////////////////////////////////////////////////////

rm(list=ls())
dat <- readRDS('./simData.RDS') # Assumes working directory is set correctly.
source("functions.R")

# Bootstrapping requires the dataframe of AR status `AR.df` is defined.
AR.df <- subset(dat, select = c("ssid", "final_ar_status"))

# The workhorse function is fixITB.
args(fixITB)

# Ensure that formatData is tailored to use case (i.e. to the variable names in
# data).
args(formatData)

# By default this performs the 'base case' analysis, here owing to the simulated
# data having different distribution of event times, alter it to be 'capped' at
# 5 years.
fit1 <- fixITB(dat, truncation.time = 365.25 * 5,
               rmsts = seq(1, 5, 1) * 365.25,
               verbose = TRUE)
# Print the output, S3 method
fit1

# Can make a quick plot via S3 method, too
args(plot.fixITB)
plot(fit1, 'KM', legend.pos = 'bottom')
# Plot cumulative incidence
par(mfrow = c(1,2))
plot(fit1, 'CI')
# -"- of the competing risks, too
plot(fit1, 'CR')
par(mfrow=c(1,1))

# Note you can get _any_ S(t) you want for each arm from the fitted objects,
# since the original survival objects are stored in the list.
summary(fit1$noCR$S.surgery, times = 365.25 * 4) # if you were interested in e.g. 4 year survival
summary(fit1$noCR$S.control, times = 365.25 * 4) # Shows just S(t|arm)
summary(fit1$CR$sCR, times = 365.25 * 4)         # Shows Probabilities of each state for each arm.


# If we're only interested in the composite endpoint (or only competing risks),
# we can specify what output we do want
fit2 <- fixITB(dat, outputs = 'composite', , truncation.time = 365.25 * 5,
               rmsts = seq(2, 4, .5) * 365.25)
fit2 # Just prints the composite survival outcome stuff

fit3 <- fixITB(dat, outputs = 'CR', truncation.time = 365.25 * 5,
               rmsts = seq(2, 4, .5) * 365.25)
fit3
# And we can change arguments to the S3 plotting method, e.g. of an ugly plot:
plot(fit3, surgery.col = 'magenta', control.col = 'tomato',
     aneurysm.lty = 5, GP.lty = 3, GP.col = 'grey20')

