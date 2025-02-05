#
# Custom code for "Dose-Dependent Effects of Niacin Intake on Mortality: A U-Shaped Curve"
# 
# Author: Marc Ferrell
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Set working directory to Project Folder
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("path/to/data")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Libraries
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
library(readr)
library(tidyverse)
library(survival)
library(MASS)
library(weights)
library(rms)


###################################
#
# Custom Functions
#
###################################

##### Spline model of Cox residuals to adjust for known confounders

## Returns Spline curves, CIs, and X coordinates with HR significantly different from 1

cox_spline <- function(t, v, x, covar, knots = 3, w = rep(1,length(x)),
                       ref = median(x, na.rm = TRUE)){
  require(rms)
  xx <- x
  dd <- datadist(xx)
  options(datadist='dd')
  
  # Unadjusted
  spl_mod_un <- cph(Surv(t, v) ~ rcs(xx, knots), weights = w)
  print(anova(spl_mod_un))
  
  # Move HRs to make ref have a HR of 1
  p_un <- Predict(spl_mod_un, xx, np = 1e4)
  ref_un <- p_un$yhat[which.min(abs(p_un$xx - ref))]
  p_un$yhat <- p_un$yhat - ref_un
  p_un$lower <- p_un$lower - ref_un
  p_un$upper <- p_un$upper - ref_un
  
  ## Risk regions -- yhat in ln(HR)
  
  hi <- pos_neg_Regions(p_un$xx, p_un$lower)
  hi <- hi[hi$sign == "+",]
  lo <- pos_neg_Regions(p_un$xx, p_un$upper)
  lo <- lo[lo$sign == "-",]
  
  r_un <- rbind(hi[,1:2], lo[,1:2])
  r_un$Risk <- c(rep("High", nrow(hi)), rep("Low", nrow(lo)))
  r_un <- r_un[order(r_un$start),]
  
  # Adjusted
  
  # Calculate Residuals from covariates
  cox_mod <- coxph(Surv(t, v) ~ covar, weights = w)
  
  # Spline model of residual risk
  res_mod <- glm(cox_mod$residuals ~ rcs(xx, knots), weights = w,
                 x=TRUE, y=TRUE)
  
  print(summary(res_mod))
  print(summary(aov(res_mod)))
  
  # Extract data points of adjusted curve
  p <- predict(res_mod, newdata = data.frame(xx), se.fit = TRUE)
  
  # Move HRs to make ref have a HR of 1
  ref_ad <- p$fit[which.min(abs(xx - ref))]
  
  p_ad <- data.frame(x = xx, y = p$fit - ref_ad, 
                     lower = p$fit - 1.96 * p$se.fit - ref_ad,
                     upper = p$fit + 1.96 * p$se.fit - ref_ad)
  
  p_ad <- p_ad[order(p_ad$x),]
  
  ## Risk regions -- yhat in ln(HR)
  
  hi <- pos_neg_Regions(p_ad$x, p_ad$lower)
  hi <- hi[hi$sign == "+",]
  lo <- pos_neg_Regions(p_ad$x, p_ad$upper)
  lo <- lo[lo$sign == "-",]
  
  r_ad <- rbind(hi[,1:2], lo[,1:2])
  r_ad$Risk <- c(rep("High", nrow(hi)), rep("Low", nrow(lo)))
  r_ad <- r_ad[order(r_ad$start),]
  
  (list(
    Unadjusted =      p_un,
    Risk_Unadjusted = r_un,
    Adjusted   =      p_ad,
    Risk_Adjusted =   r_ad
  ))
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Load Data
#   * dat - all observations
#   * covar1 - Matrix of covariates for Model 1
#   * covar3 - Matrix of covariates for Model 2
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("nhanes-ndi-imputed-092924.RData")
attach(nhanes_ndi_merge_imputed)

dat$DSQTNIACcat100 <- ifelse(dat$DSQTNIAC == 0, 0, 
                             ifelse(dat$DSQTNIAC < 100, 1, 2))

# write.table(dat[,c("UNNIAC", "UNNIACcat", "DSQTNIAC", "DSQTNIACcat", "TOTNIAC", 
#                 "TOTNIACcat", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", 
#                 "HTN","EverSmoker", "ALQ130", "INDFMPIR", "HXCVD", "HXDiabetes",
#                 "CKD", "permth_int", "mortstat")],
#             file = "NHANES-NDI-29627.csv", quote=FALSE, row.names=FALSE, sep=",")

###################################
#
# Proportional Hazard Assumption Testing
# Plotted lines are flat
#
###################################

## Covariates: Age, Sex, Ethnicity, Education, Hypertension, Smoking, Alcohol Use,
##             Income:Poverty Ratio, History of: CVD, Diabetes, CKD

plot(
  cox.zph(
    robcov(
      coxph(Surv(dat$permth_int, dat$mortstat) ~ dat$RIDAGEYR + dat$RIAGENDR + 
              dat$RIDRETH1 + dat$DMDEDUC2 + dat$HTN + dat$EverSmoker + dat$ALQ130 + 
              dat$INDFMPIR + dat$HXCVD + dat$HXDiabetes + dat$CKD, 
            weights = dat$WTNDI10YR)
    )
  )
)

## Predictors: Niacin from Food (Contiuous  and Categorical), 
##             Niacin from Supplements (Contiuous  and Categorical),
##             Total Niacin Intake (Contiuous  and Categorical)

for(predictor in list(dat$UNNIAC, dat$UNIACcat, dat$DSQTNIAC, dat$DSQTNIACcat, dat$TOTNIAC, dat$TOTNIACcat)){
  
  plot(cox.zph(robcov(
      coxph(Surv(dat$permth_int, dat$mortstat) ~ dat$predictor, 
            weights = dat$WTNDI10YR)
    )))
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 1
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### F1A - Niacin Intake in Supplement USers and Non-users

F1A <- ggplot(dat, aes(x = as.factor(User), y = UNNIAC)) + 
  geom_boxplot(outlier.shape = NA, aes(weight = WTNDI10YR)) +
  coord_cartesian(ylim = c(0,100))

F1Aa <- ggplot(dat, aes(x = as.factor(User), y = TOTNIAC)) + 
  geom_boxplot(outlier.shape = NA, aes(weight = WTNDI10YR)) +
  coord_cartesian(ylim = c(0,100))

print(F1A)
print(F1Aa)

wilcox.test(UNNIAC ~ User, data = dat)$p.value # P < 1.1e-24

nrow(dat[dat$User == 1,])


#### F1B - Niacin deficiency and >UL rates in supplement users and non-users

### 31.5% Users

wpct(dat$User == 1, weight = dat$WTNDI10YR)

### <EAR Incidence 29.0% Nonuser / 3.4% User
### >UL  Incidence 14.1% Nonuser / 64.1% User

wpct(dat$TOTNIAC[dat$User == 0] < 15, weight = dat$WTNDI10YR[dat$User == 0])
wpct(dat$TOTNIAC[dat$User == 1] < 15, weight = dat$WTNDI10YR[dat$User == 1])

wpct(dat$TOTNIAC[dat$User == 0] > 35, weight = dat$WTNDI10YR[dat$User == 0])
wpct(dat$TOTNIAC[dat$User == 1] > 35, weight = dat$WTNDI10YR[dat$User == 1])

### Chi Square tests

wtd.chi.sq(as.numeric(dat$TOTNIAC < 15), 
           as.numeric(dat$User == 0),
           weight = dat$WTNDI10YR)

wtd.chi.sq(as.numeric(dat$TOTNIAC > 35), 
           as.numeric(dat$User == 0),
           weight = dat$WTNDI10YR)


## F1C - Diet / Supp Distribution within Total Niacin

F1B <- ggplot(dat[dat$User == 1,], aes(y=Niac_supp_percent, x= TOTNIACcat)) +
  geom_boxplot(outlier.shape = NA, aes(weight = WTNDI10YR))

kruskal.test(Niac_supp_percent ~ TOTNIACcat, data = dat)$p.value # P < 2.2e-16

nrow(dat[dat$User == 1 & dat$TOTNIACcat == ">100",])

## F1D - Diet / Supp Distribution within Total Niacin

### Model 1

#### All intake groups

summary(robcov(coxph(
  Surv(dat$permth_int, dat$mortstat) ~ dat$User + covar1, 
  weights = dat$WTNDI10YR
)))

summary(robcov(coxph(
  Surv(dat$permth_int, dat$mortstat) ~ dat$User + covar3, 
  weights = dat$WTNDI10YR
)))



#### Intake from Food < EAR (15mg/d), >UL (35mg/d), Ref

for(intake in c("Low", "Med", "High")){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$User + covar1, 
    subset = dat$UNNIACcat == intake, weights = dat$WTNDI10YR
  )))
}

### Model 2

#### All intake groups

summary(robcov(coxph(
  Surv(dat$permth_int, dat$mortstat) ~ dat$User + covar3, 
  weights = dat$WTNDI10YR
)))

#### Intake from Food < EAR (15mg/d), >UL (35mg/d), Ref

for(intake in c("Low", "Med", "High")){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$User + covar3, 
    subset = dat$UNNIACcat == intake, weights = dat$WTNDI10YR
  )))
}

### Trend Test (covar3)

# Agresti A. Categorical data analysis. 2nd ed. New York: John Wiley & Sons; 2002.

my_coef = c(-0.0972081, -0.1051189, -0.1508957)

w = tapply(dat$TOTNIAC, dat$UNNIACcat, length)

# Ref group 15-35 coded as zero for coxph function 
summary(lm(my_coef ~ c(1, 0, 2), weights = w))
# P = 0.622

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 2
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## F2A - Total Niacin Intake for All-Comers: Cox HRs

for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$TOTNIACcat + model, 
         weights = dat$WTNDI10YR
  )))
}

### Trend Test (covar3)

# Agresti A. Categorical data analysis. 2nd ed. New York: John Wiley & Sons; 2002.

my_coef = c(0, 0.1723613, -0.0917559, -0.1834883)

w = tapply(dat$TOTNIAC, dat$TOTNIACcat, length)

# Ref group 15-35 coded as zero for coxph function 
summary(lm(my_coef ~ c(1, 0, 2, 3), weights = w))


for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$TOTNIACcat + model, 
    weights = dat$WTNDI10YR
  )))
}

## F2B - Total Niacin Intake for All-Comers: Spline Analysis (See custom fxns)

x <- dat[, "TOTNIAC"]
xx <- x
dd <- datadist(xx)
options(datadist = "dd")

f2b <- cox_spline(t = dat$permth_int, v = dat$mortstat, x = dat$TOTNIAC,
                  w = dat$WTNDI10YR,
                  covar = covar3, ref = 15)
print(
ggplot(f2b$Adjusted, 
       aes(x=x, y=exp(y))) + 
  geom_line() +
  geom_line( aes(y=exp(lower), linetype = "dashed")) +
  geom_line( aes(y=exp(upper), linetype = "dashed"))
)

## F2C - Total Niacin Intake for Non-users: Cox HRs

for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$TOTNIACcat + model, 
    weights = dat$WTNDI10YR, subset = dat$User == 0
  )))
}

### Trend Test (covar3)

# Agresti A. Categorical data analysis. 2nd ed. New York: John Wiley & Sons; 2002.

my_coef = c(0, 0.1946102, -0.0005053, 0.2852331)

w = tapply(dat$TOTNIAC, dat$TOTNIACcat, length)

# Ref group 15-35 coded as zero for coxph function 
summary(lm(my_coef ~ c(1, 0, 2, 3), weights = w))

## F2D - Total Niacin Intake for Non-users: Spline Analysis (See custom fxns)
x <- dat[dat$User == 0, "TOTNIAC"]
xx <- x
dd <- datadist(xx)
options(datadist = "dd")

f2d <- cox_spline(t = dat$permth_int[dat$User == 0], v = dat$mortstat[dat$User == 0],
                  w = dat$WTNDI10YR[dat$User == 0],
                  x = dat$TOTNIAC[dat$User == 0], covar = covar3[dat$User == 0,], ref = 15)

print(
  ggplot(f2d$Adjusted, 
         aes(x=x, y=exp(y))) + 
    geom_line() +
    geom_line( aes(y=exp(lower), linetype = "dashed")) +
    geom_line( aes(y=exp(upper), linetype = "dashed"))
)

## F2E - Total Niacin Intake for Users: Cox HRs

for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$TOTNIACcat + model, 
    weights = dat$WTNDI10YR, subset = dat$User == 1
  )))
}

### Trend Test (covar3)

# Agresti A. Categorical data analysis. 2nd ed. New York: John Wiley & Sons; 2002.

my_coef = c(0, 2.062e-01, -1.577e-01, -2.317e-01)

w = tapply(dat$TOTNIAC, dat$TOTNIACcat, length)

# Ref group 15-35 coded as zero for coxph function 
summary(lm(my_coef ~ c(1, 0, 2, 3), weights = w))

## F2F - Total Niacin Intake for Non-users: Spline Analysis (See custom fxns)

x <- dat[dat$User == 1, "TOTNIAC"]
xx <- x
dd <- datadist(xx)
options(datadist = "dd")

f2f <- cox_spline(t = dat$permth_int[dat$User == 1], v = dat$mortstat[dat$User == 1],
                  w = dat$WTNDI10YR[dat$User == 1],
                  x = dat$TOTNIAC[dat$User == 1], covar = covar3[dat$User == 1,], ref = 15)

print(
  ggplot(f2f$Adjusted, 
         aes(x=x, y=exp(y))) + 
    geom_line() +
    geom_line( aes(y=exp(lower), linetype = "dashed")) +
    geom_line( aes(y=exp(upper), linetype = "dashed"))
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Supplemental Figure 2
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## SF2A - Niacin Intake and Follow-up time for those alive and deceased at last follow-up

SF2A <- ggplot(dat[dat$User == 0,], 
              aes(x=TOTNIAC, y = permth_int, color = as.factor(mortstat))) + 
  geom_point() +
  facet_grid(mortstat ~ .)

print(SF2A)

## SF2B - Total Niacin Intake for Non-users: Spline Analysis (See custom fxns)

####* With high niacin intake outliers (>75 mg/d) excluded

x <- dat[dat$User == 0 & dat$TOTNIAC < 75, "TOTNIAC"]
xx <- x
dd <- datadist(xx)
options(datadist = "dd")

sf2b <- cox_spline(t = dat$permth_int[dat$User == 0 & dat$TOTNIAC <= 75], 
                  v = dat$mortstat[dat$User == 0 & dat$TOTNIAC <= 75],
                  x = dat$TOTNIAC[dat$User == 0 & dat$TOTNIAC <= 75], 
                  w = dat$WTNDI10YR[dat$User == 0 & dat$TOTNIAC < 75],
                  covar = covar3[dat$User == 0 & dat$TOTNIAC <= 75,], ref = 15)


ggplot(sf2b$Adjusted, 
       aes(x=x, y=exp(y))) + 
  geom_line() +
  geom_line( aes(y=exp(lower), linetype = "dashed")) +
  geom_line( aes(y=exp(upper), linetype = "dashed"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## F3A - Histogram of NHANES Participant Niacin Supplement Doses

# Nonusers, 21421
table(dat$User)

# Highest bar, 7905
nrow(dat[dat$User == 1 & dat$DSQTNIAC < 100,])

my_hist <- ggplot(dat[dat$User == 1,], 
                  aes(x=DSQTNIAC, y=after_stat(ncount) * 7905, weight = WTNDI10YR)) + 
  geom_histogram(binwidth = 100) 

print(my_hist)

# Highest Bar, 2734
nrow(dat[dat$User == 1 & dat$DSQTNIAC == 20,])

my_hist_inset <- ggplot(dat[dat$User == 1,], 
                        aes(x=DSQTNIAC, y=after_stat(ncount) * 2734, weight = WTNDI10YR)) + 
  geom_histogram(breaks = 1:100) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0, 3000))

print(my_hist_inset)

## F3B - Stratification by Supplement Dose Categories

for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ dat$DSQTNIACcat + model, 
    weights = dat$WTNDI10YR
  )))
}

for(model in list(covar1, covar3)){
  summary(robcov(coxph(
    Surv(dat$permth_int, dat$mortstat) ~ factor(dat$DSQTNIACcat100) + model, 
    weights = dat$WTNDI10YR
  )))
}

### Trend Test (covar3)

# Agresti A. Categorical data analysis. 2nd ed. New York: John Wiley & Sons; 2002.

my_coef = c(0, -0.2738666, -0.0073975, -0.1061163, -0.2896009)

w = tapply(dat$TOTNIAC, dat$DSQTNIACcat, length)

# Ref group 15-35 coded as zero for coxph function 
summary(lm(my_coef ~ c(0:4), weights = w))

## F3C - Niacin Supplement Dose: Spline Analysis (See custom fxns)

x <- dat[, "DSQTNIAC"]
xx <- x
dd <- datadist(xx)
options(datadist = "dd")

f3c <- cox_spline(t = dat$permth_int, v = dat$mortstat, w = dat$WTNDI10YR,
                  x = dat$DSQTNIAC, covar = covar3, ref = 0)

ggplot(f3c$Adjusted, 
       aes(x=x, y=exp(y))) + 
  geom_line() +
  geom_line( aes(y=exp(lower), linetype = "dashed")) +
  geom_line( aes(y=exp(upper), linetype = "dashed"))

## F3D - Low dose Niacin versus Non-users, Stratified by Chronic disease
for(i in c(1:2)){
  for(disease in list("CKD", "Dyslipid", "HXCVD", "HXDiabetes")){
    summary(robcov(coxph(
      Surv(dat$permth_int, dat$mortstat) ~ dat$DSQTNIACcat + covar3, 
      weights = dat$WTNDI10YR, subset = dat[,disease] == i
    )))
  }
}

# Having any of the above: Yes "Metabolic Disease" 
summary(robcov(coxph(
  Surv(dat$permth_int, dat$mortstat) ~ dat$DSQTNIACcat + covar3, 
  weights = dat$WTNDI10YR, subset = dat$CKD == 1 | dat$HXCVD == 1 | dat$HXDiabetes == 1
)))

# Having none of the above: No "Metabolic Disease" 
summary(robcov(coxph(
  Surv(dat$permth_int, dat$mortstat) ~ dat$DSQTNIACcat + covar3, 
  weights = dat$WTNDI10YR, subset = dat$CKD == 0 & dat$HXCVD == 0 & dat$HXDiabetes == 0
)))
  



