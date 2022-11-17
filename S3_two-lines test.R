source("http://webstimate.org/twolines/twolines.R")

library(tidyverse)
library(scales)

#load data
#friendship
friend <- read.csv('friend_baseline.csv')
#mental health
mh <- read.csv('cbcl_baseline.csv')
#cognition
cog <- read.csv('tbss_baseline.csv')
#covariates
covs <- read.csv('cov_beh_baseline.csv')
covs$eth <- as.factor(covs$eth)
covs$site <- as.factor(covs$site)

#the number of friends is log-transformed
friend$cf_log <- log10(friend$totalCloseFriend+1)
cfD <- subset(friend, select = c("cf_log"))

#select significant nonlinear DVs
result.QR <- read.csv("Result_cf_mh&cog_quadraticR_baseline.csv")
DV_raw <- cbind(mh_s,cog_s)
DV_s <- DV_raw[,which(result.QR$p_quadratic<0.05/60)]
DVn <- colnames(DV_s)
Dcf <- cbind(cfD,DV_s,covs) 

#two-lines test, find the break point
#formula
tl.fmla <- list()
for (i in 1:length(DVn)){
  tl.fmla[i] = paste0(DVn[i],"~cf_log+sex+age+income+edu+bmi+pbs+eth+site")
}

tl.cf <- list()
for (i in 1:length(DVn)) {
  tl.cf[[i]] <- twolines(as.formula(tl.fmla[[i]]), data = Dcf, graph=0)
}

#break point
cf.bp <- sapply(tl.cf,function(x) c(breakPoint = x[["xc"]]))
Result.cf.bp <- data.frame(DV = c(DVn),
                           BreakPoint = c(cf.bp),
                           BreakPoint_raw = c(10^cf.bp-1))
rownames(Result.cf.bp) <- NULL

#run lm seperately
result_2lm <- list()
for (i in 1:length(DVn)){
  #group 1
  Dg1 <- subset(Dcf, cf_log <= Result.cf.bp[i,2])
  lm.g1 <- lm(as.formula(tl.fmla[[i]]),data = Dg1)
  n.g1 <- nrow(Dg1)
  b1 <- summary(lm.g1)[["coefficients"]][2,1]
  p1 <- summary(lm.g1)[["coefficients"]][2,4]
  ci.g1 <- confint(lm.g1)
  b1.lower <- ci.g1[2,1]
  b1.upper <- ci.g1[2,2]
  #group 2
  Dg2 <- subset(Dcf, cf_log > Result.cf.bp[i,2])
  lm.g2 <- lm(as.formula(tl.fmla[[i]]),data = Dg2)
  n.g2 <- nrow(Dg2)
  b2 <- summary(lm.g2)[["coefficients"]][2,1]
  p2 <- summary(lm.g2)[["coefficients"]][2,4]
  ci.g2 <- confint(lm.g2)
  b2.lower <- ci.g2[2,1]
  b2.upper <- ci.g2[2,2]
  #result table
  result <- data.frame(DV = c(DVn[i]),coef_g1 = b1,coef_lower_g1 = b1.lower,coef_upper_g1 = b1.upper,pval_g1 = p1,n_g1 = c(n.g1),
                       coef_g2 = b2, coef_lower_g2 = b2.lower,coef_upper_g2 = b2.upper,pval_g2 = p2,n_g2 = c(n.g2))
  result_2lm <- rbind(result_2lm,result)
}

result_2lt <- cbind(Result.cf.bp,result_2lm[,2:12])
write.csv(result_2lt,file = "Result_cf_mh&cog_2linesTest_baseline.csv")

