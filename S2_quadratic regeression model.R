###########load data###########
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

#Data scale: the number of friends is log-transformed and mean centered
friend$cf_log <- log10(friend$totalCloseFriend+1)
friend$cf_log_c <- scale(friend$cf_log)

#prepare data frame for quadratic analysis
cfD <- subset(friend, select = c("cf_log_c"))
Dcf <- cbind(cfD,mh,cog,covs)
DVn <- colnames(Dcf[,2:31])

#quadratic regression function
quaR <- function(y,Data){
  lm(y ~ cf_log_c + I(cf_log_c^2) + 
       sex+age+income+edu+bmi+pbs+eth+site,data = Data)
}

M.qua <- apply(Dcf[,2:31], 2, function(x) quaR(x,Dcf))

#create result table: beta,95% cI,t,p
sumResult <- function(LMobject){
  LM.sum = lapply(LMobject, summary)
  LM.ci = lapply(LMobject, function(x) confint(x))
  coef.linear <- sapply(LM.sum,function(x) x[["coefficients"]][2,c(1,3,4)])
  rownames(coef.linear) <- c('beta_linear','t_linear','p_linear')
  coef.quadratic <- sapply(LM.sum,function(x) x[["coefficients"]][3,c(1,3,4)])
  rownames(coef.quadratic) <- c('beta_quadratic','t_quadratic','p_quadratic')
  coef.CI.linear <- sapply(LM.ci,function(x) c(linear.lowerCI = x[2,1], linear.upperCI = x[2,2]))
  coef.CI.quadratic <- sapply(LM.ci,function(x) c(quadratic.lowerCI = x[3,1], quadratic.upperCI = x[3,2]))
  tbl.result <- cbind(t(coef.linear),t(coef.quadratic),t(coef.CI.linear),t(coef.CI.quadratic))
  return(tbl.result)
}

result.qua.coef <- sumResult(M.qua)

#calculate delta quadratic R^2
#linear regression function define
lnR <- function(y,Data){
  lm(y ~ cf_log_c + 
       sex+age+income+edu+bmi+pbs+eth+site,data = Data)
}
M.ln <- apply(Dcf[,2:31], 2, function(x) lnR(x,Dcf))

covR <- function(y,Data){
  lm(y ~ sex+age+income+edu+bmi+pbs+eth+site,data = Data)
}
M.cov <- apply(Dcf[,2:31], 2, function(x) covR(x,Dcf))
                            
#create result table: r square, delta quadratic r square
Rsquare <- function(Mqua,Mln){
  Mqua.sum = lapply(Mqua, summary)
  Mln.sum = lapply(Mln, summary)
  Mqua.rsquares <- sapply(Mqua.sum,function(x) c(qua_r_sq = x$r.squared, qua_adj_r_sq = x$adj.r.squared))
  Mln.rsquares <- sapply(Mln.sum,function(x) c(ln_r_sq = x$r.squared, ln_adj_r_sq = x$adj.r.squared))
  delta.rsquares <- t(Mqua.rsquares)[,1] - t(Mln.rsquares)[,1]
  delta.adj.rsquares <- t(Mqua.rsquares)[,2] - t(Mln.rsquares)[,2]
  tbl.result <- cbind(t(Mqua.rsquares),t(Mln.rsquares),delta.rsquares,delta.adj.rsquares)
  return(tbl.result)
}

result.Rsquare.quadratic <- Rsquare(M.qua,M.ln)
result.Rsquare.linear <- Rsquare(M.ln,M.cov)
result.Rsquare <- data.frame(delta.adjR2.linear = result.Rsquare.linear[,6],
                             delta.adjR2.quadratic = result.Rsquare.quadratic[,6])

#save results                         
result.QR <- cbind.data.frame(result.qua.coef,result.Rsquare)
write.csv(result.baseline,file = "Result_beh_QR.csv")

#model comparison
F.modelComp <- list()
for (i in 1:length(DVn)){
  fit.qua = M.qua[[i]]
  fit.ln = M.ln[[i]]
  F.modelComp[[i]] = anova(fit.ln, fit.qua, test="F")
}

Fval <- sapply(F.modelComp, function(x) x[["F"]][2])
Pval <- sapply(F.modelComp, function(x) x[["Pr(>F)"]][2])
R.modelComp <- data.frame(DV = DVn,
                          Fvalue = Fval,
                          Pvalue = Pval,
                          bonf = as.factor(Pval<0.05/30))
write.csv(R.modelComp,"Result_cf_beh_modelComparison_F.csv")
