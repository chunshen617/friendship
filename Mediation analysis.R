library(mediation)

######Load data
#friendship
IDV <- read.csv("med_IDV_less5.csv")
#covariates
COV <- read.csv("med_cov_less5.csv")
COV$eth <- as.factor(COV$eth)
COV$site <- as.factor(COV$site)
COV$hand <- as.factor(COV$hand)
COV$mri <- as.factor(COV$mri)
#significant brain area
Area <- read.csv("med_sigArea_less5.csv")
AreaN <- colnames(Area)

#standardized
totArea <- data.frame(totArea = scale(rowSums(Area)))
df.tot <- cbind(scale(IDV),totArea,COV)

##########Model 1: close friend --> total area --> ADHD
model1.m <- lm(totArea ~ cf_log+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model1.m)
model1.y <- lm(ADHD_syn ~ cf_log+totArea+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model1.y)
model1.tot <- lm(ADHD_syn ~ cf_log+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model1.tot)
medM1_result <- mediate(model.m = model1.m,
                        model.y = model1.y,
                        sims = 5000,
                        boot = TRUE,
                        boot.ci.type = "bca",
                        mediator = "totArea",
                        treat = "cf_log")
t1 <- (medM1_result)

##########Model 2: close friend --> total area --> crystalized intelligence
model2.m <- lm(totArea ~ cf_log+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model2.m)
model2.y <- lm(CrysInt ~ cf_log+totArea+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model2.y)
model2.tot <- lm(CrysInt ~ cf_log++sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df.tot)
summary(model2.tot)
medM2_result <- mediate(model.m = model3.m,
                        model.y = model3.y,
                        sims = 5000,
                        boot = TRUE,
                        boot.ci.type = "bca",
                        mediator = "totArea",
                        treat = "cf_log")
t2 <- summary(medM2_result)

############Model1 of each region###########
model1.area.Mprop <- numeric(0)
model1.area.Mp <- numeric(0)
for (i in 1:ncol(Area)){
  AreaS <- data.frame(AreaS = scale(Area[,i]))
  df <- cbind(scale(IDV),AreaS,COV)
  model_m1 <- lm(AreaS ~ cf_log+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df)
  model_y1 <- lm(ADHD_syn ~ cf_log+AreaS+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df)
  med_result1 <- mediate(model.m = model_m1,
                         model.y = model_y1,
                         sims = 5000,
                         boot = TRUE,
                         boot.ci.type = "bca",
                         mediator = "AreaS",
                         treat = "cf_log")
  s <- summary(med_result1)
  model1.area.Mprop[i] <- s[["n0"]]
  model1.area.Mp[i] <- s[["n0.p"]]
}
sigArea_med_less5 <- data.frame(region = AreaN,
                                mediatedP = model1.area.Mprop,
                                pval = model1.area.Mp,
                                fdrp = p.adjust(model1.area.Mp, method = "fdr", n = length(model1.area.Mp)))
write.csv(sigArea_med_less5,"sigArea_medP_less5.csv")

############Model2 of each region###########
model2.area.Mprop <- numeric(0)
model2.area.Mp <- numeric(0)
for (i in 1:ncol(Area)){
  AreaS <- data.frame(AreaS = scale(Area[,i]))
  df <- cbind(scale(IDV),AreaS,COV)
  model_m3 <- lm(AreaS ~ cf_log+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df)
  model_y3 <- lm(CrysInt ~ cf_log+AreaS+sex+age+income+edu+bmi+pbs+eth+site+motion+hand+mri,data=df)
  med_result3 <- mediate(model.m = model_m3,
                         model.y = model_y3,
                         sims = 5000,
                         boot = TRUE,
                         boot.ci.type = "bca",
                         mediator = "AreaS",
                         treat = "cf_log")
  s <- summary(med_result3)
  model3.area.Mprop[i] <- s[["n0"]]
  model3.area.Mp[i] <- s[["n0.p"]]
}
sigArea_med3_less5 <- data.frame(region = AreaN,
                                 mediatedP = model3.area.Mprop,
                                 pval = model3.area.Mp,
                                 fdrp = p.adjust(model3.area.Mp, method = "fdr", n = length(model3.area.Mp)))
write.csv(sigArea_med3_less5,"sigArea_medP_cfTOcryInt_less5.csv")
