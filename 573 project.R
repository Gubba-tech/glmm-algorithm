library(ggplot2)
library(GGally)
library(reshape2)
library(Matrix)
library(lme4)
library(compiler)
library(parallel)
library(dplyr)
library(boot)
library(lattice)
library(viridisLite)
library(digest)
library(viridis)
library(sjPlot)
library(sjmisc)
library(glmmTMB)
library(caret)
library(MASS)
library(caTools)

hdp <- read.csv("https://stats.idre.ucla.edu/stat/data/hdp.csv")
str(hdp)
hdp <- within(hdp,{
  Married <- factor(Married, levels = 0:1, labels = c("no","yes"))
  DID <- factor(DID)
  HID <- factor(HID)
})
ggpairs(hdp[, c("IL6", "ntumors", "Experience")])
ggplot(hdp, aes(x = CancerStage, y = ntumors)) +
  stat_sum(aes(size = ..n.., group = 1)) +
  scale_size_area(max_size=10, name = "Population")


errbar_lims <- group_by(tmp, CancerStage) %>% 
  summarize(mean=mean(value), se=sd(value)/sqrt(nrow(hdp)), 
            upper=mean+(2*se), lower=mean-(2*se))

tmp <- melt(hdp[, c("CancerStage", "IL6")], id.vars="CancerStage")
ggplot(tmp, aes(x = CancerStage, y = value)) +
  geom_jitter(alpha = .1) +
  geom_point(data=errbar_lims, aes(x=CancerStage, y=mean)) +
  geom_violin(data = tmp, aes(x = CancerStage, y = value, 
                              fill = CancerStage, 
                              color=CancerStage), alpha = .75) +
  geom_boxplot(width=.1, outlier.shape=NA) +
  theme_minimal()+
  facet_grid(variable ~ .) +
  scale_y_sqrt()


tmp1 <- melt(hdp[, c("remission", "IL6", "ntumors", "Experience")],
            id.vars="remission")
ggplot(tmp1, aes(factor(remission), y = value, fill=factor(remission))) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free_y")

hdp <- hdp[1:500,]
accu1 <- rep(0, 10)
accu2 <- rep(0, 10)


for(i in 1:10) {
  set.seed(i+100)
  split <- sample.split(hdp$remission, SplitRatio = 0.8)
  hdp_train <- subset(hdp, split == TRUE)
  hdp_test <- subset(hdp, split == FALSE)
  
  fit2 <- glm(remission ~ IL6 + CancerStage + ntumors + Experience + DID, 
              data = hdp_train,
              family = "binomial")
  pred2 <- predict(fit2, newdata = hdp_test, type = "response")
  t1 = table(hdp_test$remission, pred2 > 0.5)
  accu1[i] = (t1[1,1] + t1[2,2])/sum(t1)
  
  fit1 <- glmer(remission ~ IL6 + CancerStage + ntumors + Experience +
                  (1 | DID), data = hdp_train, family = "binomial", 
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10)
  pred1 <- predict(fit1, newdata = hdp_test, type = "response")
  t2 <- table(hdp_test$remission, pred1 > 0.5)
  accu2[i] = (t2[1,1]+t2[2,2])/sum(t2)
  print(i)
}
sum(accu1)
sum(accu2)




plot_model(fit1, type = "re")
lattice::dotplot(ranef(fit1, which = "DID", condVar = TRUE), 
                 scales = list(y = list(alternating = 0)))
vcov(fit1)
print(fit1, corr = FALSE)
se <- sqrt(diag(vcov(fit1)))
tab <- cbind(Est = fixef(fit1), LL = fixef(fit1) - 1.96 * se, UL = fixef(fit1) + 1.96 *
               se)
#step_fit1 <- stepAIC(fit1, trace = FALSE)

tmpdat <- hdp[, c("IL6", "CancerStage", "ntumors", "Experience",
                  "DID")]
summary(hdp$ntumors)
jvalues <- with(hdp, seq(from = min(ntumors), to = max(ntumors), length.out = 100))
pp <- lapply(jvalues, function(j) {
  tmpdat$ntumors <- j
  predict(fit1, newdata = tmpdat, type = "response")
})

plotdat <- t(sapply(pp, function(x) {
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))
plotdat <- as.data.frame(cbind(plotdat, jvalues))
colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "ntumors")
head(plotdat)
ggplot(plotdat, aes(x = ntumors, y = PredictedProbability)) + geom_line() +
  ylim(c(0, 1))
ggplot(plotdat, aes(x = ntumors, y = PredictedProbability)) + 
  geom_linerange(aes(ymin = Lower, ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))


pp_stage <- lapply(levels(hdp$CancerStage), function(can_stage){
  tmpdat$CancerStage <- can_stage
  lapply(jvalues, function(j){
    tmpdat$ntumors <- j
    predict(fit1, newdata = tmpdat, type = "response")
  })
})

## pp of cancer stage and length of date 
plotdat2 <- (lapply(pp_stage, function(y){
  temp <- t(sapply(y, function(x){
    c(M = mean(x), quantile(x, c(0.25, 0.75)))
  }))
  temp <- as.data.frame(cbind(temp, jvalues))
  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "ntumors")
  return(temp)
}))
plotdat2 <- do.call(rbind, plotdat2)
plotdat2$CancerStage <- factor(rep(levels(hdp$CancerStage), each = length(jvalues)))
head(plotdat2)
ggplot(plotdat2, aes(x = ntumors, y = PredictedProbability)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = CancerStage), alpha = .15) +
  geom_line(aes(colour = CancerStage), size = 2) +
  ylim(c(0, 1)) + facet_wrap(~  CancerStage)


## pp of cancer stage only
pp_stageonly <- lapply(levels(hdp$CancerStage), function(stage1){
  tmpdat$CancerStage <- stage1
  predict(fit1, newdata = tmpdat, type = "response")
})
plotdat3 <- t(sapply(pp_stageonly, function(x){
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))
plotdat3 <- as.data.frame(cbind(plotdat3, levels(hdp$CancerStage)))
colnames(plotdat3) <- c("PredictedProbability", "Lower", "Upper", "CancerStage")
head(plotdat3)
ggplot(plotdat3, aes(x = CancerStage, y = PredictedProbability)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PredictedProbability), alpha = .25) +
  geom_boxplot() 



##
m3b <- glmer(remission ~  ntumors  + IL6  + CancerStage +
               Experience + (1 + ntumors | DID) + (1 | HID), data = hdp, family = binomial,
             nAGQ = 1)

lattice::dotplot(ranef(m3b, which = "DID", condVar = TRUE), 
                 scales = list(y = list(alternating = 0)))




