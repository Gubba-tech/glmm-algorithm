Doc1 <- hdp[hdp[,22]==1,]
Doc2 <- hdp[hdp[,22]==2,]
Doc3 <- hdp[hdp[,22]==3,]
Doc4 <- hdp[hdp[,22]==4,]
Doc5 <- hdp[hdp[,22]==5,]
Doc6 <- hdp[hdp[,22]==6,]
Doc7 <- hdp[hdp[,22]==7,]
Doc8 <- hdp[hdp[,22]==8,]
Doc11 <- hdp[hdp[,22]==11,]
Doc14 <- hdp[hdp[,22]==14,]
library(sjPlot)
library(sjmisc)
library(ggeffects)

fit_doc1 <- glm(remission ~IL6 + CRP + CancerStage + LengthofStay + Experience, 
                family = "binomial", data = Doc1)
plot(fit_doc1)

hdp_test1 <- rbind(Doc1, Doc2, Doc3, Doc4, Doc5, Doc6, Doc7, Doc8, Doc11, Doc14)
m = glmer(remission ~ IL6 + CancerStage + ntumors + Experience + (1 + ntumors | DID),
          family = binomial, 
          control = glmerControl(optimizer = "bobyqa"), 
          data=hdp_test1)
plot_model(m, type = "re")
lattice::dotplot(ranef(m, which = "DID", condVar = TRUE), 
                 scales = list(y = list(alternating = 0)))

data(efc)
fit <- glm(remission ~IL6 + CRP + CancerStage + LengthofStay + Experience, 
           data = Doc11, 
           family = binomial(link = "logit"))

plot(fit)

ggplot(Doc2) + plot_model(aes(x = CRP, y = remission, color = DID ))
