library(LogisticDx)
library(MASS)
library(ggplot2)
library(latex2exp)
library(MKmisc)
library(ResourceSelection)
library(car)
library(gridExtra)


################Correlation Analysis########################
cat("Before going into complex model building, looking at data relation is a sensible step to understand how your different variable interact together. Correlation look at trends shared between two variables
    
    Pearson coefficient -> it is the covariance of the two variable divided by the product of their variance, it is scaled between 1 (for a perfect positive correlation) to -1 (for a perfect negative correlation), 0 would be complete randomness.
    ")

fset <- read.csv(paste(getwd(),"/fullset.csv",sep = ""),header=T)

subset = fset[-c(22,23,24,40,41,42)]
subset2 = na.omit(subset)
purcalc = subset2[-c(1:3)]
#, fig.height = 100, fig.width = 100}

ggcorr(purcalc, label = TRUE, label_alpha = TRUE)

cat("Further analysis is needed for :
    Strong correlation between FSW and TPW (= 0.9). Only one variable should be included in the model.
    FNL is positively/negatively correlated with ST1 and ST2. 
    ")

# Plotting "FSW" (First Serve won) versus "TPW" (Total points won) 

qplot(FSW.1, 
      TPW.1, 
      data = purcalc, 
      geom = c("point", "smooth"), 
      method = "lm", 
      alpha = I(1 / 5), 
      se = FALSE)
# Correlation Coefficient
cor(purcalc$FSW.1, purcalc$TPW.1)
## [1] 0.8611495
cat("We see that there is a strong correlation betwen FSW and TPW. Henc we can remove TPW from the model.
    ")
# Plotting FNL1, ST1.1 and ST2.1

ggpairs(purcalc, 
        columns = c("FNL1", "ST1.1", "ST2.1"), 
        upper = list(continuous = wrap("cor", 
                                       size = 10)), 
        lower = list(continuous = "smooth"))

cat("We see that there is a strong correlation between FNL & ST1 and FNL & ST2 . Hence we are excluding FNL from the model fit. 
    Based on similar analysis for other variable correlations, we have decided to exclude the following variables from model fit.
    FNL (Final number of set won by player)
    TPW (Total points won by player)
    SSP (Second Serve percentage that was playable)
    NPA (Net Points Attempted by player)
    BPC (Break Points Created by player)
    NPW (Net Points Won by player)
    BPW (Break Points Won by player)
    ")

################Determining model with step procedure########################

###### Why Stepwise 

# Define full and null models and do step procedure
model.null = glm(Result~1,purcalc,family = binomial)
model.full = glm(Result~ST1.1+ST2.1+FSP.1+FSW.1+SSW.1+ST1.2+ST2.2+FSP.2+FSW.2+SSW.2+ACE.1+DBF.1+WNR.1+UFE.1+ACE.2+DBF.2+WNR.2+UFE.2,purcalc,family = binomial) 

u = model.full
par(mfrow=c(2,2))
cat("\n P-value for deviance tested \n ",1 - pchisq(u$null.deviance-u$deviance,u$df.null-u$df.residual),"\n")
cat("\n Deviance / df \n ",u$deviance/u$df.residual,"\n")
hoslem.test(u$y,  fitted(u),g=10)


cat("-----Forward----")
step(model.null,scope = list(upper=model.full),direction="forward",test="Chisq",data=purcalc, trace  = FALSE)
cat("-----Backward----")
step(model.full,scope = list(upper=model.null),direction="backward",test="Chisq",data=purcalc, trace  = FALSE)
cat("-----Both----")
step(model.null,scope = list(upper=model.full),direction="both",test="Chisq",data=purcalc, trace  = FALSE)

cat("After perfomring the forward, backward and stepwise regression, we found that forward regression is not efficient as the
    data has many variables. Backward and Stepwise regression has given a total of 9 same variables. 
    So, We chose backward regression for variable selection
    and built our final model using variables from backward regression.")
cat("--------How the backward Regression worked-------")
cat("Began with the full least squares model containing all n predictors, and then iteratively removes the least useful predictor, one-at-a-time.
    Start with all variables in the model.
    Remove the variable with the largest p-value(Pr(>Chi)) ==> the variable that is the least statistically significant.
    The new (n - 1)-variable model is t, and the variable with the largest p-value is removed.
    Continued until a stopping rule is reached. Finally all remaining variables have a significant p-value above significance threshold i.e,0.1.")

cat("But in this step process, there was one variable ACE.2 which has a p value sightly greater than 0.1 
    but wasn't removed during the process and was included in the model.
    ")

sub_fnl = glm(formula = Result ~ UFE.2 + SSW.2 + WNR.1 + FSW.2 + FSW.1 + 
                UFE.1 + WNR.2 + SSW.1 + ACE.2, family = binomial, data = purcalc)
final = glm(formula = Result ~ UFE.2 + SSW.2 + WNR.1 + FSW.2 + FSW.1 + 
              UFE.1 + WNR.2 + SSW.1, family = binomial, data = purcalc)

## $$H_0: \beta_{\text {ACE}_2} = 0 \text {vs} H_A:\beta_{\text {ACE}_2} \not = 0$$

par(mfrow=c(2,2))
cat("\n P-value for deviance tested \n ",1 - pchisq(sub_fnl$null.deviance-sub_fnl$deviance,sub_fnl$df.null-sub_fnl$df.residual),"\n")
cat("\n Deviance / df \n ",sub_fnl$deviance/sub_fnl$df.residual,"\n")
hoslem.test(sub_fnl$y,  fitted(sub_fnl),g=10)

cat("Deviace of model with ACE2",sub_fnl$deviance,"\n")
cat("Deviance of model without ACE2",final$deviance)

## After the test test stat = dev(final) - dev(sub_final) = 2.52 which not greater than chi sq, thus not significant.

CI.OR=exp(cbind(coef(sub_fnl) ,confint(sub_fnl, level=0.95)))######## CI FOR ODDS RATIO OF THE 
CI.OR
confint(sub_fnl)

################Model Adequecy######################
table = dx(final)
lp = log(table$P/(1-table$P))

p1 = ggplot(data = table ,aes(x = P, y = dr)) + geom_point(aes(color = y)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x = TeX("$\\hat{ \\pi}$"), y = TeX("$\\d_i$"), title = "Probability vs Deviance residual")

p2 = ggplot(table, aes(sample = dr))+stat_qq()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs( title = "Normal Probability of Deviance Residual", x = TeX("$\\d_i$"))

p3 = ggplot(data = table ,aes(x = lp, y = Pr)) + geom_point(aes(color = y)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                   panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x = TeX("$\\hat{ \\eta}$"), y = TeX("Pearson residual"), title = "Linear predictor vs Pearson residual") + geom_hline(yintercept = 3) + geom_hline(yintercept = -3)

grid.arrange(p1,p2,p3,ncol = 3)

newmod = cbind(purcalc,table)
subnew = subset(newmod, Pr < 3)
subnew = subset(subnew,Pr > -3)
xyz = subnew[-c(34:55)]

new_final = glm(formula = Result ~ UFE.2 + SSW.2 + WNR.1 + FSW.2 + FSW.1 + 
                  UFE.1 + WNR.2 + SSW.1, family = binomial,data = xyz ) # Model after dropping outliers

par(mfrow=c(2,2))
cat("\n P-value for deviance tested \n ",1 - pchisq(new_final$null.deviance-new_final$deviance,new_final$df.null-new_final$df.residual),"\n")
cat("\n Deviance / df \n ",new_final$deviance/new_final$df.residual,"\n")
hoslem.test(new_final$y,  fitted(new_final),g=10)


Anova(new_final, type="II", test="Wald")