install.packages("survival")
library("survival")
install.packages("rms")
library(rms)
install.packages("survminer")
library("survminer")

aids = read.csv("C:/UCONN/Semester 2/Data Analytics - R/Project/Dataset_1.csv")
View(aids)
#_________________________________________________________________________________________
# Data Cleaning

# removing rows with txgrp = 3 and txgrp = 4
aids = aids[aids$txgrp == 1 || aids$txgrp == 2]
View(aids)

# recoding variables
#1,2,3,4,5,6 in Raceth to ("WNH", "BNH", "H", "AI", "A", "U"). 
aids$raceth[aids$raceth == 1] = "WNH" 
aids$raceth[aids$raceth == 2] = "BNH" 
aids$raceth[aids$raceth == 3] = "H" 
aids$raceth[aids$raceth == 4] = "AI" 
aids$raceth[aids$raceth == 5] = "A" 
aids$raceth[aids$raceth == 6] = "U" 

#_________________________________________________________________________________________________

#Data Visualization
data1 = aids
mean(data1$time)
median(data1$time)

# Number of censored observations
sum(data1$censor)
# number of uncensored observations
length(data1$censor) - sum(data1$censor)

# CD4 stratus at screening > 50
sum(data1$strat2)
length(data1$strat2) - sum(data1$strat2)

#ivdrug
table(aids$ivdrug)
abline(h = mean(data1$age), col = "pink")
abline(v = mean(data1$time), col = "orange")
#plot between age and time
plot(data1$time, data1$age)
#________________________________________________________________________________________________

# Modeling

#kaplan meier estimation for tx
Kms1 = survfit(Surv(aids$time,aids$censor)~aids$tx, data = aids)
summary(Kms1)
ggsurvplot(Kms1,ylim = c(0.8,1), censor = F, legend = "right", title = "Kaplan Meier Curves", legend.title = "Treatment", legend.labs = c("Without IDV", "With IDV"))

#Showing that the two groups have a statistically significant differernce since P value < 0.05
survdiff(Surv(aids$time,aids$censor)~aids$tx, data = aids)

# building cox model with all variables
cox1 = coxph(Surv(aids$time,aids$censor)~tx+factor(sex)+factor(raceth)+factor(ivdrug)+
               hemophil+factor(karnof)+cd4+priorzdv+age, data = aids)
summary(cox1)

#Considering retained covariates - tx, karnoff, cd4, age for further models
cox2 = coxph(Surv(aids$time,aids$censor)~tx+factor(karnof)+cd4+age, data = aids)
summary(cox2)

#Checking if both models cox1 and cox2 are different with statistical significan : result :NO
anova(cox1, cox2)

#Plotting survival functions for both treatements seperately
ggsurvplot(survfit(cox2),ylim = c(0.9,1), censor = F, legend = "right", 
           title = "Cox Proportional Hazards", legend.title = "Treatment",
           legend.labs = c("Without IDV", "With IDV"))

#checking for non proportional hazards.

h = cox.zph(cox2)
h # testing with p values if p value < 0.05 there is a chance of non proportional haxards or
#interaction of co variates with time.

#Graphically checking for non prorortional hazards or interaction of covariates with time
ggcoxzph(h)

#Checking for non-linearity only for continuous variables using martingales residuals
res = resid(cox2,type='martingale')
aids2 = cbind(aids,res)
res_plot_cd4 = ggplot(aids2, aes(cd4,res))+geom_point(color = "orange")+geom_smooth(method = lm)+
  labs(title = "Martingale Residuals - CD4")
res_plot_age = ggplot(aids2, aes(age,res))+geom_point(color = "orange")+geom_smooth(method = lm)+
  labs(title = "Martingale Residuals - Age")
library(gridExtra)
grid.arrange(res_plot_age, res_plot_cd4, ncol = 2)

#both lines are straight and and a smooth fit aligns with 0 which suggests that there is 
#no non-linearty

# Checking for influential observations
ggcoxdiagnostics(cox2, type = "dfbeta",linear.predictions = FALSE)
i = 2/sqrt(length(aids$time))
i
#exclude if the dfbeta value is greater than i (2/sqrt(n)) for small datasets
#Comparing both the approaches based grouped by tx (three drug treatement and two drug treatment)
#Comparing two treatments 
cox3 = coxph(Surv(aids$time,aids$censor)~strata(tx)+karnof+cd4+age, data = aids)
summary(cox3)
coxdf = surv_summary(survfit(cox3))
Kms2 = survfit(Surv(aids$time,aids$censor)~strata(tx), data = aids)
summary(Kms2)
k1 = summary(Kms2)
kmsdf <- as.data.frame(k1[c("strata", "time", "n.risk", "n.event", "surv", "std.err", "lower",
                            "upper")])
#Building the plot to compare all the four curves. 
g = ggplot() 
coxgraph = g + geom_step(data = coxdf,  aes(x = time, y = surv, color = coxdf$strata, linetype = coxdf$strata), size = 0.75)

combined = coxgraph + geom_step(data = kmsdf, aes(x = time , y = surv, color = kmsdf$strata, linetype = kmsdf$strata),size = 0.75)

final = combined +labs(title = "COX Model vs Kaplan Meier\n", x = "Time", y = "Survival Probability") +
  scale_color_manual("Legend\n", labels = c("KM - without IDV ", "KM - with IDV", "Cox - without IDV", "Cox - with IDV"), values = c("blue", "red", "Green", "Orange")) +
  scale_linetype_manual("Legend\n", labels = c("KM - without IDV ", "KM - with IDV", "Cox - without IDV", "Cox - with IDV"), values = c(2,2,1,1)) +
  theme_grey() +
  theme(panel.background = element_rect(fill='lavender', colour='black'))+
  theme(axis.text.x=element_text(size=10), axis.title.x=element_text(size=14),
        axis.text.y=element_text(size=10), axis.title.y =element_text(size=14),
        plot.title=element_text(size=16, face="bold", color="darkblue"))
final
