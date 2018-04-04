######################
### Power analysis ###
######################

### script setup ####

# we want to try and do a power analysis for our various analyses, 
# to indicate whether lack of interactions between effects is likely a power issue or a genuine lack of effect

## set up the script

rm(list=ls())

# make a list of required packages, check they're all installed, then load them up

j <- c("rstudioapi","ggplot2","pbapply")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://star-www.st-andrews.ac.uk/cran/")

lapply(j, require, character.only = TRUE)

devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc)
library(lme4)

# also source my own version of the sim.glmm function which I've tried to repair
source("sim.glmm.R")

# set working directory to current script location and load up data (from a folder named "Data" in the same location)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### power analysis of pollination success (binomial) ####


# import the true data to check parameters against
dframe1 <- read.csv("Data\\SeedSetBinom.csv")

dframe1$fPlot <- factor(dframe1$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1$fRound <- factor(dframe1$Round)
dframe1$fPlantNo <- factor(dframe1$PlantNo)
dframe1$fDistance <- factor(dframe1$Distance) # make "Distance" available as a factor if reqd

dframe1$Regime <- relevel(dframe1$Regime,"Control") # relevel variables related to control level
dframe1$Pollinators <- relevel(dframe1$Pollinators,"Control")
dframe1$LitUnlit <- relevel(dframe1$LitUnlit, "Lit")

# usage variable - first sort dataframe by fPlantNo
dframe1 <- dframe1[order(dframe1$fPlantNo),]
# now run down counting down sequences of fPlantNo
dframe1$usage <- sequence(rle(as.character(dframe1$fPlantNo))$lengths)

summary(dframe1)


with(dframe1,table(Plot,Round)) # plants per plot and round - all plants included as they should be


#### first, simulate a single data set - just the basic experimental structure

# to do this we construct the model we want to test, based on the real data (i.e. the actual model!)

model1 <- glmer(cbind(Successes,Failures) ~ Treatment + Pollinators + Distance + usage
                +(1|fRound),
                family = binomial,
                data = dframe1)

summary(model1)
drop1(model1, test = "Chi")


# then feed this model into sim.glmm - which pulls out all the effect sizes it needs from the model to simulate the output

testdata <- sim.glmm(model1)

# also check it works with my tweaked version

testdata_CJM <- sim.glmm_CJM(model1) # it does...


# plot lighting effect against variance factors
plot(response/n ~ Treatment, data = testdata)
plot(response/n ~ Pollinators, data = testdata)
plot(response/n ~ Distance, data = testdata)
plot(response/n ~ usage, data = testdata)
plot(response/n ~ fRound, data = testdata)


boxplot(response/n ~ Treatment, data = testdata, border = "white")
points(response/n ~ Treatment, data = testdata, pch = 1)

boxplot(response/n ~ Pollinators, data = testdata, border = "white")
points(response/n ~ Pollinators, data = testdata, pch = 1)



# now test for variation in pollination between lit & unlit

fit <- glmer(cbind(response, n - response) ~ 
               Treatment + Pollinators + Distance + usage +
               (1|fRound),
             family = binomial,
             data = testdata)

summary(fit)
drop1(fit, test = "Chi")


# note that although this model is based on data directly simulated from the previous one, it doesn't produce the same significance level
# (which is exactly what we wanted to happen!)


### now functionalise this simulation

# to actually feed this straight into the analysis...
fit <- glmer(cbind(response, n - response) ~ 
               Treatment + Pollinators + Distance + usage +
               (1|fRound),
             family = binomial,
             data = sim.glmm(model1))

summary(fit)
drop1(fit, test = "Chi")



# ...and now the function
sim.test.pval <- function(...){
  fit <- glmer(cbind(response, n - response) ~ 
                 Treatment + Pollinators + Distance + usage +
                 (1|fRound),
               family = binomial,
               data = sim.glmm(model1))
  
  drop1(fit, test= "Chi")["Treatment","Pr(Chi)"]
}

sim.test.pval()



# now we can simulate this several times to estimate the power
# in the final reckoning we'll want to do 1000 sims
sim.pvals <- sapply(1:20, sim.test.pval)
mean(sim.pvals < 0.05)  # 0.8 means we can detect 80% of real significant effects of this size

# show confints - will be wide because number of sims is low
binom.test(table(factor(sim.pvals < 0.05, c(T, F))))$conf.int



### right, that's the basis!
## now we want to (a) set it up to simulate different models (including different combinations of interaction terms)
## (b) do it over 1000 simulations for each

### (a) set up the various functions

# (i) basic model (as above) - no interactions

sim.test.pval.i <- function(...){
  fit <- glmer(cbind(response, n - response) ~ 
                 Treatment + Pollinators + Distance + usage +
                 (1|fRound),
               family = binomial,
               data = sim.glmm(model1))
  
  drop1(fit, test= "Chi")["Treatment","Pr(Chi)"]
}

sim.test.pval.i()

# (ii) both two-way interactions simultaneously

model2 <- glmer(cbind(Successes,Failures) ~ Treatment * Pollinators + Treatment * Distance + usage
                +(1|fRound),
                family = binomial,
                data = dframe1)

summary(model2)
drop1(model2, test = "Chi")
drop1(model2, test = "Chi")[c(3:4),4]


sim.test.pval.ii <- function(...){
  fit <- glmer(cbind(response, n - response) ~ 
                 Treatment * Pollinators + Treatment * Distance + usage +
                 (1|fRound),
               family = binomial,
               data = sim.glmm(model2))
  
  drop1(fit, test= "Chi")[c(3:4),4]
}

sim.test.pval.ii()

# (iii) Treatment:Pollinators interaction

model3 <- glmer(cbind(Successes,Failures) ~ Treatment * Pollinators + Distance + usage
                +(1|fRound),
                family = binomial,
                data = dframe1)

summary(model3)
drop1(model3, test = "Chi")
drop1(model3, test = "Chi")[4,4]


sim.test.pval.iii <- function(...){
  fit <- glmer(cbind(response, n - response) ~ 
                 Treatment * Pollinators + Distance + usage +
                 (1|fRound),
               family = binomial,
               data = sim.glmm(model3))
  
  drop1(fit, test= "Chi")[4,4]
}

sim.test.pval.iii()


# (iv) Treatment:Distance interaction

model4 <- glmer(cbind(Successes,Failures) ~ Treatment * Distance + Pollinators + usage
                +(1|fRound),
                family = binomial,
                data = dframe1)

summary(model4)
drop1(model4, test = "Chi")
drop1(model4, test = "Chi")[4,4]


sim.test.pval.iv <- function(...){
  fit <- glmer(cbind(response, n - response) ~ 
                 Treatment * Distance + Pollinators + usage +
                 (1|fRound),
               family = binomial,
               data = sim.glmm(model4))
  
  drop1(fit, test= "Chi")[4,4]
}

sim.test.pval.iv()







# (b) - set them all off running

# nsims <- 1000    ### use this for real - blanked out to save running time!
nsims <- 1

sim.pvals.i <- pbsapply(1:nsims, sim.test.pval.i)

sim.pvals.ii <- pbsapply(1:nsims, sim.test.pval.ii)

sim.pvals.iii <- pbsapply(1:nsims, sim.test.pval.iii)

sim.pvals.iv <- pbsapply(1:nsims, sim.test.pval.iv)





# (c) - extract power for each model and make a dataframe

model <- c("i", "ii.TP", "ii.TD", "iii", "iv")

sim.pvals.ii <- t(sim.pvals.ii)
colnames(sim.pvals.ii) <- c("Treatment:Pollinators","Treatment:Distance")


# first change any failures to converge into 1's
sim.pvals.i[which(sim.pvals.i=="Fail")] <- 1
sim.pvals.ii[which(sim.pvals.ii=="Fail")] <- 1
sim.pvals.iii[which(sim.pvals.iii=="Fail")] <- 1
sim.pvals.iv[which(sim.pvals.iv=="Fail")] <- 1

# and force them to be numerics

sim.pvals.i <- as.numeric(as.character(sim.pvals.i))
sim.pvals.ii[,1] <- as.numeric(as.character(sim.pvals.ii[,1]))
sim.pvals.ii[,2] <- as.numeric(as.character(sim.pvals.ii[,2]))
sim.pvals.iii <- as.numeric(as.character(sim.pvals.iii))
sim.pvals.iv <- as.numeric(as.character(sim.pvals.iv))

mean.i <- mean(sim.pvals.i < 0.05)
mean.ii.TP <- mean(sim.pvals.ii[,1] < 0.05)
mean.ii.TD <- mean(sim.pvals.ii[,2] < 0.05)
mean.iii <- mean(sim.pvals.iii < 0.05)
mean.iv <- mean(sim.pvals.iv < 0.05)

power <- c(mean.i,mean.ii.TP,mean.ii.TD,mean.iii,mean.iv)


# show confints - will be wide because number of sims is low
lci.i <- binom.test(table(factor(sim.pvals.i < 0.05, c(T, F))))$conf.int[1]
uci.i <- binom.test(table(factor(sim.pvals.i < 0.05, c(T, F))))$conf.int[2]

lci.ii.TP <- binom.test(table(factor(sim.pvals.ii[,1] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TP <- binom.test(table(factor(sim.pvals.ii[,1] < 0.05, c(T, F))))$conf.int[2]

lci.ii.TD <- binom.test(table(factor(sim.pvals.ii[,2] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TD <- binom.test(table(factor(sim.pvals.ii[,2] < 0.05, c(T, F))))$conf.int[2]

lci.iii <- binom.test(table(factor(sim.pvals.iii < 0.05, c(T, F))))$conf.int[1]
uci.iii <- binom.test(table(factor(sim.pvals.iii < 0.05, c(T, F))))$conf.int[2]

lci.iv <- binom.test(table(factor(sim.pvals.iv < 0.05, c(T, F))))$conf.int[1]
uci.iv <- binom.test(table(factor(sim.pvals.iv < 0.05, c(T, F))))$conf.int[2]



lci <- c(lci.i,lci.ii.TP,lci.ii.TD,lci.iii,lci.iv)
uci <- c(uci.i,uci.ii.TP,uci.ii.TD,uci.iii,uci.iv)

# put it all in a table

power.table <- data.frame(cbind(model,power,lci,uci))
power.table



### this has worked great! Now we need to repeat it for seed count and seed weight models
# the results might differ slightly due to the different error family and the smaller dataset

#### power analysis of seed count and weight ####

# import the true data to check parameters against

### Seed count & weight

dframe1s <- read.csv("Data\\SeedCount.csv")

names(dframe1s)

dframe1s$fPlot <- factor(dframe1s$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1s$fRound <- factor(dframe1s$Round)
dframe1s$fPlantNo <- factor(dframe1s$PlantNo)
dframe1s$fDistance <- factor(dframe1s$Distance) # make "Distance" available as a factor if reqd

dframe1s$Regime <- relevel(dframe1s$Regime,"Control") # relevel variables related to control level
dframe1s$Pollinators <- relevel(dframe1s$Pollinators,"Control")
dframe1s$LitUnlit <- relevel(dframe1s$LitUnlit,"Lit")

dframe1s$SeedWeight.mg <- round(dframe1s$SeedWeight*1000)

summary(dframe1s)

with(dframe1s,table(Plot,Round)) # plants per plot and round - all plants included as they should be

dframe1s$combo <- paste(dframe1s$PlantNo,dframe1s$Round,sep="-")


# extract usage from dframe1 and tie it into this one
dframe1$combo <- paste(dframe1$PlantNo,dframe1$Round,sep="-")
usages <- dframe1[,c(17:18)]
dframe1s <- merge(dframe1s,usages)


########################## FAO Paul Johnson ####

# set up that final model

model1s <- glmer.nb(SeedCount ~ Pollinators + Treatment + Distance + usage
                       + (1|fRound),
                       data = dframe1s)

summary(model1s)
drop1(model1s, test = "Chi")


# then feed this model into sim.glmm - which pulls out all the effect sizes it needs from the model to simulate the output

testdata <- sim.glmm(model1s)

### strangely, this is not actually simulating the model, but just outputting the exact input data again...?

# however, I've tracked down the problem and fixed it (hence the different function being read in at the start!)

testdata_CJM <- sim.glmm_CJM(model1s)


# plot lighting effect against variance factors
plot(response ~ Treatment, data = testdata_CJM)
plot(response ~ Pollinators, data = testdata_CJM)
plot(response ~ Distance, data = testdata_CJM)
plot(response ~ usage, data = testdata_CJM)
plot(response ~ fRound, data = testdata_CJM)



# now test for variation in pollination between lit & unlit

fit <- glmer.nb(response ~ 
                  Treatment + Pollinators + Distance + usage 
                  + (1|fRound),
                data = testdata_CJM)

summary(fit)
drop1(fit, test = "Chi")


### repeat for seed weight


# final model

model1w <- glmer.nb(SeedWeight.mg ~ Pollinators + Treatment + Distance + usage
                    + (1|fRound),
                    data = dframe1s)

summary(model1w)
drop1(model1w, test = "Chi")

# then feed this model into sim.glmm - which pulls out all the effect sizes it needs from the model to simulate the output

testdata_CJM <- sim.glmm_CJM(model1w)


# plot lighting effect against variance factors
plot(response ~ Treatment, data = testdata_CJM)
plot(response ~ Pollinators, data = testdata_CJM)
plot(response ~ Distance, data = testdata_CJM)
plot(response ~ usage, data = testdata_CJM)
plot(response ~ fRound, data = testdata_CJM)


testdata_CJM <- sim.glmm_CJM(model1w)

# now test for variation in pollination between lit & unlit

fit <- glmer.nb(response ~ 
                  Treatment + Pollinators + Distance + usage 
                + (1|fRound),
                data = testdata_CJM)

summary(fit)
drop1(fit, test = "Chi")


### (a) now set up the various functions

## seed count

# (i) basic model (as above) - no interactions

sim.test.pval.i.sc <- function(...){
  fit <- tryCatch(glmer.nb(response ~ 
                      Treatment + Pollinators + Distance + usage 
                      + (1|fRound),
                      data = sim.glmm_CJM(model1s)),
                  error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")["Treatment","Pr(Chi)"], error = function(e) print("Fail"))
}

sim.test.pval.i.sc()





# (ii) both two-way interactions simultaneously

model2s <- glmer.nb(SeedCount ~ Pollinators * Treatment + Distance * Treatment + usage
                    + (1|fRound),
                    data = dframe1s)

summary(model2s)
drop1(model2s, test = "Chi")
drop1(model2s, test = "Chi")[c(3:4),4]


sim.test.pval.ii.sc <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators * Treatment + Distance * Treatment + usage
               + (1|fRound),
               data = sim.glmm_CJM(model2s)),
               error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[c(3:4),4], error = function(e) print("Fail"))
}

sim.test.pval.ii.sc()

# (iii) Treatment:Pollinators interaction

model3s <- glmer.nb(SeedCount ~ Pollinators * Treatment + Distance + usage
                 + (1|fRound),
                data = dframe1s)

summary(model3s)
drop1(model3s, test = "Chi")
drop1(model3s, test = "Chi")[4,4]


sim.test.pval.iii.sc <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators * Treatment + Distance + usage
               + (1|fRound),
               data = sim.glmm_CJM(model3s)),
               error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[4,4], error = function(e) print("Fail"))
}

sim.test.pval.iii.sc()


# (iv) Treatment:Distance interaction

model4s <- glmer.nb(SeedCount ~ Pollinators + Treatment * Distance + usage
                    + (1|fRound),
                data = dframe1s)

summary(model4s)
drop1(model4s, test = "Chi")
drop1(model4s, test = "Chi")[4,4]


sim.test.pval.iv.sc <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators + Treatment * Distance + usage
               + (1|fRound),
               data = sim.glmm_CJM(model4s)),
               error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[4,4], error = function(e) print("Fail"))
}

sim.test.pval.iv.sc()


## seed weight

# (i) basic model (as above) - no interactions

sim.test.pval.i.sw <- function(...){
  fit <- tryCatch(glmer.nb(response ~ 
                    Treatment + Pollinators + Distance + usage 
                  + (1|fRound),
                  data = sim.glmm_CJM(model1w)),
                  error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")["Treatment","Pr(Chi)"], error = function(e) print("Fail"))
}

sim.test.pval.i.sw()

# (ii) both two-way interactions simultaneously

model2w <- glmer.nb(SeedWeight.mg ~ Pollinators * Treatment + Distance * Treatment + usage
                    + (1|fRound),
                    data = dframe1s)

summary(model2w)
drop1(model2w, test = "Chi")
drop1(model2w, test = "Chi")[c(3:4),4]


sim.test.pval.ii.sw <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators * Treatment + Distance * Treatment + usage
                  + (1|fRound),
                  data = sim.glmm_CJM(model2w)),
                  error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[c(3:4),4], error = function(e) print("Fail"))
}

sim.test.pval.ii.sw()

# (iii) Treatment:Pollinators interaction

model3w <- glmer.nb(SeedWeight.mg ~ Pollinators * Treatment + Distance + usage
                    + (1|fRound),
                    data = dframe1s)

summary(model3w)
drop1(model3w, test = "Chi")
drop1(model3w, test = "Chi")[4,4]


sim.test.pval.iii.sw <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators * Treatment + Distance + usage
                  + (1|fRound),
                  data = sim.glmm_CJM(model3w)),
                  error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[4,4], error = function(e) print("Fail"))
}

sim.test.pval.iii.sw()


# (iv) Treatment:Distance interaction

model4w <- glmer.nb(SeedWeight.mg ~ Pollinators + Treatment * Distance + usage
                    + (1|fRound),
                    data = dframe1s)

summary(model4w)
drop1(model4w, test = "Chi")
drop1(model4w, test = "Chi")[4,4]


sim.test.pval.iv.sw <- function(...){
  fit <- tryCatch(glmer.nb(response ~ Pollinators + Treatment * Distance + usage
                  + (1|fRound),
                  data = sim.glmm_CJM(model4w)),
                  error = function(e) print("Fail"))
  
  tryCatch(drop1(fit, test= "Chi")[4,4], error = function(e) print("Fail"))
}

sim.test.pval.iv.sw()


### in theory this should work



## seed count

# (b) - set them all off running

nsims <- 1000
# nsims <- 1000    ### use this for real - blanked out to save running time!


suppressMessages(sim.pvals.i.sc <- pbsapply(1:nsims, sim.test.pval.i.sc))

suppressMessages(sim.pvals.ii.sc <- pbsapply(1:nsims, sim.test.pval.ii.sc))
sim.pvals.ii.sc <- do.call("rbind",sim.pvals.ii.sc)
colnames(sim.pvals.ii.sc) <- c("Treatment:Pollinators","Treatment:Distance")

suppressMessages(sim.pvals.iii.sc <- pbsapply(1:nsims, sim.test.pval.iii.sc))

suppressMessages(sim.pvals.iv.sc <- pbsapply(1:nsims, sim.test.pval.iv.sc))





# (c) - extract power for each model and make a dataframe

# first change any failures to converge into 1's
sim.pvals.i.sc[which(sim.pvals.i.sc=="Fail")] <- 1
sim.pvals.ii.sc[which(sim.pvals.ii.sc=="Fail")] <- 1
sim.pvals.iii.sc[which(sim.pvals.iii.sc=="Fail")] <- 1
sim.pvals.iv.sc[which(sim.pvals.iv.sc=="Fail")] <- 1

# and force them to be numerics

sim.pvals.i.sc <- as.numeric(as.character(sim.pvals.i.sc))
sim.pvals.ii.sc[,1] <- as.numeric(as.character(sim.pvals.ii.sc[,1]))
sim.pvals.ii.sc[,2] <- as.numeric(as.character(sim.pvals.ii.sc[,2]))
sim.pvals.iii.sc <- as.numeric(as.character(sim.pvals.iii.sc))
sim.pvals.iv.sc <- as.numeric(as.character(sim.pvals.iv.sc))

# then calculate the powers

mean.i.sc <- mean(sim.pvals.i.sc < 0.05)
mean.ii.TP.sc <- mean(sim.pvals.ii.sc[,1] < 0.05)
mean.ii.TD.sc <- mean(sim.pvals.ii.sc[,2] < 0.05)
mean.iii.sc <- mean(sim.pvals.iii.sc < 0.05)
mean.iv.sc <- mean(sim.pvals.iv.sc < 0.05)

power.sc <- c(mean.i.sc,mean.ii.TP.sc,mean.ii.TD.sc,mean.iii.sc,mean.iv.sc)


# show confints - will be wide because number of sims is low
lci.i.sc <- binom.test(table(factor(sim.pvals.i.sc < 0.05, c(T, F))))$conf.int[1]
uci.i.sc <- binom.test(table(factor(sim.pvals.i.sc < 0.05, c(T, F))))$conf.int[2]

lci.ii.TP.sc <- binom.test(table(factor(sim.pvals.ii.sc[,1] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TP.sc <- binom.test(table(factor(sim.pvals.ii.sc[,1] < 0.05, c(T, F))))$conf.int[2]

lci.ii.TD.sc <- binom.test(table(factor(sim.pvals.ii.sc[,2] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TD.sc <- binom.test(table(factor(sim.pvals.ii.sc[,2] < 0.05, c(T, F))))$conf.int[2]

lci.iii.sc <- binom.test(table(factor(sim.pvals.iii.sc < 0.05, c(T, F))))$conf.int[1]
uci.iii.sc <- binom.test(table(factor(sim.pvals.iii.sc < 0.05, c(T, F))))$conf.int[2]

lci.iv.sc <- binom.test(table(factor(sim.pvals.iv.sc < 0.05, c(T, F))))$conf.int[1]
uci.iv.sc <- binom.test(table(factor(sim.pvals.iv.sc < 0.05, c(T, F))))$conf.int[2]



lci.sc <- c(lci.i.sc,lci.ii.TP.sc,lci.ii.TD.sc,lci.iii.sc,lci.iv.sc)
uci.sc <- c(uci.i.sc,uci.ii.TP.sc,uci.ii.TD.sc,uci.iii.sc,uci.iv.sc)

# put it all in a table

power.table.sc <- data.frame(cbind(model,power.sc,lci.sc,uci.sc))
power.table.sc




## seed weight

# (b) - set them all off running

suppressMessages(sim.pvals.i.sw <- pbsapply(1:nsims, sim.test.pval.i.sw))

suppressMessages(sim.pvals.ii.sw <- pbsapply(1:nsims, sim.test.pval.ii.sw))
sim.pvals.ii.sw <- do.call("rbind",sim.pvals.ii.sw)
colnames(sim.pvals.ii.sw) <- c("Treatment:Pollinators","Treatment:Distance")

suppressMessages(sim.pvals.iii.sw <- pbsapply(1:nsims, sim.test.pval.iii.sw))

suppressMessages(sim.pvals.iv.sw <- pbsapply(1:nsims, sim.test.pval.iv.sw))





# (c) - extract power for each model and make a dataframe

# first change any failures to converge into 1's
sim.pvals.i.sw[which(sim.pvals.i.sw=="Fail")] <- 1
sim.pvals.ii.sw[which(sim.pvals.ii.sw=="Fail")] <- 1
sim.pvals.iii.sw[which(sim.pvals.iii.sw=="Fail")] <- 1
sim.pvals.iv.sw[which(sim.pvals.iv.sw=="Fail")] <- 1

# and force them to be numerics

sim.pvals.i.sw <- as.numeric(as.character(sim.pvals.i.sw))
sim.pvals.ii.sw[,1] <- as.numeric(as.character(sim.pvals.ii.sw[,1]))
sim.pvals.ii.sw[,2] <- as.numeric(as.character(sim.pvals.ii.sw[,2]))
sim.pvals.iii.sw <- as.numeric(as.character(sim.pvals.iii.sw))
sim.pvals.iv.sw <- as.numeric(as.character(sim.pvals.iv.sw))


# then calculate the powers

mean.i.sw <- mean(sim.pvals.i.sw < 0.05)
mean.ii.TP.sw <- mean(sim.pvals.ii.sw[,1] < 0.05)
mean.ii.TD.sw <- mean(sim.pvals.ii.sw[,2] < 0.05)
mean.iii.sw <- mean(sim.pvals.iii.sw < 0.05)
mean.iv.sw <- mean(sim.pvals.iv.sw < 0.05)

power.sw <- c(mean.i.sw,mean.ii.TP.sw,mean.ii.TD.sw,mean.iii.sw,mean.iv.sw)


# show confints - will be wide because number of sims is low
lci.i.sw <- binom.test(table(factor(sim.pvals.i.sw < 0.05, c(T, F))))$conf.int[1]
uci.i.sw <- binom.test(table(factor(sim.pvals.i.sw < 0.05, c(T, F))))$conf.int[2]

lci.ii.TP.sw <- binom.test(table(factor(sim.pvals.ii.sw[,1] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TP.sw <- binom.test(table(factor(sim.pvals.ii.sw[,1] < 0.05, c(T, F))))$conf.int[2]

lci.ii.TD.sw <- binom.test(table(factor(sim.pvals.ii.sw[,2] < 0.05, c(T, F))))$conf.int[1]
uci.ii.TD.sw <- binom.test(table(factor(sim.pvals.ii.sw[,2] < 0.05, c(T, F))))$conf.int[2]

lci.iii.sw <- binom.test(table(factor(sim.pvals.iii.sw < 0.05, c(T, F))))$conf.int[1]
uci.iii.sw <- binom.test(table(factor(sim.pvals.iii.sw < 0.05, c(T, F))))$conf.int[2]

lci.iv.sw <- binom.test(table(factor(sim.pvals.iv.sw < 0.05, c(T, F))))$conf.int[1]
uci.iv.sw <- binom.test(table(factor(sim.pvals.iv.sw < 0.05, c(T, F))))$conf.int[2]



lci.sw <- c(lci.i.sw,lci.ii.TP.sw,lci.ii.TD.sw,lci.iii.sw,lci.iv.sw)
uci.sw <- c(uci.i.sw,uci.ii.TP.sw,uci.ii.TD.sw,uci.iii.sw,uci.iv.sw)

# put it all in a table

power.table.sw <- data.frame(cbind(model,power.sw,lci.sw,uci.sw))
power.table.sw




