#############################################
#  Binomial logistic regression of seed set #
#############################################

rm(list=ls())

# make a list of required packages, check they're all installed, then load them up

j <- c("rstudioapi","plyr","lme4","MASS","car","multcomp","AICcmodavg","effects","ggplot2","scales")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://star-www.st-andrews.ac.uk/cran/")

lapply(j, require, character.only = TRUE)

# source required functions
# a function for plotting glmer residuals (deviance and Pearson)
# a function for panel plots in ggplot2 - see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# a function to check whether model convergence is sufficient - see https://stats.stackexchange.com/questions/110004/how-scared-should-we-be-about-convergence-warnings-in-lme4

f <- c("CheckResidsFunction.R","MultiplotFunction.R","CheckConvergenceFunction.R")
lapply(f, source)


# set working directory to current script location and load up data (from a folder named "Data" in the same location)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dframe1 <- read.csv("Data\\SeedSetBinom.csv")
# dframe1 <- read.csv(file.choose())

names(dframe1)


dframe1$fPlot <- factor(dframe1$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1$fRound <- factor(dframe1$Round)
dframe1$fPlantNo <- factor(dframe1$PlantNo)
dframe1$fDistance <- factor(dframe1$Distance) # make "Distance" available as a factor if reqd

dframe1$Regime <- relevel(dframe1$Regime,"Control") # relevel variables related to control level
dframe1$Pollinators <- relevel(dframe1$Pollinators,"Control")
dframe1$LitUnlit <- relevel(dframe1$LitUnlit, "Lit")

summary(dframe1)

with(dframe1,table(Plot,Round)) # plants per plot and round - all plants included as they should be

# finally we want to label all plants by whether it's the first, second, third or fourth time they're used

# take *only* plants with more than one use
# use dframe1 duplicate to determine which these are
dframeduplicate <- dframe1
dframeduplicate$count <- 1
dframecollapse <- ddply(dframeduplicate, .(fPlantNo), numcolwise(sum))

dframereused <- dframecollapse[which(dframecollapse$count>1),]

reused <- droplevels(dframereused$fPlantNo)

# now we want a variable for usage of each plant - is there a change from 1st to 2nd to 3rd?

# first sort dataframe by fPlantNo
dframe1 <- dframe1[order(dframe1$fPlantNo),]
# now run down counting down sequences of fPlantNo
dframe1$usage <- sequence(rle(as.character(dframe1$fPlantNo))$lengths)




# also need alternatively coded dataframe

dframeYN <- read.csv("Data\\SeedSetYN.csv")
# dframeYN <- read.csv(file.choose())

names(dframeYN)

dframeYN$fPlot <- factor(dframeYN$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframeYN$fRound <- factor(dframeYN$Round)
dframeYN$fPlantNo <- factor(dframeYN$PlantNo)
dframeYN$fDistance <- factor(dframeYN$Distance) # make "Distance" available as a factor if reqd

dframeYN$Regime <- relevel(dframeYN$Regime,"Control") # relevel variables related to control level
dframeYN$Pollinators <- relevel(dframeYN$Pollinators,"Control")
dframeYN$LitUnlit <- relevel(dframeYN$LitUnlit, "Unlit")

summary(dframeYN)

with(dframeYN,table(Plot,Round)) # seedheads per plot and round - fairly even spread, no obvious patterns, good


# extract usage from dframe1 and tie it into this one
dframe1$combo <- paste(dframe1$PlantNo,dframe1$Round,sep="-")
usages <- dframe1[,c(17:18)]

dframeYN$combo <- paste(dframeYN$PlantNo,dframeYN$Round,sep="-")
dframeYN <- merge(dframeYN,usages)




###### Question 1 - complementarity and redundancy of pollination

### Chance of pollination

summary(dframe1)

model1 <- glmer(cbind(Successes,Failures) ~ Pollinators + LitUnlit + Distance + usage
                + (1|fPlantNo),
                family = binomial,
                data = dframe1)


summary(model1)
drop1(model1, test = "Chi")

# set up a post-hoc comparison between levels
summary(glht(model1, mcp(Pollinators="Tukey")))



# look only within unlit treatment

dframe1c <- subset(dframe1,Light=="CON")

# drop out the 'control' pollinator treatment too as it has only failures & so makes the model too complex
dframe1cres <- subset(dframe1c, Pollinators!="Control")

dframe1c$Pollinators <- droplevels(dframe1c$Pollinators)

model1c <- glmer(cbind(Successes,Failures) ~ Pollinators + Distance + usage
                +(1|fPlantNo) + (1|fRound),
                family = binomial,
                data = dframe1cres)
summary(model1c)
drop1(model1c, test = "Chi")

summary(glht(model1c, mcp(Pollinators="Tukey")))


# now, to work out observed proportions let's bring in the YN coded dataset...

summary(dframeYN)

# we can use this to get information about the observed pollination rate under each treatment
# that's the mean value for SeedSetYN
summary(dframeYN[which(dframeYN$Pollinators=="All"),])
summary(dframeYN[which(dframeYN$Pollinators=="Diurnal"),])
summary(dframeYN[which(dframeYN$Pollinators=="Nocturnal"),])
summary(dframeYN[which(dframeYN$Pollinators=="Control"),])


# repeat for CON only

dframeYNc <- subset(dframeYN,Light=="CON")
summary(dframeYNc)

# we can use this to get information about the observed pollination rate under each treatment
summary(dframeYNc[which(dframeYNc$Pollinators=="All"),])
summary(dframeYNc[which(dframeYNc$Pollinators=="Diurnal"),])
summary(dframeYNc[which(dframeYNc$Pollinators=="Nocturnal"),])
summary(dframeYNc[which(dframeYNc$Pollinators=="Control"),])


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

summary(dframe1s)

with(dframe1s,table(Plot,Round)) # plants per plot and round - all plants included as they should be

dframe1s$combo <- paste(dframe1s$PlantNo,dframe1s$Round,sep="-")
dframe1s <- merge(dframe1s,usages)


# test

hist(dframe1s$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1s$SeedCount ~ dframe1s$Pollinators)
boxplot(dframe1s$SeedCount ~ dframe1s$Distance)

model1s <- glmer(SeedCount ~ Pollinators + LitUnlit + Distance + usage
                 +(1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1s)
summary(model1s)
drop1(model1s, test = "Chi")

# possible failure to converge - check this out further
chkconv(model1s)


summary(glht(model1s, mcp(Pollinators="Tukey")))


# get observed mean & se for each grouping
dframe1sa <- dframe1s[which(dframe1s$Pollinators=="All"),]
mean(dframe1sa$SeedCount)
sd(dframe1sa$SeedCount)/nrow(dframe1sa)
# get no. flowers and no. plants
nrow(dframe1sa)
droplevels(dframe1sa$fPlantNo)

dframe1sd <- dframe1s[which(dframe1s$Pollinators=="Diurnal"),]
mean(dframe1sd$SeedCount)
sd(dframe1sd$SeedCount)/nrow(dframe1sd)
# get no. flowers and no. plants
nrow(dframe1sd)
droplevels(dframe1sd$fPlantNo)

dframe1sn <- dframe1s[which(dframe1s$Pollinators=="Nocturnal"),]
mean(dframe1sn$SeedCount)
sd(dframe1sn$SeedCount)/nrow(dframe1sn)
# get no. flowers and no. plants
nrow(dframe1sn)
droplevels(dframe1sn$fPlantNo)

dframe1scon <- dframe1s[which(dframe1s$Pollinators=="Control"),]
mean(dframe1scon$SeedCount)
sd(dframe1scon$SeedCount)/nrow(dframe1scon)
# get no. flowers and no. plants
nrow(dframe1scon)
droplevels(dframe1scon$fPlantNo)

#Restrict to unlit
dframe1sc <- subset(dframe1s,Light=="CON")
summary(dframe1sc)

hist(dframe1sc$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1sc$SeedCount ~ dframe1sc$Pollinators)

model1sc <- glmer(SeedCount ~ Pollinators + Distance + usage
                 +(1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1sc)
summary(model1sc)
drop1(model1sc, test = "Chi")

summary(glht(model1sc, mcp(Pollinators="Tukey")))


# get observed mean & se for each grouping
dframe1sca <- dframe1sc[which(dframe1sc$Pollinators=="All"),]
mean(dframe1sca$SeedCount)
sd(dframe1sca$SeedCount)/nrow(dframe1sca)
# get no flowers & plants
nrow(dframe1sca)
droplevels(dframe1sca$fPlantNo)

dframe1scd <- dframe1sc[which(dframe1sc$Pollinators=="Diurnal"),]
mean(dframe1scd$SeedCount)
sd(dframe1scd$SeedCount)/nrow(dframe1scd)
# get no flowers & plants
nrow(dframe1scd)
droplevels(dframe1scd$fPlantNo)

dframe1scn <- dframe1sc[which(dframe1sc$Pollinators=="Nocturnal"),]
mean(dframe1scn$SeedCount)
sd(dframe1scn$SeedCount)/nrow(dframe1scn)
# get no flowers & plants
nrow(dframe1scn)
droplevels(dframe1scn$fPlantNo)

dframe1scc <- dframe1sc[which(dframe1sc$Pollinators=="Control"),]
mean(dframe1scc$SeedCount)
sd(dframe1scc$SeedCount)/nrow(dframe1scc)
# get no flowers & plants
nrow(dframe1scc)
droplevels(dframe1scc$fPlantNo)

### Seed weight

hist(dframe1s$SeedWeight)
hist(log(dframe1s$SeedWeight),10)

dframe1s$lSeedWeight <- log(dframe1s$SeedWeight,10)
plot(dframe1s$lSeedWeight ~ dframe1s$Pollinators)

model1w <- lmer(lSeedWeight ~ Pollinators + LitUnlit + Distance + usage
                +(1|fPlantNo) + (1|fRound),
                data = dframe1s)

summary(model1w)
drop1(model1w, test = "Chi")

summary(glht(model1w, mcp(Pollinators="Tukey")))


# restrict to unlit

hist(dframe1sc$SeedWeight)
hist(log(dframe1sc$SeedWeight),10)

dframe1sc$lSeedWeight <- log(dframe1sc$SeedWeight,10)
plot(dframe1sc$lSeedWeight ~ dframe1sc$Pollinators)

model1wc <- lmer(lSeedWeight ~ Pollinators + Distance + usage
                 +(1|fPlantNo) + (1|fRound),
                 data = dframe1sc)

summary(model1wc)
drop1(model1wc, test = "Chi")

summary(glht(model1w, mcp(Pollinators="Tukey")))

# for completeness, let's also work out the mean + se weights

# get observed mean & se for each grouping
dframe1sa <- dframe1s[which(dframe1s$Pollinators=="All"),]
mean(dframe1sa$SeedWeight)
sd(dframe1sa$SeedWeight)/nrow(dframe1sa)

dframe1sd <- dframe1s[which(dframe1s$Pollinators=="Diurnal"),]
mean(dframe1sd$SeedWeight)
sd(dframe1sd$SeedWeight)/nrow(dframe1sd)

dframe1sn <- dframe1s[which(dframe1s$Pollinators=="Nocturnal"),]
mean(dframe1sn$SeedWeight)
sd(dframe1sn$SeedWeight)/nrow(dframe1sn)

dframe1scon <- dframe1s[which(dframe1s$Pollinators=="Control"),]
mean(dframe1scon$SeedWeight)
sd(dframe1scon$SeedWeight)/nrow(dframe1scon)

# get observed mean & se for each grouping - unlit
dframe1sca <- dframe1sc[which(dframe1sc$Pollinators=="All"),]
mean(dframe1sca$SeedWeight)
sd(dframe1sca$SeedWeight)/nrow(dframe1sca)

dframe1scd <- dframe1sc[which(dframe1sc$Pollinators=="Diurnal"),]
mean(dframe1scd$SeedWeight)
sd(dframe1scd$SeedWeight)/nrow(dframe1scd)

dframe1scn <- dframe1sc[which(dframe1sc$Pollinators=="Nocturnal"),]
mean(dframe1scn$SeedWeight)
sd(dframe1scn$SeedWeight)/nrow(dframe1scn)

dframe1scc <- dframe1sc[which(dframe1sc$Pollinators=="Control"),]
mean(dframe1scc$SeedWeight)
sd(dframe1scc$SeedWeight)/nrow(dframe1scc)



### figures

#PollinatorsYN

modelYNf <- glmer(SeedSetYN ~ Pollinators + LitUnlit + usage + (1|fPlot),
                 family=binomial,
                  data = dframeYN)

summary(modelYNf)
drop1(modelYNf, test = "Chi")



newdata1<-expand.grid(Pollinators=(c("Control","All","Diurnal","Nocturnal")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errors <- data.frame(effect(c("Pollinators"),modelYNf))

newdata1 <- merge(newdata1,errors)

newdata1$Pollinators <- revalue(newdata1$Pollinators, c("Control"="Caged","All"="Open"))
newdata1  

newdata1$Pollinators<-relevel(newdata1$Pollinators,ref="Open")

#Plot


g1 <- ggplot(newdata1,
              aes(x=Pollinators, y=fit, fill=Pollinators))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray30","white","gray70","gray50"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x="Pollinator treatment", y="Pollination success")

g1





###### Question 2 - effect of lighting (fully lit vs unlit)

dframe1t <- dframe1[which(dframe1$Regime!="Midnight"),]

model2 <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance + usage
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1t)

summary(model2)
drop1(model2, test = "Chi")

chkres(model2)

# break down into pollinator classes
# don't bother with unlit control as too little pollination occurred

# nocturnal only

dframe1tn <- subset(dframe1t,Pollinators=="Nocturnal")

model2n <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1tn)

summary(model2n)
drop1(model2n, test = "Chi")

chkres(model2n)



# diurnal only

dframe1td <- subset(dframe1t,Pollinators=="Diurnal")

model2d <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1td)

summary(model2d)
drop1(model2d, test = "Chi")


# open only

dframe1to <- subset(dframe1t,Pollinators=="All")

model2o <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1to)

summary(model2o)
drop1(model2o, test = "Chi")



### Seed Count

dframe1st <- dframe1s[which(dframe1s$Regime!="Midnight"),]
summary(dframe1st)

dframe1st$LitUnlit <- relevel(dframe1st$LitUnlit, "Unlit")

hist(dframe1st$SeedCount)
plot(dframe1st$SeedCount ~ interaction(dframe1st$LitUnlit,dframe1st$Distance))

model2a <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance + usage
                + (1|fPlantNo) + (1|fRound),
                family = poisson (link = "log"),
                data = dframe1st)

summary(model2a)
drop1(model2a, test = "Chi")

chkres(model2a)

# open pollination

dframe1sto <- subset(dframe1st,Pollinators=="All")

hist(dframe1sto$SeedCount)
plot(dframe1sto$SeedCount ~ dframe1sto$LitUnlit)

model2ao <- glmer(SeedCount ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1sto)

summary(model2ao)
drop1(model2ao, test = "Chi")

chkres(model2ao)


# diurnal pollination

dframe1std <- subset(dframe1st,Pollinators=="Diurnal")

hist(dframe1std$SeedCount)
plot(dframe1std$SeedCount ~ dframe1std$LitUnlit)

model2ad <- glmer(SeedCount ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1std)

summary(model2ad)
drop1(model2ad, test = "Chi")

chkres(model2ad)


# nocturnal pollination

dframe1stn <- subset(dframe1st,Pollinators=="Nocturnal")

hist(dframe1stn$SeedCount)
plot(dframe1stn$SeedCount ~ dframe1stn$LitUnlit)

model2an <- glmer(SeedCount ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1stn)

summary(model2an)
drop1(model2an, test = "Chi")

chkres(model2an)


# caged, no pollination

dframe1stc <- subset(dframe1st,Pollinators=="Control")

hist(dframe1stc$SeedCount)
plot(dframe1stc$SeedCount ~ dframe1stc$LitUnlit)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1st$SeedWeight)

# this will be easier if I convert seed weight from grams to milligrams
# rounding to nearest whole milligram will allow Poisson to be fitted
# that's not unreasonable as it's about the level of accuracy of the balance

dframe1st$iSeedWeight <- round(dframe1st$SeedWeight*1000)

hist(dframe1st$iSeedWeight)
plot(dframe1st$iSeedWeight ~ dframe1st$LitUnlit)

# Poisson

model2w <- glmer(iSeedWeight ~ LitUnlit + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson(link="log"),
                 data = dframe1st)

summary(model2w)
drop1(model2w, test = "Chi")

chkres(model2w)

# check convergence
chkconv(model2w)

# residuals are not 100% perfect
# try Gaussian + log transform

hist(log(dframe1st$SeedWeight,10))
dframe1st$lSeedWeight <- log(dframe1st$SeedWeight,10)

hist(dframe1st$lSeedWeight)
plot(dframe1st$lSeedWeight ~ dframe1st$LitUnlit)

model2wa <- lmer(lSeedWeight ~ LitUnlit + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1st)

summary(model2wa)
drop1(model2wa, test = "Chi")

chkres(model2wa)  # much worse!

# try quasi-poisson on untransformed data

model2wb <- glmmPQL(SeedWeight ~ LitUnlit + Pollinators + Distance + usage,
                 random = list(~1|fPlantNo, ~1|fRound),
                 family = quasipoisson (link = "log"),
                 data = dframe1st)

chkres.PQL(model2wb)  # still worse

summary(model2wa)
Anova(model2wa, type = "III")

# Poisson weren't perfect but best fit of the lot nonetheless
summary(model2w)
drop1(model2w, test = "Chi")


# open pollination

dframe1sto <- subset(dframe1st,Pollinators=="All")

hist(dframe1sto$SeedWeight)
hist(dframe1sto$iSeedWeight)
plot(dframe1sto$iSeedWeight ~ dframe1sto$LitUnlit)

model2wo <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1sto)

chkres(model2wo)

summary(model2wo)
drop1(model2wo, test = "Chi")

# diurnal pollination

dframe1std <- subset(dframe1st,Pollinators=="Diurnal")

hist(dframe1std$SeedWeight)
hist(dframe1std$iSeedWeight)
plot(dframe1std$iSeedWeight ~ dframe1std$LitUnlit)

model2wd <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1std)

chkres(model2wd)

summary(model2wd)
drop1(model2wd, test = "Chi")


# nocturnal pollination

dframe1stn <- subset(dframe1st,Pollinators=="Nocturnal")

hist(dframe1stn$SeedWeight)
hist(dframe1stn$iSeedWeight)
plot(dframe1stn$iSeedWeight ~ dframe1stn$LitUnlit)

model2wn <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1stn)

chkres(model2wn)

summary(model2wn)
drop1(model2wn, test = "Chi")


# caged, no pollination

dframe1stc <- subset(dframe1st,Pollinators=="Control")

hist(dframe1stc$SeedWeight)
hist(dframe1stc$iSeedWeight)
plot(dframe1stc$iSeedWeight ~ dframe1stc$LitUnlit)

# no data points!






###### Question 2 - effect of lighting cont. (part-lit vs unlit)

dframe1p <- dframe1[which(dframe1$Regime!="AllNight"),]
summary(dframe1p)

model2p <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance + usage
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1p)

summary(model2p)
drop1(model2p, test = "Chi")

chkres(model2p)

chkconv(model2p) # convergence unacceptable but this full model not for use in final analysis anyway...

# nocturnal only

dframe1pn <- subset(dframe1p,Pollinators=="Nocturnal")

model2pn <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1pn)

summary(model2pn)
drop1(model2pn, test = "Chi")

chkres(model2pn)


# diurnal only

dframe1pd <- subset(dframe1p,Pollinators=="Diurnal")

model2pd <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1pd)

summary(model2pd)
drop1(model2pd, test = "Chi")


# open only

dframe1po <- subset(dframe1p,Pollinators=="All")

model2po <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1po)

summary(model2po)
drop1(model2po, test = "Chi")



### Seed Count

dframe1sp <- dframe1s[which(dframe1s$Regime!="AllNight"),]
summary(dframe1sp)

dframe1sp$LitUnlit <- relevel(dframe1sp$LitUnlit, "Unlit")

hist(dframe1sp$SeedCount)
plot(dframe1sp$SeedCount ~ interaction(dframe1sp$LitUnlit,dframe1sp$Distance))

model2pa <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1sp)

summary(model2pa)
drop1(model2pa, test = "Chi")

chkres(model2pa)

chkconv(model2pa)

# open pollination

dframe1spo <- subset(dframe1sp,Pollinators=="All")

hist(dframe1spo$SeedCount)
plot(dframe1spo$SeedCount ~ dframe1spo$LitUnlit)

model2pao <- glmer(SeedCount ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1spo)

summary(model2pao)
drop1(model2pao, test = "Chi")

chkres(model2pao)


# diurnal pollination

dframe1spd <- subset(dframe1sp,Pollinators=="Diurnal")

hist(dframe1spd$SeedCount)
plot(dframe1spd$SeedCount ~ dframe1spd$LitUnlit)

model2pad <- glmer(SeedCount ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1spd)

summary(model2pad)
drop1(model2pad, test = "Chi")

chkres(model2pad)


# nocturnal pollination

dframe1spn <- subset(dframe1sp,Pollinators=="Nocturnal")

hist(dframe1spn$SeedCount)
plot(dframe1spn$SeedCount ~ dframe1spn$LitUnlit)

model2pan <- glmer(SeedCount ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1spn)

summary(model2pan)
drop1(model2pan, test = "Chi")

chkres(model2pan)


# caged, no pollination

dframe1spc <- subset(dframe1sp,Pollinators=="Control")

hist(dframe1spc$SeedCount)
plot(dframe1spc$SeedCount ~ dframe1spc$LitUnlit)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1sp$SeedWeight)

# again convert seed weight from grams to milligrams, rounding to nearest whole milligram

dframe1sp$iSeedWeight <- round(dframe1sp$SeedWeight*1000)

hist(dframe1sp$iSeedWeight)
plot(dframe1sp$iSeedWeight ~ dframe1sp$LitUnlit)

# Poisson

model2pw <- glmer(iSeedWeight ~ LitUnlit + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson(link="log"),
                 data = dframe1sp)

summary(model2pw)
drop1(model2pw, test = "Chi")

chkres(model2pw)

chkconv(model2pw)

# residuals are not 100% perfect
# try Gaussian + log transform

hist(log(dframe1sp$SeedWeight,10))
dframe1sp$lSeedWeight <- log(dframe1sp$SeedWeight,10)

hist(dframe1sp$lSeedWeight)
plot(dframe1sp$lSeedWeight ~ dframe1sp$LitUnlit)

model2pwa <- lmer(lSeedWeight ~ LitUnlit + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1sp)

summary(model2pwa)
drop1(model2pwa, test = "Chi")

chkres(model2pwa)  # much worse!

# try quasi-poisson on untransformed data

model2pwb <- glmmPQL(SeedWeight ~ LitUnlit + Pollinators + Distance + usage,
                    random = list(~1|fPlantNo, ~1|fRound),
                    family = quasipoisson (link = "log"),
                    data = dframe1sp)

chkres.PQL(model2pwb)  # still worse


# Poisson weren't perfect but best of the lot nonetheless, as for full-night
summary(model2pw)
drop1(model2pw, test = "Chi")


# open pollination

dframe1spo <- subset(dframe1sp,Pollinators=="All")

hist(dframe1spo$SeedWeight)
hist(dframe1spo$iSeedWeight)
plot(dframe1spo$iSeedWeight ~ dframe1spo$LitUnlit)

model2pwo <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1spo)

chkres(model2pwo)

summary(model2pwo)
drop1(model2pwo, test = "Chi")

# diurnal pollination

dframe1spd <- subset(dframe1sp,Pollinators=="Diurnal")

hist(dframe1spd$SeedWeight)
hist(dframe1spd$iSeedWeight)
plot(dframe1spd$iSeedWeight ~ dframe1spd$LitUnlit)

model2pwd <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1spd)

chkres(model2pwd)

summary(model2pwd)
drop1(model2pwd, test = "Chi")


# nocturnal pollination

dframe1spn <- subset(dframe1sp,Pollinators=="Nocturnal")

hist(dframe1spn$SeedWeight)
hist(dframe1spn$iSeedWeight)
plot(dframe1spn$iSeedWeight ~ dframe1spn$LitUnlit)

model2pwn <- glmer(iSeedWeight ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1spn)

chkres(model2pwn)

summary(model2pwn)
drop1(model2pwn, test = "Chi")


# caged, no pollination

dframe1spc <- subset(dframe1sp,Pollinators=="Control")

hist(dframe1spc$SeedWeight)
hist(dframe1spc$iSeedWeight)
plot(dframe1spc$iSeedWeight ~ dframe1spc$LitUnlit)

# no data points!






### figures

# FN, PR
summary(model2n)
summary(model2o)
summary(model2d)

# We need equivalent models to these based on the binomial (1/0) data frame
dframeYNt <- dframeYN[which(dframeYN$Regime!="Midnight"),]
summary(dframeYNt)

# nocturnal only
dframeYNtn <- subset(dframeYNt,Pollinators=="Nocturnal")

model2nx <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframeYNtn)

summary(model2nx)
drop1(model2nx, test = "Chi")

# this has produced identical results, as it should, so let's just press ahead with the others

# diurnal only

dframeYNtd <- subset(dframeYNt,Pollinators=="Diurnal")

model2dx <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframeYNtd)

summary(model2dx)


# open only

dframeYNto <- subset(dframeYNt,Pollinators=="All")

model2ox <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframeYNto)

summary(model2ox)

# prepare the data
#open

newdata2.1o<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2ox))

newdata2.1o <- merge(newdata2.1o,se)

#nocturnal

newdata2.1n<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2nx))
newdata2.1n <- merge(newdata2.1n,se)

#diurnal

newdata2.1d<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2dx))
newdata2.1d <- merge(newdata2.1d,se)


# stitch together

newdata2.1 <- rbind(newdata2.1o,newdata2.1d,newdata2.1n)
newdata2.1

newdata2.1$Pollinators <- revalue(newdata2.1$Pollinators, c("All"="Open"))
newdata2.1  

newdata2.1$LitUnlit <- relevel(newdata2.1$LitUnlit, "Lit")

#Plot


g2.1 <- ggplot(newdata2.1,
              aes(x=Pollinators, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1))+
  xlab(" ")+ ylab("Pollination success\n ")+ 
  ggtitle("Full-night lighting")+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

g2.1




## FN, SC

# pollination treatments separately
#open
summary(model2ao)

newdata2ao<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2ao<-model.matrix(terms(model2ao),newdata2ao)
newdata2ao$SeedCount = mm2ao %*% fixef(model2ao)
pvar2ao <- diag(mm2ao %*% tcrossprod(vcov(model2ao),mm2ao))
newdata2ao$SeedCounte <- exp(newdata2ao$SeedCount)
newdata2ao <- data.frame(
  newdata2ao
  , plo = exp(newdata2ao$SeedCount-1.96*sqrt(pvar2ao))
  , phi = exp(newdata2ao$SeedCount+1.96*sqrt(pvar2ao))
)
newdata2ao 

#nocturnal
summary(model2an)

newdata2an<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2an<-model.matrix(terms(model2an),newdata2an)
newdata2an$SeedCount = mm2an %*% fixef(model2an)
pvar2an <- diag(mm2an %*% tcrossprod(vcov(model2an),mm2an))
newdata2an$SeedCounte <- exp(newdata2an$SeedCount)
newdata2an <- data.frame(
  newdata2an
  , plo = exp(newdata2an$SeedCount-1.96*sqrt(pvar2an))
  , phi = exp(newdata2an$SeedCount+1.96*sqrt(pvar2an))
)
newdata2an

#diurnal
summary(model2ad)

newdata2ad<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2ad<-model.matrix(terms(model2ad),newdata2ad)
newdata2ad$SeedCount = mm2ad %*% fixef(model2ad)
pvar2ad <- diag(mm2ad %*% tcrossprod(vcov(model2ad),mm2ad))
newdata2ad$SeedCounte <- exp(newdata2ad$SeedCount)
newdata2ad <- data.frame(
  newdata2ad
  , plo = exp(newdata2ad$SeedCount-1.96*sqrt(pvar2ad))
  , phi = exp(newdata2ad$SeedCount+1.96*sqrt(pvar2ad))
)
newdata2ad


# stitch together

newdata2.2 <- rbind(newdata2ao,newdata2ad,newdata2an)
newdata2.2

newdata2.2$Pollinators <- revalue(newdata2.2$Pollinators, c("All"="Open"))
newdata2.2  

newdata2.2$LitUnlit <- relevel(newdata2.2$LitUnlit, "Lit")


#Plot


g2.2 <- ggplot(newdata2.2,
             aes(x=Pollinators, y=SeedCounte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,500), breaks = c(0,100,200,300,400,500))+
  xlab(" ")+ ylab("Seed count per\nseed capsule")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.2


## FN, SW

# pollination treatments separately
#open
summary(model2wo)

newdata2wo<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wo<-model.matrix(terms(model2wo),newdata2wo)
newdata2wo$iSeedWeight = mm2wo %*% fixef(model2wo)
pvar2wo <- diag(mm2wo %*% tcrossprod(vcov(model2wo),mm2wo))
newdata2wo$iSeedWeighte <- exp(newdata2wo$iSeedWeight)
newdata2wo <- data.frame(
  newdata2wo
  , plo = exp(newdata2wo$iSeedWeight-1.96*sqrt(pvar2wo))
  , phi = exp(newdata2wo$iSeedWeight+1.96*sqrt(pvar2wo))
)
newdata2wo 

#nocturnal
summary(model2wn)

newdata2wn<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wn<-model.matrix(terms(model2wn),newdata2wn)
newdata2wn$iSeedWeight = mm2wn %*% fixef(model2wn)
pvar2wn <- diag(mm2wn %*% tcrossprod(vcov(model2wn),mm2wn))
newdata2wn$iSeedWeighte <- exp(newdata2wn$iSeedWeight)
newdata2wn <- data.frame(
  newdata2wn
  , plo = exp(newdata2wn$iSeedWeight-1.96*sqrt(pvar2wn))
  , phi = exp(newdata2wn$iSeedWeight+1.96*sqrt(pvar2wn))
)

newdata2wn

#diurnal
summary(model2wd)

newdata2wd<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wd<-model.matrix(terms(model2wd),newdata2wd)
newdata2wd$iSeedWeight = mm2wd %*% fixef(model2wd)
pvar2wd <- diag(mm2wd %*% tcrossprod(vcov(model2wd),mm2wd))
newdata2wd$iSeedWeighte <- exp(newdata2wd$iSeedWeight)
newdata2wd <- data.frame(
  newdata2wd
  , plo = exp(newdata2wd$iSeedWeight-1.96*sqrt(pvar2wd))
  , phi = exp(newdata2wd$iSeedWeight+1.96*sqrt(pvar2wd))
)

newdata2wd


# stitch together

newdata2.3 <- rbind(newdata2wo,newdata2wd,newdata2wn)
newdata2.3

newdata2.3$Pollinators <- revalue(newdata2.3$Pollinators, c("All"="Open"))
newdata2.3  

newdata2.3$LitUnlit <- relevel(newdata2.3$LitUnlit, "Lit")

#Plot


g2.3 <- ggplot(newdata2.3,
               aes(x=Pollinators, y=iSeedWeighte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits=c(0,100),
                     breaks=c(0,20,40,60,80,100))+
  xlab("Pollinator treatment")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.3



# PN, PR
summary(model2pn)
summary(model2po)
summary(model2pd)

# We need equivalent models to these based on the binomial (1/0) data frame
dframeYNp <- dframeYN[which(dframeYN$Regime!="AllNight"),]
summary(dframeYNp)

# nocturnal only
dframeYNpn <- subset(dframeYNp,Pollinators=="Nocturnal")

model2pnx <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNpn)

summary(model2pnx)
drop1(model2pnx, test = "Chi")

# diurnal only

dframeYNpd <- subset(dframeYNp,Pollinators=="Diurnal")

model2pdx <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNpd)

summary(model2pdx)


# open only

dframeYNpo <- subset(dframeYNp,Pollinators=="All")

model2pox <- glmer(SeedSetYN ~ LitUnlit + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNpo)

summary(model2pox)

# prepare the data
#open

newdata2.4o<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2pox))

newdata2.4o <- merge(newdata2.4o,se)


#nocturnal

newdata2.4n<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2pnx))

newdata2.4n <- merge(newdata2.4n,se)

#diurnal

newdata2.4d<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2pdx))

newdata2.4d <- merge(newdata2.4d,se)


# stitch together

newdata2.4 <- rbind(newdata2.4o,newdata2.4d,newdata2.4n)
newdata2.4

newdata2.4$Pollinators <- revalue(newdata2.4$Pollinators, c("All"="Open"))
newdata2.4  

newdata2.4$LitUnlit <- relevel(newdata2.4$LitUnlit, "Lit")

#Plot


g2.4 <- ggplot(newdata2.4,
               aes(x=Pollinators, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1), oob=squish)+
  xlab(" ")+ ylab(" \n ")+ 
  ggtitle("Part-night lighting")+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
g2.4




## PN, SC

# pollination treatments separately
#open
summary(model2pao)

newdata2pao<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pao<-model.matrix(terms(model2pao),newdata2pao)
newdata2pao$SeedCount = mm2pao %*% fixef(model2pao)
pvar2pao <- diag(mm2pao %*% tcrossprod(vcov(model2pao),mm2pao))
newdata2pao$SeedCounte <- exp(newdata2pao$SeedCount)
newdata2pao <- data.frame(
  newdata2pao
  , plo = exp(newdata2pao$SeedCount-1.96*sqrt(pvar2pao))
  , phi = exp(newdata2pao$SeedCount+1.96*sqrt(pvar2pao))
)
newdata2pao 

#nocturnal
summary(model2pan)

newdata2pan<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pan<-model.matrix(terms(model2pan),newdata2pan)
newdata2pan$SeedCount = mm2pan %*% fixef(model2pan)
pvar2pan <- diag(mm2pan %*% tcrossprod(vcov(model2pan),mm2pan))
newdata2pan$SeedCounte <- exp(newdata2pan$SeedCount)
newdata2pan <- data.frame(
  newdata2pan
  , plo = exp(newdata2pan$SeedCount-1.96*sqrt(pvar2pan))
  , phi = exp(newdata2pan$SeedCount+1.96*sqrt(pvar2pan))
)
newdata2pan

#diurnal
summary(model2pad)

newdata2pad<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,usage=1,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pad<-model.matrix(terms(model2pad),newdata2pad)
newdata2pad$SeedCount = mm2pad %*% fixef(model2pad)
pvar2pad <- diag(mm2pad %*% tcrossprod(vcov(model2pad),mm2pad))
newdata2pad$SeedCounte <- exp(newdata2pad$SeedCount)
newdata2pad <- data.frame(
  newdata2pad
  , plo = exp(newdata2pad$SeedCount-1.96*sqrt(pvar2pad))
  , phi = exp(newdata2pad$SeedCount+1.96*sqrt(pvar2pad))
)
newdata2pad


# stitch together

newdata2.5 <- rbind(newdata2pao,newdata2pad,newdata2pan)
newdata2.5

newdata2.5$Pollinators <- revalue(newdata2.5$Pollinators, c("All"="Open"))
newdata2.5  

newdata2.5$LitUnlit <- relevel(newdata2.5$LitUnlit, "Lit")

#Plot


g2.5 <- ggplot(newdata2.5,
               aes(x=Pollinators, y=SeedCounte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,500), breaks = c(0,100,200,300,400,500))+
  xlab(" ")+ ylab(" \n ")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.5


## PN, SW

# pollination treatments separately
#open
summary(model2pwo)

newdata2pwo<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwo<-model.matrix(terms(model2pwo),newdata2pwo)
newdata2pwo$iSeedWeight = mm2pwo %*% fixef(model2pwo)
pvar2pwo <- diag(mm2pwo %*% tcrossprod(vcov(model2pwo),mm2pwo))
newdata2pwo$iSeedWeighte <- exp(newdata2pwo$iSeedWeight)
newdata2pwo <- data.frame(
  newdata2pwo
  , plo = exp(newdata2pwo$iSeedWeight-1.96*sqrt(pvar2pwo))
  , phi = exp(newdata2pwo$iSeedWeight+1.96*sqrt(pvar2pwo))
)
newdata2pwo 



#nocturnal
summary(model2pwn)

newdata2pwn<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwn<-model.matrix(terms(model2pwn),newdata2pwn)
newdata2pwn$iSeedWeight = mm2pwn %*% fixef(model2pwn)
pvar2pwn <- diag(mm2pwn %*% tcrossprod(vcov(model2pwn),mm2pwn))
newdata2pwn$iSeedWeighte <- exp(newdata2pwn$iSeedWeight)
newdata2pwn <- data.frame(
  newdata2pwn
  , plo = exp(newdata2pwn$iSeedWeight-1.96*sqrt(pvar2pwn))
  , phi = exp(newdata2pwn$iSeedWeight+1.96*sqrt(pvar2pwn))
)
newdata2pwn

#diurnal
summary(model2pwd)

newdata2pwd<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,usage=1,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwd<-model.matrix(terms(model2pwd),newdata2pwd)
newdata2pwd$iSeedWeight = mm2pwd %*% fixef(model2pwd)
pvar2pwd <- diag(mm2pwd %*% tcrossprod(vcov(model2pwd),mm2pwd))
newdata2pwd$iSeedWeighte <- exp(newdata2pwd$iSeedWeight)
newdata2pwd <- data.frame(
  newdata2pwd
  , plo = exp(newdata2pwd$iSeedWeight-1.96*sqrt(pvar2pwd))
  , phi = exp(newdata2pwd$iSeedWeight+1.96*sqrt(pvar2pwd))
)
newdata2pwd


# stitch together

newdata2.6 <- rbind(newdata2pwo,newdata2pwd,newdata2pwn)
newdata2.6

newdata2.6$Pollinators <- revalue(newdata2.6$Pollinators, c("All"="Open"))
newdata2.6  

newdata2.6$LitUnlit <- relevel(newdata2.6$LitUnlit, "Lit")

#Plot


g2.6 <- ggplot(newdata2.6,
               aes(x=Pollinators, y=iSeedWeighte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  xlab("Pollinator treatment")+ ylab(" \n ")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(0,100),
                     breaks=c(0,20,40,60,80,100))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.6



multiplot(g2.1,g2.2,g2.3,g2.4,g2.5,g2.6,cols=2)













###### Question 3 - effect of lighting treatments (LED vs HPS & FN vs PN, crossed)

### chance of pollination

# overall

dframe1a <- subset(dframe1,Regime!="Control") # Midnight + AllNight
summary(dframe1a)

model3 <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance + usage
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1a)

summary(model3)
drop1(model3, test = "Chi")

chkres(model3) # fine

chkconv(model3)

# open only

dframe1ao <- subset(dframe1a,Pollinators=="All")

model3o <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ao)

summary(model3o)
drop1(model3o, test = "Chi")



# diurnal only

dframe1ad <- subset(dframe1a,Pollinators=="Diurnal")

model3d <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ad)

summary(model3d)
drop1(model3d, test = "Chi")


# nocturnal only

dframe1an <- subset(dframe1a,Pollinators=="Nocturnal")

model3n <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1an)

summary(model3n)
drop1(model3n, test = "Chi")

chkres(model3n)


# control only

dframe1ac <- subset(dframe1a,Pollinators=="Control")  

model3c <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ac)

summary(model3c)
drop1(model3c, test = "Chi")


### Seed Count

dframe1as <- subset(dframe1s,Regime!="Control") # Midnight + AllNight
summary(dframe1as)

hist(dframe1as$SeedCount)
plot(dframe1as$SeedCount ~ dframe1as$Light)
plot(dframe1as$SeedCount ~ dframe1as$Regime)

model3a <- glmer(SeedCount ~ Light + Regime + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1as)

summary(model3a)
drop1(model3a, test = "Chi")

chkres(model3a)


# open pollination

dframe1aso <- subset(dframe1as,Pollinators=="All")

hist(dframe1aso$SeedCount)
plot(dframe1aso$SeedCount ~ dframe1aso$Light)
plot(dframe1aso$SeedCount ~ dframe1aso$Regime)

model3ao <- glmer(SeedCount ~ Light + Regime + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1aso)

summary(model3ao)
drop1(model3ao, test = "Chi")

chkres(model3ao)  

chkconv(model3ao)

# diurnal pollination

dframe1asd <- subset(dframe1as,Pollinators=="Diurnal")

hist(dframe1asd$SeedCount)
plot(dframe1asd$SeedCount ~ dframe1asd$Light)
plot(dframe1asd$SeedCount ~ dframe1asd$Regime)

model3ad <- glmer(SeedCount ~ Light + Regime + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asd)

summary(model3ad)
drop1(model3ad, test = "Chi")

chkres(model3ad)

# nocturnal pollination

dframe1asn <- subset(dframe1as,Pollinators=="Nocturnal")

hist(dframe1asn$SeedCount)
plot(dframe1asn$SeedCount ~ dframe1asn$Light)
plot(dframe1asn$SeedCount ~ dframe1asn$Regime)

model3an <- glmer(SeedCount ~ Light + Regime + Distance + usage
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asn)

summary(model3an)
drop1(model3an, test = "Chi")

chkres(model3an)


# caged, no pollination

dframe1asc <- subset(dframe1as,Pollinators=="Control")

hist(dframe1asc$SeedCount)
plot(dframe1asc$SeedCount ~ dframe1asc$Light)
plot(dframe1asc$SeedCount ~ dframe1asc$Regime)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1as$SeedWeight)

dframe1as$iSeedWeight <- round(dframe1as$SeedWeight*1000)

hist(dframe1as$iSeedWeight)

plot(dframe1as$iSeedWeight ~ dframe1as$Light)
plot(dframe1as$iSeedWeight ~ dframe1as$Regime)

model3w <- glmer(iSeedWeight ~ Light + Regime + Pollinators + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family=poisson(link="log"),
                 data = dframe1as)

summary(model3w)
drop1(model3w, test = "Chi")

chkres(model3w)  
chkconv(model3w)

# open pollination

dframe1aso <- subset(dframe1as,Pollinators=="All")

hist(dframe1aso$iSeedWeight)
plot(dframe1aso$iSeedWeight ~ dframe1aso$Light)
plot(dframe1aso$iSeedWeight ~ dframe1aso$Regime)

model3wo <- glmer(iSeedWeight ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family=poisson(link="log"),
                 data = dframe1aso)

chkres(model3wo)

summary(model3wo)
drop1(model3wo, test = "Chi")

chkconv(model3wo)

# diurnal pollination

dframe1asd <- subset(dframe1as,Pollinators=="Diurnal")

hist(dframe1asd$iSeedWeight)
plot(dframe1asd$iSeedWeight ~ dframe1asd$Light)
plot(dframe1asd$iSeedWeight ~ dframe1asd$Regime)

model3wd <- glmer(iSeedWeight ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family=poisson(link="log"),
                 data = dframe1asd)

chkres(model3wd)

summary(model3wd)
drop1(model3wd, test = "Chi")


# nocturnal pollination

dframe1asn <- subset(dframe1as,Pollinators=="Nocturnal")

hist(dframe1asn$iSeedWeight)
plot(dframe1asn$iSeedWeight ~ dframe1asn$Light)
plot(dframe1asn$iSeedWeight ~ dframe1asn$Regime)

model3wn <- glmer(iSeedWeight ~ Light + Regime + Distance + usage
                 + (1|fPlantNo) + (1|fRound),
                 family=poisson(link="log"),
                 data = dframe1asn)

chkres(model3wn)

summary(model3wn)
drop1(model3wn, test = "Chi")

# caged, no pollination

dframe1asc <- subset(dframe1as,Pollinators=="Control")

hist(dframe1asc$iSeedWeight)
plot(dframe1asc$iSeedWeight ~ dframe1asc$Light)
plot(dframe1asc$iSeedWeight ~ dframe1asc$Regime)

# only 2 data points so no point going further!




### figures

# Light + Regime

# nocturnal only, pollination rate (the only analysis with a sig. effect)

# need alternative coding of data

dframeYNa <- subset(dframeYN,Regime!="Control") # Midnight + AllNight
summary(dframeYNa)

dframeYNa$Regime <- relevel(dframeYNa$Regime,"Midnight")
dframeYNa$Light <- relevel(dframeYNa$Light,"LED")

dframeYNan <- subset(dframeYNa,Pollinators=="Nocturnal")
summary(dframeYNan)

#revised model

model3YN <- glmer(SeedSetYN ~ Light + Regime + Distance + usage
                   + (1|fPlantNo) + (1|fRound),
                   family = binomial (link = "logit"),
                   data = dframeYNan)

summary(model3YN)
drop1(model3YN, test = "Chi")   # very similar outputs so fine to proceed

# prep data
newdata3<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
se <- data.frame(effect(c("Light","Regime"),model3YN))
newdata3 <- merge(newdata3,se)

newdata3

newdata3$Regime <- revalue(newdata3$Regime, c("AllNight"="Full-night","Midnight"="Part-night"))
newdata3

#Plot


g3 <- ggplot(newdata3,
             aes(x=Regime, y=fit, fill=Light))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  xlab("Lighting regime")+ ylab("Pollination success")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3



############################ re-analysis #################################

# we need to dig down into whether re-using some plants in multiple experimental runs
# had any impact on the results



## create a dataframe that omits all 're-used' plants (i.e. includes only 'first use')

dframe1fu <- ddply(dframe1, .(PlantNo), function(x) head(x,1))

## we need to do the same for dframeYN - this is not so simple as we may want to retain multiple rows per plant

# first generate a vector containing all combos of PlantNo & Round to be used
combos <- paste(dframe1fu$PlantNo,dframe1fu$Round,sep = "-")

# use same procedure to label each row in dframeYN
dframeYN$combo <- paste(dframeYN$PlantNo,dframeYN$Round,sep = "-")

# now select out only rows in dframeYN that are included in the vector
dframeYNfu <- dframeYN[which(dframeYN$combo %in% combos),]


# and do the same for dframe1s
dframe1s$combo <- paste(dframe1s$PlantNo,dframe1s$Round,sep = "-")
dframe1sfu <- dframe1s[which(dframe1s$combo %in% combos),]

## now we have the suitable dataframes to reanalyse everything, 
## checking that there is no bias caused by including second-, third- and even fourth-use plants in random treatments


## what follows is a full replication of the above analysis using this restricted dataset


###### Question 1 - complementarity and redundancy of pollination

### Chance of pollination

summary(dframe1fu)

model1fu <- glmer(cbind(Successes,Failures) ~ Pollinators + LitUnlit + Distance
                + (1|fRound),
                family = binomial,
                data = dframe1fu)


summary(model1fu)
drop1(model1fu, test = "Chi")

# set up a post-hoc comparison between levels
summary(glht(model1fu, mcp(Pollinators="Tukey")))



# look only within unlit treatment

dframe1cfu <- subset(dframe1fu,Light=="CON")
summary(dframe1cfu)

model1cfu <- glmer(cbind(Successes,Failures) ~ Pollinators + Distance
                 + (1|fRound),
                 family = binomial,
                 data = dframe1cfu)
summary(model1cfu)
drop1(model1cfu, test = "Chi")

summary(glht(model1cfu, mcp(Pollinators="Tukey")))

# check how many plants retained under each treatment
summary(dframe1fu[which(dframe1fu$Pollinators=="All"),])
summary(dframe1fu[which(dframe1fu$Pollinators=="Diurnal"),])
summary(dframe1fu[which(dframe1fu$Pollinators=="Nocturnal"),])
summary(dframe1fu[which(dframe1fu$Pollinators=="Control"),])

summary(dframe1cfu[which(dframe1cfu$Pollinators=="All"),])
summary(dframe1cfu[which(dframe1cfu$Pollinators=="Diurnal"),])
summary(dframe1cfu[which(dframe1cfu$Pollinators=="Nocturnal"),])
summary(dframe1cfu[which(dframe1cfu$Pollinators=="Control"),])



# now, to work out observed proportions let's bring in the YN coded dataset...

summary(dframeYNfu)

# we can use this to get information about the observed pollination rate under each treatment
summary(dframeYNfu[which(dframeYNfu$Pollinators=="All"),])
summary(dframeYNfu[which(dframeYNfu$Pollinators=="Diurnal"),])
summary(dframeYNfu[which(dframeYNfu$Pollinators=="Nocturnal"),])
summary(dframeYNfu[which(dframeYNfu$Pollinators=="Control"),])


# repeat for CON only

dframeYNcfu <- subset(dframeYNfu,Light=="CON")
summary(dframeYNcfu)

# we can use this to get information about the observed pollination rate under each treatment
summary(dframeYNcfu[which(dframeYNcfu$Pollinators=="All"),])
summary(dframeYNcfu[which(dframeYNcfu$Pollinators=="Diurnal"),])
summary(dframeYNcfu[which(dframeYNcfu$Pollinators=="Nocturnal"),])
summary(dframeYNcfu[which(dframeYNcfu$Pollinators=="Control"),])


### Seed count & weight

hist(dframe1sfu$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1sfu$SeedCount ~ dframe1sfu$Pollinators)
boxplot(dframe1sfu$SeedCount ~ dframe1sfu$Distance)



model1sfu <- glmer(SeedCount ~ Pollinators + LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1sfu)
summary(model1sfu)
drop1(model1sfu, test = "Chi")


summary(glht(model1sfu, mcp(Pollinators="Tukey")))



# get observed mean & se for each grouping
dframe1safu <- dframe1sfu[which(dframe1sfu$Pollinators=="All"),]
mean(dframe1safu$SeedCount)
sd(dframe1safu$SeedCount)/nrow(dframe1safu)
# get no. flowers and no. plants
nrow(dframe1safu)
droplevels(dframe1safu$fPlantNo)

dframe1sdfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Diurnal"),]
mean(dframe1sdfu$SeedCount)
sd(dframe1sdfu$SeedCount)/nrow(dframe1sdfu)
# get no. flowers and no. plants
nrow(dframe1sdfu)
droplevels(dframe1sdfu$fPlantNo)

dframe1snfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Nocturnal"),]
mean(dframe1snfu$SeedCount)
sd(dframe1snfu$SeedCount)/nrow(dframe1snfu)
# get no. flowers and no. plants
nrow(dframe1snfu)
droplevels(dframe1snfu$fPlantNo)

dframe1sconfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Control"),]
mean(dframe1sconfu$SeedCount)
sd(dframe1sconfu$SeedCount)/nrow(dframe1sconfu)
# get no. flowers and no. plants
nrow(dframe1sconfu)
droplevels(dframe1sconfu$fPlantNo)



#Restrict to unlit
dframe1scfu <- subset(dframe1sfu,Light=="CON")
summary(dframe1scfu)

hist(dframe1scfu$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1scfu$SeedCount ~ dframe1scfu$Pollinators)

model1scfu <- glmer(SeedCount ~ Pollinators + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1scfu)
summary(model1scfu)
drop1(model1scfu, test = "Chi")

summary(glht(model1scfu, mcp(Pollinators="Tukey")))


# get observed mean & se for each grouping
dframe1scafu <- dframe1scfu[which(dframe1scfu$Pollinators=="All"),]
mean(dframe1scafu$SeedCount)
sd(dframe1scafu$SeedCount)/nrow(dframe1scafu)
# get no. flowers & plants
nrow(dframe1scafu)
droplevels(dframe1scafu$fPlantNo)

dframe1scdfu <- dframe1scfu[which(dframe1scfu$Pollinators=="Diurnal"),]
mean(dframe1scdfu$SeedCount)
sd(dframe1scdfu$SeedCount)/nrow(dframe1scdfu)
# get no. flowers & plants
nrow(dframe1scdfu)
droplevels(dframe1scdfu$fPlantNo)

dframe1scnfu <- dframe1scfu[which(dframe1scfu$Pollinators=="Nocturnal"),]
mean(dframe1scnfu$SeedCount)
sd(dframe1scnfu$SeedCount)/nrow(dframe1scnfu)
# get no. flowers & plants
nrow(dframe1scnfu)
droplevels(dframe1scnfu$fPlantNo)

dframe1sccfu <- dframe1scfu[which(dframe1scfu$Pollinators=="Control"),]
mean(dframe1sccfu$SeedCount)
sd(dframe1sccfu$SeedCount)/nrow(dframe1sccfu)
# get no. flowers & plants
nrow(dframe1sccfu)
droplevels(dframe1sccfu$fPlantNo)


### Seed weight

hist(dframe1sfu$SeedWeight)
hist(log(dframe1sfu$SeedWeight),10)

dframe1sfu$lSeedWeight <- log(dframe1sfu$SeedWeight,10)
plot(dframe1sfu$lSeedWeight ~ dframe1sfu$Pollinators)

model1wfu <- lmer(lSeedWeight ~ Pollinators + LitUnlit + Distance
                +(1|fPlantNo) + (1|fRound),
                data = dframe1sfu)

summary(model1wfu)
drop1(model1wfu, test = "Chi")

summary(glht(model1wfu, mcp(Pollinators="Tukey")))


# restrict to unlit

hist(dframe1scfu$SeedWeight)
hist(log(dframe1scfu$SeedWeight),10)

dframe1scfu$lSeedWeight <- log(dframe1scfu$SeedWeight,10)
plot(dframe1scfu$lSeedWeight ~ dframe1scfu$Pollinators)

model1wcfu <- lmer(lSeedWeight ~ Pollinators + Distance
                 +(1|fPlantNo) + (1|fRound),
                 data = dframe1scfu)

summary(model1wcfu)
drop1(model1wcfu, test = "Chi")

summary(glht(model1wfu, mcp(Pollinators="Tukey")))

# for completeness, let's also work out the mean + se weights

# get observed mean & se for each grouping
dframe1safu <- dframe1sfu[which(dframe1sfu$Pollinators=="All"),]
mean(dframe1safu$SeedWeight)
sd(dframe1safu$SeedWeight)/nrow(dframe1safu)

dframe1sdfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Diurnal"),]
mean(dframe1sdfu$SeedWeight)
sd(dframe1sdfu$SeedWeight)/nrow(dframe1sdfu)

dframe1snfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Nocturnal"),]
mean(dframe1snfu$SeedWeight)
sd(dframe1snfu$SeedWeight)/nrow(dframe1snfu)

dframe1sconfu <- dframe1sfu[which(dframe1sfu$Pollinators=="Control"),]
mean(dframe1sconfu$SeedWeight)
sd(dframe1sconfu$SeedWeight)/nrow(dframe1sconfu)



### figures
#PollinatorsYN

modelYNffu <- glmer(SeedSetYN ~ Pollinators + LitUnlit + (1|fPlot),
                  family=binomial,
                  data = dframeYNfu)

summary(modelYNffu)
drop1(modelYNffu, test = "Chi")



newdata1fu<-expand.grid(Pollinators=(c("Control","All","Diurnal","Nocturnal")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errorsfu <- data.frame(effect(c("Pollinators"),modelYNffu))

newdata1fu <- merge(newdata1fu,errorsfu)

newdata1fu$Pollinators <- revalue(newdata1fu$Pollinators, c("Control"="Caged","All"="Open"))
newdata1fu 

newdata1fu$Pollinators<-relevel(newdata1fu$Pollinators,ref="Open")

#Plot


g1fu <- ggplot(newdata1fu,
             aes(x=Pollinators, y=fit, fill=Pollinators))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray30","white","gray70","gray50"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x="Pollinator treatment", y="Pollination success")

g1fu



###### Question 2 - effect of lighting (fully lit vs unlit)

dframe1tfu <- dframe1fu[which(dframe1fu$Regime!="Midnight"),]

model2fu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance
                + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1tfu)

summary(model2fu)
drop1(model2fu, test = "Chi")

chkres(model2fu)

# nocturnal only

dframe1tnfu <- subset(dframe1tfu,Pollinators=="Nocturnal")

model2nfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1tnfu)

summary(model2nfu)
drop1(model2nfu, test = "Chi")

chkres(model2nfu)

# model failed to converge with fPlot included, but fPlot has very low variance, so removed


# diurnal only

dframe1tdfu <- subset(dframe1tfu,Pollinators=="Diurnal")

model2dfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1tdfu)

summary(model2dfu)
drop1(model2dfu, test = "Chi")


# open only

dframe1tofu <- subset(dframe1tfu,Pollinators=="All")

model2ofu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1tofu)

summary(model2ofu)
drop1(model2ofu, test = "Chi")

# control only

dframe1tcfu <- subset(dframe1tfu,Pollinators=="Control")  

model2cfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1tcfu)

summary(model2cfu)
drop1(model2cfu, test = "Chi")



### Seed Count

dframe1stfu <- dframe1sfu[which(dframe1sfu$Regime!="Midnight"),]
summary(dframe1stfu)

dframe1stfu$LitUnlit <- relevel(dframe1stfu$LitUnlit, "Unlit")

hist(dframe1stfu$SeedCount)
plot(dframe1stfu$SeedCount ~ interaction(dframe1stfu$LitUnlit,dframe1stfu$Distance))

model2afu <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1stfu)

summary(model2afu)
drop1(model2afu, test = "Chi")

chkres(model2afu)

# open pollination

dframe1stofu <- subset(dframe1stfu,Pollinators=="All")

hist(dframe1stofu$SeedCount)
plot(dframe1stofu$SeedCount ~ dframe1stofu$LitUnlit)

model2aofu <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1stofu)

summary(model2aofu)
drop1(model2aofu, test = "Chi")

chkres(model2aofu)


# diurnal pollination

dframe1stdfu <- subset(dframe1stfu,Pollinators=="Diurnal")

hist(dframe1stdfu$SeedCount)
plot(dframe1stdfu$SeedCount ~ dframe1stdfu$LitUnlit)

model2adfu <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1stdfu)

summary(model2adfu)
drop1(model2adfu, test = "Chi")

chkres(model2adfu)


# nocturnal pollination

dframe1stnfu <- subset(dframe1stfu,Pollinators=="Nocturnal")

hist(dframe1stnfu$SeedCount)
plot(dframe1stnfu$SeedCount ~ dframe1stnfu$LitUnlit)

model2anfu <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1stnfu)

summary(model2anfu)
drop1(model2anfu, test = "Chi")

chkres(model2anfu)


# caged, no pollination

dframe1stcfu <- subset(dframe1stfu,Pollinators=="Control")

hist(dframe1stcfu$SeedCount)
plot(dframe1stcfu$SeedCount ~ dframe1stcfu$LitUnlit)

# no data points so no point going further!


### Seed Weight

hist(dframe1stfu$SeedWeight)

# this will be easier if I convert seed weight from grams to milligrams, rounding to nearest whole milligram

dframe1stfu$iSeedWeight <- round(dframe1stfu$SeedWeight*1000)

hist(dframe1stfu$iSeedWeight)
plot(dframe1stfu$iSeedWeight ~ dframe1stfu$LitUnlit)

# Poisson

model2wfu <- glmer(iSeedWeight ~ LitUnlit + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson(link="log"),
                 data = dframe1stfu)

summary(model2wfu)
drop1(model2wfu, test = "Chi")

chkres(model2wfu)

# residuals are not 100% perfect
# try Gaussian + log transform

hist(log(dframe1stfu$SeedWeight,10))
dframe1stfu$lSeedWeight <- log(dframe1stfu$SeedWeight,10)

hist(dframe1stfu$lSeedWeight)
plot(dframe1stfu$lSeedWeight ~ dframe1stfu$LitUnlit)

model2wafu <- lmer(lSeedWeight ~ LitUnlit + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1stfu)

summary(model2wafu)
drop1(model2wafu, test = "Chi")

chkres(model2wafu)  # much worse!

# try quasi-poisson on untransformed data

model2wbfu <- glmmPQL(SeedWeight ~ LitUnlit + Pollinators + Distance,
                    random = list(~1|fPlantNo, ~1|fRound),
                    family = quasipoisson (link = "log"),
                    data = dframe1stfu)

chkres.PQL(model2wbfu)  # still worse

summary(model2wafu)
Anova(model2wafu, type = "III")

# Poisson weren't perfect but best of the lot nonetheless
summary(model2wfu)
drop1(model2wfu, test = "Chi")


# open pollination

dframe1stofu <- subset(dframe1stfu,Pollinators=="All")

hist(dframe1stofu$SeedWeight)
hist(dframe1stofu$iSeedWeight)
plot(dframe1stofu$iSeedWeight ~ dframe1stofu$LitUnlit)

model2wofu <- glmer(iSeedWeight ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1stofu)

chkres(model2wofu)

summary(model2wofu)
drop1(model2wofu, test = "Chi")

# diurnal pollination

dframe1stdfu <- subset(dframe1stfu,Pollinators=="Diurnal")

hist(dframe1stdfu$SeedWeight)
hist(dframe1stdfu$iSeedWeight)
plot(dframe1stdfu$iSeedWeight ~ dframe1stdfu$LitUnlit)

model2wdfu <- glmer(iSeedWeight ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1stdfu)

chkres(model2wdfu)

summary(model2wdfu)
drop1(model2wdfu, test = "Chi")


# nocturnal pollination

dframe1stnfu <- subset(dframe1stfu,Pollinators=="Nocturnal")

hist(dframe1stnfu$SeedWeight)
hist(dframe1stnfu$iSeedWeight)
plot(dframe1stnfu$iSeedWeight ~ dframe1stnfu$LitUnlit)

model2wnfu <- glmer(iSeedWeight ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1stnfu)

chkres(model2wnfu)

summary(model2wnfu)
drop1(model2wnfu, test = "Chi")


# caged, no pollination

dframe1stcfu <- subset(dframe1stfu,Pollinators=="Control")

hist(dframe1stcfu$SeedWeight)
hist(dframe1stcfu$iSeedWeight)
plot(dframe1stcfu$iSeedWeight ~ dframe1stcfu$LitUnlit)

# no data points!






###### Question 2 - effect of lighting cont. (part-lit vs unlit)

dframe1pfu <- dframe1fu[which(dframe1fu$Regime!="AllNight"),]
summary(dframe1pfu)

model2pfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1pfu)

summary(model2pfu)
drop1(model2pfu, test = "Chi")

chkres(model2pfu)

# nocturnal only

dframe1pnfu <- subset(dframe1pfu,Pollinators=="Nocturnal")

model2pnfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframe1pnfu)

summary(model2pnfu)
drop1(model2pnfu, test = "Chi")

chkres(model2pnfu)

# model failed to converge with fPlot included, but fPlot has very low variance, so removed


# diurnal only

dframe1pdfu <- subset(dframe1pfu,Pollinators=="Diurnal")

model2pdfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframe1pdfu)

summary(model2pdfu)
drop1(model2pdfu, test = "Chi")


# open only

dframe1pofu <- subset(dframe1pfu,Pollinators=="All")

model2pofu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframe1pofu)

summary(model2pofu)
drop1(model2pofu, test = "Chi")

# control only

dframe1pcfu <- subset(dframe1pfu,Pollinators=="Control")  

model2pcfu <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframe1pcfu)

summary(model2pcfu)
drop1(model2pcfu, test = "Chi")



### Seed Count

dframe1spfu <- dframe1sfu[which(dframe1sfu$Regime!="AllNight"),]
summary(dframe1spfu)

dframe1spfu$LitUnlit <- relevel(dframe1spfu$LitUnlit, "Unlit")

hist(dframe1spfu$SeedCount)
plot(dframe1spfu$SeedCount ~ interaction(dframe1spfu$LitUnlit,dframe1spfu$Distance))

model2pafu <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1spfu)

summary(model2pafu)
drop1(model2pafu, test = "Chi")

chkres(model2pafu)

# open pollination

dframe1spofu <- subset(dframe1spfu,Pollinators=="All")

hist(dframe1spofu$SeedCount)
plot(dframe1spofu$SeedCount ~ dframe1spofu$LitUnlit)

model2paofu <- glmer(SeedCount ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson (link = "log"),
                   data = dframe1spofu)

summary(model2paofu)
drop1(model2paofu, test = "Chi")

chkres(model2paofu)


# diurnal pollination

dframe1spdfu <- subset(dframe1spfu,Pollinators=="Diurnal")

hist(dframe1spdfu$SeedCount)
plot(dframe1spdfu$SeedCount ~ dframe1spdfu$LitUnlit)

model2padfu <- glmer(SeedCount ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson (link = "log"),
                   data = dframe1spdfu)

summary(model2padfu)
drop1(model2padfu, test = "Chi")

chkres(model2padfu)


# nocturnal pollination

dframe1spnfu <- subset(dframe1spfu,Pollinators=="Nocturnal")

hist(dframe1spnfu$SeedCount)
plot(dframe1spnfu$SeedCount ~ dframe1spnfu$LitUnlit)

model2panfu <- glmer(SeedCount ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson (link = "log"),
                   data = dframe1spnfu)

summary(model2panfu)
drop1(model2panfu, test = "Chi")

chkres(model2panfu)


# caged, no pollination

dframe1spcfu <- subset(dframe1spfu,Pollinators=="Control")

hist(dframe1spcfu$SeedCount)
plot(dframe1spcfu$SeedCount ~ dframe1spcfu$LitUnlit)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1spfu$SeedWeight)

# this will be easier if I convert seed weight from grams to milligrams, rounding to nearest whole milligram

dframe1spfu$iSeedWeight <- round(dframe1spfu$SeedWeight*1000)

hist(dframe1spfu$iSeedWeight)
plot(dframe1spfu$iSeedWeight ~ dframe1spfu$LitUnlit)

# Poisson

model2pwfu <- glmer(iSeedWeight ~ LitUnlit + Pollinators + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson(link="log"),
                  data = dframe1spfu)

summary(model2pwfu)
drop1(model2pwfu, test = "Chi")

chkres(model2pwfu)



# open pollination

dframe1spofu <- subset(dframe1spfu,Pollinators=="All")

hist(dframe1spofu$SeedWeight)
hist(dframe1spofu$iSeedWeight)
plot(dframe1spofu$iSeedWeight ~ dframe1spofu$LitUnlit)

model2pwofu <- glmer(iSeedWeight ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson(link="log"),
                   data = dframe1spofu)

chkres(model2pwofu)

summary(model2pwofu)
drop1(model2pwofu, test = "Chi")

# diurnal pollination

dframe1spdfu <- subset(dframe1spfu,Pollinators=="Diurnal")

hist(dframe1spdfu$SeedWeight)
hist(dframe1spdfu$iSeedWeight)
plot(dframe1spdfu$iSeedWeight ~ dframe1spdfu$LitUnlit)

model2pwdfu <- glmer(iSeedWeight ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson(link="log"),
                   data = dframe1spdfu)

chkres(model2pwdfu)

summary(model2pwdfu)
drop1(model2pwdfu, test = "Chi")


# nocturnal pollination

dframe1spnfu <- subset(dframe1spfu,Pollinators=="Nocturnal")

hist(dframe1spnfu$SeedWeight)
hist(dframe1spnfu$iSeedWeight)
plot(dframe1spnfu$iSeedWeight ~ dframe1spnfu$LitUnlit)

model2pwnfu <- glmer(iSeedWeight ~ LitUnlit + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = poisson(link="log"),
                   data = dframe1spnfu)

chkres(model2pwnfu)

summary(model2pwnfu)
drop1(model2pwnfu, test = "Chi")


# caged, no pollination

dframe1spcfu <- subset(dframe1spfu,Pollinators=="Control")

hist(dframe1spcfu$SeedWeight)
hist(dframe1spcfu$iSeedWeight)
plot(dframe1spcfu$iSeedWeight ~ dframe1spcfu$LitUnlit)

# no data points!






### figures

# FN, PR
summary(model2nfu)
summary(model2ofu)
summary(model2dfu)

# We need equivalent models to these based on the binomial (1/0) data frame
dframeYNtfu <- dframeYNfu[which(dframeYNfu$Regime!="Midnight"),]
summary(dframeYNtfu)

# nocturnal only
dframeYNtnfu <- subset(dframeYNtfu,Pollinators=="Nocturnal")

model2nxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNtnfu)

summary(model2nxfu)
drop1(model2nxfu, test = "Chi")

# this has produced identical results, as it should, so let's just press ahead with the others

# diurnal only

dframeYNtdfu <- subset(dframeYNtfu,Pollinators=="Diurnal")

model2dxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNtdfu)

summary(model2dxfu)


# open only

dframeYNtofu <- subset(dframeYNtfu,Pollinators=="All")

model2oxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNtofu)

summary(model2oxfu)

# prepare the data
#open

newdata2.1ofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2oxfu))

newdata2.1ofu <- merge(newdata2.1ofu,se)

#nocturnal

newdata2.1nfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2nxfu))
newdata2.1nfu <- merge(newdata2.1nfu,se)

#diurnal

newdata2.1dfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2dxfu))
newdata2.1dfu <- merge(newdata2.1dfu,se)


# stitch together

newdata2.1fu <- rbind(newdata2.1ofu,newdata2.1dfu,newdata2.1nfu)
newdata2.1fu

newdata2.1fu$Pollinators <- revalue(newdata2.1fu$Pollinators, c("All"="Open"))
newdata2.1fu  

newdata2.1fu$LitUnlit <- relevel(newdata2.1fu$LitUnlit, "Lit")

#Plot


g2.1fu <- ggplot(newdata2.1fu,
               aes(x=Pollinators, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1))+
  xlab(" ")+ ylab("Pollination success\n ")+ 
  ggtitle("Full-night lighting")+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

g2.1fu




## FN, SC

# pollination treatments separately
#open
summary(model2aofu)

newdata2aofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2aofu<-model.matrix(terms(model2aofu),newdata2aofu)
newdata2aofu$SeedCount = mm2aofu %*% fixef(model2aofu)
pvar2aofu <- diag(mm2aofu %*% tcrossprod(vcov(model2aofu),mm2aofu))
newdata2aofu$SeedCounte <- exp(newdata2aofu$SeedCount)
newdata2aofu <- data.frame(
  newdata2aofu
  , plo = exp(newdata2aofu$SeedCount-1.96*sqrt(pvar2aofu))
  , phi = exp(newdata2aofu$SeedCount+1.96*sqrt(pvar2aofu))
)
newdata2aofu

#nocturnal
summary(model2anfu)

newdata2anfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2anfu<-model.matrix(terms(model2anfu),newdata2anfu)
newdata2anfu$SeedCount = mm2anfu %*% fixef(model2anfu)
pvar2anfu <- diag(mm2anfu %*% tcrossprod(vcov(model2anfu),mm2anfu))
newdata2anfu$SeedCounte <- exp(newdata2anfu$SeedCount)
newdata2anfu <- data.frame(
  newdata2anfu
  , plo = exp(newdata2anfu$SeedCount-1.96*sqrt(pvar2anfu))
  , phi = exp(newdata2anfu$SeedCount+1.96*sqrt(pvar2anfu))
)
newdata2anfu

#diurnal
summary(model2adfu)

newdata2adfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2adfu<-model.matrix(terms(model2adfu),newdata2adfu)
newdata2adfu$SeedCount = mm2adfu %*% fixef(model2adfu)
pvar2adfu <- diag(mm2adfu %*% tcrossprod(vcov(model2adfu),mm2adfu))
newdata2adfu$SeedCounte <- exp(newdata2adfu$SeedCount)
newdata2adfu <- data.frame(
  newdata2adfu
  , plo = exp(newdata2adfu$SeedCount-1.96*sqrt(pvar2adfu))
  , phi = exp(newdata2adfu$SeedCount+1.96*sqrt(pvar2adfu))
)
newdata2adfu


# stitch together

newdata2.2fu <- rbind(newdata2aofu,newdata2adfu,newdata2anfu)
newdata2.2fu

newdata2.2fu$Pollinators <- revalue(newdata2.2fu$Pollinators, c("All"="Open"))
newdata2.2fu

newdata2.2fu$LitUnlit <- relevel(newdata2.2fu$LitUnlit, "Lit")


#Plot


g2.2fu <- ggplot(newdata2.2fu,
               aes(x=Pollinators, y=SeedCounte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,300), breaks = c(0,50,100,150,200,250,300))+
  xlab(" ")+ ylab("Seed count per\nseed capsule")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.2fu


## FN, SW

# pollination treatments separately
#open
summary(model2wofu)

newdata2wofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wofu<-model.matrix(terms(model2wofu),newdata2wofu)
newdata2wofu$iSeedWeight = mm2wofu %*% fixef(model2wofu)
pvar2wofu <- diag(mm2wofu %*% tcrossprod(vcov(model2wofu),mm2wofu))
newdata2wofu$iSeedWeighte <- exp(newdata2wofu$iSeedWeight)
newdata2wofu <- data.frame(
  newdata2wofu
  , plo = exp(newdata2wofu$iSeedWeight-1.96*sqrt(pvar2wofu))
  , phi = exp(newdata2wofu$iSeedWeight+1.96*sqrt(pvar2wofu))
)
newdata2wofu 

#nocturnal
summary(model2wnfu)

newdata2wnfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wnfu<-model.matrix(terms(model2wnfu),newdata2wnfu)
newdata2wnfu$iSeedWeight = mm2wnfu %*% fixef(model2wnfu)
pvar2wnfu <- diag(mm2wnfu %*% tcrossprod(vcov(model2wnfu),mm2wnfu))
newdata2wnfu$iSeedWeighte <- exp(newdata2wnfu$iSeedWeight)
newdata2wnfu <- data.frame(
  newdata2wnfu
  , plo = exp(newdata2wnfu$iSeedWeight-1.96*sqrt(pvar2wnfu))
  , phi = exp(newdata2wnfu$iSeedWeight+1.96*sqrt(pvar2wnfu))
)

newdata2wnfu

#diurnal
summary(model2wdfu)

newdata2wdfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2wdfu<-model.matrix(terms(model2wdfu),newdata2wdfu)
newdata2wdfu$iSeedWeight = mm2wdfu %*% fixef(model2wdfu)
pvar2wdfu <- diag(mm2wdfu %*% tcrossprod(vcov(model2wdfu),mm2wdfu))
newdata2wdfu$iSeedWeighte <- exp(newdata2wdfu$iSeedWeight)
newdata2wdfu <- data.frame(
  newdata2wdfu
  , plo = exp(newdata2wdfu$iSeedWeight-1.96*sqrt(pvar2wdfu))
  , phi = exp(newdata2wdfu$iSeedWeight+1.96*sqrt(pvar2wdfu))
)

newdata2wdfu


# stitch together

newdata2.3fu <- rbind(newdata2wofu,newdata2wdfu,newdata2wnfu)
newdata2.3fu

newdata2.3fu$Pollinators <- revalue(newdata2.3fu$Pollinators, c("All"="Open"))
newdata2.3fu

newdata2.3fu$LitUnlit <- relevel(newdata2.3fu$LitUnlit, "Lit")

#Plot


g2.3fu <- ggplot(newdata2.3fu,
               aes(x=Pollinators, y=iSeedWeighte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits=c(0,120),
                     breaks=c(0,20,40,60,80,100,120))+
  xlab("Pollinator treatment")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.3fu



# PN, PR
summary(model2pnfu)
summary(model2pofu)
summary(model2pdfu)

# We need equivalent models to these based on the binomial (1/0) data frame
dframeYNpfu <- dframeYNfu[which(dframeYNfu$Regime!="AllNight"),]
summary(dframeYNpfu)

# nocturnal only
dframeYNpnfu <- subset(dframeYNpfu,Pollinators=="Nocturnal")

model2pnxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                   + (1|fRound),
                   family = binomial (link = "logit"),
                   data = dframeYNpnfu)

summary(model2pnxfu)
drop1(model2pnxfu, test = "Chi")

# diurnal only

dframeYNpdfu <- subset(dframeYNpfu,Pollinators=="Diurnal")

model2pdxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                   + (1|fRound),
                   family = binomial (link = "logit"),
                   data = dframeYNpdfu)

summary(model2pdxfu)


# open only

dframeYNpofu <- subset(dframeYNpfu,Pollinators=="All")

model2poxfu <- glmer(SeedSetYN ~ LitUnlit + Distance
                   + (1|fRound),
                   family = binomial (link = "logit"),
                   data = dframeYNpofu)

summary(model2poxfu)

# prepare the data
#open

newdata2.4ofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2poxfu))

newdata2.4ofu <- merge(newdata2.4ofu,se)


#nocturnal

newdata2.4nfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2pnxfu))

newdata2.4nfu <- merge(newdata2.4nfu,se)

#diurnal

newdata2.4dfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(effect(c("LitUnlit"),model2pdxfu))

newdata2.4dfu <- merge(newdata2.4dfu,se)


# stitch together

newdata2.4fu <- rbind(newdata2.4ofu,newdata2.4dfu,newdata2.4nfu)
newdata2.4fu

newdata2.4fu$Pollinators <- revalue(newdata2.4fu$Pollinators, c("All"="Open"))
newdata2.4fu  

newdata2.4fu$LitUnlit <- relevel(newdata2.4fu$LitUnlit, "Lit")

#Plot


g2.4fu <- ggplot(newdata2.4fu,
               aes(x=Pollinators, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1), oob=squish)+
  xlab(" ")+ ylab(" \n ")+ 
  ggtitle("Part-night lighting")+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
g2.4fu




## PN, SC

# pollination treatments separately
#open
summary(model2paofu)

newdata2paofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2paofu<-model.matrix(terms(model2paofu),newdata2paofu)
newdata2paofu$SeedCount = mm2paofu %*% fixef(model2paofu)
pvar2paofu <- diag(mm2paofu %*% tcrossprod(vcov(model2paofu),mm2paofu))
newdata2paofu$SeedCounte <- exp(newdata2paofu$SeedCount)
newdata2paofu <- data.frame(
  newdata2paofu
  , plo = exp(newdata2paofu$SeedCount-1.96*sqrt(pvar2paofu))
  , phi = exp(newdata2paofu$SeedCount+1.96*sqrt(pvar2paofu))
)
newdata2paofu

#nocturnal
summary(model2panfu)

newdata2panfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2panfu<-model.matrix(terms(model2panfu),newdata2panfu)
newdata2panfu$SeedCount = mm2panfu %*% fixef(model2panfu)
pvar2panfu <- diag(mm2panfu %*% tcrossprod(vcov(model2panfu),mm2panfu))
newdata2panfu$SeedCounte <- exp(newdata2panfu$SeedCount)
newdata2panfu <- data.frame(
  newdata2panfu
  , plo = exp(newdata2panfu$SeedCount-1.96*sqrt(pvar2panfu))
  , phi = exp(newdata2panfu$SeedCount+1.96*sqrt(pvar2panfu))
)
newdata2panfu

#diurnal
summary(model2padfu)

newdata2padfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2padfu<-model.matrix(terms(model2padfu),newdata2padfu)
newdata2padfu$SeedCount = mm2padfu %*% fixef(model2padfu)
pvar2padfu <- diag(mm2padfu %*% tcrossprod(vcov(model2padfu),mm2padfu))
newdata2padfu$SeedCounte <- exp(newdata2padfu$SeedCount)
newdata2padfu <- data.frame(
  newdata2padfu
  , plo = exp(newdata2padfu$SeedCount-1.96*sqrt(pvar2padfu))
  , phi = exp(newdata2padfu$SeedCount+1.96*sqrt(pvar2padfu))
)
newdata2padfu


# stitch together

newdata2.5fu <- rbind(newdata2paofu,newdata2padfu,newdata2panfu)
newdata2.5fu

newdata2.5fu$Pollinators <- revalue(newdata2.5fu$Pollinators, c("All"="Open"))
newdata2.5fu  

newdata2.5fu$LitUnlit <- relevel(newdata2.5fu$LitUnlit, "Lit")

#Plot


g2.5fu <- ggplot(newdata2.5fu,
               aes(x=Pollinators, y=SeedCounte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,300), breaks = c(0,50,100,150,200,250,300))+
  xlab(" ")+ ylab(" \n ")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.5fu


## PN, SW

# pollination treatments separately
#open
summary(model2pwofu)

newdata2pwofu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="All",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwofu<-model.matrix(terms(model2pwofu),newdata2pwofu)
newdata2pwofu$iSeedWeight = mm2pwofu %*% fixef(model2pwofu)
pvar2pwofu <- diag(mm2pwofu %*% tcrossprod(vcov(model2pwofu),mm2pwofu))
newdata2pwofu$iSeedWeighte <- exp(newdata2pwofu$iSeedWeight)
newdata2pwofu <- data.frame(
  newdata2pwofu
  , plo = exp(newdata2pwofu$iSeedWeight-1.96*sqrt(pvar2pwofu))
  , phi = exp(newdata2pwofu$iSeedWeight+1.96*sqrt(pvar2pwofu))
)
newdata2pwofu



#nocturnal
summary(model2pwnfu)

newdata2pwnfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Nocturnal",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwnfu<-model.matrix(terms(model2pwnfu),newdata2pwnfu)
newdata2pwnfu$iSeedWeight = mm2pwnfu %*% fixef(model2pwnfu)
pvar2pwnfu <- diag(mm2pwnfu %*% tcrossprod(vcov(model2pwnfu),mm2pwnfu))
newdata2pwnfu$iSeedWeighte <- exp(newdata2pwnfu$iSeedWeight)
newdata2pwnfu <- data.frame(
  newdata2pwnfu
  , plo = exp(newdata2pwnfu$iSeedWeight-1.96*sqrt(pvar2pwnfu))
  , phi = exp(newdata2pwnfu$iSeedWeight+1.96*sqrt(pvar2pwnfu))
)
newdata2pwnfu

#diurnal
summary(model2pwdfu)

newdata2pwdfu<-expand.grid(LitUnlit=(c("Unlit","Lit")),Pollinators="Diurnal",Distance=0,iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2pwdfu<-model.matrix(terms(model2pwdfu),newdata2pwdfu)
newdata2pwdfu$iSeedWeight = mm2pwdfu %*% fixef(model2pwdfu)
pvar2pwdfu <- diag(mm2pwdfu %*% tcrossprod(vcov(model2pwdfu),mm2pwdfu))
newdata2pwdfu$iSeedWeighte <- exp(newdata2pwdfu$iSeedWeight)
newdata2pwdfu <- data.frame(
  newdata2pwdfu
  , plo = exp(newdata2pwdfu$iSeedWeight-1.96*sqrt(pvar2pwdfu))
  , phi = exp(newdata2pwdfu$iSeedWeight+1.96*sqrt(pvar2pwdfu))
)
newdata2pwdfu


# stitch together

newdata2.6fu <- rbind(newdata2pwofu,newdata2pwdfu,newdata2pwnfu)
newdata2.6fu

newdata2.6fu$Pollinators <- revalue(newdata2.6fu$Pollinators, c("All"="Open"))
newdata2.6fu  

newdata2.6fu$LitUnlit <- relevel(newdata2.6fu$LitUnlit, "Lit")

#Plot


g2.6fu <- ggplot(newdata2.6fu,
               aes(x=Pollinators, y=iSeedWeighte, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  xlab("Pollinator treatment")+ ylab(" \n ")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(0,100),
                     breaks=c(0,20,40,60,80,100))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2.6fu



multiplot(g2.1fu,g2.2fu,g2.3fu,g2.4fu,g2.5fu,g2.6fu,cols=2)













###### Question 3 - effect of lighting treatments (LED vs HPS & FN vs PN)

### chance of pollination

# overall

dframe1afu <- subset(dframe1fu,Regime!="Control") # Midnight + AllNight
summary(dframe1afu)

model3fu <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1afu)

summary(model3fu)
drop1(model3fu, test = "Chi")

chkres(model3fu) # fine


# open only

dframe1aofu <- subset(dframe1afu,Pollinators=="All")

model3ofu <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1aofu)

summary(model3ofu)
drop1(model3ofu, test = "Chi")



# diurnal only

dframe1adfu <- subset(dframe1afu,Pollinators=="Diurnal")

model3dfu <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1adfu)

summary(model3dfu)
drop1(model3dfu, test = "Chi")


# nocturnal only

dframe1anfu <- subset(dframe1afu,Pollinators=="Nocturnal")

model3nfu <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1anfu)

summary(model3nfu)
drop1(model3nfu, test = "Chi")

chkres(model3nfu)


# control only

dframe1acfu <- subset(dframe1afu,Pollinators=="Control")  

model3cfu <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1acfu)

summary(model3cfu)
drop1(model3cfu, test = "Chi")


### Seed Count

dframe1asfu <- subset(dframe1sfu,Regime!="Control") # Midnight + AllNight
summary(dframe1asfu)

hist(dframe1asfu$SeedCount)
plot(dframe1asfu$SeedCount ~ dframe1asfu$Light)
plot(dframe1asfu$SeedCount ~ dframe1asfu$Regime)

model3afu <- glmer(SeedCount ~ Light + Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1asfu)

summary(model3afu)
drop1(model3afu, test = "Chi")

chkres(model3afu)  # possible overdispersion - try quasi-poisson

source("OverdispersalFunction.R")
overdisp_fun(model3afu)


# open pollination

dframe1asofu <- subset(dframe1asfu,Pollinators=="All")

hist(dframe1asofu$SeedCount)
plot(dframe1asofu$SeedCount ~ dframe1asofu$Light)
plot(dframe1asofu$SeedCount ~ dframe1asofu$Regime)

model3aofu <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asofu)

summary(model3aofu)
drop1(model3aofu, test = "Chi")

chkres(model3aofu)


# diurnal pollination

dframe1asdfu <- subset(dframe1asfu,Pollinators=="Diurnal")

hist(dframe1asdfu$SeedCount)
plot(dframe1asdfu$SeedCount ~ dframe1asdfu$Light)
plot(dframe1asdfu$SeedCount ~ dframe1asdfu$Regime)

model3adfu <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asdfu)

summary(model3adfu)
drop1(model3adfu, test = "Chi")

chkres(model3adfu)

# nocturnal pollination

dframe1asnfu <- subset(dframe1asfu,Pollinators=="Nocturnal")

hist(dframe1asnfu$SeedCount)
plot(dframe1asnfu$SeedCount ~ dframe1asnfu$Light)
plot(dframe1asnfu$SeedCount ~ dframe1asnfu$Regime)

model3anfu <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asnfu)

summary(model3anfu)
drop1(model3anfu, test = "Chi")

chkres(model3anfu)


# caged, no pollination

dframe1ascfu <- subset(dframe1asfu,Pollinators=="Control")

hist(dframe1ascfu$SeedCount)
plot(dframe1ascfu$SeedCount ~ dframe1ascfu$Light)
plot(dframe1ascfu$SeedCount ~ dframe1ascfu$Regime)

# no data points so no point going further!


### Seed Weight

hist(dframe1asfu$SeedWeight)

dframe1asfu$iSeedWeight <- round(dframe1asfu$SeedWeight*1000)

hist(dframe1asfu$iSeedWeight)

plot(dframe1asfu$iSeedWeight ~ dframe1asfu$Light)
plot(dframe1asfu$iSeedWeight ~ dframe1asfu$Regime)

model3wfu <- glmer(iSeedWeight ~ Light + Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family=poisson(link="log"),
                 data = dframe1asfu)

summary(model3wfu)
drop1(model3wfu, test = "Chi")

chkres(model3wfu)  
checkConv(model3wfu)

# open pollination

dframe1asofu <- subset(dframe1asfu,Pollinators=="All")

hist(dframe1asofu$iSeedWeight)
plot(dframe1asofu$iSeedWeight ~ dframe1asofu$Light)
plot(dframe1asofu$iSeedWeight ~ dframe1asofu$Regime)

model3wofu <- glmer(iSeedWeight ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family=poisson(link="log"),
                  data = dframe1asofu)

chkres(model3wofu)

summary(model3wofu)
drop1(model3wofu, test = "Chi")


# diurnal pollination

dframe1asdfu <- subset(dframe1asfu,Pollinators=="Diurnal")

hist(dframe1asdfu$iSeedWeight)
plot(dframe1asdfu$iSeedWeight ~ dframe1asdfu$Light)
plot(dframe1asdfu$iSeedWeight ~ dframe1asdfu$Regime)

model3wdfu <- glmer(iSeedWeight ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family=poisson(link="log"),
                  data = dframe1asdfu)

chkres(model3wdfu)

summary(model3wdfu)
drop1(model3wdfu, test = "Chi")


# nocturnal pollination

dframe1asnfu <- subset(dframe1asfu,Pollinators=="Nocturnal")

hist(dframe1asnfu$iSeedWeight)
plot(dframe1asnfu$iSeedWeight ~ dframe1asnfu$Light)
plot(dframe1asnfu$iSeedWeight ~ dframe1asnfu$Regime)

model3wnfu <- glmer(iSeedWeight ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family=poisson(link="log"),
                  data = dframe1asnfu)

chkres(model3wnfu)

summary(model3wnfu)
drop1(model3wnfu, test = "Chi")

# caged, no pollination

dframe1ascfu <- subset(dframe1asfu,Pollinators=="Control")

hist(dframe1ascfu$iSeedWeight)
plot(dframe1ascfu$iSeedWeight ~ dframe1ascfu$Light)
plot(dframe1ascfu$iSeedWeight ~ dframe1ascfu$Regime)

# no data points so no point going further!




### figures

# Light + Regime

# nocturnal only, pollination rate (the only analysis with a sig. effect)

# need alternative coding of data

dframeYNafu <- subset(dframeYNfu,Regime!="Control") # Midnight + AllNight
summary(dframeYNafu)

dframeYNafu$Regime <- relevel(dframeYNafu$Regime,"Midnight")
dframeYNafu$Light <- relevel(dframeYNafu$Light,"LED")

dframeYNanfu <- subset(dframeYNafu,Pollinators=="Nocturnal")
summary(dframeYNanfu)

#revised model

model3YNfu <- glmer(SeedSetYN ~ Light + Regime + Distance
                  + (1|fRound),
                  family = binomial (link = "logit"),
                  data = dframeYNanfu)

summary(model3YNfu)
drop1(model3YNfu, test = "Chi")   # very similar outputs so fine to proceed

# prep data
newdata3fu<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
se <- data.frame(effect(c("Light","Regime"),model3YNfu))
newdata3fu <- merge(newdata3fu,se)

newdata3fu

newdata3fu$Regime <- revalue(newdata3fu$Regime, c("AllNight"="Full-night","Midnight"="Part-night"))
newdata3fu

#Plot


g3fu <- ggplot(newdata3fu,
             aes(x=Regime, y=fit, fill=Light))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  xlab("Lighting regime")+ ylab("Pollination success")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3fu



############ more development ##########################

# as a very last step:
# look at relationship of usage to variables

# models

model1reused <- glmer(cbind(Successes,Failures) ~ usage
                      + (1|Round) + (1|fPlantNo),
                      family = binomial (link = "logit"),
                      data = dframe1)

summary(model1reused)
drop1(model1reused, test = "Chi")

# do this with YN for comparison

modelYNreused <- glmer(SeedSetYN ~ usage
                      + (1|Round) + (1|fPlantNo),
                      family = binomial (link = "logit"),
                      data = dframeYN)

summary(modelYNreused)
drop1(modelYNreused, test = "Chi")


# now for seed count

plot(dframe1s$SeedCount ~ dframe1s$usage)

model1sreused <- glmer(SeedCount ~ usage
                       + (1|Round) + (1|fPlantNo),
                       family = poisson (link = "log"),
                       data = dframe1s)

summary(model1sreused)
drop1(model1sreused, test="Chi")



# now for seed weight

dframe1s$iSeedWeight <- round(dframe1s$SeedWeight*1000)

plot(dframe1s$iSeedWeight ~ dframe1s$usage)

# Poisson

model1wreused <- glmer(iSeedWeight ~ usage
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson(link="log"),
                 data = dframe1s)

summary(model1wreused)
drop1(model1wreused, test = "Chi")





### figures
newdataYNreused <- expand.grid(usage=c(1:4),SeedSetYN=1)
se <- data.frame(effect(c("usage"),modelYNreused))
newdataYNreused <- merge(newdataYNreused,se)


#Plot


g4.1 <- ggplot(newdataYNreused,
               aes(x=usage, y=fit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1))+
  xlab(" ")+ ylab("Pollination success")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

g4.1




## FN, SC

# pollination treatments separately
#open
summary(model1sreused)

newdata1sreused<-expand.grid(usage=(c(1:4)),SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm1sreused<-model.matrix(terms(model1sreused),newdata1sreused)
newdata1sreused$SeedCount = mm1sreused %*% fixef(model1sreused)
pvar1sreused <- diag(mm1sreused %*% tcrossprod(vcov(model1sreused),mm1sreused))
newdata1sreused$SeedCounte <- exp(newdata1sreused$SeedCount)
newdata1sreused <- data.frame(
  newdata1sreused
  , plo = exp(newdata1sreused$SeedCount-1.96*sqrt(pvar1sreused))
  , phi = exp(newdata1sreused$SeedCount+1.96*sqrt(pvar1sreused))
)
newdata1sreused 


#Plot


g4.2 <- ggplot(newdata1sreused,
               aes(x=usage, y=SeedCounte))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  xlab("Plant usage no.")+ ylab("Seed count per\nseed capsule")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g4.2


## FN, SW

# pollination treatments separately
#open
summary(model1wreused)

newdata1wreused<-expand.grid(usage=c(1:4),iSeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm1wreused<-model.matrix(terms(model1wreused),newdata1wreused)
newdata1wreused$iSeedWeight = mm1wreused %*% fixef(model1wreused)
pvar1wreused <- diag(mm1wreused %*% tcrossprod(vcov(model1wreused),mm1wreused))
newdata1wreused$iSeedWeighte <- exp(newdata1wreused$iSeedWeight)
newdata1wreused <- data.frame(
  newdata1wreused
  , plo = exp(newdata1wreused$iSeedWeight-1.96*sqrt(pvar1wreused))
  , phi = exp(newdata1wreused$iSeedWeight+1.96*sqrt(pvar1wreused))
)
newdata1wreused 


#Plot


g4.3 <- ggplot(newdata1wreused,
               aes(x=usage, y=iSeedWeighte))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits=c(0,40),
                     breaks=c(0,10,20,30,40))+
  xlab(" ")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g4.3


multiplot(g4.1,g4.2,g4.3, cols=3)





############ further checks #########################

# we want to pull out records of seed parasitism and have a look at them
# because we felt it was invalid to investigate this, it's not coded in the dataframe...
# but it's retrievable - any flower with pollination in dframeYN but no entry in dframe 1s was parasitised

dframeYN$ExactSeed <- paste(dframeYN$combo,dframeYN$SeedNo,sep="-")
dframe1s$ExactSeed <- paste(dframe1s$fPlantNo,dframe1s$Round,dframe1s$SeedNo,sep="-")

dframepara <- merge(dframeYN, dframe1s, all = T)

# now assign it

dframepara$Parasitism <- ifelse(dframepara$SeedSetYN==1 & is.na(dframepara$SeedCount), 1, 0)

# first extract how many seedheads (of total) were parasitised
summary(dframepara$Parasitism) #11.03%

# now how many of pollinated were parasitised
dframeparaY <- dframepara[which(dframepara$SeedSetYN==1),]
summary(dframeparaY$Parasitism) #23.71%

# now look at some aspects of the parasitised flowers:
dframeparaonly <- dframepara[which(dframepara$Parasitism==1),]

summary(dframeparaonly)

# happened across *all* pollinator treatments

dframenopara <- dframepara[which(dframepara$Parasitism==0),]


# test whether parasitism rates are significantly related to pollinator treatments
# i.e. were they most likely to occur when exposed in the field at night (i.e. probably field parasitism)
# or were they random (i.e. probably glasshouse parasitism)

tbl <- table(dframeparaY$Pollinators, dframeparaY$Parasitism) 
tbl

chisq.test(tbl) 

# a simple chi-squared test suggests that parasitized flowers are drawn randomly
# from the pool of pollinated flowers

# confirm this with an lm

model1 <- lm(Parasitism ~ Pollinators,
             data = dframeparaY)
summary(model1)
drop1(model1, test = "Chi")

# yes, no difference!