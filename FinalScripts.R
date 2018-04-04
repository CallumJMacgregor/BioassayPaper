#############################################
#  Binomial logistic regression of seed set #
#############################################

rm(list=ls())

# make a list of required packages, check they're all installed, then load them up

j <- c("rstudioapi","plyr","lme4","MASS","car","multcomp","AICcmodavg","effects","ggplot2","scales","gridExtra")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://star-www.st-andrews.ac.uk/cran/")

lapply(j, require, character.only = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source required functions
# a function for plotting glmer residuals (deviance and Pearson)
# a function for panel plots in ggplot2 - see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# a function to check whether model convergence is sufficient - see https://stats.stackexchange.com/questions/110004/how-scared-should-we-be-about-convergence-warnings-in-lme4

f <- c("CheckResidsFunction.R","MultiplotFunction.R","CheckConvergenceFunction.R","OverdispersalFunction.R")
lapply(f, source)


theta <- function(model){
  summary(model)$AICtab[4] / summary(model)$AICtab[5]
}
  

# set working directory to current script location and load up data (from a folder named "Data" in the same location)

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




###### Question 1 - complementarity and redundancy of pollination; and effect of lighting (just five-way treatment)

### Chance of pollination

summary(dframe1)

model1i <- glmer(cbind(Successes,Failures) ~ Pollinators * Treatment + Treatment * Distance + usage
                + (1|fRound),
                family = binomial,
                data = dframe1)

chkconv(model1i)

summary(model1i)
drop1(model1i, test = "Chi")

# non-sig interaction terms


model1ia <- glmer(cbind(Successes,Failures) ~ Pollinators * Treatment + Distance + usage
                 + (1|fRound),
                 family = binomial,
                 data = dframe1)

chkconv(model1ia)

summary(model1ia)
drop1(model1ia, test = "Chi")

##
model1ib <- glmer(cbind(Successes,Failures) ~ Pollinators + Treatment * Distance + usage
                 + (1|fRound),
                 family = binomial,
                 data = dframe1)

chkconv(model1ib)

summary(model1ib)
drop1(model1ib, test = "Chi")


# all non-sig!

model1 <- glmer(cbind(Successes,Failures) ~ Pollinators + Treatment + Distance + usage
                 + (1|fRound),
                 family = binomial,
                 data = dframe1)

chkconv(model1)

summary(model1)
drop1(model1, test = "Chi")


# set up a post-hoc comparison between levels
# for pollinator treatments
summary(glht(model1, mcp(Pollinators="Tukey")))

# for lighting treatments
summary(glht(model1, mcp(Treatment="Tukey")))



# now, to work out observed proportions let's bring in the YN coded dataset...

summary(dframeYN)

# we can use this to get information about the observed pollination rate under each treatment
# that's the mean value for SeedSetYN
summary(dframeYN[which(dframeYN$Pollinators=="All"),])
summary(dframeYN[which(dframeYN$Pollinators=="Diurnal"),])
summary(dframeYN[which(dframeYN$Pollinators=="Nocturnal"),])
summary(dframeYN[which(dframeYN$Pollinators=="Control"),])



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

dframe1s$SeedWeight.mg <- round(dframe1s$SeedWeight*1000) # turn SeedWeight into milligrams (from grams) and make it an integer

summary(dframe1s)

with(dframe1s,table(Plot,Round)) # plants per plot and round - all plants included as they should be

dframe1s$combo <- paste(dframe1s$PlantNo,dframe1s$Round,sep="-")
dframe1s <- merge(dframe1s,usages)


# test

hist(dframe1s$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1s$SeedCount ~ dframe1s$Pollinators)
boxplot(dframe1s$SeedCount ~ dframe1s$Distance)

model1si <- glmer(SeedCount ~ Pollinators + Treatment * Distance + usage
                 + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1s)
summary(model1si)
drop1(model1si, test = "Chi")


# possible failure to converge - check this out further
chkconv(model1si)
chkres(model1si)

# possible overdispersion - check this out further
overdisp_fun(model1si)

theta(model1si)


# definitely overdispersed - theta > 15 so try neg.binom
model1si.nb <- glmer.nb(SeedCount ~ Pollinators + Treatment * Distance + usage
                       + (1|fRound),
                       data = dframe1s)

chkconv(model1si.nb)
chkres(model1si.nb)

summary(model1si.nb)
drop1(model1si.nb, test = "Chi")

# non-sig interaction term

model1s.nb <- glmer.nb(SeedCount ~ Pollinators + Treatment + Distance + usage
                        + (1|fRound),
                        data = dframe1s)

chkconv(model1s.nb)
chkres(model1s.nb)


summary(model1s.nb)
drop1(model1s.nb, test = "Chi")


summary(glht(model1s.nb, mcp(Pollinators="Tukey")))
summary(glht(model1s.nb, mcp(Treatment="Tukey")))


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




### Seed weight

hist(dframe1s$SeedWeight.mg)
hist(log(dframe1s$SeedWeight.mg),10)

dframe1s$lSeedWeight.mg <- log(dframe1s$SeedWeight.mg,10)
plot(dframe1s$lSeedWeight.mg ~ dframe1s$Pollinators)

# data is non-integer and roughly Poisson
# logged data is roughly Gaussian

model1w <- lmer(lSeedWeight.mg ~ Pollinators + Treatment * Distance + usage
                + (1|fRound),
                data = dframe1s)

summary(model1w)
drop1(model1w, test = "Chi")

chkres(model1w)

# residuals aren't great - try an NB on untransformed data


model1wi.nb <- glmer.nb(SeedWeight.mg ~ Pollinators + Treatment * Distance + usage
                       + (1|fRound),
                 data = dframe1s)

chkconv(model1wi.nb)

summary(model1wi.nb)
drop1(model1wi.nb, test = "Chi")

chkres(model1wi.nb)

# non-sig interaction

model1w.nb <- glmer.nb(SeedWeight.mg ~ Pollinators + Treatment + Distance + usage
                        + (1|fRound),
                        data = dframe1s)

chkconv(model1w.nb)
chkres(model1w.nb)

summary(model1w.nb)
drop1(model1w.nb, test = "Chi")

summary(glht(model1w.nb, mcp(Pollinators="Tukey")))
summary(glht(model1w.nb, mcp(Treatment="Tukey")))



### figures

#PollinatorsYN

modelYNf <- glmer(SeedSetYN ~ Pollinators + Treatment + Distance + usage + (1|fRound),
                 family=binomial,
                  data = dframeYN)


summary(modelYNf)
drop1(modelYNf, test = "Chi")



newdata1<-expand.grid(Pollinators=(c("Control","All","Diurnal","Nocturnal")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errors <- data.frame(Effect(c("Pollinators"),modelYNf))


newdata1 <- merge(newdata1,errors)

newdata1$Pollinators <- revalue(newdata1$Pollinators, c("Control"="Caged","All"="Open"))
newdata1  

newdata1$Pollinators<-relevel(newdata1$Pollinators,ref="Open")



# Plot


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
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x="Pollinator treatment", y="Pollination success")

g1

ggsave("NewFig2.svg", plot = g1, device = "svg", path = "Plots/", width = 15, height = 15, units = "cm")




# Lighting treatments (5-way)

summary(modelYNf)


newdata2<-expand.grid(Treatment=(c("CONControl","HPSAllNight","HPSMidnight","LEDAllNight","LEDMidnight")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errors <- data.frame(Effect(c("Treatment"),modelYNf))

newdata2 <- merge(newdata2,errors)

newdata2$Treatment <- revalue(newdata2$Treatment, c("CONControl"="Unlit",
                                                      "HPSAllNight" = "HPS - FN","HPSMidnight" = "HPS - PN",
                                                      "LEDAllNight" = "LED - FN","LEDMidnight" = "LED - PN"))
newdata2  

newdata2$Treatment<-relevel(newdata2$Treatment,ref="Unlit")

newdata2$Regime <- c("Unlit","Full-night","Part-night","Full-night","Part-night")
newdata2  


# Plot


g2 <- ggplot(newdata2,
             aes(x=Treatment, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70","gray30"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill = F)+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x=" ", y="Pollination success")

g2


## try plotting these for sc and sw as well..

# sc
summary(model1s.nb)

newdata2.s<-expand.grid(Treatment=(c("CONControl","HPSAllNight","HPSMidnight","LEDAllNight","LEDMidnight")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errors.s <- data.frame(Effect(c("Treatment"),model1s.nb))

newdata2.s <- merge(newdata2.s,errors.s)

newdata2.s$Treatment <- revalue(newdata2.s$Treatment, c("CONControl"="Unlit",
                                                    "HPSAllNight" = "HPS - FN","HPSMidnight" = "HPS - PN",
                                                    "LEDAllNight" = "LED - FN","LEDMidnight" = "LED - PN"))
newdata2.s  

newdata2.s$Treatment<-relevel(newdata2.s$Treatment,ref="Unlit")

newdata2.s$Regime <- c("Unlit","Full-night","Part-night","Full-night","Part-night")
newdata2.s  


# Plot


g2.s <- ggplot(newdata2.s,
             aes(x=Treatment, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70","gray30"))+
  scale_y_continuous(limits = c(0,150))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill = F)+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x=" ", y="Seed count per\nseed capsule")

g2.s

# sw
summary(model1w.nb)

newdata2.w<-expand.grid(Treatment=(c("CONControl","HPSAllNight","HPSMidnight","LEDAllNight","LEDMidnight")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix

errors.w <- data.frame(Effect(c("Treatment"),model1w.nb))

newdata2.w <- merge(newdata2.w,errors.w)

newdata2.w$Treatment <- revalue(newdata2.w$Treatment, c("CONControl"="Unlit",
                                                        "HPSAllNight" = "HPS - FN","HPSMidnight" = "HPS - PN",
                                                        "LEDAllNight" = "LED - FN","LEDMidnight" = "LED - PN"))
newdata2.w  

newdata2.w$Treatment<-relevel(newdata2.w$Treatment,ref="Unlit")

newdata2.w$Regime <- c("Unlit","Full-night","Part-night","Full-night","Part-night")

newdata2.w  


# Plot


g2.w <- ggplot(newdata2.w,
               aes(x=Treatment, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70","gray30"))+
  scale_y_continuous(limits = c(0,50))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x="Lighting treatment", y="Dry mass of seeds per\nseed capsule (mg)")

g2.w


## multiplot


# grab the legend from g2.w and tack it on...

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g2.w)

m <- grid.arrange(arrangeGrob(g2,g2.s,
                              g2.w + theme(legend.position="none"),
                              nrow=3, ncol=1),
                    mylegend, ncol =2, widths=c(10,2))

m


## export

ggsave("NewFig3.svg", plot = m, device = "svg", path = "Plots/", width = 20, height = 40, units = "cm")





###### Question 2 - effect of lighting (full-night vs unlit)

### pollination success

dframe1FN <- dframe1[which(!dframe1$Regime == "Midnight"), ]

summary(dframe1FN)


model2i <- glmer(cbind(Successes,Failures) ~ LitUnlit + Pollinators * Distance + usage
                 + (1|fRound),
                 family = binomial,
                 data = dframe1FN)

chkconv(model2i)

summary(model2i)
drop1(model2i, test = "Chi")

# non-sig interaction

model2 <- glmer(cbind(Successes,Failures) ~ LitUnlit + Pollinators + Distance + usage
                + (1|fRound),
                family = binomial,
                data = dframe1FN)

chkconv(model2)

summary(model2)
drop1(model2, test = "Chi")




### seed count

dframe1sFN <- dframe1s[which(!dframe1s$Regime == "Midnight"), ]

summary(dframe1sFN)


model2si <- glmer(SeedCount ~ LitUnlit *  Distance + Pollinators + usage
                 + (1|fRound),
                 family = poisson(link = "log"),
                 data = dframe1sFN)

chkconv(model2si)
theta(model2si)

summary(model2si)
drop1(model2si, test = "Chi")

# overdispersed - try nb

model2si.nb <- glmer.nb(SeedCount ~ LitUnlit * Distance + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sFN)

chkconv(model2si.nb)
chkres(model2si.nb)
summary(model2si.nb)
drop1(model2si.nb, test = "Chi")


# non-sig interaction

model2s.nb <- glmer.nb(SeedCount ~ LitUnlit + Distance + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sFN)

chkconv(model2s.nb)
chkres(model2s.nb)
summary(model2s.nb)
drop1(model2s.nb, test = "Chi")



### seed weight

model2wi.nb <- glmer.nb(SeedWeight.mg ~ LitUnlit * Distance + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sFN)


chkconv(model2wi.nb)
chkres(model2wi.nb)
summary(model2wi.nb)
drop1(model2wi.nb, test = "Chi")


# non-sig interaction

model2w.nb <- glmer.nb(SeedWeight.mg ~ LitUnlit + Distance + Pollinators + usage
                        + (1|fRound),
                        data = dframe1sFN)


chkconv(model2w.nb)
chkres(model2w.nb)
summary(model2w.nb)
drop1(model2w.nb, test = "Chi")


###### Question 3 - effect of lighting (part-night vs unlit)

### pollination success

dframe1PN <- dframe1[which(!dframe1$Regime == "AllNight"), ]

summary(dframe1PN)


model3i <- glmer(cbind(Successes,Failures) ~ LitUnlit * Distance + Pollinators + usage
                + (1|fRound),
                family = binomial,
                data = dframe1PN)

chkconv(model3i)

summary(model3i)
drop1(model3i, test = "Chi")

# non-sig interaction

model3 <- glmer(cbind(Successes,Failures) ~ LitUnlit + Distance + Pollinators + usage
                + (1|fRound),
                family = binomial,
                data = dframe1PN)

chkconv(model3)

summary(model3)
drop1(model3, test = "Chi")

### seed count

dframe1sPN <- dframe1s[which(!dframe1s$Regime == "AllNight"), ]

summary(dframe1sPN)


model3si <- glmer(SeedCount ~ LitUnlit * Distance + Pollinators + usage
                 + (1|fRound),
                 family = poisson(link = "log"),
                 data = dframe1sPN)

chkconv(model3si)
theta(model3si)

summary(model3si)
drop1(model3si, test = "Chi")

# overdispersed - try nb

model3si.nb <- glmer.nb(SeedCount ~ LitUnlit * Distance + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sPN)

chkconv(model3si.nb)
chkres(model3si.nb)
summary(model3si.nb)
drop1(model3si.nb, test = "Chi")

# non-sig interaction

model3s.nb <- glmer.nb(SeedCount ~ LitUnlit + Pollinators + Distance + usage
                       + (1|fRound),
                       data = dframe1sPN)

chkconv(model3s.nb)
chkres(model3s.nb)
summary(model3s.nb)
drop1(model3s.nb, test = "Chi")

### seed weight

model3wi.nb <- glmer.nb(SeedWeight.mg ~ LitUnlit * Distance + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sPN)


chkconv(model3wi.nb)
chkres(model3wi.nb)
summary(model3wi.nb)
drop1(model3wi.nb, test = "Chi")

# non-sig interaction

model3w.nb <- glmer.nb(SeedWeight.mg ~ LitUnlit + Pollinators + Distance + usage
                       + (1|fRound),
                       data = dframe1sPN)


chkconv(model3w.nb)
chkres(model3w.nb)
summary(model3w.nb)
drop1(model3w.nb, test = "Chi")



###### Question 4 - effect of light types (full-night vs part-night & LED vs HPS)

### pollination success

dframe1Lit <- dframe1[which(!dframe1$Regime == "Control"), ]

summary(dframe1Lit)


model4i <- glmer(cbind(Successes,Failures) ~ Regime * Distance + Light * Distance + Regime * Light + Pollinators + usage
                + (1|fRound),
                family = binomial,
                data = dframe1Lit)

chkconv(model4i)

summary(model4i)
drop1(model4i, test = "Chi")

# non-sig interaction

model4 <- glmer(cbind(Successes,Failures) ~ Regime + Light + Pollinators + Distance + usage
                + (1|fRound),
                family = binomial,
                data = dframe1Lit)

chkconv(model4)

summary(model4)
drop1(model4, test = "Chi")



### seed count

dframe1sLit <- dframe1s[which(!dframe1s$Regime == "Control"), ]

summary(dframe1sLit)


model4si <- glmer(SeedCount ~ Regime * Distance + Light * Distance + Regime * Light + Pollinators + usage
                 + (1|fRound),
                 family = poisson(link = "log"),
                 data = dframe1sLit)

chkconv(model4si)
theta(model4si)


# overdispersed - try nb

model4si.nb <- glmer.nb(SeedCount ~ Regime * Distance + Light * Distance + Regime * Light + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sLit)

chkconv(model4si.nb)
chkres(model4si.nb)
summary(model4si.nb)
drop1(model4si.nb, test = "Chi")


# non-sig interactions dropped (one retained)

model4s.nb <- glmer.nb(SeedCount ~ Regime * Light + Pollinators + Distance + usage
                       + (1|fRound),
                       data = dframe1sLit)

chkconv(model4s.nb)
chkres(model4s.nb)
summary(model4s.nb)
drop1(model4s.nb, test = "Chi")

### seed weight

model4wi.nb <- glmer.nb(SeedWeight.mg ~ Regime * Distance + Light * Distance + Regime * Light + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sLit)


chkconv(model4wi.nb)
chkres(model4wi.nb)
summary(model4wi.nb)
drop1(model4wi.nb, test = "Chi")

# non-sig interaction (one retained)


model4w.nb <- glmer.nb(SeedWeight.mg ~ Light * Distance + Regime + Pollinators + usage
                       + (1|fRound),
                       data = dframe1sLit)


chkconv(model4w.nb)
chkres(model4w.nb)
summary(model4w.nb)
drop1(model4w.nb, test = "Chi")


##### Figures

# FN, PR

summary(model2)

# we need an equivalent model based on the binomial (1/0) data frame
dframeYNFN <- dframeYN[which(dframeYN$Regime!="Midnight"),]

model2YN <- glmer(SeedSetYN ~ LitUnlit + Pollinators + Distance + usage
                + (1|fRound),
                family = binomial,
                data = dframeYNFN)

summary(model2YN)
drop1(model2YN, test = "Chi")

# prepare model predictions

newdata3.1 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedSetYN=1)
effects3.1 <- data.frame(Effect(c("LitUnlit"),model2YN))

newdata3.1 <- merge(newdata3.1,effects3.1)
newdata3.1$LitUnlit <- relevel(newdata3.1$LitUnlit, "Lit")

newdata3.1

# Plot

g3.1 <- ggplot(newdata3.1,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
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
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
  

g3.1


## FN, SC

summary(model2s.nb)

newdata3.2 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedCount=1)
effects3.2 <- data.frame(Effect(c("LitUnlit"),model2s.nb))
  
newdata3.2 <- merge(newdata3.2,effects3.2)
newdata3.2$LitUnlit <- relevel(newdata3.2$LitUnlit, "Lit")
newdata3.2

# Plot

g3.2 <- ggplot(newdata3.2,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  xlab(" ")+ ylab("Seed count per\nseed capsule")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3.2



## FN, SW

summary(model2w.nb)

newdata3.3 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedWeight=1)
effects3.3 <- data.frame(Effect(c("LitUnlit"),model2w.nb))

newdata3.3 <- merge(newdata3.3,effects3.3)
newdata3.3$LitUnlit <- relevel(newdata3.3$LitUnlit, "Lit")
newdata3.3



# Plot

g3.3 <- ggplot(newdata3.3,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,50))+
  xlab("Pollinator treatment")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+  
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3.3



# PN, PR

summary(model3)

# we need an equivalent model based on the binomial (1/0) data frame
dframeYNPN <- dframeYN[which(dframeYN$Regime!="AllNight"),]

model3YN <- glmer(SeedSetYN ~ LitUnlit + Pollinators + Distance + usage
                  + (1|fRound),
                  family = binomial,
                  data = dframeYNPN)

summary(model3YN)
drop1(model3YN, test = "Chi")

# prepare model predictions

newdata3.4 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedSetYN=1)
effects3.4 <- data.frame(Effect(c("LitUnlit"),model3YN))

newdata3.4 <- merge(newdata3.4,effects3.4)
newdata3.4$LitUnlit <- relevel(newdata3.4$LitUnlit, "Lit")

newdata3.4

# Plot

g3.4 <- ggplot(newdata3.4,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray70","gray30"))+
  scale_y_continuous(labels=percent_format(),limits = c(0,1))+
  xlab(" ")+ ylab(" ")+
  ggtitle("Part-night lighting")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))


g3.4


## PN, SC

summary(model3s.nb)

newdata3.5 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedCount=1)
effects3.5 <- data.frame(Effect(c("LitUnlit"),model3s.nb))

newdata3.5 <- merge(newdata3.5,effects3.5)
newdata3.5$LitUnlit <- relevel(newdata3.5$LitUnlit, "Lit")
newdata3.5

# Plot

g3.5 <- ggplot(newdata3.5,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray70","gray30"))+
  scale_y_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  xlab(" ")+ ylab(" ")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3.5



## PN, SW

summary(model3w.nb)

newdata3.6 <- expand.grid(LitUnlit=(c("Unlit","Lit")),SeedWeight=1)
effects3.6 <- data.frame(Effect(c("LitUnlit"),model3w.nb))

newdata3.6 <- merge(newdata3.6,effects3.6)
newdata3.6$LitUnlit <- relevel(newdata3.6$LitUnlit, "Lit")
newdata3.6


# Plot

g3.6 <- ggplot(newdata3.6,
               aes(x=LitUnlit, y=fit, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray70","gray30"))+
  scale_y_continuous(limits = c(0,50))+
  xlab("Pollinator treatment")+ ylab(" ")+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3.6


m1 <- grid.arrange(g3.1,g3.4,g3.2,g3.5,g3.3,g3.6,ncol=2)

m1


# grab the legend from g2 and tack it onto these...

# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(g2.w)

m1a <- grid.arrange(arrangeGrob(g3.1,g3.4,
                                g3.2,g3.5,
                                g3.3,g3.6 + theme(legend.position="none"),
                                nrow=3, ncol=2),
                    mylegend, ncol =2, widths=c(10,2))

m1a



ggsave("NewFig4.svg", plot = m1a, device = "svg", path = "Plots/", width = 25, height = 30, units = "cm")


### now do something similar for the crossed analysis (question 4)

## PR

summary(model4)

# need alternative coding

dframeYNLit <- dframeYN[which(dframeYN$LitUnlit!="Unlit"),]

model4YN <- glmer(SeedSetYN ~ Regime + Light + Pollinators + Distance + usage
                  + (1|fRound),
                  family = binomial,
                  data = dframeYNLit)
  
summary(model4YN)
drop1(model4YN, test = "Chi")  
  

# prepare model predictions
newdata4.1 <- expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),SeedSetYN=1)
effects4.1 <- data.frame(Effect(c("Light","Regime"),model4YN))

newdata4.1 <- merge(newdata4.1,effects4.1)

newdata4.1$Regime <- revalue(newdata4.1$Regime, c("AllNight"="Full-night","Midnight"="Part-night"))
newdata4.1

# Plot

g4.1 <- ggplot(newdata4.1,
               aes(x=Light, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  xlab("Light type")+ ylab("Pollination success\n ")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g4.1



## SC


summary(model4s.nb)
drop1(model4s.nb, test = "Chi")

# prepare model predictions
newdata4.2 <- expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),SeedCount=1)
effects4.2 <- data.frame(Effect(c("Light","Regime"),model4s.nb))

newdata4.2 <- merge(newdata4.2,effects4.2)

newdata4.2$Regime <- revalue(newdata4.2$Regime, c("AllNight"="Full-night","Midnight"="Part-night"))
newdata4.2

# Plot

g4.2 <- ggplot(newdata4.2,
               aes(x=Light, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70"))+
  scale_y_continuous(limits = c(0,150), breaks = c(0,25,50,75,100,125,150))+
  xlab("Light type")+ ylab("Seed count per\nseed capsule")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill = FALSE)+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g4.2



## SW


summary(model4w.nb)
drop1(model4w.nb, test = "Chi")

# prepare model predictions
newdata4.3 <- expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),SeedCount=1)
effects4.3 <- data.frame(Effect(c("Light","Regime"),model4w.nb))

newdata4.3 <- merge(newdata4.3,effects4.3)

newdata4.3$Regime <- revalue(newdata4.3$Regime, c("AllNight"="Full-night","Midnight"="Part-night"))
newdata4.3


# Plot

g4.3 <- ggplot(newdata4.3,
               aes(x=Light, y=fit, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray70"))+
  scale_y_continuous(limits = c(0,50))+
  xlab("Light type")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill = FALSE)+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g4.3


# arrange into a multiplot

m2 <- grid.arrange(g4.1,g4.2,g4.3,ncol=1)

m2


# grab the legend from g4.1 and tack it onto these...

mylegend1 <- g_legend(g4.1)

m2a <- grid.arrange(arrangeGrob(g4.1 + theme(legend.position="none"),
                                g4.2,g4.3,
                                nrow=3, ncol=1),
                    mylegend1, ncol =2, widths=c(10,2))

m2a



ggsave("NewFig5.svg", plot = m2a, device = "svg", path = "Plots/", width = 20, height = 30, units = "cm")



### finally plot a figure to show interaction between light type and distance for seed weight

## SW


summary(model4w.nb)
drop1(model4w.nb, test = "Chi")

# prepare model predictions
newdata5 <- expand.grid(Light=(c("HPS","LED")),Distance=(c(0,10,20)),SeedCount=1)
effects5 <- data.frame(Effect(c("Light","Distance"),model4w.nb))

newdata5 <- merge(newdata5,effects5)

newdata5


# Plot

g5 <- ggplot(newdata5,
               aes(x=Distance, y=fit, ymin = lower, ymax = upper, colour=Light, fill=Light))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(colour="black",size=0.8,width=3,
                position=position_dodge(9))+
  scale_colour_manual(values = c("darkorange2","royalblue"),
                      labels = c("HPS","LED"))+
  scale_fill_manual(values = c("yellow","lightsteelblue2"),
                    labels = c("HPS","LED"))+
  scale_y_continuous(limits = c(0,60))+
  xlab("Distance from light (m)")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())



g5



ggsave("NewFig6.svg", plot = g5, device = "svg", path = "Plots/", width = 20, height = 15, units = "cm")







####################################


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

model1sreused <- glmer.nb(SeedCount ~ usage
                       + (1|Round) + (1|fPlantNo),
                       data = dframe1s)

chkconv(model1sreused)

summary(model1sreused)
drop1(model1sreused, test="Chi")



# now for seed weight

plot(dframe1s$SeedWeight.mg ~ dframe1s$usage)

# Poisson

model1wreused <- glmer.nb(SeedWeight.mg ~ usage
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1s)

chkconv(model1wreused)

summary(model1wreused)
drop1(model1wreused, test = "Chi")





### figures
newdataYNreused <- expand.grid(usage=c(1:4),SeedSetYN=1)
se <- data.frame(Effect(c("usage"),modelYNreused))
newdataYNreused <- merge(newdataYNreused,se)


#Plot


g5.1 <- ggplot(newdataYNreused,
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
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

g5.1




## FN, SC

# pollination treatments separately
#open
summary(model1sreused)

newdata1sreused<-expand.grid(usage=(c(1:4)),SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(Effect(c("usage"),model1sreused))
newdata1sreused <- merge(newdata1sreused,se)


#Plot


g5.2 <- ggplot(newdata1sreused,
               aes(x=usage, y=fit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  xlab("Number of uses")+ ylab("Seed count per\nseed capsule")+ 
  guides(fill=FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g5.2


## FN, SW

# pollination treatments separately
#open
summary(model1wreused)

newdata1wreused<-expand.grid(usage=c(1:4),SeedWeight=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
se <- data.frame(Effect(c("usage"),model1wreused))
newdata1wreused <- merge(newdata1wreused,se)



#Plot


g5.3 <- ggplot(newdata1wreused,
               aes(x=usage, y=fit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits=c(0,50))+
  xlab(" ")+ ylab("Dry mass of seeds per\nseed capsule (mg)")+ 
  geom_errorbar(aes(ymin = lower, ymax = upper),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  guides(fill=FALSE) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g5.3


m3 <- grid.arrange(g5.1,g5.2,g5.3,
             nrow = 1)

ggsave("NewFigS1.svg", plot = m3, device = "svg", path = "Plots/", width = 30, height = 10, units = "cm")




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

