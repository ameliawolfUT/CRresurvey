#R code for Realistic biodiversity loss mediates ecosystem function recovery following drought

#Amelia A. Wolf, Erika S. Zavaleta and Jennifer L. Funk

## data also available in this Github repository (ameliawolfUT/CRresurvey)

#read in data
CRbiomass=read.csv(file=paste("//Users//acacia//Dropbox//CR resurvey 2017//CRresurveyR4Nov2025.csv", sep=""), stringsAsFactors=F)
#CRbiomass=read.csv("CRresurveyR16May2025.csv", sep="", stringsAsFactors=F)
CRbiomass$block<-as.character(CRbiomass$block)


library(MASS)
library(lme4)
library(lmerTest)
library(nlme)

################
#Examining the effects of the previous experimental treatments on post-drought recovery plant communities
#########
#RICHNESS ANALYSES (NATIVE, INVASIVE, TOTAL)
m1 <- lmer(sqrt(native.richness+1) ~ spnum*treat*depth+(1|block), data=CRbiomass)
summary(m1)
anova(m1)

m2 <- lmer((invasive.richness) ~ spnum*treat*depth+(1|block), data=CRbiomass) #not presented in manuscript because there is very low (<3) invasive richness
summary(m2)
anova(m2)

m3 <- lmer(sqrt(Total.richness+2) ~ spnum*treat*depth+(1|block), data=CRbiomass)
summary(m3)
anova(m3)

#########
#ENS (Effective Number of Species, e^Shannon) ANALYSES (TOTAL only)
m3a <- lmer(log(ENS) ~ spnum*treat*depth+(1|block), data=CRbiomass)
summary(m3a)
anova(m3a)

#########
#PERCENT COVER ANALYSES (NATIVE, INVASIVE, TOTAL, NFixer)
m4 <- lmer(sqrt(scaled.native.cover) ~ spnum*treat*depth+(1|block), data=CRbiomass)
summary(m4)
anova(m4)

m5 <- lmer((scaled.invasive.cover) ~ spnum*treat*depth+(1|block), data=CRbiomass)
summary(m5)
anova(m5)

m6 <- lmer((total.scaled.cover)^2 ~ spnum*treat*depth+(1|block), data=CRbiomass)#m6 has an "isSingular" warning - rerun model without the random effect of block below
summary(m6)
anova(m6) 
m6lm <- lm((total.scaled.cover)^2 ~ spnum*treat*depth, data=CRbiomass)
summary(m6lm)
anova(m6lm)

# Model with Nfixer cover needs to be fit using Negative Binomial distribution - using the raw cover data rather than scaled because negative binomial requires integer values
m7_negbin <- glmer.nb((raw.Nfixer.cover+1) ~ spnum * treat * depth + (1 | block), 
                       data = CRbiomass,
                       control = glmerControl(optimizer = "bobyqa"))
summary(m7_negbin)
Anova(m7_negbin, type = "III")

################
#PRODUCTIVITY ANALYSES (NATIVE, INVASIVE, TOTAL)
m8 <- lmer(log(totnatprod+2) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m8)
anova(m8)

m9 <- lmer(log(totprod) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m9)
anova(m9)

m10 <- lmer(log(invprod) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)#m10 has an "isSingular" warning - rerun model without the random effect of block below
summary(m10)
anova(m10) 
m10lm <- lm(log(invprod) ~ spnum*treat*depth, data=CRbiomass,na.action=na.omit)
summary(m10lm)
anova(m10lm)

#################
#NITROGEN USE EFFICIENCY (NUE) ANALYSES (NATIVE, INVASIVE, TOTAL)
m11 <- lmer(log(natNUE) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m11)
anova(m11)

m12 <- lmer((invNUE) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m12)
anova(m12)

m13 <- lmer((totNUE) ~ spnum*treat*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m13)
anova(m13)

##################
#Models examining the relationship between post-drought recovery richness and Post-drought recovery productivity/cover/NUE

############
#Post-drought cover (Native, Invasive, Total) ~ Post-drought native richness
m14<-lmer(sqrt(scaled.native.cover)~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m14)
anova(m14)

m15<-lmer(sqrt(scaled.invasive.cover)~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m15)
anova(m15)

m16<-lmer((total.scaled.cover)^2~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)#m16 has an "isSingular" warning - rerun model without the random effect of block below
summary(m16)
anova(m16)
m16lm<-lm((total.scaled.cover)^2~sqrt(native.richness+1)*depth, data=CRbiomass,na.action=na.omit)
summary(m16lm)
anova(m16lm)

#Post-drought cover (Native, Invasive, Total) ~ Post-drought total richness
m17<-lmer(sqrt(scaled.native.cover)~(sqrt(Total.richness+2))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m17)
anova(m17)

m18<-lmer(sqrt(scaled.invasive.cover)~(sqrt(Total.richness+2))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m18)
anova(m18)

m19<-lmer(total.scaled.cover^2~(sqrt(Total.richness+2))*depth+(1|block), data=CRbiomass,na.action=na.omit)#m19 has an "isSingular" warning - rerun model without the random effect of block below
summary(m19)
anova(m19)
m19lm<-lm(total.scaled.cover^2~(sqrt(Total.richness+2))*depth, data=CRbiomass,na.action=na.omit)
summary(m19lm)
anova(m19lm)

#Post-drought cover (Native, Invasive, Total) ~ Post-drought ENS
m17a<-lmer(sqrt(scaled.native.cover)~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m17a)
anova(m17a)

m18a<-lmer(sqrt(scaled.invasive.cover)~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m18a)
anova(m18a)
m18alm<-lm(sqrt(scaled.invasive.cover)~(log(ENS))*depth, data=CRbiomass,na.action=na.omit)
summary(m18alm)
anova(m18alm)

m19a<-lmer(total.scaled.cover^2~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)#m19 has an "isSingular" warning - rerun model without the random effect of block below
summary(m19a)
anova(m19a)
m19alm<-lm(total.scaled.cover^2~(log(ENS))*depth, data=CRbiomass,na.action=na.omit)
summary(m19alm)
anova(m19alm)

#Post-drought productivity (Native, Invasive, Total) ~ Post-drought native richness
m20<-lmer(log(totnatprod+2)~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m20)
anova(m20)

m21<-lmer(log(invprod)~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m21)
anova(m21)

m22<-lmer((totprod)~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m22)
anova(m22) 

#Post-drought productivity (Native, Invasive, Total) ~ Post-drought total richness
m23<-lmer(sqrt(totnatprod+1)~(sqrt(Total.richness+1))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m23)
anova(m23)

m24<-lmer(log(invprod)~(sqrt(Total.richness+1))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m24)
anova(m24)

m25<-lmer((totprod)~(sqrt(Total.richness+2))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m25)
anova(m25)

#Post-drought productivity (Native, Invasive, Total) ~ Post-drought ENS
m23a<-lmer(sqrt(totnatprod+1)~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m23a)
anova(m23a)

m24a<-lmer(log(invprod)~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m24a)
anova(m24a)

m25a<-lmer((totprod)~(log(ENS))*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m25a)
anova(m25a)

#Post-drought NUE (Native, Invasive, Total) ~ Post-drought native richness
m26 <- lmer(log(natNUE) ~ sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m26)
anova(m26)

m27<-lmer(invNUE~sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m27)
anova(m27)

m28 <- lmer((totNUE) ~ sqrt(native.richness+1)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m28)
anova(m28)

#Post-drought NUE (Native, Invasive, Total) ~ Post-drought total richness
m29 <- lmer(log(natNUE) ~ sqrt(Total.richness+2)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m29)
anova(m29)

m30<-lmer(invNUE~sqrt(Total.richness+2)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m30)
anova(m30)

m31 <- lmer((totNUE) ~ sqrt(Total.richness+2)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m31)
anova(m31)

#Post-drought NUE (Native, Invasive, Total) ~ Post-drought ENS
m29 <- lmer(log(natNUE) ~ sqrt(ENS)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m29)
anova(m29)

m30<-lmer(invNUE~sqrt(ENS)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m30)
anova(m30)

m31 <- lmer((totNUE) ~ sqrt(ENS)*depth+(1|block), data=CRbiomass,na.action=na.omit)
summary(m31)
anova(m31)
########################

#Subset data by depth and loss scenario
CRbiomassOrd <- subset(CRbiomass, treat=="Order")
CRbioOrdSha <- subset(CRbiomassOrd, depth=="Shallow")
CRbioOrdMed <- subset(CRbiomassOrd, depth=="Medium")
CRbioOrdDee <- subset(CRbiomassOrd, depth=="Deep")
CRbiomassRan <- subset(CRbiomass, treat=="Random")
CRbioRanSha <- subset(CRbiomassRan, depth=="Shallow")
CRbioRanMed <- subset(CRbiomassRan, depth=="Medium")
CRbioRanDee <- subset(CRbiomassRan, depth=="Deep")
CRbiomassSha <- subset(CRbiomass, depth=="Shallow")
CRbiomassMed <- subset(CRbiomass, depth=="Medium")
CRbiomassDee <- subset(CRbiomass, depth=="Deep")

#Order depths in proper order
CRbiomassforbox <- CRbiomass
CRbiomassforbox$depth <- factor(CRbiomassforbox$depth, c("Shallow", "Medium", "Deep"))
CRbiomassOrdbox <- CRbiomassOrd
CRbiomassOrdbox$depth <- factor(CRbiomassOrdbox$depth, c("Shallow", "Medium", "Deep"))
CRbiomassRanbox <- CRbiomassRan
CRbiomassRanbox$depth <- factor(CRbiomassRanbox$depth, c("Shallow", "Medium", "Deep"))


#Splitting out the effects of species loss scenario on post-drought recovery native richness
mord <- lmer(sqrt(native.richness+1) ~ spnum*depth+(1|block), data=CRbiomassOrd)
summary(mord)
anova(mord)
mran <- lmer(sqrt(native.richness+1) ~ spnum*depth+(1|block), data=CRbiomassRan) #mran has an "isSingular" warning - rerun model without the random effect of block below
summary(mran)
anova(mran)
mranlm <- lm(sqrt(native.richness+1) ~ spnum*depth, data=CRbiomassRan)
summary(mranlm)
anova(mranlm)

#Plots examining the significant treat*spnum interaction on post-drought recovery native richness (Figure 2)
plot(jitter(CRbiomassRan$spnum, factor = .4), CRbiomassRan$native.richness, main="Random species loss", xlab="Experimental phase target richness", ylab="Post-drought recovery native richness", ylim=c(0, 9), pch=19)
abline(lm(CRbiomassRan$native.richness~CRbiomassRan$spnum), lty = 2)
plot(jitter(CRbiomassOrd$spnum, factor = .4), (CRbiomassOrd$native.richness), main="Realistic species loss", xlab="Experimental phase target richness", ylab="Post-drought recovery native richness", ylim=c(0, 9), pch=19)
abline(lm(CRbiomassOrd$native.richness~CRbiomassOrd$spnum))

#Plots examining the significant treat*spnum interaction on post-drought recovery ENS (Figure 2)
plot(jitter(CRbiomassRan$spnum, factor = .4), CRbiomassRan$ENS, main="Random species loss", xlab="Experimental phase target richness", ylab="Post-drought effective number of species (ENS)", ylim=c(0, 5), pch=19)
abline(lm(CRbiomassRan$ENS~CRbiomassRan$spnum), lty = 2)
boxplot(ENS~spnum, data=CRbiomassRan, notch=F, xlab="Target richness", ylab="ENS")
plot(jitter(CRbiomassOrd$spnum, factor = .4), (CRbiomassOrd$ENS), main="Realistic species loss", xlab="Experimental phase target richness", ylab="Post-drought effective number of species (ENS)", ylim=c(0, 5), pch=19)
abline(lm(CRbiomassOrd$ENS~CRbiomassOrd$spnum))
boxplot(ENS~spnum, data=CRbiomassOrd, notch=F, xlab="Target richness", ylab="ENS")


#Plots for invasive productivity (Figure 3)
par(mar = c(6, 6, 4, 2) + 0.1)
plot(CRbioOrdSha$spnum, CRbioOrdSha$invprod, main="Realistic species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioOrdSha$invprod~CRbioOrdSha$spnum))
plot(CRbioOrdMed$spnum, CRbioOrdMed$invprod, main="Realistic species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioOrdMed$invprod~CRbioOrdMed$spnum), lty=2)
plot(CRbioOrdDee$spnum, CRbioOrdDee$invprod, main="Realistic species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioOrdDee$invprod~CRbioOrdDee$spnum))

par(mar = c(6, 6, 4, 2) + 0.1)
plot(CRbioRanSha$spnum, CRbioRanSha$invprod, main="Random species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioRanSha$invprod~CRbioRanSha$spnum), lty=2)
plot(CRbioRanMed$spnum, CRbioRanMed$invprod, main="Random species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioRanMed$invprod~CRbioRanMed$spnum), lty=2)
plot(CRbioRanDee$spnum, CRbioRanDee$invprod, main="Random species loss", xlab="Experimental phase target richnesss", ylab=expression("Post-recovery Invasive Productivity (g m"^-2*" yr"^-1*")"), ylim=c(150,900), pch=19)
abline(lm(CRbioRanDee$invprod~CRbioRanDee$spnum), lty=2)


#plots of Post-drought recovery richness (native, total) effects on productivity, cover, NUE (Including plots included in Figure 4 and Figure S4)
plot(CRbiomass$Total.richness, CRbiomass$totprod, xlab="Post-drought recovery total species richness", ylab="Post-drought recovery total ANPP (g m-2 yr-1)", xlim=c(0,10), ylim=c(0,1000), pch=19)
abline(lm(CRbiomass$totprod~CRbiomass$Total.richness))
plot(CRbiomass$native.richness, CRbiomass$totnatprod, xlab="Post-drought recovery native species richness", ylab="Post-drought recovery native ANPP (g m-2 yr-1)", xlim=c(0,10), ylim=c(0,200), pch=19)
abline(lm(CRbiomass$totnatprod~CRbiomass$native.richness))
plot(CRbiomass$native.richness, CRbiomass$scaled.native.cover, xlab="Post-drought recovery native species richness", ylab="Post-drought recovery native cover", pch=19)
abline(lm(CRbiomass$scaled.native.cover~CRbiomass$native.richness))
plot(CRbiomass$native.richness, CRbiomass$invprod, xlab="Post-drought recovery native species richness", ylab="Post-drought recovery invasive ANPP (g m-2 yr-1)", xlim=c(0,10), ylim=c(0,900), pch=19)
abline(lm(CRbiomass$invprod~CRbiomass$native.richness))
plot(CRbiomass$total.scaled.cover, CRbiomass$totprod, xlab="Post-drought recovery total species cover", ylab="Post-drought recovery total ANPP", pch=19)
abline(lm(CRbiomass$totprod~CRbiomass$total.scaled.cover))
plot(CRbiomass$Total.richness, CRbiomass$invNUE, xlab="Post-drought recovery total richness", ylab="Post-drought recovery invasive NUE", pch=19)
abline(lm(CRbiomass$invNUE~CRbiomass$Total.richness))
plot(CRbiomass$native.richness, CRbiomass$invNUE, xlab="Post-drought recovery native richness", ylab="Post-drought recovery invasive NUE", pch=19)
abline(lm(CRbiomass$invNUE~CRbiomass$native.richness))

#boxplot of post-drought recovery native richness as a factor of soil depth (Figure S1)
boxplot(native.richness~depth, data=CRbiomassforbox,na.action=na.omit, notch=F, xlab="Soil Depth", ylab="Post-drought recovery native richness")

#graphs of significant treat*spnum interaction effects on post-drought Nfixer cover (Figure S2)
plot(jitter(CRbiomassRan$spnum, factor = .4), CRbiomassRan$raw.Nfixer.cover, main="Random species loss", xlab="Experimental phase target diversity", ylab="Post-drought recovery N fixer cover (%)", ylim=c(0, 18), pch=19)
abline(lm(CRbiomassRan$raw.Nfixer.cover~CRbiomassRan$spnum), lty = 2)

plot(jitter(CRbiomassOrd$spnum, factor = .4), CRbiomassOrd$raw.Nfixer.cover, main="Realistic species loss", xlab="Experimental phase target diversity", ylab="Post-drought recovery N fixer cover (%)", ylim=c(0, 18), pch=19)
abline(lm(CRbiomassOrd$raw.Nfixer.cover~CRbiomassOrd$spnum), lty = 2)

#boxplots of Total NUE/visualizing significant Treatment*Depth interaction (Figure S3)
boxplot(totNUE~depth, data=CRbiomassOrdbox,na.action=na.omit, notch=F, main="Realistic species loss", xlab="Soil Depth", ylab="Post-recovery total NUE", ylim=c(30,100))
boxplot(totNUE~depth, data=CRbiomassRanbox,na.action=na.omit, notch=F, main="Randomized species loss", xlab="Soil Depth", ylab="Post-recovery total NUE", ylim=c(30,100))
