rm(list=ls())

# TO DO-----

 # turn modelling code into functions
     # 1) tps
     # 2) costs
     # 3) qol
       # 4) run model: calls 1-3 and extracts key outputs

# --------



#### SET UP ######
# Packages ------
library(dplyr)
library(survHE)
library(ggplot2)
library(cowplot)
library(scales)
library(mice)
library(survival)
library(flexsurv)
library(rms)
library(car)

# Data ------
load("Y:/CPRD for analysis.RData")
# for average characteristics

load("Y:/multiple imputated data knee/models.RData")
# models (for each mi dataset) 

load("Y:/proms with MI.RData")
# get predicted post-op qol 

#load("Y:/CPRD for analysis.RData")


# Functions -----


#### AVERAGE CHARACTERISTICS ------
# patient profiles -----

Mode <- function(x, na.rm = T) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


tkr.average.characteristics<-cprd_patients %>% 
  filter(tkr.1==1) %>% 
  summarise(median_tkr_age=round(median(tkr.1.age, na.rm=T)),
            mode_tkr_diagnosis=as.character(Mode(diagnosis_knee)),
            mode_tkr_gender=as.character(Mode(gender)),
            mode_tkr_IMD_2004_quintiles=as.character(Mode(IMD_2004_quintiles)),
            mode_tkr_RCS=as.character(Mode(tkr.1.RCS.charlson.full)),
            median_tkr_BMI=round(median(tkr.1.BMI, na.rm=T)),
            mode_tkr_smoke=as.character(Mode(tkr.1.smoke)),
            median_tkr_q1_eq5d=median(tkr.1.q1_eq5d_index, na.rm=T),
            median_tkr_year=median(tkr.1.year, na.rm=T))

median_tkr_age<-tkr.average.characteristics$median_tkr_age
mode_tkr_diagnosis<-tkr.average.characteristics$mode_tkr_diagnosis
mode_tkr_gender<-tkr.average.characteristics$mode_tkr_gender
mode_tkr_IMD_2004_quintiles<-tkr.average.characteristics$mode_tkr_IMD_2004_quintiles
mode_tkr_RCS<-tkr.average.characteristics$mode_tkr_RCS
mode_tkr_smoke<-tkr.average.characteristics$mode_tkr_smoke
median_tkr_BMI<-tkr.average.characteristics$median_tkr_BMI
median_tkr_q1_eq5d<-tkr.average.characteristics$median_tkr_q1_eq5d
median_tkr_year<-tkr.average.characteristics$median_tkr_year

 characteristics<-
  rbind(
  expand.grid(age=median_tkr_age, #c(60,70),
              gender=mode_tkr_gender,
              diagnosis=mode_tkr_diagnosis,
              IMD=as.character(mode_tkr_IMD_2004_quintiles),
              RCS=mode_tkr_RCS,
              BMI=median_tkr_BMI,
              smoke=mode_tkr_smoke,
              q1_eq5d=median_tkr_q1_eq5d,
              year=median_tkr_year))

characteristics$group<-seq(1, length(characteristics$age))

characteristics$diagnosis<-as.character(characteristics$diagnosis)
characteristics$gender<-as.character(characteristics$gender)
characteristics$RCS<-as.character(characteristics$RCS)
characteristics$IMD<-as.character(characteristics$IMD)
characteristics$smoke<-as.character(characteristics$smoke)











# relative improvements------
# for many and combined 
qol.improvement<-seq(1,1.05, by=0.0025)     # relative improvement in post-op qol
revision.reduction<-seq(0, 1, by=0.025) # relative reduction in risk of revision

scenarios<-expand.grid(qol.improvement=qol.improvement,
            revision.reduction=revision.reduction)

# first scenario: usual care 
# *important* reference when calculaing incremental effects
scenarios<-scenarios %>% 
  arrange(desc(revision.reduction),qol.improvement)

scenarios$scenario.id<-1:length(scenarios[,1])

# run model -------
# number of mi -----
n.mi<-1#2#5#0
# number of sims ----
sims<-10#0   # number of bootstraps


##### GET TRANSITION PROBABILITIES -----
# survival probabilities ------
# get tkr survival curves ------
  # risk of rev - from surgery
  # risk of death - from surgery
  # risk of death- from diagnosis




# bootstrap survival curves  -------

bs.tkr_revision.risk<-list()
for(i in 1:n.mi) {
bs.tkr_revision.risk[[i]]<-normboot.flexsurvreg(mi.models.tkr_revision[[i]], 
            B=sims, 
            newdata = data.frame(age=characteristics$age, 
                             diagnosis=characteristics$diagnosis, 
                             gender=characteristics$gender,
                             RCS=characteristics$RCS, 
                             IMD=characteristics$IMD, 
                             BMI=characteristics$BMI, 
                             smoke=characteristics$smoke))}

bs.model.mortality.tkr<-list()
for(i in 1:n.mi) {
bs.model.mortality.tkr[[i]]<-normboot.flexsurvreg(mi.models.tkr.to.death[[i]], 
            B=sims, 
            newdata = data.frame(age=characteristics$age, 
                             diagnosis=characteristics$diagnosis, 
                             gender=characteristics$gender,
                             RCS=characteristics$RCS, 
                             IMD=characteristics$IMD, 
                             BMI=characteristics$BMI, 
                             smoke=characteristics$smoke))}

bs.model.mortality.knee<-list()
for(i in 1:n.mi) {
bs.model.mortality.knee[[i]]<-normboot.flexsurvreg(mi.models.diagnosis.to.death[[i]], 
            B=sims, 
            newdata = data.frame(age=characteristics$age, 
                             diagnosis_knee=characteristics$diagnosis, 
                             gender=characteristics$gender,
                             charlson=characteristics$RCS, 
                             IMD=characteristics$IMD, 
                             bmi_diagnosis_knee=characteristics$BMI, 
                             smoke_diagnosis_knee=characteristics$smoke))}


# deterministic probs  -----
all.m.probs<-list()
#using.mi<-1
for(using.mi in 1:n.mi) {
# unrevised_revised probability 
survival.tkr_revision<-summary(mi.models.tkr_revision[[using.mi]],
                               t = seq(0,50), 
        newdata = data.frame(age=characteristics$age, 
                             diagnosis=characteristics$diagnosis, 
                             gender=characteristics$gender,
                             RCS=characteristics$RCS, 
                             IMD=characteristics$IMD, 
                             BMI=characteristics$BMI, 
                             smoke=characteristics$smoke), 
                             ci=FALSE)

names(survival.tkr_revision)<-characteristics$group
survival.tkr_revision <-mapply(`[<-`, survival.tkr_revision, 
                               'group', value = names(survival.tkr_revision), SIMPLIFY = FALSE)
survival.tkr_revision <-bind_rows(survival.tkr_revision) %>% 
  mutate(group=as.integer(group))
# we now have the survival curves for each 
# set of characteristics

# need to transform survival curves to
# transition probabilities
tkr.revision.risk<-NULL
for (g in 1:max(characteristics$group)) {
  working.group<-g
for (cycle in seq(1:50)) {
  t.using<-cycle
  tp.current<-1-(survival.tkr_revision[which(survival.tkr_revision$group==working.group & survival.tkr_revision$time==t.using), 2]/
                   survival.tkr_revision[which(survival.tkr_revision$group==working.group & survival.tkr_revision$time==t.using-1), 2])
  tkr.revision.risk<-rbind(tkr.revision.risk, 
                           data.frame(time=t.using, 
                                      prob.revision=tp.current, 
                                      group=working.group))
  }}
rm(survival.tkr_revision)
# now have transition probabilities for 
# unrevised to revision

# only want them up to age 99
tkr.revision.risk<-tkr.revision.risk %>% 
  left_join(characteristics %>% select(group, age), by=c("group"))

# only up to age 99
tkr.revision.risk$working.age<-tkr.revision.risk$age+(tkr.revision.risk$time-1)
tkr.revision.risk<-tkr.revision.risk %>% 
  filter(working.age<100) %>% 
  select(-working.age, -age)





# unrevised_dead probability 
 # first 10 years based on tkr to death
 # subsequent years based on diagnosis to death

# 10 year post op mortality risk
survival.tkr_mortality<-summary(mi.models.tkr.to.death[[using.mi]], 
                                t = seq(0,10), 
                               newdata = data.frame(age=characteristics$age, 
                                                    diagnosis=characteristics$diagnosis, 
                                                    gender=characteristics$gender,
                                                    RCS=characteristics$RCS, 
                                                    IMD=characteristics$IMD, 
                                                    BMI=characteristics$BMI, 
                                                    smoke=characteristics$smoke), ci=FALSE)
names(survival.tkr_mortality)<-characteristics$group
survival.tkr_mortality <-mapply(`[<-`, survival.tkr_mortality, 
                               'group', value = names(survival.tkr_mortality), SIMPLIFY = FALSE)
survival.tkr_mortality <-bind_rows(survival.tkr_mortality) %>% 
  mutate(group=as.integer(group))

# to transition probabilities
tkr.mortality.risk<-NULL
for (g in 1:max(characteristics$group)) {
  working.group<-g
  for (cycle in seq(1:10)) {
    t.using<-cycle
    tp.current<-1-(survival.tkr_mortality[which(survival.tkr_mortality$group==working.group & survival.tkr_mortality$time==t.using), 2]/
                     survival.tkr_mortality[which(survival.tkr_mortality$group==working.group & survival.tkr_mortality$time==t.using-1), 2])
    tkr.mortality.risk<-rbind(tkr.mortality.risk, 
                             data.frame(time=t.using, 
                                        prob.mortality=tp.current, 
                                        group=working.group)) }}
rm(survival.tkr_mortality)
# now we have survival curve for 10 years following op
# want to add subsequent risk of death

# subsequent risk of death 
# n.b this model uses age as the time scale
# starts from age of 60
survival.diagnosis_mortality<-summary(mi.models.diagnosis.to.death[[using.mi]],  
                                      t = seq(0,50), #t=0 is age 60 
                                newdata = data.frame(#age=characteristics$age, #nb not used for calculation
                                                     diagnosis_knee=characteristics$diagnosis, 
                                                     gender=characteristics$gender,
                                                     charlson=characteristics$RCS, 
                                                     IMD=characteristics$IMD, 
                                                     bmi_diagnosis_knee=characteristics$BMI, 
                                                     smoke_diagnosis_knee=characteristics$smoke), ci=FALSE)
names(survival.diagnosis_mortality)<-characteristics$group
survival.diagnosis_mortality <-mapply(`[<-`, survival.diagnosis_mortality, 
                               'group', value = names(survival.diagnosis_mortality), SIMPLIFY = FALSE)
survival.diagnosis_mortality <-bind_rows(survival.diagnosis_mortality) %>% 
  mutate(group=as.integer(group))
# to transition probabilities
diagnosis.mortality.risk<-NULL
for (g in 1:max(characteristics$group)) {
  working.group<-g
  for (cycle in seq(1:50)) {
    t.using<-cycle
    tp.current<-1-(survival.diagnosis_mortality[which(survival.diagnosis_mortality$group==working.group & survival.diagnosis_mortality$time==t.using), 2]/
                     survival.diagnosis_mortality[which(survival.diagnosis_mortality$group==working.group & survival.diagnosis_mortality$time==t.using-1), 2])
    diagnosis.mortality.risk<-rbind(diagnosis.mortality.risk, 
                              data.frame(time=t.using, 
                                         prob.mortality=tp.current, 
                                         group=working.group)) }}
rm(survival.diagnosis_mortality)
# we have transition probabilities for background risk of death
# these are based on age (rather than time since op)
# time 1 is risk from age 60 to 61

# want to keep those from start age plus 10
# i.e time point where we switch to 
# background mortality

diagnosis.mortality.risk$surv.age<-diagnosis.mortality.risk$time+59

diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
  left_join(characteristics %>% select(age, group),
            by=c("group")) 

diagnosis.mortality.risk$t.x<-(diagnosis.mortality.risk$surv.age)-diagnosis.mortality.risk$age
# want from t.x=11
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
             filter(t.x>=11)
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
             mutate(time=t.x) %>% 
             select(-c(surv.age,t.x))
# now t for diagnosis_mortality is time since op, starting from 11 years post-op

# up to age 99
diagnosis.mortality.risk$working.age<-diagnosis.mortality.risk$age+(diagnosis.mortality.risk$time-1)
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
  filter(working.age<100) %>% 
  select(-working.age, -age)

# combine risks of death
tkr.mortality.risk<-tkr.mortality.risk %>% rbind(diagnosis.mortality.risk)
tkr.mortality.risk<-tkr.mortality.risk %>% arrange(group)
rm(diagnosis.mortality.risk)
  


###
# unrevised_unrevised probability 
m.probs<- left_join(tkr.revision.risk,tkr.mortality.risk,
              by=c("time","group"))
m.probs<-m.probs %>% 
  rename(prob.unrevised_revision=prob.revision)
m.probs<-m.probs %>% 
  rename(prob.unrevised_dead=prob.mortality)

m.probs$mi<-using.mi

name<-paste0(" n.mi:",using.mi)
all.m.probs[[name]]<-m.probs
}

#as one long dataframe
m.probs.deterministic<- plyr::ldply(all.m.probs, 
                      data.frame, .id=NULL)

m.probs.deterministic$sim<-"deterministic"
rm(all.m.probs)
# probabilistic probs -----
# risk of revision

tp.tkr.revision.probabilistic<-list()
for(using.mi in 1:n.mi) {


survival.tkr_revision<-NULL
for (j in 1:max(characteristics$group)) {
for (t in 0:50) {
working.group<-j
current.time<-t

if (length(characteristics$group)==1) {
  g<-as.data.frame(bs.tkr_revision.risk[[using.mi]])
}

if (length(characteristics$group)>1) {
g<-as.data.frame(bs.tkr_revision.risk[[using.mi]][[working.group]]) # 1st is mi, 2nd is group
}



current.survival<-data.frame(time=current.time,
                             survival.tkr_revision=psurvspline(current.time,  
            gamma = matrix(c(g[,1],  #gamma0
                      g[,2],  #gamma1
                      g[,3], #gamma2
                      g[,4]), #gamma3
                      nrow = sims, ncol = 4, byrow = F), 
            knots = mi.models.tkr_revision[[i]]$knots,
            lower.tail = FALSE),
            group=working.group)
current.survival$sim<-1:length(current.survival$time)

survival.tkr_revision<-rbind(survival.tkr_revision, current.survival)
}}
#survival.tkr_revision


#to transition probabilities
tkr.revision.risk<-survival.tkr_revision %>% 
  group_by(group, sim) %>% 
  arrange(sim,group,time) %>% 
  mutate(prob.revision = 1-(survival.tkr_revision/ lag(survival.tkr_revision))) %>% 
  select(time,prob.revision, group, sim) %>% 
  filter(time!=0) %>% 
  ungroup()
  
rm(survival.tkr_revision)
tkr.revision.risk<-tkr.revision.risk %>%  left_join(characteristics %>% select(age, group)
                                                    , by=c("group"))

# only up to age 99
tkr.revision.risk$working.age<-tkr.revision.risk$age+(tkr.revision.risk$time-1)
tkr.revision.risk<-tkr.revision.risk %>% 
  filter(working.age<100) %>% 
  select(-working.age, age)

tkr.revision.risk$mi<-using.mi
tp.tkr.revision.probabilistic[[using.mi]]<-tkr.revision.risk
}

tp.tkr.revision.probabilistic<- plyr::ldply(
  tp.tkr.revision.probabilistic, 
                      data.frame, .id=NULL)
tp.tkr.revision.probabilistic<-tp.tkr.revision.probabilistic %>% 
  select(-age)

#risk of death
# 10 year post op mortality risk

tp.tkr.death.probabilistic<-list()
for(using.mi in 1:using.mi) {

survival.tkr_mortality<-NULL
for (j in 1:max(characteristics$group)) {
for (t in 0:10) {
working.group<-j
current.time<-t


if (length(characteristics$group)==1) {
  g<-as.data.frame(bs.model.mortality.tkr[[using.mi]])
}

if (length(characteristics$group)>1) {
  g<-as.data.frame(bs.model.mortality.tkr[[using.mi]][[working.group]])
}


current.survival<-data.frame(time=current.time,
                             survival.tkr_mortality=psurvspline(current.time,  
            gamma = matrix(c(g[,1],  #gamma0
                      g[,2],  #gamma1
                      g[,3],  #gamma2
                      g[,4]), #gamma3
                      nrow = sims, ncol = 4, byrow = F), 
            knots = mi.models.tkr.to.death[[i]]$knots,
            lower.tail = FALSE),
            group=working.group)
current.survival$sim<-1:length(current.survival$time)

survival.tkr_mortality<-rbind(survival.tkr_mortality, current.survival)
}}

#to transition probabilities
tkr.mortality.risk<-survival.tkr_mortality %>% 
  group_by(group, sim) %>% 
  arrange(sim,group,time) %>% 
  mutate(prob.mortality = 1-(survival.tkr_mortality/ lag(survival.tkr_mortality))) %>% 
  select(time,prob.mortality, group, sim) %>% 
  filter(time!=0) %>% 
  ungroup()
  
rm(survival.tkr_mortality)
tkr.mortality.risk<-tkr.mortality.risk %>% 
  left_join(characteristics %>% select(group,age), by=c("group"))

# subsequent risk of death 
survival.diagnosis_mortality<-NULL
for (j in 1:max(characteristics$group)) {
for (t in 0:50) {
working.group<-j
current.time<-t

if (length(characteristics$group)==1) {
  g<-as.data.frame(bs.model.mortality.knee[[using.mi]])
}

if (length(characteristics$group)>1) {
  g<-as.data.frame(bs.model.mortality.knee[[using.mi]][[working.group]])
}



current.survival<-data.frame(time=current.time,
                             survival.tkr_mortality=pgompertz(current.time,  
            shape=g[,1], 
            rate= g[,2], 
            lower.tail = FALSE),
            group=working.group)
current.survival$sim<-1:length(current.survival$time)

survival.diagnosis_mortality<-rbind(survival.diagnosis_mortality, current.survival)
}}

# to transition probabilities
diagnosis.mortality.risk<-survival.diagnosis_mortality %>% 
  group_by(group, sim) %>% 
  arrange(sim,group,time) %>% 
  mutate(prob.mortality = 1-(survival.tkr_mortality/ lag(survival.tkr_mortality))) %>% 
  select(time,prob.mortality, group, sim) %>% 
  filter(time!=0) %>% 
  ungroup()

rm(survival.diagnosis_mortality)
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
  left_join(characteristics %>% select(group,age), by=c("group"))

# combine mortality for each group 
# want first ten years from tkr mortality
# then want next ten years- choice of 'time' depends on starting age
# eg for 70 year old we want from time 20 onwards (e.g. from age 80....)

diagnosis.mortality.risk$surv.age<-diagnosis.mortality.risk$time+59
diagnosis.mortality.risk$t.x<-(diagnosis.mortality.risk$surv.age)-diagnosis.mortality.risk$age
# want from t.x=11
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
             filter(t.x>=11)
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
             mutate(time=t.x) %>% 
             select(-c(surv.age,t.x))
# now t for diagnosis_mortality is time since op, starting from 11 years post-op


# up to age 99
diagnosis.mortality.risk$working.age<-diagnosis.mortality.risk$age+(diagnosis.mortality.risk$time-1)
diagnosis.mortality.risk<-diagnosis.mortality.risk %>% 
  filter(working.age<100) %>% 
  select(-working.age, age)


# combine risks of death
tkr.mortality.risk<-tkr.mortality.risk %>% rbind(diagnosis.mortality.risk)
tkr.mortality.risk<-tkr.mortality.risk %>% arrange(group, sim)
rm(diagnosis.mortality.risk)

tkr.mortality.risk$mi<-using.mi

tp.tkr.death.probabilistic[[using.mi]]<-tkr.mortality.risk
}


tp.tkr.death.probabilistic<- plyr::ldply(
  tp.tkr.death.probabilistic, 
                      data.frame, .id=NULL)

tp.tkr.death.probabilistic<-tp.tkr.death.probabilistic %>% 
  select(-age)


###
# unrevised_unrevised probability 
m.probs<- left_join(tp.tkr.revision.probabilistic,
                    tp.tkr.death.probabilistic,
              by=c("time","group", "sim", "mi"))
m.probs<-m.probs %>% 
  rename(prob.unrevised_revision=prob.revision)
m.probs<-m.probs %>% 
  rename(prob.unrevised_dead=prob.mortality)



m.probs.probabilistic<-m.probs
rm(m.probs)

# combine----
# names(m.probs.deterministic)
# names(m.probs.probabilistic)
m.probs.probabilistic<-m.probs.probabilistic %>% 
  select(time,prob.unrevised_revision,group,                 
         prob.unrevised_dead,mi, sim )

m.probs<-rbind(m.probs.deterministic,
      m.probs.probabilistic)
rm(m.probs.deterministic,
      m.probs.probabilistic)


### GET QOL  ----
# for models
dd <- datadist(tkr.proms)
options(datadist="dd")

# eq5d 
preds<-NULL
for(i in 1:n.mi) {
using.mi<-i
using.tkr.proms.imp<-complete(tkr.proms.imp, using.mi)
using.tkr.proms.imp$tkr.1.smoke<-as.character(using.tkr.proms.imp$tkr.1.smoke)

# expected

m<-ols(tkr.1.q2_eq5d_index ~ diagnosis + 
    rcs(tkr.1.age, 3) + gender + tkr.1.RCS.charlson.ra.omitted + 
    IMD_2004_quintiles + rcs(tkr.1.q1_eq5d_index, 3) + tkr.1.BMI + 
    tkr.1.smoke,
    x=TRUE, y=TRUE,
    data=using.tkr.proms.imp)


working.expected.preds<-characteristics %>% 
  mutate(diagnosis=as.character(ifelse(diagnosis=="k_ost", "OA", "RA"))) %>% 
  rename(tkr.1.age=age,
         tkr.1.RCS.charlson.ra.omitted=RCS,
         IMD_2004_quintiles=IMD,
         tkr.1.q1_eq5d_index=q1_eq5d,
         tkr.1.BMI=BMI,
         tkr.1.smoke=smoke)
working.expected.preds<-working.expected.preds %>% 
  mutate(pred.tkr.1.q2_eq5d_index=predict(m,newdata=working.expected.preds)) %>% 
  mutate(mi=using.mi) %>% 
  mutate(sim="deterministic") %>% 
  select(group,pred.tkr.1.q2_eq5d_index,mi, sim)


# bstraps
working.bstrap.preds<-characteristics %>% 
  mutate(diagnosis=as.character(ifelse(diagnosis=="k_ost", "OA", "RA"))) %>% 
  rename(tkr.1.age=age,
         tkr.1.RCS.charlson.ra.omitted=RCS,
         IMD_2004_quintiles=IMD,
         tkr.1.q1_eq5d_index=q1_eq5d,
         tkr.1.BMI=BMI,
         tkr.1.smoke=smoke)

regressAndPredict <- function( dat ) {
  model<-ols(tkr.1.q2_eq5d_index ~ diagnosis + 
    rcs(tkr.1.age, 3) + gender + tkr.1.RCS.charlson.ra.omitted + 
    IMD_2004_quintiles + rcs(tkr.1.q1_eq5d_index, 3) + tkr.1.BMI + 
    tkr.1.smoke,
    x=TRUE, y=TRUE,
    data=dat)
  
    data.frame(pred.tkr.1.q2_eq5d_index= predict(m, working.bstrap.preds),
               group=working.bstrap.preds$group,
               mi=using.mi)

      
}

regressAndPredict(using.tkr.proms.imp)
working.bstrap.preds<-plyr::rdply(sims, 
             regressAndPredict(using.tkr.proms.imp[sample(seq(length(using.tkr.proms.imp[,1])),
                                          replace=TRUE) ,]),
             .id="sim")

working.bstrap.preds<-working.bstrap.preds %>% 
                      select(pred.tkr.1.q2_eq5d_index,
                             group, mi, sim)

preds<-rbind(preds,
  rbind(working.expected.preds,working.bstrap.preds))

}

qol<-preds %>% mutate(group=as.numeric(group))
qol<-qol %>% left_join(characteristics, by="group")





### GET COSTS ----
## tkr primary cost ------
# data ----
tkr.costs<-cprd_patients %>% 
 filter(tkr.1==1) %>% 
 filter(!is.na(tkr.1.reimbursement.cost))  %>% 
 mutate(diagnosis=ifelse(diagnosis_knee=="rheum", "RA", "OA")) %>% 
 select(tkr.1.age, gender, 
         diagnosis,
          tkr.1.RCS.charlson.ra.omitted,
          tkr.1.BMI,
          tkr.1.smoke,
          IMD_2004_quintiles,
          tkr.1.reference.cost)
# reduce rcs categories  -----
table(tkr.costs$tkr.1.RCS.charlson.ra.omitted)
# going to combine 1, 2 and 3+ 
tkr.costs$tkr.1.RCS.charlson.ra.omitted<-ifelse(tkr.costs$tkr.1.RCS.charlson.ra.omitted=="0", "0", "1+")  
table(tkr.costs$tkr.1.RCS.charlson.ra.omitted)

# knots for models ####
knots<-c(0,3,4,5,6,7,8)

# when zero, linear relationship
# will get lowest aic and bic for each event 

# age knots for tkr -----
tkr.costs %>% 
   ggplot()+
  geom_point(aes(tkr.1.reference.cost, tkr.1.age))+ 
  geom_smooth(aes(tkr.1.reference.cost, tkr.1.age),
              method="lm")

tkr.age.knots<-NULL # going to add to this
for (i in 1:length(knots)) {
current.age.knots<-knots[i]

if (current.age.knots==0){ 
glm<-glm(tkr.1.reference.cost~
            tkr.1.age, 
                      family = Gamma(link="log"),
                      data=subset(cprd_patients, tkr.1==1))

} else if (current.age.knots>0) {

glm<-glm(tkr.1.reference.cost~
           rcs(tkr.1.age,current.age.knots), 
                      family = Gamma(link="log"),
                      data=subset(cprd_patients, tkr.1==1))

} 

aic<-glm$aic
bic<-BIC(glm)
tkr.age.knots<-rbind(tkr.age.knots,
                     data.frame(age.knots=current.age.knots,
                                   aic=aic, bic=bic))

}
tkr.age.knots

#by aic
tkr.aic.age.knots<-tkr.age.knots[which.min(tkr.age.knots$aic),1] #age knots
tkr.aic.age.knots
tkr.aic.year.knots<-tkr.age.knots[which.min(tkr.age.knots$aic),2] #tkr.1.year knots
#by bic
tkr.bic.age.knots<-tkr.age.knots[which.min(tkr.age.knots$bic),1] #age knots
tkr.bic.age.knots
tkr.bic.year.knots<-tkr.age.knots[which.min(tkr.age.knots$bic),1] #dmidate_year knots








# BMI knots for tkr -----
tkr.costs %>% 
   ggplot()+
  geom_point(aes(tkr.1.reference.cost, tkr.1.BMI))+ 
  geom_smooth(aes(tkr.1.reference.cost, tkr.1.BMI),
              method="lm")

tkr.BMI.knots<-NULL # going to add to this
for (i in 1:length(knots)) {
current.BMI.knots<-knots[i]

if (current.BMI.knots==0){ 
glm<-glm(tkr.1.reference.cost~tkr.1.BMI, 
                      family = Gamma(link="log"),
                      data=subset(cprd_patients, tkr.1==1))

} else if (current.BMI.knots>0) {

glm<-glm(tkr.1.reference.cost~
           rcs(tkr.1.BMI,current.BMI.knots), 
                      family = Gamma(link="log"),
                      data=subset(cprd_patients, tkr.1==1))

} 

aic<-glm$aic
bic<-BIC(glm)
tkr.BMI.knots<-rbind(tkr.BMI.knots,
                     data.frame(BMI.knots=current.BMI.knots,
                                   aic=aic, bic=bic))

}
tkr.BMI.knots

#by aic
tkr.aic.BMI.knots<-tkr.BMI.knots[which.min(tkr.BMI.knots$aic),1] #BMI knots
tkr.aic.BMI.knots
tkr.aic.year.knots<-tkr.BMI.knots[which.min(tkr.BMI.knots$aic),2] #tkr.1.year knots
#by bic
tkr.bic.BMI.knots<-tkr.BMI.knots[which.min(tkr.BMI.knots$bic),1] #BMI knots
tkr.bic.BMI.knots
tkr.bic.year.knots<-tkr.BMI.knots[which.min(tkr.BMI.knots$bic),1] #dmidate_year knots









# model -----
# knots based on bic
tkr.bic.age.knots
tkr.bic.BMI.knots

# # simple model (not including bmi or smoking) -----
# m<-glm(tkr.1.reference.cost~diagnosis+
#                           tkr.1.age+gender+
#                         tkr.1.RCS.charlson.ra.omitted + IMD_2004_quintiles, 
#                       family = Gamma(link="log"),
#                       data=tkr.costs)
# #plot(m)
# #exp(glm.tkr.reference_cost$coefficients)
# #exp(confint(glm.tkr.reference_cost))
# 
# 
# 
# 
# 
# multiple imputation: tkr ------
#n.mi<-5#0
#predictor matrix
p<-matrix(0, length(names(tkr.costs)), 
          length(names(tkr.costs)))
dimnames(p)<-list(from = names(tkr.costs), 
                       to = names(tkr.costs))
names(tkr.costs)

# predictors 
# going to use all
predictor.vars<-c("tkr.1.age", "gender",                                 
   "diagnosis",                              
  "tkr.1.RCS.charlson.ra.omitted", "tkr.1.BMI",                            
  "tkr.1.smoke", "IMD_2004_quintiles",                     
  "tkr.1.reference.cost")

# vars to impute (characteristics and pre-op scores)
to.impute.vars<-c( "tkr.1.BMI", "tkr.1.smoke","IMD_2004_quintiles")


for(i in 1:length(to.impute.vars)) {
p[to.impute.vars[i],                                              # variable being imputed
   predictor.vars[predictor.vars != to.impute.vars[i]]]<-1
}

# impute
tkr.costs.imp <- mice(tkr.costs, predictorMatrix = p,
                     #method = method.imp,
                      #ridge = 1e-04,
                     # donors = 10L,
                          m = n.mi, maxit = 5, seed=123,
                      pri = FALSE) # increase m for more imputations
#contrasts(tkr.costs.imp$data$IMD_2004_quintiles) <- contrasts(tkr.costs.imp$IMD_2004_quintiles)
#contrasts(tkr.costs.imp$data$tkr.1.smoke) <- contrasts(tkr.proms.imp$tkr.1.smoke)


# check imputations
tkr.costs.imp$meth # these are the methods used
plot(tkr.costs.imp) # this shows convergence

densityplot(tkr.costs.imp,~tkr.1.BMI)
#stripplot(tkr.costs.imp, tkr.1.BMI, pch = 20, cex = 1.2)



# to check that no values below zero
# n.b. have used pmm so should be fine- we're using 'donors'
# a<-complete(tkr.proms.imp, action="long")
# min(a$tkr.1.BMI)












# mi models ----
# combined.tkr.mi.model<-with(tkr.costs.imp,
#   glm(tkr.1.reference.cost~diagnosis+
#                           tkr.1.age+gender+
#                         tkr.1.RCS.charlson.ra.omitted + IMD_2004_quintiles+
#                         tkr.1.BMI+tkr.1.smoke, 
#                       family = Gamma(link="log")))
# pool(combined.tkr.mi.model)

#each mi model
each.tkr.mi.model<-list()
for(i in 1:n.mi) {
each.tkr.mi.model[[i]]<-glm(tkr.1.reference.cost~diagnosis+
                          tkr.1.age+gender+
                          tkr.1.RCS.charlson.ra.omitted + 
                          IMD_2004_quintiles+
                          tkr.1.BMI+tkr.1.smoke, 
                      family = Gamma(link="log"),
                      data=complete(tkr.costs.imp, i))
}





# get predictions -----

# get predicted costs
using.characteristics<-
  characteristics %>% 
  rename(tkr.1.year=year,
         tkr.1.age=age,
         tkr.1.RCS.charlson.ra.omitted=RCS,
         IMD_2004_quintiles=IMD,
         tkr.1.BMI=BMI,
         tkr.1.smoke=smoke) %>% 
 mutate(diagnosis=ifelse(diagnosis=="rheum", "RA", "OA"))


tkr.primary.costs<-NULL
for(i in 1:n.mi) {
using.mi<-i

deterministic.costs<-data.frame(predicted.tkr.1.cost=exp(predict(each.tkr.mi.model[[using.mi]],
                          using.characteristics)),
           group=using.characteristics$group,
           mi=using.mi,
           sim="deterministic")
                       
# bstrapped
each.tkr.mi.model[[using.mi]]$call
#each.tkr.mi.model[[using.mi]]$data
using.data<-each.tkr.mi.model[[using.mi]]$data

regressAndPredict <- function( dat ) {
  model <- glm(tkr.1.reference.cost~diagnosis +
                 tkr.1.age + 
                 gender + 
                 tkr.1.RCS.charlson.ra.omitted + 
                 IMD_2004_quintiles + 
                 tkr.1.BMI + tkr.1.smoke, 
                 family = Gamma(link="log"),
                 data=dat)
     
cbind(data.frame(predicted.tkr.1.cost=
                   exp(predict(model,
                               using.characteristics))),
           group=using.characteristics$group)
      
}

regressAndPredict(using.data)
prob.costs<-plyr::rdply(sims, 
             regressAndPredict(using.data[sample(seq(length(using.data$diagnosis)),
                                          replace=TRUE) ,]),
             .id="sim")
prob.costs<-prob.costs %>% 
  mutate(sim=as.character(sim), 
         mi=using.mi) %>% 
  select(predicted.tkr.1.cost,group,mi,sim)

working.costs<-rbind(deterministic.costs, prob.costs )
tkr.primary.costs<-rbind(tkr.primary.costs,working.costs)
}

##
##### ####
# tkr revision cost----
#load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/model tkr_revision costs.RData")
# data ----
#load("Y:/CPRD for analysis.RData")

tkr_revision.costs<-cprd_patients %>% 
 filter(tkr_revision.post.diagnosis.1==1) %>% 
 filter(!is.na(tkr_revision.post.diagnosis.1.reimbursement.cost))  %>% 
 rename(tkr_revision.1.age=tkr_revision.post.diagnosis.1.age,
        tkr_revision.1.RCS.charlson.ra.omitted=tkr_revision.post.diagnosis.1.RCS.charlson.ra.omitted,
          tkr_revision.1.reference.cost=tkr_revision.post.diagnosis.1.reference.cost)  %>% 
 mutate(diagnosis=ifelse(diagnosis_knee=="rheum", "RA", "OA")) %>% 
 select(tkr_revision.1.age, gender, 
         diagnosis,
          tkr_revision.1.RCS.charlson.ra.omitted,
          IMD_2004_quintiles,
          tkr_revision.1.reference.cost)
# reduce rcs categories  -----
table(tkr_revision.costs$tkr_revision.1.RCS.charlson.ra.omitted)
# going to combine 1, 2 and 3+ 
tkr_revision.costs$tkr_revision.1.RCS.charlson.ra.omitted<-ifelse(tkr_revision.costs$tkr_revision.1.RCS.charlson.ra.omitted=="0", "0", "1+")  
table(tkr_revision.costs$tkr_revision.1.RCS.charlson.ra.omitted)

# knots for models ####
knots<-c(0,3,4,5,6,7,8)

# when zero, linear relationship
# will get lowest aic and bic for each event 

# age knots for tkr_revision -----
tkr_revision.costs %>% 
   ggplot()+
  geom_point(aes(tkr_revision.1.reference.cost, tkr_revision.1.age))+ 
  geom_smooth(aes(tkr_revision.1.reference.cost, tkr_revision.1.age),
              method="lm")

tkr_revision.age.knots<-NULL # going to add to this
for (i in 1:length(knots)) {
current.age.knots<-knots[i]

if (current.age.knots==0){ 
glm<-glm(tkr_revision.1.reference.cost~scale(tkr_revision.1.age), 
                      family = Gamma(link="log"),
                      data=tkr_revision.costs)

} else if (current.age.knots>0) {

glm<-glm(tkr_revision.1.reference.cost~
           rcs(tkr_revision.1.age,current.age.knots), 
                      family = Gamma(link="log"),
                      data=tkr_revision.costs)

} 

aic<-glm$aic
bic<-BIC(glm)
tkr_revision.age.knots<-rbind(tkr_revision.age.knots,
                     data.frame(age.knots=current.age.knots,
                                   aic=aic, bic=bic))

}
tkr_revision.age.knots

#by aic
tkr_revision.aic.age.knots<-tkr_revision.age.knots[which.min(tkr_revision.age.knots$aic),1] #age knots
tkr_revision.aic.age.knots
tkr_revision.aic.year.knots<-tkr_revision.age.knots[which.min(tkr_revision.age.knots$aic),2] #tkr_revision.1.year knots
#by bic
tkr_revision.bic.age.knots<-tkr_revision.age.knots[which.min(tkr_revision.age.knots$bic),1] #age knots
tkr_revision.bic.age.knots
tkr_revision.bic.year.knots<-tkr_revision.age.knots[which.min(tkr_revision.age.knots$bic),1] #dmidate_year knots








# model -----
# knots based on bic
tkr_revision.bic.age.knots

#  model 
glm.tkr_revision.reference_cost<-glm(tkr_revision.1.reference.cost~diagnosis+
                          tkr_revision.1.age+gender+
                        tkr_revision.1.RCS.charlson.ra.omitted + 
                          IMD_2004_quintiles, 
                      family = Gamma(link="log"),
                      data=tkr_revision.costs)
#plot(glm.tkr_revision.reference_cost)
#exp(glm.tkr_revision.reference_cost$coefficients)
#exp(confint(glm.tkr_revision.reference_cost))







# predictions -----
# need for all ages
using.characteristics<-
  characteristics 

# need to expand to go up to age 100 for each
time<-NULL
#for(using.mi in 1:n.mi) {
for(g in 1:max(characteristics$group)) {
working.group<-g
a<-data.frame(time=1:max(m.probs$time[m.probs$group==working.group]), 
           group=working.group)
time<-rbind(time,a)
}#}


using.characteristics<-using.characteristics %>% 
  left_join(time, by="group")

using.characteristics<-using.characteristics %>% 
  mutate(working.age=age+time)



# get predicted costs
# nb revision cost models not mi
using.characteristics<-
  using.characteristics %>% 
  rename(tkr_revision.1.year=year,
         tkr_revision.1.age=working.age,
         tkr_revision.1.RCS.charlson.ra.omitted=RCS,
         IMD_2004_quintiles=IMD,
         tkr_revision.1.BMI=BMI,
         tkr_revision.1.smoke=smoke) %>% 
 mutate(diagnosis=ifelse(diagnosis=="rheum", "RA", "OA"))




deterministic.costs<-data.frame(
                       exp(data.frame(
                          predicted.tkr_revision.cost=predict(glm.tkr_revision.reference_cost,
                          using.characteristics))),
                          group=using.characteristics$group,
                          time=using.characteristics$time,
                          sim="deterministic")


# bstrapped
glm.tkr_revision.reference_cost$call
using.data<-glm.tkr_revision.reference_cost$data

regressAndPredict <- function( dat ) {
  model<-glm(tkr_revision.1.reference.cost ~ diagnosis + tkr_revision.1.age + 
    gender + tkr_revision.1.RCS.charlson.ra.omitted + IMD_2004_quintiles, 
                      family = Gamma(link="log"),
                      data=dat)
  
     cbind(data.frame(predicted.tkr_revision.cost=exp(predict(model,using.characteristics))),
           group=using.characteristics$group,
           time=using.characteristics$time)
      
}

regressAndPredict(using.data)
prob.costs<-plyr::rdply(sims, 
             regressAndPredict(using.data[sample(seq(length(using.data$diagnosis)),
                                          replace=TRUE) ,]),
             .id="sim")
prob.costs<-prob.costs %>% 
  mutate(sim=as.character(sim)) %>% 
  select(predicted.tkr_revision.cost,
           group,time,sim)

tkr.revision.costs<-rbind(
 deterministic.costs,
   prob.costs )
rm(deterministic.costs,prob.costs)

##### RUN MODEL FOR ALL SCENARIOS -----
          


# for every scenario
 # get tranisition matrix
 # get qol estimeates 
 # get cost estimates



all.scenarios.proportion_revised<-NULL
all.scenarios.markov_trace<-NULL


#using.scenario<-1
for(using.scenario in 1:length(scenarios[,1])) {
working.qol.improvement<-scenarios$qol.improvement[using.scenario]
working.revision.reduction<-scenarios$revision.reduction[using.scenario]
working.scenario.id<-scenarios$scenario.id[using.scenario]


### #####
# deterministic transitions ----
markov_trace<-NULL
proportion_revised<-NULL

for(using.mi in 1:n.mi) {
for(g in 1:max(characteristics$group)) {
  
working.group<-g

using.m.probs<-m.probs %>% 
  filter(group==working.group,
         mi==using.mi)

# reduce risk of revision given working scenario
using.m.probs<-using.m.probs %>%
  mutate(prob.unrevised_revision=prob.unrevised_revision*working.revision.reduction)


# transition probabilities
#from unrevised
p.U.Rn <- using.m.probs$prob.unrevised_revision
   # probability unrevised to revision
p.U.D <- using.m.probs$prob.unrevised_dead
   # probability unrevised to dead
p.U.U <- 1-p.U.Rn-p.U.D                                        
   # probability unrevised to unrevised

#from revision
p.Rn.D<-p.U.D
# probability revision to death 
p.Rn.Rd<-1-p.Rn.D
# probability revision to revised 

#from revised
p.Rd.D<-p.U.D
# probability revision to death 
p.Rd.Rd<-1-p.Rd.D
# probability revised to revised 


# set up Markov trace -
# states
state.names <- c("unrevised", "revision",
                 "revised","dead")         # state names
number.states<- length(state.names)                     # number of states
# starting distribution
start.distribution<-c(1,0,0,0)                      # all unrevised
# cycles
number.cycles <- max(m.probs$time[m.probs$group==working.group])


# set up Markov trace -
m.TR <- matrix(NA, nrow = number.cycles+1, ncol = number.states, 
               dimnames = list(1:(number.cycles+1), state.names))  
#head(m.TR)
m.TR[1,] <- start.distribution               # initialize Markov trace
#head(m.TR)

revisions<-NULL
# Transition matrix 
for(t in 1:(number.cycles)) 
{
# transition probabilities   
m.P.t <- rbind(c(p.U.U[t],   # unrevised-> unrevised
                 p.U.Rn[t],   # unrevised->revision 
                 0,          # unrevised->revised 
                 p.U.D[t]),  # unrevised->dead 

             c(0, # revision-> unrevised
               0, # revision-> revision
               p.Rn.Rd[t], # revision-> revised
               p.Rn.D[t]),  # revision-> dead   
             
             c(0, # revised-> unrevised
               0, # revised-> revision
               p.Rd.Rd[t], # revised-> revised
               p.Rd.D[t]),  # revised-> dead  
            
              c(0, # dead-> unrevised
                0, # dead-> revision
                0, # dead-> revised
                1)  # dead-> dead 
             )     


#transition
m.TR[t+1, ] <- m.TR[t,] %*% m.P.t

 revision.count<-p.U.Rn[t] * m.TR[t]
 revisions<-rbind(revisions, revision.count)
}
m.TR<-data.frame(m.TR)
m.TR$time<-0:number.cycles
m.TR$group<-working.group
m.TR<- m.TR %>% left_join(characteristics, by="group")
m.TR$mi<-using.mi
m.TR$sim<-"deterministic"

end.revised<-data.frame(prop.revised=sum(revisions),
                        group=working.group) # sum of transitions to revision
end.revised<- end.revised %>% left_join(characteristics, by="group")
end.revised$mi<-using.mi
end.revised$sim<-"deterministic"


markov_trace<-rbind(markov_trace,m.TR)
proportion_revised<-rbind(proportion_revised,end.revised)

}}

  



# probabilistic transitions ----
#markov_trace<-NULL
#proportion_revised<-NULL

for (using.sim in 1:sims) {
for(using.mi in 1:n.mi) {
for(g in 1:max(characteristics$group)) {
  
working.group<-g

using.m.probs<-m.probs %>% 
  filter(group==working.group,
         mi==using.mi, 
         sim==using.sim)


# reduce risk of revision given working scenario
using.m.probs<-using.m.probs %>%
  mutate(prob.unrevised_revision=prob.unrevised_revision*working.revision.reduction)



# transition probabilities
#from unrevised
p.U.Rn <- using.m.probs$prob.unrevised_revision
   # probability unrevised to revision
p.U.D <- using.m.probs$prob.unrevised_dead
   # probability unrevised to dead
p.U.U <- 1-p.U.Rn-p.U.D                                        
   # probability unrevised to unrevised

#from revision
p.Rn.D<-p.U.D
# probability revision to death 
p.Rn.Rd<-1-p.Rn.D
# probability revision to revised 

#from revised
p.Rd.D<-p.U.D
# probability revision to death 
p.Rd.Rd<-1-p.Rd.D
# probability revised to revised 


# set up Markov trace -
# states
state.names <- c("unrevised", "revision",
                 "revised","dead")         # state names
number.states<- length(state.names)                     # number of states
# starting distribution
start.distribution<-c(1,0,0,0)                      # all unrevised
# cycles
number.cycles <- max(m.probs$time[m.probs$group==working.group])


# set up Markov trace -
m.TR <- matrix(NA, nrow = number.cycles+1, ncol = number.states, 
               dimnames = list(1:(number.cycles+1), state.names))  
#head(m.TR)
m.TR[1,] <- start.distribution               # initialize Markov trace
#head(m.TR)

revisions<-NULL
# Transition matrix 
for(t in 1:(number.cycles)) {
# transition probabilities   
m.P.t <- rbind(c(p.U.U[t],   # unrevised-> unrevised
                 p.U.Rn[t],   # unrevised->revision 
                 0,          # unrevised->revised 
                 p.U.D[t]),  # unrevised->dead 

             c(0, # revision-> unrevised
               0, # revision-> revision
               p.Rn.Rd[t], # revision-> revised
               p.Rn.D[t]),  # revision-> dead   
             
             c(0, # revised-> unrevised
               0, # revised-> revision
               p.Rd.Rd[t], # revised-> revised
               p.Rd.D[t]),  # revised-> dead  
            
              c(0, # dead-> unrevised
                0, # dead-> revision
                0, # dead-> revised
                1)  # dead-> dead 
             )     


#transition
m.TR[t+1, ] <- m.TR[t,] %*% m.P.t

 revision.count<-p.U.Rn[t] * m.TR[t]
 revisions<-rbind(revisions, revision.count)
}
m.TR<-data.frame(m.TR)
m.TR$time<-0:number.cycles
m.TR$group<-working.group
m.TR<- m.TR %>% left_join(characteristics, by="group")
m.TR$mi<-using.mi
m.TR$sim<-using.sim


end.revised<-data.frame(prop.revised=sum(revisions),
                        group=working.group) # sum of transitions to revision
end.revised<- end.revised %>% left_join(characteristics, by="group")
end.revised$mi<-using.mi
end.revised$sim<-using.sim


markov_trace<-rbind(markov_trace,m.TR)
proportion_revised<-rbind(proportion_revised,end.revised)

}}}

  



# drop time zero ----
# time zero is just the starting state
# time 1 is the first year in the model
markov_trace<-markov_trace %>% 
  filter(time!=0)
# add qol to markov trace -----
a<-markov_trace %>% 
  select(time, group, mi, sim,q1_eq5d)

a<-a %>% left_join((
  qol %>% 
  select(group, mi, sim, pred.tkr.1.q2_eq5d_index)),
  by=c("group", "mi", "sim"))


# add improvements
a<-a %>% 
  mutate(unrevised.qol=
           pred.tkr.1.q2_eq5d_index*working.qol.improvement)



# year with primary surgery
a<-a %>% 
  mutate(unrevised.y1.qol=((((q1_eq5d+unrevised.qol)/2)
                        *6)+
                        (unrevised.qol*6))/12) 

# years onwards no surgery
a<-a %>% 
  mutate(unrevised.y2.plus.qol=unrevised.qol)

# revision and revised as 75% of 
# unevised with surgery and
# unrevised without surgery
# BUT without any improvement 
# (i.e. intervention only improves qol if unrevised)
a<-a %>% 
  mutate(revision.qol=(((((q1_eq5d+pred.tkr.1.q2_eq5d_index)/2)
                        *6)+
                        (pred.tkr.1.q2_eq5d_index*6))/12) * 0.75)
a<-a %>% 
  mutate(revised.qol=pred.tkr.1.q2_eq5d_index * 0.75)

# if revised is also improved by intervention
# a<-a %>% 
#   mutate(revision.qol=unrevised.y1.qol*0.75,
#          revised.qol=unrevised.y2.plus.qol*0.75)



a<-markov_trace %>% 
  select(time, group, mi,
         sim,
         unrevised, revision, revised, dead) %>% 
  left_join(a,
            by=c( "group","time", "mi", "sim"))


a<-a %>% 
  mutate(qol=(ifelse(time==1, 
                    unrevised*unrevised.y1.qol,
                    unrevised*unrevised.y2.plus.qol)) # unrevised 
             + 
             (revision*revision.qol)  # revision 
             +
             (revised*revised.qol)    # revised 
         )
a<-a %>% 
  select(group,time, mi,sim, unrevised.qol,qol)

markov_trace<-markov_trace %>% 
  left_join(a,
  by=c("group","time", "mi","sim"))

# discount qol -----
markov_trace<-markov_trace %>% 
   mutate(discounted.qol=(qol)/
                         ((1.035)^time))







# add primary cost: only in first year-----
#incurred by all
a<-markov_trace %>% 
  filter(time==1) %>% 
  left_join(tkr.primary.costs,
            by=c("group", "sim", "mi")) %>% 
  select(group, sim, mi,time, predicted.tkr.1.cost)
 # need to add in mi when the model has bmi and smoke

markov_trace<-markov_trace %>% 
  left_join(a, by=c("group", "sim", "mi","time"))
  
#zero if not first year
markov_trace$predicted.tkr.1.cost<-
  ifelse(is.na(markov_trace$predicted.tkr.1.cost),
         0, markov_trace$predicted.tkr.1.cost)


# add revision costs to trace -----
markov_trace<-markov_trace %>% 
  left_join(tkr.revision.costs,
            by=c("group", "sim", "time")) 
 # need to add in mi when the model has bmi and smoke


# only for those in the revision state
markov_trace<-markov_trace %>% 
              mutate(predicted.tkr_revision.cost=
                       predicted.tkr_revision.cost*revision)

markov_trace<-markov_trace %>% 
              mutate(cost=predicted.tkr.1.cost+
                          predicted.tkr_revision.cost)

markov_trace<-markov_trace %>% 
              select(-predicted.tkr.1.cost,
                     -predicted.tkr_revision.cost)



# discounted cost ----
markov_trace<-markov_trace %>% 
   mutate(discounted.cost=(cost)/
                         ((1.035)^time))




markov_trace$scenario.id<-working.scenario.id
markov_trace$qol.improvement<-working.qol.improvement
markov_trace$revision.reduction<-working.revision.reduction

proportion_revised$scenario.id<-working.scenario.id
proportion_revised$qol.improvement<-working.qol.improvement
proportion_revised$revision.reduction<-working.revision.reduction








all.scenarios.proportion_revised<-rbind(all.scenarios.proportion_revised,
      proportion_revised)
all.scenarios.markov_trace<-rbind(all.scenarios.markov_trace,
      markov_trace)



}



### TO KEEP -----
rm(list= ls()[!(ls() %in% c('characteristics',
                            'all.scenarios.markov_trace',
                            'all.scenarios.proportion_revised'))])






## SUMMARISE MODEL OUTPUT -----
# get proportion revised ------
deterministic.proportion_revised<-
  all.scenarios.proportion_revised %>% 
  filter(sim=="deterministic") %>% 
  group_by(group, scenario.id,qol.improvement,revision.reduction) %>% 
  summarise(mean=mean(prop.revised))   # mean across mi

probabilistic.proportion_revised<-
  all.scenarios.proportion_revised %>% 
  filter(sim!="deterministic") %>% 
  group_by(group, scenario.id,qol.improvement,revision.reduction) %>% 
  summarise(low.ci=quantile(prop.revised, probs = c(0.025), na.rm=TRUE),
            high.ci=quantile(prop.revised, probs = c(0.975), na.rm=TRUE))

proportion_revised<-deterministic.proportion_revised %>% 
  left_join(probabilistic.proportion_revised,
            by=c("group", "scenario.id",
                 "qol.improvement","revision.reduction"))


proportion_revised<-
  proportion_revised %>% 
  mutate(proportion_revised=
           paste0(sprintf("%.1f",mean*100), "%",
                  " (",
                  sprintf("%.1f",low.ci*100), "%",
                  " to ",
                  sprintf("%.1f",high.ci*100), "%",
                  ")")) %>% 
  select(-mean, -low.ci, -high.ci)
# get qol of unrevised -----
deterministic.qol_unrevised<-
all.scenarios.markov_trace %>% 
  filter(sim=="deterministic") %>% 
  filter(time==1) %>%
  group_by(group, scenario.id,qol.improvement,revision.reduction) %>% 
  summarise(mean=mean(unrevised.qol))

probabilistic.qol_unrevised<-
  all.scenarios.markov_trace %>% 
  filter(sim!="deterministic") %>% 
  filter(time==1) %>%
  group_by(group, scenario.id,qol.improvement,revision.reduction) %>% 
  summarise(low.ci=quantile(unrevised.qol, probs = c(0.025), na.rm=TRUE),
            high.ci=quantile(unrevised.qol, probs = c(0.975), na.rm=TRUE))

qol_unrevised<-deterministic.qol_unrevised %>% 
  left_join(probabilistic.qol_unrevised,
            by=c("group", "scenario.id",
                 "qol.improvement","revision.reduction"))

qol_unrevised<-
  qol_unrevised %>% 
  mutate(qol_unrevised=
           paste0(sprintf("%.2f",mean), 
                  " (",
                  sprintf("%.2f",low.ci), 
                  " to ",
                  sprintf("%.2f",high.ci), 
                  ")")) %>% 
  select(-mean, -low.ci, -high.ci)
# get QALYs and costs -----
# and net benefit and incremental net benefit
wtp<-20000


## deterministic
#calculate for each mi
deterministic.QALYs.costs<-all.scenarios.markov_trace %>% 
  filter(sim=="deterministic") %>% 
  group_by(group, scenario.id, qol.improvement,revision.reduction, 
           mi) %>% 
  summarise(QALYs=sum(discounted.qol),
            Costs=sum(discounted.cost))

deterministic.QALYs.costs$nmb<-wtp*deterministic.QALYs.costs$QALYs-
                        deterministic.QALYs.costs$Costs

# and average across mi
deterministic.QALYs.costs<-
  deterministic.QALYs.costs %>% 
  group_by(group, scenario.id, qol.improvement,revision.reduction) %>% 
  summarise(mean.QALYs=mean(QALYs),
            mean.Costs=mean(Costs),
            mean.nmb=mean(nmb))

# add incremental net benefit relative to scenario 1
deterministic.QALYs.costs<-deterministic.QALYs.costs %>% 
  arrange(qol.improvement, desc(revision.reduction))

deterministic.QALYs.costs$mean.inc.nmb<-0


#by group and scenario
for (working.screnario.id in 
     2:(max(deterministic.QALYs.costs$scenario.id))) {
for (working.group in 
     1:(max(deterministic.QALYs.costs$group))) {
  deterministic.QALYs.costs$mean.inc.nmb[
    deterministic.QALYs.costs$group==working.group &
    deterministic.QALYs.costs$scenario.id== working.screnario.id] <-
  deterministic.QALYs.costs$mean.nmb[
    deterministic.QALYs.costs$group==working.group &
    deterministic.QALYs.costs$scenario.id== working.screnario.id] -
  deterministic.QALYs.costs$mean.nmb[
     deterministic.QALYs.costs$group==working.group &
     deterministic.QALYs.costs$scenario.id== 1] 
}}
  




## probabilistic
probabilistic.QALYs.costs<-
  all.scenarios.markov_trace %>% 
  filter(sim!="deterministic") %>% 
  group_by(group, scenario.id, qol.improvement,revision.reduction, 
           mi, sim) %>% 
  summarise(QALYs=sum(discounted.qol),
            Costs=sum(discounted.cost))

probabilistic.QALYs.costs$nmb<-wtp*probabilistic.QALYs.costs$QALYs-
                        probabilistic.QALYs.costs$Costs

probabilistic.QALYs.costs<-probabilistic.QALYs.costs %>% 
  arrange(qol.improvement, desc(revision.reduction))
probabilistic.QALYs.costs$inc.nmb<-0
# get inc nmb by scenario, group, sim and mi 


for (working.screnario.id in 
     2:(max(probabilistic.QALYs.costs$scenario.id))) {
for (working.group in 
     1:(max(probabilistic.QALYs.costs$group))) {
for (using.mi in 
     1: max(as.numeric(probabilistic.QALYs.costs$mi))){  
for (using.sim in 
     1: max(as.numeric(probabilistic.QALYs.costs$sim))) {  
  
  probabilistic.QALYs.costs$inc.nmb[
    probabilistic.QALYs.costs$group==working.group &
    probabilistic.QALYs.costs$scenario.id== working.screnario.id & 
    probabilistic.QALYs.costs$mi == using.mi &
    probabilistic.QALYs.costs$sim ==using.sim] <-
 #working nmb
    probabilistic.QALYs.costs$nmb[
    probabilistic.QALYs.costs$group==working.group &
    probabilistic.QALYs.costs$scenario.id== working.screnario.id & 
    probabilistic.QALYs.costs$mi == using.mi &
    probabilistic.QALYs.costs$sim ==using.sim] -
    # nmb for base case scenario
   probabilistic.QALYs.costs$nmb[
    probabilistic.QALYs.costs$group==working.group &
    probabilistic.QALYs.costs$scenario.id== 1 & 
    probabilistic.QALYs.costs$mi == using.mi &
    probabilistic.QALYs.costs$sim ==using.sim]
 }}}}


probabilistic.QALYs.costs<-
  probabilistic.QALYs.costs %>% 
    group_by(group, scenario.id, 
             qol.improvement,revision.reduction) %>% 
  summarise(low.ci.QALYs=quantile(QALYs,
                      probs = c(0.025), na.rm=TRUE),            
            high.ci.QALYs=quantile(QALYs,
                      probs = c(0.975), na.rm=TRUE), 
            low.ci.Costs=quantile(Costs,
                      probs = c(0.025), na.rm=TRUE),            
            high.ci.Costs=quantile(Costs,
                      probs = c(0.975), na.rm=TRUE),
            low.ci.nmb=quantile(nmb,
                      probs = c(0.025), na.rm=TRUE),            
            high.ci.nmb=quantile(nmb,
                      probs = c(0.975), na.rm=TRUE),
            low.ci.inc.nmb=quantile(inc.nmb,
                      probs = c(0.025), na.rm=TRUE),            
            high.ci.inc.nmb=quantile(inc.nmb,
                      probs = c(0.975), na.rm=TRUE))




QALYs.costs<-deterministic.QALYs.costs %>% 
  left_join(probabilistic.QALYs.costs,
            by=c("group", "scenario.id",
                 "qol.improvement","revision.reduction"))



QALYs.costs<-
  QALYs.costs %>% 
  mutate(QALYs=paste0(
           sprintf("%.1f",mean.QALYs), 
           " (",
           sprintf("%.1f",low.ci.QALYs),
           " to ",
           sprintf("%.1f",high.ci.QALYs),
           ")")) %>% 
   mutate(Costs=paste0(
           sprintf("%.0f",mean.Costs), 
           " (",
           sprintf("%.0f",low.ci.Costs),
           " to ",
           sprintf("%.0f",high.ci.Costs),
           ")")) %>% 
  mutate(threshold.price=paste0(
           sprintf("%.0f",mean.inc.nmb), 
           " (",
           sprintf("%.0f",low.ci.inc.nmb),
           " to ",
           sprintf("%.0f",high.ci.inc.nmb),
           ")")) %>% 
  mutate(threshold.price=ifelse(scenario.id==1, NA, 
                                threshold.price)) %>% 
    select(group, scenario.id, qol.improvement,revision.reduction, 
           mean.QALYs, low.ci.QALYs, high.ci.QALYs,
           QALYs, 
           mean.Costs,low.ci.Costs, high.ci.Costs,
           Costs, 
           mean.inc.nmb, low.ci.inc.nmb, high.ci.inc.nmb,
           threshold.price)









         
         
# combine   -----      

summary<-proportion_revised %>% 
  left_join(qol_unrevised, 
            by=c("group", "scenario.id",
                 "qol.improvement","revision.reduction"))
summary<-summary %>% 
  left_join(QALYs.costs, 
            by=c("group", "scenario.id",
                 "qol.improvement","revision.reduction"))


summary<-characteristics %>% 
   right_join(summary,
             by=c("group"))
  

# add mapped oks from eq5d ------    
#OKS 
oks.beta<- 0.0224412
oks.intercept<- -0.0404485

summary<-summary %>% 
   mutate(q1_eq5d_mapped_to_OKS=
     ifelse(((summary$q1_eq5d-oks.intercept)/oks.beta)<0, 0, 
     ifelse(((summary$q1_eq5d-oks.intercept)/oks.beta)>48, 48, 
              (summary$q1_eq5d-oks.intercept)/oks.beta)))

summary<-summary %>% 
   mutate(qol_unrevised_mapped_to_OKS=
     ifelse((((summary$q1_eq5d*summary$qol.improvement)-oks.intercept)/oks.beta)<0, 0, 
     ifelse((((summary$q1_eq5d*summary$qol.improvement)-oks.intercept)/oks.beta)>48, 48, 
            ((summary$q1_eq5d*summary$qol.improvement)-oks.intercept)/oks.beta)))
  

# save summary tables -----         
# write.csv(summary,
#           "tkr.threshold_prices.csv",
#             row.names = F, 
#           na="")

### threshold prices grid
# (alternative to heat map)
summary.varying.improvement<-QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction %in%  seq(0.5, 1, 0.1)) %>% 
  filter(qol.improvement %in%  seq(1, 1.05, 0.01)) %>% 
  mutate(revision.reduction=
  paste0((1-revision.reduction)*100,
        "%")) %>% 
  mutate(qol.improvement=
  paste0((qol.improvement-1)*100,
        "%")) %>% 
  select(revision.reduction, qol.improvement, threshold.price) %>% 
  tidyr::spread(qol.improvement, threshold.price)

summary.varying.improvement<-summary.varying.improvement %>% 
  rename(`qol.improvement: 0%`=`0%`)

# write.csv(summary.varying.improvement,
#           "tkr.threshold_prices.varying.improvement.csv",
#             row.names = F, 
#           na="")

# save output for plotting -----
# only need QALYs.costs

# add characteristics 
QALYs.costs<-QALYs.costs %>% 
  left_join(characteristics,
            by="group")

tkr.QALYs.costs<-QALYs.costs
rm(list= ls()[!(ls() %in% c('tkr.QALYs.costs'))])

#save.image("tkr.QALYs.costs.average.characteristics.RData")


# save output --------