rm(list=ls())

#### SET UP ######
# Packages ------
library(dplyr)
library(survHE)
library(ggplot2)
library(cowplot)
library(scales)
#library(mice)
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


#cost models
#tkr
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/mi models tkr costs.RData")
tkr.cost.mi.model<-each.tkr.mi.model
rm(each.tkr.mi.model)
# bootstrapped 
# (for first mi dataset only!)
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/tkr.cost.model.bstrap.mi_1.RData")


#tkr_revision
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/tkr_revision.cost.model.RData")

# bstrapped
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/tkr_revision.cost.model.bstrap.RData")




## uk life tables 
life.table<-read.csv("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/diagnosis replacement and mortality/uk life table.csv")
life.table<-life.table %>% 
  rename(Male=male,
         Female=female)


# average characteristics
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



# Functions -----


# get.tps
# returns transition probabilities
# deterministic and probabilistic
# for a given flexsurvspline model 
# (see OA lifetime risk code for extracting for 
# other types of models)
# for a given set of times
# for a given patient profile
# and, for a given number of simulations
get.tps <- function(model,
                    times,
                    patient.profile,
                    n.sim) {
#deterministic
det<-summary(model,
        t = times, 
        ci=F,
        newdata = patient.profile)
det<-plyr::ldply(det, data.frame, .id = NULL)

det<-det %>% 
  select(time,est) %>% 
  rename(surv=est)

det<-det %>% 
  mutate(tp =1-(surv/ lag(surv)))

det<-det %>% 
  mutate(sim="deterministic")

#probabilistic
sims <- normboot.flexsurvreg(model, 
            B=n.sim, 
            newdata = patient.profile)
prob<-NULL
# model is flexsurv
for (using.time in 0:tail(times,1)){
prob<-rbind(prob,
           data.frame(time=using.time,
           surv=psurvspline(using.time,  
                  gamma = sims, 
                  knots = model$knots,
                  lower.tail = FALSE),
           sim=1:n.sim)) }
#to transition probabilities
prob<-prob %>% 
  group_by(sim) %>% 
  arrange(sim, time) %>% 
  mutate(tp =  1-(surv/lag(surv))) %>% 
  ungroup()

tps<-rbind(det, prob)
}
# a<-get.tps(model=mi.models.tkr_revision[[1]],
#         times=seq(0,50),
#         patient.profile=characteristics[1,],
#         n.sim=10)
# b<-get.tps(model=mi.models.tkr.to.death[[1]],
#         times=seq(0,10),
#         patient.profile=characteristics[1,],
#         n.sim=10)



# get.transitions
# returns markov trace for 
# patient profiles
# all revision reductions
# and all sims
get.transitions <- function(patient.profiles,
                            revision.reductions,
                            n.sim) {
# transition probabilities 
   # risk of rev - from surgery- extrpolated to age 100
  # risk of death - from surgery- first 10 years- then lifetable to 100
tps<-NULL
for(using.group in 1:max(patient.profiles$group)){
#using.group<-1

#1) get start age 
start.age<-patient.profiles$age[using.group]

# 2) risk of revision
tkr.to.tkr_revision.tps<-get.tps(model=mi.models.tkr_revision[[1]],
         times=seq(0, (100-start.age)),
         patient.profile=patient.profiles %>% 
                     filter(group==!!using.group),
         n.sim=n.sim)
tkr.to.tkr_revision.tps<-tkr.to.tkr_revision.tps %>% 
  mutate(group=!!using.group)

#3) risk of death 
# first ten years based on the survival model
tkr.to.death.tps<-
      get.tps(
           model = mi.models.tkr.to.death[[1]],
           times=seq(0, 10),
           patient.profile=patient.profiles %>% 
                     filter(group==!!using.group),
           n.sim=n.sim)
tkr.to.death.tps<-tkr.to.death.tps %>% 
  mutate(group=!!using.group)
# subsequent years will be based on lifetables
# add lifetables for subsequent mortality
using.lifetable<-life.table %>% 
  select(age,
         !!(as.character(
           patient.profiles$gender[using.group]))) %>% 
  rename(tp=
           !!(as.character(
           patient.profiles$gender[using.group])))

using.lifetable<-using.lifetable %>% 
    filter(age>!!start.age+10)

# add to prevoious 10 years 
# duplicate for deterministic and 
# each simulation
using.lifetable<-using.lifetable %>% 
  mutate(time=age-!!start.age)

using.lifetable<-using.lifetable %>% 
  mutate(surv=NA) %>% 
  select(time, surv, tp)

using.lifetable<- 
  cbind(
       using.lifetable[rep(1:nrow(using.lifetable),
                       times = n.sim+1),],
   rbind(data.frame(sim=rep("deterministic", 
                         length(using.lifetable[,1]))),
      data.frame(sim=as.character(rep(1:n.sim,
                         length(using.lifetable[,1])))) %>% 
      arrange(sim)))

using.lifetable$group<-using.group

tkr.to.death.tps<-
  rbind(tkr.to.death.tps,
    using.lifetable) %>%  
  arrange(sim, time)

working.tps<-
  rbind(
  tkr.to.tkr_revision.tps %>% mutate(transition="tkr.to.tkr_revision"),
  tkr.to.death.tps %>% mutate(transition="tkr.to.death"))

#working.tps$mi<-1
tps<-rbind(tps, #previous
          working.tps)
}


# transitions
# for all relative improvements
markov.trace<-NULL
proportion.tkr_revision<-NULL

for(using.group in 1:max(patient.profiles$group)){
for(using.revision.reduction in 1:length(revision.reductions)){
using.revision.reduction<-revision.reductions[using.revision.reduction]
# get start age 
start.age<-patient.profiles$age[using.group]
for(using.sim in 1:(n.sim+1)){
  
if(using.sim==1){
working.sim<-"deterministic"}

if(using.sim>1){
working.sim<-using.sim-1}
  
  
# transition probabilities
#from tkr
p.tkr.tkr_revision <- tps %>% 
              filter(sim==working.sim,
                     group==using.group) %>% 
              filter(time>0) %>% 
              filter(transition=="tkr.to.tkr_revision") %>% 
              select(tp)
   # probability tkr to tkr_revision
# add revision reduction
p.tkr.tkr_revision<-p.tkr.tkr_revision*using.revision.reduction

p.tkr.death <- tps %>% 
              filter(sim==working.sim,
                     group==using.group) %>% 
              filter(time>0) %>% 
              filter(transition=="tkr.to.death") %>% 
              select(tp)
   # probability tkr to death
p.tkr.tkr <- 1-p.tkr.tkr_revision-p.tkr.death                                        
   # probability tkr to tkr


#from revision
p.tkr_revision.death<-p.tkr.death
# probability revision to death 
p.tkr_revision.revised<-1-p.tkr.death
# probability revision to revised 

#from revised
p.revised.death<-p.tkr.death
# probability revision to death 
p.revised.revised<-1-p.revised.death
# probability revised to revised 



  
# set up Markov trace 
# states
state.names <- c("tkr", 
                 "tkr_revision",
                 "revised",
                 "death")         # state names
number.states<- length(state.names)                     
    # number of states
# starting distribution
start.distribution<-c(1,0,0,0)                      
   # from diagnosis
# cycles
number.cycles <- 100-start.age

m.TR <- matrix(NA, 
               nrow = number.cycles+1, 
               ncol = number.states, 
               dimnames = list(1:(number.cycles+1), 
                               state.names))  
#head(m.TR)
m.TR[1,] <- start.distribution  
 # initialize Markov trace
tkr_revision.count<-NULL



# Transition matrix 
for(t in 1:(number.cycles)) {
# transition probabilities   
 m.P.t <- rbind(c(p.tkr.tkr[t,],   # unrevised-> unrevised
                 p.tkr.tkr_revision[t,],   # unrevised->revision 
                 0,          # unrevised->revised 
                 p.tkr.death[t,]),  # unrevised->dead 

             c(0, # revision-> unrevised
               0, # revision-> revision
               p.tkr_revision.revised[t,], # revision-> revised
               p.tkr_revision.death[t,]),  # revision-> dead   
             
             c(0, # revised-> unrevised
               0, # revised-> revision
               p.revised.revised[t,], # revised-> revised
               p.revised.death[t,]),  # revised-> dead  
            
              c(0, # dead-> unrevised
                0, # dead-> revision
                0, # dead-> revised
                1)  # dead-> dead 
             )   

#transition
m.TR[t+1, ] <- m.TR[t,] %*% m.P.t

working.tkr_revision.count<-p.tkr.tkr_revision[t,] * m.TR[t,1]
tkr_revision.count<-rbind(tkr_revision.count, 
                          working.tkr_revision.count)

}

m.TR<-data.frame(m.TR)
m.TR$time<-0:number.cycles
m.TR$group<-using.group
m.TR<- m.TR %>% 
  left_join(patient.profiles, by="group")
#m.TR$mi<-using.mi
m.TR$sim<-working.sim
m.TR$revision.reduction<-using.revision.reduction


end.tkr_revision<-data.frame(prop.revised=sum(tkr_revision.count),
                        group=using.group) # sum of transitions to revision
end.tkr_revision<- end.tkr_revision %>% 
     left_join(patient.profiles, by="group")
#end.tkr_revision$mi<-using.mi
end.tkr_revision$sim<-working.sim
end.tkr_revision$revision.reduction<-using.revision.reduction

markov.trace<-rbind(markov.trace,m.TR)
proportion.tkr_revision<-rbind(proportion.tkr_revision,
                               end.tkr_revision)

}}
}

#add id for each markov trace
markov.trace<-markov.trace %>% 
  mutate(mtr.id=
  as.character(markov.trace %>% 
  group_indices(group, sim, revision.reduction))) %>% 
  arrange(group, sim, revision.reduction)


tps<<-tps
markov.trace<<-markov.trace
proportion.tkr_revision<<-proportion.tkr_revision}

#patient.profiles<-characteristics
#n.sim<-10
#revision.reductions<-seq(0, 1, by=0.025) 
     # relative reduction in risk of revision
     # 0: no reduction

#get.costs
# primary costs for all
# revision costs for proportion who get one
get.costs <- function() {

  markov.trace$working.age<-
  markov.trace$age +
  markov.trace$time

n.sim<-as.numeric(markov.trace %>% 
  filter(sim!="deterministic") %>% 
  summarise(max(as.numeric(sim))))

## tkr ------
# depends on
# sim 
# and patient profile

#n.sim<-2
#using.sim<-2
tkr_costs<-NULL
for(using.sim in 1:(n.sim+1)){
  
if(using.sim==1){
working.sim<-"deterministic"}
if(using.sim>1){
working.sim<-using.sim-1}

a<-markov.trace %>% 
  filter(sim==working.sim) %>% 
  filter(time=="0")  

if(working.sim=="deterministic"){
a$tkr.cost<-
  exp(predict(tkr.cost.mi.model[[1]],
            a %>% 
            mutate(diagnosis=
                     ifelse(diagnosis=="k_ost",
                            "OA", "RA"),
        tkr.1.age=working.age,
        gender= gender,#"Male",
        tkr.1.RCS.charlson.ra.omitted= RCS ,#"0",
        IMD_2004_quintiles= IMD ,#"2",
        tkr.1.BMI=BMI,#30,
        tkr.1.smoke= smoke  #"Ex"
            )))}

if(working.sim!="deterministic"){
a$tkr.cost<-
  exp(predict(tkr.cost.model.bstrap.mi_1[[working.sim]],
            a %>% 
            mutate(diagnosis=
                     ifelse(diagnosis=="k_ost",
                            "OA", "RA"),
        tkr.1.age=working.age,
        gender= gender,#"Male",
        tkr.1.RCS.charlson.ra.omitted= RCS ,#"0",
        IMD_2004_quintiles= IMD ,#"2",
        tkr.1.BMI=BMI,#30,
        tkr.1.smoke= smoke  #"Ex"
            )))
  }

#for entire cohort at start
a$tkr.cost<-a$tkr*a$tkr.cost
# no need to discount

tkr_costs<-rbind(tkr_costs,a)
}

markov.trace<-markov.trace %>% 
  left_join(tkr_costs %>% 
            select("mtr.id","time",
                   "tkr.cost"),
            by=c("mtr.id","time"))
rm(tkr_costs)

markov.trace<-markov.trace %>% 
  mutate(tkr.cost=ifelse(is.na(tkr.cost), 0, tkr.cost))




## tkr_revision ------
# depends on
# sim 
# and patient profile

#n.sim<-2
#using.sim<-2
tkr_revision_costs<-NULL

for(using.sim in 1:(n.sim+1)){
  
if(using.sim==1){
working.sim<-"deterministic"}
if(using.sim>1){
working.sim<-using.sim-1}

a<-markov.trace %>% 
  filter(sim==working.sim)

if(working.sim=="deterministic"){
a$tkr_revision.cost<-
  exp(predict(glm.tkr_revision.reference_cost,
            a %>% 
            mutate(diagnosis=
                     ifelse(diagnosis=="k_ost",
                            "OA", "RA"),
        tkr_revision.1.age=working.age,
        gender= gender,#"Male",
        tkr_revision.1.RCS.charlson.ra.omitted= RCS ,#"0",
        IMD_2004_quintiles= IMD ,#"2",
        tkr_revision.1.BMI=BMI,#30,
        tkr_revision.1.smoke= smoke  #"Ex"
            )))}

if(working.sim!="deterministic"){
a$tkr_revision.cost<-
  exp(predict(tkr_revision.cost.model.bstrap.mi_1[[working.sim]],
            a %>% 
            mutate(diagnosis=
                     ifelse(diagnosis=="k_ost",
                            "OA", "RA"),
        tkr_revision.1.age=working.age,
        gender= gender,#"Male",
        tkr_revision.1.RCS.charlson.ra.omitted= RCS ,#"0",
        IMD_2004_quintiles= IMD ,#"2",
        tkr_revision.1.BMI=BMI,#30,
        tkr_revision.1.smoke= smoke  #"Ex"
            )))
  }

#by proportion in tkr_revision
a$tkr_revision.cost<-a$tkr_revision*a$tkr_revision.cost

tkr_revision_costs<-rbind(tkr_revision_costs,a)
}

markov.trace<-markov.trace %>% 
  left_join(tkr_revision_costs %>% 
            select("mtr.id","time",
                   "tkr_revision.cost"),
            by=c("mtr.id","time"))
rm(tkr_revision_costs)


# total cost ------
markov.trace<-markov.trace %>% 
  mutate(cost=tkr.cost+tkr_revision.cost,
         discounted.cost=(cost)/
                         ((1.035)^time))

markov.trace<<-markov.trace
}


# get.qol
get.qol<-
  function(qol.improvements) {
  
# n.sim used for model
n.sim<-as.numeric(markov.trace %>% 
  filter(sim!="deterministic") %>% 
  summarise(max(as.numeric(sim))))  

# need transition matricies for each qol improvements
markov.trace_with_qol<-NULL
for (i in 1:length(qol.improvements)) {
using.markov.trace<-markov.trace %>% 
    mutate(qol.improvement=qol.improvements[i])
markov.trace_with_qol<-rbind(markov.trace_with_qol,
                             using.markov.trace)

}

# update matrix id
markov.trace_with_qol<-
  markov.trace_with_qol %>% 
  mutate(mtr.id=
  as.character(markov.trace_with_qol %>% 
  group_indices(group, sim, revision.reduction,
                qol.improvement))) %>% 
  arrange(group, sim, revision.reduction,
          qol.improvement)
#table(markov.trace_with_qol$mtr.id)

# get post-op qol
# depends on 
# patient characteristics
# and sim

# model
# expected
m<-ols(tkr.1.q2_eq5d_index ~ diagnosis + 
    rcs(tkr.1.age, 3) + 
      gender + tkr.1.RCS.charlson.ra.omitted + 
    IMD_2004_quintiles + rcs(tkr.1.q1_eq5d_index, 3) + tkr.1.BMI + 
    tkr.1.smoke,
    x=TRUE, y=TRUE,
    data=mice::complete(tkr.proms.imp, 1)) #1st mi dataset

# bootstrapped
set.seed(5)
m.bstrap<-list()
for(i in 1:n.sim) {
using.data<-mice::complete(tkr.proms.imp, 1)

sample<-using.data[sample(seq(length(using.data$diagnosis)),
                                          replace=TRUE) ,]

m.bstrap[[i]]<-ols(tkr.1.q2_eq5d_index ~ diagnosis + 
    rcs(tkr.1.age, 3) + 
      gender + tkr.1.RCS.charlson.ra.omitted + 
    IMD_2004_quintiles + rcs(tkr.1.q1_eq5d_index, 3) + tkr.1.BMI + 
    tkr.1.smoke,
    x=TRUE, y=TRUE,
    data=sample)
}


# predict
# given profile 
# and sim

qol.unrevised<-NULL
#using.sim<-"deterministic"
for(using.sim in 1:(n.sim+1)){

if(using.sim==1){
working.sim<-"deterministic"}
if(using.sim>1){
working.sim<-using.sim-1}

using.data<-markov.trace_with_qol %>% 
  filter(sim==!!working.sim) %>% 
  filter(time==0) %>% #just need for first- will merge with all
  mutate(diagnosis=as.character(
    ifelse(diagnosis=="k_ost", "OA", "RA"))) %>% 
  rename(tkr.1.age=age,
         tkr.1.RCS.charlson.ra.omitted=RCS,
         IMD_2004_quintiles=IMD,
         tkr.1.q1_eq5d_index=q1_eq5d,
         tkr.1.BMI=BMI,
         tkr.1.smoke=smoke)


if(working.sim=="deterministic"){
working.qol.unrevised<-
  data.frame(mtr.id=using.data$mtr.id,
             sim=as.character(working.sim),
             pred.qol.unrevised=predict(m, using.data),
             stringsAsFactors = F)
}

if(working.sim!="deterministic"){
working.qol.unrevised<-
  data.frame(mtr.id=using.data$mtr.id,
             sim=as.character(working.sim),
             pred.qol.unrevised=predict(m.bstrap[[working.sim]],
                                   using.data),
             stringsAsFactors = F)
  }

qol.unrevised<-rbind(qol.unrevised, working.qol.unrevised)
}

markov.trace_with_qol<-markov.trace_with_qol %>% 
  left_join(qol.unrevised,
            by=c("mtr.id", "sim"))
#same for each time

####
# qol unrevised
# incorporate relative markov.trace_with_qol
markov.trace_with_qol<-markov.trace_with_qol %>% 
  mutate(unrevised.qol=
           pred.qol.unrevised*
           qol.improvement)


# add qol for states
#tkr
# depends on if in first year (i.e. year with surgery) 
# or after
tkr.qol.t0<-markov.trace_with_qol %>% 
  filter(time==0) %>% 
  mutate(tkr.qol=((((q1_eq5d+unrevised.qol)/2)
                        *6)+
                        (unrevised.qol*6))/12) %>% 
  select(mtr.id, time, tkr.qol)

tkr.qol.t1.plus<-markov.trace_with_qol %>% 
  filter(time!=0) %>% 
  mutate(tkr.qol=unrevised.qol) %>% 
  select(mtr.id, time, tkr.qol)
tkr.qol<-rbind(tkr.qol.t0, tkr.qol.t1.plus)
rm(tkr.qol.t0, tkr.qol.t1.plus)

markov.trace_with_qol<-markov.trace_with_qol %>% 
  left_join(tkr.qol,
            by=c("mtr.id", "time"))




# revision and revised as 75% of 
# unevised with surgery and
# unrevised without surgery
# BUT without any improvement 
# (i.e. intervention only improves qol if unrevised)

# tkr state always includes surgery
tkr_revision.qol<-markov.trace_with_qol %>% 
#  filter(time==0) %>% 
  mutate(tkr_revision.qol=(((((q1_eq5d+unrevised.qol)/2)
                        *6)+
                        (unrevised.qol*6))/12)*0.75) %>% 
  select(mtr.id, time, tkr_revision.qol)

markov.trace_with_qol<-markov.trace_with_qol %>% 
  left_join(tkr_revision.qol,
            by=c("mtr.id", "time"))

# if revised is also improved by intervention
# could just use 
# markov.trace_with_qol<-markov.trace_with_qol %>% 
#          mutate(tkr.qol*0.75)

# revised qol is same throughout
# i.e. never includes surgery
markov.trace_with_qol<-markov.trace_with_qol %>% 
  mutate(revised.qol=pred.qol.unrevised * 0.75)



# get qol for distribution across states
markov.trace_with_qol<-markov.trace_with_qol %>% 
  mutate(qol=(tkr*tkr.qol) # unrevised 
             + 
             (tkr_revision*tkr_revision.qol)  # revision 
             +
             (revised*revised.qol)    # revised 
         )

# discount qol
markov.trace_with_qol<-markov.trace_with_qol %>% 
   mutate(discounted.qol=(qol)/
                         ((1.035)^time))


markov.trace<<-markov.trace_with_qol }

#get.QALYs.costs
get.QALYs.costs <- function(wtp){


# by group

## deterministic
#calculate for each mi
deterministic.QALYs.costs<-markov.trace %>% 
  filter(sim=="deterministic") %>% 
  group_by(group, 
           revision.reduction,
           qol.improvement) %>% 
  summarise(mean.QALYs=sum(discounted.qol),
            mean.Costs=sum(discounted.cost))

deterministic.QALYs.costs$nmb<-wtp*
               deterministic.QALYs.costs$mean.QALYs-
               deterministic.QALYs.costs$mean.Costs

# add incremental net benefit relative to scenario 1
# by group
det.inc.nmb<-NULL
for (i in 1:max(deterministic.QALYs.costs$group)) {
  using.deterministic.QALYs.costs<-
  deterministic.QALYs.costs %>% 
    filter(group==!!i)
  
using.deterministic.QALYs.costs<-
  using.deterministic.QALYs.costs %>% 
  arrange(qol.improvement, desc(revision.reduction))

#deterministic.QALYs.costs$mean.inc.nmb<-0

using.deterministic.QALYs.costs$scenario.id<-
  1:nrow(using.deterministic.QALYs.costs)

# inc.nmb relative to scenario 1
using.deterministic.QALYs.costs$mean.inc.nmb<-
using.deterministic.QALYs.costs$nmb-
using.deterministic.QALYs.costs$nmb[
  using.deterministic.QALYs.costs$scenario.id==1]

det.inc.nmb<-rbind(det.inc.nmb,
                   using.deterministic.QALYs.costs)
}
deterministic.QALYs.costs<-det.inc.nmb
rm(det.inc.nmb)


## probabilistic
probabilistic.QALYs.costs<-
  markov.trace %>% 
  filter(sim!="deterministic") %>% 
  group_by(group, 
           qol.improvement,revision.reduction, 
           sim) %>% 
  summarise(QALYs=sum(discounted.qol),
            Costs=sum(discounted.cost))

probabilistic.QALYs.costs$nmb<-wtp*probabilistic.QALYs.costs$QALYs-
                        probabilistic.QALYs.costs$Costs

probabilistic.QALYs.costs<-probabilistic.QALYs.costs %>% 
  arrange(qol.improvement, desc(revision.reduction))

# get inc.nmb by simulation
prob.inc.nmb<-NULL
for (i in  1: max(as.numeric(probabilistic.QALYs.costs$sim))) {
for (j in 1:max(deterministic.QALYs.costs$group)) {

using.probabilistic.QALYs.costs<-  
  probabilistic.QALYs.costs %>% 
  filter(sim==!!i) %>% 
  filter(group==!!j)
       
using.probabilistic.QALYs.costs$scenario.id<-
  1:nrow(using.probabilistic.QALYs.costs)

# inc.nmb relative to scenario 1
using.probabilistic.QALYs.costs$inc.nmb<-
using.probabilistic.QALYs.costs$nmb-
using.probabilistic.QALYs.costs$nmb[
  using.probabilistic.QALYs.costs$scenario.id==1]
  
prob.inc.nmb<-
  rbind(prob.inc.nmb,
  using.probabilistic.QALYs.costs)
}}

probabilistic.QALYs.costs<-
  prob.inc.nmb
rm(prob.inc.nmb)

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
            by=c("group",
                 "scenario.id",
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
    select(group, 
           scenario.id, qol.improvement,revision.reduction, 
           mean.QALYs, low.ci.QALYs, high.ci.QALYs,
           QALYs, 
           mean.Costs,low.ci.Costs, high.ci.Costs,
           Costs, 
           mean.inc.nmb, low.ci.inc.nmb, high.ci.inc.nmb,
           threshold.price)



# add prop revised
prop_revised<-
  proportion.tkr_revision %>%
  filter(sim=="deterministic") %>%
  group_by(group,
           revision.reduction) %>%
  summarise(mean.prop.revised=mean(prop.revised)) %>%
  left_join(
  proportion.tkr_revision %>%
  filter(sim!="deterministic") %>%
  group_by(group,
           revision.reduction) %>%
  summarise(low.ci.prop.revised=quantile(prop.revised,
                      probs = c(0.025), na.rm=TRUE),
            high.ci.prop.revised=quantile(prop.revised,
                      probs = c(0.975), na.rm=TRUE)),
  by=c("group", "revision.reduction"))
    
QALYs.costs<-QALYs.costs %>% 
  left_join(prop_revised,
            by=c("group", "revision.reduction"))

# add estimated qol if unrevised
qol.improvement<-markov.trace %>%
  filter(sim=="deterministic") %>% 
  filter(time==0) %>%
  group_by(group,
           qol.improvement) %>%
  summarise(mean.unrevised.qol=mean(unrevised.qol)) %>%
  left_join(
markov.trace %>%
  filter(sim!="deterministic") %>% 
  filter(time==0) %>%
  group_by(group,
           qol.improvement) %>%
  summarise(low.ci.unrevised.qol=quantile(unrevised.qol,
                      probs = c(0.025), na.rm=TRUE),
            high.ci.unrevised.qol=quantile(unrevised.qol,
                      probs = c(0.975), na.rm=TRUE)),
  by=c("group", "qol.improvement"))



QALYs.costs<-QALYs.costs %>% 
  left_join(qol.improvement,
            by=c("group", "qol.improvement"))


# add characteristics
QALYs.costs<-QALYs.costs %>% 
  left_join(characteristics,
            by="group")






QALYs.costs<<-QALYs.costs}


# AVERAGE CHARACTERISTICS -----
# 1) patient profiles -----
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





# 2) get transitions -------
# with all relative risk reductions
get.transitions(patient.profiles=characteristics,
                revision.reductions=seq(0, 1, by=0.025),#seq(0, 1, by=0.025),
                 # relative reduction in risk of revision
                 # 0: no reduction
                n.sim=15)

a<- markov.trace %>%
  mutate(check=tkr+tkr_revision+revised+death)
table(a$check)
rm(a)

# 3) get costs -----
get.costs()

# 4) get qol -----
# with all relative improvements
get.qol(qol.improvements=seq(1,1.05, 0.0025))
# 5) summarise results -----
# and net benefit and incremental net benefit
get.QALYs.costs(wtp=20000)


    
# example plot
# qol improvements 
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1) %>% 
  mutate(qol.improvement.name=qol.improvement-1) %>% 
  ggplot(aes(x=qol.improvement.name))+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)
# 6) save -----
save("QALYs.costs",
 file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/tkr.average.characteristics.QALYs.costs.RData")

# AGE AND GENDER ------
# 1) patient profiles -----

 characteristics<-
  rbind(
  expand.grid(age=seq(50,80, by=5),
              gender=c("Male", "Female"),
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





# 2) get transitions -------
# with all relative risk reductions
get.transitions(patient.profiles=characteristics,
                revision.reductions=seq(0, 1, by=0.025),#seq(0, 1, by=0.025),
                 # relative reduction in risk of revision
                 # 0: no reduction
                n.sim=15)

# a<- markov.trace %>%
#   mutate(check=tkr+tkr_revision+revised+death)
# table(a$check)
# rm(a)

# 3) get costs -----
get.costs()

# 4) get qol -----
# with all relative improvements
get.qol(qol.improvements=seq(1,1.05, 0.0025))
# 5) summarise results -----
# and net benefit and incremental net benefit
get.QALYs.costs(wtp=20000)


# example plot
# qol improvements 
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1) %>% 
  mutate(qol.improvement.name=qol.improvement-1) %>% 
  ggplot(aes(x=qol.improvement.name, group=group))+
  geom_line(aes(y=mean.inc.nmb,
                colour=group))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb,
                  fill=group),
              alpha=0.1)

# 6) save -----
save("QALYs.costs",
file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/tkr.age.gender.QALYs.costs.RData")