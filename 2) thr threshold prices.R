rm(list=ls())

######### SET UP ######
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

library(parallel)

# Data ------
load("Y:/CPRD for analysis.RData")
# for average characteristics


#transition probabilities
# models (for each mi dataset) 
load("Y:/cprd hes/working data/thr.to.revision_model.RData")
load("Y:/cprd hes/working data/thr.to.death_model.RData")


#load("Y:/proms with MI.RData")
load(file = "Y:/PROMs/models/thr_prom_model.RData")
load(file = "Y:/PROMs/models/thr_prom_model.bstrap.RData")


# get predicted post-op qol 


#cost models
#thr
load("Y:/reference cost models/first mi model thr costs.RData")
# bootstrapped 
# (for first mi dataset only!)
#load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/thr.cost.model.bstrap.mi_1.RData")
load("Y:/reference cost models/thr.cost.model.bstrap.mi_1.RData")

#thr_revision
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/thr_revision.cost.model.RData")
# bstrapped
load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/costs and los/reference cost models/thr_revision.cost.model.bstrap.RData")
# save(thr.cost.model.bstrap.mi_1, 
#      file = "Y:/reference cost models/thr.cost.model.bstrap.mi_1.RData")




## uk life tables 
life.table<-read.csv("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/diagnosis replacement and mortality/uk life table.csv")
life.table<-life.table %>% 
  rename(Male=male,
         Female=female)

load(file="Y:/CPRD HES/working data/thr_cohort.average.characteristics.RData")
thr.average.characteristics.risks<-thr_cohort.average.characteristics

# average characteristics
Mode <- function(x, na.rm = T) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


thr.average.characteristics<-cprd_patients %>% 
  filter(thr.1==1) %>% 
  summarise(median_thr_age=round(median(thr.1.age, na.rm=T)),
            mode_thr_diagnosis=as.character(Mode(diagnosis_hip)),
            mode_thr_gender=as.character(Mode(gender)),
            mode_thr_IMD_2004_quintiles=as.character(Mode(IMD_2004_quintiles)),
            mode_thr_RCS=as.character(Mode(thr.1.RCS.charlson.full)),
            median_thr_BMI=round(median(thr.1.BMI, na.rm=T)),
            mode_thr_smoke=as.character(Mode(thr.1.smoke)),
            median_thr_q1_eq5d=median(thr.1.q1_eq5d_index, na.rm=T),
            median_thr_year=median(thr.1.year, na.rm=T))

median_thr_age<-thr.average.characteristics$median_thr_age
mode_thr_diagnosis<-thr.average.characteristics$mode_thr_diagnosis
mode_thr_gender<-thr.average.characteristics$mode_thr_gender
mode_thr_IMD_2004_quintiles<-thr.average.characteristics$mode_thr_IMD_2004_quintiles
mode_thr_RCS<-thr.average.characteristics$mode_thr_RCS
mode_thr_smoke<-thr.average.characteristics$mode_thr_smoke
median_thr_BMI<-thr.average.characteristics$median_thr_BMI
median_thr_q1_eq5d<-thr.average.characteristics$median_thr_q1_eq5d
median_thr_year<-thr.average.characteristics$median_thr_year





# n.sim ----
n.sim<-1000



# set up cluster -----
# Calculate the number of cores
cl <- makeCluster(detectCores() - 1)
#leaves a core available for other tasks....
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(rms))

#clusterExport(cl, "mi.models.thr_revision")
clusterExport(cl, "life.table")
#clusterExport(cl,"prom.m")
#clusterExport(cl,"prom.m.bstrap")

########## ------

# Functions ------ 
run.model<-function(){
  #profvis::profvis({
  # first, a few things for general set up
  newdata<<-characteristics
  # list- one for each profile
  newdata.list<<-split(newdata, seq(nrow(newdata)))
  
  min.age<<-min(newdata$age)
  times<<-seq(0, (100-min.age))
  # up to at least 100
  
  
  ## GET TPs  
  print("getting tps")
  start<-Sys.time()
  ## deterministic tps revision 
  
  # tp to revision
  
  # function to get deterministic surviaval curves
  # for revision for each profile
  get.det.surv.revision<-function(df){
    
    det.tps<-summary(mi.models.thr.to.revision[[1]], #1st MI
                     t = times, 
                     ci=F,
                     newdata = df)
    
    
    # as data.frame 
    det.tps<-purrr::flatten_dfr(det.tps)
    #group 
    det.tps$group<-df$group
    #add group info
    det.tps <-  det.tps %>% 
      left_join(newdata,by="group")
    # add sim info
    det.tps$sim<-"deterministic"
    
    # models only up to age 100
    # remove tps if above 100 
    # delete rows where age is above 100
    det.tps$working.age<-det.tps$age+det.tps$time
    det.tps<-det.tps %>% filter(working.age<=100)
    # return dataframe
    det.tps
  }
  
  
  # get deterministic survival curves 
  # for thr for each profile
  clusterExport(cl, "times")
  clusterExport(cl, "newdata")
  #  clusterExport(cl, "mi.models.thr_revision")
  
  det.tps.revision<-lapply(newdata.list, get.det.surv.revision)
  
  # function to get tps from survival curve
  # specify survival data
  # and name for probability
  # drop time zero (for which there is no tp)
  # drop survival curve too? (not doing so at the moment)
  get.tp.from.surv<-function(df, name) {
    df <-  df %>% 
      mutate(!!name:=1-(est/ lag(est))) %>% #survival curve
      filter(time>0) %>% # drop time zero
      select(-est) # drop survival curve
  }
  
  # get thr tps from survival curve
  det.tps.revision<-parLapply(cl, det.tps.revision, 
                              get.tp.from.surv,
                              "p.revision") # name for probability
  
  ## deterministic tps mortality
  # 10 years from model
  # revert to lifetable after that
  
  # 1) first 10 years
  det.tps.mortality.model<-summary(mi.models.thr.to.death[[1]],
                                   t = seq(0,10), # 1st 10 years
                                   ci=F,
                                   newdata = newdata)
  names(det.tps.mortality.model)<-1:length(det.tps.mortality.model)
  
  # get  tps from survival curve
  det.tps.mortality.model<-parLapply(cl,det.tps.mortality.model, 
                                     get.tp.from.surv,
                                     "p.mortality") # name for probability
  
  
  # add group info to each set of tps
  for (i in 1:length(det.tps.mortality.model)) {
    det.tps.mortality.model[[i]]<- det.tps.mortality.model[[i]] %>% 
      mutate(group=i)
  }
  
  # add patient characteristics to tp thr
  get.characteristics<- function(df) {
    df <-  df %>% 
      left_join(newdata, 
                by="group")}
  det.tps.mortality.model<-parLapply(cl,det.tps.mortality.model, 
                                     get.characteristics)
  
  # add sim info
  det.tps.mortality.model<-parLapply(cl,det.tps.mortality.model, 
                                     function(df){
                                       df$sim<-"deterministic"
                                       df
                                     })
  
  # models only up to age 100
  # remove tps if above 100 
  # delete rows where age is above 100
  #  clusterExport(cl, "life.table")
  det.tps.mortality.model<-parLapply(cl,det.tps.mortality.model, 
                                     function(df){
                                       df$working.age<-df$age+df$time
                                       df<-df %>% 
                                         filter(working.age<=100)
                                       df
                                     })
  
  
  # 2) subsequent mortality from lifetables
  # depends on start age and gender for each group
  det.tps.mortality.lifetable<-parLapply(cl,det.tps.mortality.model, 
                                         function(df){
                                           # add lifetables for subsequent mortality
                                           using.lifetable<-life.table %>% 
                                             select(age,
                                                    !!(as.character(
                                                      df$gender[1]))) %>% 
                                             rename(p.mortality=
                                                      !!(as.character(
                                                        df$gender[1])))
                                           
                                           using.lifetable<-using.lifetable %>% 
                                             filter(age>!!df$age[1]+10)
                                           
                                           # add to prevoious 10 years 
                                           using.lifetable<-using.lifetable %>% 
                                             mutate(time=age-!!df$age[1]) %>% 
                                             mutate(group=df$group[1]) %>% 
                                             select(time,  p.mortality, group)
                                           
                                           using.lifetable
                                         })
  
  
  # add patient characteristics to tp thr
  det.tps.mortality.lifetable<-parLapply(cl,det.tps.mortality.lifetable,
                                         get.characteristics)
  
  # add sim info
  det.tps.mortality.lifetable<-parLapply(cl,det.tps.mortality.lifetable,
                                         function(df){
                                           df$sim<-"deterministic"
                                           df
                                         })
  
  # add working age
  det.tps.mortality.lifetable<-parLapply(cl,det.tps.mortality.lifetable,
                                         function(df){
                                           df$working.age<-df$age+df$time
                                           df
                                         })
  
  #combine
  det.tps.mortality<-rbind(
    plyr::ldply(det.tps.mortality.model, data.frame, .id=NULL),
    plyr::ldply(det.tps.mortality.lifetable, data.frame, .id=NULL)) %>% 
    arrange(group, time) # so that data.frane is in right order
  
  det.tps.mortality<-det.tps.mortality %>% 
    split(det.tps.mortality, 
          f = det.tps.mortality$group) # to list (dataframe by group)
  
  #  rm(det.tps.mortality.model, det.tps.mortality.lifetable)     
  
  ## probabilistic tps revision 
  sims.revision <- normboot.flexsurvreg(
    mi.models.thr.to.revision[[1]], 
    B=n.sim, 
    newdata = newdata)
  
  prob.revision<-NULL
  
  #mi.models.thr_revision[[1]]$call
  #spline
  if(length(newdata$group)==1){
    for (using.time in 0:tail(times,1)){
      for(using.group in 1:max(newdata$group)) {
        prob.revision<-rbind(prob.revision,
                             data.frame(time=using.time,
                                        surv=psurvspline(using.time,  
                                                         gamma = sims.revision, 
                                                         knots = mi.models.thr.to.revision[[1]]$knots,
                                                         lower.tail = FALSE),
                                        sim=1:n.sim,
                                        group=using.group))}}}
  
  
  if(length(newdata$group)>1){
    for (using.time in 0:tail(times,1)){
      for(using.group in 1:max(newdata$group)) {
        prob.revision<-rbind(prob.revision,
                             data.frame(time=using.time,
                                        surv=psurvspline(using.time,  
                                                         gamma = sims.revision[[using.group]], 
                                                         knots = mi.models.thr.to.revision[[1]]$knots,
                                                         lower.tail = FALSE),
                                        sim=1:n.sim,
                                        group=using.group))}}}
  
  
  #to transition probabilities
  prob.revision<-prob.revision %>% 
    group_by(sim, group) %>% 
    arrange(sim, group, time) %>% 
    mutate(est =  1-(surv/lag(surv))) %>% 
    filter(time>0) %>% # drop time zero
    select(-surv) %>% 
    ungroup()
  
  
  
  prob.revision<-split(prob.revision, 
                       prob.revision[,c('sim','group')])
  # to list (dataframe by sim and group)
  
  # add patient characteristics to tp thr
  get.characteristics<- function(df) {
    df <-  df %>% 
      left_join(newdata, 
                by="group")}
  
  prob.revision<-parLapply(cl,prob.revision, 
                           get.characteristics)
  
  names(prob.revision)
  
  # models only up to age 100
  # remove tps if above 100 
  # delete rows where age is above 100
  prob.revision<-parLapply(cl,prob.revision, 
                           function(df){
                             df$working.age<-df$age+df$time
                             df<-df %>% filter(working.age<=100)
                             df
                           })
  
  
  
  ## probabilistic tps mortality 
  #mortality
  # 1) 10 years from model
  sims.death <- normboot.flexsurvreg(
    mi.models.thr.to.death[[1]], 
    B=n.sim, 
    newdata = newdata)
  
  prob.death<-NULL
  #mi.models.thr.to.death[[1]]$call
  #spline
  if(length(newdata$group)==1){
    for (using.time in 0:10){
      for(using.group in 1:max(newdata$group)) {
        prob.death<-rbind(prob.death,
                          data.frame(time=using.time,
                                     surv=psurvspline(using.time,  
                                                      gamma = sims.death, 
                                                      knots = mi.models.thr.to.death[[1]]$knots,
                                                      lower.tail = FALSE),
                                     sim=1:n.sim,
                                     group=using.group))}}}
  
  if(length(newdata$group)>1){
    for (using.time in 0:10){
      for(using.group in 1:max(newdata$group)) {
        prob.death<-rbind(prob.death,
                          data.frame(time=using.time,
                                     surv=psurvspline(using.time,  
                                                      gamma = sims.death[[using.group]], 
                                                      knots = mi.models.thr.to.death[[1]]$knots,
                                                      lower.tail = FALSE),
                                     sim=1:n.sim,
                                     group=using.group))}}}
  
  
  # get tps from survival curve
  #to transition probabilities
  prob.death<-prob.death %>% 
    group_by(sim, group) %>% 
    arrange(sim, group, time) %>% 
    mutate(est =  1-(surv/lag(surv))) %>% 
    filter(time>0) %>% # drop time zero
    select(-surv) %>% 
    ungroup()
  
  prob.death<-split(prob.death, 
                    prob.death[,c('sim','group')])
  # to list (dataframe by sim and group)
  
  
  # add patient characteristics
  prob.death<-parLapply(cl,prob.death, 
                        get.characteristics)
  
  # models only up to age 100
  # remove tps if above 100 
  # delete rows where age is above 100
  prob.death<-parLapply(cl,prob.death, 
                        function(df){
                          df$working.age<-df$age+df$time
                          df<-df %>% 
                            filter(working.age<=100)
                          df
                        })
  
  
  # 2) subsequent mortality from lifetables
  # depends on start age and gender for each group
  # no uncertainty
  
  prob.tps.death.lifetable<-parLapply(cl,prob.death, 
                                      function(df){
                                        # add lifetables for subsequent mortality
                                        using.lifetable<-life.table %>% 
                                          select(age,
                                                 !!(as.character(
                                                   df$gender[1]))) %>% 
                                          rename(p.mortality=
                                                   !!(as.character(
                                                     df$gender[1])))
                                        
                                        using.lifetable<-using.lifetable %>% 
                                          filter(age>!!df$age[1]+10)
                                        
                                        # add to prevoious 10 years 
                                        using.lifetable<-using.lifetable %>% 
                                          mutate(time=age-!!df$age[1]) %>% 
                                          mutate(group=df$group[1]) %>% 
                                          select(time, p.mortality, group)
                                        
                                        using.lifetable
                                      })
  
  
  # add patient characteristics to tp thr
  prob.tps.death.lifetable<-parLapply(cl,prob.tps.death.lifetable,
                                      get.characteristics)
  
  # add sim info
  # same as for corresponding first 10 years
  for (i in 1:length(prob.tps.death.lifetable)) {
    prob.tps.death.lifetable[[i]]$sim<- 
      prob.death[[i]]$sim[1]
  }
  
  # reorder to match
  prob.tps.death.lifetable<-
    parLapply(cl,prob.tps.death.lifetable, 
              function(df){
                df$working.age<-df$age+df$time
                df<-df %>% 
                  select(time,sim,group,         
                         p.mortality,age, gender,diagnosis, 
                         IMD, RCS,charlson, BMI, smoke, smoking_status, q1_eq5d,
                         year,age.centered, BMI.centered,working.age)
                df
              })
  
  
  prob.death<-rbind(
    plyr::ldply(prob.death, 
                data.frame,
                .id=NULL) %>% 
      rename(p.mortality=est),
    plyr::ldply(prob.tps.death.lifetable,
                data.frame,
                .id=NULL)) %>% 
    arrange( group, sim, time) # so that data.frane is in right order
  
  prob.death<-split(prob.death, 
                    prob.death[,c('sim','group')])
  
  
  prob.tps.mortality<-prob.death
  # rm(prob.death,prob.tps.death.lifetable)     
  
  # combine det and prob revision 
  
  names(det.tps.revision[[1]])
  names(prob.tps.mortality[[1]])
  
  for (i in 1:length(det.tps.revision)){
    det.tps.revision[[i]]<-det.tps.revision[[i]] %>%
      select(time, group, sim, p.revision) }
  
  for (i in 1:length(prob.revision)){
    prob.revision[[i]]<-prob.revision[[i]] %>%
      rename(p.revision=est) %>% 
      select(time, group, sim, p.revision) }
  
  tps.revision<-append(det.tps.revision, 
                       prob.revision)
  
  names(tps.revision)
  names(tps.revision)<-1:length(tps.revision)
  
  # combine det and prob mortality
  
  names(det.tps.mortality[[1]])
  names(prob.tps.mortality[[1]])
  
  
  for (i in 1:length(det.tps.mortality)){
    det.tps.mortality[[i]]<-det.tps.mortality[[i]] %>%
      select(time, group, sim, p.mortality) }
  
  for (i in 1:length(prob.tps.mortality)){
    prob.tps.mortality[[i]]<-prob.tps.mortality[[i]] %>%
      select(time, group, sim, p.mortality) }
  
  tps.mortality<-append(det.tps.mortality, 
                        prob.tps.mortality)
  
  
  names(tps.mortality)
  names(tps.mortality)<-1:length(tps.mortality)
  
  
  # combine revision and revision tps 
  
  tps<-vector("list", length(tps.revision))
  for (i in 1:length(tps.revision)){
    tps[[i]]<-tps.revision[[i]] %>% 
      left_join(tps.mortality[[i]], 
                by=c("time", "group", "sim"))
  }
  
  
  tps<-parLapply(cl,tps,
                 get.characteristics)
  
  #  tps<<-tps
  
  print(Sys.time()-start)
  
  
  
  print("Adding revision reductions")
  start<-Sys.time()
  
  revision.reductions<-seq(0, 1, by=0.025)
  
  # add tps with reductions
  tps.with.reductions<-NULL
  for(i in 1:length(revision.reductions)){
    using.revision.reduction<-revision.reductions[i]
    
    working.tps.with.reductions<-
      lapply(tps, 
             function(df){
               df$p.revision<-df$p.revision*using.revision.reduction
               df$revision.reduction<-using.revision.reduction
               df
             })
    
    tps.with.reductions<-append(tps.with.reductions,
                                working.tps.with.reductions)
    
    tps.with.reductions
  }
  tps<-tps.with.reductions
  # tps<<-tps.with.reductions
  
  print(Sys.time()-start)
  
  
  # GET MARKOV TRACE AND PROP REVISED
  print("Getting Markov trace and proportion revised")
  start<-Sys.time()
  
  # set up Markov trace 
  # states
  state.names <<- c("thr", 
                    "thr_revision",
                    "revised",
                    "death")         # state names
  number.states<<- length(state.names)                     
  # number of states
  # starting distribution
  start.distribution<<-c(1,0,0, 0)                      
  # from "at risk"
  
  #df<-tps[[1]]
  clusterExport(cl, "state.names")
  clusterExport(cl, "number.states")
  clusterExport(cl, "start.distribution")
  
  
  # df is transition probabilities
  #trans<-parLapply(cl,tps,function(df) {
  
  trans<-parLapply(cl, 
                   tps,
                   function(df) {
                     
                     # cycles
                     number.cycles <- 99-  min(df$age)
                     
                     
                     m.TR <- matrix(NA, 
                                    nrow = number.cycles+1, 
                                    ncol = number.states, 
                                    dimnames = list(1:(number.cycles+1), 
                                                    state.names))  
                     m.TR[1,] <- start.distribution   
                     
                     
                     thr_revision.count<-NULL
                     
                     
                     for(t in 1:number.cycles) {   
                       m.P.t <- rbind(c(1-df$p.revision[t]-
                                          df$p.mortality[t],   # unrevised-> unrevised
                                        df$p.revision[t],   # unrevised->revision 
                                        0,          # unrevised->revised 
                                        df$p.mortality[t]),  # unrevised->dead 
                                      
                                      c(0, # revision-> unrevised
                                        0, # revision-> revision
                                        1-df$p.mortality[t], # revision-> revised
                                        df$p.mortality[t]),  # revision-> dead   
                                      
                                      c(0, # revised-> unrevised
                                        0, # revised-> revision
                                        1-df$p.mortality[t], # revised-> revised
                                        df$p.mortality[t]),  # revised-> dead  
                                      
                                      c(0, # dead-> unrevised
                                        0, # dead-> revision
                                        0, # dead-> revised
                                        1)  # dead-> dead 
                       )   
                       
                       
                       
                       
                       #transition
                       m.TR[t+1, ] <- m.TR[t,] %*% m.P.t
                       
                       working.thr_revision.count<-df$p.revision[t] * m.TR[t,1]
                       thr_revision.count<-rbind(thr_revision.count, 
                                                 working.thr_revision.count)
                       
                       
                     }
                     m.TR<-as.data.frame(m.TR)
                     m.TR$time<-1:nrow(m.TR)
                     m.TR$sim<-df$sim[1]
                     m.TR$group<-df$group[1]
                     m.TR$revision.reduction<-df$revision.reduction[1]
                     
                     
                     end.thr_revision<-data.frame(prop.revised=sum(thr_revision.count)) # sum of transitions to revision
                     end.thr_revision$group<-df$group[1]
                     end.thr_revision$sim<-df$sim[1]
                     end.thr_revision$revision.reduction<-df$revision.reduction[1]
                     
                     
                     # markov.trace<-rbind(markov.trace,m.TR)
                     # proportion.thr_revision<-rbind(proportion.thr_revision,
                     #                                end.thr_revision)
                     # 
                     
                     output <- list("markov.trace" = m.TR, 
                                    "proportion.thr_revision" = end.thr_revision)
                     output
                   })
  
  length(trans)
  proportion.thr_revision<-NULL
  markov.trace<-NULL
  for (i in 1:length(trans)) {
    proportion.thr_revision[[i]]<-trans[[i]]$proportion.thr_revision
    markov.trace[[i]]<-trans[[i]]$markov.trace
  }
  # rm(trans)
  
  get.characteristics<- function(df) {
    df <-  df %>% 
      left_join(newdata, 
                by="group")}
  
  
  proportion.thr_revision<-parLapply(cl,proportion.thr_revision, 
                                     get.characteristics)
  markov.trace<-parLapply(cl,markov.trace, 
                          get.characteristics)
  
  proportion.thr_revision<-proportion.thr_revision
  markov.trace<-markov.trace
  
  
  print(Sys.time()-start)
  
  
  ## GET MARKOV TRACE WITH QOL IMPROVEMENTS
  print("Adding qol improvements")
  start<-Sys.time()
  
  qol.improvements<-seq(1,1.05, 0.005) 
  
  
  # add markov traces qol with improvements
  markov.trace.with.qol.improvements<-NULL
  for(i in 1:length(qol.improvements)){
    using.qol.improvements<-qol.improvements[i]
    
    working.markov.trace.with.qol.improvements<-
      lapply(markov.trace, 
             function(df){
               df$qol.improvement<-using.qol.improvements
               df
             })
    
    markov.trace.with.qol.improvements<-append(markov.trace.with.qol.improvements,
                                               working.markov.trace.with.qol.improvements)
    
  }
  markov.trace<-markov.trace.with.qol.improvements
  # rm(markov.trace.with.qol.improvements)
  #markov.trace<-markov.trace
  print(Sys.time()-start)
  
  
  #get.markov.trace.with.qol  
  print("Getting qol for model states")
  start<-Sys.time()
  
  
  # predict post op qol for each profile and sim
  pred.qol<-expand.grid(group=newdata$group,
                        sim=c("deterministic",
                              1:n.sim))
  pred.qol$sim<-as.character(pred.qol$sim)
  pred.qol<-split(pred.qol, seq(nrow(pred.qol)))
  
  get.characteristics<- function(df) {
    df <-  df %>% 
      left_join(newdata, 
                by="group")}
  
  pred.qol<-lapply(pred.qol, 
                   get.characteristics)
  
  pred.qol<-lapply(pred.qol,
                   function(df) {
                     
                     if(df$sim[1]=="deterministic"){
                       # post-op qol
                       df$pred.qol.unrevised<-  predict(prom.m, 
                                                        df[1,] %>% 
                                                          mutate(diagnosis=as.character(
                                                            ifelse(diagnosis=="h_ost", "OA", "RA"))) %>% 
                                                          rename(thr.1.age=age,
                                                                 thr.1.RCS.charlson.ra.omitted=RCS,
                                                                 IMD_2004_quintiles=IMD,
                                                                 thr.1.q1_eq5d_index=q1_eq5d,
                                                                 thr.1.BMI=BMI,
                                                                 thr.1.smoke=smoke))
                     }            
                     
                     if(df$sim[1]!="deterministic"){
                       # post-op qol
                       df$pred.qol.unrevised<-  predict(prom.m.bstrap[[as.numeric(df$sim)]], 
                                                        df[1,] %>% 
                                                          mutate(diagnosis=as.character(
                                                            ifelse(diagnosis=="h_ost", "OA", "RA"))) %>% 
                                                          rename(thr.1.age=age,
                                                                 thr.1.RCS.charlson.ra.omitted=RCS,
                                                                 IMD_2004_quintiles=IMD,
                                                                 thr.1.q1_eq5d_index=q1_eq5d,
                                                                 thr.1.BMI=BMI,
                                                                 thr.1.smoke=smoke))
                     }      
                     df })
  
  pred.qol<-plyr::ldply(pred.qol, 
                        data.frame, .id=NULL) %>% 
    select(group, sim, pred.qol.unrevised)
  
  # add predicted post op qol to markov trace
  # markov trace as dataframe
  # (no need to be a list any more- will just merge in predicted post-op qol)
  markov.trace<-plyr::ldply(markov.trace, 
                            data.frame, .id=NULL)
  
  markov.trace<-markov.trace %>%
    mutate(sim=as.character(sim)) %>%
    inner_join(pred.qol, by=c("group", "sim"))
  
  # incorporate relative improvement
  markov.trace$unrevised.qol<-markov.trace$pred.qol.unrevised*markov.trace$qol.improvement                      
  
  
  # add qol for states
  #thr
  # depends on if in first year (i.e. year with surgery) 
  # or after
  markov.trace$thr.qol<-NA
  
  markov.trace<-markov.trace %>% 
    #first year
    mutate(thr.qol=
             ifelse(time==1,
                    (((((q1_eq5d+unrevised.qol)/2)
                       *6)+
                        (unrevised.qol*6))/12) , 
                    # subsequent years
                    unrevised.qol
             ))
  
  # revision and revised as 75% of 
  # unevised with surgery and
  # unrevised without surgery
  # BUT without any improvement 
  # (i.e. intervention only improves qol if unrevised)
  
  # thr state always includes surgery
  
  markov.trace<-markov.trace %>% 
    mutate(thr_revision.qol=(((((q1_eq5d+pred.qol.unrevised)/2)
                               *6)+
                                (pred.qol.unrevised*6))/12)*0.75) 
  
  # if revised is also improved by intervention
  # could just use 
  # markov.trace_with_qol<-markov.trace_with_qol %>% 
  #          mutate(thr.qol*0.75)
  
  # revised qol is same throughout
  # i.e. never includes surgery
  markov.trace<-markov.trace %>% 
    mutate(revised.qol=pred.qol.unrevised*0.75)
  
  
  # get qol for distribution across states
  markov.trace<-markov.trace %>% 
    mutate(qol=(thr*thr.qol) # unrevised 
           + 
             (thr_revision*thr_revision.qol)  # revision 
           +
             (revised*revised.qol)    # revised 
    )
  
  # discount qol
  markov.trace<-markov.trace %>% 
    mutate(discounted.qol=(qol)/
             ((1.035)^(time-1))) #i.e. year 1 undiscounted
  
  markov.trace<-markov.trace
  
  
  # will also make tps and proportion.thr_revision
  # dataframes
  tps<-plyr::ldply(tps, 
                   data.frame, .id=NULL)
  # tps<<-tps
  proportion.thr_revision<-plyr::ldply(proportion.thr_revision, 
                                       data.frame, .id=NULL)
  #  proportion.thr_revision<<-proportion.thr_revision
  print(Sys.time()-start)
  
  
  
  #get.costs 
  print("Getting costs for model states")
  start<-Sys.time()
  
  # costs depend on
  # patient profile 
  # sim
  # and, for revision, time in model (e.g. time at which revision happened)
  # predict post op qol for each profile and sim
  
  
  
  # thr
  # only incurred at time=1
  # and incurred by all
  # so depends on profile and time
  pred.cost<-expand.grid(group=newdata$group,
                         sim=c("deterministic",
                               1:n.sim),
                         time=times[times>0])
  pred.cost$sim<-as.character(pred.cost$sim)
  #add characteristics
  pred.cost<-pred.cost %>% 
    left_join(newdata,by="group")
  
  
  
  #as list
  pred.cost<-split(pred.cost, seq(nrow(pred.cost)))
  
  #get costs for each
  #df<-pred.cost[[1]]
  pred.cost<-lapply(pred.cost,
                    function(df) {
                      #thr cost
                      df$thr.cost<-NA
                      if(df$sim[1]=="deterministic"){
                        df$thr.cost<-exp(predict(thr.cost.mi.model,
                                                 df %>% 
                                                   mutate(diagnosis=
                                                            ifelse(diagnosis=="h_ost",
                                                                   "OA", "RA"),
                                                          thr.1.age=age,
                                                          gender= gender,#"Male",
                                                          thr.1.RCS.charlson.ra.omitted= RCS ,#"0",
                                                          IMD_2004_quintiles= IMD ,#"2",
                                                          thr.1.BMI=BMI,#30,
                                                          thr.1.smoke= smoke  #"Ex"
                                                   )))}
                      
                      if(df$sim[1]!="deterministic"){
                        df$thr.cost<-
                          exp(predict(thr.cost.model.bstrap.mi_1[[as.numeric(df$sim)]],
                                      df %>% 
                                        mutate(diagnosis=
                                                 ifelse(diagnosis=="h_ost",
                                                        "OA", "RA"),
                                               thr.1.age=age,
                                               gender= gender,#"Male",
                                               thr.1.RCS.charlson.ra.omitted= RCS ,#"0",
                                               IMD_2004_quintiles= IMD ,#"2",
                                               thr.1.BMI=BMI,#30,
                                               thr.1.smoke= smoke  #"Ex"
                                        )))}
                      
                      
                      #thr_revision cost
                      df$working.age<-  df$age + (df$time-1)
                      
                      df$thr_revision.cost<-NA
                      if(df$sim=="deterministic"){
                        df$thr_revision.cost<-exp(predict(glm.thr_revision.reference_cost,
                                                          df %>% 
                                                            mutate(diagnosis=
                                                                     ifelse(diagnosis=="h_ost",
                                                                            "OA", "RA"),
                                                                   thr_revision.1.age=working.age,
                                                                   gender= gender,#"Male",
                                                                   thr_revision.1.RCS.charlson.ra.omitted= RCS ,#"0",
                                                                   IMD_2004_quintiles= IMD ,#"2",
                                                                   thr_revision.1.BMI=BMI,#30,
                                                                   thr_revision.1.smoke= smoke  #"Ex"
                                                            )))}
                      
                      
                      if(df$sim!="deterministic"){
                        df$thr_revision.cost<-exp(predict(thr_revision.cost.model.bstrap[[as.numeric(df$sim)]],
                                                          df %>% 
                                                            mutate(diagnosis=
                                                                     ifelse(diagnosis=="h_ost",
                                                                            "OA", "RA"),
                                                                   thr_revision.1.age=working.age, #dependent on time in model
                                                                   gender= gender,#"Male",
                                                                   thr_revision.1.RCS.charlson.ra.omitted= RCS ,#"0",
                                                                   IMD_2004_quintiles= IMD ,#"2",
                                                                   thr_revision.1.BMI=BMI,#30,
                                                                   thr_revision.1.smoke= smoke  #"Ex"
                                                            )))}
                      
                      
                      
                      
                      
                      
                      df
                    })
  # as data.frame
  pred.cost<-plyr::ldply(pred.cost, 
                         data.frame, .id=NULL) %>% 
    select(group, sim, time, thr.cost, thr_revision.cost)
  
  # merge in with markove trace 
  markov.trace<-markov.trace %>%
    mutate(sim=as.character(sim)) %>%
    left_join(pred.cost, by=c("group", "sim", "time"))
  
  
  # thr cost- incurred by all only at time 1
  head(markov.trace$thr.cost)
  markov.trace<-markov.trace %>% 
    mutate(thr.cost=ifelse(time==1, thr.cost, NA))
  
  
  # thr_revision cost 
  # multiply by number progressing to revision state
  markov.trace$thr_revision.cost<-markov.trace$thr_revision.cost*markov.trace$thr_revision
  
  
  
  # total cost
  # missing is zero
  markov.trace<-markov.trace %>% 
    mutate(cost=rowSums(cbind(thr.cost,thr_revision.cost), na.rm=TRUE)) 
  #discounted cost
  markov.trace<-markov.trace %>% 
    mutate(discounted.cost=(cost)/
             ((1.035)^(time-1)))
  
  # markov.trace<<-markov.trace
  print(Sys.time()-start)
  
  
  
  #get.QALYs.costs
  print("Getting summary QALYs and costs")
  start<-Sys.time()
  
  
  wtp<-20000
  # ce threshold
  
  
  get.characteristics<- function(df) {
    df <-  df %>% 
      left_join(newdata, 
                by="group")}
  
  
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
  # rm(det.inc.nmb)
  
  
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
  # rm(prob.inc.nmb)
  
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
    proportion.thr_revision %>%
    filter(sim=="deterministic") %>%
    group_by(group,
             revision.reduction) %>%
    summarise(mean.prop.revised=mean(prop.revised)) %>%
    left_join(
      proportion.thr_revision %>%
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
    filter(time==1) %>%
    group_by(group,
             qol.improvement) %>%
    summarise(mean.unrevised.qol=mean(unrevised.qol)) %>%
    left_join(
      markov.trace %>%
        filter(sim!="deterministic") %>% 
        filter(time==1) %>%
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
  
  # translate eq5d to oks/ohs
  coef.okhs<- 0.0222
  cons.okhs<-	-0.0697
  # Based on Dakin et al 2012 - OKS to EQ-5D calculator
  # paper: https://www.ncbi.nlm.nih.gov/pubmed/22555470
  # AND
  
  # paper Mapping the Oxford hip score onto the EQ-5D utility index
  # https://link.springer.com/article/10.1007%2Fs11136-012-0174-y
  
  
  
  QALYs.costs$okhs<-(QALYs.costs$mean.unrevised.qol-cons.okhs)/coef.okhs
  
  
  
  # add characteristics
  QALYs.costs<-QALYs.costs %>% 
    left_join(newdata,
              by="group")
  
  print(Sys.time()-start)
  
  
  
  
  tps<<-tps
  markov.trace<<-markov.trace
  proportion.thr_revision<<-proportion.thr_revision
  QALYs.costs<<-QALYs.costs
  
}

########## -----
# Average profile ----
# patient characteristics -----
characteristics<-
  rbind(
    expand.grid(age=median_thr_age,
                gender=mode_thr_gender,
                diagnosis=mode_thr_diagnosis,
                IMD=as.character(mode_thr_IMD_2004_quintiles),
                RCS=mode_thr_RCS,
                charlson=thr.average.characteristics.risks$mode_charlson,
                BMI=median_thr_BMI,
                smoke=mode_thr_smoke,
                smoking_status=thr.average.characteristics.risks$mode_smoking_status,
                q1_eq5d=median_thr_q1_eq5d,
                year=median_thr_year))

characteristics<-characteristics%>% 
  mutate(age.centered=age-median_thr_age,
         BMI.centered=BMI-median_thr_BMI)

#thr.average.characteristics.risks$mode_smoking_status

characteristics$group<-seq(1, length(characteristics$age))

characteristics$diagnosis<-as.character(characteristics$diagnosis)
characteristics$gender<-as.character(characteristics$gender)
characteristics$RCS<-as.character(characteristics$RCS)
characteristics$IMD<-as.character(characteristics$IMD)
characteristics$smoke<-as.character(characteristics$smoke)




# run model ----
run.model()


# save -----
save("tps",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.average.characteristics.tps.RData")
save("markov.trace",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.average.characteristics.markov.trace.RData")
save("QALYs.costs",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.average.characteristics.QALYs.costs.RData")

########## ------
# Age and gender ----
# patient characteristics -----
characteristics<-
  rbind(
    expand.grid(age=seq(50,80, by=5),
                gender=c("Male", "Female"),
                diagnosis=mode_thr_diagnosis,
                IMD=as.character(mode_thr_IMD_2004_quintiles),
                RCS=mode_thr_RCS,
                charlson=thr.average.characteristics.risks$mode_charlson,
                BMI=median_thr_BMI,
                smoke=mode_thr_smoke,
                smoking_status=thr.average.characteristics.risks$mode_smoking_status,
                q1_eq5d=median_thr_q1_eq5d,
                year=median_thr_year))

characteristics<-characteristics%>% 
  mutate(age.centered=age-median_thr_age,
         BMI.centered=BMI-median_thr_BMI)

#thr.average.characteristics.risks$mode_smoking_status

characteristics$group<-seq(1, length(characteristics$age))

characteristics$diagnosis<-as.character(characteristics$diagnosis)
characteristics$gender<-as.character(characteristics$gender)
characteristics$RCS<-as.character(characteristics$RCS)
characteristics$IMD<-as.character(characteristics$IMD)
characteristics$smoke<-as.character(characteristics$smoke)





# run model ----
run.model()


# save -----
save("tps",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.age.gender.tps.RData")
save("markov.trace",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.age.gender.markov.trace.RData")
save("QALYs.costs",
     file="C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.age.gender.QALYs.costs.RData")

########## -----

# stop cluster----
stopCluster(cl)
