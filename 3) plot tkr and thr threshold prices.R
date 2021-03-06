rm(list=ls())

#packages ----
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(ggrepel)

#1 for average characteristics -----
# data -----
rm(list=ls())

load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/tkr.average.characteristics.QALYs.costs.RData")
tkr.QALYs.costs<-QALYs.costs
rm(QALYs.costs)

load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.average.characteristics.QALYs.costs.RData")
thr.QALYs.costs<-QALYs.costs
rm(QALYs.costs)

QALYs.costs<-rbind(tkr.QALYs.costs %>% 
                   mutate(op="Knee replacement"), 
                   thr.QALYs.costs %>% 
                   mutate(op="Hip replacement")) %>% 
  mutate(op=factor(op,levels=c("Knee replacement", 
                               "Hip replacement")))



# specific estimates -----
#lifetime risk of revision
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
           qol.improvement==1 ) %>%
  mutate(est=paste0(sprintf("%.2f",mean.prop.revised*100),
                    "% (",
                    sprintf("%.2f",low.ci.prop.revised*100),
                    "% to ",
                    sprintf("%.2f",high.ci.prop.revised*100),
                    "%)")) %>% 
  select(op, est)

#qalys
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
           qol.improvement==1 ) %>%
  select(op,QALYs)


#qalys
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
           qol.improvement==1 ) %>%
  select(op,Costs)

#5% improvement in qol
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
           qol.improvement==1.05 ) %>%
  select(op,QALYs)
# threshold price 
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
           qol.improvement==1.05 ) %>%
  mutate(est=threshold.price) %>% 
  select(op, est)

# 50% reducion in revision
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1 ) %>%
  select(op,QALYs)
#cost
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1 ) %>%
  select(op,Costs)
# threshold price
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1 ) %>%
  mutate(est=threshold.price) %>% 
  select(op, est)


#5% improvement in qol, and
# 50% reducion in revision
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1.05 ) %>%
  select(op,QALYs)
#cost
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1.05 ) %>%
  select(op,Costs)
# threshold price
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1.05 ) %>%
  mutate(est=threshold.price) %>% 
  select(op, est)


# lifetime risk of revision estimate ------
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1 ) %>% 
  select(op, mean.prop.revised,
         low.ci.prop.revised,
         high.ci.prop.revised)
# heatmap -----
# shows the mean threshold price
# for combined improvements in revision risk and qol
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1) %>% 
  filter(qol.improvement==1.05) %>% 
  select(op,threshold.price)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5) %>% 
  filter(qol.improvement==1) %>% 
  select(op,threshold.price)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5) %>% 
  filter(qol.improvement==1.05) %>% 
  select(op,threshold.price)

QALYs.costs %>% 
  arrange(desc(revision.reduction),qol.improvement) %>% 
  ggplot()+
  geom_tile(aes(x=abs(revision.reduction-1),
                y=qol.improvement-1,
                fill=mean.inc.nmb)) +
   # scale_x_reverse()+
  xlab("Relative reduction in revision risk")+
  ylab("Relative improvement in quality of life if unrevised")+
  labs(fill = "Threshold price (\u00A3)")+
  scale_fill_gradient(low = "white", high = "black")+
  theme_bw(base_size = 18)+ 
  scale_x_continuous(labels=percent)+
  scale_y_continuous(labels=percent)+
  facet_grid(. ~ op)


#colour
QALYs.costs %>% 
  # filter(revision.reduction %in% seq(0,1, 0.25)) %>% # less granular
  # filter(qol.improvement %in% seq(1,1.05, 0.015)) %>% # less granular
  filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
  arrange(desc(revision.reduction),qol.improvement) %>% 
  ggplot()+
  geom_tile(aes(x=abs(revision.reduction-1),
                y=qol.improvement-1,
                fill=mean.inc.nmb)) +
  # scale_x_reverse()+
  xlab("Relative reduction in revision risk")+
  ylab("Relative improvement\nin quality of life if unrevised")+
  labs(fill = "Threshold price (\u00A3)")+
  #scale_fill_gradient(low = "white", high = "black")+ 
  scale_fill_gradientn(colours = 
                         c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))+
  theme_bw(base_size = 24)+ 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand=c(0,0))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand=c(0,0))+
  facet_grid(. ~ op)+theme(panel.spacing = unit(2, "lines"))


 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.heatmap.tiff",
          dpi=300,
              width = 14, height = 6)
 
 # with labels
 QALYs.costs %>% 
   # filter(revision.reduction %in% seq(0,1, 0.25)) %>% # less granular
   # filter(qol.improvement %in% seq(1,1.05, 0.015)) %>% # less granular
   filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
   mutate(lab=NA) %>% 
   mutate(lab=ifelse((revision.reduction==0.5 & qol.improvement==1) |
                       (revision.reduction==0.5 & qol.improvement==1.05) |
                       (revision.reduction==1 & qol.improvement==1) |
                       (revision.reduction==1 & qol.improvement==1.05),
                     paste0("Lifetime risk of\nrevision: ",
                            format(round(mean.prop.revised*100,1),nsmall = 1), 
                            "%",
                            "\n", 
                            paste0("EQ-5D-3L index: ",
                                   format(round(mean.unrevised.qol,2),nsmall = 1), 
                                   "\nOKHS: ",
                                   format(round(okhs,1),nsmall = 1))),
                     lab
   )) %>% 
   arrange(desc(revision.reduction),qol.improvement) %>% 
   ggplot(aes(x=abs(revision.reduction-1),
              y=qol.improvement-1,
              label=lab))+
   geom_tile(aes(fill=mean.inc.nmb)) +
   geom_point(colour=NA)+
   xlab("Relative reduction in revision risk")+
   ylab("Relative improvement\nin quality of life if unrevised")+
   labs(fill = "Threshold price (\u00A3)")+
   #scale_fill_gradient(low = "white", high = "black")+ 
   scale_fill_gradientn(colours = 
                          c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))+
   theme_bw(base_size = 24)+ 
   scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                      expand=c(0,0))+
   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                      expand=c(0,0))+
   facet_grid(. ~ op)+theme(panel.spacing = unit(2, "lines"))+
   geom_label_repel(segment.color=NA) 
 
 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.heatmap.label.tiff",
        dpi=300,
        width = 14, height = 6)
 
 
# improvement in qol assumption ------
 a<-rbind(
   QALYs.costs %>% 
     filter(revision.reduction==1 &
              qol.improvement==1),
   QALYs.costs %>% 
     filter(revision.reduction==1 &
              qol.improvement==1.05)) %>% 
   mutate(type=ifelse(qol.improvement==1,
                      "Current",
                      "Improved \nby 5%"))
 
 
 
 a<-rbind(
   a %>% 
     ungroup() %>% 
     select(op,type, age, gender, q1_eq5d) %>% 
     rename(mean.qol=q1_eq5d) %>% 
     mutate(low.ci.qol=mean.qol,
            high.ci.qol=mean.qol,
            time=0),
   
   a %>% 
     ungroup() %>% 
     select(op, type, age, gender, mean.unrevised.qol,
            low.ci.unrevised.qol,
            high.ci.unrevised.qol) %>% 
     rename(mean.qol=mean.unrevised.qol,
            low.ci.qol=low.ci.unrevised.qol,
            high.ci.qol=high.ci.unrevised.qol) %>% 
     mutate(time=6),
   
   a %>% 
     ungroup() %>% 
     select(op, type, age, gender, mean.unrevised.qol,
            low.ci.unrevised.qol,
            high.ci.unrevised.qol) %>% 
     rename(mean.qol=mean.unrevised.qol,
            low.ci.qol=low.ci.unrevised.qol,
            high.ci.qol=high.ci.unrevised.qol) %>% 
     mutate(time=12) 
 )
 
 
 a %>% 
   ggplot(aes(x=time, group=type))+
   facet_grid(. ~ op)+
   geom_line(aes(y=mean.qol, 
                 linetype=type))+  
   geom_ribbon(aes(ymin=low.ci.qol,
                   ymax=high.ci.qol),
               alpha=0.1)+
   scale_linetype_manual(name = "Post-operative \nquality of life",
                         values=c("solid", "longdash"))+
   ylab("EQ-5D-3L over year of\nprimary surgery (if unrevised)")+
   xlab("Time (months)")+
   scale_x_continuous(limits = c(0, 12),
                      breaks= c(0,3,6,9,12))+
   theme_bw(base_size = 24)
 
 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/CPRD HES lifetime risk of revision/plots/tkr and thr markov model/average_characteristics.qol.improvement.tiff",
        dpi=300,
        width = 12, height = 6)
 
 
 
 
 
 
 
 
 
 
 
# qol improvements -----
max.y<-as.numeric(QALYs.costs %>% 
   filter(revision.reduction==1) %>% 
   mutate(qol.improvement.name=qol.improvement-1) %>% 
   summarise(max.y=max(high.ci.inc.nmb)) %>% 
     ungroup() %>% 
     select(max.y))
 
 QALYs.costs %>% 
  filter(revision.reduction==1) %>% 
  mutate(qol.improvement.name=qol.improvement-1) %>% 
  ggplot(aes(x=qol.improvement.name))+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative improvement in quality of life if unrevised")+
   scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
   scale_y_continuous(limits=c(0, max.y))+
   theme_bw(base_size = 24)+ 
  theme(legend.title=element_blank())+
  facet_grid(. ~ op)

 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.qol.improvement.tiff",
          dpi=300,
              width = 12, height = 7)
 
 #with labels
 QALYs.costs %>% 
   filter(revision.reduction==1) %>% 
   mutate(qol.improvement.name=qol.improvement-1) %>% 
   mutate(lab=ifelse(qol.improvement==1 | qol.improvement==1.05,
                     paste0("EQ-5D-3L index: ",
                            format(round(mean.unrevised.qol,2),nsmall = 1), 
                            "\nOKHS: ",
                            format(round(okhs,1),nsmall = 1)),
                     NA)) %>% 
   ggplot(aes(x=qol.improvement.name, y=mean.inc.nmb,
              label=lab))+
   geom_line()+
   geom_point(colour=NA)+
   geom_ribbon(aes(ymin=low.ci.inc.nmb,
                   ymax=high.ci.inc.nmb),
               alpha=0.1)+
   scale_linetype_manual(values=c("solid", "longdash"))+
   ylab("Threshold price (\u00A3)")+
   xlab("Relative improvement in quality of life if unrevised")+
   scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
   scale_y_continuous(limits=c(0, max.y))+
   theme_bw(base_size = 24)+ 
   theme(legend.title=element_blank())+
   facet_grid(. ~ op) +
   geom_label_repel() 

 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.qol.improvement.label.tiff",
        dpi=300,
        width = 12, height = 7)
 
 
# revision reduction assumption -----
 # data 
 load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/tkr.average.characteristics.tps.RData")
 tkr.tps<-tps
 rm(tps)
 
 tkr.tps<- tkr.tps %>% 
   filter(revision.reduction==1 |
            revision.reduction==0.5) %>% 
   mutate(revision.reduction=as.character(revision.reduction))
 
 
 load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.average.characteristics.tps.RData")
 thr.tps<-tps
 rm(tps)
 
 thr.tps<- thr.tps %>% 
   filter(revision.reduction==1 |
            revision.reduction==0.5) %>% 
   mutate(revision.reduction=as.character(revision.reduction))
 
 
 
 
 # summarise tkr transition probabilities 
 # get confidence intervals
 tkr.summary.prob.tp<-tkr.tps %>% 
   filter(sim!="deterministic") %>% 
   filter(time!=0) %>%
   group_by(time, group, revision.reduction) %>% 
   summarise(low.ci.tp=quantile(p.revision, 
                                probs = c(0.025), na.rm=TRUE),
             high.ci.tp=quantile(p.revision, 
                                 probs = c(0.975), na.rm=TRUE))
 
 # get deterministic (average across mi) 
 tkr.summary.det.tp<-tkr.tps %>% 
   filter(sim=="deterministic") %>%
   filter(time!=0) %>%
   group_by(time, group, revision.reduction)  %>%
   summarise(det.tp=mean(p.revision))
 
 tkr.summary.tp<-tkr.summary.det.tp %>% 
   left_join(tkr.summary.prob.tp,
             by=c("time", "group", "revision.reduction"))
 rm(tkr.summary.det.tp,tkr.summary.prob.tp)    
 
 # tkr.summary.tp<-tkr.summary.tp %>% 
 #   left_join(tkr.patient.profiles, by="group")
 
 
 # summarise thr transition probabilities 
 # get confidence intervals
 thr.summary.prob.tp<-thr.tps %>% 
   filter(sim!="deterministic") %>% 
   filter(time!=0) %>%
   group_by(time, group, revision.reduction) %>% 
   summarise(low.ci.tp=quantile(p.revision, 
                                probs = c(0.025), na.rm=TRUE),
             high.ci.tp=quantile(p.revision, 
                                 probs = c(0.975), na.rm=TRUE))
 
 # get deterministic (average across mi) 
 thr.summary.det.tp<-thr.tps %>% 
   filter(sim=="deterministic") %>%
   filter(time!=0) %>%
   group_by(time, group, revision.reduction)  %>%
   summarise(det.tp=mean(p.revision))
 
 thr.summary.tp<-thr.summary.det.tp %>% 
   left_join(thr.summary.prob.tp,
             by=c("time", "group", "revision.reduction"))
 rm(thr.summary.det.tp,thr.summary.prob.tp)    
 
 # thr.summary.tp<-thr.summary.tp %>% 
 #   left_join(thr.patient.profiles, by="group")
 
 
 
 # plot 
 # plot transition probability for tkr
 rbind(tkr.summary.tp %>% 
         mutate(op="Knee replacement"), 
       thr.summary.tp %>% 
         mutate(op="Hip replacement")) %>% 
   mutate(op=factor(op,levels=c("Knee replacement", 
                                "Hip replacement"))) %>% 
   mutate(revision.reduction=ifelse(revision.reduction=="1", 
                                    "Current", "Reduced\nby 50%")) %>% 
   mutate(revision.reduction=factor(revision.reduction,
                                    levels=c("Current", 
                                             "Reduced\nby 50%"))) %>% 
   ggplot()+
   facet_grid(.~op)+
   geom_point(aes(time, det.tp,
                  shape=revision.reduction),
              position=position_dodge(0.5)) +
   geom_errorbar(aes(time, 
                     ymin=low.ci.tp,
                     ymax=high.ci.tp,
                     linetype=revision.reduction),
                 colour="black", 
                 width=0,
                 position=position_dodge(0.5))+
   scale_y_continuous(labels = scales::percent)+ 
   scale_linetype(name="Risk of\nrevision")+
   scale_shape(name="Risk of\nrevision")+
   xlab("Years since surgery")+
   ylab("Transition probability")+
   theme_bw(base_size = 24)
 
 
 ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/CPRD HES lifetime risk of revision/plots/tkr and thr markov model/average_characteristics.revision.improvement.tiff",
        dpi=300,
        width = 12, height = 6)
 rm(tkr.summary.tp, tkr.tps,
    thr.summary.tp, thr.tps)
 
 
# revision reduction -----
 max.y # from qol improvement- so they have the same y axis
 
 QALYs.costs %>% 
   filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
  filter(qol.improvement==1) %>% 
  mutate(revision.reduction.name=1-revision.reduction) %>% 
  ggplot(aes(x=revision.reduction.name))+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative reduction in risk of revision (%)")+
   scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
   scale_y_continuous(limits=c(0, max.y))+
   theme_bw(base_size = 24)+ 
  theme(legend.title=element_blank())+
  facet_grid(. ~ op)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.revision.reduction.tiff",
          dpi=300,
              width = 12, height = 7)


#with label
QALYs.costs %>% 
  filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
  filter(qol.improvement==1) %>% 
  mutate(revision.reduction.name=1-revision.reduction) %>% 
  mutate(lab=ifelse(revision.reduction==0.5 | revision.reduction==1,
                    paste0("Lifetime risk of\nrevision: ",
                           format(round(mean.prop.revised*100,1),nsmall = 1), "%"),
                    NA)) %>% 
  ggplot(aes(x=revision.reduction.name,
             y=mean.inc.nmb,
             label=lab))+
  geom_point(colour=NA) +
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative reduction in risk of revision (%)")+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_y_continuous(limits=c(0, max.y))+
  theme_bw(base_size = 24)+ 
  theme(legend.title=element_blank())+
  facet_grid(. ~ op)+
  geom_label_repel(nudge_y = 1000, 
                   segment.color=NA) 

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.revision.reduction.tiff",
       dpi=300,
       width = 12, height = 7)



# # absolute risk 
# QALYs.costs %>% 
#   filter(qol.improvement==1) %>% 
#   ggplot(aes(x=mean.prop.revised))+
#   geom_line(aes(y=mean.inc.nmb))+
#   geom_ribbon(aes(ymin=low.ci.inc.nmb,
#                   ymax=high.ci.inc.nmb),
#               alpha=0.1)+
#   scale_linetype_manual(values=c("solid", "longdash"))+
#   ylab("Threshold price (\u00A3)")+
#   xlab("Absolute risk of revision (%)")+
#   scale_x_continuous(labels=percent)+
#  theme_bw(base_size = 18)+ 
#   theme(legend.title=element_blank())+
#   facet_grid(. ~ op)
# 
# # absolute reduction
# risk.revised<-QALYs.costs %>% 
#   group_by(op) %>% 
#   filter(qol.improvement==1) %>% 
#   filter(revision.reduction==1) %>% 
#   ungroup() %>% 
#   select(op, mean.prop.revised) %>% 
#   rename(risk.revised=mean.prop.revised)
# 
# QALYs.costs<-QALYs.costs %>% 
#   left_join(risk.revised,
#             by="op") %>% 
#   mutate(absolute.reduction=risk.revised-mean.prop.revised)
#   
# 
# QALYs.costs %>% 
#   filter(qol.improvement==1) %>% 
#  # mutate(revision.reduction.name=1-revision.reduction) %>% 
#   ggplot(aes(x=absolute.reduction))+
#   geom_line(aes(y=mean.inc.nmb))+
#   geom_ribbon(aes(ymin=low.ci.inc.nmb,
#                   ymax=high.ci.inc.nmb),
#               alpha=0.1)+
#   scale_linetype_manual(values=c("solid", "longdash"))+
#   ylab("Threshold price (\u00A3)")+
#   xlab("Absolute reduction in risk of revision (%)")+
#   scale_x_continuous(labels=percent)+
#  theme_bw(base_size = 18)+ 
#   theme(legend.title=element_blank())+
#   facet_grid(. ~ op)


# together -----
plot_grid(
  QALYs.costs %>% 
    filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
    filter(qol.improvement==1) %>% 
    mutate(revision.reduction.name=1-revision.reduction) %>% 
    ggplot(aes(x=revision.reduction.name))+
    geom_line(aes(y=mean.inc.nmb))+
    geom_ribbon(aes(ymin=low.ci.inc.nmb,
                    ymax=high.ci.inc.nmb),
                alpha=0.1)+
    scale_linetype_manual(values=c("solid", "longdash"))+
    scale_y_continuous(limits = c(0,11750))+
    ylab("Threshold price (\u00A3)")+
    xlab("Relative reduction in risk of revision (%)")+
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
    theme_bw(base_size = 24)+ 
    theme(legend.title=element_blank())+
    facet_grid(. ~ op)+ ggtitle("a) Improvements in risk of revision"),
  
  QALYs.costs %>% 
    filter(revision.reduction==1) %>% 
    mutate(qol.improvement.name=qol.improvement-1) %>% 
    ggplot(aes(x=qol.improvement.name))+
    geom_line(aes(y=mean.inc.nmb))+
    geom_ribbon(aes(ymin=low.ci.inc.nmb,
                    ymax=high.ci.inc.nmb),
                alpha=0.1)+
    scale_linetype_manual(values=c("solid", "longdash"))+
    scale_y_continuous(limits = c(0,11750))+
    ylab("Threshold price (\u00A3)")+
    xlab("Relative improvement in quality of life if unrevised")+
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
    theme_bw(base_size = 24)+ 
    theme(legend.title=element_blank())+
    facet_grid(. ~ op)+ ggtitle("b) Improvements in quality of life if unrevised ")
  , nrow=2
)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.revision.reduction.and.qol.tiff",
       dpi=300,
       width = 12, height = 14)

#with labels
plot_grid(
  QALYs.costs %>% 
    filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
    filter(qol.improvement==1) %>% 
    mutate(revision.reduction.name=1-revision.reduction) %>% 
    mutate(lab=ifelse(revision.reduction==0.5 | revision.reduction==1,
                      paste0("Lifetime risk of\nrevision: ",
                             format(round(mean.prop.revised*100,1),nsmall = 1), "%"),
                      NA)) %>% 
    ggplot(aes(x=revision.reduction.name,
               y=mean.inc.nmb,
               label=lab))+
    geom_point(colour=NA) +
    geom_line(aes(y=mean.inc.nmb))+
    geom_ribbon(aes(ymin=low.ci.inc.nmb,
                    ymax=high.ci.inc.nmb),
                alpha=0.1)+
    scale_linetype_manual(values=c("solid", "longdash"))+
    ylab("Threshold price (\u00A3)")+
    xlab("Relative reduction in risk of revision (%)")+
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
    scale_y_continuous(limits=c(0, max.y))+
    theme_bw(base_size = 24)+ 
    theme(legend.title=element_blank())+
    facet_grid(. ~ op)+
    geom_label_repel(nudge_y = 1000, 
                     segment.color=NA,
                     size = 5.5) + 
    ggtitle("a) Improvements in risk of revision"),
  
  QALYs.costs %>% 
    filter(revision.reduction==1) %>% 
    mutate(qol.improvement.name=qol.improvement-1) %>% 
    mutate(lab=ifelse(qol.improvement==1 | qol.improvement==1.05,
                      paste0("EQ-5D-3L index: ",
                             format(round(mean.unrevised.qol,2),nsmall = 1), 
                             "\nOKHS: ",
                             format(round(okhs,1),nsmall = 1)),
                      NA)) %>% 
    ggplot(aes(x=qol.improvement.name, y=mean.inc.nmb,
               label=lab))+
    geom_line()+
    geom_point(colour=NA)+
    geom_ribbon(aes(ymin=low.ci.inc.nmb,
                    ymax=high.ci.inc.nmb),
                alpha=0.1)+
    scale_linetype_manual(values=c("solid", "longdash"))+
    ylab("Threshold price (\u00A3)")+
    xlab("Relative improvement in quality of life if unrevised")+
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
    scale_y_continuous(limits=c(0, max.y))+
    theme_bw(base_size = 24)+ 
    theme(legend.title=element_blank())+
    facet_grid(. ~ op) +
    geom_label_repel(size = 5.5) +
    ggtitle("b) Improvements in quality of life if unrevised ")
  , nrow=2
)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.revision.reduction.and.qol.label.tiff",
       dpi=300,
       width = 12, height = 14)

#2 Age and gender -----
# data -----
rm(list=ls())

load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/tkr.age.gender.QALYs.costs.RData")
tkr.QALYs.costs<-QALYs.costs
rm(QALYs.costs)

load("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/model output/thr.age.gender.QALYs.costs.RData")
thr.QALYs.costs<-QALYs.costs
rm(QALYs.costs)

QALYs.costs<-rbind(tkr.QALYs.costs %>% 
                   mutate(op="Knee replacement"), 
                   thr.QALYs.costs %>% 
                   mutate(op="Hip replacement")) %>% 
  mutate(op=factor(op,levels=c("Knee replacement", 
                               "Hip replacement"))) %>% 
  mutate(gender=factor(gender, levels=c("Male", "Female")))

# lifetime risk of revision ------
a<-rbind(
QALYs.costs %>% 
  filter(revision.reduction==1 &
           qol.improvement==1),
QALYs.costs %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1)) %>% 
  mutate(type=ifelse(revision.reduction==1,
                     "Current",
                     "Reduced \nby 50%"))


  a %>% 
  ggplot(aes(x=age, group=type))+  #, group=gender
  facet_grid(gender ~ op)+
  geom_line(aes(y=mean.prop.revised,
                linetype=type))+  #,linetype=gender
  geom_ribbon(aes(ymin=low.ci.prop.revised,
                  ymax=high.ci.prop.revised),
              alpha=0.1)+
  scale_linetype_manual(name = "Revision risk",
                        values=c("solid", "longdash"))+
      scale_y_continuous(labels = scales::percent)+
  ylab("Lifetime risk of revision")+
    xlab("Age at surgery")+
  theme_bw(base_size = 18)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.revision.risk.tiff",
         dpi=300,width = 12, height = 7)


# change in qol------
a<-rbind(
QALYs.costs %>% 
  filter(revision.reduction==1 &
           qol.improvement==1),
QALYs.costs %>% 
  filter(revision.reduction==1 &
           qol.improvement==1.05)) %>% 
  mutate(type=ifelse(qol.improvement==1,
                     "Current",
                     "Improved \nby 5%"))



a<-rbind(
a %>% 
  ungroup() %>% 
  select(op,type, age, gender, q1_eq5d) %>% 
  rename(mean.qol=q1_eq5d) %>% 
  mutate(low.ci.qol=mean.qol,
         high.ci.qol=mean.qol,
         time=0),
  
a %>% 
  ungroup() %>% 
  select(op, type, age, gender, mean.unrevised.qol,
         low.ci.unrevised.qol,
         high.ci.unrevised.qol) %>% 
  rename(mean.qol=mean.unrevised.qol,
         low.ci.qol=low.ci.unrevised.qol,
         high.ci.qol=high.ci.unrevised.qol) %>% 
  mutate(time=6),
  
 a %>% 
  ungroup() %>% 
  select(op, type, age, gender, mean.unrevised.qol,
         low.ci.unrevised.qol,
         high.ci.unrevised.qol) %>% 
  rename(mean.qol=mean.unrevised.qol,
         low.ci.qol=low.ci.unrevised.qol,
         high.ci.qol=high.ci.unrevised.qol) %>% 
  mutate(time=12) 
)
 
  
a %>% 
    filter(age==50| age==65| age==80) %>% 
  ggplot(aes(x=time, group=type))+
  facet_grid(age+ gender ~ op)+
  geom_line(aes(y=mean.qol, 
            linetype=type))+  
  geom_ribbon(aes(ymin=low.ci.qol,
                  ymax=high.ci.qol),
              alpha=0.1)+
  scale_linetype_manual(name = "Post-operative \nquality of life",
                        values=c("solid", "longdash"))+
  ylab("EQ-5D-3L over year of primary surgery (if unrevised)")+
  xlab("Time (months)")+
  scale_x_continuous(limits = c(0, 12),
                   breaks= c(0,3,6,9,12))+
  theme_bw(base_size = 18)


ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.qol.improvement.tiff",
         dpi=300,
             width = 8.27, height = 11.7)









# heatmap -----
# shows the mean threshold price
# for combined improvements in revision risk and qol
QALYs.costs %>% 
  filter(revision.reduction>=0.5) %>%  #up to a 50% reduction
  filter(age==50| age==65| age==80) %>% 
 # filter(gender=="Female") %>% 
  arrange(desc(revision.reduction),qol.improvement) %>% 
  ggplot()+
  geom_tile(aes(x=abs(revision.reduction-1),
                y=qol.improvement-1,
                fill=mean.inc.nmb)) +
   # scale_x_reverse()+
  xlab("Relative reduction in revision risk")+
  ylab("Relative improvement in quality of life if unrevised")+
  labs(fill = "Threshold \nprice (\u00A3)")+
  scale_fill_gradientn(colours = 
                         c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))+
  theme_bw(base_size = 18)+ 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand=c(0,0))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand=c(0,0))+
  facet_grid(age+gender ~ op)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.heatmap.tiff",
         dpi=300,
             width = 8.27, height = 11.7)


# qol improvements -----
QALYs.costs %>% 
  filter(revision.reduction==1) %>% 
  filter(age==55| age==65 |age==75) %>% 
  mutate(age=as.character(age)) %>% 
 # filter(gender=="Female") %>% 
  mutate(qol.improvement.name=qol.improvement-1) %>% 
  ggplot(aes(x=qol.improvement.name, 
             group=age,
             #colour=age,
             #fill=age,
             linetype=age))+
    facet_grid(gender ~ op)+
  #facet_grid(age+gender ~ op)+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
 # scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative improvement in quality of life if unrevised")+
  scale_x_continuous(labels=percent)+
  theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.qol.improvement.tiff",
         dpi=300,
         width = 12, height = 7)




# age and qol improvements -----


QALYs.costs %>% 
  filter(revision.reduction==1) %>% 
  filter(qol.improvement==1.05) %>% 
  ggplot(aes(x=age))+  #, group=gender
  facet_grid(. ~ op)+
  geom_line(aes(y=mean.inc.nmb))+  #,linetype=gender
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3) if quality of life \nif unrevised were increased by 5%")+
  xlab("Age at surgery")+
  theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(gender ~ op)

# ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.qol.improvement.5percent.tiff",
#          dpi=300,
#              width = 12, height = 7)


# revision reduction -----
QALYs.costs %>% 
  filter(qol.improvement==1) %>% 
  filter(age==55| age==65 |age==75) %>% 
  mutate(age=as.character(age)) %>% 
  mutate(revision.reduction.name=1-revision.reduction) %>% 
  ggplot(aes(x=revision.reduction.name, 
             group=age,
             #colour=age,
             #fill=age,
             linetype=age))+
    facet_grid(gender ~ op)+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative reduction in risk of revision (%)")+
  scale_x_continuous(labels=percent)+
 theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.revision.reduction.tiff",
         dpi=300,
         width = 12, height = 7)



QALYs.costs %>% 
  filter(qol.improvement==1) %>% 
    filter(age==55| age==65 |age==75) %>% 
  mutate(revision.reduction.name=1-revision.reduction) %>% 
    mutate(age=as.character(age)) %>% 
  ggplot(aes(x=revision.reduction.name,
             group=age,
             linetype=age))+
  geom_line(aes(y=mean.inc.nmb))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
 # scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3)")+
  xlab("Relative reduction in risk of revision (%)")+
  scale_x_continuous(labels=percent)+
 theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(gender ~ op)


# age and revision reduction -----

QALYs.costs %>% 
  filter(revision.reduction==0.5) %>% 
  filter(qol.improvement==1) %>% 
  ggplot(aes(x=age, group=gender))+
  facet_grid(. ~ op)+
  geom_line(aes(y=mean.inc.nmb))+  #,linetype=gender
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
 # scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (\u00A3) if \nrevision risk were reduced by 50%")+
  xlab("Age at surgery")+
  theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(gender ~ op)

# ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.revision.reduction.50percent.tiff",
#          dpi=300,
#              width = 12, height = 7)


# 5% qol improvement and 50% revision reduction----

# QALYs.costs<-rbind(
#   tkr.QALYs.costs %>% 
#   mutate(op="TKR"),
#   thr.QALYs.costs %>% 
#   mutate(op="THR")) %>% 
#   mutate(op=factor(op,levels=c("TKR", "THR")))

QALYs.costs<-QALYs.costs %>% 
   mutate(type=
            ifelse(qol.improvement==1.05, 
                   "Quality of life\nif unrevised\nimproved by 5%\n", 
            ifelse(revision.reduction==0.5, 
                   "Risk of revision\nreduced by 50%\n",  NA      
                   )))



rbind(QALYs.costs %>% 
  filter(revision.reduction==0.5 &
           qol.improvement==1), 
  QALYs.costs %>% 
  filter(revision.reduction==1 &
           qol.improvement==1.05)) %>% 
  ggplot(aes(x=age, 
             group=type))+
  geom_line(aes(y=mean.inc.nmb, 
             linetype=type))+
  geom_ribbon(aes(ymin=low.ci.inc.nmb,
                  ymax=high.ci.inc.nmb),
              alpha=0.1)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  ylab("Threshold price (£)")+
  xlab("Age at surgery")+
 theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(gender ~ op)

  
ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.revision.reduction.50percent.qol.5percent.tiff",
         dpi=300,width = 12, height = 7)

# specific estimates -----

# all for 65 year old woman
# lifetime risk
QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1 ) %>%
  filter(age=="65", gender=="Female") %>% 
  mutate(est=paste0(sprintf("%.2f",mean.prop.revised*100),
                    "% (",
                   sprintf("%.2f",low.ci.prop.revised*100),
                   "% to ",
                   sprintf("%.2f",high.ci.prop.revised*100),
                   "%)")) %>% 
    select(op, est)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1 ) %>%
  filter(age=="50", gender=="Female") %>% 
  mutate(est=paste0(sprintf("%.2f",mean.prop.revised*100),
                    "% (",
                   sprintf("%.2f",low.ci.prop.revised*100),
                   "% to ",
                   sprintf("%.2f",high.ci.prop.revised*100),
                   "%)")) %>% 
    select(op, est)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1 ) %>%
  filter(age=="65", gender=="Female") %>% 
  mutate(est=paste0(sprintf("%.2f",mean.unrevised.qol),
                    " (",
                   sprintf("%.2f",low.ci.unrevised.qol),
                   " to ",
                   sprintf("%.2f",high.ci.unrevised.qol),
                   ")")) %>% 
    select(op, est)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1.05 ) %>%
  filter(age=="65", gender=="Female") %>% 
  mutate(est=paste0(sprintf("%.2f",mean.unrevised.qol),
                    " (",
                   sprintf("%.2f",low.ci.unrevised.qol),
                   " to ",
                   sprintf("%.2f",high.ci.unrevised.qol),
                   ")")) %>% 
    select(op, est)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==0.5 &
          qol.improvement==1 ) %>%
  filter(age=="50"| age=="65"|age=="80") %>% 
  filter(gender=="Female") %>% 
  mutate(est=threshold.price) %>% 
    select(op, age,est)

QALYs.costs %>% 
  ungroup() %>% 
  filter(revision.reduction==1 &
          qol.improvement==1.05 ) %>%
  filter(age=="50"| age=="65"|age=="80") %>% 
  filter(gender=="Female") %>% 
  mutate(est=threshold.price) %>% 
    select(op, age,est)
