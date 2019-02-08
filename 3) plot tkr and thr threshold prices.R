rm(list=ls())

library(dplyr)
library(ggplot2)
library(scales)

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



# heatmap -----
# shows the mean threshold price
# for combined improvements in revision risk and qol
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

# ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.heatmap.tiff",
#          dpi=300,
#              width = 12, height = 7)
# qol improvements -----
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
  scale_x_continuous(labels=percent)+
  theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(. ~ op)

# ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.qol.improvement.tiff",
#          dpi=300,
#              width = 12, height = 7)

# revision reduction -----
QALYs.costs %>% 
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
  scale_x_continuous(labels=percent)+
 theme_bw(base_size = 18)+ 
  theme(legend.title=element_blank())+
  facet_grid(. ~ op)

# ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/average.characteristics.revision.reduction.tiff",
#          dpi=300,
#              width = 12, height = 7)


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
  ylab("EQ-5D over year of primary surgery (if unrevised)")+
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
  scale_fill_gradient(low = "white", high = "black")+
  theme_bw(base_size = 18)+ 
  scale_x_continuous(labels=percent)+
  scale_y_continuous(labels=percent)+
  facet_grid(age+gender ~ op)

ggsave("C:/Users/Ed/Dropbox/DPhil data cprd hes analysis/threshold prices tkr thr improved effectiveness/plots/age.gender.heatmap.tiff",
         dpi=300,
             width = 8.27, height = 11.7)


# qol improvements -----
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
