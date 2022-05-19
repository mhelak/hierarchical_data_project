
# general data inquiry and result summarization from simulations


setwd("Documents/work/coursera_stats/2_hierarchical_data/")
library(tidyverse)
library(nlme)
library(Rmisc)

pilot<-read.csv("project/G6.pilot.data.csv") %>%
  mutate(garden=as.factor(garden),
         species=as.factor(species),
         rater=as.factor(rater),
         combi=paste0("species",species,"_garden",garden,"_rater",rater))

ggplot(pilot, aes(x=tot.vase.days))+geom_density() #the distribution seems bimodal, hence simply modelling a lambda+1 might generate lower variance
summary(pilot$tot.vase.days)

ggplot(pilot, aes(x=combi, y=tot.vase.days))+geom_violin()

pilot %>%
  group_by(species) %>%
  summarise(n=mean(tot.vase.days))


######################
# GLM sample size ####
######################


simlist_GLM<-readRDS("project/glm_sim_out_effsize2.rds")

simmodel_out<-bind_rows(simlist_GLM, .id = "model") %>%
  mutate(model=gsub("glm_Hochberg", "glm_Benjamini_Hochberg", model),
         model=factor(model, levels=c("glm_Bonferroni", "glm_Holm", "glm_Benjamini_Hochberg")),
         power=round(power, digits=2)) 

summaryMin<-simmodel_out %>% filter(power>=0.85) %>% group_by(model) %>% slice_min(n) #248 per group.
summaryMin2<-simmodel_out %>% filter(power>=0.80) %>% group_by(model) %>% slice_min(n) #248 per group.
summaryMin3<-simmodel_out %>% filter(power>=0.9) %>% group_by(model) %>% slice_min(n) #248 per group.
summaryMin4<-simmodel_out %>% filter(power>=0.95) %>% group_by(model) %>% slice_min(n) #248 per group.

ggplot(simmodel_out %>% filter(n<225),aes(x=n, y=power, group=model))+
  geom_hline(yintercept = 0.95, linetype="dashed")+
  geom_hline(yintercept = 0.8, linetype="dashed")+geom_point(aes(color=model))+
  geom_hline(yintercept = 0.9, linetype="dashed")+geom_point(aes(color=model))+
  geom_point(aes(color=model))+geom_line(aes( color=model)) +
  geom_label(aes(label=n), data=summaryMin3, nudge_y = 0.015, color="black") +
  geom_label(aes(label=n), data=summaryMin2, nudge_y = 0.015, color="black") +
  geom_label(aes(label=n), data=summaryMin4, nudge_y = 0.015, color="black")+
  scale_color_manual(values=c("#59A5D8","#386FA4","#133C55"))+theme_bw()

ggsave("project/additive_model_sim_BH_effsize2.png", height=4.5, width=8)


# sample size calculations with n=96/compound
min=3
compounds=15
n=96*compounds
n_perSpecies = n/2
minPP = min* n_perSpecies / 6
minPP / 60 #13.3h working day


summaryMin2<-simmodel_out %>% filter(power>0.8) %>% group_by(model) %>% slice_min(n)

simmodel_out %>% filter(power>0.78) %>% group_by(model) %>% slice_min(n)




#######################
# GLMM sample size ####
#######################

glmm_outlist<-list.files("project/glmm_sepNpowsim", full.names = T)

simmodel_out_GLMM<-bind_rows(lapply(glmm_outlist, function(x){
  readRDS(x) %>% 
    bind_rows %>% 
    mutate(file=x)})) %>%
  mutate(model=ifelse(grepl("Bonf",file),"glmm_Bonferroni","glmm_Holm"),
         model=factor(model, levels=c("glmm_Bonferroni", "glmm_Holm")),
         power=round(power, digits=2)) 

summaryMin_glmm<-simmodel_out_GLMM %>% filter(power>=0.85) %>% group_by(model) %>% slice_min(n) #248 per group.
summaryMin_glmm2<-simmodel_out_GLMM %>% filter(power>=0.8 & power < 0.85) %>% group_by(model) %>% slice_min(n) #248 per group


ggplot(simmodel_out_GLMM,aes(x=n, y=power, group=model))+
  geom_hline(yintercept = 0.85, linetype="dashed")+geom_point(aes(color=model))+
  geom_hline(yintercept = 0.8, linetype="dashed")+geom_point(aes(color=model))+
  geom_line(aes( color=model)) + #+theme_bw()
  geom_label(aes(label=n), data=summaryMin_glmm, nudge_y = 0.02, color="black")+
  geom_label(aes(label=n), data=summaryMin_glmm2, nudge_y = 0.02, color="black") +
  ylim(c(0.55,NA))+
  scale_color_manual(values=c("#59A5D8","#386FA4"))+theme_bw()
  

ggsave("project/mixed_model_sim_BH.png", height=4.5, width=7)






# ok we decided on n = 80 per compound; how to distribute these across the gardens, raters, subplots and bushes?
n_perCompoound=96
n_compounds=15
n_perSpeciesFinal=n_perCompoound/2

n_total=n_perSpeciesFinal*compounds

n_raters=6
bushes=n_total/n_compounds #dooh..
bushes_perGarden=bushes/2


bushes_perSubplot=bushes_perGarden/8
bushes_perRater=bushes/n_raters


summary_perSpecies<-data.frame(garden=rep(1,))

40/6


96/2
24/8


