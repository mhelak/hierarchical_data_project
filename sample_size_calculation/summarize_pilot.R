
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




simlist_GLM<-readRDS("project/glm_sim_out.rds")

simmodel_out<-bind_rows(simlist_GLM, .id = "model") %>%
  mutate(model=gsub("glm_Hochberg", "glm_Benjamini_Hochberg", model),
         model=factor(model, levels=c("glm_Bonferroni", "glm_Holm", "glm_Benjamini_Hochberg")))

summaryMin<-simmodel_out %>% filter(power>0.85) %>% group_by(model) %>% slice_min(n) #248 per group.

ggplot(simmodel_out,aes(x=n, y=power, group=model))+
  geom_hline(yintercept = 0.85, linetype="dashed")+geom_point(aes(color=model))+geom_line(aes( color=model)) + #+theme_bw()
  geom_label(aes(label=n), data=summaryMin, nudge_y = 0.02, color="black")

ggsave("project/additive_model_sim_BH.png", height=4.5, width=7)


# sample size calculations with n=248
min=3
compounds=15
n=232*compounds
n_perSpecies = n/2
minPP = min* n_perSpecies / 7
minPP / 60 #13.3h working day


summaryMin2<-simmodel_out %>% filter(power>0.8) %>% group_by(model) %>% slice_min(n)

