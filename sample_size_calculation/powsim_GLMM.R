

library(lme4)
library(tidyverse)

setwd("/home/zane/RNA_seq_scripts/random_learning/")

pilot<-read.csv("G6.pilot.data.csv") %>%
  mutate(garden=as.factor(garden),
         species=as.factor(species),
         rater=as.factor(rater))


#version without garden
glm_out_noGard <- glmer(tot.vase.days ~  species  + (1|rater), family=poisson(link = "log"), data=pilot)
summary(glm_out_noGard)


#this one allows to pull out the intercepts and use these directly for the simulations
b_rater_noG <- ranef(glm_out_noGard)$rater %>%
  as_tibble(rownames = 'r') %>%
  dplyr::rename(int = `(Intercept)`)

b0=fixef(glm_out_noGard)[1]
b1=fixef(glm_out_noGard)[2]






simulate_glmm<-function(n=408, alpha=0.05, initial_pow=1, power_stop=0.6, effect_size=1, N=1000, method="Holm", n_comparisons=14){
  
  sim_list<-list()
  pow=initial_pow
  set.seed(6247555)
  
  while (pow>power_stop){
    values<-data.frame(power=NA, n=NA) 
    pval <- matrix(0, nrow = N, ncol = n_comparisons)
    n=n-8
    for(i in 1:N) {
      
      #null model
      Y0_s1_r1<-rpois(n/4, lambda=exp(b0  + b_rater_noG$int[b_rater_noG$r==1]))
      Y0_s1_r2<-rpois(n/4, lambda=exp(b0  + b_rater_noG$int[b_rater_noG$r==2]))
      
      Y0_s2_r1<-rpois(n/4, lambda=exp(b0 + b1 + b_rater_noG$int[b_rater_noG$r==1]))
      Y0_s2_r2<-rpois(n/4, lambda=exp(b0 + b1 + b_rater_noG$int[b_rater_noG$r==2]))
      
      
      #lambda=lambda+1
      treatment_list<-list()
      for (c in 1:14){
        
        Y1_s1_r1<-rpois(n/4, lambda=exp(b0  + b_rater_noG$int[b_rater_noG$r==1])+effect_size)
        Y1_s1_r2<-rpois(n/4, lambda=exp(b0  + b_rater_noG$int[b_rater_noG$r==2])+effect_size)
        
        Y1_s2_r1<-rpois(n/4, lambda=exp(b0 + b1 + b_rater_noG$int[b_rater_noG$r==1])+effect_size)
        Y1_s2_r2<-rpois(n/4, lambda=exp(b0 + b1 + b_rater_noG$int[b_rater_noG$r==2])+effect_size)
        
        
        treatment_list[[paste0("compound",c)]]<-data.frame(species = as.factor(c(rep(1, n/2), rep(2, n/2))), 
                                                           rater = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 2)),
                                                           tot.vase.days=c(Y1_s1_r1, Y1_s1_r2, Y1_s2_r1, Y1_s2_r2),
                                                           type="sim", condition=as.factor(c(rep(paste0("compound",c), n))))
        
      }
      
      compounds_df<-bind_rows(treatment_list, .id="compound") %>%
        bind_rows(data.frame(species = as.factor(c(rep(1, n/2), rep(2, n/2))), 
                             rater = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 2)),
                             tot.vase.days=c(Y0_s1_r1, Y0_s1_r2, Y0_s2_r1, Y0_s2_r2),
                             type="sim", condition=as.factor(c(rep("ctrl", n)))))
      
      compounds_df$condition<-relevel(compounds_df$condition, "ctrl")
      
      glm_out_cond <- summary(glmer(tot.vase.days ~ condition + species + (1|rater), family=poisson, data=compounds_df))
      
      pval[i,] <- glm_out_cond$coeff[2:(n_comparisons+1), "Pr(>|z|)"]/2
      
      # Multiplicity adjustment(s) 
      if (method == "Bonferroni"){
        pval[i,] = p.adjust(pval[i,], method = "bonferroni")
      }
      if (method == "Holm"){
        pval[i,] = p.adjust(pval[i,], method = "holm")
      }
      if (method == "Benjamini-Hochberg"){
        pval[i,] = p.adjust(pval[i,], method = "BH")
      }
    }
    # 4. Estimate power
    pow<- mean(apply(pval, 2,FUN = function(x){sum(x < alpha)}))/N
    print(paste(pow, n))
    sim_list[[paste0("n",n)]]<-data.frame(power=pow, n=n) 
  }
  return(bind_rows(sim_list))
}


glmm_Holm<-simulate_glmm()
glmm_Bonferroni<-simulate_glmm(method = "Bonferroni")

saveRDS(list(glmm_Holm=glmm_Holm, glmm_Bonferroni=glmm_Bonferroni), "glmm_sim_out.rds")



