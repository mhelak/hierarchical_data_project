
library(tidyverse)
library(nlme)
library(Rmisc)


setwd("Documents/work/coursera_stats/2_hierarchical_data/")
pilot<-read.csv("project/G6.pilot.data.csv") %>%
  mutate(garden=as.factor(garden),
         species=as.factor(species),
         rater=as.factor(rater))



#############################################
# GET AN ESTIMATE OF PREDICTOR VARIABLES ####
#############################################


glm_out <- glm(tot.vase.days ~ species + rater + garden, family=poisson(link = "log"), data=pilot)
summary(glm_out)


glm_out_qp <- glm(tot.vase.days ~ species + rater + garden, family=quasipoisson(link = "log"), data=pilot)
summary(glm_out_qp)

#variance= scaler*mean

b0=glm_out$coefficients[1]
b1_s=glm_out$coefficients[2]
b2_r=glm_out$coefficients[3]
b3_g=glm_out$coefficients[4]


simulate_glm<-function(n=408, alpha=0.05, initial_pow=1, power_stop=0.6, effect_size=1, N=1000, method="Holm", n_comparisons=14){
  
  sim_list<-list()
  pow=initial_pow
  set.seed(6247555)
  
  while (pow>power_stop){
    values<-data.frame(power=NA, n=NA) 
    pval <- matrix(0, nrow = N, ncol = n_comparisons)
    n=n-8
    for(i in 1:N) {
      
      #null model
      Y0_s1_r1_g1<-rpois(n/8, lambda=exp(b0))
      Y0_s1_r1_g2<-rpois(n/8, lambda=exp(b0 + b3_g))
      
      Y0_s2_r1_g1<-rpois(n/8, lambda=exp(b0 + b1_s))
      Y0_s2_r1_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b3_g))
      
      Y0_s1_r2_g1<-rpois(n/8, lambda=exp(b0 + b2_r))
      Y0_s1_r2_g2<-rpois(n/8, lambda=exp(b0 + b2_r + b3_g))
      
      Y0_s2_r2_g1<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r))
      Y0_s2_r2_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r + b3_g))
      
      
      #lambda=lambda+1
      treatment_list<-list()
      for (c in 1:14){
        
        Y1_s1_r1_g1<-rpois(n/8, lambda=exp(b0)+effect_size)
        Y1_s1_r1_g2<-rpois(n/8, lambda=exp(b0 + b3_g)+effect_size)
        
        Y1_s2_r1_g1<-rpois(n/8, lambda=exp(b0 + b1_s)+effect_size)
        Y1_s2_r1_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b3_g)+effect_size)
        
        Y1_s1_r2_g1<-rpois(n/8, lambda=exp(b0 + b2_r)+effect_size)
        Y1_s1_r2_g2<-rpois(n/8, lambda=exp(b0 + b2_r + b3_g)+effect_size)
        
        Y1_s2_r2_g1<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r)+effect_size)
        Y1_s2_r2_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r + b3_g)+effect_size)
        
        
        treatment_list[[paste0("compound",c)]]<-data.frame(species = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 2)), 
                                                           rater = as.factor(rep(c(rep(1, n/2), rep(2, n/2)), 1)),
                                                           garden = as.factor(rep(c(rep(1, n/8), rep(2, n/8)), 4)),
                                                           tot.vase.days=c(Y1_s1_r1_g1, 
                                                                           Y1_s1_r1_g2, Y1_s2_r1_g1, Y1_s2_r1_g2, 
                                                                           Y1_s1_r2_g1, Y1_s1_r2_g2,
                                                                           Y1_s2_r2_g1, Y1_s2_r2_g2),
                                                           type="sim", condition=as.factor(c(rep(paste0("compound",c), n))))
        
      }
      
      compounds_df<-bind_rows(treatment_list, .id="compound") %>%
        bind_rows(data.frame(species = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 2)), 
                             rater = as.factor(rep(c(rep(1, n/2), rep(2, n/2)), 1)),
                             garden = as.factor(rep(c(rep(1, n/8), rep(2, n/8)), 4)),
                             tot.vase.days=c(Y0_s1_r1_g1, Y0_s1_r1_g2, Y0_s2_r1_g1, 
                                             Y0_s2_r1_g2, Y0_s1_r2_g1, Y0_s1_r2_g2, 
                                             Y0_s2_r2_g1, Y0_s2_r2_g2),
                             type="sim", condition=as.factor(c(rep("ctrl", n)))))
      
      compounds_df$condition<-relevel(compounds_df$condition, "ctrl")
      
      
      glm_out1 <- summary(glm(tot.vase.days ~ condition + species + rater + garden, family=poisson(link = "log"), data=compounds_df))
      
      
      pval[i,] <- glm_out1$coeff[2:(n_comparisons+1), "Pr(>|z|)"]/2
      
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
    # 4. Estimate power: for each compound (column) calculate proportion of times of significance (power for each compound) and take the mean across compounds (mean power across compounds)
    pow<- mean(apply(pval, 2,FUN = function(x){sum(x < alpha)}))/N
    print(paste(pow, n))
    sim_list[[paste0("n",n)]]<-data.frame(power=pow, n=n) 
  }
  return(bind_rows(sim_list))
}


glm_Holm<-simulate_glm()
glm_Bonferroni<-simulate_glm(method = "Bonferroni")
glm_Hochberg<-simulate_glm(method = "Benjamini-Hochberg")


saveRDS(list(glm_Holm=glm_Holm, glm_Bonferroni=glm_Bonferroni, glm_Hochberg=glm_Hochberg), "project/glm_sim_out.rds")



