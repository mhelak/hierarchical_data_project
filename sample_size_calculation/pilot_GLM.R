
# checkup pilot data
# this was my first attempt, there is no balancing or random effects. Maybe with only balancing this could be already a good way (although of course in theory, the closer we get to the final model the better)

library(nlme)
library(Rmisc)

pilot<-read.csv("project/G6.pilot.data.csv") %>%
  mutate(garden=as.factor(garden),
         species=as.factor(species),
         rater=as.factor(rater))

library(tidyverse)
ggplot(pilot, aes(x=tot.vase.days))+geom_density() #the distribution seems bimodal, hence simply modelling a lambda+1 might generate lower variance
summary(pilot$tot.vase.days)

# lambda 10.78

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

# don't know yet what are the distributions, but its binary and looks approx randomly distributed
table(pilot$garden, pilot$species, pilot$rater)



##########################
# SIMULATE AN EXAMPLE ####
##########################


# Use rpois to make balanced draws from a poisson distribution
n=200
effect_size=1


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
Y1_s1_r1_g1<-rpois(n/8, lambda=exp(b0)+effect_size)
Y1_s1_r1_g2<-rpois(n/8, lambda=exp(b0 + b3_g)+effect_size)

Y1_s2_r1_g1<-rpois(n/8, lambda=exp(b0 + b1_s)+effect_size)
Y1_s2_r1_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b3_g)+effect_size)

Y1_s1_r2_g1<-rpois(n/8, lambda=exp(b0 + b2_r)+effect_size)
Y1_s1_r2_g2<-rpois(n/8, lambda=exp(b0 + b2_r + b3_g)+effect_size)

Y1_s2_r2_g1<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r)+effect_size)
Y1_s2_r2_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r + b3_g)+effect_size)



new_df_balanced_eff<-data.frame(species = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 4)), 
                                rater = as.factor(rep(c(rep(1, n/2), rep(2, n/2)), 2)),
                                garden = as.factor(rep(c(rep(1, n/8), rep(2, n/8)), 8)),
                                tot.vase.days=c(Y0_s1_r1_g1, Y0_s1_r1_g2, Y0_s2_r1_g1, 
                                                Y0_s2_r1_g2, Y0_s1_r2_g1, Y0_s1_r2_g2, 
                                                Y0_s2_r2_g1, Y0_s2_r2_g2,Y1_s1_r1_g1, Y1_s1_r1_g2,
                                                Y1_s2_r1_g1, Y1_s2_r1_g2, Y1_s1_r2_g1, Y1_s1_r2_g2,
                                                Y1_s2_r2_g1, Y1_s2_r2_g2),
                                type="sim", condition=c(rep("ctrl", n), rep("compound1", n))) %>% 
  mutate(condition=factor(condition, levels=c("ctrl","compound1")))

ggplot(new_df_balanced_eff, aes(color=condition, x=tot.vase.days))+geom_density()

glm_out1 <- glm(tot.vase.days ~ condition + species + rater + garden, family=poisson(link = "log"), data=new_df_balanced_eff)

summary(glm_out1)
exp(glm_out1$coefficients[2])
summary(glm_out1)$coefficients[2,4]

# confint
exp(confint(glm_out1, level = 0.90))
exp(confint(glm_out1, level = 0.90))[2,1]>=1 #check how often this is true?



pglm_out1 <- glm(tot.vase.days ~ condition + species + rater + garden, family=quasipoisson(link = "log"), data=new_df_balanced_eff)
summary(pglm_out1)




###########################
# DO POWER SIMULATIONS ####
###########################


set.seed(6247555)

n=400
effect_size=1
N=1000
pow=1
n_comparisons=14
alpha=0.05
alpha_corrected=alpha/n_comparisons

sim_list<-list()

while (pow>0.6){
  values<-data.frame(pval=NA, beta_1=NA, CI_lower=NA, power=NA, CI_prop=NA) 
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
    Y1_s1_r1_g1<-rpois(n/8, lambda=exp(b0)+effect_size)
    Y1_s1_r1_g2<-rpois(n/8, lambda=exp(b0 + b3_g)+effect_size)
    
    Y1_s2_r1_g1<-rpois(n/8, lambda=exp(b0 + b1_s)+effect_size)
    Y1_s2_r1_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b3_g)+effect_size)
    
    Y1_s1_r2_g1<-rpois(n/8, lambda=exp(b0 + b2_r)+effect_size)
    Y1_s1_r2_g2<-rpois(n/8, lambda=exp(b0 + b2_r + b3_g)+effect_size)
    
    Y1_s2_r2_g1<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r)+effect_size)
    Y1_s2_r2_g2<-rpois(n/8, lambda=exp(b0 + b1_s + b2_r + b3_g)+effect_size)
    
    
    new_df_balanced_eff<-data.frame(species = as.factor(rep(c(rep(1, n/4), rep(2, n/4)), 4)), 
                                    rater = as.factor(rep(c(rep(1, n/2), rep(2, n/2)), 2)),
                                    garden = as.factor(rep(c(rep(1, n/8), rep(2, n/8)), 8)),
                                    tot.vase.days=c(Y0_s1_r1_g1, Y0_s1_r1_g2, Y0_s2_r1_g1, 
                                                    Y0_s2_r1_g2, Y0_s1_r2_g1, Y0_s1_r2_g2, 
                                                    Y0_s2_r2_g1, Y0_s2_r2_g2,Y1_s1_r1_g1, 
                                                    Y1_s1_r1_g2, Y1_s2_r1_g1, Y1_s2_r1_g2, 
                                                    Y1_s1_r2_g1, Y1_s1_r2_g2,
                                                    Y1_s2_r2_g1, Y1_s2_r2_g2),
                                    type="sim", condition=c(rep("ctrl", n), rep("compound1", n))) %>% 
      mutate(condition=factor(condition, levels=c("ctrl","compound1")))
    
    
    glm_out1 <- glm(tot.vase.days ~ condition + species + rater + garden, family=poisson(link = "log"), data=new_df_balanced_eff)
    
    values[i,2]<-exp(glm_out1$coefficients[2]) #getting beta 1 (condition)
    values[i,1]<-summary(glm_out1)$coefficients[2,4]/2 #getting pvalue 
    values[i,3]<-exp(confint(glm_out1, level = 0.90))[2,1]

  } 
  signif=values %>% filter(pval<alpha_corrected) 
  pow = dim(signif)[1]/dim(values)[1]
  values$power<-pow
  
  largeB=values %>% filter(CI_lower>=1) 
  betaConf = dim(largeB)[1]/dim(values)[1]
  values$CI_prop<-betaConf
  
  sim_list[[paste0("n",n)]]<-values
  
}

sim_df<-bind_rows(sim_list, .id="n")
ggplot(sim_df, aes(x=beta_1))+geom_density()


#for each power level, do 20x, then calculate mean power and CI of the power estimate, plot and connect
# we also want to know at which point the upper CI of the beta_1 (lambda) is above 1
# and we want to add multiple testing correction

sim_summarized<-sim_df %>%
  group_by(n) %>%
  dplyr::summarise(mean_beta_1=mean(beta_1),
            lowerCI_above1=CI_prop[1],
            power=power[1],
            n=n[1]) %>%
  mutate(n=as.integer(gsub("n","",n)))
  
sim_summarized %>% filter(power>0.85) %>% slice_min(n) #160 per group.

cowplot::plot_grid(ggplot(sim_summarized,aes(x=n, y=power))+geom_hline(yintercept = 0.85)+geom_point(),
                   ggplot(sim_summarized,aes(x=n, y=lowerCI_above1))+geom_vline(xintercept = 328)+geom_point()+ylab("% runs with lower 5% CI bound >=1")
)

ggsave("project/additive_model_sim.png", height=4.5, width=10)







