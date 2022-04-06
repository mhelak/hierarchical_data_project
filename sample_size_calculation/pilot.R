
# checkup pilot data

pilot<-read.csv("Documents/work/coursera_stats/2_hierarchical_data/project/G6.pilot.data.csv")

library(tidyverse)
ggplot(pilot, aes(x=tot.vase.days))+geom_density()
summary(pilot$tot.vase.days)

# lambda 10.78

#############################################
# GET AN ESTIMATE OF PREDICTOR VARIABLES ####
#############################################


glm_out <- glm(tot.vase.days ~ garden + species + rater, family=poisson(link = "log"), data=pilot)
summary(glm_out)


b0=glm_out$coefficients[1]
b1=glm_out$coefficients[2]
b2=glm_out$coefficients[3]
b3=glm_out$coefficients[4]

# don't know yet what are the distributions, but its binary and looks approx randomly distributed
table(pilot$garden, pilot$species, pilot$rater)



##########################
# SIMULATE AN EXAMPLE ####
##########################


# Use rpois to make random draws from a poisson distribution
n=200
effect_size=1

new_df<-data.frame(garden=sample(pilot$garden, size = n, replace = T),
                   species=sample(pilot$species, size = n, replace = T),
                   rater=sample(pilot$rater, size = n, replace = T))


Y0 <- rpois(n/2, lambda=exp(b0 + b1*new_df$garden[1:n/2] + b2*new_df$species[1:n/2] + b3*new_df$rater[1:n/2]))

summary(Y0)

Y1 <- rpois(n/2, lambda=exp(b0 + b1*new_df$garden[n/2+1:nrow(new_df)] + b2*new_df$species[n/2+1:nrow(new_df)] + b3*new_df$rater[n/2+1:nrow(new_df)])+effect_size)

summary(Y1)


combined<-new_df %>% 
  mutate(condition=as.factor(c(rep("water",n/2), rep("compound1",n/2))),
         tot.vase.days=c(Y0, Y1))

combined$condition<-relevel(combined$condition, "water")

glm_out1 <- glm(tot.vase.days ~ condition + garden + species + rater, family=poisson(link = "log"), data=combined)

summary(glm_out1)
exp(glm_out1$coefficients[2])
summary(glm_out1)$coefficients[2,4]




###########################
# DO POWER SIMULATIONS ####
###########################


set.seed(6247555)

n=1000
effect_size=1
N=100
pow=1

sim_list<-list()

while (pow>0.8){
  values<-data.frame(pval=NA, beta_1=NA, power=NA) 
  n=n-10
  for(i in 1:N) {
    
    #simulate a larger sample size of predictors
    new_df<-data.frame(garden=sample(pilot$garden, size = n, replace = T),
                       species=sample(pilot$species, size = n, replace = T),
                       rater=sample(pilot$rater, size = n, replace = T))
    
    #generate a larger control set
    Y0 <- rpois(n/2, lambda=exp(b0 + b1*new_df$garden[1:n/2] + b2*new_df$species[1:n/2] + b3*new_df$rater[1:n/2]))
    
    #generate a condition set
    Y1 <- rpois(n/2, lambda=exp(b0 + b1*new_df$garden[n/2+1:nrow(new_df)] + b2*new_df$species[n/2+1:nrow(new_df)] + b3*new_df$rater[n/2+1:nrow(new_df)])+effect_size)
    
    combined<-new_df %>% 
      mutate(condition=as.factor(c(rep("water",n/2), rep("compound1",n/2))),
             tot.vase.days=c(Y0, Y1))
    combined$condition<-relevel(combined$condition, "water")
    
    glm_out1 <- glm(tot.vase.days ~ condition + garden + species + rater, 
                    family=poisson(link = "log"), data=combined)
    
    values[i,2]<-exp(glm_out1$coefficients[2]) #getting beta 1 (condition)
    values[i,1]<-summary(glm_out1)$coefficients[2,4] #getting pvalue 

  } 
  signif=values %>% filter(pval<0.05) 
  pow = dim(signif)[1]/dim(values)[1]
  values$power<-pow
  sim_list[[paste0("n",n)]]<-values
  
}

sim_df<-bind_rows(sim_list, .id="n")
ggplot(sim_df, aes(x=beta_1))+geom_density()


# missing:
# garden, rater as random effects
# multiple testing correction




ggplot(sim_df %>% filter(n=="n6"), aes(x=beta_1, y=n, color=pval<0.05))+geom_point()+xlab("Estimated effect size (beta 1)")+ylab("minimal sample size (6)")


