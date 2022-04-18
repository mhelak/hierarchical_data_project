
# checkup pilot data
# this was my first attempt, there is no balancing or random effects. Maybe with only balancing this could be already a good way (although of course in theory, the closer we get to the final model the better)

library(nlme)

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


#shoot, I only later on added the lines 8-10 where i recode the predictors as categorical. therefore the next code won't work (with predictors as continuous it did..) Will necessarily have to code it up similarly to the other script, i.e. f.e. in case of species1 and garden1 and rater1, i guess the lambda should be just the intercept? basically, do it in a balanced way and per combination of predictors.. I will do it the following days, don't have time right now.

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

n=400
effect_size=1
N=1000
pow=1

sim_list<-list()

while (pow>0.6){
  values<-data.frame(pval=NA, beta_1=NA, power=NA) 
  n=n-4
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
    values[i,1]<-summary(glm_out1)$coefficients[2,4]/2 #getting pvalue 

  } 
  signif=values %>% filter(pval<0.05) 
  pow = dim(signif)[1]/dim(values)[1]
  values$power<-pow
  sim_list[[paste0("n",n)]]<-values
  
}

sim_df<-bind_rows(sim_list, .id="n")
ggplot(sim_df, aes(x=beta_1))+geom_density()




sim_summarized<-sim_df %>%
  group_by(n) %>%
  summarise(mean_beta_1=mean(beta_1),
            power=power[1],
            n=n[1]) %>%
  mutate(n=as.integer(gsub("n","",n)))
  
ggplot(sim_summarized,aes(x=n, y=power))+geom_hline(yintercept = 0.85)+geom_point()

sim_summarized %>% filter(power>0.85) %>% slice_min(n) #332: 166 per group.






# missing:
# garden, rater as random effects
# multiple testing correction
# CI for the beta?



