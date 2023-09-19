### Bayesian  model for PI curve fitting
### Created on 9/7/2023
### Created by Nyssa Silbiger and Hollie Putnam
### https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-Source

### Load Libraries ##########
library(tidyverse)
library(tidybayes)
library(brms)
library(posterior)
library(rstan)
library(here)
library(ggthemes)
library(purrr)
library(reshape2)

##### Load data ###########
Data<-read_csv(here("output","pi_curve_extracted_rates.csv"))
length(unique(Data$colony_id))

#load metadata
md <- read_csv("data/1_pi_curves/coral_metadata.csv")
md <- md  %>%filter(species!="Blank")

All <- left_join(Data, md, by="colony_id")

#identify light and oxygen flux rates
Data$PAR <- as.numeric(Data$Light_Value.x)
Data$Pc <- as.numeric(Data$micromol.cm2.h)
Data_sub <- Data 
length(unique(Data$colony_id))

ggplot(Data_sub, aes(x = PAR, y=Pc))+
  geom_point()

#set priors
prior1 <- c(set_prior("normal(0, 5)", nlpar = "Am", lb = 0),
          set_prior("normal(0, 1)", nlpar = "AQY", lb = 0),
          set_prior("normal(0, 3)", nlpar = "Rd", lb = 0))

# #model 
# #Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd)
# 
# fit <- brm(bf(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), Am ~ 1, AQY ~ 1, Rd ~ 1, nl = TRUE), 
#               data = Data_sub_Dlab_A1, family = gaussian(),
#               prior = prior1,
#               chains = 4, iter = 2000, seed = 333)
# 
# #model summary
# summary(fit)
# 
# #plot model fit
# plot(fit)
# 
# #posterior predictive checks
# pp_check(fit)
# 
# #leave-one-out cross validation (LOO) method. The LOO assesses the predictive ability of posterior distributions 
# #(a little like the pp_check function). It is a good way to assess the fit of your model. 
# #You should look at the elpd estimate for each model, the higher value the better the fit. 
# #By adding compare = TRUE, we get a comparison already done for us at the bottom of the summary. 
# #The value with an elpd of 0 should appear, that’s the model that shows the best fit to our data.
# loo(fit, compare = TRUE)
# 
# model_fit <- Data_sub_Dlab %>%
#     add_predicted_draws(fit) %>%  # adding the posterior distribution
#     ggplot(aes(x = PAR, y = Pc)) +  
#     stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
#                     alpha = 0.5, colour = "black") +
#     geom_point(data = Data_sub_Dlab, colour = "darkseagreen4", size = 3) +   # raw data
#     scale_fill_brewer(palette = "Greys") +
#     ylab("Oxygen Flux") +  # latin name for red knot
#     xlab("PAR µmol m-2 s-1") +
#     theme_bw() +
#     theme(legend.title = element_blank(),
#           legend.position = c(0.15, 0.85))
# model_fit
# 
# #Am = Pmax
# #AQY = alpha = photochemical efficiency
# #Rd_Intercept = Rdark
# #NEED TO ADD THESE CALCULATIONS TO THE PARAMETER LIST
# #Ik=Am/AQY
# #Ic=(Am*Rd)/(AQY*(sqrt(Am^2-Rd^2))))
# fixef(fit)
# 
# 
# # Extract posterior values for each parameter
# samples1 <- posterior_samples(fit, "^b")
# head(samples1)
# 
# # get the predicted draws from the model
# pred_draws<-fit %>% 
#   epred_draws(newdata = expand_grid(PAR = seq(1,800, by = 100)), 
#               re_formula = NA)
# 
# # plot the fits
# ggplot(pred_draws, 
#        aes(x = PAR, y = .epred)) +
#   stat_lineribbon() +
#   geom_point(data =pred_draws, aes(x = PAR, y =.epred ) )+
#   scale_fill_brewer(palette = "Reds") +
#   labs(x = "PAR", y = "Rate",
#        fill = "Credible interval") +
#   theme_clean() +
#   theme(legend.position = "bottom")

#fit many models with a for loop

names <- unique(Data_sub$colony_id)
names 
fits <- setNames(vector("list", length(names)), names)

for (i in names) {
  data_model<-Data_sub %>% filter(colony_id==i)
  fits[[i]] <- fixef(brm(bf(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), Am ~ 1, AQY ~ 1, Rd ~ 1, nl = TRUE), 
                   data = data_model, family = gaussian(),
                   prior = prior1,
                   chains = 4, iter = 2000, seed = 333))
}
fits

Estimates <- as.data.frame(fits)
Estimates$parameter <- row.names(Estimates)
param.estimates <- melt(Estimates, id="parameter")
param.estimates <- param.estimates %>% separate(variable,  into=c("colony_id", "metric"), sep=7, remove = FALSE)
param.estimates <- param.estimates %>% separate(colony_id, into = c('species', 'colony_id'), sep=4, remove = FALSE)
param.estimates$colony_id <- paste0(param.estimates$species, param.estimates$colony_id)
param.estimates$colony_id <- gsub("\\.","-", param.estimates$colony_id)

#join parameters and metadata
Data.out <- left_join(param.estimates, md, by="colony_id")

#Plot all of the estimates
estimate.plots <-Data.out %>% filter(metric==".Estimate") %>%
  ggplot(aes(Temp.Cat, value, color=species.x)) +
  geom_point(size = 2) +
  #geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~parameter*species.x)

estimate.plots

