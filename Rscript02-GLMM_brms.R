################################################################################
# This code was written to perform GLMs for:
# Habitat associations and abundance estimates for two least-disturbed Redside 
# Dace (Clinostomus elongatus) populations in tributaries of Lake Huron 
# (Ontario, Canada). Canadian Journal of Fisheries and Aquatic Sciences
# Authors: Karl A. Lamothe, D. Andrew R. Drake
# Date revised: 10/1/2024; R Version 4.4.0. 
################################################################################
rm(list=ls())

# load libraries
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(patchwork)   # combine ggplots
p_load(brms)        # glms
p_load(posterior)   # tools for working with posterior and prior distributions
p_load(reshape2)    # manipulate tables
p_load(tidyverse)

# Personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=11, family="sans", colour="black"),
                  axis.text.x  = element_text(size=10, family="sans", colour="black"),
                  axis.text.y  = element_text(size=10, family="sans", colour="black"),
                  legend.title = element_text(size=10, family="sans", colour="black"),
                  legend.text  = element_text(size=10, family="sans", colour="black"),
                  plot.title   = element_text(size=11, family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black"),
                  legend.background = element_blank(),
                  legend.key = element_blank()))

# function to control decimal points in ggplots
scaleFUN <- function(x) sprintf("%.2f", x)

################################################################################
# read in abundance and habitat data
site.covs  <- read.csv("Site_Covariates.csv", header=T)
RSD.counts <- read.csv("Adult_Counts.csv",    header=T)

#  Separate adult counts
Mature.RSD.Detections <- RSD.counts[,c(4:6)]
Mature.RSD.Detections <- rowSums(Mature.RSD.Detections)

# average abundance across sites
sum(Mature.RSD.Detections)/43

# add in site abundance to site covariates
site.covs$RSD <- Mature.RSD.Detections

# quick plot to look at this monstrosity of a distribution
plot(site.covs$RSD, las=1, pch=20, 
     ylab="Counts", xlab="Site", main="Redside Dace")
plot(site.covs$RSD[site.covs$Waterbody=='Gully Creek'], 
     las=1, pch=20, ylab="Counts", xlab="Site", col="red")
ggplot(site.covs, aes(x=RSD, fill=Waterbody))+
  geom_histogram(color='white')

################################################################################
# Consider habitat covariates
################################################################################
colnames(site.covs)

# scaled and centered to allow comparisons of effects between variables
site.covs2<-data.frame(Depth.max     = scale(site.covs$Depth_max, center=T),
                       HH            = scale(site.covs$HH_mean,   center=T),
                       Temp          = scale(site.covs$Wtemp,     center=T),
                       P.vol         = scale(site.covs$Pool_vol,  center=T),
                       Co.occurring  = scale(site.covs$Co.occurring,  center=T),
                       Woody         = scale(site.covs$Woody_Debris,  center=T),
                       UnderCut.prop = scale(site.covs$UnderCut.Prop, center=T),
                       Nesters       = scale(site.covs$Nesters, center=T),
                       Waterbody     = site.covs$Waterbody,
                       Year          = site.covs$Year,
                       RSD           = Mature.RSD.Detections,
                       RelRSD        = Mature.RSD.Detections/site.covs$Pool_area,
                       ID            = site.covs$Field_Number)

# Convert to factors
# Waterbody
site.covs2$wateryear <- paste0(site.covs2$Waterbody, " ", site.covs2$Year)

# Gully Creek only
site.covs3 <- site.covs2[site.covs2$Waterbody=='Gully Creek',]                  
site.covs3$Year[site.covs3$Year=="2019"]<-0
site.covs3$Year[site.covs3$Year=="2020"]<-1
site.covs3$Year<-as.factor(site.covs3$Year)

# full data
site.covs2$Waterbody[site.covs2$Waterbody=="Gully Creek"]<-1
site.covs2$Waterbody[site.covs2$Waterbody=="Unknown Stan J"]<-0
site.covs2$Waterbody<-as.factor(site.covs2$Waterbody)

# Year
site.covs2$Year[site.covs2$Year=="2019"]<-0
site.covs2$Year[site.covs2$Year=="2020"]<-1
site.covs2$Year<-as.factor(site.covs2$Year)

# Site ID
site.covs2$ID<-as.factor(site.covs2$ID)

site.covs$densm2 <- site.covs$RSD/site.covs$Pool_area
site.covs$wateryear <- paste0(site.covs$Waterbody, " ", site.covs$Year)

# summaries
aggregate(site.covs2$RSD, list(site.covs2$wateryear), median)
aggregate(site.covs2$RSD, list(site.covs2$wateryear), sum)
aggregate(site.covs2$RSD, list(site.covs2$wateryear), min)
aggregate(site.covs2$RSD, list(site.covs2$wateryear), max)

################################################################################
# Fit GLMs
################################################################################
mod.form.1  <- RSD ~ 1
mod.form.2  <- RSD ~ Depth.max
mod.form.3  <- RSD ~ Depth.max + Temp
mod.form.4  <- RSD ~ Depth.max + UnderCut.prop
mod.form.5  <- RSD ~ Depth.max + Nesters
mod.form.6  <- RSD ~ Depth.max + wateryear
mod.form.7  <- RSD ~ Temp
mod.form.8  <- RSD ~ Temp + UnderCut.prop
mod.form.9  <- RSD ~ Temp + Nesters
mod.form.10 <- RSD ~ Temp + wateryear
mod.form.11 <- RSD ~ UnderCut.prop
mod.form.12 <- RSD ~ UnderCut.prop + Nesters
mod.form.13 <- RSD ~ UnderCut.prop + wateryear
mod.form.14 <- RSD ~ Nesters
mod.form.15 <- RSD ~ Nesters + wateryear
mod.form.16 <- RSD ~ wateryear

n.chains = 3; n.iter = 50000; n.burn = 10000; n.thin = 20

# identify prior for intercept model
mod.prior1 <- c(prior(normal(0, 2.72),   class = Intercept),
                prior(gamma(0.01, 0.01), class= shape))

# run model
GLM.mod1 <- brm(mod.form.1, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)

# identify priors for rest of models
mod.prior <- c(prior(normal(0, 2.72),  class = Intercept),
               prior(normal(0, 2.72),  class = b),
               prior(gamma(0.01, 0.01), class= shape))

## run models. split up in batches of 5 when running to avoid crashing
GLM.mod2 <- brm(mod.form.2, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod3 <- brm(mod.form.3, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod4 <- brm(mod.form.4, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod5 <- brm(mod.form.5, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod6 <- brm(mod.form.6, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)

# batch
GLM.mod7 <- brm(mod.form.7, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod8 <- brm(mod.form.8, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod9 <- brm(mod.form.9, family = "negbinomial", data = site.covs2,
                chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod10 <- brm(mod.form.10, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod11 <- brm(mod.form.11, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)

# batch
GLM.mod12 <- brm(mod.form.12, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod13 <- brm(mod.form.13, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod14 <- brm(mod.form.14, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod15 <- brm(mod.form.15, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)
GLM.mod16 <- brm(mod.form.16, family = "negbinomial", data = site.covs2,
                 chains = n.chains, iter = n.iter, warmup = n.burn, thin = n.thin,
                 save_pars = save_pars(all=TRUE), prior = mod.prior)

################################################################################
# traceplots
################################################################################

plot(GLM.mod1)
plot(GLM.mod2)
plot(GLM.mod3)
plot(GLM.mod4)
plot(GLM.mod5)
plot(GLM.mod6)
plot(GLM.mod7)
plot(GLM.mod8)
plot(GLM.mod9)
plot(GLM.mod10)
plot(GLM.mod11)
plot(GLM.mod12)
plot(GLM.mod13)
plot(GLM.mod14)
plot(GLM.mod15)
plot(GLM.mod16)

################################################################################
# Look at model diagnostics and estimates
################################################################################

#launch_shinystan(GLM.mod1)
#launch_shinystan(GLM.mod2)
#launch_shinystan(GLM.mod3)
#launch_shinystan(GLM.mod4)
#launch_shinystan(GLM.mod5)
#launch_shinystan(GLM.mod6)
#launch_shinystan(GLM.mod7)
#launch_shinystan(GLM.mod8)
#launch_shinystan(GLM.mod9)
#launch_shinystan(GLM.mod10)
#launch_shinystan(GLM.mod11)
#launch_shinystan(GLM.mod12)
#launch_shinystan(GLM.mod13)
#launch_shinystan(GLM.mod14)
#launch_shinystan(GLM.mod15)
#launch_shinystan(GLM.mod16)

################################################################################
# add loo-CV for each model for model comparisons
################################################################################
GLM.mod1  <- add_criterion(GLM.mod1,  criterion = c("loo"), moment_match=T)
GLM.mod2  <- add_criterion(GLM.mod2,  criterion = c("loo"), moment_match=T)
GLM.mod3  <- add_criterion(GLM.mod3,  criterion = c("loo"), moment_match=T)
GLM.mod4  <- add_criterion(GLM.mod4,  criterion = c("loo"), moment_match=T)
GLM.mod5  <- add_criterion(GLM.mod5,  criterion = c("loo"), moment_match=T)
GLM.mod6  <- add_criterion(GLM.mod6,  criterion = c("loo"), moment_match=T)
GLM.mod7  <- add_criterion(GLM.mod7,  criterion = c("loo"), moment_match=T)
GLM.mod8  <- add_criterion(GLM.mod8,  criterion = c("loo"), moment_match=T)
GLM.mod9  <- add_criterion(GLM.mod9,  criterion = c("loo"), moment_match=T)
GLM.mod10 <- add_criterion(GLM.mod10, criterion = c("loo"), moment_match=T)
GLM.mod11 <- add_criterion(GLM.mod11, criterion = c("loo"), moment_match=T)
GLM.mod12 <- add_criterion(GLM.mod12, criterion = c("loo"), moment_match=T)
GLM.mod13 <- add_criterion(GLM.mod13, criterion = c("loo"), moment_match=T)
GLM.mod14 <- add_criterion(GLM.mod14, criterion = c("loo"), moment_match=T)
GLM.mod15 <- add_criterion(GLM.mod15, criterion = c("loo"), moment_match=T)
GLM.mod16 <- add_criterion(GLM.mod16, criterion = c("loo"), moment_match=T)

# model comparisons LOO
comp.loo<-loo_compare(GLM.mod1,  GLM.mod2,  GLM.mod3,  GLM.mod4,  GLM.mod5,  
                      GLM.mod6,  GLM.mod7,  GLM.mod8,  GLM.mod9,  GLM.mod10, 
                      GLM.mod11, GLM.mod12, GLM.mod13, GLM.mod14, GLM.mod15, 
                      GLM.mod16)
print(comp.loo, digits = 3, simplify=FALSE)
comp.loo.df <- cbind.data.frame(comp.loo)

################################################################################
# Bayesian R2 for each model
################################################################################
bayes_R2_res <- function(fit) {
  y <- rstanarm::get_y(fit)
  ypred <- rstanarm::posterior_epred(fit)
  if (family(fit)$family == "binomial" && NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1]
    ypred <- ypred %*% diag(trials)
  }
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}

round(median(bayesR2res<-bayes_R2_res(GLM.mod1)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod2)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod3)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod4)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod5)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod6)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod7)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod8)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod9)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod10)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod11)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod12)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod13)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod14)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod15)), 2)
round(median(bayesR2res<-bayes_R2_res(GLM.mod16)), 2)

################################################################################
## R hat
################################################################################
rhat.df<-cbind.data.frame(Rhat = c(brms::rhat(GLM.mod1),  brms::rhat(GLM.mod2),  brms::rhat(GLM.mod3),
                                   brms::rhat(GLM.mod4),  brms::rhat(GLM.mod5),  brms::rhat(GLM.mod6),
                                   brms::rhat(GLM.mod7),  brms::rhat(GLM.mod8),  brms::rhat(GLM.mod9),
                                   brms::rhat(GLM.mod10), brms::rhat(GLM.mod11), brms::rhat(GLM.mod12),
                                   brms::rhat(GLM.mod13), brms::rhat(GLM.mod14), brms::rhat(GLM.mod15),
                                   brms::rhat(GLM.mod16)),
                          Model = c(rep("1",  length(brms::rhat(GLM.mod1))),  rep("2",  length(brms::rhat(GLM.mod2))),
                                    rep("3",  length(brms::rhat(GLM.mod3))),  rep("4",  length(brms::rhat(GLM.mod4))),
                                    rep("5",  length(brms::rhat(GLM.mod5))),  rep("6",  length(brms::rhat(GLM.mod6))),
                                    rep("7",  length(brms::rhat(GLM.mod7))),  rep("8",  length(brms::rhat(GLM.mod8))),
                                    rep("9",  length(brms::rhat(GLM.mod9))),  rep("10", length(brms::rhat(GLM.mod10))),
                                    rep("11", length(brms::rhat(GLM.mod11))), rep("12", length(brms::rhat(GLM.mod12))),
                                    rep("13", length(brms::rhat(GLM.mod13))), rep("14", length(brms::rhat(GLM.mod14))),
                                    rep("15", length(brms::rhat(GLM.mod15))), rep("16", length(brms::rhat(GLM.mod16)))))
rhat.df$Model <- factor(rhat.df$Model, levels=c("1","2","3","4","5","6","7","8",
                                                "9","10","11","12","13",'14',
                                                '15','16'))

# plot
ggplot(rhat.df, aes(y=Rhat, x=Model))+
  geom_point() +
  coord_flip() +
  labs(y = expression(''*italic(widehat(R))))

# perform contrasts between wateryear across models that include the variable
q2 <- c(q2 = "wateryearGullyCreek2020 + wateryearUnknownStanJ2020 = 0")
plot(hypothesis(GLM.mod6, q2))
plot(hypothesis(GLM.mod10, q2))
plot(hypothesis(GLM.mod13, q2))
plot(hypothesis(GLM.mod15, q2))
plot(hypothesis(GLM.mod16, q2))

################################################################################
# extract model coefficients and plot
################################################################################

GLM.fixedef.res <- rbind.data.frame(fixef(GLM.mod1),  fixef(GLM.mod2),  fixef(GLM.mod3),
                                    fixef(GLM.mod4),  fixef(GLM.mod5),  fixef(GLM.mod6),
                                    fixef(GLM.mod7),  fixef(GLM.mod8),  fixef(GLM.mod9),
                                    fixef(GLM.mod10), fixef(GLM.mod11), fixef(GLM.mod12),
                                    fixef(GLM.mod13), fixef(GLM.mod14), fixef(GLM.mod15),
                                    fixef(GLM.mod16))
GLM.fixedef.res$ModNo <- c(1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,7,7,8,8,8,
                           9,9,9,10,10,10,10,11,11,12,12,12,13,13,13,13,
                           14,14,15,15,15,15,16,16,16)
GLM.fixedef.res$Variables <- c('Intercept',
                               'Intercept','Max depth',
                               'Intercept','Max depth','Water temp',
                               'Intercept','Max depth','Undercut banks',
                               'Intercept','Max depth','Nest-builders',
                               'Intercept','Max depth','Gully Creek 2020','USJ Tributary 2020',
                               'Intercept','Water temp',
                               'Intercept','Water temp','Undercut banks',
                               'Intercept','Water temp','Nest-builders',
                               'Intercept','Water temp','Gully Creek 2020','USJ Tributary 2020',
                               'Intercept','Undercut banks',
                               'Intercept','Undercut banks','Nest-builders',
                               'Intercept','Undercut banks','Gully Creek 2020','USJ Tributary 2020',
                               'Intercept','Nest-builders',
                               'Intercept','Nest-builders','Gully Creek 2020','USJ Tributary 2020',
                               'Intercept','Gully Creek 2020','USJ Tributary 2020')
GLM.fixedef.res$ModNo <- as.character(GLM.fixedef.res$ModNo)

# remove intercept values for visual purposes
GLM.fixedef.res.plt   <- GLM.fixedef.res[!GLM.fixedef.res$Variables=="Intercept",]

# order the variables for visual purposes
GLM.fixedef.res.plt$Variables <- factor(GLM.fixedef.res.plt$Variables,
                                        levels = c('USJ Tributary 2020', 'Gully Creek 2020','Nest-builders',
                                                   'Undercut banks', 'Water temp', 'Max depth'))

modcoef.gg<-ggplot(GLM.fixedef.res.plt, aes(y=Estimate, x=Variables, group=ModNo))+
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), size=0.25, 
                  position = position_dodge2(width=0.7, preserve="single"),
                  show.legend = FALSE, linewidth=0.25)+
  geom_point(position = position_dodge2(width=0.7, preserve="single"), size=0.25)+
  coord_flip()+
  labs(y = "Posterior estimate", color = "Model") +
  guides(color=guide_legend(ncol=2))+
  theme(axis.title.y = element_blank(),
        legend.position = 'none') 

########################################
# compute model-averaged fitted values #
########################################
# this code makes it easier to get the raw weights
stacks <- model_weights(GLM.mod1,  GLM.mod2,  GLM.mod3,  GLM.mod4,  GLM.mod5,  
                        GLM.mod6,  GLM.mod7,  GLM.mod8,  GLM.mod9,  GLM.mod10, 
                        GLM.mod11, GLM.mod12, GLM.mod13, GLM.mod14, GLM.mod15, 
                        GLM.mod16,
                        weights = "stacking", # this is loo2
                        method = 'fitted',
                        newdata = site.covs2,
                        seed=09812)
stacks2 <- cbind.data.frame(stacks)

# this code provides site-level predictions
pp_a<-pp_average(GLM.mod1,  GLM.mod2,  GLM.mod3,  GLM.mod4,  GLM.mod5,  
                 GLM.mod6,  GLM.mod7,  GLM.mod8,  GLM.mod9,  GLM.mod10, 
                 GLM.mod11, GLM.mod12, GLM.mod13, GLM.mod14, GLM.mod15, 
                 GLM.mod16,
                 weights = "stacking", # this is loo2
                 method  = 'fitted',
                 #newdata = site.covs2,
                 summary = F,
                 seed    = 09812)

# plot model weights
mw.manual <- cbind.data.frame(model =  seq(1:16),
                              weight = round(stacks2$stacks,2),
                              model2 = c('Intercept',
                                         'Max depth',
                                         'Max depth + Water temp',
                                         'Max depth + Undercut banks',
                                         'Max depth + Nest-builders',
                                         'Max depth + Period',
                                         'Water temp',
                                         'Water temp + Undercut banks',
                                         'Water temp + Nest-builders',
                                         'Water temp + Period',
                                         'Undercut banks',
                                         'Undercut banks + Nest-builders',
                                         'Undercut banks + Period',
                                         'Nest-builders',
                                         'Nest-builders + Period',
                                         'Period'))

# Bayesian stacking with greater than 0.05 weight
mw.manual[mw.manual$weight > 0.05,]

modweight.gg<-ggplot(mw.manual, aes(x=reorder(model2,weight), y=weight))+
  geom_point(size=1.5) +
  labs(y="Model weight", x="Model") +
  coord_flip() +
  ylim(0,0.5) +
  theme(axis.title.y=element_blank())

modcoef.gg+modweight.gg+
  plot_annotation(tag_levels='a', tag_suffix = ")")

################################################################################
################################################################################
# model averaged site predictions
################################################################################
################################################################################

sitepreds <- cbind.data.frame(post   = c(pp_a),
                              site   = rep(seq(1:43), each = 6000),
                              Area   = rep(site.covs$Pool_area, each = 6000),
                              Volume = rep(site.covs$Pool_vol,  each = 6000),
                              Length = rep(site.covs$P_Length, each = 6000),
                              Strahler = rep(site.covs$Strahler,  each = 6000),
                              Waterbody = rep(c(rep('Gully Creek', 16),  
                                                rep('USJ Tributary', 9),
                                                rep('Gully Creek', 18)), each = 6000),
                              Year = rep(c(rep('2019', 16), 
                                           rep('2020', 9),
                                           rep('2020', 18)), each = 6000),
                              Period = rep(c(rep('Gully Creek 2019', 16), 
                                             rep('USJ Tributary 2020', 9),
                                             rep('Gully Creek 2020', 18)), each = 6000))
sitepreds$Strahler2 <- sitepreds$Strahler

sitepreds$Strahler2[sitepreds$Waterbody == 'Gully Creek' & sitepreds$Strahler2 == 2] <- 15.19
sitepreds$Strahler2[sitepreds$Waterbody == 'Gully Creek' & sitepreds$Strahler2 == 3] <- 16.75
sitepreds$Strahler2[sitepreds$Waterbody == 'Gully Creek' & sitepreds$Strahler2 == 4] <- 33.54
sitepreds$Strahler2[sitepreds$Waterbody == 'USJ Tributary' & sitepreds$Strahler2 == 2] <- 15.19
sitepreds$Strahler2[sitepreds$Waterbody == 'USJ Tributary' & sitepreds$Strahler2 == 3] <- 13.15
sitepreds$Strahler2[sitepreds$Waterbody == 'USJ Tributary' & sitepreds$Strahler2 == 4] <- 20.74

# divide the site-predicted abundance by length and area
sitepreds$site      <- as.character(sitepreds$site)
sitepreds$Strahler  <- as.character(sitepreds$Strahler)
sitepreds$fishm2    <- sitepreds$post/sitepreds$Area
sitepreds$fishm     <- sitepreds$post/sitepreds$Length
sitepreds$fishreach <- sitepreds$post/sitepreds$Strahler2

fish.m.df <- cbind.data.frame(Fish = c(sitepreds$fishm, sitepreds$fishm2),
                              Variable = rep(c("Redside Dace per m",'Redside Dace per m^2'), each = length(sitepreds$post)),
                              Strahler = c(sitepreds$Strahler, sitepreds$Strahler),
                              Waterbody= c(sitepreds$Waterbody, sitepreds$Waterbody))

ggplot()+
  geom_density(data = fish.m.df, aes(x=log(Fish), fill=Strahler), alpha=0.5)+
  labs(x='log(Abundance)', y='Probability Density Function')+
  scale_x_continuous(labels=scaleFUN)+
  scale_fill_manual(values=c("gold","grey","pink"))+
  facet_grid(Variable~Waterbody)+
  guides(fill=guide_legend(position='inside')) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.4,0.75))

# only sample estimates from gully sites - this excludes USJ tributary sites
sitepreds.g <- sitepreds[!c(sitepreds$site =='17' | sitepreds$site =='18' | 
                              sitepreds$site =='19' | sitepreds$site =='20' | 
                              sitepreds$site =='21' | sitepreds$site =='22' | 
                              sitepreds$site =='23' | sitepreds$site =='24' | 
                              sitepreds$site =='25' ),]

# only sample estimates from USJ Tributary
sitepreds.usj<-sitepreds[c(sitepreds$site =='17' | sitepreds$site =='18' | 
                             sitepreds$site =='19' | sitepreds$site =='20' | 
                             sitepreds$site =='21' | sitepreds$site =='22' | 
                             sitepreds$site =='23' | sitepreds$site =='24' | 
                             sitepreds$site =='25' ),]

G2 <- sitepreds.g$post[sitepreds.g$Strahler=='2']
G3 <- sitepreds.g$post[sitepreds.g$Strahler=='3']
G4 <- sitepreds.g$post[sitepreds.g$Strahler=='4']
U3 <- sitepreds.usj$post[sitepreds.usj$Strahler=='3']
U4 <- sitepreds.usj$post[sitepreds.usj$Strahler=='4']

RSD.km.df <- cbind.data.frame(Posterior = c(G2,G3,G4,U3,U4),
                              Waterbody = c(rep("Gully Creek", 204000),
                                            rep("USJ Tributary",54000)),
                              Strahler = c(rep("2",60000), rep("3",48000),
                                           rep("4",96000), rep("3", 6000),
                                           rep("4", 48000)),
                              RSD.km = c((686/10.414)*G2,
                                         (125/2.091)*G3,
                                         (191/6.396)*G4,
                                         (173/2.279)*U3,
                                         (261/5.407)*U4))


fish.m.df1 <- cbind.data.frame(Fish = c(fish.m.df$Fish, RSD.km.df$RSD.km),
                               Waterbody = c(fish.m.df$Waterbody, RSD.km.df$Waterbody),
                               Strahler = c(fish.m.df$Strahler, RSD.km.df$Strahler),
                               Variable = c(fish.m.df$Variable, rep("Redside Dace per km", length(RSD.km.df$Posterior))))

meds<-aggregate(fish.m.df1$Fish, list(fish.m.df1$Variable, 
                                      fish.m.df1$Strahler, 
                                      fish.m.df1$Waterbody), median)
meds2<-aggregate(fish.m.df1$Fish, list(fish.m.df1$Variable, 
                                       fish.m.df1$Strahler, 
                                       fish.m.df1$Waterbody), quantile, probs=c(0.025, 0.975))
round(meds[4], 2)

fish.per.gg<-ggplot()+
  geom_density(data = fish.m.df1, aes(x=log(Fish), fill=Strahler), alpha=0.5, lwd=0.25)+
  labs(x='log(Abundance)', y='Probability Density Function')+
  scale_fill_manual(values=c("gold","grey","pink"))+
  facet_grid(Waterbody~Variable, scales = 'free_y')+
  guides(fill=guide_legend(position='inside')) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.07,0.25))
fish.per.gg

# site level abundance estimates by waterbody
aggregate(sitepreds$post, list(sitepreds$Waterbody), median)
aggregate(sitepreds$post, list(sitepreds$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/m of pool habitat 
aggregate(sitepreds$fishm, list(sitepreds$Waterbody), median) 
aggregate(sitepreds$fishm, list(sitepreds$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/m2 of pool habitat per waterbody
aggregate(sitepreds$fishm2, list(sitepreds$Waterbody), median)
aggregate(sitepreds$fishm2, list(sitepreds$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/m of reach habitat per waterbody
aggregate(sitepreds$fishreach, list(sitepreds$Strahler, sitepreds$Waterbody), median)
aggregate(sitepreds$fishreach, list(sitepreds$Waterbody), median)

# number of fish per kilometer = median number of fish per pool * number of pools per km
# waterbody level
aggregate(sitepreds$post, list(sitepreds$Strahler, sitepreds$Waterbody), median)

Period.abundgg<-ggplot()+
  geom_density(data = sitepreds, aes(x=log(post), fill=Period), alpha=0.5)+
  labs(x='log(Abundance)', y='Probability Density Function')+
  scale_x_continuous(labels=scaleFUN, limits=c(-3,12))+
  guides(fill=guide_legend(position='inside')) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.75,0.75),
        legend.title = element_blank())
Period.abundgg

#===============================================================================
#===============================================================================
# Sample the posterior distribution of site abundance to generate estimates of 
# the abundance of RSD across the entire system
# In this iteration separate posterior distributions are used based on the 
# waterbody and straher order are being considered. 
# Gully Creek 2019 and 2020 posteriors are combined
#===============================================================================
#===============================================================================

reps <- 20000 # number of replicates
G.pools2 <-0
G.pools3 <-0
G.pools4 <-0
U.pools3 <-0
U.pools4 <-0
Abund.Gull.2 <-0
Abund.Gull.3 <-0
Abund.Gull.4 <-0
Abund.Stan.3 <-0
Abund.Stan.4 <-0

###############
# ~~~~~~~~~~~ #
# ~~~~~~~~~~~ #
###############
# this takes a couple minutes to run
set.seed(0528)
for(i in 1:reps){
  # number of pools
  # Gully Creek
  G.pools2[i]  <- floor(runif(1, min = 482, max = 1186)) #Strahler 2
  G.pools3[i]  <- floor(runif(1, min = 88,  max = 216))  #Strahler 3
  G.pools4[i]  <- floor(runif(1, min = 143, max = 284))  #Strahler 4
  
  # USJ Tributary
  U.pools3[i] <- floor(runif(1, min = 113, max = 374)) #Strahler 3
  U.pools4[i] <- floor(runif(1, min = 207, max = 351)) #Strahler 4 
  
  # Calculate site abundance
  # Gully Creek
  Abund.Gull.2[i] <- sample(sitepreds.g$post[sitepreds.g$Strahler=='2'], size=1) * G.pools2[i]
  Abund.Gull.3[i] <- sample(sitepreds.g$post[sitepreds.g$Strahler=='3'], size=1) * G.pools3[i]
  Abund.Gull.4[i] <- sample(sitepreds.g$post[sitepreds.g$Strahler=='4'], size=1) * G.pools4[i]
  
  # USJ Tributary
  Abund.Stan.3[i] <- sample(sitepreds.usj$post[sitepreds.usj$Strahler=='3'], size=1) * U.pools3[i]
  Abund.Stan.4[i] <- sample(sitepreds.usj$post[sitepreds.usj$Strahler=='4'], size=1) * U.pools4[i]
}

# Combine data into a df for summarizing
Abund.post2 <- cbind.data.frame(Strahler = rep(c(2,3,4,3,4), each=reps),
                                Site.abund = c(Abund.Gull.2, Abund.Gull.3, Abund.Gull.4,
                                               Abund.Stan.3, Abund.Stan.4),
                                Waterbody  = c(rep("Gully Creek",   reps*3), 
                                               rep("USJ Tributary", reps*2)),
                                Pools      = c(G.pools2, G.pools3, G.pools4,
                                               U.pools3, U.pools4))
Abund.post2$Strahler <- as.character(Abund.post2$Strahler)
Abund.post2$fish.pool <- Abund.post2$Site.abund/Abund.post2$Pools

Diff.post3<-ggplot()+
  geom_density(data = Abund.post2, aes(x=log(Site.abund), fill=Waterbody), alpha=0.5) +
  labs(x='log(Abundance)', y='Probability Density Function') +
  scale_fill_manual(values=c("pink","darkblue"))+
  guides(fill=guide_legend(position='inside'))+
  geom_vline(xintercept = log(7791), linetype='dashed')+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.79,0.75),
        legend.title = element_blank())
Diff.post3

# plot together
Period.abundgg + Diff.post3+
  plot_annotation(tag_levels='a', tag_suffix = ")")

# Index of abundance 
aggregate(Abund.post2$Site.abund, list(Abund.post2$Waterbody), median)
aggregate(Abund.post2$Site.abund, list(Abund.post2$Waterbody), "quantile", probs=c(0.05, 0.95))

# probability greater than MVP
sum(Abund.post2$Site.abund[Abund.post2$Waterbody=="USJ Tributary"] >= 7791) / length(Abund.post2$Site.abund[Abund.post2$Waterbody=="USJ Tributary"])
sum(Abund.post2$Site.abund[Abund.post2$Waterbody=="Gully Creek"] >= 7791) / length(Abund.post2$Site.abund[Abund.post2$Waterbody=="Gully Creek"])

# plot
ggplot()+
  geom_density(data = Abund.post2, aes(x=log(fish.pool), fill=Strahler), adjust=2, alpha=0.5) +
  labs(x='log(Abundance)', y='Probability Density Function') +
  guides(fill=guide_legend(position='inside'))+
  facet_wrap(~Waterbody)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.4,0.75))

