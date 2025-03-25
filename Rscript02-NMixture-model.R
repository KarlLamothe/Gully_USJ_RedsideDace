#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# N-mixture model for Redside Dace in USJ Tributary and Gully Creek. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(pacman)    # download and load packages
p_load(rjags)      # Bayesian model
p_load(ggplot2)    # plotting
p_load(jagsUI)     # Bayesian model
p_load(patchwork)  # plotting

# personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=10, family="sans", colour="black"),
                  axis.text.x  = element_text(size=10, family="sans", colour="black"),
                  axis.text.y  = element_text(size=10, family="sans", colour="black"),
                  legend.title = element_text(size=10, family="sans", colour="black"),
                  legend.text  = element_text(size=10, family="sans", colour="black"),
                  plot.title   = element_text(size=10, family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black"),
                  legend.background = element_blank(),
                  legend.key = element_blank()))

################################################################################
# read in habitat data and RSD mature data
################################################################################
site.covs  <- read.csv("Site_covariates.csv", header=T)
RSD.counts <- read.csv("Adult_Counts.csv",    header=T)

#  Separate adult counts
Mature.RSD <- RSD.counts[,c(4:6)]
Mature.RSD.Detections <- rowSums(Mature.RSD)

# add in site abundance to site covariates
site.covs$RSD <- Mature.RSD.Detections
site.covs$RSD.01 <- site.covs$RSD # presence absence
site.covs$RSD.01[site.covs$RSD.01>0] <-1

# factors
site.covs$Strahler <- as.character(site.covs$Strahler)

# scaled and centered to allow comparisons of effects between variables
site.covs2<-data.frame(Depth.max     = scale(site.covs$Depth_max, center=T),
                       HH            = scale(site.covs$HH_mean,   center=T),
                       Temp          = scale(site.covs$Wtemp,     center=T),
                       Area          = scale(site.covs$Pool_area, center=T),
                       Woody         = scale(site.covs$Woody_Debris,  center=T),
                       UnderCut.prop = scale(site.covs$UnderCut.Prop, center=T),
                       Nesters       = scale(site.covs$Nesters, center=T),
                       Waterbody     = site.covs$Waterbody,
                       Year          = site.covs$Year,
                       RSD           = Mature.RSD.Detections,
                       ID            = site.covs$Field_Number)

# Site ID
site.covs2$ID<-as.factor(site.covs2$ID)

# Redside Dace per m2
site.covs$densm2 <- site.covs$RSD/site.covs$Pool_area

# Convert to dummy variables
site.covs2$Event    <- paste0(site.covs2$Waterbody, " ", site.covs2$Year)
site.covs2$GC_2020  <- ifelse(site.covs2$Event == "Gully Creek 2020", 1, 0)
site.covs2$USJ_2020 <- ifelse(site.covs2$Event == "Unknown Stan J 2020", 1, 0)

################################################################################
# Prepare modelling
################################################################################
#Define array dimensions
nsite <- length(Mature.RSD$Haul1) # Number of sites
nrep <- 3                         # Number of depletion passes per site

################################################################################
################################################################################
# Write model
# Hashtagged out are the variables that were not included in the final model
# and the inclusion parameters that were used to determine the final model.
################################################################################
################################################################################
modelFilename <- "Poisson.model.txt"
cat('model{
       # Priors for abundance model
       beta0 ~ dnorm(0, 0.01)   # Prior for intercept
       beta1 ~ dnorm(0, 0.01)   # Prior for slope of GC_2020 dummy
       beta2 ~ dnorm(0, 0.01)   # Prior for slope of USJ_2020 dummy
       #beta3 ~ dnorm(0, 0.01)   # Prior for slope of Nesters (not included)
       #beta4 ~ dnorm(0, 0.01)   # Prior for slope of maximum depth (not included)
       #beta5 ~ dnorm(0, 0.01)   # Prior for slope of undercut banks (not included)
       #beta6 ~ dnorm(0, 0.01)   # Prior for slope of water temp (not included)
       #beta7 ~ dnorm(0, 0.01)   # Prior for slope of hydraulic head (not included)
       #beta8 ~ dnorm(0, 0.01)   # Prior for slope of pool area (not included)
       phi   ~ dunif(0.01, 100) # Prior for overdispersion parameter (site-specific variation)

       # Inclusion parameters for abundance model
       #z.beta1 ~ dbern(0.5)
       #z.beta2 ~ dbern(0.5)
       #z.beta3 ~ dbern(0.5)
       #z.beta4 ~ dbern(0.5)
       #z.beta5 ~ dbern(0.5)
       #z.beta6 ~ dbern(0.5)
       #z.beta7 ~ dbern(0.5)
       #z.beta8 ~ dbern(0.5)

       # Priors for detection model
       for (i in 1:nsite) {
         q0[i] ~ dunif(0, 1) # initial detection probability
         a[i]  ~ dunif(0, 1) # detection pr decay parameter
       }
       alpha0 ~ dnorm(0, 0.01)    # Prior for intercept
       #alpha1 ~ dnorm(0, 0.01)    # Prior for HH (not included)
       #alpha2 ~ dnorm(0, 0.01)    # Prior for depth (not included)
       
       # Inclusion parameters for detection model
       #z.alpha1 ~ dbern(0.5)
       #z.alpha2 ~ dbern(0.5)
       
       # Likelihood
       for (i in 1:nsite) {
         eta[i] ~ dgamma(phi, phi) # Prior for Gamma latent variable (extra-site specific variation)
         
          # Ecological model for true abundance
         N[i,1]    <- N.total[i]
         
         #####################################
         # Poisson model
         #####################################
         N.total[i] ~ dpois(lambda[i])
         lambda[i]  <- mu[i] * eta[i]
         log(mu[i]) <- beta0 + beta2 * USJ_2020[i] + beta1 * GC_2020[i]
                       #z.beta1 * beta1 * GC_2020[i] +
                       #z.beta2 * beta2 * USJ_2020[i] +
                       #z.beta3 * beta3 * Nesters[i] +
                       #z.beta4 * beta4 * Depth.max[i] +
                       #z.beta5 * beta5 * UnderCut.prop[i] +
                       #z.beta6 * beta6 * Temp[i] +
                       #z.beta7 * beta7 * HH[i] +                       
                       #z.beta8 * beta8 * Area[i] 

         
         for (j in 1:nrep) {
           # Observation model for removal count data
           counts.multi[i,j] ~ dbin(q[i,j], N[i,j])
           N[i,j + 1] <- N[i,j] - counts.multi[i,j]
         }
       }
       
       # Detection model
       for (i in 1:nsite) {
         for (j in 1:nrep) {
           q[i, j] <- q1[i] + (q0[i] - q1[i]) * (1 - pow(a[i], (j - 1))) 
         }
         
         # intercept model because no variables were supported
         logit(q1[i]) <- alpha0 # + 
                         #z.alpha1 * alpha1 * HH[i] +
                         #z.alpha2 * alpha2 * Depth.max[i]
       }
     }', fill = TRUE, file = modelFilename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Prepare model run
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Initial abundance values
Nst <- Mature.RSD.Detections + 1

# initial values
jags.inits <- function() {
  list(
    N.total = Nst, 
    alpha0  = rnorm(1, 0, 1), 
    #alpha1  = rnorm(1, 0, 1), alpha2 = rnorm(1, 0, 1), 
    beta0   = rnorm(1, 0, 1),
    beta1   = rnorm(1, 0, 1),
    beta2   = rnorm(1, 0, 1),
    #beta3   = rnorm(1, 0, 1),beta4 = rnorm(1, 0, 1),
    #beta5   = rnorm(1, 0, 1),beta6 = rnorm(1, 0, 1),
    #beta7   = rnorm(1, 0, 1),beta8 = rnorm(1, 0, 1),
    phi     = runif(1, 0.01, 100),
    q0      = runif(nsite, 0, 1),
    a       = runif(nsite, 0, 1)
    #z.beta1 = rbinom(1, 1, 0.5), z.beta2 = rbinom(1, 1, 0.5),
    #z.beta3 = rbinom(1, 1, 0.5), z.beta4 = rbinom(1, 1, 0.5),
    #z.beta5 = rbinom(1, 1, 0.5), z.beta6 = rbinom(1, 1, 0.5),
    #z.beta7 = rbinom(1, 1, 0.5), z.beta8 = rbinom(1, 1, 0.5)
    #z.alpha1 = rbinom(1, 1, 0.5),z.alpha2 = rbinom(1, 1, 0.5)
  )
}

# Bundle data
jags.data <- list(counts.multi  = Mature.RSD, 
                  totalC        = Nst, 
                  nsite         = nsite, 
                  nrep          = nrep, 
                  GC_2020       = site.covs2$GC_2020,
                  USJ_2020      = site.covs2$USJ_2020)#,
                  #Nesters       = site.covs2$Nesters, 
                  #HH            = site.covs2$HH,
                  #Depth.max     = site.covs2$Depth.max,
                  #UnderCut.prop = site.covs2$UnderCut.prop, 
                  #Temp          = site.covs2$Temp,
                  #Area          = site.covs2$Area)

# Parameters monitored
jags.params <- c("N.total", "q", "q1","q0","a", "beta0","beta1","beta2","alpha0")#,
                 #"alpha1", "alpha2","z.alpha1", "z.alpha2",
                 #"beta3", "beta4", "beta5", "beta6","beta7", "beta8",
                 #"z.beta1", "z.beta2", "z.beta3", "z.beta4", "z.beta5", 
                 #"z.beta6","z.beta7", "z.beta8")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Run model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
set.seed(1234)
jagsfit2 <- jags(data = jags.data, init = jags.inits, 
                 parameters.to.save = jags.params, model.file = modelFilename,
                 n.chains = 3, n.iter = 500000, n.burnin = 250000, n.thin = 20,
                 parallel = TRUE, n.cores = 4, DIC = TRUE, verbose = TRUE)
jagsfit2

# Rhat - check if convergence occurred
min(unlist(jagsfit2$Rhat))
max(unlist(jagsfit2$Rhat))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# inclusion parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#unlist(jagsfit2$mean$z.beta1) # retained to improve model convergence
#unlist(jagsfit2$mean$z.beta2) # retained as important variable for abundance
#unlist(jagsfit2$mean$z.beta3) # not retained
#unlist(jagsfit2$mean$z.beta4) # not retained
#unlist(jagsfit2$mean$z.beta5) # not retained
#unlist(jagsfit2$mean$z.beta6) # not retained
#unlist(jagsfit2$mean$z.beta7) # not retained
#unlist(jagsfit2$mean$z.beta8) # not retained
#unlist(jagsfit2$mean$z.alpha1) # not retained
#unlist(jagsfit2$mean$z.alpha2) # not retained

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# traceplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
traceplot(jagsfit2, parameters = "beta0")
traceplot(jagsfit2, parameters = "alpha0")
traceplot(jagsfit2, parameters = "beta1")
traceplot(jagsfit2, parameters = "beta2")
traceplot(jagsfit2, parameters = "N.total")

#density plot
densityplot(jagsfit2, parameters = "beta0")
densityplot(jagsfit2, parameters = "alpha0")
densityplot(jagsfit2, parameters = "beta1")
densityplot(jagsfit2, parameters = "beta2")
densityplot(jagsfit2, parameters = "N.total")
densityplot(jagsfit2, parameters = "q1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Covariate effects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# mean and 95% CI 
whiskerplot(jagsfit2,"alpha0",  las = 1)
whiskerplot(jagsfit2,"beta0",   las = 1)
whiskerplot(jagsfit2,"beta1",   las = 1)
whiskerplot(jagsfit2,"beta2",   las = 1)
whiskerplot(jagsfit2,"q",       las = 1)
whiskerplot(jagsfit2,"N.total", las = 1)
whiskerplot(jagsfit2,"a",       las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# detection probability in second and third hauls #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
min(jagsfit2$mean$q[,2])
max(jagsfit2$mean$q[,2])
min(jagsfit2$mean$q[,3])
max(jagsfit2$mean$q[,3])

sum(jagsfit2$mean$q[,2] > max(jagsfit2$mean$q[,1]))
sum(jagsfit2$mean$q[,3] > max(jagsfit2$mean$q[,1]))
38/43

# ~~~~~~~~~~~~~~~~~~~~~ #
# USJ2020 
# ~~~~~~~~~~~~~~~~~~~~~ #
jagsfit2$mean$beta2
jagsfit2$sd$beta2
jagsfit2$q97.5$beta2
jagsfit2$q2.5$beta2
exp(jagsfit2$mean$beta2) # back transform from log

# ~~~~~~~~~~~~~~~~~~~~~ #
# Gully Creek 2020
# ~~~~~~~~~~~~~~~~~~~~~ #
jagsfit2$mean$beta1
jagsfit2$sd$beta1
jagsfit2$q97.5$beta1
jagsfit2$q2.5$beta1
exp(jagsfit2$mean$beta1) # back transform from log
22.98025/exp(jagsfit2$mean$beta1) # how many times greater of an effect?

# ~~~~~~~~~~~~~~~~~~~~~ #
# abundance intercept
# ~~~~~~~~~~~~~~~~~~~~~ #
jagsfit2$mean$beta0
jagsfit2$sd$beta0
jagsfit2$q97.5$beta0
jagsfit2$q2.5$beta0
exp(jagsfit2$mean$beta0) # back transform from log
22.98025/exp(jagsfit2$mean$beta0) # how many times greater of an effect?

# ~~~~~~~~~~~~~~~~~~~~~ #
# detection intercept
# ~~~~~~~~~~~~~~~~~~~~~ #
jagsfit2$mean$alpha0
jagsfit2$sd$alpha0
jagsfit2$q97.5$alpha0
jagsfit2$q2.5$alpha0
exp(jagsfit2$mean$alpha0)/(1+exp(jagsfit2$mean$alpha0)) #inverse logit

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# site abundance plot #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
abund.plot2 <- data.frame(Site      = site.covs2$ID, 
                          Lower95   = unlist(jagsfit2$q2.5$N.total), 
                          site_est  = unlist(jagsfit2$mean$N.total), 
                          Upper95   = unlist(jagsfit2$q97.5$N.total),
                          Waterbody = c(rep('Gully Creek', 16),  
                                        rep('USJ Tributary', 9),
                                        rep('Gully Creek', 18)),
                          Year      = c(rep('2019', 16), 
                                        rep('2020', 9),
                                        rep('2020', 18)),
                          Observed  = Mature.RSD.Detections)
abund.plot2$Event   <- paste0(abund.plot2$Waterbody, " ", abund.plot2$Year)

# plot
site.abundgg<-ggplot(abund.plot2, aes(x=Site, y=site_est, color=Event)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .3) +
  geom_point(aes(x=Site, y=Observed), color='red', size=1, shape=17)+
  scale_color_manual(values=c("black","grey","blue"))+
  guides(color=guide_legend(position = 'inside'))+
  annotate("text",x=1.5, y=250,label="b)")+
  labs(x="Site", y="Abundance") +
  theme(axis.text.x=element_blank(),
        legend.position.inside = c(0.2,0.65),
        legend.title = element_blank())
site.abundgg

min(abund.plot2$site_est[abund.plot2$Waterbody=='Gully Creek' & abund.plot2$Year=='2019'])
max(abund.plot2$site_est[abund.plot2$Waterbody=='Gully Creek' & abund.plot2$Year=='2019'])

min(abund.plot2$site_est[abund.plot2$Waterbody=='Gully Creek' & abund.plot2$Year=='2020'])
max(abund.plot2$site_est[abund.plot2$Waterbody=='Gully Creek' & abund.plot2$Year=='2020'])

min(abund.plot2$site_est[abund.plot2$Waterbody=='USJ Tributary'])
max(abund.plot2$site_est[abund.plot2$Waterbody=='USJ Tributary'])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# baseline detection probability #
# The detection probability that would be observed if there were no 
# change over replicates.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
det.plot2 <- data.frame(Site      = site.covs2$ID, 
                        Lower95   = unlist(jagsfit2$q2.5$q1), 
                        detection = unlist(jagsfit2$mean$q1), 
                        Upper95   = unlist(jagsfit2$q97.5$q1))
ggplot(det.plot2, aes(Site, detection)) +
  geom_point() +
  ylim(0,1)+
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .3) +
  labs(x="Site", y="Baseline detection probability") +
  theme(axis.text.x=element_blank())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Site specific detection probability
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
Pop.df2 <- cbind.data.frame(mean = c(unlist(jagsfit2$mean$q)[,1],
                                     unlist(jagsfit2$mean$q)[,2],
                                     unlist(jagsfit2$mean$q)[,3]),
                            Site = rep(seq(1:43)),
                            Haul = rep(seq(1:3), each = 43),
                            Period = abund.plot2$Event)

#plot
detgg<-ggplot(Pop.df2, aes(x=Haul, y=mean, color=Period))+
  geom_line(aes(group=as.factor(Site)), lwd=0.5)+
  geom_point(aes(group=as.factor(Site)), pch=20)+
  ylim(0,1)+
  scale_color_manual(values=c("black","grey","blue"))+
  scale_x_continuous(breaks=c(1, 2, 3))+
  annotate('text', x=1.1, y=1, label="a)")+
  labs(color="Site", y = "Detection probability") +
  theme(legend.position = "none")
detgg

detgg + site.abundgg +
  plot_layout(ncol = 2, widths = c(0.25,0.75))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# visualize site counts across iterations and sites
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
posterior.list <- cbind.data.frame(posterior=c(jagsfit2$sims.list$N.total))
ggplot(posterior.list, aes(x=posterior))+
  geom_histogram(bins=20, color='black')+
  labs(y="Count", x="Site abundance")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# separate posteriors by site
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
sitepreds2 <- cbind.data.frame(post      = c(jagsfit2$sims.list$N.total),
                               site      = rep(seq(1:43),           each = 37500),
                               Area      = rep(site.covs$Pool_area, each = 37500),
                               Length    = rep(site.covs$P_Length,  each = 37500),
                               Strahler  = rep(site.covs$Strahler,  each = 37500),
                               Waterbody = rep(c(rep('Gully Creek',  16), rep('USJ Tributary', 9),
                                                 rep('Gully Creek',  18)), each = 37500),
                               Year = rep(c(rep('2019', 16), rep('2020', 9),
                                            rep('2020', 18)), each = 37500),
                               Period = rep(c(rep('Gully Creek 2019', 16), 
                                              rep('USJ Tributary 2020', 9),
                                              rep('Gully Creek 2020', 18)), each = 37500))
sitepreds2$ReachLength <- sitepreds2$Strahler

# modify ReachLength to represent the mean ReachLength of strahler segment per waterbody
sitepreds2$ReachLength[sitepreds2$Waterbody == 'Gully Creek' & sitepreds2$ReachLength == 2] <- 15.19
sitepreds2$ReachLength[sitepreds2$Waterbody == 'Gully Creek' & sitepreds2$ReachLength == 3] <- 16.75
sitepreds2$ReachLength[sitepreds2$Waterbody == 'Gully Creek' & sitepreds2$ReachLength == 4] <- 33.54
sitepreds2$ReachLength[sitepreds2$Waterbody == 'USJ Tributary' & sitepreds2$ReachLength == 2] <- 15.19
sitepreds2$ReachLength[sitepreds2$Waterbody == 'USJ Tributary' & sitepreds2$ReachLength == 3] <- 13.15
sitepreds2$ReachLength[sitepreds2$Waterbody == 'USJ Tributary' & sitepreds2$ReachLength == 4] <- 20.74

# make site and strahler characters
sitepreds2$site      <- as.character(sitepreds2$site)
sitepreds2$Strahler  <- as.character(sitepreds2$Strahler)
sitepreds2$ReachLength <- as.numeric(sitepreds2$ReachLength)

# divide the site-predicted abundance by length and area
sitepreds2$fishm2    <- sitepreds2$post/sitepreds2$Area
sitepreds2$fishm     <- sitepreds2$post/sitepreds2$Length
sitepreds2$fishreach <- sitepreds2$post/sitepreds2$ReachLength
sitepreds2$fishkm    <- sitepreds2$fishreach * 1000

# site level fish of pool habitat 
aggregate(sitepreds2$post, list(sitepreds2$Waterbody), median) 
aggregate(sitepreds2$post, list(sitepreds2$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/m of pool habitat 
aggregate(sitepreds2$fishm, list(sitepreds2$Waterbody), median) 
aggregate(sitepreds2$fishm, list(sitepreds2$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/m2 of pool habitat per waterbody
aggregate(sitepreds2$fishm2, list(sitepreds2$Waterbody), median)
aggregate(sitepreds2$fishm2, list(sitepreds2$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/km of pool habitat per waterbody
aggregate(sitepreds2$fishkm, list(sitepreds2$Waterbody), median)
aggregate(sitepreds2$fishkm, list(sitepreds2$Waterbody), quantile, probs=c(0.025, 0.975))

# site level fish/km of pool habitat per waterbody
aggregate(sitepreds2$fishkm, list(sitepreds2$Waterbody, sitepreds2$Strahler), median)
aggregate(sitepreds2$fishkm, list(sitepreds2$Waterbody, sitepreds2$Strahler), quantile, probs=c(0.025, 0.975))

# only sample estimates from gully sites - this excludes USJ tributary sites
sitepreds2.g <- sitepreds2[!c(sitepreds2$site =='17' | sitepreds2$site =='18' | 
                                sitepreds2$site =='19' | sitepreds2$site =='20' | 
                                sitepreds2$site =='21' | sitepreds2$site =='22' | 
                                sitepreds2$site =='23' | sitepreds2$site =='24' | 
                                sitepreds2$site =='25' ),]

# only sample estimates from USJ Tributary
sitepreds2.usj<-sitepreds2[c(sitepreds2$site =='17' | sitepreds2$site =='18' | 
                               sitepreds2$site =='19' | sitepreds2$site =='20' | 
                               sitepreds2$site =='21' | sitepreds2$site =='22' | 
                               sitepreds2$site =='23' | sitepreds2$site =='24' | 
                               sitepreds2$site =='25' ),]

#===============================================================================
#===============================================================================
# Sample the posterior distribution of site abundance to generate estimates of 
# the abundance of RSD across the entire system
# In this iteration separate posterior distributions are used based on the 
# waterbody and strahler order are being considered. 
# Gully Creek 2019 and 2020 posteriors are combined
#===============================================================================
#===============================================================================
reps <- 10000 # number of replicates
set.seed(0528)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Function to sample and sum Gully Creek Strahler 2 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
G.pools2 <- floor(runif(reps, min = 482, max = 1186))
sample_and_sum_St2G <- function(n) {
  sum(sample(sitepreds2.g$post[sitepreds2.g$Strahler == '2'], n, replace = TRUE))
}
Gull.2<- sapply(G.pools2, sample_and_sum_St2G)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Function to sample and sum Gully Creek Strahler 3 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
G.pools3 <- floor(runif(reps, min = 88, max = 216))
sample_and_sum_St3G <- function(n) {
  sum(sample(sitepreds2.g$post[sitepreds2.g$Strahler == '3'], n, replace = TRUE))
}
Gull.3<- sapply(G.pools3, sample_and_sum_St3G)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Function to sample and sum Gully Creek Strahler 4 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
G.pools4 <- floor(runif(reps, min = 143, max = 284))
sample_and_sum_St4G <- function(n) {
  sum(sample(sitepreds2.g$post[sitepreds2.g$Strahler == '4'], n, replace = TRUE))
}
Gull.4<- sapply(G.pools4, sample_and_sum_St4G)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Function to sample and sum USJ Tributary Strahler 3 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
U.pools3 <- floor(runif(reps, min = 113, max = 374))
sample_and_sum_St3U <- function(n) {
  sum(sample(sitepreds2.usj$post[sitepreds2.usj$Strahler == '3'], n, replace = TRUE))
}
USJ.3<- sapply(U.pools3, sample_and_sum_St3U)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Function to sample and sum USJ Tributary Strahler 4 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
U.pools4 <- floor(runif(reps, min = 207, max = 351))
sample_and_sum_St4U <- function(n) {
  sum(sample(sitepreds2.usj$post[sitepreds2.usj$Strahler == '4'], n, replace = TRUE))
}
USJ.4<- sapply(U.pools4, sample_and_sum_St4U)

################################################################################
################################################################################
# Combine data into a df for summarizing
Abund.post2 <- cbind.data.frame(Strahler  = rep(c(2,3,4,3,4), each=reps),
                                Abundance = c(Gull.2, Gull.3, Gull.4, USJ.3, USJ.4),
                                Waterbody = c(rep("Gully Creek",reps*3), rep("USJ Tributary",reps*2)),
                                Pools     = c(G.pools2, G.pools3, G.pools4, U.pools3, U.pools4))
Abund.post2$Strahler  <- as.character(Abund.post2$Strahler)
Abund.post2$fish.pool <- Abund.post2$Abundance/Abund.post2$Pools

# sum abundance across strahler orders
Gull.Pop.Abund<-Gull.2 + Gull.3 + Gull.4
USJ.Pop.Abund<-USJ.3 + USJ.4

# create data frame for plotting
Total.pop <- cbind.data.frame(Pop.Abund = c(Gull.Pop.Abund, USJ.Pop.Abund),
                              Waterbody = rep(c("Gully Creek", "USJ Tributary"), each=10000))

# probability that abundance is above MVP
probability_within_range <- sum(Total.pop$Pop.Abund[Total.pop$Waterbody=="USJ Tributary"] >= 7791) / length(Total.pop$Pop.Abund[Total.pop$Waterbody=="USJ Tributary"])
probability_within_range

probability_within_range <- sum(Total.pop$Pop.Abund[Total.pop$Waterbody=="Gully Creek"] >= 7791) / length(Total.pop$Pop.Abund[Total.pop$Waterbody=="Gully Creek"])
probability_within_range

# fish per pool per waterbody
aggregate(Abund.post2$fish.pool, list(Abund.post2$Waterbody), median)
aggregate(Abund.post2$fish.pool, list(Abund.post2$Waterbody), "quantile", probs=c(0.05, 0.95))

ggplot()+
  geom_histogram(data = Abund.post2, aes(x=fish.pool, fill=Strahler), 
                 alpha=0.5, color="black") +
  labs(x='Abundance', y='Probability Density Function') +
  guides(fill=guide_legend(position='inside'))+
  facet_wrap(~Waterbody, scales='free')+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position.inside = c(0.4,0.75))

ggplot()+
  geom_density(data = Total.pop, aes(x=Pop.Abund, fill=Waterbody), alpha=0.5, color="black") +
  labs(x = "Adult population abundance", y = "Probability density function")+
  scale_fill_manual(values=c("red","grey"))+
  guides(fill=guide_legend(position='inside'))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position.inside = c(0.85,0.85),
        legend.title = element_blank())
