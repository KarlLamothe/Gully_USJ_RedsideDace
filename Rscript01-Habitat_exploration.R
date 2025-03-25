################################################################################
# This code was written to perform habitat exploration analysis for:
# Habitat associations and abundance estimates for two least-disturbed Redside 
# Dace (Clinostomus elongatus) populations in tributaries of Lake Huron 
# (Ontario, Canada). Canadian Journal of Fisheries and Aquatic Sciences
# Authors: Karl A. Lamothe, D. Andrew R. Drake
# Date revised: 03/25/2025; R Version 4.4.0. 
################################################################################
# load libraries
library(pacman)
p_load(ggplot2)
p_load(ggcorrplot)   
p_load(vegan)

# Personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=9, family="sans", colour="black"),
                  axis.text.x  = element_text(size=8, family="sans", colour="black"),
                  axis.text.y  = element_text(size=8, family="sans", colour="black"),
                  legend.title = element_text(size=9, family="sans", colour="black"),
                  legend.text  = element_text(size=8, family="sans", colour="black"),
                  plot.title   = element_text(size=10, family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black")))

################################################################################
# read in habitat data and RSD mature data
################################################################################
site.covs  <- read.csv("Site_covariates.csv", header=T)
RSD.counts <- read.csv("Adult_Counts.csv",    header=T)

#  Separate adult counts
Mature.RSD.Detections <- RSD.counts[,c(4:6)]
Mature.RSD.Detections <- rowSums(Mature.RSD.Detections)

# average abundance across sites
sum(Mature.RSD.Detections)/43

# add in site abundance to site covariates
site.covs$RSD <- Mature.RSD.Detections
site.covs$RSD.01 <- site.covs$RSD # presence absence
site.covs$RSD.01[site.covs$RSD.01>0] <-1

# factors
site.covs$Strahler <- as.character(site.covs$Strahler)

################################################################################
# Site covariates
################################################################################
colnames(site.covs)

# Water quality variables
Water.var <- cbind.data.frame(Value    = c(site.covs$Wtemp, site.covs$Turb.ntu,
                                           site.covs$Cond, site.covs$DO),
                              Variable = rep(c("Water temp.", "Turbidity",
                                               "Conductivity", "Dissolved oxygen"), each=43),
                              RSD      = rep(site.covs$RSD, 4),
                              RSD.01   = rep(site.covs$RSD.01, 4),
                              Waterbody= rep(site.covs$Waterbody, 4))

ggplot(Water.var, aes(x=Value, y=RSD, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside Dace count") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

ggplot(Water.var, aes(x=Value, y=RSD.01, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside detections") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

# Site characteristics
Site.var <- cbind.data.frame(Value    = c(site.covs$Width_mean, site.covs$P_Length,
                                          site.covs$Pool_area, site.covs$Pool_vol,
                                          site.covs$UnderCut.Prop,site.covs$Woody_Debris),
                             Variable = rep(c("Mean pool width", "Pool length",
                                              "Pool area", "Pool volume",
                                              "Prop. undercut banks","Prop. woody debris"), each=43),
                             RSD      = rep(site.covs$RSD, 6),
                             RSD.01   = rep(site.covs$RSD.01, 6),
                             Waterbody= rep(site.covs$Waterbody, 6))

ggplot(Site.var, aes(x=Value, y=RSD, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside Dace count") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

ggplot(Site.var, aes(x=Value, y=RSD.01, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside detections") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

# Depth velocity
Depth.HH.var <- cbind.data.frame(Value = c(site.covs$Depth_mean, site.covs$Depth_max,
                                           site.covs$Depth_sd, site.covs$HH_mean,
                                           site.covs$HH_max, site.covs$HH_sd),
                                 Variable = rep(c("Mean depth", "Max depth",
                                                  "Depth sd", "Mean HH",
                                                  "Max HH", "HH sd"), each=43),
                                 RSD      = rep(site.covs$RSD, 6),
                                 RSD.01   = rep(site.covs$RSD.01, 6),
                                 Waterbody= rep(site.covs$Waterbody, 6))

ggplot(Depth.HH.var, aes(x=Value, y=RSD, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside Dace count") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

ggplot(Depth.HH.var, aes(x=Value, y=RSD.01, color=Waterbody))+
  geom_point()+
  facet_wrap(~Variable, scales='free') +
  labs(y = "Adult Redside detections") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

################################################################################
# Boxplots
################################################################################
boxplot.data <- cbind.data.frame(Value=c(site.covs$HH_mean,
                                         site.covs$Depth_max,
                                         site.covs$UnderCut.Prop, 
                                         site.covs$Woody_Debris, 
                                         site.covs$Wtemp,
                                         site.covs$Pool_area),
                                 Variable=c(rep("Hydraulic head (mm)",43),  rep("Maximum depth (m)",43),
                                            rep("Under-cut banks (prop.)",43), rep("Woody debris (prop.)",43),
                                            rep("Temperature (C)",43), rep("Pool area (m2)", 43)),
                                 Waterbody = site.covs$Waterbody,
                                 RSD.01    = as.factor(site.covs$RSD.01),
                                 Year      = as.factor(site.covs$Year))

str(boxplot.data)

ggplot(boxplot.data, aes(x=Waterbody, y=Value, color=Year))+
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  facet_wrap(.~Variable, scales='free_y', ncol=3)+
  scale_color_manual(values=c("#648FFF","#FE6100"))+
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(hjust=0.75)) +
  labs(color="Year")

ggplot(boxplot.data, aes(x=Waterbody, y=Value, color=RSD.01))+
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  facet_wrap(.~Variable, scales='free_y', ncol=3)+
  scale_color_manual(values=c("#648FFF","#FE6100"))+
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(hjust=0.75)) +
  labs(color="RSD")

################################################################################
# Consider correlations
################################################################################
colnames(site.covs)
site.covs2<-data.frame(Max.depth     = site.covs$Depth_max,
                       Mean.HH       = site.covs$HH_mean,
                       Water.temp    = site.covs$Wtemp,
                       Nest.builders = site.covs$Co.occurring,
                       Undercut.bank.prop. = site.covs$UnderCut.Prop,
                       Woody.debris.prop   = site.covs$Woody_Debris,
                       Pool.area     = site.covs$Pool_area,
                       Waterbody     = site.covs$Waterbody,
                       Event         = paste(site.covs$Waterbody, site.covs$Year)
)
str(site.covs2)

site.covs2$Nest.builders <- as.character(site.covs2$Nest.builders)
site.covs2$Nest.builders <- as.numeric(site.covs2$Nest.builders)

# Generate correlation plot with pearson correlations
ggcorrplot(cor(site.covs2[1:7],method="spearman"), 
                            hc.order = TRUE, outline.col = "black",
                            type = "lower", lab=TRUE,
                            insig = "blank")+
  theme(panel.border = element_rect(colour = "black", fill=NA),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

################################################################################
# Perform permanova
################################################################################
ord.vars <- cbind.data.frame(Max.depth     = scale(site.covs$Depth_max, center=T),
                             HH            = scale(site.covs$HH_mean, center=T),
                             Temp          = scale(site.covs$Wtemp, center=T),
                             Cooccurring   = scale(site.covs$Co.occurring, center=T),
                             Undercut.prop = scale(site.covs$UnderCut.Prop, center=T),
                             Woody         = scale(site.covs$Woody_Debris, center=T),
                             Pool.area     = scale(site.covs$P_Length, center=T))

# based on sampling period
# analysis of variance
event.perm<-adonis2(ord.vars~site.covs2$Event, permutations=999, method="euclidean")
event.perm

#adonis2(dist(scale(ord.vars,center=T), method = "euclidean")
#        ~site.covs$wateryear, permutations=999, method="euclidean")
pairwise.adonis(ord.vars,site.covs2$Event, sim.method = "euclidean")

# homogeneity of variance
mod3 <- betadisper(dist(scale(ord.vars,center=T), method = "euclidean"), group=site.covs2$Event)

# extract distances to centroid
Event_distance <- mod3$distances
Event_distance <- cbind.data.frame(Distance=Event_distance, Event=site.covs2$Event)

# plot
ggplot(Event_distance, aes(x=Event, y=Distance))+
  geom_boxplot() +  
  labs(y="Distance") +
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))

# extract vectors for plotting
mod3vectors<-data.frame(mod3$vectors)
mod3vectors$Event <- site.covs2$Event

# fit environmental vectors onto an Ordination
modfit<-envfit(mod3, ord.vars, perm=999)
habvec <- scores(modfit, "vectors")

# plot
Eventplot<-ggplot(mod3vectors, aes(x=PCoA1,y=PCoA2,color=Event))+
  geom_hline(yintercept=0, lty="dashed")+
  geom_vline(xintercept=0, lty="dashed")+
  geom_point(size=1)+ 
  scale_color_manual(values=c("black","grey","blue"))+
  stat_ellipse(type='norm', lwd=0.4)+
  geom_segment(data=habvec, aes(x=0,y=0,xend=PCoA1*2, yend=PCoA2*2), 
               arrow=arrow(length=unit(0.2,"cm")),
               lwd =0.4, colour = "black") + 
  geom_text(data=habvec, aes(x=PCoA1*2, y=PCoA2*2),
            label=rownames(habvec),
            inherit.aes = F, 
            nudge_y = ifelse(habvec[,2]*2 > 0, 0.1, -0.1),
            nudge_x = ifelse(habvec[,1]*2 > 0, 0.1, -0.1),
            size=3)+
  xlab("PCA Axis 1") + ylab("PCA Axis 2")+
  coord_fixed()+
  labs(color="Year and Waterbody")+
  theme(legend.background = element_blank(),
        legend.key=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Eventplot

################################################################################
# based on redside date detections
################################################################################
# not significant
adonis2(ord.vars~as.factor(site.covs$RSD.01), permutations=999, method="euclidean")

# not significant
mod4<-betadisper(dist(scale(ord.vars,center=T), method = "euclidean"), group=site.covs$RSD.01)
permutest(mod4, pairwise = T, permutations = 999)

# extract distances to centroid
RSD.01_distance <- mod4$distances
RSD.01_distance <- cbind.data.frame(Distance=RSD.01_distance, RSD.01=site.covs$RSD.01)

# plot
ggplot(RSD.01_distance, aes(x=factor(RSD.01), y=Distance))+
  geom_boxplot() +  
  labs(y="Distance") +
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))

mod4vectors<-data.frame(mod4$vectors)
mod4vectors$RSD.01 <- site.covs$RSD.01
mod4vectors$RSD.01 <- ifelse(mod4vectors$RSD.01 == 1, "RSD Detected", "RSD Not detected")

multiellipsegg<-ggplot(mod3vectors, aes(x=PCoA1,y=PCoA2,color=Event))+
  geom_hline(yintercept=0, lty="dashed", lwd=0.25)+
  geom_vline(xintercept=0, lty="dashed", lwd=0.25)+
  geom_point(size=1)+ 
  scale_color_manual(values=c("red","purple","black","grey","blue"))+
  stat_ellipse(type='norm', lwd=0.4)+
  geom_segment(data=habvec, aes(x=0,y=0,xend=PCoA1*2, yend=PCoA2*2), 
               arrow=arrow(length=unit(0.2,"cm")),
               lwd =0.4, colour = "black") + 
  geom_text(data=habvec, aes(x=PCoA1*2, y=PCoA2*2),
            label=rownames(habvec),
            inherit.aes = F, 
            nudge_y = ifelse(habvec[,2]*2 > 0, 0.1, -0.1),
            nudge_x = ifelse(habvec[,1]*2 > 0, 0.1, -0.1),
            size=3)+
  stat_ellipse(data = mod4vectors, aes(x=PCoA1,y=PCoA2, color=factor(RSD.01)), 
               linetype="dashed")+
  geom_point(data = mod4vectors, aes(x=PCoA1,y=PCoA2,color=factor(RSD.01)))+ 
  xlab("PCA Axis 1") + ylab("PCA Axis 2")+
  coord_fixed()+
  labs(color="Year and Waterbody")+
  theme(legend.background = element_blank(),
        legend.key=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
multiellipsegg

#~~~~~~~~~~~~~~~~~~~~~~#
# test between sampling periods
aggregate(site.covs2$Mean.HH, list(site.covs2$Event), mean)
aggregate(site.covs2$Mean.HH, list(site.covs2$Event), sd)

kruskal.test(Mean.HH ~ Event, data = site.covs2)
pairwise.wilcox.test(site.covs2$Mean.HH, site.covs2$Event, p.adjust.method = "BH")

aggregate(site.covs2$Pool.area, list(site.covs2$Event), mean)
aggregate(site.covs2$Pool.area, list(site.covs2$Event), sd)

kruskal.test(Pool.area ~ Event, data = site.covs2)
pairwise.wilcox.test(site.covs2$Pool.area, site.covs2$Event, p.adjust.method = "BH")
