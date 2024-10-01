################################################################################
# This code was written to perform habitat exploration analysis for:
# Habitat associations and abundance estimates for two least-disturbed Redside 
# Dace (Clinostomus elongatus) populations in tributaries of Lake Huron 
# (Ontario, Canada). Canadian Journal of Fisheries and Aquatic Sciences
# Authors: Karl A. Lamothe, D. Andrew R. Drake
# Date revised: 10/1/2024; R Version 4.4.0. 
################################################################################
rm(list=ls())

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

# Boxplots
boxplot.data <- cbind.data.frame(Value=c(site.covs$HH_mean,
                                         site.covs$Depth_max,
                                         site.covs$UnderCut.Prop, 
                                         site.covs$Woody_Debris, 
                                         site.covs$Wtemp,
                                         site.covs$P_Length,
                                         site.covs$Width_mean),
                                 Variable=c(rep("Hydraulic head (mm)",43),  rep("Maximum depth (m)",43),
                                            rep("Under-cut banks (prop)",43), rep("Woody debris (prop)",43),
                                            rep("Temperature (C)",43), 
                                            rep("Pool length (m)", 43),
                                            rep("Pool width (m)", 43)),
                                 Waterbody = site.covs$Waterbody,
                                 RSD.01    = as.factor(site.covs$RSD.01),
                                 Year      = as.factor(site.covs$Year))
boxplot.data$Variable <- factor(boxplot.data$Variable,
                                levels = c("Hydraulic head (mm)", "Maximum depth (m)",
                                           "Under-cut banks (prop)", "Woody debris (prop)",
                                           "Temperature (C)","Pool length (m)",
                                           "Pool width (m)"),
                                labels = c("Hydraulic head (mm)", "Maximum depth (m)",
                                           "Undercut banks (prop)", "Woody debris (prop)",
                                           "Temperature (\u00b0C)","Pool length (m)",
                                           "Pool width (m)"))

boxplot.data$Waterbody[boxplot.data$Waterbody=="Unknown Stan J Trib a"]<-"Unknown Stan J"
boxplot.data$Waterbody[boxplot.data$Waterbody=="Unknown Stan J"]<-"USJ Tributary"

# view
ggplot(boxplot.data, aes(x=Waterbody, y=Value, color=Year))+
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  facet_wrap(.~Variable, scales='free_y', ncol=4)+
  scale_color_manual(values=c("#648FFF","#FE6100"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8))+
  theme(axis.title.y  = element_blank(),
        axis.title.x  = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.title  = element_text(hjust=0.75)) +
  labs(color="Year")

################################################################################
# Consider correlations
colnames(site.covs)
site.covs2<-data.frame(Max.depth     = site.covs$Depth_max,
                       HH            = site.covs$HH_mean,
                       Temp          = site.covs$Wtemp,
                       Nest.builders = site.covs$Co.occurring,
                       Undercut.bank = site.covs$UnderCut.Prop,
                       Woody.debris  = site.covs$Woody_Debris,
                       Pool.length   = site.covs$P_Length,
                       Pool.width    = site.covs$Width_mean)

# Generate correlation plot with pearson correlations
correlationplot<-ggcorrplot(cor(site.covs2,method="spearman"), 
                            hc.order = TRUE, outline.col = "black",
                            type = "lower", lab=TRUE,
                            #ggtheme = ggplot2::theme_gray,
                            insig = "blank")+
  theme(panel.border = element_rect(colour = "black", fill=NA),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
correlationplot

#~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~#
site.covs$wateryear <- paste0(site.covs$Waterbody," ",site.covs$Year)
site.covs$wateryear[site.covs$wateryear=="Unknown Stan J Trib a 2020"]<-"Unknown Stan J 2020"
site.covs$wateryear[site.covs$wateryear=="Unknown Stan J 2020"]<-"USJ Tributary 2020"

# kruskal wallis test
kruskal.test(HH_mean ~ wateryear, data = site.covs)
kruskal.test(Width_mean ~ wateryear, data = site.covs)

# pairwise tests
pairwise.wilcox.test(site.covs$HH_mean, site.covs$wateryear, p.adjust.method = "BH")
pairwise.wilcox.test(site.covs$Width_mean, site.covs$wateryear, p.adjust.method = "BH")

#summary
aggregate(site.covs$Width_mean, list(site.covs$wateryear), mean)
aggregate(site.covs$Width_mean, list(site.covs$wateryear), sd)
aggregate(site.covs$HH_mean, list(site.covs$wateryear), length)
aggregate(site.covs$HH_mean, list(site.covs$wateryear), mean)
aggregate(site.covs$HH_mean, list(site.covs$wateryear), sd)

ord.vars <- cbind.data.frame(Max.dept      = scale(site.covs$Depth_max, center=T),
                             HH            = scale(site.covs$HH_mean, center=T),
                             Temp          = scale(site.covs$Wtemp, center=T),
                             Cooccurring   = scale(site.covs$Co.occurring, center=T),
                             Undercut.prop = scale(site.covs$UnderCut.Prop, center=T),
                             Woody         = scale(site.covs$Woody_Debris, center=T),
                             Pool.length   = scale(site.covs$P_Length, center=T),
                             Pool.width    = scale(site.covs$Width_mean, center=T))

# PERMANOVA
adonis2(ord.vars~site.covs$wateryear, permutations=999, method="euclidean")
adonis2(dist(scale(ord.vars,center=T), method = "euclidean")
        ~site.covs$wateryear, permutations=999, method="euclidean")

# variance
mod3<-betadisper(dist(scale(ord.vars,center=T), method = "euclidean"), group=site.covs$wateryear)
permutest(mod3)

# extract distances to centroid
wateryear_distance <- mod3$distances
wateryear_distance <- cbind.data.frame(Distance=wateryear_distance, wateryear=site.covs$wateryear)
mod3$eig

# plot
mod3_gg<-ggplot(wateryear_distance, aes(x=wateryear, y=Distance))+
  geom_boxplot() +  
  labs(y="Distance") +
  geom_jitter(size=1.5, width=0.1, alpha=0.5) +
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))
mod3_gg

set.seed(876)
permutest(mod3, pairwise = T, permutations = 999)
TukeyHSD(mod3)

# plot
mod3vectors<-data.frame(mod3$vectors)
mod3vectors$wateryear <- site.covs$wateryear

modfit<-envfit(mod3, ord.vars, perm=999)
habvec <- scores(modfit, "vectors")

wateryearplot<-ggplot(mod3vectors, aes(x=PCoA1,y=PCoA2,color=wateryear))+
  geom_hline(yintercept=0, lty="dashed")+
  geom_vline(xintercept=0, lty="dashed")+
  geom_point(size=1)+ 
  scale_color_manual(values=c("black","gold","pink"))+
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
        legend.title = element_blank())
wateryearplot
