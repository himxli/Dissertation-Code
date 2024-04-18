rm(list=ls())

#set working directory here


#install packages
library(ggplot2)
library(vegan)
library(ape)
library(dplyr)

#Code used for RQ1

#Latitude vs richness model 
latitude<- read.csv ("latitude.csv", row.names = 1)

colours <- c( "Richness" = "darkorchid")
shapes <- c("Richness" = 16)
ggplot() +
  
  geom_point(data=latitude, aes(x=Latitude, y=Richness, color="Richness", shape="Richness")) +
 
  scale_color_manual(name = "",
                     labels = c("Richness"),
                     values =c( "darkorchid")) +
  scale_shape_manual(name = "",
                     labels = c("Richness"),
                     values =c(16)) +
 
  geom_smooth(data=latitude, method="lm", se=FALSE, aes(x= Latitude, y=Richness), color = "darkorchid", fullrange= TRUE )+

  
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Latitude") +
  ylab("Shannon's Diversity Index")


#Spearman's rank
corr_richness <- cor.test(x=latitude$Latitude, y=latitude$Richness, method = 'spearman')
corr_richness

#Latitude vs evenness model
colours <- c( "Evenness" = "gold")
shapes <- c("Evenness" = 16)
ggplot() +
  
  geom_point(data=latitude, aes(x=Latitude, y=Evenness, color="Richness", shape="Richness")) +
  
  scale_color_manual(name = "",
                     labels = c("Evenness"),
                     values =c( "gold")) +
  scale_shape_manual(name = "",
                     labels = c("Evenness"),
                     values =c(16)) +
  
  geom_smooth(data=latitude, method="lm", se=FALSE, aes(x= Latitude, y=Evenness), color = "gold", fullrange= TRUE )+
  
  
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Latitude") +
  ylab("Pielou's Evenness Index")

#Spearman's rank
corr_evenness <- cor.test(x=latitude$Latitude, y=latitude$Evenness, method = 'spearman')
corr_evenness


#Code used for RQ2

#boxplot for diversity
highlowland<- read.csv("highlandlowland.csv", row.names = 1)

(diversity_boxplot <- ggplot(highlowland, aes(location, Richness)) + 
    geom_boxplot(aes(fill = location)) +
    theme_bw() +
    scale_fill_manual(values = c("darkseagreen", "cornflowerblue")) +               
    scale_colour_manual(values = c("darkseagreen", "cornflowerblue")) +            
    ylab("Shannon's Diversity Index\n") +                             
    xlab("\n Location")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(),                                                
          legend.position = "none"))       
#t-test
shapiro.test(highlowland$Richness)
t.test(Richness ~ location, data = highlowland)


#boxplot for evenness

(diversity_boxplot <- ggplot(highlowland, aes(location, Evenness)) + 
    geom_boxplot(aes(fill = location)) +
    theme_bw() +
    scale_fill_manual(values = c("darkseagreen", "cornflowerblue")) +               
    scale_colour_manual(values = c("darkseagreen", "cornflowerblue")) +             
    ylab("Pielou's Evenness Index\n") +                             
    xlab("\n Location")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(),                                
          legend.position = "none"))      

#t test for evenness
shapiro.test(highlowland$Evenness)
t.test(Evenness ~ location, data = highlowland)


#NMDS plot
nmds_rq2 <- read.csv ("NMDS_attempt.csv", row.names = 1)

dist <- vegdist(nmds_rq2,  method = "bray")

NMDS.scree <- function(nmds_rq2) 
  plot(rep(1, 10), replicate(10, metaMDS(nmds_rq2, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(nmds_rq2, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(dist)

set.seed(2)

NMDS3 <- metaMDS(nmds_rq2, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
plot(NMDS3)

colors = c(rep("cornflowerblue", 5), rep("darkseagreen", 7))

ordiplot(NMDS3, type = "n")
for(i in unique(group)) {
  ordihull(NMDS3$point[grep(i, group),], draw="polygon",
           groups = group[group == i],col = colors[grep(i,group)],label=F) } 


orditorp(NMDS3, display = "sites", col = c(rep("black",5),
                                         rep("black", 7)), air = 0.01, cex = 1.25)
#ANOSIM 
anosim_attempt = read.csv("NMDS_attempt.csv", header= TRUE)
com = anosim_attempt[,3:ncol(anosim_attempt)]
m_com = as.matrix(com)
ano = anosim(m_com, anosim_attempt$Location, distance = "bray", permutations = 9999)
ano


#Code used for RQ3

eastwest<- read.csv("east_west.csv", row.names = 1) #SPREADSHEET FOR BOXPLOT
eastwest2<- read.csv("east_west2.csv", header = TRUE) #SPREADSHEET FOR ANOSIM

#ANOSIM
com = eastwest2[,3:ncol(eastwest2)]
n_com = as.matrix(com)
ano = anosim(n_com, eastwest2$Location, distance = "bray", permutations = 9999)
ano

#richness boxplot
(diversity_boxplot <- ggplot(eastwest, aes(location, Richness)) + 
    geom_boxplot(aes(fill = location)) +
    theme_bw() +
    scale_fill_manual(values = c("darksalmon", "darkgoldenrod2")) +               
    scale_colour_manual(values = c("darksalmon", "darkgoldenrod2")) +             
    ylab("Shannon's Diversity Index\n") +                             
    xlab("\n Location")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(),                                                
          legend.position = "none"))    

#t-test for richness
t.test(Richness ~ location, data = eastwest)
shapiro.test(eastwest$Richness)

#evenness boxplot
(evenness2_boxplot <- ggplot(eastwest, aes(location, Evenness)) + 
    geom_boxplot(aes(fill = location)) +
    theme_bw() +
    scale_fill_manual(values = c("darksalmon", "darkgoldenrod2")) +               
    scale_colour_manual(values = c("darksalmon", "darkgoldenrod2")) +             
    ylab("Pielou's Evenness Index\n") +                             
    xlab("\n Location")  +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),                     
          panel.grid = element_blank(),                                                
          legend.position = "none")) 

#t-test for evenness
t.test(Evenness ~ location, data = eastwest)
shapiro.test(eastwest$Evenness)
