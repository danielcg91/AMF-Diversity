#Load all necessary libraries
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(phyloseq)
library(tidyr)
library(gridExtra)
library(dplyr)
library(FSA)


#Choose all setting to plots
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5, family = "sans" ),
  #axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
  axis.text.y = element_text(color="black",size=12, vjust = 1.5, family = "sans"),
  axis.title = element_text(color="black",size=12, face = "bold", family = "sans"),
  axis.title.x.bottom = element_blank(),
  panel.border =element_rect(color = "black", fill = NA),#element_blank(),
  strip.text.x = element_text(size=12, color="black",face="bold", family = "sans"),
  strip.text.y = element_text(size=12, color="black",face="bold", family = "sans"),
  strip.placement = "outside", strip.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
  panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
  legend.position = "right", legend.text = element_text(color = "black",size=12, family = "sans"), 
  legend.direction = "vertical", legend.title = element_text(color = "black",size=12, face = "bold", family = "sans"),
  legend.key.size = unit(0.4,"cm"))

# Set PARAMETERS
data=c("daniel", "victor", "all")[1]
amplicon=c("ITS2","16S")[2]
base=c("uniteall", "silva")[2]
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F","#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","dimgrey"))
pc=c("soil", "rhizosphere", "root endosphere", "spores", "T0", "T4", "T8", "T12", 
     "Agave.tequilana", "Agave.salmiana", "Myrtillocactus.geometrizans")


#Upload data R Base
if (amplicon == "ITS2"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agaveITS/ITS2_1/") 
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$AIM.sample.id
  met=met[order(met$AIM.sample.id),]
  met=met[order(met$AIM.sample.id),]
  met=met[colnames(otu),]
  #Subset of data to process only at 12 months and soil
  met= subset(met, met$treatment != "T4")
  met= subset(met, met$treatment != "T8")
  otu=otu[rownames(met)]
} else if ( amplicon == "16S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$AIM.sample.id
  met=met[order(met$AIM.sample.id),]
  met=met[order(met$AIM.sample.id),]
  met=met[colnames(otu),]
  #Subset of data to process only at 12 months and soil
  met= subset(met, met$treatment != "T4")
  met= subset(met, met$treatment != "T8")
  otu=otu[rownames(met)]
}


#Rarefy otu table
set.seed(23171341)
#all fungi
ric=(rarefy(t(otu), sample = min(colSums(otu))))
#all fungi
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
#Checklist
colSums(otu.s)
shan=diversity(t(otu.s),index="shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan)
alpha=alpha[met$AIM.sample.id,]
alpha=cbind(alpha, met[,c( "plant.species","plant.compartment","treatment","AIM.sample.id")])
alpha$treatment=factor(alpha$treatment, level = pc)
alpha$plant.species=factor(alpha$plant.species, level = pc)
alpha$plant.compartment=factor(alpha$plant.compartment, level = pc)

#KW by plant.compartment
kruskal.test(richness~plant.compartment, data = alpha) ; dunnTest(richness~plant.compartment, data = alpha, method = "bh")
kruskal.test(shannon~plant.compartment, data = alpha) ; dunnTest(shannon~plant.compartment, data = alpha, method = "bh")

#KW by treatment
kruskal.test(richness~treatment, data = alpha) ; dunnTest(richness~treatment, data = alpha, method = "bh")
kruskal.test(shannon~treatment, data = alpha) ; dunnTest(shannon~treatment, data = alpha, method = "bh")

#KW by plant.species
kruskal.test(richness~plant.species, data = alpha) ; dunnTest(richness~plant.species, data = alpha, method = "bh")
kruskal.test(shannon~plant.species, data = alpha) ; dunnTest(shannon~plant.species, data = alpha, method = "bh")

#ANOVA by plant.compartment
test=aov(richness~plant.compartment, data = alpha) ; print(test) ; TukeyHSD(test)
test=aov(shannon~plant.compartment, data = alpha) ; print(test) ; TukeyHSD(test)

#Select the folder to save the figures
if (amplicon == "ITS2"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURAS-ITS/")
} else if ( amplicon == "16S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURAS-16S/")
} 

mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "plant.compartment")
mean.r$plant.compartment=factor(mean.r$plant.compartment, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "plant.compartment")
mean.s$plant.compartment=factor(mean.s$plant.compartment, level = pc)

#RICHENESS BY PLANT.COMPARTMENT
pdf(paste(amplicon, file = "plantcompartment.richenessplot.pdf"), colormodel = "cmyk", width = 8.5, height =8.5, compress = F)
p=ggplot(data=mean.r, aes(x=plant.compartment, y=richness))+
  geom_bar(stat="identity", fill=c("gray100","gray80","gray50", "gray30"), colour="black",size=1)+
  geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:8])+tema
p
dev.off()
#SHANNON BY PLANT.COMPARTMENT
pdf(paste(amplicon, file = "plantcompartment.shannonplot.pdf"), colormodel = "cmyk", width = 8.5, height = 8.5, compress = F)
p1=ggplot(data=mean.s, aes(x=plant.compartment, y=shannon))+
  geom_bar(stat="identity",fill=c("gray100","gray80","gray50", "gray30"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()


####NMDS#####

#Beta diversity
otu.n=as.data.frame(t(t(otu.s)/colSums(otu.s)*100)) #relative abundance
otu.n=log10(otu.n+1)
otu.n=sqrt(otu.n)
#calculate dimensions
scaling=vegdist(t(otu.n), method = "bray", binary = T, na.rm = FALSE) #calculate distance
#scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
#scaling2=metaMDS(otu.n, distance = "bray") ; scaling2$stress
scaling2=monoMDS(scaling); scaling2$stress
scaling3=data.frame(scaling2$points) #select cordenates
scaling3=cbind(scaling3,alpha)

#plot
#plant.compartment
pdf(paste(amplicon, file = "NMDS.plantcompartment.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
n=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=plant.compartment, shape=plant.species))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#FF876F","#AC4EA7","#F05D8E","#3F4DAA"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:18))+
  labs(x = "NMDS1",y = "NMDS2")+tema
n
dev.off()

#plant.compartment
pdf(paste(amplicon,file = "NMDS.plantcompartment-1.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
m=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=plant.compartment, shape=plant.species))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#FF876F","#AC4EA7","#F05D8E","#3F4DAA"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:18))+
  labs(x = "NMDS1",y = "NMDS2")+tema+
  theme(legend.position = "none")
print(m)
dev.off()

#Plant species
pdf(paste(amplicon, file = "ITSNMDS.plantspecies.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
a=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=plant.species, shape=plant.compartment))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#FF876F","#AC4EA7","#F05D8E","#3F4DAA"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:18))+
  labs(x = "NMDS1",y = "NMDS2")+tema
print(a)
dev.off()


mod = with(met, betadisper(scaling, plant.compartment))
alpha$distance_to_centroid=mod$distances
mean.c=summarySE(data = alpha, measurevar = "distance_to_centroid", groupvars = "plant.compartment")
mean.c$plant.compartment=factor(mean.c$plant.compartment, level = pc)
######
pdf(paste(amplicon, file = "plantcompartment.distances.pdf"),colormodel = "srgb", width = 8.5, height = 8.5, compress = F)
d=ggplot(data=mean.c, aes(x=plant.compartment, y=distance_to_centroid, fill=plant.compartment))+
  geom_bar(stat="identity",colour="black",size=1, fill=c("gray100","gray80","gray50", "gray30"))+
  geom_errorbar(size=1,aes(ymin=distance_to_centroid-sd, ymax=distance_to_centroid+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values = my_pal(12)[c(1,4,7,9)])+
  scale_y_continuous(limits = c(0,0.6), breaks = seq(0,1,0.2))+tema
d
dev.off()
TukeyHSD(mod)
#KW by distance_to_centroid
kruskal.test(distance_to_centroid~plant.compartment, data = alpha) ; dunnTest(distance_to_centroid~plant.compartment, data = alpha, method = "bh")
kruskal.test(distance_to_centroid~plant.species, data = alpha) ; dunnTest(distance_to_centroid~plant.species, data = alpha, method = "bh")
#kruskal.test(distance_to_centroid~treatment, data = alpha) ; dunnTest(distance_to_centroid~treatment, data = alpha, method = "bh")

plot(mod)

###PERMANOVA ANALYSIS
#Beta diversity
otu.n=as.data.frame(t(t(otu.s)/colSums(otu.s)*100)) #relative abundance
otu.n=log10(otu.n+1)
otu.n=sqrt(otu.n)
#
p=t(otu.n)
set.seed(173612)
scaling=vegdist(p, method = "bray", binary = T) 
adonis2(scaling~plant.compartment, data = met, permutations = 1000)
adonis2(scaling~plant.species, data = met, permutations = 1000)
adonis2(scaling~treatment, data = met, permutations = 1000)
adonis2(scaling~plant.compartment*plant.species*treatment, data = met, permutations = 1000)
