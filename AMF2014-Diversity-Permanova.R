setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")
#Load neccesary libraries
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

#All features to graph the data
tema=theme(#axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5),
  axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
  axis.text.y = element_text(color="black",size=12, vjust = 1.5),
  axis.title = element_text(color="black",size=12, face = "bold"),
  axis.title.x.bottom = element_blank(),
  panel.border =element_rect(color = "black", fill = NA),#element_blank(),
  strip.text.x = element_text(size=12, color="black",face="bold"),
  strip.text.y = element_text(size=12, color="black",face="bold"),
  strip.placement = "outside", strip.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
  panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
  legend.position = "right", legend.text = element_text(color = "black",size=12), legend.direction = "vertical",
  legend.title = element_text(color = "black",size=14, face = "bold"),
  legend.key.size = unit(0.3,"cm"))

#Factores
ranks=c("domain","phylum","class","order","family","genus","otu.id")
data.dir=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
colores=c("#2C9854","#85C33B","#C3B63B","#C33B3B","#C33BB4","#3B44C3", "#3B9AC3", "grey80","grey20")
pc=c("soil", "leaf.endosphere", "phyllosphere", "rhizosphere", "roo.zone.soil" , "root.endosphere", "dry", "rainy", 
     "Amatitan","Penjamo","Agave.Hill","Boyd.Ridge","Pinyon.Flat","Magueyal","SanFelipe",
     "Agave.tequilana","Agave.deserti","Agave.salmiana","M.geometrizans","O.robusta")
pc1=c("SanFelipe", "Magueyal", "Boyd.Ridge", "Pinyon.Flat", "Agave.Hill")
pc2=c("Agave.deserti","Agave.salmiana","M.geometrizans","O.robusta")
#Load data
#potu=read.delim(paste(data.dir,"base16S_agaves-cactus_measureable.txt",sep = "/"), header = T)
## Upload data R base
otu=read.delim(paste(data.dir,"baseITS_agaves-cactus_measureable.txt",sep = "/"), header = T)
otu$OTU.ID=paste("OTU",otu$OTU.ID,sep="_")
rownames(otu)= otu$OTU.ID; otu=otu[,-1]
otu = otu[,!grepl("D.le",names(otu))]
otu = otu[,!grepl("R.le",names(otu))]
otu = otu[,!grepl("D.e.",names(otu))]
otu = otu[,!grepl("R.e.",names(otu))]
otu = otu[,!grepl("At.",names(otu))]
meta=read.delim(paste(data.dir,"metadataAMF_AgaveCactus_measureable.txt",sep = "/"), header = T, stringsAsFactors = F)
meta= subset(meta, meta$Specie != "Agave.tequilana")
rownames(meta)= meta$SampleID
tax=read.delim(paste(data.dir,"TaxaITS_agaves-cactus_measureable.txt",sep = "/"),header = T, stringsAsFactors = F)
tax$OTU.ID=paste("OTU",tax$OTU.ID,sep="_")
rownames(tax)= tax$OTU.ID; tax=tax[,-1]
#write.table(otu, file = paste(paste("taxasITS.table","txt",sep = "."),sep = "/"),
            #row.names = T, quote = F, sep = "\t", col.names = T)

#Rarefy otu table
set.seed(23171341)
#AMF
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
#Checklist
colSums(otu.s)
#Subset AMF
tax= subset(tax, tax$phylum == "p__Glomeromycota", select= c("kingdom","phylum","class","order","family","genus"))
otu=otu[rownames(tax),]
#subset only amf otus
otu.s=as.data.frame(otu.s[rownames(tax),])
#Calculate alpha-diversity
ric=(rarefy(t(otu), sample = min(colSums(otu))))
#ric=(rarefy(t(otu), sample = 7161))
shan=diversity(t(otu.s),index = "shannon")

#Diversity analysis
alpha=data.frame(richness=ric, shannon=shan)
alpha=cbind(alpha, meta[,c("Sample", "Season","Specie", "Location")])
alpha$Specie=factor(alpha$Specie, level = pc2)

####ALPHA DIVERSITY
#KW by plant.compartment
kruskal.test(richness~Sample, data = alpha) ; dunnTest(richness~Sample, data = alpha, method = "bh")
kruskal.test(shannon~Sample, data = alpha) ; dunnTest(shannon~Sample, data = alpha, method = "bh")
#KW by plant.species
kruskal.test(richness~Specie, data = alpha) ; dunnTest(richness~Specie, data = alpha, method = "bh")
kruskal.test(shannon~Specie, data = alpha) ; dunnTest(shannon~Specie, data = alpha, method = "bh")
#KW by Location
kruskal.test(richness~Location, data = alpha) ; dunnTest(richness~Location, data = alpha, method = "bh")
kruskal.test(shannon~Location, data = alpha) ; dunnTest(shannon~Location, data = alpha, method = "bh")
#Ttest by Season
t.test(richness~Season, data = alpha)
t.test(shannon~Season, data = alpha)


####Barplot of RICHNESS AND SHANNON###
#SEASON
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "Season")
mean.r$Season=factor(mean.r$Season, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "Season")
mean.s$Season=factor(mean.s$Season, level = pc)
#RICHENESS BY SEASON
pdf(file = "AMF.season.richnessplot.pdf")
p=ggplot(data=mean.r, aes(x=Season, y=richness))+
  geom_bar(stat="identity", fill=c("lemonchiffon4","tan3"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:9])+tema
p
dev.off()
#SHANNON BY SEASON
pdf(file = "AMF.season.shannonplot.pdf")
p1=ggplot(data=mean.s, aes(x=Season, y=shannon))+
  geom_bar(stat="identity",fill=c("lemonchiffon4","tan3"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()
#SPECIE
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "Specie")
mean.r$Specie=factor(mean.r$Specie, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "Specie")
mean.s$Specie=factor(mean.s$Specie, level = pc)
#RICHENESS BY SPECIE
pdf(file = "AMF.specie.richnessplot.pdf")
p=ggplot(data=mean.r, aes(x=Specie, y=richness))+
  geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna", "lightblue"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:9])+tema
p
dev.off()
#SHANNON BY SPECIE
pdf(file = "AMF.specie.shannonplot.pdf")
p1=ggplot(data=mean.s, aes(x=Specie, y=shannon))+
  geom_bar(stat="identity",fill=c("lemonchiffon4","tan3","navajowhite4", "sienna","lightblue"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()
##SAMPLE
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "Sample")
mean.r$Sample=factor(mean.r$Sample, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "Sample")
mean.s$Sample=factor(mean.s$Sample, level = pc)
#RICHENESS BY Sample
pdf(file = "AMF.sample.richnessplot.pdf")
p=ggplot(data=mean.r, aes(x=Sample, y=richness))+
  geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:9])+tema
p
dev.off()
#SHANNON BY SAMPLE
pdf(file = "AMF.sample.shannonplot.pdf")
p1=ggplot(data=mean.s, aes(x=Sample, y=shannon))+
  geom_bar(stat="identity",fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()

##LOCATION
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "Location")
mean.r$Location=factor(mean.r$Location, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "Location")
mean.s$Location=factor(mean.s$Location, level = pc)
#RICHENESS BY Location
pdf(file = "AMF.location.richnessplot.pdf")
p=ggplot(data=mean.r, aes(x=Location, y=richness))+
  geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna", "lightblue", "red", "grey20"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:9])+tema
p
dev.off()
#SHANNON BY LOCATION
pdf(file = "AMF.location.shannonplot.pdf")
p1=ggplot(data=mean.s, aes(x=Location, y=shannon))+
  geom_bar(stat="identity",fill=c("lemonchiffon4","tan3","navajowhite4", "sienna","lightblue" ,"red", "grey20"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()

#############
###NMDS####
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")
ranks=c("domain","phylum","class","order","family","genus","otu.id")
data.dir=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")## Upload data R base
## Upload data R base
otu=read.delim(paste(data.dir,"baseITS_agaves-cactus_measureable.txt",sep = "/"), header = T)
otu$OTU.ID=paste("OTU",otu$OTU.ID,sep="_")
rownames(otu)= otu$OTU.ID; otu=otu[,-1]
otu = otu[,!grepl("D.le",names(otu))]
otu = otu[,!grepl("R.le",names(otu))]
otu = otu[,!grepl("D.e.",names(otu))]
otu = otu[,!grepl("R.e.",names(otu))]
otu = otu[,!grepl("At.",names(otu))]
meta=read.delim(paste(data.dir,"metadataAMF_AgaveCactus_measureable.txt",sep = "/"), header = T, stringsAsFactors = F)
meta= subset(meta, meta$Specie != "Agave.tequilana")
rownames(meta)= meta$SampleID
tax=read.delim(paste(data.dir,"TaxaITS_agaves-cactus_measureable.txt",sep = "/"),header = T, stringsAsFactors = F)
tax$OTU.ID=paste("OTU",tax$OTU.ID,sep="_")
rownames(tax)= tax$OTU.ID; tax=tax[,-1]

#rarefy OTU
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
#otu.s=t(rrarefy(t(otu), sample = 20000))
#Checklist
colSums(otu.s)
#Beta diversity
otu.n=as.data.frame(t(t(otu.s)/colSums(otu.s)*100)) #relative abundance
#subset only OTUsAMF
tax= subset(tax, tax$phylum == "p__Glomeromycota", select= c("kingdom","phylum","class","order","family","genus"))
otu=otu[rownames(tax),]
#separar solo los otus de AMF
otu.n=otu.n[rownames(otu),]
colSums(otu.n)
#FILTER data with 0
otu.n=otu.n[,colSums(otu.n)>0]
colSums(otu.n)
otu.n=log10(otu.n+1)
otu.n=sqrt(otu.n)
#calculate dimensions
scaling=vegdist(t(otu.n), method = "bray", binary = T) #calculate distance
#scaling2=isoMDS(scaling)  ; scaling2$stress #create NMDS 
##scaling2=metaMDS(otu.n, distance = "bray")
scaling2=metaMDS(t(scaling), distance = "bray") ; scaling2$stress
#scaling2=monoMDS(scaling)
scaling3=data.frame(scaling2$points) #select cordenates
alpha=alpha[rownames(scaling3),]
scaling3=cbind(scaling3,alpha)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/")
pdf(file = "AMFNMDS.plantcompartment.pdf")
AMF=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=Specie))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=5, aes(shape=Specie))+
  scale_shape_manual(values = c(15:19))+
  scale_colour_manual(values = c("darksalmon","tan3","navajowhite4", "sienna", "blue"))+
  labs(x = "NMDS1",y = "NMDS2")
print(AMF)
dev.off()
####ANOTHER NDMS
ggplot(data=scalingasal, aes(x=MDS1, y=MDS2, colour=Specie))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=5, aes(shape=Location))+
  scale_shape_manual(values = c(16:17))+
  scale_colour_manual(values = c("darksalmon","tan3","navajowhite4", "sienna", "blue"))+
  labs(x = "NMDS1",y = "NMDS2")
print(AMF)
dev.off()

###PERMANOVA
p=t(otu.n)
set.seed(173612)
scaling=vegdist(p, method = "bray", binary = T) 
adonis2(scaling ~ Season, data = alpha, permutations = 1000)
adonis2(scaling ~ Sample, data = alpha, permutations = 1000)
adonis2(scaling ~ Location, data = alpha, permutations = 1000)
adonis2(scaling ~ Specie, data = alpha, permutations = 1000)
adonis2(scaling ~ Season*Sample*Location*Specie, data = alpha, permutations = 1000)


#Another step we use before it is tranform the data to "log"
#transponer la base: 
#otu3 = otu.n with the subset
datos<-t(otu3)
datos<-datos[rowSums(datos) > 0,]
rowSums(datos)
#transformar los datos con distancia "log": 
datostransformados <- decostand(datos, method = "log")