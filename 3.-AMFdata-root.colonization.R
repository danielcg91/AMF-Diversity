setwd <- ("~/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL")
#Load excel
#install.packages("readxl")
#Load Packages
library(readxl)
library(Rmisc)
library(ggplot2)
library(labdsv)
library(FSA)

#theme
tema=theme(axis.text.x = element_text(color="black",size=14, angle=0,hjust=0.5,vjust=0),
           axis.text.y = element_text(color="black",size=14),
           axis.title = element_text(color="black",size=14),
           legend.text = element_text(color = "black",size=14),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=10, color="black",face="bold"),
           strip.text.y = element_text(size=10, color="black",face="bold"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

#Load data
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization.xlsx"
#Rootcol= read_excel(ruta_excel)
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization2.xlsx"
Rootcol2= read_excel(ruta_excel2)
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
pc=c("soil", "rhizosphere", "root endosphere", "spores", "T0", "T4", "T8", "T12", 
     "Ateq", "Asal", "Mgeo")

###Barplots ROOT COLONIZATION
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
mean.r=summarySE(data = Rootcol2, measurevar = "root.colonization", groupvars = c("plant.species", "treatment"))
mean.r$plant.species=factor(mean.r$plant.species, level = pc)
mean.r$treatment=factor(mean.r$treatment, level = pc)
#TREATMENT-TIME
pdf(file = "Rootcolonization by time.pdf")
z=ggplot(data=mean.r, aes(x=treatment, y=root.colonization, fill= plant.species))+
  geom_col(position = "dodge")+
  #geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=root.colonization-sd, ymax=root.colonization+sd),width=0.4,position=position_dodge(0.9))+
  scale_fill_manual(values = c("gray80","gray50", "gray30"))+tema
z
dev.off()
#PLANT.SPECIES
pdf(file = "Rootcolonization by plant.pdf")
a=ggplot(data=mean.r, aes(x=plant.species, y=root.colonization, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=root.colonization-sd, ymax=root.colonization+sd),width=0.4,position=position_dodge(0.9))+
  scale_fill_manual(values = c("gray80","gray50", "gray30"))+tema
a
dev.off()

#Shapiro wilk test
shapiro.test(Rootcol2$root.colonization)
#shapiro.test(mean.rT8$root.colonization)

mean.r=summarySE(data = Rootcol2, measurevar = "root.colonization", groupvars = c("treatment", "plant.species"))

###KRUSKAL WALLIS. TEST by time
#T4 by plant
mean.rT4=subset(Rootcol2, Rootcol2$treatment == "T4", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~plant.species, data = mean.rT4) ; dunnTest(root.colonization~plant.species, data = mean.rT4, method = "bonferroni")
#test=aov(root.colonization~plant.species, data = mean.rT4); print(test) ; TukeyHSD(test)
#T8 bt plant
mean.rT8=subset(Rootcol2, Rootcol2$treatment == "T8", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~plant.species, data = mean.rT8) ; dunnTest(root.colonization~plant.species, data = mean.rT8, method = "bonferroni")
#test=aov(root.colonization~plant.species, data = mean.rT8); print(test) ; TukeyHSD(test)
#T12 by plant
mean.rT12=subset(Rootcol2, Rootcol2$treatment == "T12", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~plant.species, data = mean.rT12) ; dunnTest(root.colonization~plant.species, data = mean.rT12, method = "bonferroni")
#test=aov(root.colonization~plant.species, data = mean.rT12); print(test) ; TukeyHSD(test)


###KRUSKAL WALLIS. TEST by plant
#Ateq
mean.rateq=subset(Rootcol2, Rootcol2$plant.species == "Ateq", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~treatment, data = mean.rateq) ; dunnTest(root.colonization~treatment, data = mean.rateq, method = "bonferroni")
#test=aov(root.colonization~treatment, data = mean.rateq); print(test) ; TukeyHSD(test)
#Asal
mean.rasal=subset(Rootcol2, Rootcol2$plant.species == "Asal", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~treatment, data = mean.rasal) ; dunnTest(root.colonization~treatment, data = mean.rasal, method = "bonferroni")
#test=aov(root.colonization~treatment, data = mean.rasal); print(test) ; TukeyHSD(test)
#Mgeo
mean.rmgeo=subset(Rootcol2, Rootcol2$plant.species == "Mgeo", select = c("plant.species", "biological.replicates", "treatment", "root.colonization"))
kruskal.test(root.colonization~treatment, data = mean.rmgeo) ; dunnTest(root.colonization~treatment, data = mean.rmgeo, method = "bonferroni")
#test=aov(root.colonization~treatment, data = mean.rmgeo); print(test) ; TukeyHSD(test)



###Barplots NUMBER OF SPORES

numberspo= "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Number of spores.xlsx"
nspores= read_excel(numberspo)

###Barplots ROOT COLONIZATION
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
mean.r=summarySE(data = nspores, measurevar = "number.spores", groupvars = c("plant.species", "treatment"))
mean.r$plant.species=factor(mean.r$plant.species, level = pc)
mean.r$treatment=factor(mean.r$treatment, level = pc)

#TREATMENT-TIME
pdf(file = "numberofspores by time.pdf")
z=ggplot(data=mean.r, aes(x=treatment, y=number.spores, fill= plant.species))+
  geom_col(position = "dodge")+
  #geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=number.spores-sd, ymax=number.spores+sd),width=0.4,position=position_dodge(0.9))+
  scale_fill_manual(values = c("gray80","gray50", "gray30"))+tema
z
dev.off()


#PLANT SPECIES-TIME
pdf(file = "numberofspores by plant.pdf")
a=ggplot(data=mean.r, aes(x=plant.species, y=number.spores, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat="identity", fill=c("lemonchiffon4","tan3","navajowhite4", "sienna"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=number.spores-sd, ymax=number.spores+sd),width=0.4,position=position_dodge(0.9))+
  scale_fill_manual(values = c("gray80","gray50", "gray30"))+tema
a
dev.off()


#Shapiro wilk test
shapiro.test(nspores$number.spores)

#ANOVA BY TREATMENT AND PLANT

#Ateq
mean.rateq=subset(nspores, nspores$plant.species == "Ateq", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~treatment, data = mean.rateq); print(test) ; TukeyHSD(test)
#Asal
mean.rasal=subset(nspores, nspores$plant.species == "Asal", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~treatment, data = mean.rasal); print(test) ; TukeyHSD(test)
#Mgeo
mean.rmgeo=subset(nspores, nspores$plant.species == "Mgeo", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~treatment, data = mean.rmgeo); print(test) ; TukeyHSD(test)

##ANOVA BY TIME 
#T4
mean.rT4=subset(nspores, nspores$treatment == "T4", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~plant.species, data = mean.rT4); print(test) ; TukeyHSD(test)
#T8
mean.rT8=subset(nspores, nspores$treatment == "T8", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~plant.species, data = mean.rT8); print(test) ; TukeyHSD(test)
#T12
mean.rT12=subset(nspores, nspores$treatment == "T12", select = c("plant.species", "biological.replicates", "treatment", "number.spores"))
test=aov(number.spores~plant.species, data = mean.rT12); print(test) ; TukeyHSD(test)



######DORUGTH DATA

setwd <- ("~/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS IRVING SEQUIA/Colonización raices_1.xlsx")
#Load excel
#install.packages("readxl")
#Load Packages
library(readxl)
library(Rmisc)
library(ggplot2)
library(labdsv)
library(FSA)

#theme
tema=theme(axis.text.x = element_text(color="black",size=14, angle=0,hjust=0.5,vjust=0),
           axis.text.y = element_text(color="black",size=14),
           axis.title = element_text(color="black",size=14),
           legend.text = element_text(color = "black",size=14),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=10, color="black",face="bold"),
           strip.text.y = element_text(size=10, color="black",face="bold"),
           panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

#Load data
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS IRVING SEQUIA/Colonización raices_2.xlsx"
Rootcol2= read_excel(ruta_excel2)
pc=c("WW", "LD", "MD", "SD", "Mock", "AMF")
Rootcol2= subset(Rootcol2, Rootcol2$rootcolonization != "0")

###Barplots ROOT COLONIZATION
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
mean.r=summarySE(data = Rootcol2, measurevar = "rootcolonization", groupvars = c("regimen", "treatment"))
mean.r$regimen=factor(mean.r$regimen, level = pc)
mean.r$treatment=factor(mean.r$treatment, level = pc)

pdf(file = "Rootcolonization-drought.pdf")
z=ggplot(data=mean.r, aes(x=regimen, y=rootcolonization))+
  #geom_col(position = "dodge")+
  geom_bar(stat="identity", fill=c("gray90","gray70","gray50", "gray30"), colour="black",size=1)+
  geom_errorbar(size=1,aes(ymin=rootcolonization-sd, ymax=rootcolonization+sd),width=0.4,position=position_dodge(0.9))+
  scale_fill_manual(values = c("gray90","gray70","gray50", "gray30"))+tema
z
dev.off()

#Shapiro wilk test
shapiro.test(Rootcol2$rootcolonization)

###KRUSKAL WALLIS. TEST by time
#T4 by plant
mean.rT4=subset(Rootcol2, Rootcol2$treatment == "AMF", select = c("rootcolonization","regimen", "treatment"))
kruskal.test(rootcolonization~regimen, data = mean.rT4) ; dunnTest(rootcolonization~regimen, data = mean.rT4, method = "bh")
test=aov(rootcolonization~regimen, data = mean.rT4); print(test) ; TukeyHSD(test)





