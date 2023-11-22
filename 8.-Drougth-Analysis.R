setwd <- ("~/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA")
#Load excel
#Load Packages
library(readxl)
library(Rmisc)
library(ggplot2)
library(labdsv)
library(FSA)
library(dplyr)
library(tidyverse)
library(rstatix)


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

###LOAD DATA PLANT TRAITS
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA/datatraitsdrought.xlsx"
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA/RWC.xlsx"
#Rootcol= read_excel(ruta_excel)
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization2.xlsx"
Rootcol2= read_excel(ruta_excel2)
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
pc=c("WW", "LD", "MD", "SD", "Mock", "AMF")
colnames(Rootcol2)

###SUMMARY OF DATA
{
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
##Fresh.weight
mean.fw=summarySE(data = Rootcol2, measurevar = "fresh.weight", groupvars = c("regimen", "treatment"))
mean.fw$regimen=factor(mean.fw$regimen, level = pc)
mean.fw$treatment=factor(mean.fw$treatment, level = pc)
##Dry.weight
Rootcol3=Rootcol2[!(is.na(Rootcol2$dry.weight)), ]
mean.dw=summarySE(data = Rootcol3, measurevar = "dry.weight", groupvars = c("regimen", "treatment"))
mean.dw$regimen=factor(mean.dw$regimen, level = pc)
mean.dw$treatment=factor(mean.dw$treatment, level = pc)
##RWC
##Dry.weight
Rootcol4=Rootcol2[!(is.na(Rootcol2$RWC)), ]
mean.rwc=summarySE(data = Rootcol3, measurevar = "RWC", groupvars = c("regimen", "treatment"))
mean.rwc$regimen=factor(mean.rwc$regimen, level = pc)
mean.rwc$treatment=factor(mean.rwc$treatment, level = pc)
##foliar.area
#mean.fa=summarySE(data = Rootcol2, measurevar = "foliar.area", groupvars = c("regimen", "treatment"))
#mean.fa$regimen=factor(mean.fa$regimen, level = pc)
#mean.fa$treatment=factor(mean.fa$treatment, level = pc)
##leaf.length
mean.ll=summarySE(data = Rootcol2, measurevar = "leaf.length", groupvars = c("regimen", "treatment"))
mean.ll$regimen=factor(mean.ll$regimen, level = pc)
mean.ll$treatment=factor(mean.ll$treatment, level = pc)
##Med.Num.roots
mean.menr=summarySE(data = Rootcol2, measurevar = "Med.Num.roots" , groupvars = c("regimen", "treatment"))
mean.menr$regimen=factor(mean.menr$regimen, level = pc)
mean.menr$treatment=factor(mean.menr$treatment, level = pc)
##Max.Num.roots
mean.maxnr=summarySE(data = Rootcol2, measurevar = "Max.Num.roots" , groupvars = c("regimen", "treatment"))
mean.maxnr$regimen=factor(mean.maxnr$regimen, level = pc)
mean.maxnr$treatment=factor(mean.maxnr$treatment, level = pc)
##Total.Root.Length
mean.trl=summarySE(data = Rootcol2, measurevar = "Total.Root.Length" , groupvars = c("regimen", "treatment"))
mean.trl$regimen=factor(mean.trl$regimen, level = pc)
mean.trl$treatment=factor(mean.trl$treatment, level = pc)
##Depth
mean.de=summarySE(data = Rootcol2, measurevar = "Depth", groupvars = c("regimen", "treatment"))
mean.de$regimen=factor(mean.de$regimen, level = pc)
mean.de$treatment=factor(mean.de$treatment, level = pc)
##Maximum.Width
mean.mw=summarySE(data = Rootcol2, measurevar = "Maximum.Width", groupvars = c("regimen", "treatment"))
mean.mw$regimen=factor(mean.mw$regimen, level = pc)
mean.mw$treatment=factor(mean.mw$treatment, level = pc)
}

######BARPLOTS OF PLANT TRAITS
#FRESH.WEIGTH
pdf(file = "fresh.weight-drought.pdf")
fw=ggplot(data=mean.fw, aes(x=regimen, y=fresh.weight, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=fresh.weight-sd, ymax=fresh.weight+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
fw
dev.off()
#DRY.WEIGHT
pdf(file = "Dry.weight-drought.pdf")
dw=ggplot(data=mean.dw, aes(x=regimen, y=dry.weight, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=dry.weight-sd, ymax=dry.weight+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
dw
dev.off()
#RWC
pdf(file = "RWC-drought.pdf")
dw=ggplot(data=mean.rwc, aes(x=regimen, y=RWC, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=RWC-sd, ymax=RWC+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
dw
dev.off()
##leaf.length
pdf(file = "leaf.length-drought.pdf")
ll=ggplot(data=mean.ll, aes(x=regimen, y=leaf.length, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=leaf.length-sd, ymax=leaf.length+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
ll
dev.off()
##Med.Num.roots
pdf(file = "Med.Num.roots-drought.pdf")
menr=ggplot(data=mean.menr, aes(x=regimen, y=Med.Num.roots, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Med.Num.roots-sd, ymax=Med.Num.roots+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
menr
dev.off()
##Max.Num.roots
pdf(file = "Max.Num.roots-drought.pdf")
maxnr=ggplot(data=mean.maxnr, aes(x=regimen, y=Max.Num.roots, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Max.Num.roots-sd, ymax=Max.Num.roots+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
maxnr
dev.off()
##Total.Root.Length
pdf(file = "Total.Root.Length-drought.pdf")
trl=ggplot(data=mean.trl, aes(x=regimen, y=Total.Root.Length, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Total.Root.Length-sd, ymax=Total.Root.Length+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
trl
dev.off()
##Depth
pdf(file = "Depth-drought.pdf")
de=ggplot(data=mean.de, aes(x=regimen, y=Depth, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Depth-sd, ymax=Depth+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
de
dev.off()
##Maximum.Width
pdf(file = "Maximum.Width-drought.pdf")
mw=ggplot(data=mean.mw, aes(x=regimen, y=Maximum.Width, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Maximum.Width-sd, ymax=Maximum.Width+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
mw
dev.off()


###STATISTICAL ANALISYS
##
#kruskal.test(fresh.weigth~regimen, data = Rootcol2) ; dunnTest(fresh.weigth~regimen, data = Rootcol2, method = "bh")
#test=aov(fresh.weight~regimen * treatment, data = Rootcol2); print(test) ; TukeyHSD(test)
#summary(test)

##Fresh weight
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(fresh.weight~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Dry.weight
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(dry.weight~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#RWC
stat.test=Rootcol3 %>%
  group_by(regimen) %>%
  t_test(RWC~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#leaf.length 
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(leaf.length~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Med.Num.roots
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Med.Num.roots~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Max.Num.roots
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Max.Num.roots~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Total.Root.Length
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Total.Root.Length~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Depth
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Depth~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Maximum.Width  
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Maximum.Width~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test



###LOAD DATA PIGMENTS
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA/datapigmentdrought.xlsx"
#Rootcol= read_excel(ruta_excel)
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization2.xlsx"
Rootcol2= read_excel(ruta_excel2)
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
pc=c("WW", "LD", "MD", "SD", "Mock", "AMF")

###SUMMARY OF DATA
{
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
  ##Chla
  mean.chla=summarySE(data = Rootcol2, measurevar = "Chla", groupvars = c("regimen", "treatment"))
  mean.chla$regimen=factor(mean.chla$regimen, level = pc)
  mean.chla$treatment=factor(mean.chla$treatment, level = pc)
  ##Chlb
  mean.chlb=summarySE(data = Rootcol2, measurevar = "Chlb", groupvars = c("regimen", "treatment"))
  mean.chlb$regimen=factor(mean.chlb$regimen, level = pc)
  mean.chlb$treatment=factor(mean.chlb$treatment, level = pc)
  ##Chl_tot
  mean.chltot=summarySE(data = Rootcol2, measurevar = "Chl_tot", groupvars = c("regimen", "treatment"))
  mean.chltot$regimen=factor(mean.chltot$regimen, level = pc)
  mean.chltot$treatment=factor(mean.chltot$treatment, level = pc)
  ##Carotenoids
  mean.car=summarySE(data = Rootcol2, measurevar = "Carotenoids", groupvars = c("regimen", "treatment"))
  mean.car$regimen=factor(mean.car$regimen, level = pc)
  mean.car$treatment=factor(mean.car$treatment, level = pc)
  #Radio_Chla/Chlab
  mean.chab=summarySE(data = Rootcol2, measurevar = "Radio_Chla_Chlb", groupvars = c("regimen", "treatment"))
  mean.chab$regimen=factor(mean.chab$regimen, level = pc)
  mean.chab$treatment=factor(mean.chab$treatment, level = pc)
  ##Radio_Chltot/Carotenoides
  mean.chcar=summarySE(data = Rootcol2, measurevar = "Radio_Chltot_Carotenoides", groupvars = c("regimen", "treatment"))
  mean.chcar$regimen=factor(mean.chcar$regimen, level = pc)
  mean.chcar$treatment=factor(mean.chcar$treatment, level = pc)
}


######BARPLOTS OF PIGMENTS
###Chla
pdf(file = "Chla-drought.pdf")
chla=ggplot(data=mean.chla, aes(x=regimen, y=Chla, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Chla-sd, ymax=Chla+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
chla
dev.off()
###Chlb
pdf(file = "Chalb-drought.pdf")
chlb=ggplot(data=mean.chlb, aes(x=regimen, y=Chlb, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Chlb-sd, ymax=Chlb+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
chlb
dev.off()
##Chl_tot
pdf(file = "Chl_tot-drought.pdf")
chtot=ggplot(data=mean.chltot, aes(x=regimen, y=Chl_tot, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Chl_tot-sd, ymax=Chl_tot+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
chtot
dev.off()
##Carotenoids
pdf(file = "Carotenoids-drought.pdf")
car=ggplot(data=mean.car, aes(x=regimen, y=Carotenoids, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Carotenoids-sd, ymax=Carotenoids+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
car
dev.off()
##Radio_Chla_Chlb
pdf(file = "Radio_Chla_Chlb-drought.pdf")
rachab=ggplot(data=mean.chab, aes(x=regimen, y=Radio_Chla_Chlb, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Radio_Chla_Chlb-sd, ymax=Radio_Chla_Chlb+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
rachab
dev.off()
##Radio_Chltot_Carotenoides
pdf(file = "Radio_Chltot_Carotenoides-drought.pdf")
rchcar=ggplot(data=mean.chcar, aes(x=regimen, y=Radio_Chltot_Carotenoides, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=Radio_Chltot_Carotenoides-sd, ymax=Radio_Chltot_Carotenoides+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
rchcar
dev.off()

###STATISTICAL ANALISYS
##
#kruskal.test(fresh.weigth~regimen, data = Rootcol2) ; dunnTest(fresh.weigth~regimen, data = Rootcol2, method = "bh")
#test=aov(fresh.weight~regimen * treatment, data = Rootcol2); print(test) ; TukeyHSD(test)
#summary(test)

##Chla
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Chla~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Chlb
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Chlb~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Chl_tot
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Chl_tot~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Carotenoids
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Carotenoids~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
##Radio_Chla_Chlb
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Radio_Chla_Chlb~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
#Radio_Chltot_Carotenoides
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(Radio_Chltot_Carotenoides~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


###LOAD DATA PROLINE
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA/dataproline.xlsx"
#Rootcol= read_excel(ruta_excel)
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization2.xlsx"
Rootcol2= read_excel(ruta_excel2)
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
pc=c("WW", "LD", "MD", "SD", "Mock", "AMF")


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
##prolin
mean.pro=summarySE(data = Rootcol2, measurevar = "data", groupvars = c("regimen", "treatment"))
mean.pro$regimen=factor(mean.pro$regimen, level = pc)
mean.pro$treatment=factor(mean.pro$treatment, level = pc)

##proline
pdf(file = "proline-drought.pdf")
chtot=ggplot(data=mean.pro, aes(x=regimen, y=data, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=data-sd, ymax=data+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
chtot
dev.off()

#Proline
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(proline~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


###LOAD DATA FOLIAR AREA
ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/DATOS-SEQUIA/foliar_area.xlsx"
#Rootcol= read_excel(ruta_excel)
#ruta_excel2 = "/Users/danielchavez/OneDrive - CINVESTAV/Doctorado AMF/Resultados/AMF-DANIEL/Root Colonization2.xlsx"
Rootcol2= read_excel(ruta_excel2)
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20"))
pc=c("WW", "LD", "MD", "SD", "Mock", "AMF")


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
##foliar_area
mean.fa=summarySE(data = Rootcol2, measurevar = "foliar_area", groupvars = c("regimen", "treatment"))
mean.fa$regimen=factor(mean.fa$regimen, level = pc)
mean.fa$treatment=factor(mean.fa$treatment, level = pc)

##foliar_area
pdf(file = "foliar_area-drought.pdf")
chtot=ggplot(data=mean.fa, aes(x=regimen, y=foliar_area, fill= treatment))+
  geom_col(position = "dodge")+
  #geom_bar(stat = "identity", fill=c("gray90", "gray70"), colour="black",size=1)+
  geom_errorbar(size=0.5,aes(ymin=foliar_area-sd, ymax=foliar_area+sd),width=0.2,position=position_dodge(0.9))+
  scale_fill_manual(values=c("gray90","gray50", "gray30"))+tema
chtot
dev.off()
#foliar_area
stat.test=Rootcol2 %>%
  group_by(regimen) %>%
  t_test(foliar_area~treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

