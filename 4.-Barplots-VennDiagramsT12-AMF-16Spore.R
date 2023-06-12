#Load all necessary libraries
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(scales)
library(phyloseq)
library(tidyr)
library(dplyr)
library(plyr)
library(ggVennDiagram)
library(RColorBrewer)

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
amplicon=c("AMF","16S-S")[2]
base=c("uniteall", "silva")[2]
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F","#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","dimgrey"))
pc=c("soil", "rhizosphere", "root endosphere", "spores", "T0", "T4", "T8", "T12", 
     "Agave.tequilana", "Agave.salmiana", "Myrtillocactus.geometrizans")

#Upload data R Base
if (amplicon == "AMF"){
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
  
  ####Total reads by PLANT SPECIE AND PLANT COMPARTMENT
  t=c("domain","phylum","class","order","family","genus","otu.id")[4]
  ##Sumarise data
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
  #Checklist
  colSums(otu.s)
  #sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]
  otu.s=as.data.frame(cbind(otu.s,tax))
  otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu))
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  for (i in met$AIM.sample.id){
    otu.t[otu.t$sample == i, c( "plant.species","plant.compartment","treatment","AIM.sample.id" )] = met[i,c( "plant.species","plant.compartment","treatment","AIM.sample.id" )]
  }
  otu.t0 = otu.t %>% group_by(plant.species,plant.compartment) %>% summarise(absab = sum(abs))
#subset only glomeromycota
  otu.t= subset(otu.t, tax$phylum =="p:Glomeromycota", 
  select= c("domain","phylum","class","order","family","genus","otu.id", "sample","abs","plant.species","plant.compartment","treatment","AIM.sample.id"))
  otu.t1 = otu.t %>% group_by(otu.t[,t],plant.species,plant.compartment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$plant.species,otu.t1$plant.compartment,sep = ".")==paste(otu.t0$plant.species,otu.t0$plant.compartment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  #Remove low abundant taxa
  orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab))
  orden1 = otu.t1 %>% group_by(rank, plant.compartment) %>% summarise(relab = mean(relab)) 
  low = orden[orden$relab < 0, "rank"]
  otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "other"
  n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
  #Reorder taxa if needed
  orden=orden[order(orden$relab,decreasing = T),]
  niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
  otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "other"))
  otu.t1$plant.species=factor(otu.t1$plant.species, level = pc)
  otu.t1$plant.compartment=factor(otu.t1$plant.compartment, level = pc)
  
} else if ( amplicon == "16S-S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  #ONLY SPORES OTUS
  otueb=read.delim(file=paste(getwd(),paste("otu.table.EB",data,"txt",sep = "."),sep = "/"))
  colnames(otueb)=gsub("X","",colnames(otueb))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$AIM.sample.id
  met=met[order(met$AIM.sample.id),]
  met=met[order(met$AIM.sample.id),]
  met=met[colnames(otu),]
  #same lenght for OTUS in spores
  tax=tax[rownames(otu),]
  taxeb=tax[rownames(otueb),]
  #Subset of data to process only at 12 months and soil
  met= subset(met, met$treatment != "T4")
  met= subset(met, met$treatment != "T8")
  otu=otu[rownames(met)]
  ####Total reads by PLANT SPECIE AND PLANT COMPARTMENT
  t=c("domain","phylum","class","order","family","genus","otu.id")[5]
  ##Sumarise data
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
  #Checklist
  colSums(otu.s)
  #sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]
  otu.s=as.data.frame(cbind(otu.s,tax))
  otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu))
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  for (i in met$AIM.sample.id){
    otu.t[otu.t$sample == i, c( "plant.species","plant.compartment","treatment","AIM.sample.id" )] = met[i,c( "plant.species","plant.compartment","treatment","AIM.sample.id" )]
  }
  otu.t0 = otu.t %>% group_by(plant.species,plant.compartment) %>% summarise(absab = sum(abs))
  #Subset only amf-spores bacteria
  rows = rownames(taxeb)
  otu.t= otu.t[(otu.t$otu.id) %in% rows, ]
  #continue
  otu.t1 = otu.t %>% group_by(otu.t[,t],plant.species,plant.compartment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$plant.species,otu.t1$plant.compartment,sep = ".")==paste(otu.t0$plant.species,otu.t0$plant.compartment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  #Remove low abundant taxa
  orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab))
  order1 = otu.t1 %>% group_by(rank, plant.compartment) %>% summarise(relab = mean(relab))
  #genus
  #low = orden[orden$relab < 0.43, "rank"]
  #otu.id
  low = orden[orden$relab < 0, "rank"]
  otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "other"
  n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
  #Reorder taxa if needed
  orden=orden[order(orden$relab,decreasing = T),]
  niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
  otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "other"))
  otu.t1$plant.species=factor(otu.t1$plant.species, level = pc)
  otu.t1$plant.compartment=factor(otu.t1$plant.compartment, level = pc)
}


#Select the folder to save the figures
if (amplicon == "AMF"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURASS-ITS/")
} else if ( amplicon == "16S-S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURASS-16S/")
}  


#FIGURES withou labels by genus
pdf(paste(amplicon, t, file = "RA BY PC-1.1.pdf"),  colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=plant.species, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+tema+
  #facet_grid(~plant.compartment, scales = "free", space = "free_x")+
  facet_wrap(~plant.compartment, nrow = 4, ncol=4, scales = "free", strip.position = "top")+
  #scale_y_continuous(limits = c(0,0.8))+
  theme(legend.position = "none")
p
dev.off()

#FIGURES with labels by genus
pdf(paste(amplicon, t, file = "RA BY PC.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=plant.species, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+
  facet_grid(~plant.compartment, scales = "free", space = "free_x")+
  #facet_wrap(~plant.compartment, nrow = 4, ncol=4, scales = "free", strip.position = "top")+
  guides(fill=guide_legend(title = "genus", ncol = 1, keyheight = 0.7, default.unit="cm",
                           label.theme = element_text(size= 12, face = "plain", family = "sans")))+tema
p
dev.off()

#FIGURES withou labels by genus by plant species
pdf(paste(amplicon, t, file = "RA BY PS-1.pdf"),  colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=plant.compartment, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+tema+
  facet_grid(~plant.species, scales = "free", space = "free_x")+
  #facet_wrap(~plant.species, nrow = 4, ncol=4, scales = "free", strip.position = "top")+
  theme(legend.position = "none")
p
dev.off()

#FIGURES with labels by genus by plant species
pdf(paste(amplicon, t,file = "RA BY PS.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=plant.compartment, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+
  facet_grid(~plant.species, scales = "free", space = "free_x")+
  #facet_wrap(~plant.species, nrow = 4, ncol=4, scales = "free", strip.position = "top")+
  guides(fill=guide_legend(title = "genus", ncol = 1, keyheight = 0.7, default.unit="cm",
  label.theme = element_text(size= 12, face = "plain", family = "sans")))+tema
p
dev.off()


####SUBSET FOR Venn diagrams ONLY FUNGI OTUS IN SPORES#######
#At
ATspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank", "plant.compartment", "plant.species", "relab"))
ATspores= subset(ATspores, ATspores$plant.species =="Agave.tequilana", select = c("rank", "plant.compartment", "plant.species", "relab"))
ATspores=subset(ATspores, ATspores$relab >0.00)
ATspores= as.character(unique(ATspores$rank))
#As
ASspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank","plant.compartment", "plant.species","relab"))
ASspores= subset(ASspores, ASspores$plant.species =="Agave.salmiana", select = c("rank","plant.compartment", "plant.species","relab"))
ASspores=subset(ASspores, ASspores$relab >0.00)
ASspores= as.character(unique(ASspores$rank))
#Mg
Mgspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank","plant.compartment", "plant.species","relab"))
Mgspores= subset(Mgspores, Mgspores$plant.species =="Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species", "relab"))
Mgspores=subset(Mgspores, Mgspores$relab >0.00)
Mgspores= as.character(unique(Mgspores$rank))
#
setlist <- list(AT=ATspores, AS=ASspores, MG=Mgspores)
#Venn diagrams of AMF OTUS in SPORES
pdf(paste(amplicon, t, file = "DV-OTUS IN SPORES.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggVennDiagram(setlist, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"), 
                label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p
dev.off()


####SUBSET FOR Venn diagrams OF FUNGI OTUS IN ROOT ENDOSPHERE
#AT
ATRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
ATRE= subset(ATRE, ATRE$plant.species =="Agave.tequilana", select = c("rank","plant.compartment", "plant.species", "relab"))
ATRE=subset(ATRE, ATRE$relab >0.00)
ATRE= as.character(unique(ATRE$rank))
#AS
ASRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
ASRE= subset(ASRE, ASRE$plant.species =="Agave.salmiana", select = c("rank","plant.compartment", "plant.species", "relab"))
ASRE=subset(ASRE, ASRE$relab >0.00)
ASRE= as.character(unique(ASRE$rank))
#MG
MGRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
MGRE= subset(MGRE, MGRE$plant.species =="Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species", "relab"))
MGRE=subset(MGRE, MGRE$relab >0.00)
MGRE= as.character(unique(MGRE$rank))

setlist <- list(AT=ATRE, AS=ASRE, MG=MGRE)

#Venn diagram of FUNGI OTUS in Root endosphere
pdf(paste(amplicon, t, file = "DV-OTUS IN ROOT ENDOSPHERE.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggVennDiagram(setlist, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"), 
                label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p
dev.off()


####SUBSET fo Venn diagrams  BY PLANT COMPARTMENT##########
#SOIL
SOIL= subset(otu.t1, otu.t1$plant.compartment == "soil", select = c("rank", "plant.compartment", "plant.species", "relab"))
SOIL=subset(SOIL, SOIL$relab >0.00)
SOIL= as.character(unique(SOIL$rank))
#RE
RE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank", "plant.compartment", "plant.species", "relab"))
RE=subset(RE, RE$relab >0.00)
RE= as.character(unique(RE$rank))
#RHI
RHI= subset(otu.t1, otu.t1$plant.compartment == "rhizosphere", select = c("rank", "plant.compartment", "plant.species", "relab"))
RHI=subset(RHI, RHI$relab >0.00)
RHI= as.character(unique(RHI$rank))
#SPO
SPO= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank", "plant.compartment", "plant.species", "relab"))
SPO=subset(SPO, SPO$relab >0.00)
SPO= as.character(unique(SPO$rank))

setlist2 <- list(S=SOIL, RHI= RHI, RE=RE, SPO=SPO)

#Venn diagrams of AMF OTUS by plant compartment
pdf(paste(amplicon, t, file = "DV-OTUS IN PLANT COMPARTMENT.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p1=ggVennDiagram(setlist2, category.names = c("soil","rhizosphere", "root endosphere", "spores"), 
                 label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black", "black"))+
  theme(legend.position = "none")
p1
dev.off()


####SUBSET for Venn diagrams by PLANT SPECIE########
#AT
AT= subset(otu.t1, otu.t1$plant.species == "Agave.tequilana", select = c("rank","plant.compartment", "plant.species","relab"))
AT=subset(AT, AT$relab >0.00)
AT= as.character(unique(AT$rank))
#AS
AS= subset(otu.t1, otu.t1$plant.species == "Agave.salmiana", select = c("rank","plant.compartment", "plant.species","relab"))
AS=subset(AS, AS$relab >0.00)
AS= as.character(unique(AS$rank))
#MG
MG= subset(otu.t1, otu.t1$plant.species == "Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species","relab"))
MG=subset(MG, MG$relab >0.00)
MG= as.character(unique(MG$rank))

setlist3 <- list(AT, AS, MG)

#Venn diagrams of AMF OTUS in plant specie
pdf(paste(amplicon,t, file = "DV-OTUS IN PLANT SPECIE.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p2=ggVennDiagram(setlist3, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"), 
                 label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p2
dev.off()