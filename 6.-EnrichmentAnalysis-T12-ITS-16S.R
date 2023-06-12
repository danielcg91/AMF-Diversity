### Performs OTU enrichment analysis between treatments and plots ###

#Directories
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/")
its=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agaveITS/ITS2_1/") 
s16=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")  

#Packages
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(Rmisc)
library(tidyr)
library(dplyr)
library(FSA)
library(ggtern)
library(ggrepel)
library(extrafont)

#Graphics
colores=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
nivel=c("soil", "rhizosphere", "root endosphere", "spores")
etiquetas=c("Agave.tequilana", "Agave.salmiana", "Myrtillocactus.geometrizans")

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
           legend.position = "right", legend.text = element_text(color = "black",size=14, family = "sans"), 
           legend.direction = "vertical", legend.title = element_text(color = "black",size=14, face = "bold", family = "sans"),
           legend.key.size = unit(0.4,"cm"))


# Set parameters
amplicon=c("AMF","16S-S")[2]
data=c("syncoms","daniel", "all")[2]  #what sub sample would you need
base=c( "uniteall","silva")[2]
t=c("domain","phylum","class","order","family","genus","otu.id")[4] #taxa level analysis
s=c("soil", "rhizosphere", "root endosphere", "spores")[1:4] #Analyze only one treatment or all
m=c("Agave.tequilana", "Agave.salmiana", "Myrtillocactus.geometrizans","soil")
test=c("kruskal")[1]

# Load metadata
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
  # Removing unwanted samples
  met= subset(met, met$treatment != "T4")
  met= subset(met, met$treatment != "T8")
  otu=otu[rownames(met)]
  # Sub sample reads 
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu)))) 
  #Checklist
  colSums(otu.s)
  #otu.r=t(t(otu.s)/colSums(otu.s)*100)
  otu.r=cbind(otu.s,tax)
  otu.t=gather(data=otu.r, key = "sample", value = "abs", colnames(otu.s))
  for (i in met$AIM.sample.id){
    otu.t[otu.t$sample == i, c( "plant.species","plant.compartment","treatment","AIM.sample.id" )] = met[i,c( "plant.species","plant.compartment","treatment","AIM.sample.id" )]
  }
  
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  #Subset to eliminate unclassified
  otu.t= subset(otu.t, otu.t$order != "Unclassified")
  otu.t= subset(otu.t, otu.t$phylum == "p:Glomeromycota", 
                select= c("domain","phylum","class","order","family","genus","otu.id", "sample", 
                          "abs","plant.species","plant.compartment","treatment","AIM.sample.id"))
  
  #Calculates relative abundance for each taxa (OTU) in each treatment
  otu.t0 = otu.t %>% group_by(sample,plant.species, plant.compartment) %>% summarise(absab = sum(abs))
  otu.t1 = otu.t %>% group_by(otu.t[,t],sample, plant.species, plant.compartment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$plant.species,sep = ".")==paste(otu.t0$plant.species,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  
  #Subset the data 
  otu.t1=otu.t1[otu.t1$sample %in% met[met$plant.compartment %in% s, "AIM.sample.id"],]
  #otu.t1=otu.t1[otu.t1$sample %in% met[met$plant.species %in% m, "AIM.sample.id"],]
  taxon=sort(unique(otu.t1$rank))
  

} else if ( amplicon == "16S-S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$AIM.sample.id
  met=met[order(met$AIM.sample.id),]
  met=met[order(met$AIM.sample.id),]
  met=met[colnames(otu),]
  # Removing unwanted samples
  met= subset(met, met$treatment != "T4")
  met= subset(met, met$treatment != "T8")
  otu=otu[rownames(met)]
  #ONLY SPORES OTUS
 #otueb=read.delim(file=paste(getwd(),paste("otu.table.EB",data,"txt",sep = "."),sep = "/"))
  # Sub sample reads 
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu)))) 
  #Checklist
  colSums(otu.s)
  #otu.r=t(t(otu.s)/colSums(otu.s)*100)
  otu.r=cbind(otu.s,tax)
  otu.t=gather(data=otu.r, key = "sample", value = "abs", colnames(otu.s))
  for (i in met$AIM.sample.id){
    otu.t[otu.t$sample == i, c( "plant.species","plant.compartment","treatment","AIM.sample.id" )] = met[i,c( "plant.species","plant.compartment","treatment","AIM.sample.id" )]
  }
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
 
  #Calculates relative abundance for each taxa (OTU) in each treatment
  otu.t0 = otu.t %>% group_by(sample,plant.species, plant.compartment) %>% summarise(absab = sum(abs))

  #Subset to eliminate unclassified
  #otu.t= subset(otu.t, otu.t$order != "Unclassified")
  otu.t= subset(otu.t, otu.t$family != "Unclassified")
  #contnue
  otu.t1 = otu.t %>% group_by(otu.t[,t],sample, plant.species, plant.compartment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$plant.species,sep = ".")==paste(otu.t0$plant.species,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  
  #Subset the data 
  otu.t1=otu.t1[otu.t1$sample %in% met[met$plant.compartment %in% s, "AIM.sample.id"],]
  #otu.t1=otu.t1[otu.t1$sample %in% met[met$plant.species %in% m, "AIM.sample.id"],]
  taxon=sort(unique(otu.t1$rank))
  
}

#Test
taxon=as.character(sort(unique(otu.t1$rank)))

if(test=="kruskal"){
kruskal=data.frame(row.names = taxon, rank=taxon, p.value=rep(1, length(taxon)[1]))

  for (id in taxon){
    sub1=otu.t1[otu.t1$rank==id,]
    kruskal[id,"p.value"]=kruskal.test(relab ~ plant.compartment, data = sub1)$p.value
    kruskal$p.value=format(kruskal$p.value, scientific = F)
  }

  if (sum(kruskal$p.value<=0.05, na.rm = T)>0){
    
    dunn=data.frame(row.names = taxon)
    
    for (h in taxon){
      sub1=otu.t1[otu.t1$rank==h,]
    
  
      
      if (sum(sub1$relab)>0){
        test1=dunnTest(relab ~ plant.compartment, data = sub1, method = "bh")
        p.value=test1$res$P.adj
        dunn[h,1:choose(length(unique(sub1$plant.compartment)),2)]=p.value
      } else if (sum(sub1$relab)==0){
        p.value=rep(1,choose(length(unique(sub1$plant.compartment)),2))
        dunn[h,1:choose(length(unique(sub1$plant.compartment)),2)]=p.value
      }
  }
}
      
  colnames(dunn)=test1$res$Comparison
  dunn[is.na(dunn)]=1
  pv=dunn
  pv=select(pv, "rhizosphere - spores","root endosphere - spores","soil - spores")
  #pv=dunn[,grep("Agave.",colnames(dunn))]
  colnames(pv)=gsub(" - spores","",colnames(pv))
  #colnames(pv)=gsub(" - spores","",colnames(pv))
  #colnames(pv)=gsub(" - Agave.salmiana - ","",colnames(pv))
  
  for (h in colnames(pv)){
    print(h)
    print(sum(pv[,h]<=0.05))
  }
}

  
#Calculate log fold for plotting 
otu.d=summarySE(otu.t1, measurevar = "relab", groupvars = c("rank","plant.compartment")) #mean abundance of each OTU
en=data.frame(row.names=taxon)
for (i in colnames(pv)){
  logfold=log2(otu.d[otu.d$plant.compartment==i,"relab"]/otu.d[otu.d$plant.compartment=="spores","relab"])
  en=cbind(en, logfold)
}
colnames(en)=colnames(pv)
en=en[rownames(pv), colnames(pv)]
en[is.infinite(as.matrix(en))]=0
en[is.nan(as.matrix(en))]=0

#Select any colors 
colores=c("#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#F56AB0")

#Select the folder to save the figures
if (amplicon == "AMF"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURASS-ITS/")
} else if ( amplicon == "16S-S"){
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/Figuras/FIGURASS-16S/")
}  


# Generates a plot per each comparison 

for (w in colnames(pv)){
  
  # threshold for enrichment
  w="rhizosphere"
  #w="root endosphere"
  thf=1.5
  thp=0.05
  lab="order"

  #ab=otu.d[otu.d$plant.species == "Agave.salmiana",]
  ab=otu.d[otu.d$plant.compartment == w,]
  plot.data=data.frame(row.names = rownames(pv), pvalue=pv[,w], logf=en[,w], relab=ab$relab)

  
  for (i in rownames(plot.data)){
    plot.data[i, c("domain","phylum","class","order","family","genus","otu.id")]=unique(tax[tax[,t] %in% i, c("domain","phylum","class","order","family","genus","otu.id")])
  }


sig=plot.data[plot.data$pvalue > 0.05, "order"]
  
 
  #plot.data$order="differential"
  plot.data[plot.data$order %in% sig,"order"] ="no differencial"
  #plot.data$order=factor(plot.data$order, levels = c(niveles,"Unclassified","Low.abundant"))
  niveles=unique(plot.data[plot.data$order %in% plot.data$order, "order"])
  col=data.frame(color=c(colores[1:length(niveles)],"grey50"),taxa=c(niveles,"no differencial"))

  plot.data$label=plot.data[,lab]
  #plot.data[plot.data$pvalue>thp | abs(plot.data$logf)<thf,"label"]="_"
  
}
  
pdf(paste(amplicon, t,w, file = "enrichment.pdf"),  colormodel = "srgb", width = 11, height = 8.5, compress = F)
p=ggplot(data=plot.data, aes(x=logf, y=-log10(pvalue), color=order))+
    geom_point(size=4)+
    scale_color_manual(values=c(col[col$taxa %in% unique(plot.data$order),"color"]))+
    geom_text_repel(aes(label=label),size=3.5, color="black", max.overlaps = 42, force_pull = 2, force = 2)+
    geom_hline(yintercept = -log10(thp), linetype=2)+
    geom_vline(xintercept = c(-thf,thf), linetype=2)+#tema+
    ggtitle(paste("spores","vs",w))+
    scale_x_continuous(limits = c(-12,5), breaks = c(seq(-10,5,2.5)))
p
dev.off()

pdf(paste(amplicon, t,w, file = "enrichment-1.pdf"),  colormodel = "srgb", width = 11, height = 8.5, compress = F)
p=ggplot(data=plot.data, aes(x=logf, y=-log10(pvalue), color=order))+
  geom_point(size=4)+
  scale_color_manual(values=c(col[col$taxa %in% unique(plot.data$order),"color"]))+
  geom_text_repel(aes(label=label),size=3.5, color="black", max.overlaps = 42, force_pull = 2, force = 2)+
  geom_hline(yintercept = -log10(thp), linetype=2)+
  geom_vline(xintercept = c(-thf,thf), linetype=2)+tema+
  ggtitle(paste("spores","vs",w))+
  scale_x_continuous(limits = c(-12,5), breaks = c(seq(-10,5,2.5)))+
  theme(legend.position = "none")
p
dev.off()
  