setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")

library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(igraph)
library(tidyr)
library(dplyr)
library(phyloseq)
library(FSA)


#All features to graph the data
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5),
  #axis.text.x = element_text(color="black",size=12, angle=90,hjust=0,vjust=0.5),
  axis.text.y = element_text(color="black",size=12, vjust = 1.5),
  axis.title = element_text(color="black",size=12, face = "bold"),
  panel.border =element_rect(color = "black", fill = NA),#element_blank(),
  strip.text.x = element_text(size=12, color="black",face="bold"),
  strip.text.y = element_text(size=12, color="black",face="bold"),
  strip.placement = "outside", strip.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
  panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
  legend.position = "bottom", legend.text = element_text(color = "black",size=12), legend.direction = "horizontal",
  legend.title = element_text(color = "black",size=14, face = "bold"),
  legend.key.size = unit(0.3,"cm"))

pc1=c("At.Am.D.s.1","At.Am.R.s.1","At.Am.D.s.2","At.Am.R.s.2","At.Am.R.s.3","At.Pe.D.s.1","At.Pe.R.s.1","At.Pe.D.s.2","At.Pe.R.s.2","At.Pe.D.s.3",
      "At.Pe.R.s.3","At.Am.D.rzs.1","At.Am.R.rzs.1","At.Am.D.rzs.2","At.Am.R.rzs.2","At.Am.D.rzs.3","At.Am.R.rzs.3","At.Pe.D.rzs.1","At.Pe.R.rzs.1",
      "At.Pe.D.rzs.2","At.Pe.R.rzs.2","At.Pe.D.rzs.3","At.Pe.R.rzs.3","At.Am.D.rz.1","At.Am.R.rz.1","At.Am.D.rz.2","At.Am.R.rz.2" ,"At.Am.D.rz.3",
      "At.Am.R.rz.3","At.Pe.D.rz.1","At.Pe.R.rz.1","At.Pe.D.rz.2","At.Pe.R.rz.2","At.Pe.D.rz.3","At.Pe.R.rz.3","At.Am.D.re.1","At.Am.R.re.1","At.Am.D.re.2",
      "At.Am.R.re.2","At.Am.D.re.3","At.Am.R.re.3","At.Pe.R.re.1","At.Pe.D.re.2","At.Pe.R.re.2","At.Pe.R.re.3",
      "Ad.Ah.D.s.1","Ad.Ah.R.s.1","Ad.Ah.D.s.2","Ad.Ah.R.s.2","Ad.Ah.D.s.3","Ad.Ah.R.s.3","Ad.Br.D.s.1","Ad.Br.R.s.1","Ad.Br.D.s.2","Ad.Br.R.s.2","Ad.Br.D.s.3","Ad.Br.R.s.3",
      "Ad.Pf.D.s.1","Ad.Pf.R.s.1","Ad.Pf.D.s.2","Ad.Pf.R.s.2","Ad.Pf.D.s.3","Ad.Pf.R.s.3","Ad.Ah.D.rzs.1","Ad.Ah.R.rzs.1","Ad.Ah.D.rzs.2","Ad.Ah.R.rzs.2","Ad.Ah.D.rzs.3","Ad.Ah.R.rzs.3",
      "Ad.Br.D.rzs.1","Ad.Br.R.rzs.1","Ad.Br.D.rzs.2","Ad.Br.R.rzs.2","Ad.Br.D.rzs.3","Ad.Br.R.rzs.3","Ad.Pf.D.rzs.1","Ad.Pf.R.rzs.1","Ad.Pf.D.rzs.2","Ad.Pf.R.rzs.2","Ad.Pf.D.rzs.3",
      "Ad.Pf.R.rzs.3","Ad.Ah.D.rz.1","Ad.Ah.R.rz.1","Ad.Ah.R.rz.2","Ad.Ah.D.rz.3","Ad.Ah.R.rz.3","Ad.Br.D.rz.1","Ad.Br.R.rz.1","Ad.Br.D.rz.2","Ad.Br.R.rz.2","Ad.Br.D.rz.3","Ad.Br.R.rz.3",
      "Ad.Pf.D.rz.1","Ad.Pf.R.rz.1","Ad.Pf.D.rz.2","Ad.Pf.R.rz.2","Ad.Pf.D.rz.3","Ad.Pf.R.rz.3","Ad.Ah.D.re.1","Ad.Ah.R.re.1","Ad.Ah.D.re.2","Ad.Ah.R.re.2","Ad.Ah.D.re.3","Ad.Ah.R.re.3",
      "Ad.Br.D.re.1","Ad.Br.R.re.1","Ad.Br.D.re.2","Ad.Br.D.re.3","Ad.Br.R.re.3","Ad.Pf.D.re.1","Ad.Pf.R.re.1","Ad.Pf.R.re.2","Ad.Pf.D.re.3","Ad.Pf.R.re.3",
      "As.Ma.D.s.1","As.Ma.R.s.1","As.Ma.D.s.2","As.Ma.R.s.2","As.Ma.D.s.3","As.Ma.R.s.3","As.Sf.D.s.1","As.Sf.R.s.1","As.Sf.D.s.2","As.Sf.R.s.2","As.Sf.D.s.3",
      "As.Sf.R.s.3", "As.Ma.D.rzs.1","As.Ma.R.rzs.1","As.Ma.D.rzs.2","As.Ma.R.rzs.2","As.Ma.D.rzs.3","As.Ma.R.rzs.3","As.Sf.D.rzs.1","As.Sf.R.rzs.1","As.Sf.D.rzs.2",
      "As.Sf.R.rzs.2","As.Sf.D.rzs.3","As.Sf.R.rzs.3","As.Ma.D.rz.1","As.Ma.R.rz.1","As.Ma.D.rz.2","As.Ma.R.rz.2","As.Ma.D.rz.3","As.Ma.R.rz.3","As.Sf.D.rz.1",
      "As.Sf.R.rz.1","As.Sf.D.rz.2","As.Sf.R.rz.2","As.Sf.D.rz.3","As.Sf.R.rz.3","As.Ma.D.re.1","As.Ma.R.re.1","As.Ma.D.re.2","As.Ma.R.re.2","As.Ma.D.re.3","As.Ma.R.re.3",
      "As.Sf.D.re.1","As.Sf.R.re.1","As.Sf.D.re.2","As.Sf.R.re.2","As.Sf.D.re.3","As.Sf.R.re.3",
      "Mg.Ma.D.s.1","Mg.Sf.D.s.1","Mg.Ma.R.s.1","Mg.Sf.R.s.1", "Mg.Ma.D.rzs.1","Mg.Sf.D.rzs.1","Mg.Ma.R.rzs.1","Mg.Sf.R.rzs.1","Mg.Ma.D.rz.1","Mg.Sf.D.rz.1","Mg.Ma.R.rz.1",
      "Mg.Sf.R.rz.1","Mg.Ma.D.re.1","Mg.Sf.D.re.1","Mg.Ma.R.re.1","Mg.Sf.R.re.1",
      "Or.Ma.D.s.1","Or.Sf.D.s.1","Or.Ma.R.s.1","Or.Sf.R.s.1", "Or.Ma.D.rzs.1","Or.Sf.D.rzs.1","Or.Ma.R.rzs.1","Or.Sf.R.rzs.1","Or.Ma.D.rz.1","Or.Sf.D.rz.1","Or.Ma.R.rz.1",
      "Or.Sf.R.rz.1","Or.Ma.D.re.1","Or.Sf.D.re.1","Or.Ma.R.re.1","Or.Sf.R.re.1")

pc2=c("At.Am.D.s.1","At.Am.R.s.1","At.Am.D.s.2","At.Am.R.s.2","At.Am.R.s.3","At.Am.D.rzs.1","At.Am.R.rzs.1","At.Am.D.rzs.2","At.Am.R.rzs.2","At.Am.D.rzs.3","At.Am.R.rzs.3",
      "At.Am.D.rz.1","At.Am.R.rz.1","At.Am.D.rz.2","At.Am.R.rz.2","At.Am.D.rz.3","At.Am.R.rz.3","At.Am.D.re.1","At.Am.R.re.1","At.Am.D.re.2","At.Am.R.re.2","At.Am.D.re.3","At.Am.R.re.3",
      "At.Pe.D.s.1","At.Pe.R.s.1","At.Pe.D.s.2","At.Pe.R.s.2","At.Pe.D.s.3","At.Pe.R.s.3","At.Pe.D.rzs.1","At.Pe.R.rzs.1","At.Pe.D.rzs.2","At.Pe.R.rzs.2","At.Pe.D.rzs.3","At.Pe.R.rzs.3",
      "At.Pe.D.rz.1","At.Pe.R.rz.1","At.Pe.D.rz.2","At.Pe.R.rz.2","At.Pe.D.rz.3","At.Pe.R.rz.3","At.Pe.R.re.1","At.Pe.D.re.2","At.Pe.R.re.2","At.Pe.R.re.3",
      "Ad.Ah.D.s.1","Ad.Ah.R.s.1","Ad.Ah.D.s.2","Ad.Ah.R.s.2","Ad.Ah.D.s.3","Ad.Ah.R.s.3","Ad.Ah.D.rzs.1","Ad.Ah.R.rzs.1","Ad.Ah.D.rzs.2","Ad.Ah.R.rzs.2","Ad.Ah.D.rzs.3","Ad.Ah.R.rzs.3",
      "Ad.Ah.D.rz.1","Ad.Ah.R.rz.1","Ad.Ah.R.rz.2","Ad.Ah.D.rz.3","Ad.Ah.R.rz.3","Ad.Ah.D.re.1","Ad.Ah.R.re.1","Ad.Ah.D.re.2","Ad.Ah.R.re.2","Ad.Ah.D.re.3","Ad.Ah.R.re.3",
      "Ad.Br.D.s.1","Ad.Br.R.s.1","Ad.Br.D.s.2","Ad.Br.R.s.2","Ad.Br.D.s.3","Ad.Br.R.s.3","Ad.Br.D.rzs.1","Ad.Br.R.rzs.1","Ad.Br.D.rzs.2","Ad.Br.R.rzs.2","Ad.Br.D.rzs.3","Ad.Br.R.rzs.3",
      "Ad.Br.D.rz.1","Ad.Br.R.rz.1","Ad.Br.D.rz.2","Ad.Br.R.rz.2","Ad.Br.D.rz.3","Ad.Br.R.rz.3", "Ad.Br.D.re.1","Ad.Br.R.re.1","Ad.Br.D.re.2","Ad.Br.D.re.3","Ad.Br.R.re.3",
      "Ad.Pf.D.s.1","Ad.Pf.R.s.1","Ad.Pf.D.s.2","Ad.Pf.R.s.2","Ad.Pf.D.s.3","Ad.Pf.R.s.3","Ad.Pf.D.rzs.1","Ad.Pf.R.rzs.1","Ad.Pf.D.rzs.2","Ad.Pf.R.rzs.2","Ad.Pf.D.rzs.3","Ad.Pf.R.rzs.3",
      "Ad.Pf.D.rz.1","Ad.Pf.R.rz.1","Ad.Pf.D.rz.2","Ad.Pf.R.rz.2","Ad.Pf.D.rz.3","Ad.Pf.R.rz.3","Ad.Pf.D.re.1","Ad.Pf.R.re.1","Ad.Pf.R.re.2","Ad.Pf.D.re.3","Ad.Pf.R.re.3",
      "As.Ma.D.s.1","As.Ma.R.s.1","As.Ma.D.s.2","As.Ma.R.s.2","As.Ma.D.s.3","As.Ma.R.s.3","Mg.Ma.D.s.1","Mg.Ma.R.s.1","Or.Ma.D.s.1", "Or.Ma.R.s.1",
      "As.Ma.D.rzs.1","As.Ma.R.rzs.1","As.Ma.D.rzs.2","As.Ma.R.rzs.2","As.Ma.D.rzs.3","As.Ma.R.rzs.3","Mg.Ma.D.rzs.1","Mg.Ma.R.rzs.1","Or.Ma.D.rzs.1","Or.Ma.R.rzs.1",
      "As.Ma.D.rz.1","As.Ma.R.rz.1","As.Ma.D.rz.2","As.Ma.R.rz.2","As.Ma.D.rz.3","As.Ma.R.rz.3","As.Sf.D.rz.1","Mg.Ma.D.rz.1","Mg.Ma.R.rz.1","Or.Ma.D.rz.1","Or.Ma.R.rz.1",
      "As.Ma.D.re.1","As.Ma.R.re.1","As.Ma.D.re.2","As.Ma.R.re.2","As.Ma.D.re.3","As.Ma.R.re.3","Mg.Ma.D.re.1","Mg.Ma.R.re.1", "Or.Ma.D.re.1","Or.Ma.R.re.1",
      "As.Sf.D.s.1","As.Sf.R.s.1","As.Sf.D.s.2","As.Sf.R.s.2","As.Sf.D.s.3","As.Sf.R.s.3","Mg.Sf.D.s.1","Mg.Sf.R.s.1","Or.Sf.D.s.1","Or.Sf.R.s.1",
      "As.Sf.D.rzs.1","As.Sf.R.rzs.1","As.Sf.D.rzs.2","As.Sf.R.rzs.2","As.Sf.D.rzs.3","As.Sf.R.rzs.3","Mg.Sf.D.rzs.1","Mg.Sf.R.rzs.1","Or.Sf.D.rzs.1","Or.Sf.R.rzs.1",
      "As.Sf.R.rz.1","As.Sf.D.rz.2","As.Sf.R.rz.2","As.Sf.D.rz.3","As.Sf.R.rz.3","Mg.Sf.D.rz.1","Mg.Sf.R.rz.1","Or.Sf.D.rz.1","Or.Sf.R.rz.1",
      "As.Sf.D.re.1","As.Sf.R.re.1","As.Sf.D.re.2","As.Sf.R.re.2","As.Sf.D.re.3","As.Sf.R.re.3","Mg.Sf.D.re.1","Mg.Sf.R.re.1","Or.Sf.D.re.1","Or.Sf.R.re.1")

##Factors
{
ranks=c("domain","phylum","class","order","family","genus","otu.id")
data.dir=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "khaki","grey20"))
pc=c("soil", "leaf.endosphere", "phyllosphere", "rhizosphere", "roo.zone.soil" , "root.endosphere", "dry", "rainy", 
     "Amatitan","Penjamo","Agave.Hill","Boyd.Ridge","Pinyon.Flat", "Magueyal","SanFelipe",
     "Agave.tequilana","Agave.deserti","Agave.salmiana","M.geometrizans","O.robusta")
pc=c("SanFelipe", "Magueyal", "Boyd.Ridge", "Pinyon.Flat", "Agave.Hill")


## Upload data R base
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
meta=read.delim(paste(data.dir,"metadataAMF_AgaveCactus_measureable.txt",sep = "/"), header = T, stringsAsFactors = F)
rownames(meta)= meta$SampleID
tax=read.delim(paste(data.dir,"TaxaITS_agaves-cactus_measureable.txt",sep = "/"),header = T, stringsAsFactors = F)
tax$OTU.ID=paste("OTU",tax$OTU.ID,sep="_")
rownames(tax)= tax$OTU.ID; tax=tax[,-1]
#####Modified the data set eliminating the first letters and all data that was not classified put the name as Unclassified
tax$kingdom=gsub("k__", "", tax$kingdom)
tax$kingdom=gsub("__", "Unclassified", tax$kingdom)
tax$phylum=gsub("p__", "", tax$phylum)
tax$phylum=gsub("__", "Unclassified", tax$phylum)
tax$class=gsub("c__", "", tax$class)
tax$class=gsub("__", "Unclassified", tax$class)
tax$order=gsub("o__", "", tax$order)
tax$order=gsub("__", "Unclassified", tax$order)
tax$family=gsub("f__", "", tax$family)
tax$family=gsub("__", "Unclassified", tax$family)
tax$genus=gsub("g__", "", tax$genus)
tax$genus=gsub("__", "Unclassified", tax$genus)

#Rarefy otu table only AMF
set.seed(23171341)
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
#Checklist
colSums(otu.s)
{
print("Kingdom")
paste(round(sum(is.na(tax$domain))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$domain)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$domain)[!is.na(unique(tax$domain))]),"unique domains")

print("phyla")
paste(round(sum(is.na(tax$phylum))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$phylum)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$phylum)[!is.na(unique(tax$phylum))]),"unique phyla")

print("class")
paste(round(sum(is.na(tax$class))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$class)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$class)[!is.na(unique(tax$class))]),"unique classes")

print("order")
paste(round(sum(is.na(tax$order))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$order)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$order)[!is.na(unique(tax$order))]),"unique orders")

print("family")
paste(round(sum(is.na(tax$family))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$family)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$family)[!is.na(unique(tax$family))]),"unique families")

print("genus")
paste(round(sum(is.na(tax$genus))/dim(tax)[1]*100,1),"%","Unclassified OTUs")
paste(round((dim(tax)[1]-sum(is.na(tax$genus)))/dim(tax)[1]*100,1),"%","Classified OTUs")
paste(length(unique(tax$genus)[!is.na(unique(tax$genus))]),"unique genera")
}

#######BARPLOTS GENUS########
otu.r=t(t(otu.s)/colSums(otu.s)*100)
sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]
#subset only glomeromycota
tax= subset(tax, tax$phylum =="Glomeromycota", select= c("kingdom","phylum","class","order","family","genus"))
otu=otu[rownames(tax),]
otu.r=otu.r[rownames(tax),]
otu.r=otu.r[,colSums(otu.r)>0]
colSums(otu.r)
otu=otu[,colnames(otu.r)]
otu.r=cbind(otu.r,tax)
otu.t=gather(data=otu.r, key = "sample", value = "relab", colnames(otu))
#otu.t$genus[is.na(otu.t$domain)]="Unclassified"

for (i in meta$SampleID){
  otu.t[otu.t$sample == i, c("Sample","Season","Location","Specie")] = meta[i,c("Sample","Season","Location","Specie")]
}
}
##BY GENUS
#Remove low abundant taxa by genus
orden = otu.t %>% group_by(genus,sample) %>% summarise(relab = sum(relab))  %>% summarise (relab = mean(relab)) 
low = orden[orden$relab < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Refamily taxa if needed
orden=orden[order(orden$relab,decreasing = T),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
otu.t$genus=factor(otu.t$genus, level = c(niveles[niveles!="Unclassified"],"Unclassified", "low abundant"))
#otu.t$Specie=factor(otu.t$Specie, level = pc)
otu.t$sample=factor(otu.t$sample, level = pc2)


#PLOTS BY sample
pdf(file = "RA BY EACHSAMPLE-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=otu.t, aes(y=relab, x=sample, fill=genus))+
  geom_bar(stat="identity",width = 0.9)+
  #scale_fill_manual(values = my_pal(n))+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #facet_grid(Specie~Season, scales = "free_x", space = "free")+
  #facet_wrap(Season~Specie, nrow = 2, scales = "free",strip.position = "top")+
  labs(x = "Sample",y = "relative abundance (%)")+
  ggtitle("Relative abudance across samples")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()

##subset
otu.t.SF=subset(otu.t, otu.t$Location == "SanFelipe")
otu.t.SF$sample=factor(otu.t.SF$sample, level = pc2)

#PLOTS BY SANFELIPE
pdf(file = "RA BY EACHSAMPLE-genus", colormodel = "cmyk", width = 8.5, height = 11, compress = F)
a=ggplot(data=otu.t.SF, aes(y=relab, x=sample, fill=genus))+
  geom_bar(stat="identity",width = 0.9)+
  scale_fill_manual(values = my_pal(n))+
  #facet_grid(Specie~Season, scales = "free_x", space = "free")+
  #facet_wrap(Season~Specie, nrow = 2, scales = "free",strip.position = "top")+
  labs(x = "Sample",y = "relative abundance (%)")+
  ggtitle("Relative abudance across samples in the soil ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()



####DATOS CON RIQUEZA###
ranks=c("domain","phylum","class","order","family","genus","otu.id")
data.dir=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave-cactus2014/Datos-agaves-cactus/")
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "khaki","grey20"))
pc=c("soil", "leaf.endosphere", "phyllosphere", "roo.zone.soil","rhizosphere","root.endosphere", "dry", "rainy", 
     "Amatitan","Penjamo","Agave.Hill","Boyd.Ridge","Pinyon.Flat","Magueyal","SanFelipe",
     "Agave.tequilana","Agave.salmiana","Agave.deserti","M.geometrizans","O.robusta")
pc=c("SanFelipe", "Magueyal", "Boyd.Ridge", "Pinyon.Flat", "Agave.Hill", "Amatitan","Penjamo")

## Upload data R base
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
meta=read.delim(paste(data.dir,"metadataAMF_AgaveCactus_measureable.txt",sep = "/"), header = T, stringsAsFactors = F)
rownames(meta)= meta$SampleID
tax=read.delim(paste(data.dir,"TaxaITS_agaves-cactus_measureable.txt",sep = "/"),header = T, stringsAsFactors = F)
tax$OTU.ID=paste("OTU",tax$OTU.ID,sep="_")
rownames(tax)= tax$OTU.ID; tax=tax[,-1]
#####Modified the data set eliminating the first lettets and all data that were nof calsified put the name as Unclassified
tax$kingdom=gsub("k__", "", tax$kingdom); tax$kingdom=gsub("__", "Unclassified", tax$kingdom)
tax$phylum=gsub("p__", "", tax$phylum); tax$phylum=gsub("__", "Unclassified", tax$phylum)
tax$class=gsub("c__", "", tax$class); tax$class=gsub("__", "Unclassified", tax$class)
tax$order=gsub("o__", "", tax$order) ; tax$order=gsub("__", "Unclassified", tax$order)
tax$family=gsub("f__", "", tax$family) ; tax$family=gsub("__", "Unclassified", tax$family)
tax$genus=gsub("g__", "", tax$genus);tax$genus=gsub("__", "Unclassified", tax$genus)
tax$OTU.ID=rownames(tax)

#Subset AMF
tax= subset(tax, tax$phylum == "Glomeromycota", select= c("kingdom","phylum","class","order","family","genus", "OTU.ID"))
otu=otu[rownames(tax),]

#RICHNESS
otu.r=otu
otu.r=otu.r[rownames(tax),]
otu.r=otu.r[,colSums(otu.r)>0]
colSums(otu.r)
otu=otu[,colnames(otu.r)]
otu.r$ric=rowSums(otu.r)
otu.r=otu.r %>% 
  mutate(ric = ifelse(ric >0, 1, ric))
otu.r=cbind(otu.r,tax)

#Diversity analysis
otu.t=gather(data=otu.r, key = "sampleID", value = "ric", colnames(otu))
otu.t=otu.t %>% 
  mutate(ric = ifelse(ric >0, 1, ric))

for (i in meta$SampleID){
  otu.t[otu.t$sample == i, c("Sample","Season","Location","Specie")] = meta[i,c("Sample","Season","Location","Specie")]
}
#Sample
orden = otu.t %>% group_by(genus,ric) %>% summarise(ric = sum(ric))  %>% summarise (ric = mean(ric)) 
low = orden[orden$ric < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Refamily taxa if needed
orden=orden[order(orden$ric,decreasing = T),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
otu.t$genus=factor(otu.t$genus, level = c(niveles[niveles!="Unclassified"],"Unclassified", "low abundant"))
otu.t$sample=factor(otu.t$sample, level = pc2)

pdf(file = "Richness BY EACHSAMPLE-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=otu.t, aes(y=ric, x=sample, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Sample ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()

#Location
orden = otu.t %>% group_by(genus, Location, ric) %>% summarise(OTU.ID = unique(OTU.ID))  %>% summarise(OTU.ID = unique(OTU.ID))
low = orden[orden$ric < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Reordeer taxa if needed
orden=orden[order(orden$genus, rev(orden$genus), decreasing = T),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
#otu.t$Location=factor(otu.t$Location, level = pc)
orden$Location=factor(orden$Location, level = pc)
orden=subset(orden, orden$Location != "Amatitan", select=c("genus", "Location", "ric","OTU.ID"))
orden=subset(orden, orden$Location != "Penjamo", select=c("genus", "Location", "ric","OTU.ID"))

mean.l=summarySE(data = orden, measurevar = "ric", groupvars = "Location")
kruskal.test(ric~Location, data = orden) ; dunnTest(ric~Location, data = orden, method = "bh")

pdf(file = "Richness BY LOCATION-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=orden, aes(y=ric, x=Location, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  #scale_fill_manual(values = c("darkblue","darkred"))+
  scale_fill_manual(values = c("#4A5DCB","#CE0937"))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Location ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()


##Summarise TOTAL OTUS BY LOCATION
TOTALOTUS= orden %>%
  group_by(Location) %>%
  summarise(ric = sum (ric))
PHOSPORUS = as.data.frame(c(87, 42.17, 18.11, 74,4.53))
colnames(PHOSPORUS)= "phosphorus"
TOTALOTUS=cbind(TOTALOTUS,PHOSPORUS)

TOTALOTUS=TOTALOTUS[order(TOTALOTUS$phosphorus),]

write.table(TOTALOTUS, file = paste(paste("Linear reggresion","txt",sep = "."),sep = "/"),
row.names = T, quote = F, sep = "\t", col.names = T)


###REGRESION LINEAL

rownames(TOTALOTUS)=TOTALOTUS$Location; TOTALOTUS=TOTALOTUS[,-1]
cor(TOTALOTUS)
regresion=lm(phosphorus~ric, data=TOTALOTUS)
summary(regresion)

TOTALOTUS$LOC=c("SF","Ma","BR","PF","AG")

pdf(file = "Regresion-richeness.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
r=ggplot(data=TOTALOTUS, aes(y=ric, x=phosphorus))+
  geom_point(size=6, aes(shape=LOC))+
  geom_smooth(method = glm)+
  scale_shape_manual(values = c(15:19))+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous(breaks = c(0,30,60,90))+
  annotate("text", x = 14, y = 5,
           label = "paste(italic(R) ^ 2, \" = 0.70\")", parse = TRUE)+tema
print(r)
dev.off()

summary(regresion)
shapiro.test(resid(regresion))
anova(regresion)


#Specie
orden = otu.t %>% group_by(genus, Specie, ric) %>% summarise(OTU.ID = unique(OTU.ID))  %>% summarise(OTU.ID = unique(OTU.ID))
low = orden[orden$ric < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Reordeer taxa if needed
orden=orden[order(orden$genus, rev(orden$genus), decreasing = F),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
#otu.t$Location=factor(otu.t$Location, level = pc)
orden$Specie=factor(orden$Specie, level = pc)

pdf(file = "Richness BY SPECIE-genus.pdf", colormodel = "cmyk", width = 8.5, height = 11, compress = F)
a=ggplot(data=orden, aes(y=ric, x=Specie, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Specie ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()

#Compartment
orden = otu.t %>% group_by(genus, Sample, ric) %>% summarise(OTU.ID = unique(OTU.ID))  %>% summarise(OTU.ID = unique(OTU.ID))
low = orden[orden$ric < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Reordeer taxa if needed
orden=orden[order(orden$genus, rev(orden$genus), decreasing = F),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
#otu.t$Location=factor(otu.t$Location, level = pc)
orden$Sample=factor(orden$Sample, level = pc)

pdf(file = "Richness BY SAMPLE-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=orden, aes(y=ric, x=Sample, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Compartment ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()


#Season
orden = otu.t %>% group_by(genus, Season, ric) %>% summarise(OTU.ID = unique(OTU.ID))  %>% summarise(OTU.ID = unique(OTU.ID))
low = orden[orden$ric < 0,"genus"]
otu.t[otu.t$genus %in% low$genus, "genus" ] = "low abundant"
n=dim(orden)[1]-length(low$genus)+1 ; print(c(n,"taxa")) 

#Reordeer taxa if needed
orden=orden[order(orden$genus, rev(orden$genus), decreasing = F),]
niveles=orden[!orden$genus %in% low$genus, "genus"]$genus
#otu.t$Location=factor(otu.t$Location, level = pc)
orden$Season=factor(orden$Season, level = pc)

pdf(file = "Richness BY SEASON-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=orden, aes(y=ric, x=Season, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Season ")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()


#####SUBSET FOR RICHNESS
#LOCATION
otu.t.SF=subset(otu.t, otu.t$Location == "SanFelipe")
otu.t.Ma=subset(otu.t, otu.t$Location == "Magueyal")
otu.t.SF$sampleID=factor(otu.t.SF$sampleID, level = pc1)
otu.t.Ma$sampleID=factor(otu.t.Ma$sampleID, level = pc1)

#Plot SANFELIPE
pdf(file = "Richness BY Location-SF-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
b=ggplot(data=otu.t.SF, aes(y=ric, x=sampleID, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Location in San Felipe")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(b)
dev.off()

#Plot SANFELIPE
pdf(file = "Richness BY Location-Ma-genus.pdf", colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a=ggplot(data=otu.t.Ma, aes(y=ric, x=sampleID, fill=genus))+
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = c("darkblue","darkred"))+
  #scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "richness")+
  ggtitle("Richness by Location in El Magueyal")+
  guides(fill=guide_legend(title = "genus", ncol = 1))+tema
print(a)
dev.off()
















