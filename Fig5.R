# Figure5.G
library(ProjecTILs)
CD8T[["RNA3"]] <- as(object = CD8T[["RNA"]], Class = "Assay")
ref <- load.reference.map(ref = "CD8T_human_ref_v1.rds")
DimPlot(ref, label = T)
DefaultAssay(CD8T) <- "RNA3"

query.projected <- Run.ProjecTILs(CD8T, ref = ref, filter.cells=TRUE,split.by="Treatment")
naive<-subset(query.projected,subset=Treatment=="naive")
NACT<-subset(query.projected,subset=Treatment=="NACT") 
           
plot.projection(ref,NACT, linesize = 0.5, pointsize = 0.5)
plot.projection(ref,naive, linesize = 0.5, pointsize = 0.5)

# Figure5.H: 
library(ggplot2)
library(ggalluvial)
#naive subset
group_counts <- table(naive@meta.data$functional.cluster)
total_count <- sum(group_counts)
percentage_by_group <- (group_counts / total_count) * 100
percentage_naive <- data.frame(
   Population = names(percentage_by_group), 
   Percentage = percentage_by_group
  )

#NACT subset
group_counts <- table(NACT@meta.data$functional.cluster)
total_count <- sum(group_counts)
percentage_by_group <- (group_counts / total_count) * 100
percentage_NACT <- data.frame(
   Population = names(percentage_by_group), 
   Percentage = percentage_by_group
  )

# Combine naive+NACT
coldata<-data.frame("Treatment"=c(rep("naive",6),rep("NACT",7)))
per<-rbind(percentage_naive,percentage_NACT)
data<-data.frame(per,coldata$Treatment)
data$coldata.Treatment <- factor(data$coldata.Treatment, levels=c("naive", "NACT"))
data$Population <- factor(data$Population, levels=c("CD8.EM", "CD8.CM","CD8.MAIT","CD8.NaiveLike","CD8.TEMRA","CD8.TEX","CD8.TPEX"))

ggplot(data,
       aes(x = coldata.Treatment, stratum = Population, alluvium = Population,y=Percentage.Freq,
           fill = Population, label = Population)) +
  scale_fill_manual(values = c("darkseagreen","#F8766D","lightpink1","#d1cfcc","tan","#edbe2a","#A58AFF")) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme(legend.position = "bottom")+theme_classic()


# Figure5.I
genes <- as.data.frame(GetAssayData(object =CD8T, assay = "RNA", slot="data"))
b<-genes%>%filter(row.names(genes) %in% c( "TCF7", "CCR7", "IL7R", "LMNA", "FGFBP2",
    "XCL1", "CD200", "CRTAM", "TOX", "PDCD1", "HAVCR2", "GNLY",
    "KLRB1"))
b<-t(b)
b<-cbind(b,CD8T@meta.data$Treatment)
b<-data.frame(b)
colnames(b)[14]="Treatment"
#% of CD8 possitive cells. One txt for each gene
counts<-b%>%group_by(Treatment)%>%summarise(test=sum(TCF7>1))
total.counts<-table(b$Treatment)
tgd<-data.frame(total.counts,counts)
percentage<-round(tgd$test/tgd$Freq*100,1)
FOXP3<-data.frame(total.counts,counts,percentage)
print(FOXP3)
write.table(FOXP3,file="FOXP3.txt",sep="\t",col.names=NA)
names(b)

#combine files and put to one file
filelist <- list.files(pattern = ".*.txt")
library(readr)
a<-read_delim(filelist,id="name")
a$name=sub(".txt","", a$name)
View(a)
a<-data.table(a)

#Circle.plot

a$name <- factor(a$name, c("TCF7", "CCR7", "IL7R", "LMNA", "FGFBP2",
    "XCL1", "CD200", "CRTAM", "TOX", "PDCD1", "HAVCR2", "GNLY",
    "KLRB1"))
ggplot(a, aes(x = name, y = percentage, fill = Treatment)) +
  geom_col(position = "dodge")+scale_fill_manual(values=c("red","blue"))+
  coord_polar()+ylim(-10,60)+theme_bw()


#Supplemental Figure 7
genes4radar <- c("CD4", "CD8A", "TCF7", "CCR7", "IL7R", "LMNA", "GZMA", "GZMK", "FGFBP2",
    "XCL1", "CD200", "CRTAM", "TOX", "PDCD1", "HAVCR2", "PRF1", "GNLY", "GZMB", "TRAV1-2",
    "KLRB1", "FOXP3")
plot.states.radar(ref = ref, query = naive,query.assay='RNA3', genes4radar = genes4radar,
    min.cells = 20)



