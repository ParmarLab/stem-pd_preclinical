# Libraries
library(Seurat)
library(ggplot2)
library(ArchR)
col.ls <- ArchRPalettes[16:30]

# Read data
graft.seq.processed <- readRDS("raw_data/graft.seq.processed.20230602.rds")

# Calculate cycling cells

graft.seq.processed <- CellCycleScoring(
  object = graft.seq.processed,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)



# A - Predicted cycling boxplot
m1 <- reshape2::melt(prop.table(table(graft.seq.processed$orig.ident,graft.seq.processed$S_G2M.score>.4),1)[,-1])
m1$age <- factor(gsub("m","",stringr::str_match(rownames(m1),pattern = "[0-9]+m")),levels =c("3","6","9","12"))
m1<-data.frame(table(graft.seq.processed$orig.ident),m1)

p1<-ggplot(m1[m1$Freq>500,],aes(y=value*100,x=age,fill=value))+
  geom_boxplot(fill=col.ls$greyMagma[1])+ylab("% cells predicted to\nbe cycling")+xlab("Age [months]")+theme_cowplot()
p1
ggsave("plots/predicted_cycling.pdf",w=8,h=5)



# B -  KI67+ boxplot
m.ki<-data.frame(prop.table(table(GetAssayData(graft.seq.processed)["MKI67",]>0,graft.seq.processed$orig.ident),2)[2,]*100)
colnames(m.ki)<-"value"
m.ki$age <- factor(gsub("m","",stringr::str_match(rownames(m.ki),pattern = "[0-9]+m")),levels =c("3","6","9","12"))
m.ki<-data.frame(table(graft.seq.processed$orig.ident),m.ki)

p2<-ggplot(m.ki[m.ki$Freq>0,],aes(y=value,x=age,fill=value))+
  geom_boxplot(fill=col.ls$greyMagma[1])+ylab("% KI67+ cells")+xlab("Age [months]")+theme_cowplot()
p2
ggsave("plots/ki67_boxplot.with_dots.pdf",w=8,h=5)



