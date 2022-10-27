if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("clusterProfiler")
library(clusterProfiler)
setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//TargetGene//")
path="//Users//houxiangqing//Desktop//temp//Datasets_PPMI//TargetGene//"
data <- read.csv("search_result.csv")
data1 <- read.csv("search_result1.csv")
data2 <- read.csv("search_result2.csv")
data3 <- read.csv("search_result3.csv")
data4 <- read.csv("search_result4.csv")
data5 <- read.csv("search_result5.csv")
data6 <- read.csv("search_result6.csv")
data7 <- read.csv("search_result7.csv")
data8 <- read.csv("search_result8.csv")
data9 <- read.csv("search_result9.csv")

total<-rbind(data,data1,data2,data3,data4,data5,data6,data7,data8,data9)
total1<-total[c('ID','miRNA','Target')]
total1$Target[total1$Target=='']<-NA
total2<-na.omit(total1)
total3<-total2[!duplicated(total2$Target), ]
library(org.Hs.eg.db)
gene.df <- bitr(total3$Target,fromType="SYMBOL",toType=c
                ("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db)
gene.df1<-gene.df[!duplicated(gene.df$ENTREZID), ]

de<-gene.df1$ENTREZID
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="MF", readable=TRUE)
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
#ego3 <- mutate(ego, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))
library(ggplot2) 
library(forcats)
###################################Go analysis
#ggplot(ego3, showCategory = 10, aes(richFactor,
  fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) + 
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"), 
  trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) + 
  theme_dose(font.size = 12) +
  xlab("Rich Factor") +lab(NULL) + 
  #ggtitle("Biological Processes")

dotplot(ego2,title="Cellular Component")
plotGOgraph(ego2)
#########################KEGG analysis
library(clusterProfiler)
kk <- enrichKEGG(gene = gene.df1$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
dim(kk)
dotplot(kk,title="Enrichment KEGG_dot")

