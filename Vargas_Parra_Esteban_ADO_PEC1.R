##------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, 
                      comment = NA, prompt = TRUE, tidy = FALSE, 
                      fig.width = 7, fig.height = 7, fig_caption = TRUE,
                      cache=FALSE)
##------------------------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_TIME", "C")

  if(!(require(printr))) {
    install.packages(
      'printr',
      type = 'source',
      repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
    )
  }
##------------------------------------------------------------------------------------------------------------------------------------------
setwd(".")
dir.create("data")
dir.create("results")
#if(!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()
#install.packages("knitr")
#install.packages("colorspace")
#install.packages("gplots")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("htmlTable")
#install.packages("prettydoc")
#install.packages("devtools")
#install.packages("BiocManager")
#BiocManager::install("oligo")
#BiocManager::install("pd.clariom.s.human")
#BiocManager::install("arrayQualityMetrics")
#BiocManager::install("pvca")
# NOT NEEDED UNTIL ANALYSES ARE PERFORMED
#BiocManager::install("limma")
#BiocManager::install("genefilter")
#BiocManager::install("clariomshumantranscriptcluster.db")
#BiocManager::install("annotate")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("ReactomePA")
##------------------------------------------------------------------------------------------------------------------------------------------
targets<- read.csv2("./data/targets.csv", header = TRUE, sep = ";")
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Content of the targets file used for the current analysis')
##------------------------------------------------------------------------------------------------------------------------------------------
library(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)
library(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)
##------------------------------------------------------------------------------------------------------------------------------------------
my.targets@data$ShortName->rownames(pData(rawData))
colnames(rawData) <-rownames(pData(rawData))
head(rawData)
##------------------------------------------------------------------------------------------------------------------------------------------
library(arrayQualityMetrics)
#arrayQualityMetrics(rawData, outdir = file.path("./results", "raw.Data_quality"), force = TRUE)                         
##------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}
##------------------------------------------------------------------------------------------------------------------------------------------
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow", "darkviolet", "cyan"))
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/PCA_RawData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow", "darkviolet", "cyan"))
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
boxplot(rawData, cex.axis=0.5, las=2,  which="all", 
        col = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
        main="Distribution of raw intensity values")
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/Intensity_RawData.tiff", res = 200, width = 4, height = 4, units = 'in')
boxplot(rawData, cex.axis=0.5, las=2,  which="all", 
        col = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
        main="Distribution of raw intensity values")
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
eset_rma <- rma(rawData)
#arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"), force = TRUE)
##------------------------------------------------------------------------------------------------------------------------------------------
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow", "darkviolet", "cyan"))
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/PCA_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow", "darkviolet", "cyan"))
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
boxplot(eset_rma, cex.axis=0.5, las=2, which="all",
        col = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
        main="Boxplot for arrays intensity: Normalized Data")
##------------------------------------------------------------------------------------------------------------------------------------------  
tiff("figures/Intensity_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')  
boxplot(eset_rma, cex.axis=0.5, las=2, which="all",
        col = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
        main="Boxplot for arrays intensity: Normalized Data")
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
library(pvca)
pData(eset_rma)<-targets
pct_threshold<-0.6
batch.factors<-c("CellLine", "Response")
pvcaObj <- pvcaBatchAssess (eset_rma, batch.factors, pct_threshold)
##------------------------------------------------------------------------------------------------------------------------------------------
bp <- barplot(pvcaObj$dat, xlab = "Effects",
              ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("mediumorchid"), las=2,
              main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.55, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.7)
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/PVCAplot.tiff", res = 150, width = 5, height = 5, units = 'in')
bp <- barplot(pvcaObj$dat, xlab = "Effects",
              ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("mediumorchid"), las=2,
              main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.55, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.7)
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
sds <- apply (exprs(eset_rma), 1, sd)
sdsO<- sort(sds)
head(sdsO)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/SDplot.tiff", res = 150, width = 5, height = 5, units = 'in')
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
library(genefilter)
library(clariomshumantranscriptcluster.db)
annotation(eset_rma) <- "clariomshumantranscriptcluster.db"
filtered <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")
##------------------------------------------------------------------------------------------------------------------------------------------
names(filtered)
class(filtered$eset)
##------------------------------------------------------------------------------------------------------------------------------------------
print(filtered$filter.log)
eset_filtered <-filtered$eset
##------------------------------------------------------------------------------------------------------------------------------------------
write.csv(exprs(eset_rma), file="./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file="./results/normalized.Filtered.Data.csv")
save(eset_rma, eset_filtered, file="./results/normalized.Data.Rda")
##------------------------------------------------------------------------------------------------------------------------------------------
if (!exists("eset_filtered")) load (file="./results/normalized.Data.Rda")
##------------------------------------------------------------------------------------------------------------------------------------------
library(limma)
designMat<-model.matrix(~0+Group, pData(eset_filtered))
colnames(designMat) <- c("NCIH520.PR", "NCIH520.PS", "RH41.PR", "RH41.PS", "SJCRH30.PR", "SJCRH30.PS")
row.names(designMat) <- c("RH41.PS.1", "RH41.PS.2", "RH41.PS.3", 
                          "RH41.PR.1", "RH41.PR.2", "RH41.PR.3", "RH41.PR.4", "RH41.PR.5", 
                          "NCIH520.PS.1", "NCIH520.PS.2", "NCIH520.PS.3",
                          "NCIH520.PR.1", "NCIH520.PR.2", "NCIH520.PR.3",
                          "SJCRH30.PS.1", "SJCRH30.PS.2", "SJCRH30.PS.3",
                          "SJCRH30.PR.1", "SJCRH30.PR.2", "SJCRH30.PR.3")
print(designMat)
##------------------------------------------------------------------------------------------------------------------------------------------
cont.matrix <- makeContrasts (NCIH520.PRvsNCIH520.PS = NCIH520.PR-NCIH520.PS,
                              RH41.PRvsRH41.PS = RH41.PR-RH41.PS,
                              SJCRH30.PRvsSJCRH30.PS = SJCRH30.PR - SJCRH30.PS,
                              INT = (NCIH520.PR-NCIH520.PS) - (RH41.PR-RH41.PS) - (SJCRH30.PR - SJCRH30.PS),
                              levels=designMat)
print(cont.matrix)
##------------------------------------------------------------------------------------------------------------------------------------------
library(limma)
fit<- lmFit(eset_filtered, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
##------------------------------------------------------------------------------------------------------------------------------------------
topTab_NCIH520.PRvsNCIH520.PS <- topTable (fit.main, number=nrow(fit.main), coef="NCIH520.PRvsNCIH520.PS", adjust="fdr") 
head(topTab_NCIH520.PRvsNCIH520.PS)
##------------------------------------------------------------------------------------------------------------------------------------------
topTab_RH41.PRvsRH41.PS <- topTable (fit.main, number=nrow(fit.main), coef="RH41.PRvsRH41.PS", adjust="fdr") 
head(topTab_RH41.PRvsRH41.PS)
##------------------------------------------------------------------------------------------------------------------------------------------
topTab_SJCRH30.PRvsSJCRH30.PS <- topTable (fit.main, number=nrow(fit.main), coef="SJCRH30.PRvsSJCRH30.PS", adjust="fdr") 
head(topTab_SJCRH30.PRvsSJCRH30.PS)
##------------------------------------------------------------------------------------------------------------------------------------------
topTab_INT  <- topTable (fit.main, number=nrow(fit.main), coef="INT", adjust="fdr") 
head(topTab_INT)
##------------------------------------------------------------------------------------------------------------------------------------------
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}
##------------------------------------------------------------------------------------------------------------------------------------------
topAnnotated_NCIH520.PRvsNCIH520.PS <- annotatedTopTable(topTab_NCIH520.PRvsNCIH520.PS, anotPackage="clariomshumantranscriptcluster.db")
topAnnotated_RH41.PRvsRH41.PS <- annotatedTopTable(topTab_RH41.PRvsRH41.PS, anotPackage="clariomshumantranscriptcluster.db")
topAnnotated_SJCRH30.PRvsSJCRH30.PS <- annotatedTopTable(topTab_SJCRH30.PRvsSJCRH30.PS, anotPackage="clariomshumantranscriptcluster.db")
topAnnotated_INT <- annotatedTopTable(topTab_INT, anotPackage="clariomshumantranscriptcluster.db")
write.csv(topAnnotated_NCIH520.PRvsNCIH520.PS, file="./results/topAnnotated_NCIH520.PRvsNCIH520.PS.csv")
write.csv(topAnnotated_RH41.PRvsRH41.PS, file="./results/topAnnotated_RH41.PRvsRH41.PS.csv")
write.csv(topAnnotated_SJCRH30.PRvsSJCRH30.PS, file="./results/topAnnotated_SJCRH30.PRvsSJCRH30.PS.csv")
write.csv(topAnnotated_INT, file="./results/topAnnotated_INT.csv")
##------------------------------------------------------------------------------------------------------------------------------------------
short<- head(topAnnotated_NCIH520.PRvsNCIH520.PS[1:5,1:4])
library(kableExtra)
knitr::kable(
  short, booktabs = TRUE,
  caption = 'Annotations added to results "topTable" for the comparison "NCIH520.PRvsNCIH520.PS"')
show(short)
##------------------------------------------------------------------------------------------------------------------------------------------
library(clariomshumantranscriptcluster.db)
geneSymbols <- select(clariomshumantranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
par(mfrow=c(2,2))
for (i in colnames(cont.matrix)){
  volcanoplot(fit.main, coef=i, highlight=10, names=SYMBOLS,
              main=paste("Differentially expressed genes",i, sep="\n"))
  abline(v=c(-1,1))
}
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/VolcanoPlot.tiff", res = 150, width = 5, height = 5, units = 'in')
volcanoplot(fit.main, coef=1, highlight=4, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n")) 
abline(v=c(-1,1))
dev.off()

pdf("figures/Volcanos.pdf")
for (i in colnames(cont.matrix)){
  volcanoplot(fit.main, coef=i, highlight=4, names=SYMBOLS,
              main=paste("Differentially expressed genes",i, sep="\n"))
  abline(v=c(-1,1))
}
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
library(limma)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)
##------------------------------------------------------------------------------------------------------------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))
##------------------------------------------------------------------------------------------------------------------------------------------
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/VennPlot.tiff", res = 150, width = 5.5, height = 5.5, units = 'in')
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]

geneSymbols <- select(clariomshumantranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
write.csv(HMdata, file = file.path("./results/data4Heatmap.csv"))
##------------------------------------------------------------------------------------------------------------------------------------------
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
library(gplots)

heatmap.2(HMdata,
          Rowv = FALSE,
          Colv = FALSE,
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
          tracecol = NULL,
          dendrogram = "none",
          srtCol = 30)
##------------------------------------------------------------------------------------------------------------------------------------------
heatmap.2(HMdata,
          Rowv = TRUE,
          Colv = TRUE,
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
          tracecol = NULL,
          dendrogram = "both",
          srtCol = 30)
##------------------------------------------------------------------------------------------------------------------------------------------
tiff("figures/Heatmap1.tiff", res = 150, width = 5.5, height = 5.5, units = 'in')
heatmap.2(HMdata,
          Rowv = FALSE,
          Colv = FALSE,
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
          tracecol = NULL,
          dendrogram = "none",
          srtCol = 30)
dev.off()

tiff("figures/Heatmap2.tiff", res = 150, width = 5.5, height = 5.5, units = 'in')
heatmap.2(HMdata,
          Rowv = TRUE,
          Colv = TRUE,
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red", 3), rep("blue", 5), rep("green", 3), rep("yellow", 3), rep("darkviolet", 3), rep("cyan", 3)),
          tracecol = NULL,
          dendrogram = "both",
          srtCol = 30)
dev.off()
##------------------------------------------------------------------------------------------------------------------------------------------
listOfTables <- list(NCIH520.PRvsNCIH520.PS = topTab_NCIH520.PRvsNCIH520.PS,
                     RH41.PRvsRH41.PS = topTab_RH41.PRvsRH41.PS,
                     SJCRH30.PRvsSJCRH30.PS = topTab_SJCRH30.PRvsSJCRH30.PS,
                     INT = topTab_INT)
listOfSelected <- list()
for (i in 1:length(listOfTables)){
  topTab <- listOfTables[[i]]
  whichGenes<-topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  EntrezIDs<- select(clariomshumantranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)
##------------------------------------------------------------------------------------------------------------------------------------------
mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)
##------------------------------------------------------------------------------------------------------------------------------------------
library(ReactomePA)

listOfData <- listOfSelected[1:3]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "human",
                                 universe = universe)
  
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))
  
  if (length(rownames(enrich.result@result)) != 0) {
    write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
              row.names = FALSE)
    
    pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 7, 
                  title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
    dev.off()
    
    pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
                   vertex.label.cex = 0.75))
    dev.off()
  }
}
##------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
dev.off()
cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
         vertex.label.cex = 0.2)
##------------------------------------------------------------------------------------------------------------------------------------------
Tab.react <- read.csv2(file.path("./results/ReactomePA.Results.NCIH520.PRvsNCIH520.PS.csv"), 
                       sep = ",", header = TRUE, row.names = 1)

Tab.react <- Tab.react[1:4, 1:5]
knitr::kable(Tab.react, booktabs = TRUE, caption = "First rows and columns for Reactome results on KOvsWT.RT.csv comparison")
##------------------------------------------------------------------------------------------------------------------------------------------
listOfFiles <- dir("./results/") 
knitr::kable(
  listOfFiles, booktabs = TRUE,
  caption = 'List of files generated in the analysis',
  col.names="List_of_Files"
)