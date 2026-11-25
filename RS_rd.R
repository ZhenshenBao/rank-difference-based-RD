##############################pre-disease state detection###########
library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(KEGGgraph)
library(openxlsx)
library(philentropy)
library(GEOquery)
library(readxl)
library(fgsea)


rdps <- function(tt1, same, edgelsits, hours){
    rs <- NULL
    rexps1 <- NULL
    rexps2 <- NULL
    rownames(tt1) <- tt1[,2]
    tt1 <- tt1[hours,]
    samsone <- tt1[,3]
    
    bexp <-same[,1]
    bfc <- rank(bexp)
    
    rbexp1 <- bfc[as.vector(t(edgelsits[,1]))]
    rexps1 <- cbind(rexps1,rbexp1)
    rbexp2 <- bfc[as.vector(t(edgelsits[,2]))]
    rexps2 <- cbind(rexps2,rbexp2)
    rbexp <- abs(rbexp2-rbexp1)
    
    rs <- cbind(rs,rbexp)
    
    for(j in 2:16){
      jexp <- same[,j]
      names(jexp) <- rownames(allgvs)
      jexp <- jexp+0.001
      jfc <- rank(jexp)
      
      rjexp1 <- jfc[as.vector(t(edgelsits[,1]))]
      rexps1 <- cbind(rexps1,rjexp1)
      rjexp2 <- jfc[as.vector(t(edgelsits[,2]))]
      rexps2 <- cbind(rexps2,rjexp2)
      rjexp <- abs(rjexp1-rjexp2)
      rs <- cbind(rs,rjexp)

      rs[is.na(rs)] <- 0
      cat("\n");
      cat(paste("Done  ", j,sep=""))
    }
    ranks1 <- rs[,1]
    sds1 <- ranks1
    
    rankscore <- NULL
    for(ii in 2:16){
      ranks2 <- rs[,ii]

      ranks <- abs(ranks2-sds1)
      rankscore <- cbind(rankscore,ranks)
      cat("\n");
      cat(paste("Done  ", ii,sep=""))
    }
    rankscore <- -apply(-rankscore,2,sort)
    rankscore <- rankscore[1:1690,]
    rankscores <- apply(rankscore,2,sum)/length(ranks)
    
    
    alls <- rankscores
    names(alls) <- hours[2:16]
    alls1 <- abs(alls-alls[1])
    alls2 <- abs(alls-alls[2])
    alls3 <- abs(alls-alls[3])
  
    nalls <- NULL
    als <- cbind(alls1, alls2, alls3)
    nalls <- c(nalls, apply(als,1,min))
    return(nalls)
}

#####################################################
FINet <- FIsInGene_122220_with_annotations[,1:2]
edgelsits <- FINet
genesy <- AnnotationDbi::select(org.Hs.eg.db,keys = rownames(allgvs),columns = c("ENSEMBL","SYMBOL","ENTREZID"),keytype = "ENTREZID")
edgelsits = edgelsits[edgelsits[,1] %in% genesy$SYMBOL, ]
edgelsits = edgelsits[edgelsits[,2] %in% genesy$SYMBOL, ]
edgelsits <- cbind(as.character(genesy$ENTREZID[match(edgelsits[,1], genesy$SYMBOL)]),as.character(genesy$ENTREZID[match(edgelsits[,2], genesy$SYMBOL)])) 


GSE30550_series_matrix <- read.delim("GSE30550_series_matrix.txt", header=FALSE)
GPL9188 <- read_excel("GPL9188.xlsx")
anot <- GPL9188
anot <- as.data.frame(anot)
colnames(anot) <- c("ID","ENTREZID")
anot <- as.data.frame(anot)
samples <- GSE30550_series_matrix[1,]
exdata <- GSE30550_series_matrix[-1,]
colnames(exdata) <- exdata[1,]
exdata <- exdata[-1,]
rownames(exdata) <- exdata$ID_REF
exdata <- exdata[,-1]
esetm = exdata[rownames(exdata) %in% anot$ID, ]
esetm$ENTREZID <- as.character(anot$ENTREZID[match(rownames(esetm), anot$ID)])
esetm<-esetm[!is.na(esetm$ENTREZID),]
esetm<-esetm[!duplicated(esetm$ENTREZID),]
rownames(esetm) <- as.character(esetm$ENTREZID)
esetm <- esetm[,1:268]
samples <- samples[-1]
circsp<-strsplit(as.character(samples),", ")
samples<-do.call(rbind,circsp)
tts <- cbind(samples,colnames(esetm))
colnames(tts) <- c("subject","hour","GSM")

dsubjects <- c("Subject 01","Subject 05","Subject 06","Subject 07","Subject 08","Subject 10","Subject 12","Subject 13","Subject 15")#???? 8 13ȱʧ
nsubjects <- c("Subject 02","Subject 03","Subject 04","Subject 09","Subject 11","Subject 14","Subject 16","Subject 17")#???? 17ȱʧ
subjects <- c(dsubjects,nsubjects)
hours<-c("Baseline","Hour 00","Hour 005","Hour 012",
         "Hour 021","Hour 029","Hour 036","Hour 045",
         "Hour 053","Hour 060","Hour 069","Hour 077",
         "Hour 084","Hour 093","Hour 101","Hour 108")
allgvs <- esetm
allgvs <- apply(allgvs,2,as.numeric)
rownames(allgvs) <- as.character(rownames(esetm))


tts<-rbind(tts,c("Subject 08","Hour 021","GSM758012"))
tts<-rbind(tts,c("Subject 13","Baseline","GSM758088"))
tts<-rbind(tts,c("Subject 13","Hour 036","GSM758092"))
tts<-rbind(tts,c("Subject 17","Hour 036","GSM758155"))
tt1 <- tts[which(tts[,1] == "Subject 01"),]

same <- allgvs[,samsone]
same <- apply(same,2,as.numeric)
rownames(same) <- rownames(allgvs)

hours<-c("Baseline","Hour 00","Hour 005","Hour 012",
         "Hour 021","Hour 029","Hour 036","Hour 045",
         "Hour 053","Hour 060","Hour 069","Hour 077",
         "Hour 084","Hour 093","Hour 101","Hour 108")

res <- rdps(tt1, same, edgelsits, hours)
#############################################################


#############################signaling pathway enrichment analysis##############
erdp <- function(edgelsits1, de, aT1, esetm_e){
  p_e <- NULL
  pSize <- NULL
  NDE <- NULL
  x4es <- NULL
  for(j in 1:length(edgelsits1)){
    edgelsits1i <- edgelsits1[[j]]
    x1 <- length(intersect(edgelsits1i,de))
    x2 <- length(de)
    x3 <- nrow(esetm_e)-length(intersect(rownames(esetm_e),edgelsits1i))
    x4 <- length(intersect(rownames(esetm_e),edgelsits1i))
    x4e <- intersect(rownames(esetm_e),edgelsits1i)
    x4e <- paste(na.omit(x4e), collapse = ",")
    x4es <- c(x4es,x4e)
    p <- phyper(q=x1-1,m=x4, n=x3, k=x2,lower.tail=FALSE)
    p_e <- c(p_e,p)
    pSize <- c(pSize,x4)
    NDE <- c(NDE,x1)
  }
  peFdr=p.adjust(p_e,"fdr")
  res<-data.frame(Name=pathwaynames, ID=pathwayIDs, pSize=pSize, NDE=NDE, Ratio= NDE/pSize, pNDE=p_e, peFdr=peFdr, edges=x4es, stringsAsFactors=FALSE)
  res<-res[order(res$peFdr),]
  rownames(res)<-NULL
  return(res)
}

##############################
library(KEGGgraph)
library(igraph)
library(limma)
library(Biobase)
library(PADOG)
library(circlize)
library(GOplot)
library(grid)
library(graphics)
library(ComplexHeatmap)
library(ggplot2)
kgml_file <- "/Users/baozhenshen/R workspace/kgml_downloads"
files <- list.files(kgml_file)
edgelsits1 <- NULL
pathwaynames <- NULL
pathwayIDs <- NULL
for(i in 97:length(files)){
  mapkG <- parseKGML2Graph(paste(kgml_file,"/",files[[i]], sep=""),expandGenes=TRUE, genesOnly = TRUE)
  mapkEdges <- KEGGgraph::edges(mapkG)
  mapkEdges <- mapkEdges[sapply(mapkEdges, length) > 0]
  KEGGpathway <- parseKGML(paste(kgml_file,"/",files[[i]], sep=""))
  pathwaynames <- c(pathwaynames, KEGGpathway@pathwayInfo@title)
  pathwayIDs <- c(pathwayIDs, KEGGpathway@pathwayInfo@number)
  if(length(mapkEdges)>0){
    res <- lapply(1:length(mapkEdges), function(t){
      name <- names(mapkEdges)[t]
      len  <- length(mapkEdges[[t]])
      do.call(rbind, lapply(1:len, function(n){
        c(name, mapkEdges[[t]][n])
      }))
    })
    result <- data.frame(do.call(rbind, res))
    circsp<-strsplit(as.character(result[,1]),":")
    res1<-do.call(rbind,circsp)
    circsp<-strsplit(as.character(result[,2]),":")
    res2<-do.call(rbind,circsp)
    edgel <- paste(res1[,2],"-",res2[,2],sep="")
    edgel <- edgel[!duplicated(edgel)]
    edgelsits1[[i-96]] <- edgel
  }
} 
names(edgelsits1) <- pathwayIDs



set = "GSE4107"
data(list = set, package = "KEGGdzPathwaysGEO")
sets <- get(set)
#write a function to extract required info into a list
getdataaslist = function(x) {
  x = get(x)
  exp = experimentData(x)
  dataset = exp@name
  disease = notes(exp)$disease
  dat.m = exprs(x)
  ano = pData(x)
  design = notes(exp)$design
  annotation = paste(x@annotation, ".db", sep = "")
  targetGeneSets = notes(exp)$targetGeneSets
  list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
  names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
                  "targetGeneSets")
  return(list)
}
dlist = getdataaslist(set)

block = dlist$ano$Block
Block = block
block = factor(Block)
group = dlist$ano$Group
G = factor(group)
force(G)
force(block)
paired = dlist$design == "Paired"
esetm = dlist$dat.m
annotation = dlist$annotation

if (paired) {
  stopifnot(length(block) == length(group))
  stopifnot(all(table(block) == 2))
}

aT1 = filteranot(esetm, group, paired, block, annotation)
esetm = esetm[rownames(esetm) %in% aT1$ID, ]
rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]
rownames(esetm) <- as.character(rownames(esetm))
esetm_e <- NULL
for(i in 1:ncol(esetm)){
  iexp <- esetm[,i]
  bfc <- rank(iexp)
  rbexp1 <- bfc[as.vector(t(edgelsits[,1]))]
  rbexp2 <- bfc[as.vector(t(edgelsits[,2]))]
  rbexp <- abs(rbexp2-rbexp1)
  esetm_e <- cbind(esetm_e,rbexp)
}
colnames(esetm_e) <- colnames(esetm)
rownames(esetm_e) <- paste(edgelsits[,1],"-",edgelsits[,2],sep="")

gg1 <- dlist$ano

normalsamples <- dlist$ano$Sample[1:10]
normale <- esetm_e[,normalsamples]
noramles <- apply(normale,1,mean)
diseasesamples <- dlist$ano$Sample[11:22]
diseasee <- esetm_e[,diseasesamples]
diseasees <- apply(diseasee,1,mean)
edv <- abs(diseasees-noramles)
edv[which(edv=="Inf")] <- 0
edv <- sort(edv,decreasing = T)
de <- names(edv[1:(length(edv)*0.05)])   
res <- erdp(edgelsits1, de, aT1, esetm_e)
#################################################

######################BRCA subtype###########################
library(org.Hs.eg.db)
library(clusterProfiler)
library(limma)
library(igraph)
library(openxlsx)
library(philentropy)
library(GEOquery)
library(readxl)
library(fgsea)
library(KEGGgraph)
library(cowplot) 
library(Rtsne)
library(ggplot2)

PANCANExp <- `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp`
SurvivalTable <- Survival_SupplementalTable_S1_20171025_xena_sp
cancernames <- SurvivalTable$cancer.type.abbreviation
cancernames <- cancernames[!duplicated(cancernames)]
samnames <- SurvivalTable$sample
samsp <- strsplit(as.character(samnames),"-")
samname <- do.call(rbind,samsp)
samnames <- cbind(samname,samnames)
samnames <- as.matrix(samnames)

normalsams <- samnames[which(samnames[,4] == "11"),]
normalsams <- normalsams[,5]
normalsams <- gsub("-",".",normalsams)
normalsams <- intersect(normalsams,colnames(PANCANExp))
normalsams <- gsub("\\.","\\-",normalsams)
cancersams <- samnames[which(samnames[,4] != "11"),]
cancersams <- cancersams[,5]

allEXP <- PANCANExp[!duplicated(PANCANExp$sample),]
rownames(allEXP) <- allEXP$sample
allEXP <- allEXP[,-1]


cancername <- "BRCA"
cancersurvival <- SurvivalTable[which(SurvivalTable$cancer.type.abbreviation == cancername), ]
TCGASubtype <- read.delim("~/Desktop/个人文件/基于秩调控强度计算/TCGASubtype.20170308.tsv")

intersams <- intersect(cancersurvival$sample,cancersams)
intersams <- gsub("-",".",intersams)
intersams <- intersect(colnames(PANCANExp),intersams)

internsams <- intersect(cancersurvival$sample,normalsams)
internsams <- gsub("-",".",internsams)
internsams <- intersect(colnames(PANCANExp),internsams)

normalsample <- gsub("\\.","-",internsams)

nexp <- allEXP[,internsams]
cEXP <- allEXP[,intersams]
BRCAEXP <- cbind(nexp,cEXP)
colnames(BRCAEXP) <- gsub("\\.", "-", colnames(BRCAEXP))

#genesy <- bitr(rownames(BRCAEXP), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
#BRCAE <- BRCAEXP[genesy[,1],]
#rownames(BRCAE)<- as.character(genesy[,2])
#BRCAEXP <- BRCAE

rownames(TCGASubtype) <- TCGASubtype[,1]
BRCAsubtype <- TCGASubtype[colnames(BRCAEXP),]
BRCAEXP1 <- BRCAEXP[,which(BRCAsubtype$Subtype_mRNA=="Normal")]
BRCAEXP1 <- BRCAEXP1[,intersect(colnames(BRCAEXP1),normalsample)]
BRCAEXP2 <- BRCAEXP[,which(BRCAsubtype$Subtype_mRNA=="LumA")]
BRCAEXP3 <- BRCAEXP[,which(BRCAsubtype$Subtype_mRNA=="LumB")]
BRCAEXP4 <- BRCAEXP[,which(BRCAsubtype$Subtype_mRNA=="Her2")]
BRCAEXP5 <- BRCAEXP[,which(BRCAsubtype$Subtype_mRNA=="Basal")]

yzname1 <- sample(colnames(BRCAEXP1), 30)
yzname2 <- sample(colnames(BRCAEXP2), 30)
yzname3 <- sample(colnames(BRCAEXP3), 30)
yzname4 <- sample(colnames(BRCAEXP4), 30)
yzname5 <- sample(colnames(BRCAEXP5), 30)

yzexp1 <- BRCAEXP[,yzname1]
yzexp2 <- BRCAEXP[,yzname2]
yzexp3 <- BRCAEXP[,yzname3]
yzexp4 <- BRCAEXP[,yzname4]
yzexp5 <- BRCAEXP[,yzname5]

xlname1 <- setdiff(colnames(BRCAEXP1),yzname1)
xlname2 <- setdiff(colnames(BRCAEXP2),yzname2)
xlname3 <- setdiff(colnames(BRCAEXP3),yzname3)
xlname4 <- setdiff(colnames(BRCAEXP4),yzname4)
xlname5 <- setdiff(colnames(BRCAEXP5),yzname5)

xlexp1 <- BRCAEXP[,xlname1]
xlexp2 <- BRCAEXP[,xlname2]
xlexp3 <- BRCAEXP[,xlname3]
xlexp4 <- BRCAEXP[,xlname4]
xlexp5 <- BRCAEXP[,xlname5]

FINet <- FIsInGene_122220_with_annotations[,1:2]
edgelsits <- FINet
genesy <- AnnotationDbi::select(org.Hs.eg.db,keys = rownames(allgvs),columns = c("ENSEMBL","SYMBOL","ENTREZID"),keytype = "ENTREZID")
edgelsits = edgelsits[edgelsits[,1] %in% genesy$SYMBOL, ]
edgelsits = edgelsits[edgelsits[,2] %in% genesy$SYMBOL, ]
edgelsits <- cbind(as.character(genesy$ENTREZID[match(edgelsits[,1], genesy$SYMBOL)]),as.character(genesy$ENTREZID[match(edgelsits[,2], genesy$SYMBOL)])) 


edgerank <- function(esetm){
  esetm_e <- NULL
  for(i in 1:ncol(esetm)){
    iexp <- esetm[,i]
    names(iexp) <- rownames(esetm)
    bfc <- rank(iexp)
    rbexp1 <- bfc[as.vector(t(edgelsits[,1]))]
    rbexp2 <- bfc[as.vector(t(edgelsits[,2]))]
    rbexp <- abs(rbexp2-rbexp1)
    esetm_e <- cbind(esetm_e,rbexp)
  }
  colnames(esetm_e) <- colnames(esetm)
  rownames(esetm_e) <- paste(edgelsits[,1],"-",edgelsits[,2],sep="")
  return(esetm_e)
}

xler1 <- edgerank(xlexp1)
xler2 <- edgerank(xlexp2)
xler3 <- edgerank(xlexp3)
xler4 <- edgerank(xlexp4)
xler5 <- edgerank(xlexp5)

yzer1 <- edgerank(yzexp1)
yzer2 <- edgerank(yzexp2)
yzer3 <- edgerank(yzexp3)
yzer4 <- edgerank(yzexp4)
yzer5 <- edgerank(yzexp5)

yzer <- cbind(yzer1,yzer2,yzer3,yzer4,yzer5)

xler <- cbind(xler2,xler3,xler4,xler5)
xlmean <- apply(xler, 1, mean)

xlmean1 <- apply(xler1,1,mean)
xle <- abs(xlmean - xlmean1)
xler <-  sort(xle,decreasing = T)
e1 <- names(xler[1:1000])

esp <- strsplit(as.character(e1),"-")
esp <- do.call(rbind,esp)
write.csv(esp,"edges_all.csv")

edges <- c(esp[,1],esp[,2])
genesy <- bitr(edges, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
genesy <- genesy[,2]
genesy <- genesy[!duplicated(genesy)]

KEGG<-enrichKEGG(genesy,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
dotplot(KEGG)


xlmean2 <- apply(xler2,1,mean)
xlmean3 <- apply(xler3,1,mean)
xlmean4 <- apply(xler4,1,mean)
xlmean5 <- apply(xler5,1,mean)

xle2 <- abs(xlmean2 - xlmean1)
xle3 <- abs(xlmean3 - xlmean1)
xle4 <- abs(xlmean4 - xlmean1)
xle5 <- abs(xlmean5 - xlmean1)

xler2 <- sort(xle2,decreasing = T)
xler3 <- sort(xle3,decreasing = T)
xler4 <- sort(xle4,decreasing = T)
xler5 <- sort(xle5,decreasing = T)

e2 <- names(xler2[1:50])
e3 <- names(xler3[1:50])
e4 <- names(xler4[1:50])
e5 <- names(xler5[1:50])

dee <- strsplit(as.character(dee),"-")
dee <- do.call(rbind,dee)
write.csv(dee,"dee.csv")

esp2 <- strsplit(as.character(e2),"-")
esp2 <- do.call(rbind,esp2)
write.csv(esp2,"edges_LumA.csv")

esp3 <- strsplit(as.character(e3),"-")
esp3 <- do.call(rbind,esp3)
write.csv(esp3,"edges_LumB.csv")

esp4 <- strsplit(as.character(e4),"-")
esp4 <- do.call(rbind,esp4)
write.csv(esp4,"edges_Her2.csv")

esp5 <- strsplit(as.character(e5),"-")
esp5 <- do.call(rbind,esp5)
write.csv(esp5,"edges_Basal.csv")

edges2 <- c(esp2[,1],esp2[,2])
genesy2 <- bitr(edges2, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
genesy2 <- genesy2[,2]
genesy2 <- genesy2[!duplicated(genesy2)]
edges3 <- c(esp3[,1],esp3[,2])
genesy3 <- bitr(edges3, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
genesy3 <- genesy3[,2]
genesy3 <- genesy3[!duplicated(genesy3)]
edges4 <- c(esp4[,1],esp4[,2])
genesy4 <- bitr(edges4, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
genesy4 <- genesy4[,2]
genesy4 <- genesy4[!duplicated(genesy4)]
edges5 <- c(esp5[,1],esp5[,2])
genesy5 <- bitr(edges5, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
genesy5 <- genesy5[,2]
genesy5 <- genesy5[!duplicated(genesy5)]

genesys <- c(genesy2,genesy3,genesy4,genesy5)
genesys <- genesys[!duplicated(genesys)]

library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(nrichplot)

KEGG_database <- 'hsa'
KEGG<-enrichKEGG(genesys,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
dotplot(KEGG)
#enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = FALSE)#circluar为指定是否环化，基因过多时建议设置为FALSE
KEGG2 <- pairwise_termsim(KEGG)
emapplot(KEGG2,showCategory =20, color = "p.adjust", layout = "kk")

des <- c(e2,e3,e4,e5)
es <- table(des)
dee <- names(es[which(es>=1)])

# BRCAall <- cbind(xlexp2,xlexp3,xlexp4,xlexp5)
# xler <- edgerank(BRCAall)
# xlmean <- apply(xler,1,mean)
# xle <- abs(xlmean - xlmean1)
# xler <- sort(xle,decreasing = T)
# ee <- names(xler[1:50])

yzertsne <- t(yzer[dee,])
yzertsne <- apply(yzertsne,2,scale)
#yzertsne <- t(yzer[dee,])
tsne_out <- Rtsne(yzertsne, dims = 2, perplexity = 30, max_iter = 1000, verbose = TRUE)
Species <- c(rep("Normal",30),rep("LumA",30),rep("LumB",30),rep("Her2",30),rep("Basal",30))
tsne_df <- data.frame(tsne_out$Y, Subtype = Species)
ggplot(tsne_df, aes(x = X1, y = X2, color = Species)) +
  geom_point() + theme_minimal() + 
  stat_ellipse(data=tsne_df,
               geom="polygon",
               aes(x=X1,y=X2,
                   group=Species,
                   fill=Species),
               alpha=0.3,
               lty="dashed",
               color="black",
               key_glyph="blank")+
  theme_bw() + theme_cowplot()