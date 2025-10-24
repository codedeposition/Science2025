# R code for mouse CNS data Figure 3
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(cowplot) 
     
# data input
WT1 <- Read10X(data.dir = "_/WT1/")
colnames(WT1) <- paste('WT1', colnames(WT1), sep = '_')
WT2 <- Read10X(data.dir = "_/WT2/")
colnames(WT2) <- paste('WT2', colnames(WT2), sep = '_')
WT3 <- Read10X(data.dir = "_/WT3/")
colnames(WT3) <- paste('WT3', colnames(WT3), sep = '_')
WT3<-WT3[, 1:4400]
FPR1KO1 <- Read10X(data.dir = "_/FPR1KO1/")
colnames(FPR1KO1) <- paste('FPR1KO1', colnames(FPR1KO1), sep = '_')
FPR1KO1<-FPR1KO1[,1:2260000]
FPR1KO2 <- Read10X(data.dir = "_/FPR1KO2/")
colnames(FPR1KO2) <- paste('FPR1KO2', colnames(FPR1KO2), sep = '_')
FPR1KO2<-FPR1KO2[,1:3700]
FPR1KO3 <- Read10X(data.dir = "_/FPR1KO3/")
colnames(FPR1KO3) <- paste('FPR1KO3', colnames(FPR1KO3), sep = '_')
FPR1KO3<-FPR1KO3[,1:3000]

scRNAlist0 <- list(WT1,WT2,WT3,FPR1KO1,FPR1KO2,FPR1KO3)
scRNAlist <- list()
# sreamline data processing
for(i in 1:length(scRNAlist0)){ sc <- scRNAlist0[[i]]
sc <- CreateSeuratObject(sc,  min.cells = 3, min.features = 300)
scRNAlist[[i]] <- sc
rm(sc)
}

for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, pattern = "^Hb[^(p)]")
  scRNAlist[[i]] <- sc
  rm(sc)
}

violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                layer = "counts",
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}

scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x,
              subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})
scRNAlist_merge <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])

group=str_split(colnames(scRNAlist_merge@assays$RNA),'_',simplify = T)[,1]
table(group)
str_detect(group,"^F")
table(str_detect(group,"^F"))
group1=ifelse(str_detect(group,"^F"),"FPR1KO","WT")

#PCA
scRNAlist_merge <- NormalizeData(scRNAlist_merge)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge)
scRNAlist_merge <- ScaleData(scRNAlist_merge, vars.to.regress = c("mt_percent"))
scRNAlist_merge <- RunPCA(scRNAlist_merge, verbose=F)

#Integrated with harmony
scRNA_harmony <- IntegrateLayers(object = scRNAlist_merge, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca", 
                                 new.reduction = "harmony",
                                 verbose = FALSE)
scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])

ElbowPlot(scRNA_harmony, ndims = 50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:30, reduction = "harmony")
scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:30, reduction = "harmony")

#DoubletFinder
library('DoubletFinder')
ddb_merged<-scRNA_harmony
sweep.res.list <- DoubletFinder::paramSweep(ddb_merged, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
mpK
dim(ddb_merged) 
nExp_poi <- round(0.075*ncol(ddb_merged))  #10000 cells:doublets rate is ～7.5%
nExp_poi
homotypic.prop <- modelHomotypic(ddb_merged$seurat_clusters) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
ddb_merged <- doubletFinder(ddb_merged, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
paste("DF.classifications_", "0.25_", mpK, "_", nExp_poi.adj, sep="")
ddb_merged$DF.classifications_0.25_0.03_1522
DimPlot(ddb_merged, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap", group.by = "DF.classifications_0.25_0.03_1522",shuffle = T)
DimPlot(ddb_merged, pt.size = 1, label = TRUE, label.size = 5, reduction = "tsne", group.by = "DF.classifications_0.25_0.03_1522",shuffle = T)

Epi_all_filter <- subset(ddb_merged, DF.classifications_0.25_0.03_1522 == "Singlet" )
Epi_all_filter
dim(Epi_all_filter)
Idents(Epi_all_filter) <- "seurat_clusters"
DimPlot(Epi_all_filter, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap",shuffle = T)
scRNA_harmony<-Epi_all_filter

# remove cell circule genes
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
s.genes<-mmus_s
g2m.genes<-mmus_g2m
scRNA_harmony <- CellCycleScoring(scRNA_harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

scRNA_harmony <- ScaleData(scRNA_harmony, vars.to.regress = c("S.Score", "G2M.Score"), 
                           features = rownames(scRNA_harmony))

#maker genes and cell anotation
merged<-scRNA_harmony
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

plot1<-FeaturePlot(merged, reduction = "umap", features = c("Ptprc", "Cd3e", "Cd4", "Cd8a", "Klrb1c", "Cd19", "Mzb1", "Itgam",   "Tmem119", "Adgre1", "Ly6g","S100a9"),cols = c("gray", "red"))
plot1
Idents(merged)<-"RNA_snn_res.4"
new.cluster.ids <- c("Mac2", "Mac2","Mac2","Mac1","Mac1","Mac2","Mac2","Mac1",
                     "Mac1","Mac2","Mac2","Mac2","Mac2","Mac1","Mac1","Mac2","Mac2","Mac2","Mac1","Mac1",
                     "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", "Mic2", 
                     "Mic1", "Mic1", "Mic1", "Mic1", "Mic1", "Mic1", "Mic1", "Mic1", 
                     "CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell",
                     "CD8 T cell", "CD8 T cell","CD8 T cell",
                     "DCs","DCs","DCs","DCs",
                     "Neutrophil","Neutrophil", "B cell", "Plasma cell", "NK cell", "NK cell")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
merged$celltype3<-merged@active.ident
merged$celltype3<-factor(merged$celltype3,levels=c("Mic1","Mic2","Mac1","Mac2","DCs","Neutrophil","CD4 T cell","CD8 T cell","B cell","Plasma cell", "NK cell"))

Idents(merged)<-merged$celltype3

merged@meta.data$sample_type <- paste(merged@meta.data$group, merged@active.ident, sep = "_")


#DimPlot for total samples
DimPlot(merged,reduction = "umap")


# cell proportion and dot plot
cell.prop<-as.data.frame(prop.table(table(Idents(merged), merged$group)))

colname(cell.prop)<-c("cluster","group","proportion")

ggplot(cell.prop,aes(x=Var2,y=Freq,fill=Var1))+
  
  geom_bar(stat="identity",position="fill")+
  
  ggtitle("")+
  
  theme_bw()+
  
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))

#VlnPlot for FPR1 expression
plot1<-VlnPlot(merged, features = c("Fpr1"),cols=c('#B6C6E0', '#E7ACA4'), split.by = "group",pt.size =0.1)
plot1 

# DEGs analysis 
highCells=colnames(subset(x = merged, subset = Fpr1 > 0, slot = 'counts'))
highORlow=ifelse(colnames(merged) %in% highCells,'Positive','Negative')
merged@meta.data$highORlow=highORlow

MicrogliaFPR1highmarker1 <- FindMarkers(merged, ident.1 = "Positive", ident.2 ="Negative",
                                        group.by = 'highORlow', subset.ident = "MicFPR1high",
                                        only.pos = F,min.pct = 0.15,logfc.threshold = 0.25)%>%mutate( gene = rownames(.) ) 

Michighvslow<-FindMarkers(merged, ident.1 = "MicFPR1high", ident.2 = "MicFPR1low",only.pos = FALSE,min.pct = 0.15,logfc.threshold = 0.25)%>% # ident 1 vs 2
  mutate( gene = rownames(.) )

merged@meta.data$sample_type <- paste(merged@meta.data$group, merged@active.ident, sep = "_")
marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "WT_MicFPR1high", ident.2 = "FPR1KO_MicFPR1high",only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.1)%>%mutate( gene = rownames(.) ) 

macrophageFPR1highmarker1 <- FindMarkers(merged, ident.1 = "Positive", ident.2 ="Negative",
                                         group.by = 'highORlow', subset.ident = "MacFPR1high",
                                         only.pos = F,min.pct = 0.15,logfc.threshold = 0.25)%>%mutate( gene = rownames(.) ) 

merged@meta.data$sample_type <- paste(merged@meta.data$group, merged@active.ident, sep = "_")
marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "WT_MacFPR1high", ident.2 = "FPR1KO_MacFPR1high",only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.15)%>%mutate( gene = rownames(.) ) 

# KEGG Go
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(gplots)
library(ggplot2)
library(mosaicData)
library(tidyverse)
library(conflicted)
library(dplyr)

Tcells_degs<-object # DEGs between groups e.g. MicrogliaFPR1highmarker1


Tcells_degs_fil = Tcells_degs %>% 
  dplyr::filter( pct.1 > 0.1 & p_val_adj < 0.05&p_val_adj > 0) %>%
  dplyr::filter( abs( avg_log2FC ) > 0.6 ) 
Tcells_degs_fil$gene <- rownames(Tcells_degs_fil)
ids=bitr(Tcells_degs_fil$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
Tcells_degs_fil=merge(Tcells_degs_fil,ids,by.x='gene',by.y='SYMBOL')
head(Tcells_degs_fil)

Tcells_degs_fil <- Tcells_degs_fil[order(Tcells_degs_fil$avg_log2FC,decreasing = T),]
Tcells_degs_list <- as.numeric(Tcells_degs_fil$avg_log2FC)
names(Tcells_degs_list) <- Tcells_degs_fil$ENTREZID
head(Tcells_degs_list)
cluster3_de <- names(Tcells_degs_list)[abs(Tcells_degs_list) > 0.6]
head(cluster3_de)
cluster3_ego <- enrichGO(cluster3_de, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
head(cluster3_ego)
dotplot(cluster3_ego, showCategory=30, title="XXX")

cluster3_ekg <- enrichKEGG(gene= cluster3_de, organism = "mmu", pvalueCutoff =1)
head(cluster3_ekg)
group=str_split(cluster3_ekg@result$Description,'- Mus',simplify = T)[,1]
cluster3_ekg@result$Description
cluster3_ekg@result$Description<-group
dotplot(cluster3_ekg, showCategory=10,title="XXX")

# cellchat analysis
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(tidydr)


Idents(merged)<-merged$celltype3
merged@meta.data$sample_type <- paste(merged@meta.data$group, merged@active.ident, sep = "_")
table(merged@meta.data$group)
#WT FPR1KO 
#13078   7487 
DimPlot(merged,reduction = "umap",label=T)
immune.combined<-merged
rm(merged)

#set cellchat file for WT

WT.object <- subset(immune.combined,group=="WT")
WT.data.input <- GetAssayData(WT.object, assay = "RNA", slot = "data")
WT.meta <- WT.object@meta.data[,c("celltype3", "group")]
WT.meta$celltype3 %<>% as.vector(.)
WT.cellchat <- createCellChat(object = WT.data.input)
WT.cellchat <- addMeta(WT.cellchat, meta = WT.meta)
WT.cellchat <- setIdent(WT.cellchat, ident.use = "celltype3")
levels(WT.cellchat@idents)
groupSize <- as.numeric(table(WT.cellchat@idents))
groupSize
WT.cellchat@DB <- CellChatDB.mouse 
showDatabaseCategory(WT.cellchat@DB)
dplyr::glimpse(CellChatDB.mouse$interaction) 
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling") 
WT.cellchat@DB <- CellChatDB.use 
WT.cellchat <- subsetData(WT.cellchat, features = NULL)
future::plan("multisession", workers = 10)
WT.cellchat <- identifyOverExpressedGenes(WT.cellchat)
WT.cellchat <- identifyOverExpressedInteractions(WT.cellchat)
WT.cellchat <- projectData(WT.cellchat, PPI.mouse)
WT.cellchat <- computeCommunProb(WT.cellchat,raw.use=F,population.size =F)
WT.cellchat <- filterCommunication(WT.cellchat, min.cells = 10)
WT.cellchat <- computeCommunProbPathway(WT.cellchat)
WT.cellchat <- aggregateNet(WT.cellchat)
WT.cellchat <- netAnalysis_computeCentrality(WT.cellchat, slot.name = "netP") 

#set cellchat file for FPR1KO

FPR1KO.object <- subset(immune.combined,group=="FPR1KO")
FPR1KO.data.input <- GetAssayData(FPR1KO.object, assay = "RNA", slot = "data")
FPR1KO.meta = FPR1KO.object@meta.data[,c("celltype3", "group")] 
FPR1KO.meta$cellType3 %<>% as.vector(.)
FPR1KO.cellchat <- createCellChat(object = FPR1KO.data.input)
FPR1KO.cellchat <- addMeta(FPR1KO.cellchat, meta = FPR1KO.meta)
FPR1KO.cellchat <- setIdent(FPR1KO.cellchat, ident.use = "celltype3")
groupSizeFPR1KO <- as.numeric(table(FPR1KO.cellchat@idents)) 
groupSizeFPR1KO
FPR1KO.cellchat@DB <- CellChatDB.mouse 
FPR1KO.cellchat <- subsetData(FPR1KO.cellchat) 
future::plan("multisession", workers = 10) 
FPR1KO.cellchat <- identifyOverExpressedGenes(FPR1KO.cellchat)
FPR1KO.cellchat <- identifyOverExpressedInteractions(FPR1KO.cellchat)
FPR1KO.cellchat <- projectData(FPR1KO.cellchat, PPI.mouse)  
FPR1KO.cellchat <- computeCommunProb(FPR1KO.cellchat, raw.use=F, population.size = F)
FPR1KO.cellchat <- filterCommunication(FPR1KO.cellchat, min.cells = 10)
FPR1KO.cellchat <- computeCommunProbPathway(FPR1KO.cellchat)
FPR1KO.cellchat <- aggregateNet(FPR1KO.cellchat)
FPR1KO.cellchat <- netAnalysis_computeCentrality(FPR1KO.cellchat, slot.name = "netP")
object.list <- list(FPR1KO = FPR1KO.cellchat, WT = WT.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat, measure = "weight")
h1
dev.off()
