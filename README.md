# R code for mouse CNS data Figure 3
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
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

 Idents(merged)<-"RNA_snn_res.1"
  new.cluster.ids <- c("MicFPR1low", 
                       "MacFPR1high", "MacFPR1low", "MicFPR1low", "MacFPR1low", "DCs", 
                       "CD4 T cell", "CD8 T cell", "DCs", "CD4 T cell", "MicFPR1high",
                       "MacFPR1high", "MacFPR1high", "CD4 T cell", "MicFPR1high", "MicFPR1high",
                       "NK cell", "MacFPR1low", "MicFPR1low", "MicFPR1high", "MicFPR1low",
                       "CD8 T cell", "CD3 T cell", "B cell", "Plasma cell", "CD3 T cell")
  names(new.cluster.ids) <- levels(merged)
  merged <- RenameIdents(merged, new.cluster.ids)
  merged$celltype<-merged@active.ident

# cell proportion and dot plot
cell.prop<-as.data.frame(prop.table(table(Idents(merged), merged$group)))

colname(cell.prop)<-c("cluster","group","proportion")

ggplot(cell.prop,aes(x=Var2,y=Freq,fill=Var1))+
  
  geom_bar(stat="identity",position="fill")+
  
  ggtitle("")+
  
  theme_bw()+
  
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))

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

  # cellchat
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(tidydr)
merged@meta.data$sample_type <- paste(merged@meta.data$group, merged@active.ident, sep = "_")
merged.list <- SplitObject(merged, split.by = "group")
merged.list <- lapply(X = merged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = merged.list, nfeatures = 2000)
merged.list <- lapply(X = merged.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})
options(future.globals.maxSize = 8000 * 1024^2)

immune.anchors <- FindIntegrationAnchors(object.list = merged.list, anchor.features = features, 
                                         reduction = "rpca", k.anchor = 20)
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- "RNA_snn_res.0.6"
DimPlot(immune.combined, reduction = "tsne" ,label=T, pt.size = 1.5 )
p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "group")
p2 <- DimPlot(immune.combined, reduction = "umap", group.by = "celltype",label = TRUE,repel = TRUE)
plot1<-FeaturePlot(immune.combined, reduction = "tsne", features = c("Ptprc", "Cd3e", "Cd4", "Cd8a", "Cd19", "Klrb1c", "Itgam",   "Tmem119", "Cx3cr1", "Adgre1", "Ly6g", "Itgax"),cols = c("gray", "red"))
plot1

plot1<-FeaturePlot(immune.combined, reduction = "tsne", features = c("Fpr1", "Cd3e", "Cd4", "Cd8a", "Cd19", "Klrb1c", "Itgam",   "Tmem119", "Cx3cr1", "Adgre1", "Ly6g", "Itgax"),cols = c("gray", "red"))
plot1
new.cluster.ids <- c( "MacFPR1high", "DCs",   "MicFPR1high", "MicFPR1high", "CD4 T cell",   
                      "MacFPR1high", "CD8 T cell",  "B cell",  "NK cell",    "MacFPR1low",  
                      "MacFPR1high", "MicFPR1low",  "MicFPR1high", "DCs",    "MacFPR1low",  
                      "CD4 T cell",   "MicFPR1low",   "CD8 T cell", "CD4 T cell",   "CD3 T cell",   
                      "Plasma cell",   "MicFPR1low")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)

immune.combined$celltype<-NULL
immune.combined$celltype<-immune.combined@active.ident
table(immune.combined$celltype)
Idents(immune.combined)
Idents(immune.combined)<-factor(Idents(immune.combined),levels=c("MicFPR1high","MicFPR1low", "MacFPR1high","MacFPR1low","CD3 T cell","CD4 T cell", "CD8 T cell","DCs","NK cell","B cell", "Plasma cell"))

immune.combined$celltype<-factor(immune.combined$celltype, levels=c("MicFPR1high","MicFPR1low", "MacFPR1high","MacFPR1low","CD3 T cell","CD4 T cell", "CD8 T cell","DCs","NK cell","B cell", "Plasma cell"))
Idents(immune.combined)<-"celltype"
table(immune.combined@active.ident)

FPR1KO.object <- subset(immune.combined,group=="FPR1KO") # same as following code to set WT.cellchat
FPR1KO.data.input <- GetAssayData(FPR1KO.object, assay = "RNA", slot = "data")
FPR1KO.meta = FPR1KO.object@meta.data[,c("celltype", "group")] 
FPR1KO.meta$cellType %<>% as.vector(.)
FPR1KO.cellchat <- createCellChat(object = FPR1KO.data.input)
FPR1KO.cellchat <- addMeta(FPR1KO.cellchat, meta = FPR1KO.meta)
FPR1KO.cellchat <- setIdent(FPR1KO.cellchat, ident.use = "celltype")
FPR1KO.cellchat@DB <- CellChatDB.mouse 
FPR1KO.cellchat <- subsetData(FPR1KO.cellchat) 
future::plan("multisession", workers = 10) 
FPR1KO.cellchat <- identifyOverExpressedGenes(FPR1KO.cellchat)
FPR1KO.cellchat <- identifyOverExpressedInteractions(FPR1KO.cellchat)
FPR1KO.cellchat <- projectData(FPR1KO.cellchat, PPI.mouse)  
FPR1KO.cellchat <- computeCommunProb(FPR1KO.cellchat, raw.use=F)
FPR1KO.cellchat <- filterCommunication(FPR1KO.cellchat, min.cells = 10)
FPR1KO.cellchat <- computeCommunProbPathway(FPR1KO.cellchat)
FPR1KO.cellchat <- aggregateNet(FPR1KO.cellchat)
FPR1KO.cellchat <- netAnalysis_computeCentrality(FPR1KO.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
group2.net <- subsetCommunication(FPR1KO.cellchat)
save (FPR1KO.cellchat, file="FPR1KO.cellchat.rds")
WT.cellchat <- readRDS("WT.cellchat.rds")
FPR1KO.cellchat <- readRDS("FPR1KO.cellchat.rds")

object.list <- list(FPR1KO = FPR1KO.cellchat, WT = WT.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


