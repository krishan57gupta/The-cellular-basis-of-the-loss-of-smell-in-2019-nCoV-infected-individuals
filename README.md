
## iCTC: The-molecular-basis-of-loss-of-smell-in-2019-nCoV-infected-individuals

Among the prominent clinical symptoms such as fatigue, shortness of breath, fever, and cough, 2019-nCoV infected individuals 
often experience hyposmia/anosmia (decrease or loss of sense of smell). Angiotensin I Converting Enzyme 2 (ACE2), a key host 
receptor has now been established as an important moiety for the entry of 2019-nCoV into the cell. A multitude of studies 
estimated the expression of ACE2 in multiple organs including heart, kidney, intestines, lungs, buccal cavity, etc. The 
ongoing medical examinations and the autopsy reports of the diseased individuals strongly corroborate these 
organ/tissue-level molecular insights. Olfactory mucosa harbors multiple functionally distinct cell types. 
Zeroing in on the cell lineages that underpin infection associated loss of olfaction may provide new leads for 
diagnostics/clinical management of 2019-nCoV infected individuals. Our pointed bioinformatic analysis of single 
expression profiles underscored selective expression of ACE2 in a subset of horizontal basal cells (HBCs) and microvillar 
cells (MVCs) of the olfactory mucosa.  Inspection of the ACE2 levels in the olfactory mucosa of 5 additional mammalian 
species revealed comparable expression patterns. In summary, our findings pinpoint the molecular rationale of loss of smell 
in 2019-nCoV infected patients.


``` r
library(Seurat)
library(ggplot2)
```

``` r
folder="~/corona_project/plots_new/"
batch_list=list("P2","P3")
batch_data_list=list("P2"=1,"P3"=1)
sce_3_1_1_5000_23_new=list("1","2")
for( i in 1:length(batch_list))
{
  print(batch_list[[i]])
  s_object=Read10X(paste("~/corona_project/Input_files/",batch_list[[i]],sep=""))
  s_object=CreateSeuratObject(counts =s_object, min.cells = 0, min.features = 400, project = "P23")
  s_object[["percent.mt"]] <- PercentageFeatureSet(s_object, pattern = "^MT-")
  s_object <- subset(s_object, subset = nFeature_RNA >100 & nFeature_RNA <8000 & percent.mt <10)
  s_object@meta.data[, "run"] <- batch_list[i]
  s_object=NormalizeData(s_object)
  batch_data_list[[i]]=FindVariableFeatures(s_object, selection.method = "dispersion", nfeatures =5000)
  run.combined <- ScaleData(batch_data_list[[i]], verbose = FALSE)
  run.combined <- RunPCA(run.combined, npcs = 30, verbose = FALSE)
  # t-SNE and Clustering
  run.combined <- RunUMAP(run.combined, reduction = "pca", dims = 1:30)
  run.combined <- FindNeighbors(run.combined, reduction = "pca", dims = 1:30)
  sce_3_1_1_5000_23_new[[i]] <- FindClusters(run.combined, resolution = 0.5)
}
p2<-RenameIdents(sce_3_1_1_5000_23_new[[1]],`0`="Fibroblasts", `1`="CD8T cells",`2`="Pericytes", `3`="Olfactory HBCs", `4`="Pericytes",
                 `5`="Vascular Smooth Muscle cells", `6`="Macrophages", `7`="Subtentacular Cells", `8`="Neuron Cells",
                 `9`="Plasma Cells", `10`="Bowman's Gland", `11`="Bowman's Gland", `12`="Respiratory HBC Cells",
                 `13`="Olfactory Ensheathing Glia", `14`="Dendritic cells", `15`="Monocytes", `16`="Subtentacular Cells", 
                 `17`="Natural Killer cells", `18`="Fibroblasts", `19`="Mast Cells", `20`="Fibroblasts", `21`="Olfactory Progenator Cells",
                 `22`="GBCs", `23`="Pericytes" )
p3<-RenameIdents(sce_3_1_1_5000_23_new[[2]],`0`="Pericytes", `1`="Fibroblasts",`2`="CD8T cells", `3`="Vascular Smooth Muscle cells",
                 `4`="Plasma Cells", `5`="Vascular Smooth Muscle cells", `6`="Pericytes", `7`="Olfactory HBCs",
                 `8`="CD4T Cells", `9`="B Cells", `10`="Bowman's Gland", `11`="Olfactory HBCs", `12`="Immature Neurons",
                 `13`="Olfactory Microvillar Cells", `14`="Natural Killer cells", `15`="Subtentacular Cells", `16`="Monocytes", 
                 `17`="Subtentacular Cells", `18`="Dendritic cells", `19`="Olfactory Ensheathing Glia", `20`="GBCs",
                 `21`="Mature Neurons", `22`="Plasma Cells", `23`="Macrophages" ,`24`="Mast Cells")




p2s=subset(p2, idents = names(table(Idents(p2)))[c(4,8,10,18,7)])
p3s=subset(p3, idents = names(table(Idents(p3)))[c(17,18,10,13,11,6,9)])

a=table(Seurat::Idents(p2s)[colnames(p2s[,as.matrix(p2s@assays$RNA["ACE2",])>0])])
a=a*100/sum(a)
df=data.frame(cell_type=names(a),count_percentage=a)
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
a=table(Seurat::Idents(p3s)[colnames(p3s[,as.matrix(p3s@assays$RNA["ACE2",])>0])])
a=a*100/sum(a)
df=data.frame(cell_type=names(a),count_percentage=a)

ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))



DimPlot(p3s)

DimPlot(p2s)


genes_marker=c("ACE2",
               "TP63","KRT5","CXCL14","SOX2","MEG3",
               "ASCL3","CFTR","HEPACAM2","FOXL1",
               "GNG8","OLIG2","EBF2","LHX2","CBX8",
               "GNG13","EBF2,","CBX8","RTP1",
               "HES6","ASCL1","CXCR4","SOX2","EZH2","NEUROD1","NEUROG1",
               "CYP2A13","CYP2J2","GPX6","ERMN","SOX2",
               "SOX9", "SOX10", "MUC5", "GPX3")
genes_index=c(1,2,7,11,16,20,27,33,37)
cells_marker=c("ACE2",names(table(Idents(p3s)))[c(1,4,3,7,6,5,2)])
for(i in 1:length(cells_marker))
{
  print(i)
  print(FeaturePlot(p2s, features = genes_marker[genes_index[i]:(genes_index[i+1]-1)],combine = FALSE,cols=c("lightgrey", "red")))
  print(FeaturePlot(p3s, features = genes_marker[genes_index[i]:(genes_index[i+1]-1)],combine = FALSE,cols=c("lightgrey", "red")))
}




p2sHBC=subset(p2, idents = names(table(Idents(p2)))[c(4)])
pos=c(rep(1,length(colnames(p2sHBC[,as.matrix(p2sHBC@assays$RNA["ACE2",])>0]))))
names(pos)=colnames(p2sHBC[,as.matrix(p2sHBC@assays$RNA["ACE2",])>0])
neg=c(rep(-1,length(colnames(p2sHBC[,as.matrix(p2sHBC@assays$RNA["ACE2",])<=0]))))
names(neg)=colnames(p2sHBC[,as.matrix(p2sHBC@assays$RNA["ACE2",])<=0])
pos_neg=c(pos,neg)
sum(pos_neg[which(names(pos_neg)%in%rownames(p2sHBC@meta.data))]==pos_neg)
p2sHBC=AddMetaData(
  object = p2sHBC,
  metadata = pos_neg,
  col.name = 'HBC'
)
Idents(p2sHBC)=pos_neg
sum(pos_neg[which(names(pos_neg)%in%names(Idents(p2sHBC)))]==pos_neg)
p2_HBC_genes=FindMarkers(p2sHBC, ident.1=1,ident.2=-1,test.use = "wilcox")
p2_HBC_genes_a=rownames(p2_HBC_genes)[p2_HBC_genes[2]>1 & p2_HBC_genes[1]<0.05]
p2_HBC_genes_b=rownames(p2_HBC_genes)[p2_HBC_genes[2]<-1 & p2_HBC_genes[1]<0.05]
p2_HBC_genes_names=c(p2_HBC_genes_a,p2_HBC_genes_b)

p2sHBC_1 <- RunPCA(p2sHBC, npcs = 30, verbose = FALSE)
p2sHBC_1 = RunUMAP(p2sHBC_1, reduction = "pca", dims = 1:30)
p2sHBC_1=FindNeighbors(p2sHBC_1, reduction = "pca", dims = 1:30)
p2sHBC_1=FindClusters(p2sHBC_1, resolution = 0.5)
DimPlot(p2sHBC_1,pt.size = 1)



p3sHBC=subset(p3, idents = names(table(Idents(p3)))[c(6)])
pos=c(rep(1,length(colnames(p3sHBC[,as.matrix(p3sHBC@assays$RNA["ACE2",])>0]))))
names(pos)=colnames(p3sHBC[,as.matrix(p3sHBC@assays$RNA["ACE2",])>0])
neg=c(rep(-1,length(colnames(p3sHBC[,as.matrix(p3sHBC@assays$RNA["ACE2",])<=0]))))
names(neg)=colnames(p3sHBC[,as.matrix(p3sHBC@assays$RNA["ACE2",])<=0])
pos_neg=c(pos,neg)
sum(pos_neg[which(names(pos_neg)%in%rownames(p3sHBC@meta.data))]==pos_neg)
p3sHBC=AddMetaData(
  object = p3sHBC,
  metadata = pos_neg,
  col.name = 'HBC'
)
Idents(p3sHBC)=pos_neg
sum(pos_neg[which(names(pos_neg)%in%names(Idents(p3sHBC)))]==pos_neg)
p3_HBC_genes=FindMarkers(p3sHBC, ident.1=1,ident.2=-1,test.use = "wilcox")
p3_HBC_genes_a=rownames(p3_HBC_genes)[p3_HBC_genes[2]>1 & p3_HBC_genes[1]<0.05]
p3_HBC_genes_b=rownames(p3_HBC_genes)[p3_HBC_genes[2]<-1 & p3_HBC_genes[1]<0.05]
p3_HBC_genes_names=c(p3_HBC_genes_a,p3_HBC_genes_b)

p3sHBC_1 <- RunPCA(p3sHBC, npcs = 30, verbose = FALSE)
p3sHBC_1 = RunUMAP(p3sHBC_1, reduction = "pca", dims = 1:30)
p3sHBC_1=FindNeighbors(p3sHBC, reduction = "pca", dims = 1:30)
p3sHBC_1=FindClusters(p3sHBC_1, resolution = 0.5)
DimPlot(p3sHBC_1,pt.size = 1)




p3sMicro=subset(p3, idents = names(table(Idents(p3)))[c(11)])
pos=c(rep(1,length(colnames(p3sMicro[,as.matrix(p3sMicro@assays$RNA["ACE2",])>0]))))
names(pos)=colnames(p3sMicro[,as.matrix(p3sMicro@assays$RNA["ACE2",])>0])
neg=c(rep(-1,length(colnames(p3sMicro[,as.matrix(p3sMicro@assays$RNA["ACE2",])<=0]))))
names(neg)=colnames(p3sMicro[,as.matrix(p3sMicro@assays$RNA["ACE2",])<=0])
pos_neg=c(pos,neg)
sum(pos_neg[which(names(pos_neg)%in%rownames(p3sMicro@meta.data))]==pos_neg)
p3sMicro=AddMetaData(
  object = p3sMicro,
  metadata = pos_neg,
  col.name = 'Micro'
)
Idents(p3sMicro)=pos_neg
sum(pos_neg[which(names(pos_neg)%in%names(Idents(p3sMicro)))]==pos_neg)
p3_Micro_genes=FindMarkers(p3sMicro, ident.1=1,ident.2=-1,test.use = "wilcox")
p3_Micro_genes_a=rownames(p3_Micro_genes)[p3_Micro_genes[2]>1 & p3_Micro_genes[1]<0.05]
p3_Micro_genes_b=rownames(p3_Micro_genes)[p3_Micro_genes[2]<-1 & p3_Micro_genes[1]<0.05]
p3_Micro_genes_names=c(p3_Micro_genes_a,p3_Micro_genes_b)

p3sMicro_1 <- RunPCA(p3sMicro, npcs = 30, verbose = FALSE)
p3sMicro_1 = RunUMAP(p3sMicro_1, reduction = "pca", dims = 1:30)
p3sMicro_1=FindNeighbors(p3sMicro, reduction = "pca", dims = 1:30)
p3sMicro_1=FindClusters(p3sMicro_1, resolution = 0.5)
DimPlot(p3sMicro_1,pt.size = 1)
```

