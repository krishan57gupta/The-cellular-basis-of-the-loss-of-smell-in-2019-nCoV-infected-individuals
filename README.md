
## The-molecular-basis-of-loss-of-smell-in-2019-nCoV-infected-individuals

Among the prominent clinical symptoms such as fatigue, shortness of
breath, fever, and cough, 2019-nCoV infected individuals often
experience hyposmia/anosmia (decrease or loss of sense of smell).
Angiotensin I Converting Enzyme 2 (ACE2), a key host receptor has now
been established as an important moiety for the entry of 2019-nCoV into
the cell. A multitude of studies estimated the expression of ACE2 in
multiple organs including heart, kidney, intestines, lungs, buccal
cavity, etc. The ongoing medical examinations and the autopsy reports of
the diseased individuals strongly corroborate these organ/tissue-level
molecular insights. Olfactory mucosa harbors multiple functionally
distinct cell types. Zeroing in on the cell lineages that underpin
infection associated loss of olfaction may provide new leads for
diagnostics/clinical management of 2019-nCoV infected individuals. Our
pointed bioinformatic analysis of single expression profiles underscored
selective expression of ACE2 in a subset of horizontal basal cells
(HBCs) and microvillar cells (MVCs) of the olfactory mucosa. Inspection
of the ACE2 levels in the olfactory mucosa of 5 additional mammalian
species revealed comparable expression patterns. In summary, our
findings pinpoint the molecular rationale of loss of smell in 2019-nCoV
infected patients.

### Required data can can be downloaded by the link:

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522>

### Libraries need to be loaded before running

``` r
library(Seurat)
library(ggplot2)
library(ggpubr)
```

### initilization

``` r
folder="~/corona_project/plots_intigrate/"
cell_type_genes=list("Bowman’s glands"=c("SOX9", "SOX10", "MUC5", "GPX3"),
                     "olfactory HBCs"=c("TP63", "KRT5", "CXCL14", "SOX2", "MEG3"),
                     "olfactory ensheathing glia"=c("S100B", "PLP1", "PMP2","MPZ", "ALX3"),
                     "olfactory microvillar cells"=c("ASCL3", "CFTR", "HEPACAM2", "FOXL1"),
                     "immature neurons"=c("GNG8", "OLIG2", "EBF2", "LHX2", "CBX8"),
                     "mature neurons"=c("GNG13", "EBF2", "CBX8", "RTP1"),
                     "GBCs"=c("HES6", "ASCL1", "CXCR4", "SOX2", "EZH2", "NEUROD1", "NEUROG1"),
                     "sustentacular cells"=c("CYP2A13", "CYP2J2", "GPX6", "ERMN", "SOX2"))
marker_genes=c("SOX9", "SOX10", "MUC5", "GPX3",
               "TP63", "KRT5", "CXCL14", "SOX2", "MEG3",
               "S100B", "PLP1", "PMP2","MPZ", "ALX3",
               "ASCL3", "CFTR", "HEPACAM2", "FOXL1",
               "GNG8", "OLIG2", "EBF2", "LHX2", "CBX8",
               "GNG13", "EBF2", "CBX8", "RTP1",
               "HES6", "ASCL1", "CXCR4", "SOX2", "EZH2", "NEUROD1", "NEUROG1",
               "CYP2A13", "CYP2J2", "GPX6", "ERMN", "SOX2")
cell_type_names=c("Bowman’s glands",
                  "olfactory HBCs",
                  "olfactory ensheathing glia",
                  "olfactory microvillar cells",
                  "immature neurons",
                  "mature neurons",
                  "GBCs",
                  "sustentacular cells")
gene_path=c("LNPEP", "ACE", "CD143","PRCP","CTSG",
            "PREP", "KLK1_2","CMA1","REN", "MME",
            "THOP1","NLN","AGTR1","AGTR2","MAS1",
            " MRGPRDCPA3","ACE2","AGT","ANPEP","ENPEP","CTSA","ATP6AP2")
```

\#\#\#intigration and batch effect removing using Seurat

``` r
batch_list=list("P2","P3")
batch_data_list=list("P2"=1,"P3"=1)
for( i in 1:length(batch_list))
{
  print(batch_list[[i]])
  s_object=Read10X(paste("~/corona_project/Input_files/",batch_list[[i]],sep=""))
  s_object=CreateSeuratObject(counts =s_object, min.cells = 0, min.features = 400, project = "P23")
  s_object[["percent.mt"]] <- PercentageFeatureSet(s_object, pattern = "^MT-")
  s_object <- subset(s_object, subset = nFeature_RNA >100 & nFeature_RNA <8000 & percent.mt <10)
  s_object@meta.data[, "run"] <- batch_list[i]
  s_object=NormalizeData(s_object)
  batch_data_list[[i]]=FindVariableFeatures(s_object, selection.method = "vst", nfeatures =5000)
}
batch_data_list
saveRDS(batch_data_list,"~/corona_project/batch_data_list_5000_23.rds")
set.seed(1)
run.anchors <- FindIntegrationAnchors(object.list = batch_data_list, dims = 1:30,anchor.features = 5000)
a=run.anchors
run.combined <- IntegrateData(anchorset = run.anchors, dims = 1:30)
b=run.combined
sce_3_1_1_before=run.combined
saveRDS(sce_3_1_1_before,"~/corona_project/sce_3_1_1_before_5000_23.rds")
sce_3_1_1_before=readRDS("~/corona_project/sce_3_1_1_before_5000_23.rds")
```

### Dim reduction and clustering

``` r
DefaultAssay(run.combined) <- "integrated"
run.combined <- ScaleData(run.combined, verbose = FALSE)
run.combined <- RunPCA(run.combined, npcs = 30, verbose = FALSE)
run.combined <- RunUMAP(run.combined, reduction = "pca", dims = 1:30)
run.combined <- FindNeighbors(run.combined, reduction = "pca", dims = 1:30)
run.combined <- FindClusters(run.combined, resolution = 0.5)
```

### cell-types annotation

``` r
for(i in 1:length(marker_genes))
{
  if(i==3)
    next
  print(i)
  pdf(paste(folder,"others/vln_gene_markers_",marker_genes[i],".pdf",sep=""))
  print(VlnPlot(run.combined, features = marker_genes[i],combine = FALSE))
  dev.off()
  svg(paste(folder,"others/vln_gene_markers_",marker_genes[i],".svg",sep=""))
  print(VlnPlot(run.combined, features = marker_genes[i],combine = FALSE))
  dev.off()
  pdf(paste(folder,"others/feat_gene_markers_",marker_genes[i],".pdf",sep=""))
  print(FeaturePlot(run.combined, features = marker_genes[i],combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order = TRUE))
  dev.off()
  svg(paste(folder,"others/feat_gene_markers_",marker_genes[i],".svg",sep=""))
  print(FeaturePlot(run.combined, features = marker_genes[i],combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order = TRUE))
  dev.off()
}
run.combined<-RenameIdents(run.combined,`13`="Immature neurons", `14`="Mature neurons", `3`="Olfactory HBCs",`19`="Olfactory microvillar cells",
                   `5`="Bowman's gland", `16`="Olfactory ensheathing glia", `9`="Sustentacular cells",`20`="GBCs" )
run.combined=subset(run.combined, idents = c("Immature neurons","Mature neurons","Olfactory HBCs","Olfactory microvillar cells","Bowman's gland",
                          "Olfactory ensheathing glia","Sustentacular cells","GBCs"))
```

### Dim plot to show clusters clearly for each cell type

``` r
pdf(paste(folder,"tSNE_plots.pdf",sep=""))
DimPlot(run.combined, reduction = "umap")
dev.off()
svg(paste(folder,"tSNe_plots.svg",sep=""))
DimPlot(run.combined, reduction = "umap")
dev.off()
```

### FeaturePlot to show the expression of ACE2 in each cell type

``` r
pdf(paste(folder,paste("gene_markers_","ACE2",".pdf",sep="")))
FeaturePlot(run.combined, features = "ACE2",combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order =TRUE)
dev.off()
svg(paste(folder,paste("gene_markers_","ACE2",".svg",sep="")))
FeaturePlot(run.combined, features = "ACE2",combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order = TRUE)
dev.off()
```

### Percentage of ACE2 expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
print(sum(a))
a=a*100/sum(a)
df=data.frame(cell_type=names(a),count_percentage=a)
svg(paste(folder,"bar_ACE2.svg",sep=""))
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dev.off()
pdf(paste(folder,"bar_ACE2.pdf",sep=""))
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dev.off()
dim(run.combined@assays$RNA)

DefaultAssay(run.combined) <- "RNA"
```

### Finding differential expressed genes (HBC+;ACE2+) vs B (HBC+;ACE2-)

``` r
HBC=subset(run.combined, idents = "Olfactory HBCs")
pos=c(rep(1,length(colnames(HBC[,as.matrix(HBC@assays$RNA["ACE2",])>0]))))
names(pos)=colnames(HBC[,as.matrix(HBC@assays$RNA["ACE2",])>0])
neg=c(rep(-1,length(colnames(HBC[,as.matrix(HBC@assays$RNA["ACE2",])<=0]))))
names(neg)=colnames(HBC[,as.matrix(HBC@assays$RNA["ACE2",])<=0])
pos_neg=c(pos,neg)
sum(pos_neg[which(names(pos_neg)%in%rownames(HBC@meta.data))]==pos_neg)
HBC=AddMetaData(
  object = HBC,
  metadata = pos_neg,
  col.name = 'HBC'
)
Idents(HBC)=pos_neg
sum(pos_neg[which(names(pos_neg)%in%names(Idents(HBC)))]==pos_neg)
HBC_genes=FindMarkers(HBC, ident.1=1,ident.2=-1,test.use = "wilcox")
HBC_genes_a=rownames(HBC_genes)[HBC_genes[2]>1 & HBC_genes[1]<0.05]
HBC_genes_b=rownames(HBC_genes)[HBC_genes[2]<-1 & HBC_genes[1]<0.05]
HBC_genes_names=c(HBC_genes_a,HBC_genes_b)
```

### Finding differential expressed genes (Sustentacular+;ACE2+) vs B (Sustentacular+;ACE2-)

``` r
Sustentacular=subset(run.combined, idents = "Sustentacular cells")
pos=c(rep(1,length(colnames(Sustentacular[,as.matrix(Sustentacular@assays$RNA["ACE2",])>0]))))
names(pos)=colnames(Sustentacular[,as.matrix(Sustentacular@assays$RNA["ACE2",])>0])
neg=c(rep(-1,length(colnames(Sustentacular[,as.matrix(Sustentacular@assays$RNA["ACE2",])<=0]))))
names(neg)=colnames(Sustentacular[,as.matrix(Sustentacular@assays$RNA["ACE2",])<=0])
pos_neg=c(pos,neg)
sum(pos_neg[which(names(pos_neg)%in%rownames(Sustentacular@meta.data))]==pos_neg)
Sustentacular=AddMetaData(
  object = Sustentacular,
  metadata = pos_neg,
  col.name = 'Sustentacular'
)
Idents(Sustentacular)=pos_neg
sum(pos_neg[which(names(pos_neg)%in%names(Idents(Sustentacular)))]==pos_neg)
Sustentacular_genes=FindMarkers(Sustentacular, ident.1=1,ident.2=-1,test.use = "wilcox")
Sustentacular_genes_a=rownames(Sustentacular_genes)[Sustentacular_genes[2]>1 & Sustentacular_genes[1]<0.05]
Sustentacular_genes_b=rownames(Sustentacular_genes)[Sustentacular_genes[2]<-1 & Sustentacular_genes[1]<0.05]
Sustentacular_genes_names=c(Sustentacular_genes_a,Sustentacular_genes_b)
```

### Saving DE genes

``` r
write.csv(HBC_genes_a,paste(folder,"HBC_genes_a.csv",sep=""))
write.csv(HBC_genes_b,paste(folder,"HBC_genes_b.csv",sep=""))
write.csv(Sustentacular_genes_a,paste(folder,"Sustentacular_genes_a.csv",sep=""))
write.csv(Sustentacular_genes_b,paste(folder,"Sustentacular_genes_b.csv",sep=""))
```
