
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
infected
patients.

### Required data can can be downloaded by the link:

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522>

### To have a look on figures with code,you can download html file link given below.

<https://github.com/krishan57gupta/The-molecular-basis-of-loss-of-smell-in-2019-nCoV-infected-individuals>

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
cell_type_genes_temp=list("Bowman’s glands"=c("SOX9", "GPX3"),
                     "olfactory HBCs"=c("CXCL14", "MEG3"),
                     "olfactory ensheathing glia"=c("S100B", "PLP1"),
                     "olfactory microvillar cells"=c( "ASCL3", "HEPACAM2"),
                     "immature neurons"=c("GNG8","OLIG2"),
                     "mature neurons"=c("GNG13","RTP1"),
                     "GBCs"=c("CXCR4", "SOX2"),
                     "sustentacular cells"=c("CYP2A13", "ERMN"))
marker_genes=c("SOX9", "SOX10", "MUC5", "GPX3",
               "TP63", "KRT5", "CXCL14", "SOX2", "MEG3",
               "S100B", "PLP1", "PMP2","MPZ", "ALX3",
               "ASCL3", "CFTR", "HEPACAM2", "FOXL1",
               "GNG8", "OLIG2", "EBF2", "LHX2", "CBX8",
               "GNG13", "EBF2", "CBX8", "RTP1",
               "HES6", "ASCL1", "CXCR4", "SOX2", "EZH2", "NEUROD1", "NEUROG1",
               "CYP2A13", "CYP2J2", "GPX6", "ERMN", "SOX2")
marker_genes_temp=c("SOX9", "GPX3","CXCL14", "MEG3","S100B", "PLP1", "ASCL3", "HEPACAM2",
                    "GNG8","OLIG2","GNG13","RTP1","CXCR4", "SOX2","CYP2A13", "ERMN")
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

### Intigration and batch effect removing using Seurat

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
sce_3_1_1_after=run.combined
saveRDS(sce_3_1_1_after,"~/corona_project/sce_3_1_1_after_5000_23.rds")
```

### Dim plot to show clusters clearly for each cell type

``` r
pdf(paste(folder,"tSNE_plots.pdf",sep=""))
DimPlot(run.combined, reduction = "umap")
dev.off()
```

### FeaturePlot to show the expression of ACE2 and TMPRSS2 in each cell type

``` r
pdf(paste(folder,paste("gene_markers_","ACE2",".pdf",sep="")))
FeaturePlot(run.combined, features = "ACE2",combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order =TRUE)
dev.off()
pdf(paste(folder,paste("gene_markers_","TMPRSS2",".pdf",sep="")))
FeaturePlot(run.combined, features = "TMPRSS2",combine = FALSE,cols=c("lightgrey", "red"),pt.size=2,order = TRUE)
dev.off()
```

### Percentage of ACE2 expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
print(sum(a))
a=a*100/sum(a)
df=data.frame(cell_type=names(a),count_percentage=a)
pdf(paste(folder,"bar_ACE2.pdf",sep=""))
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dev.off()
dim(run.combined@assays$RNA)
```

### Percentage of TMPRSS2 expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])])
print(sum(a))
a=a*100/sum(a)
df=data.frame(cell_type=names(a),count_percentage=a)
pdf(paste(folder,"bar_TMPRSS2.pdf",sep=""))
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dev.off()
dim(run.combined@assays$RNA)
```

### bar plot of union of both ACE2 and TMPRSS2

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_ACE2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])]
union_names=union(names(col_name_ACE2),names(col_name_TMPRSS2))
union_cells=table(Seurat::Idents(run.combined)[union_names])
TMPRSS2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])])
ACE2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
new_data=data.frame("cells_percenrage"=c(((TMPRSS2_cells/union_cells)*100),((ACE2_cells/union_cells)*100)),
                    genes=c(rep("TMPRSS2",length(TMPRSS2_cells)),rep("ACE2",length(ACE2_cells))),
                    cell_types=c(names(TMPRSS2_cells),names(ACE2_cells)))
pdf(paste(folder,"bar_TMPRSS2_ACE2.pdf",sep=""))
ggplot(data=new_data, aes(x=cell_types, y=cells_percenrage, fill=genes)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
```

### with common cells also

``` r
common_names=intersect(names(col_name_ACE2),names(col_name_TMPRSS2))
common_cells=table(Seurat::Idents(run.combined)[common_names])
svg(paste(folder,"bar_TMPRSS2_ACE2.svg",sep=""))
new_data=data.frame("cells_percenrage"=c(((TMPRSS2_cells/union_cells)*100),((ACE2_cells/union_cells)*100),((common_cells/union_cells)*100)),
                    genes=c(rep("TMPRSS2",length(TMPRSS2_cells)),rep("ACE2",length(ACE2_cells)),rep("common",length(common_cells))),
                    cell_types=c(names(TMPRSS2_cells),names(ACE2_cells),names(common_cells)))
pdf(paste(folder,"bar_TMPRSS2_ACE2_commom.pdf",sep=""))
ggplot(data=new_data, aes(x=cell_types, y=cells_percenrage, fill=genes)) + 
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
```

### Dot plot cluster wise with padj values of ACE2 and TMPRSS2

``` r
pdf(paste(folder,"dot_TMPRSS2_ACE2.pdf",sep=""))
DotPlot(run.combined,features = c("ACE2","TMPRSS2"),cols = c("lightgrey", "red",dot.scale = 10))
dev.off()
```

### Finding differential expressed genes (HBC+;ACE2+) vs B (HBC+;ACE2-)

``` r
DefaultAssay(run.combined) <- "RNA"
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

## Supplementary figures

### Do heat map

``` r
sce_3_1_1_after_5000_23 <- readRDS("~/corona_project/sce_3_1_1_after_5000_23.rds")
cluster.averages <- AverageExpression(sce_3_1_1_after_5000_23, return.seurat = TRUE)
DoHeatmap(sce_3_1_1_after_5000_23,features =c(marker_genes_temp),size = 3, draw.lines = FALSE )
DoHeatmap(cluster.averages,features =c(marker_genes),size = 3, draw.lines = FALSE )
```

### Dim plot to show clusters clearly for each cell type

``` r
plot=list()
for(i in 1:length(marker_genes_temp))
{
  print(i)
  if(i %in% c(1,5,9))
    plot[[i]]<-VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]]) +
       NoLegend() +   theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank(),
                                                                 axis.line.x = element_blank())
  if(i %in% c(14,15,16))
    plot[[i]]<-VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]]) + 
        NoLegend() +   theme(axis.title.y=element_blank(),
                                                                 axis.text.y=element_blank(),
                                                                 axis.ticks.y=element_blank(),
                                                                 axis.line.y = element_blank())
  if(i %in% c(13))
    plot[[i]]<-VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]]) + 
      NoLegend() 
  if(i %in% c(2,3,4,6,7,8,10,11,12))
    plot[[i]]<-VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]]) + 
      NoLegend()+   theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.line.x = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.line.y = element_blank())
  
}
cowplot::plot_grid(plotlist = plot,nrow=4,rel_heights = c(1,1,1,2))
```

### for FeaturePlot

``` r
plot=list()
plot[[1]]=DimPlot(sce_3_1_1_after_5000_23, reduction = "umap",pt.size = 1,label=TRUE)+
   labs(title = "") + NoAxes() + NoLegend()
for(i in 1:length(marker_genes_temp))
{
  print(i)
  plot[[i+1]]<-FeaturePlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]],cols=c("lightgrey", "red"),
            pt.size=1,order = TRUE) +  NoLegend() + NoAxes()

}
AT=cowplot::plot_grid(plotlist = plot[c(3,4,5)],ncol=3)
AB=cowplot::plot_grid(plotlist = plot[c(12,13,11)],ncol=3)
AR=cowplot::plot_grid(plotlist = plot[c(6,7,8,9,10)],nrow=5)
AL=cowplot::plot_grid(plotlist = plot[c(2,14,15,16,17)],nrow=5)
AC=cowplot::plot_grid(plotlist = plot[c(1)],nrow=1)
ATCB=cowplot::plot_grid(plotlist = list(AT,AC,AB),nrow=3,rel_heights = c(1,4,1))
cowplot::plot_grid(plotlist = list(AL,ATCB,AR),ncol=3,rel_widths = c(1,2,1))
```
