
## The-cellular-basis-of-the-loss-of-smell-in-2019-nCoV-infected-individuals

A prominent clinical symptom of 2019-nCoV infection is hyposmia/anosmia
(decrease or loss of sense of smell), along with general symptoms such
as fatigue, shortness of breath, fever, and cough. The identity of the
cell lineages that underpin infection associated with loss of olfaction
could be critical for the diagnostics/clinical management of 2019-nCoV
infected individuals. Angiotensin I Converting Enzyme 2 (ACE2), and
Transmembrane protease serine 2 (TMPRSS2) are emerging as associated
host receptors vital for viral entry. Accordingly, the ongoing medical
examinations and the autopsy reports of the deceased individuals
strongly corroborate with the organ/cell-type-specific expression of
ACE2, TMPRSS2, and other putative viral entry-associated genes. To
determine the cellular basis of anosmia upon 2019-nCoV infection, we
employed a targeted bioinformatic analysis of single-cell expression
profiles of human olfactory epithelium cell-types. Our results
underscored selective expression of these viral entry-associated genes
in a subset of sustentacular cells, Bowman’s gland cells, and stem cells
of the olfactory epithelium. Co-expression analysis of ACE2 and TMPRSS2
and protein-protein interaction among the host and viral proteins
elected regulatory cytoskeleton protein-enriched sustentacular cells as
the most vulnerable cell-type of the olfactory epithelium. Furthermore,
expression, structural and docking analyses of these viral-entry
moieties revealed the potential risk of olfactory dysfunction in four
additional mammalian species, revealing an evolutionarily conserved
susceptibility. In summary, our findings suggest the molecular and
cellular rationale of loss of smell in 2019-nCoV infected
patients.

### Required data can can be downloaded by the link:

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522>

### To have a look on figures with code,you can download html file link given below.

<https://github.com/krishan57gupta/The-cellular-basis-of-the-loss-of-smell-in-2019-nCoV-infected-individuals>

## Libraries need to be loaded before running

``` r
library(ggplot2)
library(ggpubr)
library(Seurat)
library(cowplot)
library(readxl)
library(corto)
library(Matrix)
library(EnhancedVolcano)
```

## Initilization

``` r
set.seed(0)
label_size=8
point_size=1
cell_type_genes=list("Bowman’s glands"=c("SOX9", "SOX10", "GPX3"),
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
marker_genes=c("SOX9", "SOX10", "GPX3",
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
```

## Loading paitent P2 and P3 matrices and formation of Seurat object

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
```

## Integration and batch effect removing using Seurat

``` r
set.seed(1)
run.anchors <- FindIntegrationAnchors(object.list = batch_data_list, dims = 1:30,anchor.features = 5000)
a=run.anchors
run.combined <- IntegrateData(anchorset = run.anchors, dims = 1:30)
b=run.combined
sce_3_1_1_before=run.combined
saveRDS(sce_3_1_1_before,"~/corona_project/sce_3_1_1_before_5000_23.rds")
```

## Dim reduction and clustering

``` r
DefaultAssay(run.combined) <- "integrated"
run.combined <- ScaleData(run.combined, verbose = FALSE)
run.combined <- RunPCA(run.combined, npcs = 30, verbose = FALSE)
run.combined <- RunUMAP(run.combined, reduction = "pca", dims = 1:30)
run.combined <- FindNeighbors(run.combined, reduction = "pca", dims = 1:30)
run.combined <- FindClusters(run.combined, resolution = 0.5)
```

## Cell-types annotation

``` r
for(i in 0:as.integer(length(marker_genes)/4))
{
  a=(i*4+1)
  b=((i+1)*4)
  if(i==as.integer(length(marker_genes)/4))
    b=length(marker_genes)
  if(a>length(marker_genes))
     next
  print(i)
  print(VlnPlot(run.combined, features = marker_genes[a:b],ncol=2))
  print(FeaturePlot(run.combined, features = marker_genes[a:b],cols=c("lightgrey", "red"),label = TRUE,pt.size=2,order = TRUE))
}
run.combined<-RenameIdents(run.combined,`13`="Immature neurons", `14`="Mature neurons", `3`="Olfactory HBCs",`19`="Olfactory microvillar cells",
                   `5`="Bowman's gland", `16`="Olfactory ensheathing glia", `9`="Sustentacular cells",`20`="GBCs" )
run.combined=subset(run.combined, idents = c("Immature neurons","Mature neurons","Olfactory HBCs","Olfactory microvillar cells","Bowman's gland",
                          "Olfactory ensheathing glia","Sustentacular cells","GBCs"))
sce_3_1_1_after=run.combined
saveRDS(sce_3_1_1_after,"~/corona_project/sce_3_1_1_after_5000_23.rds")
```

## Ploting heat map for selected gene markers on all clusters

``` r
sce_3_1_1_after_5000_23 <- readRDS("~/corona_project/sce_3_1_1_after_5000_23.rds")
DoHeatmap(sce_3_1_1_after_5000_23,features =c(marker_genes_temp),size = 3, draw.lines = FALSE )
```

## Vln plots with all selected gene markers

``` r
plot=list()
for(i in 1:length(marker_genes_temp))
{
  print(i)
  if(i %in% c(1,5,9))
    plot[[i]]<-print(VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]])) +
      NoLegend() +   theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(),
                           axis.line.x = element_blank())
  if(i %in% c(14,15,16))
    plot[[i]]<-print(VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]])) +
      NoLegend() +   theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank(),
                           axis.line.y = element_blank())
  if(i %in% c(13))
    plot[[i]]<-print(VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]])) +
      NoLegend()
  if(i %in% c(2,3,4,6,7,8,10,11,12))
    plot[[i]]<-print(VlnPlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]])) +
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

## Feature plots with all selected gene markers

``` r
plot=list()
plot[[1]]=DimPlot(sce_3_1_1_after_5000_23, reduction = "umap",pt.size = point_size,label=TRUE)+
  labs(title = "") + NoAxes() + NoLegend()
for(i in 1:length(marker_genes_temp))
{
  print(i)
  plot[[i+1]]<-print(FeaturePlot(sce_3_1_1_after_5000_23, features = marker_genes_temp[[i]],cols=c("lightgrey", "red"),
                                 pt.size = point_size,label.size = label_size,order = TRUE)) +
    NoLegend() + NoAxes()
}
AT=cowplot::plot_grid(plotlist = plot[c(3,4,5)],ncol=3)
AB=cowplot::plot_grid(plotlist = plot[c(12,13,11)],ncol=3)
AR=cowplot::plot_grid(plotlist = plot[c(6,7,8,9,10)],nrow=5)
AL=cowplot::plot_grid(plotlist = plot[c(2,14,15,16,17)],nrow=5)
AC=cowplot::plot_grid(plotlist = plot[c(1)],nrow=1)
ATCB=cowplot::plot_grid(plotlist = list(AT,AC,AB),nrow=3,rel_heights = c(1,4,1))
cowplot::plot_grid(plotlist = list(AL,ATCB,AR),ncol=3,rel_widths = c(1,2,1))
```

## Dim plot to show all cell types with different colors

``` r
DimPlot(sce_3_1_1_after_5000_23, reduction = "umap",pt.size = point_size,label=TRUE)+ labs(title = "") + NoLegend()
DimPlot(sce_3_1_1_after_5000_23, reduction = "umap",pt.size = point_size,label=FALSE)+ labs(title = "")
```

## FeaturePlot to show expression of gene markers in respective cell markers

### FeaturePlot to show the expression of ACE2 in each cell type

``` r
FeaturePlot(run.combined, features = "ACE2",combine = FALSE,cols=c("lightgrey", "red"),pt.size = point_size,label.size = label_size,order =TRUE)
```

### FeaturePlot to show the expression of TMPRSS2 in each cell type

``` r
FeaturePlot(run.combined, features = "TMPRSS2",combine = FALSE,cols=c("lightgrey", "green"),pt.size = point_size,label.size = label_size,order =TRUE)
```

### FeaturePlot to show the expression of CTSL in each cell type

``` r
FeaturePlot(run.combined, features = "CTSL",combine = FALSE,cols=c("lightgrey", "blue"),pt.size = point_size,label.size = label_size,order =TRUE)
```

### FeaturePlot to show the expression of BSG in each cell type

``` r
FeaturePlot(run.combined, features = "BSG",combine = FALSE,cols=c("lightgrey", "magenta"),pt.size = point_size,label.size = label_size,order =TRUE)
```

## Initiallizing a matrix for 24 \* 8, 8 for cells type and 24 for different cases as shown below

``` r
count_matrix=matrix(0,24,8)
rownames(count_matrix)<-c("ACE2_count","ACE2_percent","TMPRSS2_count","TMPRSS2_percent",
                          "CTSL_count","CTSL_percent","BSG_count","BSG_percent",
                          "intersect_of_ACE2_TMPRSS2_count","ACE2_contri_percent","TMPRSS2_contri_percent",
                          "intersect_of_ACE2_CTSL_count","ACE2_contri_percent","CTSL_contri_percent",
                          "intersect_of_BSG_TMPRSS2_count","BSG_contri_percent","TMPRSS2_contri_percent",
                          "intersect_of_BSG_CTSL_count","BSG_contri_percent","CSTL_contri_percent",
                          "ACE2_TMPRSS2_CTSL","ACE2_TMPRSS2_CTSB",
                          "BSB_TMPRSS2_CTSL","BSG_TMPRSS2_CTSB")
colnames(count_matrix)<-unique(Idents(run.combined))
```

### Percentage of ACE2 expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
count_matrix[1,names(a)]=a
print(sum(a))
a=a*100/sum(a)
count_matrix[2,names(a)]=a
df=data.frame(cell_type=names(a),count_percentage=a)
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dim(run.combined@assays$RNA)
```

### Percentage of TMPRSS2 expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])])
count_matrix[3,names(a)]=a
print(sum(a))
a=a*100/sum(a)
count_matrix[4,names(a)]=a
df=data.frame(cell_type=names(a),count_percentage=a)
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dim(run.combined@assays$RNA)
```

### Percentage of CTSL expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])])
count_matrix[5,names(a)]=a
print(sum(a))
a=a*100/sum(a)
count_matrix[6,names(a)]=a
df=data.frame(cell_type=names(a),count_percentage=a)
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dim(run.combined@assays$RNA)
```

### Percentage of BSG expressed’s cells in each cell type

``` r
a=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])])
count_matrix[7,names(a)]=a
print(sum(a))
a=a*100/sum(a)
count_matrix[8,names(a)]=a
df=data.frame(cell_type=names(a),count_percentage=a)
ggplot(data=df, aes(x=cell_type, y=count_percentage.Freq, fill=cell_type)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
dim(run.combined@assays$RNA)
```

### Bar plot of intersect of both ACE2 and TMPRSS2

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_ACE2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])]
intersect_names=intersect(names(col_name_ACE2),names(col_name_TMPRSS2))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[9,names(intersect_cells)]=intersect_cells
TMPRSS2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])])
ACE2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
count_matrix[10,names(ACE2_cells)]=(intersect_cells/ACE2_cells)*100
count_matrix[11,names(TMPRSS2_cells)]=(intersect_cells/TMPRSS2_cells)*100
new_data=data.frame("cells_percenrage"=c(((intersect_cells/TMPRSS2_cells)*100),((intersect_cells/ACE2_cells)*100)),
                    genes=c(rep("TMPRSS2",length(TMPRSS2_cells)),rep("ACE2",length(ACE2_cells))),
                    cell_types=c(names(TMPRSS2_cells),names(ACE2_cells)))
ggplot(data=new_data, aes(x=cell_types, y=cells_percenrage, fill=genes)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar plot of intersect of both ACE2 and CTSL

``` r
col_name_CTSL=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])]
col_name_ACE2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])]
intersect_names=intersect(names(col_name_ACE2),names(col_name_CTSL))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[12,names(intersect_cells)]=intersect_cells
CTSL_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])])
ACE2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])])
count_matrix[13,names(ACE2_cells)]=(intersect_cells/ACE2_cells)*100
count_matrix[14,names(CTSL_cells)]=(intersect_cells/CTSL_cells)*100
new_data=data.frame("cells_percenrage"=c(((intersect_cells/CTSL_cells)*100),((intersect_cells/ACE2_cells)*100)),
                    genes=c(rep("CTSL",length(CTSL_cells)),rep("ACE2",length(ACE2_cells))),
                    cell_types=c(names(CTSL_cells),names(ACE2_cells)))
ggplot(data=new_data, aes(x=cell_types, y=cells_percenrage, fill=genes)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar plot of intersect of both BSG and TMPRSS2

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_BSG=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])]
intersect_names=intersect(names(col_name_BSG),names(col_name_TMPRSS2))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[15,names(intersect_cells)]=intersect_cells
TMPRSS2_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])])
BSG_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])])
count_matrix[16,names(BSG_cells)]=(intersect_cells/BSG_cells)*100
count_matrix[17,names(TMPRSS2_cells)]=(intersect_cells/TMPRSS2_cells)*100
new_data=data.frame("cells_percenrage"=c(((intersect_cells/TMPRSS2_cells)*100),((intersect_cells/BSG_cells)*100)),
                    genes=c(rep("TMPRSS2",length(TMPRSS2_cells)),rep("BSG",length(BSG_cells))),
                    cell_types=c(names(TMPRSS2_cells),names(BSG_cells)))
ggplot(data=new_data, aes(x=cell_types, y=cells_percenrage, fill=genes)) +
  geom_bar(stat="identity") +  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar plot of intersect of both BSG and CTSL

```r
col_name_CTSL=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])]
col_name_BSG=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])]
intersect_names=intersect(names(col_name_BSG),names(col_name_CTSL))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[18,names(intersect_cells)]=intersect_cells
CTSL_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])])
BSG_cells=table(Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])])
count_matrix[19,names(BSG_cells)]=(intersect_cells/BSG_cells)*100
count_matrix[20,names(CTSL_cells)]=(intersect_cells/CTSL_cells)*100
new_data=data.frame("cells_percenrage"=c(((intersect_cells/CTSL_cells)*100),((intersect_cells/BSG_cells)*100)),
genes=c(rep("CTSL",length(CTSL_cells)),rep("BSG",length(BSG_cells))),
cell_types=c(names(CTSL_cells),names(BSG_cells))) ggplot(data=new_data,
aes(x=cell_types, y=cells_percenrage, fill=genes)) +
geom_bar(stat="identity") + theme_classic()+ theme(axis.text.x =
element_text(angle = 90, hjust = 1))
```

### Intersect of ACE2, TMPRSS2 and CTSL

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_ACE2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])]
col_name_CTSL=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])]
intersect_names=intersect(names(col_name_ACE2),names(col_name_TMPRSS2))
intersect_names=intersect(intersect_names,names(col_name_CTSL))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[21,names(intersect_cells)]=intersect_cells
```

### Intersect of ACE2, TMPRSS2 and CTSB

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_ACE2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["ACE2",])>0])]
col_name_CTSB=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSB",])>0])]
intersect_names=intersect(names(col_name_ACE2),names(col_name_TMPRSS2))
intersect_names=intersect(intersect_names,names(col_name_CTSB))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[22,names(intersect_cells)]=intersect_cells
```

### Intersect of BSG, TMPRSS2 and CTSL

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["TMPRSS2",])>0])]
col_name_BSG=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["BSG",])>0])]
col_name_CTSL=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(run.combined@assays$RNA["CTSL",])>0])]
intersect_names=intersect(names(col_name_BSG),names(col_name_TMPRSS2))
intersect_names=intersect(intersect_names,names(col_name_CTSL))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[23,names(intersect_cells)]=intersect_cells
```

### Intersect of BSG, TMPRSS2 and CTSB

``` r
col_name_TMPRSS2=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(<run.combined@assays$RNA>[“TMPRSS2”,])>0])]
col_name_BSG=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(<run.combined@assays$RNA>[“BSG”,])>0])]
col_name_CTSB=Seurat::Idents(run.combined)[colnames(run.combined[,as.matrix(<run.combined@assays$RNA>[“CTSB”,])>0])]
intersect_names=intersect(names(col_name_BSG),names(col_name_TMPRSS2))
intersect_names=intersect(intersect_names,names(col_name_CTSB))
intersect_cells=table(Seurat::Idents(run.combined)[intersect_names])
count_matrix[24,names(intersect_cells)]=intersect_cells
count_matrix
```

## Bar plots to show + expressed cells of some gene markers as shown below

### Bar stack plot of ACE2, TMPRSS2 and (ACE2 and TMPRSS2)

``` r
count=c(count_matrix["ACE2_count",],count_matrix["TMPRSS2_count",],count_matrix["intersect_of_ACE2_TMPRSS2_count",])
count[count==0]<-1
A_T_AT_df=data.frame("count"=count,"cell_markers"= names(count),
                     "gene_markers"=c(rep("ACE2",8),rep("TMPRSS2",8),rep("Intersect (ACE2 and TMPRSS2)",8)))
ggplot(data=A_T_AT_df, aes(x=cell_markers, y=count, fill=gene_markers,width=.7)) +
  geom_bar(stat="identity")+
  theme_minimal() + scale_fill_manual(values=c('red','green',"blue")) +
  scale_y_continuous(trans='log10') +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar stack plot of BSG, CTSL and (BSG and CTSL)

``` r
count=c(count_matrix["BSG_count",],count_matrix["CTSL_count",],count_matrix["intersect_of_BSG_CTSL_count",])
count[count==0]<-1
B_C_BC_df=data.frame("count"=count,"cell_markers"= names(count),
                     "gene_markers"=c(rep("BSG",8),rep("CTSL",8),rep("Intersect (BSG and CTSL)",8)))
ggplot(data=B_C_BC_df, aes(x=cell_markers, y=count, fill=gene_markers,width=.7)) +
  geom_bar(stat="identity")+
  theme_minimal() + scale_fill_manual(values=c('red','green',"blue")) +
  scale_y_continuous(trans='log10') +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar stack plot of BSG, CTSL and (BSG and CTSL)

``` r
count=c(count_matrix["ACE2_count",],count_matrix["CTSL_count",],count_matrix["intersect_of_ACE2_CTSL_count",])
count[count==0]<-1
A_C_AC_df=data.frame("count"=count,"cell_markers"= names(count),
                     "gene_markers"=c(rep("ACE2",8),rep("CTSL",8),rep("Intersect (ACE2 and CTSL)",8)))
ggplot(data=A_C_AC_df, aes(x=cell_markers, y=count, fill=gene_markers,width=.7)) +
  geom_bar(stat="identity")+
  theme_minimal() + scale_fill_manual(values=c('red','green',"blue")) +
  scale_y_continuous(trans='log10') +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### Bar stack plot of BSG, TMPRSS2 and (BSG and TMPRSS2)

``` r
count=c(count_matrix["BSG_count",],count_matrix["TMPRSS2_count",],count_matrix["intersect_of_BSG_TMPRSS2_count",])
count[count==0]<-1
B_T_BT_df=data.frame("count"=count,"cell_markers"= names(count),
                     "gene_markers"=c(rep("BSG",8),rep("TMPRSS2",8),rep("Intersect (BSG and TMPRSS2)",8)))
ggplot(data=B_T_BT_df, aes(x=cell_markers, y=count, fill=gene_markers,width=.7)) +
  geom_bar(stat="identity")+
  theme_minimal() + scale_fill_manual(values=c('red','green',"blue")) +
  scale_y_continuous(trans='log10') +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

## DE genes analysis for selected cell types

``` r
run.combined=sce_3_1_1_after_5000_23
DefaultAssay(run.combined) <- "RNA"
seurat_subset_list=list()
p_val=0.05
l_fold=1
de_methods=c("poisson")
```

### DE genes analysis on intersection of (ACE2 and TMPRSS2)+ and (ACE2 and TMPRSS2)-

``` r
gene_marker=c("ACE2","TMPRSS2")
print(sort(count_matrix[1,]))
selected_cells_type=sort(count_matrix[9,])
selected_cells_type=selected_cells_type[selected_cells_type>2]
for(k in de_methods)
  for(i in 1:length(selected_cells_type))
  {
    print(paste(k,"_",i,sep=""))
    seurat_subset=subset(run.combined, idents = names(selected_cells_type)[i])
    positive_cells_name = colnames(seurat_subset[,(as.matrix(seurat_subset@assays$RNA[gene_marker[1],])>0 & as.matrix(seurat_subset@assays$RNA[gene_marker[2],])>0)])
    negative_cells_name = colnames(seurat_subset[,(as.matrix(seurat_subset@assays$RNA[gene_marker[1],])<=0 | as.matrix(seurat_subset@assays$RNA[gene_marker[2],])<=0)])
    pos_cells=paste(gene_marker[1],"_and_",gene_marker[2],"+",sep="")
    neg_cells=paste(gene_marker[1],"_and_",gene_marker[2],"-",sep="")
    pos=c(rep(pos_cells,length(positive_cells_name)))
    names(pos)=positive_cells_name
    neg=c(rep(neg_cells,length(negative_cells_name)))
    names(neg)=negative_cells_name
    pos_neg=c(pos,neg)
    print(paste("cells in seurat object after subsetting with ",names(selected_cells_type)[i]," = ",dim(seurat_subset)[2],sep=""))
    print(paste(names(selected_cells_type)[i]," total cells = ",sum(pos_neg[which(names(pos_neg)%in%rownames(seurat_subset@meta.data))]==pos_neg),sep=""))
    print(paste(pos_cells," cells = ",sum(pos_neg==pos_cells),sep=""))
    print(paste(neg_cells," cells = ",sum(pos_neg==neg_cells),sep=""))
    seurat_subset=AddMetaData(
      object = seurat_subset,
      metadata = pos_neg,
      col.name = 'seurat_subset'
    )
    seurat_subset_list[[paste(names(selected_cells_type)[i],"_",gene_marker[1],"_and_",gene_marker[2],sep="")]]<-seurat_subset
    Idents(seurat_subset)<-pos_neg
    seurat_subset_genes=FindMarkers(seurat_subset, ident.1=pos_cells,ident.2=neg_cells,test.use = k,min.cells.group=1)
    print(seurat_subset_genes[1:2,])
    seurat_subset_genes_up=rownames(seurat_subset_genes)[seurat_subset_genes$avg_logFC>l_fold & seurat_subset_genes$p_val<p_val]
    seurat_subset_genes_down=rownames(seurat_subset_genes)[seurat_subset_genes$avg_logFC<(-l_fold) & seurat_subset_genes$p_val<p_val]
    seurat_subset_genes_names=c(seurat_subset_genes_up,seurat_subset_genes_down)
    print(paste(names(selected_cells_type)[i],"_",pos_cells," genes = ",length(seurat_subset_genes_up),
                " and ",names(selected_cells_type)[i],"_",neg_cells," genes = ",length(seurat_subset_genes_down),sep=""))
    DE_info=data.frame("log2FC"=as.numeric(seurat_subset_genes$avg_logFC),"p_val"=as.numeric(seurat_subset_genes$p_val))
    lab=rownames(seurat_subset_genes)
    print(EnhancedVolcano(DE_info,
                          lab = lab,
                          x = 'log2FC',
                          y = 'p_val',
                          xlim = c(-5, 5),
                          title = names(selected_cells_type)[i],
                          pCutoff = p_val,
                          FCcutoff = l_fold,
                          colAlpha = 1))
  }
```

## Finding Stouffer score after log transfer and zscore

### Loading seaurat object

``` r
sce_3_1_1_after_5000_23 <- readRDS("~/corona_project/sce_3_1_1_after_5000_23.rds")
mat<-sce_3_1_1_after_5000_23@assays$RNA@counts
```

### Cell filtering

``` r
cells_sum <- Matrix::colSums(mat>=3)
mat<-mat[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]
```

### Gene filtering

``` r
mat<-mat[which(Matrix::rowSums(mat > 2) > 3),]
```

### Median normalization

``` r
cells_sum<-Matrix::rowSums(t(mat))
mat<-Matrix::t(t(mat)/(cells_sum/stats::median(cells_sum)))
```

### Loading media genes file and intersecting genes and filtering matrix common genes

``` r
media_genes <- read_excel("~/corona_project/media-6.xlsx")$...3[-1]
common_genes=intersect(rownames(mat),media_genes)
mat=mat[common_genes,]
```

### Further insuring have cells with atleast 10% expressed genes

``` r
mat<-mat[,Matrix::colSums(mat>0)>(dim(mat)[1]/10)]
print(dim(mat))
seurat_RNA_mat<-mat
print(sum(apply(seurat_RNA_mat,1,function(x) sum(x)==0)))
print(sum(apply(seurat_RNA_mat,2,function(x) sum(x)==0)))
print(dim(seurat_RNA_mat))
```

### \#\#\# First adding 1 (adding 1 to only zero) then log2 then zscore

``` r
seurat_RNA_mat[seurat_RNA_mat==0]=1
seurat_RNA_mat<-log2(seurat_RNA_mat)
seurat_RNA_mat<-t(apply(seurat_RNA_mat,1,function(x) (x-mean(x))/sd(x)))
```

### Stouffer score and then boxplot for each cell type

``` r
stouffer_score<- apply(seurat_RNA_mat,2,function(x) sum(x)/sqrt(length(x)))
print(sum(is.na(stouffer_score)))
cell_types<-Idents(sce_3_1_1_after_5000_23)[names(stouffer_score)]
unique_cell_types=unique(cell_types)
Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
  geom_boxplot(position=position_dodge(.2)) +
  geom_jitter(shape=16, position=position_jitter(.1)) +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

### One sided wilcoxon test

``` r
p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
colnames(p_val_stouffer_score)=unique_cell_types
rownames(p_val_stouffer_score)=unique_cell_types
for (i in 1:dim(p_val_stouffer_score)[1])
{
  for (j in 1:dim(p_val_stouffer_score)[2])
  {
    print(paste(i,"_",j,sep=""))
    a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
    b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
    p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
  }
}
p_val_stouffer_score
```

### Scatter plot and Pearson correlation with the average vectors from patient 2 and patient 3

#### 5 average vector for patient 2 and 5 from patient 3 as “ACE2+”,“TMPRSS2+”,“BSG+”,“CTSL+”,“All\_cells”

``` r
sce_3_1_1_after_5000_23 <- readRDS("~/corona_project/sce_3_1_1_after_5000_23.rds")
run=sce_3_1_1_after_5000_23@meta.data$run
for(i in c("ACE2","TMPRSS2","BSG","CTSL","All"))
{
  print(i)
  P2=sce_3_1_1_after_5000_23@assays$RNA@counts[,run=="P2"]
  dim(P2)
  if (i!="All")
    P2=P2[,P2[i,]>0]
  P2[P2==0]<-1
  P2=t(apply(P2,1,function(x) log2(x)))
  dim(P2)
  P2=rowMeans(as.matrix(P2))
  length(P2)
  P3=sce_3_1_1_after_5000_23@assays$RNA@counts[,run=="P3"]
  dim(P3)
  if (i!="All")
    P3=P3[,P3[i,]>0]
  P3[P3==0]<-1
  P3=t(apply(P3,1,function(x) log2(x)))
  dim(P3)
  P3=rowMeans(as.matrix(P3))
  length(P3)
  corr=cor(P2,P3,method="pearson")
  df=data.frame("P2"=P2,"P3"=P3)
  print(ggplot(df, aes(x=P2, y=P3)) + 
    geom_point(shape=18, color="blue")+
    geom_smooth(method=lm,  linetype="dashed",
                color="darkred", fill="blue")+
    ggtitle(paste("pearson correlation: ",corr,sep=""))+
    xlab("Patient 2")+ylab("Patient 3")+theme_classic())
}
```
