library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(stringr)
library(readxl)
counts <- Read10X_h5("AM1.h5")
AM1 <- CreateSeuratObject(counts = counts,project = "Acral1")
counts <- Read10X_h5("AM2.h5")
AM2 <- CreateSeuratObject(counts = counts,project = "Acral2")
counts <- Read10X_h5("AM3pre.h5")
AM3pre <- CreateSeuratObject(counts = counts,project = "Acral3pre")
counts <- Read10X_h5("AM3post.h5")
AM3post <- CreateSeuratObject(counts = counts,project = "Acralpost")
counts <- Read10X_h5("AM4.h5")
AM4 <- CreateSeuratObject(counts = counts,project = "Acral4")
counts <- Read10X_h5("AM5.h5")
AM5 <- CreateSeuratObject(counts = counts,project = "Acral5")
counts <- Read10X_h5("AM6.h5")
AM6 <- CreateSeuratObject(counts = counts,project = "Acral6")
counts <- Read10X_h5("CM1.h5")
CM1 <- CreateSeuratObject(counts = counts,project = "Cutaneous1")
counts <- Read10X_h5("CM2.h5")
CM2 <- CreateSeuratObject(counts = counts,project = "Cutaneous2")
counts <- Read10X_h5("CM3.h5")
CM3 <- CreateSeuratObject(counts = counts,project = "Cutaneous3")
counts <- Read10X_h5("CM1lym.h5")
CM1lym <- CreateSeuratObject(counts = counts,project = "Cutaneous1lym")
#merge_data#
merged_seurat <- merge(AM1, y = c(AM2,AM3pre,AM3post,AM4,AM5,AM6,CM1,CM2,CM3,CM1lym ),
                       add.cell.ids = c("AM1","AM2","AM3pre","AM3post","AM4","AM5","AM6",
                                        "CM1","CM2","CM3","CM1lym"))
#count
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
rbc.genes <- c("HBA1","HBA2","HBB","HBD","HBM","HBQ1","AHSP")
merged_seurat[["percent.rbc"]] <- PercentageFeatureSet(merged_seurat,features = rbc.genes)

#creat_metadata_include_base&sample&group_information#
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "AM1"))] <- 'AM1'
metadata$sample[which(str_detect(metadata$cells, "AM2"))] <- 'AM2'
metadata$sample[which(str_detect(metadata$cells, "AM3pre"))] <- 'AM3pre'
metadata$sample[which(str_detect(metadata$cells, "AM3post"))] <- 'AM3post'
metadata$sample[which(str_detect(metadata$cells, "AM4"))] <- 'AM4'
metadata$sample[which(str_detect(metadata$cells, "AM5"))] <- 'AM5'
metadata$sample[which(str_detect(metadata$cells, "AM6"))] <- 'AM6'
metadata$sample[which(str_detect(metadata$cells, "CM1"))] <- 'CM1'
metadata$sample[which(str_detect(metadata$cells, "CM2"))] <- 'CM2'
metadata$sample[which(str_detect(metadata$cells, "CM3"))] <- 'CM3'
metadata$sample[which(str_detect(metadata$cells, "CM1lym"))] <- 'CM1lym'
metadata <- metadata %>%
  mutate(disease = case_when(
    sample %in% c("AM1","AM2","AM3pre","AM3post","AM4","AM5","AM6") ~ 'Acral',
    sample %in% c("CM1","CM2","CM3","CM1lym") ~ 'Cutaneous'
  ))
merged_seurat@meta.data <- metadata
###QC###
filtered_seurat <- subset(x = merged_seurat,
                          subset= (nUMI >= 1000) &
                            (nGene >= 200) &
                            (nGene <= 5000) &
                            (percent.rbc < 3) &
                            (mitoRatio < 0.15))
metadata <- filtered_seurat@meta.data
#QC_Visualize#
metadata %>% ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))+
  ggtitle("NCells")
metadata %>% ggplot(aes(color=sample, x=nUMI, fill= sample))+
  geom_density(alpha = 0.2) + scale_x_log10()+
  theme_classic() +ylab("Cell density")+geom_vline(xintercept = 500)
metadata %>% ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) + theme_classic() +
  scale_x_log10() + geom_vline(xintercept = 300)
metadata %>% ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +ggtitle("NCells vs NGenes")
metadata %>% ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() + scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept = 500) +geom_hline(yintercept = 250) +facet_wrap(~sample)
VlnPlot(filtered_seurat,features = c("nFeature_RNA","nCount_RNA","mitoRatio","percent.rbc"),
        group.by = "sample",pt.size = 0)

#normalize data#
scRNAlist_merge <- NormalizeData(filtered_seurat)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge,nfeatures = 2000)
scRNAlist_merge <- ScaleData(scRNAlist_merge,vars.to.regress = c('mitoRatio'))
#run_PCA#
scRNAlist_merge <- RunPCA(scRNAlist_merge, npcs = 50)
ElbowPlot(scRNAlist_merge, ndims = 50)
print(scRNAlist_merge[["pca"]], dims = 1:30, nfeatures = 5)
#visualize_check_data#
VizDimLoadings(scRNAlist_merge, dims = 1:5, reduction = "pca")
DimPlot(scRNAlist_merge, reduction = "pca")
DimHeatmap(scRNAlist_merge, dims = 1:20, cells = 500, balanced = TRUE)
#remove_Batch_effect_by_harmony#
library(harmony)
scRNA_harmony <- RunHarmony(object = scRNAlist_merge,
                            group.by.vars ="sample",
                            reduction = "pca",
                            dims.use = 1:30,
                            reduction.save = "harmony")
#joinLayers#
scRNA_harmony[['RNA']] <- JoinLayers(scRNA_harmony[['RNA']])
#chose_resolution#
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = 'harmony',dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = seq(from = 0.1,to = 1, by = 0.1))
library(clustree)
clustree(scRNA_harmony)
scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:30,reduction = 'harmony')
scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:30,reduction = 'harmony')
scRNA_harmony$RNA_snn_res.0.5
Idents(scRNA_harmony) <- 'RNA_snn_res.0.5'
DimPlot(scRNA_harmony,reduction = 'umap')
#UMAP_and_TSNE#
DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'sample',label = T)
#DimPlot(scRNA_harmony,reduction = 'tsne', label = T)
DimPlot(scRNA_harmony,reduction = 'umap', label = T)
DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'disease',label = T)

#save data

save(scRNA_harmony,file = 'scRNA_harmony2.Rdata')
########Run analysis#########
rm(list = ls())
gc()
#load packages#
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(dplyr)
###load_data###
load('scRNA_harmony2.Rdata')
scRNA_harmony$RNA_snn_res.0.5
Idents(scRNA_harmony) <- 'RNA_snn_res.0.5'
DimPlot(scRNA_harmony,reduction = 'umap')
umap_integrated_3 <- DimPlot(scRNA_harmony,reduction = 'umap', label = T,label.size = 6)+NoLegend()
umap_integrated_3
###Annotation_was_based_on_marker_genes_reported_in_the_original_publication###
current_levels <- levels(scRNA_harmony)
sorted_levels <- as.character(sort(as.numeric(current_levels)))
Idents(scRNA_harmony) <- factor(Idents(scRNA_harmony), levels = sorted_levels)
levels(scRNA_harmony)
current_levels <- levels(scRNA_harmony)
sorted_levels <- as.character(sort(as.numeric(current_levels)))
Idents(scRNA_harmony) <- factor(Idents(scRNA_harmony), levels = sorted_levels)
levels(scRNA_harmony)
new.cluster.ids <- c("Tumor_1",         #"SOX10"
                     "Tumor_2",         #"SOX10"
                     "CD8T_1",          #"CD8A"
                     "Tumor_3",         #"SOX10"
                     "Treg_1",           #"CD4"
                     "Tumor_4",         #"SOX10"
                     "Tumor_5",         #"SOX10"
                     "Tumor_6",         #"SOX10"
                     "Endothelial cell_1", #"CLDN5"
                     "Tumor_7",         #"SOX10"
                     "Fibroblast_1",         #"DCN"
                     "Tumor_8",         #"SOX10"
                     "Fibroblast_2",         #"DCN"
                     "B cell_1",         #"MS4A1"
                     "Tumor_9",         #"SOX10"
                     "Macrophage_1",         #"CD14","CD68"
                     "Tumor_10",         #"SOX10"
                     "Tumor_11",         #"SOX10"
                     "Endothelial cell_2", #"CLDN5"
                     "CD4T_1")            # "CD4"
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
scRNA_harmony<- AddMetaData(object = scRNA_harmony,
                            metadata = scRNA_harmony@active.ident,
                            col.name = "celltype")
p_umap <- DimPlot(scRNA_harmony,reduction = "umap",
                  group.by = "celltype",label = TRUE,repel = TRUE,          label.size = 7        ) +
  theme_classic() +
  theme(legend.position = "none",axis.title.x = element_text(size = 25,face = "plain"),
    axis.title.y = element_text(size = 25,face = "plain"),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.line = element_line(linewidth = 0.6),plot.title = element_blank())
p_umap
###########
Idents(scRNA_harmony) <- "celltype"
DefaultAssay(scRNA_harmony)
#DefaultAssay(scRNA_harmony) <- "RNA"
#markers_all <- FindAllMarkers(scRNA_harmony,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#library(dplyr)
#top10_markers <- markers_all %>%group_by(cluster) %>%slice_max(avg_log2FC, n = 10)
#write.csv(top10_markers,"Top10_Markers_each_cluster.csv",row.names = FALSE)
#top10_markers
##plot##
Idents(scRNA_harmony) <- "celltype"
celltypeorder <- c("Tumor_1", "Tumor_2", "Tumor_3", "Tumor_4", "Tumor_5", 
                   "Tumor_6", "Tumor_7", "Tumor_8", "Tumor_9", "Tumor_10", "Tumor_11",
                   "CD8T_1", "CD4T_1", "Treg_1", "B cell_1", "Macrophage_1",
                   "Endothelial cell_1", "Endothelial cell_2","Fibroblast_1", "Fibroblast_2")
Idents(scRNA_harmony) <- factor(Idents(scRNA_harmony),levels = celltypeorder)
cell_mark <- c("SEPT1", "SEPT2", "SEPT3", "SEPT4","SEPT5", "SEPT6", "SEPT7", "SEPT8",
               "SEPT9", "SEPT10", "SEPT11", "SEPT12","SEPT14",
               "STIM1", "STIM2","ORAI1", "ORAI2", "ORAI3")
cell_mark <- c("NCAM1",                   #NK
               "LILRA4",                  #DC
               "FOXP3",                   # Treg +CD3D
               "MS4A1",                   # B cell
               "CD3D","CD4",              # CD4T
               "CD8A",                    # CD8T
               "CD14","CD68",             # Mono/Macrophage
               "CLDN5",                   # Endothelial cell
               "DCN",                     # Fibroblast
               "SOX10",'PMEL','MLANA')                   # Tumor+NCAM1
DotPlot(scRNA_harmony,features = cell_mark,assay = "RNA",cluster.idents = F,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() +theme_bw() +labs(x = "",y = "") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 15),
    axis.text.y = element_text(size = 15),    
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)) +
  scale_color_gradient2(low = "#2166AC",mid = "white",high = "#B2182B")+
  guides(color = guide_colorbar(title = "Average\nExpression",title.position = "top"),
    size = guide_legend(title = "Percent\nExpressed",title.position = "top")) +
  theme(legend.title = element_text(hjust = 0.5,vjust = 0.5))
##cell ratio##
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
tumor_types <- c("Tumor_1", "Tumor_2", "Tumor_3", "Tumor_4","Tumor_5", 
                 "Tumor_6", "Tumor_7", "Tumor_8","Tumor_9", "Tumor_10", "Tumor_11")
tumor_seurat <- subset(scRNA_harmony,subset = celltype %in% tumor_types)
tumor_meta <- tumor_seurat@meta.data
tumor_count <- tumor_meta %>%count(sample, celltype, name = "n_cells")
tumor_count
tumor_prop <- tumor_count %>%group_by(sample) %>%
  mutate(proportion = n_cells / sum(n_cells),percentage = proportion * 100) %>%ungroup()
tumor_prop
sample_order <- c("AM1", "AM2", "AM3pre", "AM3post",
                  "AM4", "AM5", "AM6","CM1", "CM2", "CM3", "CM1lym")
tumor_prop$sample <- factor(tumor_prop$sample,levels = sample_order)
tumor_prop$celltype <- factor(tumor_prop$celltype,levels = tumor_types)
tumor_colors <- c("Tumor_1"  = "#5B88B9","Tumor_2"  = "#E6A154","Tumor_3"  = "#6AA577","Tumor_4"  = "#C96565",
                  "Tumor_5"  = "#8C78B5","Tumor_6"  = "#A9896A","Tumor_7"  = "#CF91B5","Tumor_8"  = "#9B9B9B",
                  "Tumor_9"  = "#D0BB6C","Tumor_10" = "#65B1C5","Tumor_11" = "#765B8E")
p_tumor_prop <- ggplot(tumor_prop,
                       aes(x = sample,y = proportion,fill = celltype)) +
  geom_bar(stat = "identity",width = 0.75,color = "black",linewidth = 0.3) +
  scale_fill_manual(values = tumor_colors,breaks = rev(tumor_types)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%"),expand = c(0, 0)) +
  labs(x = NULL, y = "Tumor cell proportion",fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 14),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.position = "right")
p_tumor_prop
###Just tumor###
cell_mark <- c("BANCR","COL28A1","HORMAD1","ISG15","RSAD2","IFIT3","IFIT1",
               "KIF20A","NEK2","CDC20","PLK1","TRPM1","TYRP1","FABP1","MET","DCT","PMEL")
TUMORcelltypeorder <- c("Tumor_1","Tumor_3","Tumor_8","Tumor_9","Tumor_10")
tumor_subset <- subset(scRNA_harmony,idents = TUMORcelltypeorder)
Idents(tumor_subset) <- factor(Idents(tumor_subset),levels = TUMORcelltypeorder)
DotPlot(tumor_subset,features = cell_mark,assay = "RNA",cluster.idents = FALSE,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() +theme_bw() +labs(x = "", y = "") +
  theme(axis.text.x = element_text(
      angle = 45,hjust = 1,vjust = 1,size = 15),
    axis.text.y = element_text(size = 15),    
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)) +
  scale_color_gradient2(low = "#2166AC",mid = "white",high = "#B2182B")+
  guides(color = guide_colorbar(title = "Average\nExpression",title.position = "top"),
         size = guide_legend(title = "Percent\nExpressed",title.position = "top")) +
  theme(legend.title = element_text(hjust = 0.5,vjust = 0.5))