###load packages
library(data.table)  
library(dplyr)       
library(ggplot2)     
library(pheatmap)    
library(DESeq2)      
library(edgeR)       
library(limma)       
library(tinyarray)  
library(readxl)
#remove.packages("tinyarray")
# install.packages("tinyarray")
######################################BEFORE############################################
##load data
library(readxl)
all_genesanalysis <- read_excel("HajduT_all_genes.xlsx",col_names = TRUE)
matchlist <- all_genesanalysis[, c("GeneSymbol", "GENE")]
write.csv(matchlist, "matchlist_GeneSymbol_GENES.csv", row.names = FALSE)
matchlist
#Because Gene symbols have duplicates, so they were first converted into indices#
library(dplyr)
expr_mean <- all_genesanalysis[, -1]
expr_mean <- as.data.frame(expr_mean)   
rownames(expr_mean) <- expr_mean$GENE
expr_mean <- expr_mean[ , -c(1)]
##Filter out low-expression genes##
keep <- rowSums(expr_mean >= 10) >= 3
expr_mean <- expr_mean[keep, ]
all_genesanalysis
##grouping##
expr_Melanocyte <- expr_mean[, 1:3]
expr_WM35 <- expr_mean[, 4:6]
expr_A2058 <- expr_mean[, 7:9]
expVSWM35 <- cbind(expr_Melanocyte,expr_WM35)
expVS2058 <- cbind(expr_Melanocyte,expr_A2058)
group <- c(rep('Melanocyte', ncol(expr_Melanocyte)), rep('WM35', ncol(expr_WM35)),
            rep('A2058', ncol(expr_A2058)))
group <- factor(group, levels = c("Melanocyte",'WM35', 'A2058'))
########################### Differ genes analysis prepare    ###########################
#######        DESeq2    #######
colData <- data.frame(row.names = colnames(expr_mean),group = group)
expr_meanDESeq <- apply(expr_mean, 2, as.integer)
rownames(expr_meanDESeq) <- rownames(expr_mean)
dds <- DESeqDataSetFromMatrix(countData = expr_meanDESeq, colData = colData,design = ~ group)  
dds <- DESeq(dds)
####!!!Obtain the VST matrix for PCA plotting and downstream rhythmic analysis#####
vsd <- vst(dds, blind=FALSE)
expr_for_plot <- assay(vsd)
expr_for_plot <- as.data.frame(expr_for_plot, stringsAsFactors = FALSE)
expr_for_plot$GENE <- rownames(expr_for_plot)
expr_for_plot_annot <- merge(expr_for_plot, matchlist, by = "GENE", all.x = TRUE)
expr_for_plot_annot <- expr_for_plot_annot[, c("GeneSymbol", colnames(expr_for_plot))]
expr_for_plot<-expr_for_plot_annot
write.csv(expr_for_plot, "expression_for_plot_vst.csv")
##BOXPLOT##
library(reshape2)
library(ggplot2)
plot_data <- expr_for_plot
plot_data <- plot_data[, !(colnames(plot_data) %in% c("GeneSymbol","GENE"))]
plot_data_long <- melt(as.matrix(plot_data))
colnames(plot_data_long) <- c("Gene","Sample","Expression")
plot_data_long$Group <- factor(c(rep("Melanocyte", 3),rep("WM35", 3),rep("A2058", 3))
                               [match(plot_data_long$Sample,colnames(plot_data))])
ggplot(plot_data_long,
       aes(x = Sample,y = Expression,fill = Group)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = c("Melanocyte" = "#F4B6B6","WM35" = "#B7E4C7","A2058" = "#B3D9FF"        
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1
    ),
    legend.position = "right"
  ) +
  labs(title = "Boxplot",x = "",y = "VST Expression")

###############################################PCA####################################
library(DESeq2)
library(ggplot2)
library(ggforce)
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
plotdat <- data.frame(Sample = colnames(vsd),PC1 = pca$x[,1],PC2 = pca$x[,2],group = group)
pc1_var <- percentVar[1]
pc2_var <- percentVar[2]
current_colors <- c(
  "Melanocyte" = "#F4B6B6",
  "WM35"       = "#B7E4C7",
  "A2058"      = "#B3D9FF"
)
current_treatment_shape <- 21
ggplot(plotdat, aes(x = PC1, y = PC2)) +
  ggforce::geom_mark_rect(
    aes(fill = group,color = group,group = group),
    expand = unit(3, "mm"),alpha = 0.35,size = 0.6,linetype = "solid",show.legend = FALSE,
    con.type = "none",label.buffer = unit(0, "mm")
  ) +
  geom_point(aes(fill = group),shape = current_treatment_shape,
             size = 5,colour = "black",stroke = 0.8
  ) +
  geom_hline(yintercept = 0,colour = "gray70",linetype = "dashed",linewidth = 0.5
  ) +
  geom_vline(xintercept = 0,colour = "gray70",linetype = "dashed",linewidth = 0.5
  ) +
  scale_fill_manual(values = current_colors,name = "Cell Line"
  ) +
  scale_color_manual(values = current_colors,guide = "none"
  ) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),y = paste0("PC2 (", pc2_var, "%)")
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_line(colour = "gray95",linewidth = 0.2
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14,face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = "white",color = "gray70",linewidth = 0.3))

#################################differentanalysis####################################
res_WM35 <- results(dds,contrast = c("group","WM35","Melanocyte"))
res_WM35 <- res_WM35[order(res_WM35$padj), ]
DEG_WM35 <- as.data.frame(res_WM35)
res_A2058 <- results(dds,contrast = c("group","A2058","Melanocyte"))
res_A2058 <- res_A2058[order(res_A2058$padj), ]
DEG_A2058 <- as.data.frame(res_A2058)
res_tumor <- results(dds,contrast = c("group","WM35","A2058"))
res_tumor <- res_tumor[order(res_tumor$padj), ]
DEG_tumor <- as.data.frame(res_tumor)
#
DEG_WM35$GENE <- rownames(DEG_WM35)
DEG_WM35_annot <- merge(DEG_WM35,matchlist,by = "GENE",all.x = TRUE)
DEG_WM35_annot <- DEG_WM35_annot[,c("GeneSymbol", colnames(DEG_WM35))]
DEG_A2058$GENE <- rownames(DEG_A2058)
DEG_A2058_annot <- merge(DEG_A2058,matchlist,by = "GENE",all.x = TRUE)
DEG_A2058_annot <- DEG_A2058_annot[,c("GeneSymbol", colnames(DEG_A2058))]
DEG_tumor$GENE <- rownames(DEG_tumor)
DEG_tumor_annot <- merge(DEG_tumor,matchlist,by = "GENE",all.x = TRUE)
DEG_tumor_annot <- DEG_tumor_annot[,c("GeneSymbol", colnames(DEG_tumor))]
write.csv(DEG_WM35_annot,"DEG_WM35_vs_Melanocyte.csv")
write.csv(DEG_A2058_annot,"DEG_A2058_vs_Melanocyte.csv")

####################################Gene###########################################
genelist <- c("SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN4","SEPTIN5", "SEPTIN6", "SEPTIN7", "SEPTIN8",
              "SEPTIN9", "SEPTIN10", "SEPTIN11", "SEPTIN12","SEPTIN14", 
              "STIM1", "STIM2","ORAI1", "ORAI2", "ORAI3")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
expr_gene <- expr_for_plot %>%dplyr::filter(GeneSymbol %in% genelist)
expr_gene_long <- expr_gene %>%pivot_longer(
  cols = -c(GeneSymbol, GENE),names_to = "Sample",values_to = "Expression")
expr_gene_long$Group <- dplyr::case_when(
  grepl("Melanocyte", expr_gene_long$Sample) ~ "Melanocyte",
  grepl("WM35", expr_gene_long$Sample) ~ "WM35",
  grepl("A2058", expr_gene_long$Sample) ~ "A2058")
expr_gene_long$Group <- factor(expr_gene_long$Group,
                               levels = c("Melanocyte","WM35","A2058"))
expr_gene_long$GeneSymbol <- factor(expr_gene_long$GeneSymbol,levels = genelist)
sig_genes_to_plot <- c("SEPTIN7","ORAI1","ORAI2","ORAI3","STIM1","STIM2")
sig_WM35 <- DEG_WM35_annot %>%dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%dplyr::distinct(GeneSymbol,.keep_all = TRUE) %>%dplyr::select(GeneSymbol,padj)
colnames(sig_WM35)[2] <- "padj_WM35"
sig_A2058 <- DEG_A2058_annot %>%dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%dplyr::distinct(GeneSymbol,.keep_all = TRUE) %>%
  dplyr::select(GeneSymbol,padj)
colnames(sig_A2058)[2] <- "padj_A2058"
sig_tumor <- DEG_tumor_annot %>%dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%dplyr::distinct(GeneSymbol,.keep_all = TRUE) %>%
  dplyr::select(GeneSymbol,padj)
colnames(sig_tumor)[2] <- "padj_tumor"

sig_df <- sig_WM35 %>%dplyr::left_join(
  sig_A2058,by = "GeneSymbol") %>%
  dplyr::left_join(sig_tumor,by = "GeneSymbol")
sig_df$WM35_sig <- ifelse(!is.na(sig_df$padj_WM35) &sig_df$padj_WM35 < 0.05,"*",NA)
sig_df$A2058_sig <- ifelse(!is.na(sig_df$padj_A2058) &sig_df$padj_A2058 < 0.05,"*",NA)
sig_df$tumor_sig <- ifelse(!is.na(sig_df$padj_tumor) &sig_df$padj_tumor < 0.05,"*",NA)
stat_df <- bind_rows(
  # Melanocyte vs WM35
  sig_df %>%dplyr::filter(!is.na(WM35_sig)) %>%
    dplyr::transmute(GeneSymbol,group1 = "Melanocyte",group2 = "WM35",label = "*"),
  # Melanocyte vs A2058
  sig_df %>%dplyr::filter(!is.na(A2058_sig)) %>%
    dplyr::transmute(GeneSymbol,group1 = "Melanocyte",group2 = "A2058",label = "*"),
  # WM35 vs A2058
  sig_df %>%dplyr::filter(!is.na(tumor_sig)) %>%
    dplyr::transmute(GeneSymbol,group1 = "WM35",group2 = "A2058",label = "*"))
### Y position ####
y_pos <- expr_gene_long %>%dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(ymin = min(Expression,na.rm = TRUE),ymax = max(Expression,na.rm = TRUE),.groups = "drop")
stat_df <- stat_df %>%dplyr::left_join(y_pos,by = "GeneSymbol")
stat_df <- stat_df %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::mutate(y.position = ymin +
                  (ymax - ymin) *
                  c(1.30,1.40,1.49)[seq_len(dplyr::n())]
  ) %>%dplyr::ungroup()
stat_df$GeneSymbol <- factor(stat_df$GeneSymbol,levels = genelist)
### Colors ###
mycols <- c("Melanocyte" = "#E57373","WM35"= "#66BB6A","A2058" = "#64B5F6")
y_axis_text_size <- 14
y_axis_title_size <- 18
#################################### Plot ##########################################
ggplot(expr_gene_long,aes(x = Group,y = Expression,fill = Group)
) +
  geom_violin(trim = FALSE,width = 1,alpha = 1,colour = "black",linewidth = 0.3
  ) +
  geom_boxplot(width = 0.16,fill = "white",colour = "black",linewidth = 0.45,outlier.shape = NA
  ) +
  stat_pvalue_manual(stat_df,label = "label",xmin = "group1",xmax = "group2",
                     y.position = "y.position",tip.length = 0.01,bracket.size = 0.5,size = 5
  ) +
  facet_wrap(~ GeneSymbol,scales = "free_y",ncol = 6
  ) +
  scale_fill_manual(values = mycols
  ) +
  scale_y_continuous(breaks = function(x) {
    seq(floor(min(x)),ceiling(max(x)),by = 1)},
    labels = function(x) {as.character(round(x))},
    expand = expansion(mult = c(0, 0),add = c(0, 0.5))
  ) +
  labs(x = NULL,y = "VST Expression"
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey90",colour = "grey60"),
        strip.text = element_text(size = 13),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = y_axis_text_size,colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = y_axis_title_size),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.size = unit(0.9,"cm"),
        legend.key = element_rect(fill = "white",colour = NA)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 22,size = 8,colour = "black",alpha = 1)))
###siSEPTIN7###
#expression#
head(Melanocyteexpr_for_plot, 5)
head(WM35expr_for_plot, 5)
head(A2058expr_for_plot, 5)
#DEGs#
head(DEG_Melanocyte_annot, 5)
head(DEG_WM35_annot, 5)
head(DEG_A2058_annot, 5)
genelist <- c("SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN4","SEPTIN5", "SEPTIN6", "SEPTIN7", "SEPTIN8",
              "SEPTIN9", "SEPTIN10", "SEPTIN11", "SEPTIN12","SEPTIN14", "STIM1", "STIM2","ORAI1", "ORAI2", "ORAI3")
sig_genes_to_plot <- c("SEPTIN7")
#################################### Packages #######################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(grid)
mycols <- c("Melanocyte NT"= "grey70","Melanocyte siSEPTIN7"  = "#1F4E79",
            "WM35 NT"= "grey70","WM35 siSEPTIN7"= "#1F4E79",
            "A2058 NT"= "grey70","A2058 siSEPTIN7"       = "#1F4E79")
#################################### Melanocyte ######################################
expr_Melanocyte <- Melanocyteexpr_for_plot %>%dplyr::filter(GeneSymbol %in% genelist)
expr_Melanocyte_long <- expr_Melanocyte %>%
  tidyr::pivot_longer(cols = -c(GeneSymbol, GENE),names_to = "Sample",values_to = "Expression")
expr_Melanocyte_long$Group <- dplyr::case_when(
  grepl("_NT\\(raw\\)$",expr_Melanocyte_long$Sample) ~ "Melanocyte NT",
  grepl("_siSEPT7\\(raw\\)$",expr_Melanocyte_long$Sample) ~ "Melanocyte siSEPTIN7",
  TRUE ~ NA_character_)
expr_Melanocyte_long$Group <- factor(expr_Melanocyte_long$Group,
                                     levels = c("Melanocyte NT","Melanocyte siSEPTIN7"))
expr_Melanocyte_long$GeneSymbol <- factor(as.character(expr_Melanocyte_long$GeneSymbol),levels = genelist)
expr_Melanocyte_long$GeneSymbol <- droplevels(expr_Melanocyte_long$GeneSymbol)
cat("\n================ Genes in Data ====================\n")
print(unique(as.character(expr_Melanocyte_long$GeneSymbol)))
##Significance ##
sig_Melanocyte <- DEG_Melanocyte_annot %>%dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%dplyr::distinct(GeneSymbol,.keep_all = TRUE) %>%
  dplyr::filter(!is.na(padj),padj < 0.05) %>%
  dplyr::transmute(GeneSymbol = GeneSymbol,group1 = "Melanocyte NT",group2 = "Melanocyte siSEPTIN7",label = "*")
if (nrow(sig_Melanocyte) > 0) {sig_Melanocyte$GeneSymbol <- factor(as.character(sig_Melanocyte$GeneSymbol),
                                                                   levels = genelist)
y_Melanocyte <- expr_Melanocyte_long %>%dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(ymax = max(Expression,na.rm = TRUE),ymin = min(Expression,na.rm = TRUE),.groups = "drop")
sig_Melanocyte <- sig_Melanocyte %>%dplyr::left_join(y_Melanocyte,by = "GeneSymbol") %>%
  dplyr::mutate(y.position = ymax + 0.30) %>%dplyr::select(GeneSymbol,group1,group2,label,y.position)
}
##Plot##
p_Melanocyte <- ggplot(expr_Melanocyte_long,aes(x = Group,y = Expression,fill = Group)
) +
  geom_violin(trim = FALSE,width = 1,alpha = 1,colour = "black",linewidth = 0.3
  ) +
  geom_boxplot(width = 0.16,fill = "white",colour = "black",linewidth = 0.45,outlier.shape = NA
  ) +
  {
    if (nrow(sig_Melanocyte) > 0) {
      stat_pvalue_manual(sig_Melanocyte,label = "label",xmin = "group1",xmax = "group2",
                         y.position = "y.position",tip.length = 0.01,bracket.size = 0.5,size = 5)
    }
  } +
  facet_wrap(~ GeneSymbol,scales = "free_y",ncol = 6,drop = TRUE
  ) +
  scale_fill_manual(values = mycols,breaks = c("Melanocyte NT","Melanocyte siSEPTIN7")
  ) +
  scale_y_continuous(breaks = function(x) {seq(from = min(x, na.rm = TRUE),to = max(x, na.rm = TRUE),length.out = 3)},
                     labels = function(x) {sprintf("%.1f", x)},
                     expand = expansion(mult = c(0, 0), add = c(0, 0.6))
  ) +
  labs(title = NULL,x = NULL,y = "VST Expression",fill = "Melanocyte"
  ) +
  theme_bw(base_size = 14
  ) +
  theme(plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90",colour = "grey60"),
        strip.text = element_text(
          size = 13),
        panel.grid.major = element_line(colour = "grey90",linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14,colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.key.size = unit(0.9,"cm"),
        legend.key = element_rect(fill = "white",colour = NA)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 22,size = 8,colour = "black",alpha = 1)))
p_Melanocyte
#################################### WM35 ######################################
expr_WM35 <- WM35expr_for_plot %>%dplyr::filter(GeneSymbol %in% genelist)
expr_WM35_long <- expr_WM35 %>%tidyr::pivot_longer(cols = -c(GeneSymbol, GENE),
                                                   names_to = "Sample",values_to = "Expression")
expr_WM35_long$Group <- dplyr::case_when(
  grepl("_non-targeting_[0-9]+\\(raw\\)$", expr_WM35_long$Sample) ~ "WM35 NT",
  grepl("_siSept7_[0-9]+\\(raw\\)$", expr_WM35_long$Sample) ~ "WM35 siSEPTIN7",
  TRUE ~ NA_character_)
expr_WM35_long$Group <- factor(expr_WM35_long$Group,levels = c("WM35 NT", "WM35 siSEPTIN7"))
expr_WM35_long$GeneSymbol <- factor(as.character(expr_WM35_long$GeneSymbol),levels = genelist)
expr_WM35_long$GeneSymbol <- droplevels(expr_WM35_long$GeneSymbol)
cat("\n================ WM35 Genes in Data ====================\n")
print(unique(as.character(expr_WM35_long$GeneSymbol)))
## Significance ##
sig_WM35 <- DEG_WM35_annot %>%
  dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(pvalue), pvalue < 0.05) %>%
  dplyr::transmute(
    GeneSymbol = GeneSymbol,group1 = "WM35 NT",group2 = "WM35 siSEPTIN7",label = "*")
if (nrow(sig_WM35) > 0) {
  sig_WM35$GeneSymbol <- factor(as.character(sig_WM35$GeneSymbol),levels = genelist)
  y_WM35 <- expr_WM35_long %>%dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(ymax = max(Expression, na.rm = TRUE),ymin = min(Expression, na.rm = TRUE),.groups = "drop")
  sig_WM35 <- sig_WM35 %>%
    dplyr::left_join(y_WM35, by = "GeneSymbol") %>%
    dplyr::mutate(y.position = ymax + 0.30) %>%
    dplyr::select(GeneSymbol,group1,group2,label,y.position)}
## Plot ##
p_WM35 <- ggplot(expr_WM35_long,aes(x = Group, y = Expression, fill = Group)
) +
  geom_violin(trim = FALSE,width = 1,alpha = 1,colour = "black",linewidth = 0.3
  ) +
  geom_boxplot(width = 0.16,fill = "white",colour = "black",linewidth = 0.45,outlier.shape = NA
  ) +
  {if (nrow(sig_WM35) > 0) {stat_pvalue_manual(
    sig_WM35,label = "label",xmin = "group1",xmax = "group2",y.position = "y.position",
    tip.length = 0.01,bracket.size = 0.5,size = 5)}
  } +
  facet_wrap(~ GeneSymbol,scales = "free_y",ncol = 6,drop = TRUE
  ) +
  scale_fill_manual(values = mycols,breaks = c("WM35 NT", "WM35 siSEPTIN7")
  ) +
  scale_y_continuous(breaks = function(x) {seq(from = min(x, na.rm = TRUE),to = max(x, na.rm = TRUE),length.out = 3)},
                     labels = function(x) {sprintf("%.1f", x)},
                     expand = expansion(mult = c(0, 0), add = c(0, 0.6))
  ) +
  labs(title = NULL,x = NULL,y = "VST Expression",fill = "WM35"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90",colour = "grey60"),
        strip.text = element_text(size = 13),
        panel.grid.major = element_line(colour = "grey90",linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14,colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.key.size = unit(0.9,"cm"),
        legend.key = element_rect(fill = "white",colour = NA)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 22,size = 8,colour = "black",alpha = 1)))
p_WM35
#################################### A2058 ######################################
expr_A2058 <- A2058expr_for_plot %>%dplyr::filter(GeneSymbol %in% genelist)
expr_A2058_long <- expr_A2058 %>%tidyr::pivot_longer(cols = -c(GeneSymbol, GENE),names_to = "Sample",values_to = "Expression")
expr_A2058_long$Group <- dplyr::case_when(
  grepl("_non-targeting_[0-9]+\\(raw\\)$", expr_A2058_long$Sample) ~ "A2058 NT",
  grepl("_siSept7_[0-9]+\\(raw\\)$", expr_A2058_long$Sample) ~ "A2058 siSEPTIN7",
  TRUE ~ NA_character_)
expr_A2058_long$Group <- factor(expr_A2058_long$Group,levels = c("A2058 NT", "A2058 siSEPTIN7"))
expr_A2058_long$GeneSymbol <- factor(as.character(expr_A2058_long$GeneSymbol),levels = genelist)
expr_A2058_long$GeneSymbol <- droplevels(expr_A2058_long$GeneSymbol)
cat("\n================ A2058 Genes in Data ====================\n")
print(unique(as.character(expr_A2058_long$GeneSymbol)))
## Significance ##
sig_A2058 <- DEG_A2058_annot %>%
  dplyr::filter(GeneSymbol %in% sig_genes_to_plot) %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(pvalue), pvalue < 0.05) %>%
  dplyr::transmute(GeneSymbol = GeneSymbol,group1 = "A2058 NT",group2 = "A2058 siSEPTIN7",label = "*")
if (nrow(sig_A2058) > 0) {sig_A2058$GeneSymbol <- factor(as.character(sig_A2058$GeneSymbol),levels = genelist)
y_A2058 <- expr_A2058_long %>%dplyr::group_by(GeneSymbol) %>%dplyr::summarise(
  ymax = max(Expression, na.rm = TRUE),ymin = min(Expression, na.rm = TRUE),.groups = "drop")
sig_A2058 <- sig_A2058 %>%dplyr::left_join(y_A2058,by = "GeneSymbol") %>%
  dplyr::mutate(y.position = ymax + 0.30) %>%
  dplyr::select(GeneSymbol,group1,group2,label,y.position)}
## Plot ##
p_A2058 <- ggplot(expr_A2058_long,aes(x = Group, y = Expression, fill = Group)
) +
  geom_violin(trim = FALSE,width = 1,alpha = 1,colour = "black",linewidth = 0.3
  ) +
  geom_boxplot(width = 0.16,fill = "white",colour = "black",linewidth = 0.45,outlier.shape = NA
  ) +
  {if (nrow(sig_A2058) > 0) {stat_pvalue_manual(
    sig_A2058,label = "label",xmin = "group1",xmax = "group2",y.position = "y.position",
    tip.length = 0.01,bracket.size = 0.5,size = 5)}
  } +
  facet_wrap(~ GeneSymbol,scales = "free_y",ncol = 6,drop = TRUE
  ) +
  scale_fill_manual(values = mycols,breaks = c("A2058 NT", "A2058 siSEPTIN7")
  ) +
  scale_y_continuous(breaks = function(x) {seq(from = min(x, na.rm = TRUE),to = max(x, na.rm = TRUE),length.out = 3)},
                     labels = function(x) {sprintf("%.1f", x)},
                     expand = expansion(mult = c(0, 0), add = c(0, 0.6))
  ) +
  labs(title = NULL,x = NULL,y = "VST Expression",fill = "A2058"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90",colour = "grey60"),
        strip.text = element_text(size = 13),
        panel.grid.major = element_line(colour = "grey90",linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14,colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.key.size = unit(0.9,"cm"),
        legend.key = element_rect(fill = "white",colour = NA)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 22,size = 8,colour = "black",alpha = 1)))
p_A2058
############
size <- dev.size("in")

ggsave(
  "Melanocyte.png",
  plot = p_Melanocyte,
  width = 9.319444,
  height = 6.631944,
  units = "in",
  dpi = 300
)
ggsave(
  "WM35.png",
  plot = p_WM35,
  width = 9.319444,
  height = 6.631944,
  units = "in",
  dpi = 300
)
ggsave(
  "A2058.png",
  plot = p_A2058,
  width = 9.319444,
  height = 6.631944,
  units = "in",
  dpi = 300
)