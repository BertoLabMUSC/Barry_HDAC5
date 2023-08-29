# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
library(sva)
library(viridis)
library(limma)
library(janitor)
library(UpSetR)
})


# Load data
load("dge/Sarah_Dge_Defined.RData")
load("dge/Expression_Object.RData")
load("dge/Sarah_Dge_Interaction_Posthoc.RData")

# Input for viz Plot for HDAC5c vs HDAC5s
dir.create("dge/HDAC5c_HDAC5s")
df <- HDAC5c_HDAC5s %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
           unite("Class", Genotype:Condition, na.rm = TRUE, remove = F) %>%
            filter(Class %in% c("HDAC5_COCA","HDAC5_SALINE"))

mat <- p_regressed[rownames(p)%in% top_labelled$Gene,colnames(p_regressed) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Condition = pd_filt$Condition) %>%
        pivot_longer(!Condition, names_to = "Gene", values_to="Exp")

pdf("dge/HDAC5c_HDAC5s/Boxplots_TopGenes_HDAC5c_HDAC5s.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Condition", "Exp", color = "Condition",
 palette = c("red", "black")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/HDAC5c_HDAC5s/Vulcano_Plot_HDAC5c_HDAC5s.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-1.5,+1.5)
dev.off()

# heatmap
mat <- p_regressed[rownames(p_regressed)%in% HDAC5c_HDAC5s_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Condition        <- c("red", "black")
names(Condition) <- c("COCA", "SALINE")
anno_colors <- list(Condition = Condition)
pdf("dge/HDAC5c_HDAC5s/Heatmap_HDAC5c_HDAC5s.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


# Input for viz Plot for HDAC5c vs eGFPc
dir.create("dge/HDAC5c_eGFPc")
df <- HDAC5c_eGFPc %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
           unite("Class", Genotype:Condition, na.rm = TRUE, remove = F) %>%
            filter(Class %in% c("HDAC5_COCA","EGFP_COCA"))

mat <- p_regressed[rownames(p)%in% top_labelled$Gene,colnames(p_regressed) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd_filt$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/HDAC5c_eGFPc/Boxplots_TopGenes_HDAC5c_eGFPc.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("black", "red")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/HDAC5c_eGFPc/Vulcano_Plot_HDAC5c_eGFPc.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-1.5,+1.5)
dev.off()

# heatmap
mat <- p_regressed[rownames(p_regressed)%in% HDAC5c_eGFPc_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Genotype        <- c("red", "black")
names(Genotype) <- c("HDAC5", "EGFP")
anno_colors <- list(Genotype = Genotype)
pdf("dge/HDAC5c_eGFPc/Heatmap_HDAC5c_eGFPc.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


# Input for viz Plot for HDAC5s vs eGFPs
dir.create("dge/HDAC5s_eGFPs")
df <- HDAC5s_eGFPs %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
           unite("Class", Genotype:Condition, na.rm = TRUE, remove = F) %>%
            filter(Class %in% c("HDAC5_SALINE","EGFP_SALINE"))

mat <- p_regressed[rownames(p)%in% top_labelled$Gene,colnames(p_regressed) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd_filt$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/HDAC5s_eGFPs/Boxplots_TopGenes_HDAC5s_eGFPs.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("black", "red")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/HDAC5s_eGFPs/Vulcano_Plot_HDAC5s_eGFPs.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-1.5,+1.5)
dev.off()

# heatmap
mat <- p_regressed[rownames(p_regressed)%in% HDAC5s_eGFPs_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Genotype        <- c("red", "black")
names(Genotype) <- c("HDAC5", "EGFP")
anno_colors <- list(Genotype = Genotype)
pdf("dge/HDAC5s_eGFPs/Heatmap_HDAC5s_eGFPs.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


# Input for viz Plot for eGFPc vs eGFPs
dir.create("dge/eGFPc_eGFPs")
df <- eGFPc_eGFPs %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
           unite("Class", Genotype:Condition, na.rm = TRUE, remove = F) %>%
            filter(Class %in% c("EGFP_COCA","EGFP_SALINE"))

mat <- p_regressed[rownames(p)%in% top_labelled$Gene,colnames(p_regressed) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Condition = pd_filt$Condition) %>%
        pivot_longer(!Condition, names_to = "Gene", values_to="Exp")

pdf("dge/eGFPc_eGFPs/Boxplots_TopGenes_eGFPc_eGFPs.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Condition", "Exp", color = "Condition",
 palette = c("red", "black")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/eGFPc_eGFPs/Vulcano_Plot_eGFPc_eGFPs.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-1.5,+1.5)
dev.off()

# heatmap
mat <- p_regressed[rownames(p_regressed)%in% eGFPc_eGFPs_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Condition        <- c("red", "black")
names(Condition) <- c("COCA", "SALINE")
anno_colors <- list(Condition = Condition)
pdf("dge/eGFPc_eGFPs/Heatmap_eGFPc_eGFPs.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


####################
# Cross Comparison #
####################
tmp1<- read.table("dge/HDAC5c_HDAC5s_Dge_MouseID.txt",header=T)
tmp2 <- read.table("dge/HDAC5c_eGFPc_Dge_MouseID.txt",header=T)
tmp3 <- read.table("dge/HDAC5s_eGFPs_Dge_MouseID.txt",header=T)
tmp4 <- read.table("dge/eGFPc_eGFPs_Dge_MouseID.txt",header=T)

tmp1$Class <- "HDAC5c_HDAC5s"
tmp2$Class <- "HDAC5c_eGFPc"
tmp3$Class <- "HDAC5s_eGFPs"
tmp4$Class <- "eGFPc_eGFPs"

tmp1 <- tmp1 %>%
      select(-Direction)

tmp2 <- tmp2 %>%
      select(-Direction)

tmp3 <- tmp3 %>%
      select(-Direction)

tmp4 <- tmp4 %>%
      select(-Direction)

tmp_upset <- rbind(tmp1,tmp2,tmp3,tmp4)

l <- split(as.character(tmp_upset$Gene),tmp_upset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("dge/Upset_Plot_Intersection.pdf", width = 6, height = 4)
upset(fromList(l),,nsets = 4, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(HDAC5c_HDAC5s = "#89E651", HDAC5c_eGFPc = "#DB8BE4",HDAC5s_eGFPs = "#DDA187", eGFPc_eGFPs="#DDD3DD"), 
    alpha = 0.5))))
dev.off()


# Boxplot all genes
p_regressed <- as.data.frame(p_regressed)
mat <-  p_regressed[rownames(p_regressed)%in% sign_posthoc_p$Gene,] %>% 
        rownames_to_column("Gene") %>%
        melt() %>%
        as.data.frame() 

tmp <- merge(mat,pd,by.x="variable",by.y="row.names") %>%
        mutate(Gene == as.factor(Gene)) %>%
        mutate(Condition = fct_relevel(Condition, c("COCA","SALINE")), Treatment = fct_relevel(Genotype, c("HDAC5","EGFP")))

dir.create("dge/Plot")

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- c("#454b87","#bd3106")
doPlot = function(sel_name) 
{
    df = subset(tmp, Gene == sel_name)
    PLOT= ggboxplot(df, "Condition", "value", color = "Treatment",
              palette = cols,
              outlier.shape = NA) +
      xlab("")+ 
      ylab("Gene Exp")+
        theme_classic() +
        rotate_x_text(angle = 45)

    print(PLOT)
    ggsave(sprintf("dge/Plot/%s.pdf", sel_name),width=4,height=4)
 }

lapply(unique(tmp$Gene), doPlot)



# ComplexHeatmap Fold Changes
load("dge/Sarah_Dge_Defined.RData")


HDAC5c_HDAC5s_Sign$Class <- "HDAC5c_HDAC5s"
HDAC5c_eGFPc_Sign$Class <- "HDAC5c_eGFPc"
HDAC5s_eGFPs_Sign$Class <- "HDAC5s_eGFPs"
eGFPc_eGFPs_Sign$Class <- "eGFPc_eGFPs"

dge <- rbind(HDAC5c_HDAC5s_Sign,HDAC5c_eGFPc_Sign,HDAC5s_eGFPs_Sign,eGFPc_eGFPs_Sign)


tmp1 <- dge %>% 
        group_by(Class) %>%
         mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
         mutate(Direction = case_when(logFC > 0 ~ "UpReg", logFC < -0 ~ "DownReg")) %>%
         as.data.frame()


tab<-table(tmp1$Class, tmp1$Direction)

mat<-matrix(nrow=length(unique(tmp1$Gene)),ncol=4)
rownames(mat)<-unique(tmp1$Gene)
colnames(mat)<-unique(tmp1$Class)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-tmp1[which(tmp1$Gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$Class==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$Gene),]
mat<-mat[order(noDup$Class, noDup$logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(Class=colnames(mat), 
      col= list(Class=c(
                           "HDAC5c_HDAC5s"="#89E651",
                           "HDAC5s_eGFPs"="#DB8BE4",
                           "HDAC5c_eGFPc"="#DDD3DD",
                           "eGFPc_eGFPs"="#DDA187"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,200)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-200,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)


pdf("dge/DGE_heatmap_foldCs.pdf", width=2.5, height=6)
Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(1,0,-1),c("red","white","blue")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
dev.off()



















