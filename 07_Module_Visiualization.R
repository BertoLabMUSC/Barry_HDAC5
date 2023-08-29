# WGCNA
suppressPackageStartupMessages({
library(tidyverse)
library(WGCNA)
library(cluster)
library(ggplot2)
library(ggpubr)
library(igraph)
library(reshape2)
library(clusterProfiler)
library(enrichR)
library(cowplot)
})
enableWGCNAThreads()


# Module Viz
dir.create("wgcna_output/modules_visualizations/")


# Load tables 
load("dge/Sarah_Dge_Data.RData")

## Set up input file for WGCNA
pd <- pdSv %>%
      select(Condition,Genotype) %>%
      unite("Sample", Genotype:Condition, na.rm = TRUE, remove = FALSE) %>%
      arrange(Sample)


exp <- p_regressed[,match(rownames(pd),colnames(p_regressed))]
p_regressed <- p_regressed[,match(rownames(pd),colnames(p_regressed))]

mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)



# Hubs 
hubs <- mod %>% 
            group_by(ModuleColor) %>% 
            top_n(n = 25, wt = kWithin) %>%
            as.data.frame()

hub_list <- split(hubs, hubs$ModuleColor)


tmp <- list()
adjMat <- list()
modColors <- list()
valueList <- list()
colorList <- list()
g1 <- list()
g1a <- list()
g1b <- list()
layoutFR <- list()

for(i in 1:length(hub_list)){
tmp[[i]] <- exp[rownames(exp) %in% hub_list[[i]]$Gene,]
adjMat[[i]] <- bicor(t(tmp[[i]]))
adjMat[[i]][abs(adjMat[[i]])<0.5] <- 0
adjMat[[i]] <- abs(adjMat[[i]])
adjMat[[i]] <- adjMat[[i]][match(hub_list[[i]]$Gene,rownames(adjMat[[i]])), match(hub_list[[i]]$Gene,colnames(adjMat[[i]]))]
modColors[[i]] <- as.matrix(t(hub_list[[i]]$ModuleColor))
valueList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) as.numeric(!is.na(modColors[[i]][,x])))
colorList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) modColors[[i]][,x])

g1[[i]] <- graph.adjacency(as.matrix(adjMat[[i]]),mode="undirected",weighted=TRUE,diag=FALSE)
#g1[[i]] <- delete.vertices(g1[[i]],which(degree(g1[[i]])<1))

g1a[[i]] <- graph.adjacency(as.matrix(adjMat[[i]][1:10, 1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
#g1a[[i]] <- delete.vertices(g1a[[i]],which(degree(g1a[[i]])<1))
g1b[[i]] <- graph.adjacency(as.matrix(adjMat[[i]][11:25, 11:25]),mode="undirected",weighted=TRUE,diag=FALSE)
#g1b[[i]] <- delete.vertices(g1b[[i]],which(degree(g1b[[i]])<1))

layoutFR[[i]] <- rbind(layout.circle(g1a[[i]])/2, layout.circle(g1b[[i]]))
#layoutFR[[i]] <- layout_with_lgl(g1[[i]],maxiter = 500)

pdf(paste("wgcna_output/modules_visualizations/", names(hub_list)[[i]], "_igraph.pdf",sep = ""),width=6,height=6,useDingbats=FALSE)
plot(g1[[i]], edge.color=adjustcolor("grey", alpha.f = .3),
            edge.alpha = 0.25, vertex.color = as.character(hub_list[[i]]$ModuleColor), 
             vertex.label = as.character(hub_list[[i]]$Gene), vertex.label.dist = 1.1, 
             vertex.label.degree = -pi/4, vertex.label.color = "black", 
             vertex.label.family = "Helvetica", vertex.label.font = 3, 
             vertex.label.cex = 1, vertex.frame.color = "black", 
             layout = jitter(layoutFR[[i]]), vertex.size = 6)
             
dev.off()
}




# GGnet
dir.create("wgcna_output/MODULES_CYTOSCAPE")

files = list.files(path="wgcna_output/Cyto",pattern = 'CytoEdge',full.names =TRUE)
samples <- gsub("wgcna_output/Cyto/CytoEdge", "", files )
samples <- gsub(".txt", "", samples )
myfiles = lapply(files, read.table,sep="\t",fill=TRUE,header=T)
names(myfiles) <- samples

list_mod <- within(myfiles, rm(grey)) 

list_mod_top <- list()
for(i in 1:length(list_mod)){
pdf(paste("wgcna_output/MODULES_CYTOSCAPE/", names(list_mod)[[i]], "_cytoscape.pdf",sep = ""),width=8,height=8,useDingbats=FALSE)
list_mod_top[[i]] <- list_mod[[i]] %>% 
                        top_n(n = 300, wt = weight) %>%
                        as.data.frame() %>%
                        ggnet::ggnet2(label=TRUE, size = "degree", size.cut = 20,alpha = 0.5, color=names(list_mod)[[i]]) +
                        guides(color = FALSE, size = FALSE)
print(list_mod_top[[i]]) 
dev.off()
}



# Traits
p <- read.table("wgcna_output/modTraitP.txt")
rownames(p) <- gsub("ME","",rownames(p))
cor <- read.table("wgcna_output/modTraitCor.txt")
rownames(cor) <- gsub("ME","",rownames(cor))

p_filt <- p %>% 
      rownames_to_column("Module") %>%
      dplyr::select(Module, SampleEGFP_COCA, SampleEGFP_SALINE, SampleHDAC5_COCA, SampleHDAC5_SALINE) %>%
      filter(SampleEGFP_COCA < 0.05 | SampleEGFP_SALINE < 0.05 | SampleHDAC5_COCA < 0.05 | SampleHDAC5_SALINE < 0.05)

cor_filt <- cor[rownames(cor) %in% p_filt$Module,colnames(cor)%in%colnames(p_filt)] %>% 
      rownames_to_column("Module")

p <- melt(p_filt)
cor <- melt(cor_filt)
p$cor <- -1*cor$value
p$log <- -log10(p$value)
p$Direction <- ifelse(p$cor < 0,"Pos","Neg") %>% as.factor()
#p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))
#p$variable <- gsub("_LR","",p$variable)

p$log[p$log < 1.3] <- NA



p$abs<- ifelse(is.na(p$log), p$log, p$cor) %>% abs()


pdf("wgcna_output/Module_Traits_Bubble.pdf",width=3,height=6)
ggscatter(
        p, 
                  x = "variable",
                  y = "Module",
                  size="abs",
                  color="Direction",
                  palette=c("blue","red"),
                  alpha = 0.8,
                  xlab = "",
            ylab = "",) +
                  theme_minimal() +
                  rotate_x_text(angle = 45)
dev.off()


# Module Eigengene
tab=read.table("wgcna_output/Matrix_module_correlation.txt")
df=melt(tab)
df$variable=gsub("ME","",df$variable)

coolmod=p_filt$Module
df=df[df$variable %in% coolmod,]

pdf("wgcna_output/Module_Eigengene_CoolMod.pdf",width=8,height=6,useDingbats=FALSE)
ggplot(df, aes(x=Sample, y=value, fill=Sample,group=Sample)) +
geom_pointrange(mapping = aes(x = Sample, y = value,group=Sample,colour=Sample), position=position_dodge(width=0.5),
stat = "summary",
fun.ymin = function(z) {quantile(z,0.25)},
fun.ymax = function(z) {quantile(z,0.75)},
fun.y = median)+
theme_classic()+
geom_hline(yintercept = 0,linetype="dashed")+
scale_colour_manual(values=c("blue", "orange","purple","green"))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(title="Module EigenGene",x="", y = "Eigengene")+
ylim(-0.5,+0.5) +
facet_wrap(.~variable)
dev.off()



# Gene Onto
dir.create("wgcna_output/functional_enrichment/")

p <- read.table("wgcna_output/modTraitP.txt")
rownames(p) <- gsub("ME","",rownames(p))
cor <- read.table("wgcna_output/modTraitCor.txt")
rownames(cor) <- gsub("ME","",rownames(cor))

p_filt <- p %>% 
      rownames_to_column("Module") %>%
      dplyr::select(Module, SampleEGFP_COCA, SampleEGFP_SALINE, SampleHDAC5_COCA, SampleHDAC5_SALINE) %>%
      filter(SampleEGFP_COCA < 0.05 | SampleEGFP_SALINE < 0.05 | SampleHDAC5_COCA < 0.05 | SampleHDAC5_SALINE < 0.05)


cor_filt <- cor[rownames(cor) %in% p_filt$Module,colnames(cor)%in%colnames(p_filt)] 

mod <- read.table("wgcna_output/ModuleOutput_HumanID.txt",header=T)
mod <- mod %>% mutate(ModuleColor = as.factor(ModuleColor)) %>% droplevels()


l <- split(mod, mod$ModuleColor)

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("wgcna_output/functional_enrichment/%s_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

}


cur_result <- list()
for(i in 1:length(l)){
    cur_result[[i]] <- GeneOnto[[i]] %>% as.data.frame() %>% remove_rownames()
    cur_result[[i]]$Module <- as.character(unique(mod$ModuleColor))[[i]]
    collapsed_output <- bind_rows(cur_result) %>% as.data.frame()
}


collapsed_output %>%
  write.csv(file='wgcna_output/functional_enrichment/ENRICHR_Modules_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(Module %in% rownames(cor_filt)) %>% 
    mutate(log = -log10(pvalue)) %>%
    group_by(Module) %>%
    top_n(2,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Description, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), Module = as.factor(Module))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("wgcna_output/functional_enrichment/ENRICHR_HDAC5c_SpecificMod_GO_bubblechart.pdf", width = 12, height = 5)
ggballoonplot(input_bub, x = "Module", y = "Term2",
              size = "Count", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE) + coord_flip()

  dev.off()








