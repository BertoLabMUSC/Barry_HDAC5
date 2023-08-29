
suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(enrichR)
})


#################
# HDAC5c_HDAC5s #
#################
dir.create("dge/HDAC5c_HDAC5s/functional_enrichment/")
dge <- read.table("dge/HDAC5c_HDAC5s_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

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
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/HDAC5c_HDAC5s/functional_enrichment/%s_HDAC5c_HDAC5s_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/HDAC5c_HDAC5s/functional_enrichment/%s_HDAC5c_HDAC5s_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/HDAC5c_HDAC5s/functional_enrichment/%s_HDAC5c_HDAC5s_barplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/HDAC5c_HDAC5s_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HDAC5c_HDAC5s'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/HDAC5c_HDAC5s/functional_enrichment/ENRICHR_HDAC5c_HDAC5s_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/HDAC5c_HDAC5s/functional_enrichment/ENRICHR_HDAC5c_HDAC5s_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()


#################
# HDAC5c_eGFPc #
#################
dir.create("dge/HDAC5c_eGFPc/functional_enrichment/")

dge <- read.table("dge/HDAC5c_eGFPc_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

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
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/HDAC5c_eGFPc/functional_enrichment/%s_HDAC5c_eGFPc_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/HDAC5c_eGFPc/functional_enrichment/%s_HDAC5c_eGFPc_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/HDAC5c_eGFPc/functional_enrichment/%s_HDAC5c_eGFPc_barplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/HDAC5c_eGFPc_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HDAC5c_eGFPc'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/HDAC5c_eGFPc/functional_enrichment/ENRICHR_HDAC5c_eGFPc_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/HDAC5c_eGFPc/functional_enrichment/ENRICHR_HDAC5c_eGFPc_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()


#################
# HDAC5s_eGFPs #
#################
dir.create("dge/HDAC5s_eGFPs/functional_enrichment/")

dge <- read.table("dge/HDAC5s_eGFPs_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

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
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/HDAC5s_eGFPs/functional_enrichment/%s_HDAC5s_eGFPs_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/HDAC5s_eGFPs/functional_enrichment/%s_HDAC5s_eGFPs_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/HDAC5s_eGFPs/functional_enrichment/%s_HDAC5s_eGFPs_barplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/HDAC5s_eGFPs_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HDAC5s_eGFPs'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/HDAC5s_eGFPs/functional_enrichment/ENRICHR_HDAC5s_eGFPs_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/HDAC5s_eGFPs/functional_enrichment/ENRICHR_HDAC5s_eGFPs_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()


#################
# eGFPc_eGFPs #
#################
dir.create("dge/eGFPc_eGFPs/functional_enrichment/")

dge <- read.table("dge/eGFPc_eGFPs_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

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
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/eGFPc_eGFPs/functional_enrichment/%s_eGFPc_eGFPs_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/eGFPc_eGFPs/functional_enrichment/%s_eGFPc_eGFPs_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/eGFPc_eGFPs/functional_enrichment/%s_eGFPc_eGFPs_barplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/eGFPc_eGFPs_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'eGFPc_eGFPs'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/eGFPc_eGFPs/functional_enrichment/ENRICHR_eGFPc_eGFPs_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/eGFPc_eGFPs/functional_enrichment/ENRICHR_eGFPc_eGFPs_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()


#####################
# Gene Onto for All #
#####################

tmp1<- read.table("dge/HDAC5c_HDAC5s_Dge_HumanID.txt",header=T)
tmp2 <- read.table("dge/HDAC5c_eGFPc_Dge_HumanID.txt",header=T)
tmp3 <- read.table("dge/HDAC5s_eGFPs_Dge_HumanID.txt",header=T)
tmp4 <- read.table("dge/eGFPc_eGFPs_Dge_HumanID.txt",header=T)

tmp1$Class <- "HDAC5c_HDAC5s"
tmp2$Class <- "HDAC5c_eGFPc"
tmp3$Class <- "HDAC5s_eGFPs"
tmp4$Class <- "eGFPc_eGFPs"

tmp1 <- tmp1 %>%
      select(-Direction) %>%
      distinct()

tmp2 <- tmp2 %>%
      select(-Direction) %>%
      distinct()

tmp3 <- tmp3 %>%
      select(-Direction) %>%
      distinct()

tmp4 <- tmp4 %>%
      select(-Direction) %>%
      distinct()

tmp_go <- rbind(tmp1,tmp2,tmp3,tmp4)


dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(tmp_go$Class))){
  print(cur)
  # select genes
  cur_genes <- tmp_go %>%
    subset(Class == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/ENRICHR_Together_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/ENRICHR_Together_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()







































