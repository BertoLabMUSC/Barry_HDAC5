suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
load("dge/Sarah_Dge_Defined.RData")

# Convert to human 
human = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI_1 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = HDAC5c_HDAC5s_Sign$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
MGI_2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = HDAC5c_eGFPc_Sign$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
MGI_3 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = HDAC5s_eGFPs_Sign$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
MGI_4 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = eGFPc_eGFPs_Sign$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb_1 <- merge(HDAC5c_HDAC5s_Sign,MGI_1,by.x="Gene",by.y="MGI.symbol",all=F)
signComb_2 <- merge(HDAC5c_eGFPc_Sign,MGI_2,by.x="Gene",by.y="MGI.symbol",all=F)
signComb_3 <- merge(HDAC5s_eGFPs_Sign,MGI_3,by.x="Gene",by.y="MGI.symbol",all=F)
signComb_4 <- merge(eGFPc_eGFPs_Sign,MGI_4,by.x="Gene",by.y="MGI.symbol",all=F)

#################
# HDAC5c_HDAC5s #
#################

df_1 <- signComb_1 %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp_1 <- data.frame(Gene = df_1$Gene, Direction = rep("All",nrow(df_1)))

df_1 <- rbind(df_1,tmp_1) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df_1,"dge/HDAC5c_HDAC5s_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df_1 <- HDAC5c_HDAC5s_Sign %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp_1 <- data.frame(Gene = df_1$Gene, Direction = rep("All",nrow(df_1)))

df_1 <- rbind(df_1,tmp_1) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df_1,"dge/HDAC5c_HDAC5s_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)

#################
# HDAC5c_eGFPc #
#################

df_2 <- signComb_2 %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp_2 <- data.frame(Gene = df_2$Gene, Direction = rep("All",nrow(df_2)))

df_2 <- rbind(df_2,tmp_2) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df_2,"dge/HDAC5c_eGFPc_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df_2 <- HDAC5c_eGFPc_Sign %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp_2 <- data.frame(Gene = df_2$Gene, Direction = rep("All",nrow(df_2)))

df_2 <- rbind(df_2,tmp_2) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df_2,"dge/HDAC5c_eGFPc_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)


#################
# HDAC5s_eGFPs #
#################

df_3 <- signComb_3 %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp_3 <- data.frame(Gene = df_3$Gene, Direction = rep("All",nrow(df_3)))

df_3 <- rbind(df_3,tmp_3) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df_3,"dge/HDAC5s_eGFPs_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df_3 <- HDAC5s_eGFPs_Sign %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp_3 <- data.frame(Gene = df_3$Gene, Direction = rep("All",nrow(df_3)))

df_3 <- rbind(df_3,tmp_3) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df_3,"dge/HDAC5s_eGFPs_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)


#################
# eGFPc_eGFPs #
#################


df_4 <- signComb_4 %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp_4 <- data.frame(Gene = df_4$Gene, Direction = rep("All",nrow(df_4)))

df_4 <- rbind(df_4,tmp_4) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df_4,"dge/eGFPc_eGFPs_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df_4 <- eGFPc_eGFPs_Sign %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp_4 <- data.frame(Gene = df_4$Gene, Direction = rep("All",nrow(df_4)))

df_4 <- rbind(df_4,tmp_4) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df_4,"dge/eGFPc_eGFPs_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)
