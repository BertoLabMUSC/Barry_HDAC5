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
library(emmeans)
library(broom)
library(janitor)
})

dir.create("dge")


load("futcounts/Expression_Input.RData")



pd <- data.frame(row.names=colnames(exp), 
      Condition = as.factor(c(rep("SALINE",6),rep("COCA",6))),
      Genotype = as.factor(c("EGFP","HDAC5","EGFP","HDAC5","EGFP","HDAC5","HDAC5","EGFP","HDAC5","EGFP","HDAC5","EGFP")))
write.table(pd,"pData.txt",sep="\t",quote=F)

# Expression level
exp <- exp %>%
        rownames_to_column("Gene") %>%
        group_by(Gene) %>% 
        summarise_each(list(sum)) %>%
        column_to_rownames("Gene")


# Filter the expression by condition
filter=apply(cpm, 1, function(x) (all(x[1:12] >= 0.5)))
count_filt <- exp[filter,]
cpm_filt <- cpm[filter,]

logCPM <- log2(cpm_filt+1)
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)


pdf("dge/PCA.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("red","black"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()
dev.off()

# Surrogate variable
mod <- model.matrix(~., pd)
mod0 <- model.matrix(~ 1, pd)

svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100,method="two-step")

svaobj$sv <- data.frame(svaobj$sv)
colnames(svaobj$sv) = c(paste0('SV',seq(svaobj$n.sv)))
pdSv <- cbind(pd,svaobj$sv)



# Regression Expression
pd_sva <- pdSv %>%
            dplyr::select(-Genotype, -Condition) %>% #Removing ID, Region, Diagnosis from the model
            droplevels()

betas<-future_lapply(1:nrow(p), function(x)
            {
                lm(unlist(p[x,])~., data = pd_sva)
            })

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
p_regressed <- residuals+matrix(future_apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(p_regressed)<-rownames(p)

write.table(p_regressed,"dge/expression_regressed.txt",sep="\t",quote=F)

pdf("dge/PCA_AdjustedForConfound.pdf",width=8,height=8,useDingbats=FALSE)
pca.Sample<-prcomp(t(p_regressed))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("black","red"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()  +
theme(legend.position= "none")
dev.off()

# Limma
design <- model.matrix(~Genotype + Condition + Genotype*Condition+SV1,pdSv) # describe model to be fit
fitLM <- lmFit(p, design, method ="robust")

# Check only Condition
cont_1 <- contrasts.fit(fitLM, coef = 2) # Directly test second coefficient
fitEb_1 <- eBayes(cont_1,robust = TRUE)

# Check only Genotype
cont_2 <- contrasts.fit(fitLM, coef = 3) # Directly test second coefficient
fitEb_2 <- eBayes(cont_2,robust = TRUE)

# Check only interaction
cont_3 <- contrasts.fit(fitLM, coef = 5) # Directly test second coefficient
fitEb_3 <- eBayes(cont_3,robust = TRUE)

#
Geno_FullTab = topTable(fitEb_1, coef = "GenotypeHDAC5",number=nrow(p)) %>%
             rownames_to_column("Gene")

Geno_DGE <- Geno_FullTab %>%
        mutate(Abs = abs(logFC)) %>%
        filter(adj.P.Val < 0.05, Abs > 0.3) %>%
        arrange(desc(Abs))

Geno_DGE_relaxed <- Geno_FullTab %>%
        mutate(Abs = abs(logFC)) %>%
        filter(adj.P.Val < 0.05) %>%
        arrange(desc(Abs))

# 
Cond_FullTab = topTable(fitEb_2, coef = "ConditionSALINE",number=nrow(p)) %>%
             rownames_to_column("Gene")

Cond_DGE <- Cond_FullTab %>%
        mutate(Abs = abs(logFC)) %>%
        filter(adj.P.Val < 0.05, Abs > 0.3) %>%
        arrange(desc(Abs))

Cond_DGE_relaxed <- Cond_FullTab %>%
        mutate(Abs = abs(logFC)) %>%
        filter(adj.P.Val < 0.05) %>%
        arrange(desc(Abs))


#
Inter_FullTab = topTable(fitEb_3, coef = "GenotypeHDAC5:ConditionSALINE",number=nrow(p)) %>%
             rownames_to_column("Gene")

Inter_DGE <- Inter_FullTab %>%
        mutate(Abs = abs(logFC)) %>%
        filter(adj.P.Val < 0.05, Abs > 0.3) %>%
        arrange(desc(Abs))

save(Geno_FullTab,Cond_FullTab,Inter_FullTab,Geno_DGE,Cond_DGE,Inter_DGE,pdSv,p_regressed,  file = "dge/Sarah_Dge_Data.RData")


openxlsx::write.xlsx(Geno_FullTab, 
                     file = "dge/Sarah_Genotype_FullStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Cond_FullTab, 
                     file = "dge/Sarah_Condition_FullStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Inter_FullTab, 
                     file = "dge/Sarah_Intersection_FullStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


openxlsx::write.xlsx(Geno_DGE, 
                     file = "dge/Sarah_Genotype_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Geno_DGE_relaxed, 
                     file = "dge/Sarah_Genotype_DGE_relaxed.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Cond_DGE, 
                     file = "dge/Sarah_Condition_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Cond_DGE_relaxed, 
                     file = "dge/Sarah_Condition_DGE_relaxed.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(Inter_DGE, 
                     file = "dge/Sarah_Intersection_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

#######################################

# Emmeans post-hoc


# Modeling
model1 <- 'geneExpr ~ Genotype * Condition' # Null model

# Function for fitting the two models with/without the "diagnosis".
emmeans_lm_p <- function(vectorizedExpression) {
  tmpMetaData <- cbind(pdSv, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  pval <- emmeans::emmeans(residuals1,pairwise ~ Genotype*Condition,adjust="tukey")$contrasts %>%
          broom::tidy() %>%
          select(contrast,adj.p.value) %>% 
          pivot_wider(names_from = contrast, values_from = adj.p.value) %>%
          as.data.frame()
       
}

emmeans_lm_p_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_p(vectorizedExpression))
}


posthoc_p <- future_apply(p_regressed, 1, emmeans_lm_p_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_p) <- rownames(p_regressed)

# Eff Size
emmeans_lm_eff <- function(vectorizedExpression) {
  tmpMetaData <- cbind(pdSv, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  
  effect_size <- emmeans::emmeans(residuals1, adjust="tukey",pairwise ~ Genotype*Condition)$contrasts %>%
          broom::tidy() %>%
          select(contrast,estimate) %>% 
          pivot_wider(names_from = contrast, values_from = estimate) %>%
          as.data.frame()
       
}

emmeans_lm_eff_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_eff(vectorizedExpression))
}


posthoc_eff <- future_apply(p_regressed, 1, emmeans_lm_eff_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_eff) <- rownames(p_regressed)


write.table(posthoc_p,"dge/Sarah_Dge_Interaction_Posthoc_Pvalues.txt",sep="\t",quote=F)
write.table(posthoc_eff,"dge/Sarah_Dge_Interaction_Posthoc_EffSize.txt",sep="\t",quote=F)



# Filter for significant
sign_posthoc_p <- posthoc_p %>%
                   filter(if_any(is.numeric, ~ .x < 0.05)) %>%
                   rownames_to_column("Gene")

sign_posthoc_eff <- posthoc_eff[rownames(posthoc_eff) %in% sign_posthoc_p$Gene,] %>%
                   rownames_to_column("Gene")


# Save
save(posthoc_eff,posthoc_p,sign_posthoc_p,sign_posthoc_eff,  file = "dge/Sarah_Dge_Interaction_Posthoc.RData")


openxlsx::write.xlsx(sign_posthoc_p, 
                     file = "dge/Sarah_Dge_Interaction_Posthoc_Pvalues_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


openxlsx::write.xlsx(sign_posthoc_eff, 
                     file = "dge/Sarah_Dge_Interaction_Posthoc_EffSize_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")



# Combine

HDAC5c_HDAC5s <- data.frame(Gene = rownames(posthoc_p), logFC = posthoc_eff$hdac5_coca_hdac5_saline, FDR = posthoc_p$hdac5_coca_hdac5_saline)

HDAC5c_eGFPc <- data.frame(Gene = rownames(posthoc_p), logFC = -1 * posthoc_eff$egfp_coca_hdac5_coca, FDR = posthoc_p$egfp_coca_hdac5_coca)

HDAC5s_eGFPs <- data.frame(Gene = rownames(posthoc_p), logFC = -1 * posthoc_eff$egfp_saline_hdac5_saline, FDR = posthoc_p$egfp_saline_hdac5_saline)

eGFPc_eGFPs <- data.frame(Gene = rownames(posthoc_p), logFC = posthoc_eff$egfp_coca_egfp_saline, FDR = posthoc_p$egfp_coca_egfp_saline)




HDAC5c_HDAC5s_Sign <- HDAC5c_HDAC5s %>% filter(FDR < 0.05, abs(logFC) > 0.3)
HDAC5c_eGFPc_Sign <- HDAC5c_eGFPc %>% filter(FDR < 0.05, abs(logFC) > 0.3)
HDAC5s_eGFPs_Sign <- HDAC5s_eGFPs %>% filter(FDR < 0.05, abs(logFC) > 0.3)
eGFPc_eGFPs_Sign <- eGFPc_eGFPs %>% filter(FDR < 0.05, abs(logFC) > 0.3)

save(HDAC5c_HDAC5s,HDAC5c_eGFPc,HDAC5s_eGFPs,eGFPc_eGFPs,HDAC5c_HDAC5s_Sign,HDAC5c_eGFPc_Sign,HDAC5s_eGFPs_Sign,eGFPc_eGFPs_Sign, file = "dge/Sarah_Dge_Defined.RData")


# Save input data
save(p,p_regressed,pd, file = "dge/Expression_Object.RData")
