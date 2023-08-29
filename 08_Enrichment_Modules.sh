# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19982, brain expressed = 15585)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir wgcna_output/enrichments_wgcna/

cp utils/geneset/*.RData wgcna_output/enrichments_wgcna/
cp utils/Enrichment_wgcna.r wgcna_output/enrichments_wgcna/
cp wgcna_output/ModuleOutput.txt wgcna_output/enrichments_wgcna/
cp wgcna_output/ModuleOutput_HumanID.txt wgcna_output/enrichments_wgcna/

cd wgcna_output/enrichments_wgcna/

mkdir STATS/


# scRNA healty
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_Markers -W 8 -H 4
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet_Simple.RData -p -b 15585 -o STATS/Allen_Simple_Markers -W 8 -H 3
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_scMouse.RData -p -b 15585 -o STATS/scMouse_DropViz -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_Allen_Mouse_PFC.RData -p -b 15585 -o STATS/scMouse_Allen_PFC -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_Cowan_PFC.RData -p -b 15585 -o STATS/scMouse_Cowan_PFC -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_Cowan_PFC_Simple.RData -p -b 15585 -o STATS/scMouse_Cowan_PFC_Simple -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_AdultMouse.RData -p -b 15585 -o STATS/scMouse_AdultMouse_Brain -W 8 -H 6
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_YoungMouse.RData -p -b 15585 -o STATS/scMouse_YoungMouse_Brain -W 8 -H 6
Rscript Enrichment_wgcna.r -g ModuleOutput.txt -l GeneSets_OligoMouse.RData -p -b 15585 -o STATS/scMouse_OligoMouse_Brain -W 8 -H 4


# scRNA Disorders
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 8 -H 5

# Neuropsy
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 8 -H 3
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 8 -H 5
Rscript Enrichment_wgcna.r -g ModuleOutput_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 8 -H 2

rm *.RData

