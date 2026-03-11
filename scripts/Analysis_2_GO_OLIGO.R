.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
.libPaths()

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_pseudobulked_granges.RData")

#genic annotations using chipseeker

#assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#issues with grnages objects 
#chnages chromosome naming style of grnages to USSC FORMAT
library(GenomeInfoDb)
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"
#annotation
#making a list of all samples (tutorial does it this way, helpful for visualization later on)

#create list of 4 grnages objects I made before
peak_list <- list(Bulk = bulk_gr, Opalin = opalin_gr, Plekhg1 = plekhg1_gr, OPCs = opc_gr)

#apply annotatepeak to each grnages
peakAnnoList <- lapply(peak_list, annotatePeak, TxDb = txdb,
                       tssRegion = c(-3000, 3000), verbose = FALSE)


peakAnnoList
plotAnnoBar(peakAnnoList)

#can also look at to look at the distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")


#create annotated files for reference later on 
# extract annotation data frames
bulk_anno_df    <- as.data.frame(peakAnnoList$Bulk)
opalin_anno_df  <- as.data.frame(peakAnnoList$Opalin)
plekhg1_anno_df <- as.data.frame(peakAnnoList$Plekhg1)
opc_anno_df     <- as.data.frame(peakAnnoList$OPCs)

library(dplyr)
# function to add gene symbols to annotation data frame
#add gene syblos to annotation data frame (chipseeker annoates with entrezIDs which are only numbers)
add_gene_symbols <- function(anno_df) {
  entrez <- anno_df$geneId
  
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys    = entrez,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  
  # keep only first gene name per EntrezID to avoid duplicate rows
  annotations_edb <- annotations_edb[!duplicated(annotations_edb$ENTREZID), ]
  
  anno_df %>% left_join(annotations_edb, by = c("geneId" = "ENTREZID"))
}

# apply to all four datasets
bulk_anno_df    <- add_gene_symbols(bulk_anno_df)
opalin_anno_df  <- add_gene_symbols(opalin_anno_df)
plekhg1_anno_df <- add_gene_symbols(plekhg1_anno_df)
opc_anno_df     <- add_gene_symbols(opc_anno_df)


# write to file
write.csv(bulk_anno_df,    "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/bulk_annotation.csv",    row.names = FALSE)
write.csv(opalin_anno_df,  "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/opalin_annotation.csv",  row.names = FALSE)
write.csv(plekhg1_anno_df, "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/plekhg1_annotation.csv", row.names = FALSE)
write.csv(opc_anno_df,     "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/opc_annotation.csv",     row.names = FALSE)



#check if bulk.csv has gene symbols
bulk_csv <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/bulk_annotation.csv")


#Functional enrichment (GO) using clusterprofile 
head(bulk_anno_df[, c("seqnames", "start", "end", "geneId", "GENENAME")])
library(clusterProfiler)
library(ggplot2)

# recreate overlaps
overlaps_opalin  <- findOverlaps(bulk_gr, opalin_gr)
overlaps_plekhg1 <- findOverlaps(bulk_gr, plekhg1_gr)
overlaps_opc     <- findOverlaps(bulk_gr, opc_gr)

# classify peaks
opalin_shared      <- bulk_peaks[unique(queryHits(overlaps_opalin))]
opalin_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opalin)]
opalin_unique_sc   <- opalin_peaks[!seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin)]

plekhg1_shared      <- bulk_peaks[unique(queryHits(overlaps_plekhg1))]
plekhg1_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_plekhg1)]
plekhg1_unique_sc   <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1)]

opc_shared      <- bulk_peaks[unique(queryHits(overlaps_opc))]
opc_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opc)]
opc_unique_sc   <- opc_peaks[!seq_along(opc_peaks) %in% subjectHits(overlaps_opc)]

# get bulk unique peak names
bulk_unique_df <- bulk_anno_df[paste(bulk_anno_df$seqnames, 
                                     bulk_anno_df$start, 
                                     bulk_anno_df$end, sep = "_") %in% opalin_unique_bulk, ]

# extract entrez IDs, remove NAs
bulk_unique_entrez <- unique(bulk_unique_df$geneId)
bulk_unique_entrez <- bulk_unique_entrez[!is.na(bulk_unique_entrez)]

# run GO enrichment
ego_bulk_unique <- enrichGO(gene          = bulk_unique_entrez,
                            keyType       = "ENTREZID",
                            OrgDb         = org.Hs.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)


# view results
head(ego_bulk_unique)

# extract results to data frame and write to file
bulk_unique_go <- data.frame(ego_bulk_unique)

write.csv(bulk_unique_go, "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/GO_oligo_bulk_unique.csv", row.names = FALSE)
library(ggplot2)
#dotplot, looking at top 50 sig genes 

dotplot(ego_bulk_unique, showCategory = 50)

p <- dotplot(ego_bulk_unique, showCategory = 50) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7)) +
  scale_size_continuous(range = c(1, 6))

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_GO_oligo_bulk_unique.png",
       plot = p, width = 8, height = 18, dpi = 150)

#KEGG analysis to compare. (this is another type of pathway level annotating)
ekegg_bulk_unique <- enrichKEGG(gene         = bulk_unique_entrez,
                                organism     = "hsa",
                                pvalueCutoff = 0.05)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/kegg_oligo_bulk_unique.png",
       plot = dotplot(ekegg_bulk_unique, showCategory = 30),
       width = 10, height = 12, dpi = 150)

write.csv(data.frame(ekegg_bulk_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/KEGG_oligo_bulk_unique.csv",
          row.names = FALSE)

# Repeat analysis for each of the subtype unique peaks
#Opalin+ unique peaks

# get Opalin+ unique peak names
opalin_unique_df <- opalin_anno_df[paste(gsub("chr", "", opalin_anno_df$seqnames),
                                         opalin_anno_df$start,
                                         opalin_anno_df$end, sep = "-") %in% opalin_unique_sc, ]

# extract entrez IDs, remove NAs
opalin_unique_entrez <- unique(opalin_unique_df$geneId)
opalin_unique_entrez <- opalin_unique_entrez[!is.na(opalin_unique_entrez)]

# run GO enrichment
ego_opalin_unique <- enrichGO(gene          = opalin_unique_entrez,
                              keyType       = "ENTREZID",
                              OrgDb         = org.Hs.eg.db,
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

# save and plot
write.csv(data.frame(ego_opalin_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/GO_oligo_opalin_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_GO_oligo_opalin_unique.png",
       plot = dotplot(ego_opalin_unique, showCategory = 30) +
         labs(title ="GO enrichment, Opalin+ unique peaks"),
       width = 10, height = 14, dpi = 150)


#KEGG

ekegg_opalin_unique <- enrichKEGG(gene         = opalin_unique_entrez,
                                  organism     = "hsa",
                                  pvalueCutoff = 0.05)

# save and plot
write.csv(data.frame(ekegg_opalin_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/KEGG_oligo_opalin_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_KEGG_oligo_opalin_unique.png",
       plot = dotplot(ekegg_opalin_unique, showCategory = 30) +
         labs(title = "KEGG Enrichment - Opalin+ unique peaks"),
       width = 10, height = 14, dpi = 150)

#GO for Plekhg1+
plekhg1_unique_df <- plekhg1_anno_df[paste(gsub("chr", "", plekhg1_anno_df$seqnames),
                                           plekhg1_anno_df$start,
                                           plekhg1_anno_df$end, sep = "-") %in% plekhg1_unique_sc, ]

plekhg1_unique_entrez <- unique(plekhg1_unique_df$geneId)
plekhg1_unique_entrez <- plekhg1_unique_entrez[!is.na(plekhg1_unique_entrez)]

ego_plekhg1_unique <- enrichGO(gene          = plekhg1_unique_entrez,
                               keyType       = "ENTREZID",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)

write.csv(data.frame(ego_plekhg1_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/GO_oligo_plekhg1_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_GO_oligo_plekhg1_unique.png",
       plot = dotplot(ego_plekhg1_unique, showCategory = 30) +
         labs(title = "GO Enrichment - Plekhg1+ unique peaks"),
       width = 10, height = 14, dpi = 150)

#KEGG

ekegg_plekhg1_unique <- enrichKEGG(gene         = plekhg1_unique_entrez,
                                   organism     = "hsa",
                                   pvalueCutoff = 0.05)

write.csv(data.frame(ekegg_plekhg1_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/KEGG_oligo_plekhg1_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_KEGG_oligo_plekhg1_unique.png",
       plot = dotplot(ekegg_plekhg1_unique, showCategory = 30) +
         labs(title = "KEGG Enrichment - Plekhg1+ unique peaks"),
       width = 10, height = 14, dpi = 150)


#OPCs GO

opc_unique_df <- opc_anno_df[paste(gsub("chr", "", opc_anno_df$seqnames),
                                   opc_anno_df$start,
                                   opc_anno_df$end, sep = "-") %in% opc_unique_sc, ]

opc_unique_entrez <- unique(opc_unique_df$geneId)
opc_unique_entrez <- opc_unique_entrez[!is.na(opc_unique_entrez)]

ego_opc_unique <- enrichGO(gene          = opc_unique_entrez,
                           keyType       = "ENTREZID",
                           OrgDb         = org.Hs.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)

write.csv(data.frame(ego_opc_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/GO_oligo_opc_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_GO_oligo_opc_unique.png",
       plot = dotplot(ego_opc_unique, showCategory = 30) +
         labs(title = "GO Enrichment - OPC unique peaks"),
       width = 10, height = 14, dpi = 150)

#KEGG
ekegg_opc_unique <- enrichKEGG(gene         = opc_unique_entrez,
                               organism     = "hsa",
                               pvalueCutoff = 0.05)

write.csv(data.frame(ekegg_opc_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/KEGG_oligo_opc_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_KEGG_oligo_opc_unique.png",
       plot = dotplot(ekegg_opc_unique, showCategory = 30) +
         labs(title = "KEGG Enrichment - OPC unique peaks"),
       width = 10, height = 14, dpi = 150)


#try combining the sc subsets

# combine all SC unique entrez IDs
all_sc_unique_entrez <- unique(c(opalin_unique_entrez, plekhg1_unique_entrez, opc_unique_entrez))
all_sc_unique_entrez <- all_sc_unique_entrez[!is.na(all_sc_unique_entrez)]

# GO enrichment
ego_all_sc_unique <- enrichGO(gene          = all_sc_unique_entrez,
                              keyType       = "ENTREZID",
                              OrgDb         = org.Hs.eg.db,
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

# KEGG enrichment
ekegg_all_sc_unique <- enrichKEGG(gene         = all_sc_unique_entrez,
                                  organism     = "hsa",
                                  pvalueCutoff = 0.05)

# save and plot
write.csv(data.frame(ego_all_sc_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/GO_oligo_all_sc_unique.csv",
          row.names = FALSE)

write.csv(data.frame(ekegg_all_sc_unique),
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/KEGG_oligo_all_sc_unique.csv",
          row.names = FALSE)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_GO_oligo_all_sc_unique.png",
       plot = dotplot(ego_all_sc_unique, showCategory = 30) +
         labs(title = "GO Enrichment - All SC unique peaks"),
       width = 10, height = 14, dpi = 150)

ggsave("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/dotplot_KEGG_oligo_all_sc_unique.png",
       plot = dotplot(ekegg_all_sc_unique, showCategory = 30) +
         labs(title = "KEGG Enrichment - All SC unique peaks"),
       width = 10, height = 14, dpi = 150)


