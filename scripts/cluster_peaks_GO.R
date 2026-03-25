
#GO FOR SC DATA RE-CALLED ON CLUSTERS
#chip seeker annoates each peak with genomic location and nearest gene 
#annoate peak sssigns each peak a genic location and add_gene_symbols adds
# human readable gene symbols as TxDb only gives EntrezIDs. Uses ensDb to add gene names

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(dplyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# fix chromosome naming
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

# annotate
opalin_anno_cp  <- as.data.frame(annotatePeak(opalin_gr,  tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))
plekhg1_anno_cp <- as.data.frame(annotatePeak(plekhg1_gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))
opc_anno_cp     <- as.data.frame(annotatePeak(opc_gr,     tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))

# add gene symbols
add_gene_symbols <- function(anno_df) {
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys    = anno_df$geneId,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  annotations_edb <- annotations_edb[!duplicated(annotations_edb$ENTREZID), ]
  anno_df %>% left_join(annotations_edb, by = c("geneId" = "ENTREZID"))
}

opalin_anno_cp  <- add_gene_symbols(opalin_anno_cp)
plekhg1_anno_cp <- add_gene_symbols(plekhg1_anno_cp)
opc_anno_cp     <- add_gene_symbols(opc_anno_cp)


head(opalin_anno_cp[, c("seqnames", "start", "end", "annotation", "GENENAME")])

colnames(opalin_anno_cp)
head(opalin_anno_cp[, c("seqnames", "start", "end", "annotation", "GENENAME.x")])


#GO ENRICHMENT 
library(clusterProfiler)
library(ggplot2)

out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out"

# function to get entrez IDs from unique SC peaks
get_entrez_sc <- function(anno_df, unique_peaks) {
  key <- paste(gsub("chr", "", anno_df$seqnames), anno_df$start, anno_df$end, sep = "-")
  entrez <- unique(anno_df$geneId[key %in% unique_peaks])
  entrez[!is.na(entrez)]
}

# extract entrez IDs for SC unique peaks
opalin_unique_entrez  <- get_entrez_sc(opalin_anno_cp,  opalin_unique_sc)
plekhg1_unique_entrez <- get_entrez_sc(plekhg1_anno_cp, plekhg1_unique_sc)
opc_unique_entrez     <- get_entrez_sc(opc_anno_cp,     opc_unique_sc)

cat("Opalin+ unique genes:", length(opalin_unique_entrez), "\n")
cat("Plekhg1+ unique genes:", length(plekhg1_unique_entrez), "\n")
cat("OPC unique genes:", length(opc_unique_entrez), "\n")


# run GO enrichment function
run_go <- function(entrez, label) {
  ego <- enrichGO(gene          = entrez,
                  keyType       = "ENTREZID",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  write.csv(data.frame(ego),
            file.path(out_dir, paste0("GO_clusterpeaks_", label, "_sc_unique.csv")),
            row.names = FALSE)
  
  ggsave(file.path(out_dir, paste0("dotplot_GO_clusterpeaks_", label, "_sc_unique.png")),
         plot = dotplot(ego, showCategory = 30) +
           labs(title = paste("GO -", label, "SC unique peaks")),
         width = 10, height = 14, dpi = 150)
  
  ego
}

ego_opalin  <- run_go(opalin_unique_entrez,  "opalin")
ego_plekhg1 <- run_go(plekhg1_unique_entrez, "plekhg1")
ego_opc     <- run_go(opc_unique_entrez,     "opc")

