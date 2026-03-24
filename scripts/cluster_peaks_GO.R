
#GO FOR SC DATA RE-CALLED ON CLUSTERS


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







