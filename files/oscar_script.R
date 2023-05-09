if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#my deg result volcano plot
temp_volc <- uni_genes
for (i in 1:nrow(uni_genes)) {
  temp_volc$X[i] <- gsub("\\..*","",temp_volc$X[i])
}
temp_volc$X[1] <- "Cdc5l"
temp_volc$X[2] <- "Zfp937"
temp_volc$X[3] <- "Fth1"
temp_volc$X[4] <- "Nubp2"
temp_volc$X[5] <- "3300002I08Rik"
temp_volc$X[6] <- "Zfp953"
temp_volc$X[7] <- "Selenow"
temp_volc$X[8] <- "Prdx2"
temp_volc$X[9] <- "Cpox"
temp_volc$X[10] <- "Gm42980"
temp_volc$X[11] <- "Gemin7"
temp_volc$X[12] <- "Cenps"
temp_volc$X[16] <- "Gm30400"
temp_volc$X[57] <- "Gm43674"
temp_volc$X[134] <- "Gm9312"
temp_volc$X[22] <- "Psg20"
temp_volc$X[36] <- "Plac8"
temp_volc$X[70] <- "Rab4a"
temp_volc$X[144] <- "Elovl1"
temp_volc$X[49] <- "Ddx4"

EnhancedVolcano(temp_volc,
                lab = temp_volc$X,
                x = 'log2FoldChange',
                y = 'pvalue')
#temp_genesymbols_volc <- bitr(temp_volc$X, fromType="ENSEMBL", toType=c("SYMBOL","ENSEMBL"), OrgDb="org.Mm.eg.db")

#deseq2 with bam files and count matrix from my experiment to generate volcano plot
countdata <- read.table("/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/counts.tsv", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\.[sb]am$", "", colnames(countdata))

#produce volcano plot
EnhancedVolcano(Mvh_DEGs_Down_DEG,
                lab = Mvh_DEGs_Down_DEG$geneID,
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(Mvh_DEGs_xlsx_DEGs,
                lab = Mvh_DEGs_xlsx_DEGs$geneID,
                x = 'logFC',
                y = 'PValue')

#volcano plots for 3 genes in GSC
EnhancedVolcano(RNA_Mvh_DEGs_raw_pvalue_ajust_pvalue,
                lab = RNA_Mvh_DEGs_raw_pvalue_ajust_pvalue$geneID,
                selectLab = c("Postn",  "Hnrnpk", "Ppia"),
                labSize = 4,
                x = 'logFoldchange',
                y = 'PValue',  
                drawConnectors = TRUE, maxoverlapsConnectors = Inf, arrowheads=TRUE, endsConnectors="last")

#gene enrichment analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(c("MASS", "nlme", "survival"))
BiocManager::install("DirichletMultinomial")
library(DirichletMultinomial)

BiocManager::install("clusterProfiler", force=TRUE)
BiocManager::install("pathview", force=TRUE)
BiocManager::install("enrichplot", force=TRUE)
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

df = Mvh_DEGs_xlsx_DEGs
#df <- read.csv('/gpfs/data/myajima/ylee131/alternative_splicing/GSC/Splicing/Thesis/Mvh_DEGs.xlsx-DEGs.csv', header=TRUE)
#df <- df[,1]

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$geneID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE, force = TRUE)
library(organism, character.only = TRUE)
de <- names(gene_list)

#get universe gene lst
df_uni <- Mvh_DEGs_Sheet1
df_uni <- read.csv('/gpfs/data/myajima/ylee131/alternative_splicing/GSC/Splicing/Thesis/Mvh_DEGs-Sheet1.csv', header=TRUE)
original_gene_list_uni <- df_uni$PValue
names(original_gene_list_uni) <- df_uni$geneID
gene_list_uni<-na.omit(original_gene_list_uni)
gene_list_uni = sort(gene_list_uni, decreasing = TRUE)
de_uni = names(gene_list_uni)

#enrichGO analysis: genelist = downDEG, universe = all genes in sheet 1
enrgo <- enrichGO(gene          = names(gene_list),
                  universe      = de_uni,
                  OrgDb         = organism,
                  keyType = "SYMBOL",
                  ont           = "ALL",
                  readable      = TRUE)

#enrichGO analysis with GSC candidate genes
enrgo_2 <- enrichGO(gene          = c("Postn", "Ppia", "Hnrnpk"),
                    #universe      = de_uni,
                    OrgDb         = organism,
                    keyType = "SYMBOL",
                    ont           = "ALL",
                    readable      = TRUE)
library(enrichplot)
goplot(enrgo_2)
dotplot(enrgo_2, showCategory=30, font.size=10, label_format=100)
#end of candidate gene

goplot(enrgo)
goplot(enrgo,
       showCategory = 30, 
       color = "p.adjust",
       layout = "sugiyama",
       geom = "text")

#get all genes sorted by adjusted p-value
df_all <- RNA_Mvh_DEGs_raw_pvalue_ajust_pvalue
original_gene_list_all <- df_all$adjust_Pvalue
names(original_gene_list_all) <- df_all$geneID
gene_list_all <-na.omit(original_gene_list_all)
gene_list_all = sort(gene_list_all, decreasing = TRUE)
de_all = names(gene_list_all)

#gene list for gse sorted by logFC
df_all_FC <- RNA_Mvh_DEGs_raw_pvalue_ajust_pvalue
original_gene_list_all_FC = df_all_FC[,2]
names(original_gene_list_all_FC) = as.character(df_all_FC[,1])
gene_list_all_FC <-na.omit(original_gene_list_all_FC)
gene_list_all_FC = sort(gene_list_all_FC, decreasing = TRUE)

original_gene_list_all_FC <- df_all_FC$logFoldchange
names(original_gene_list_all_FC) <- df_all_FC$geneID
gene_list_all_FC<-na.omit(original_gene_list_all_FC)
gene_list_all_FC = sort(gene_list_all_FC, decreasing = TRUE)
de_all_FC = names(gene_list_all_FC)
#this works
original_gene_list_all_FC <- df_all_FC$logFoldchange
names(original_gene_list_all_FC) <- df_all_FC$geneID
gene_list_all_FC <-na.omit(original_gene_list_all_FC)
gene_list_all_FC = sort(gene_list_all_FC, decreasing = TRUE)
de_all_FC = names(gene_list_all_FC)
x <- enrichDO(de_all_FC)

gse <- gseGO(geneList=gene_list_all, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

gse_FC <- gseGO(geneList=gene_list_all_FC, 
                ont ="ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")
#dotplot
require(DOSE)
dotplot(enrgo, showCategory=50)
goplot(enrgo)


dotplot(gse, showCategory=30, split=".sign", font.size=10, label_format=100) + facet_grid(.~.sign)
dotplot(gse, showCategory=30, font.size=10, label_format=100)

ego <- pairwise_termsim(gse)
emapplot(ego, showCategory = 30)
geneSetId
gseaplot(gse, by = "all", title = gse$Description[100], geneSetID = 100)

dotplot(gse_FC, x = "GeneRatio", showCategory=30, split=".sign", font.size=10, label_format=100) + facet_grid(.~.sign)
#selected_pathways <- sample(x$Description[1:30])
dotplot(gse_FC, showCategory=30, font.size=10, label_format=100)
dotplot(gse_FC, showCategory = 60, split=".sign") + facet_grid(.~.sign)
act_list <- array()

#for separating activated and suppressed
library(ggplot2)
library(dplyr)
library(stringr)
## count the gene number
gene_count<- gse_FC@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
## merge with the original dataframe
dot_df<- left_join(gse_FC@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
## plot
dot_df = dot_df[1:50,] ## small dataset
dot_df$type[dot_df$NES > 0] = "activated"
dot_df$type[dot_df$NES < 0] = "suppressed"
## from Tommy's code
p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment for suppressed genes")
p + facet_grid(.~type)
#end of separating code

dotplot(gse_FC, showCategory=selected_pathways, font.size=10, label_format=100)
ego <- pairwise_termsim(gse)
emapplot(ego, showCategory = 30)
geneSetId
gseaplot(gse, by = "all", title = gse$Description[100], geneSetID = 100)

#DE analysis with Eric's data
BiocManager::install("DESeq2")
library(DESeq2) 
library(readr)
library(dplyr)
library(ggplot2)
#prepping bam files
BiocManager::install("airway")
library("airway")
BiocManager::install("pasilla")
library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(file="/users/ylee131/data/ylee131/alternative_splicing/counts_byeric/Counts_w_EnsmblGeneIDs.tsv",sep="\t",row.names="gene_id"))

#Generating counts matrix from mine
BiocManager::install("Rsubread")
library("Rsubread")
library("Rsamtools")

#Started coding here
cts <- read.table('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/counts_modified.tsv', sep="\t", header=T,quote="", comment.char="", stringsAsFactors=F)
cts <- as.matrix(cts,sep="\t",row.names="GeneIDs")
cts_symbol <- cts
for (i in 1:length(cts[,1])) {
  cts_symbol[,1][i] <- gsub("\\..*","",cts[,1][i])
}
#temp_for_pca <- bitr(gene_list_to_symbol, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db", drop = FALSE)
#TO KEEP GENEIDS LATER TO COPY FROM
cts_genes <- read.table('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/counts_modified.tsv', sep="\t", header=T,quote="", comment.char="", stringsAsFactors=FALSE)
coldata <- coldata_DEG
cts_genes_symbol <- cts_genes
for (i in 1:length(cts_genes$GeneIDs)) {
  cts_symbol$SYMBOL <- gsub("\\..*","",cts_genes$GeneIDs)
}
coldata
rownames(coldata) <- coldata$...1
# cts <- cts[, rownames(coldata)]
# all(rownames(coldata) == colnames(cts))
# rownames(coldata) <- sub("fb", "", rownames(coldata))
# all(rownames(coldata) %in% colnames(cts))
# library("DESeq2")
# na.omit(cts)
# cts <- cts[complete.cases(cts), ] 
# any(is.na(cts))
# #REPLACE CHARACTER TO NUMERIC VALUES IN COLUMN OF COUNT MATRIX
# #cts[,c(0:3)] <- as.numeric(as.character(cts[,c(0:3)]))
# #cts_copy[,c(0:3)] <- as.numeric(as.character(cts[,c(0:3)]))
# cts_copy <- matrix(as.numeric(as.character(cts[,c(0:4)])),             # Duplicate vector in matrix rows
#        nrow = length(cts[,1]),
#        ncol = 4,
#        byrow = FALSE)
# sapply(cts_copy[,1], class)
# #cts[,1] <- sapply(cts[,1], as.character)
# #cts[,1] <- sapply(cts[,1], as.numeric)
# #cts <- transform(cts, cts[, c(0:3)], as.numeric)
# #ACTUAL FUNCTION TO GET DEG          
# dds <- DESeqDataSetFromMatrix(countData = cts_copy,
#                               colData = coldata,
#                               design = ~ treatment)
# #Add GeneIDs
# featureData <- data.frame(cts_genes[,1])
# featureData <- data.frame(GeneIDs=cts_genes)
# mcols(dds) <- DataFrame(mcols(dds), featureData)
# #Add IDs to count matrix
cts_genes_id <- matrix(as.numeric(as.character(cts_genes[,c(2:5)])),             # Duplicate vector in matrix rows
                       nrow = length(cts_genes[,1]),
                       ncol = 4,
                       byrow = FALSE)
cts_genes_id <- cts_genes
#switch from ENSMBL to symbol
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(clusterProfiler)
gene_list_to_symbol <- cts_genes$GeneIDs
for (i in 1:length(gene_list_to_symbol)) {
  gene_list_to_symbol[i] <- gsub("\\..*","",gene_list_to_symbol[i])
}
temp_for_pca <- bitr(gene_list_to_symbol, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db", drop = FALSE)

gene_list_to_symbol$SYMBOL <- bitr(gene_list_to_symbol, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db", drop = FALSE)$SYMBOL
#cts_genes_id$SYMBOL <- temp_for_pca
#cts_genes_id$SYMBOL <- bitr(gene_list_to_symbol, fromType="ENSEMBL", toType=c("SYMBOL","ENSEMBL"), OrgDb="org.Mm.eg.db")$SYMBOL
#set row names to gene id
rownames(cts_genes_id) <- cts_genes[,1]
#rownames(cts_genes_id) <- gene_list_to_symbol_final$SYMBOL
na.omit(cts_genes_id)
cts_genes_id <- cts_genes_id[complete.cases(cts_genes_id), ] 
any(is.na(cts_genes_id))
cts_genes_id <- cts_genes_id[-c(1) ]
dds <- DESeqDataSetFromMatrix(countData = cts_genes_id,
                              colData = coldata,
                              design = ~ treatment)
#Pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, name="treatment_WT_vs_KO")
res <- results(dds, contrast=c("treatment","WT","KO"))
#shrinking
BiocManager::install("apeglm")
library(apeglm)
resLFC <- lfcShrink(dds, coef="treatment_WT_vs_KO", type="apeglm")
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
#adjusted p-values
resOrdered <- res[order(res$pvalue),]
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")
write.csv(as.data.frame(resOrdered), 
          file="treatment_DEG_results.csv")
#volcano plot with mapped reads
vol <- EnhancedVolcano(treatment_DEG_results,
                       lab = treatment_DEG_results$"...1",
                       x = 'log2FoldChange',
                       y = 'padj')
for (i in 1:length(vol$data$...1)) {
  vol$data$...1[i] <- gsub("\\..*","",vol$data$...1[i])
  vol$data$...1[i] <- bitr(vol$data$...1[i], fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
}
vol$data$...1 <- bitr(vol$data$...1, fromType="ENSEMBL", toType=c("SYMBOL","ENSEMBL"), OrgDb="org.Mm.eg.db")

resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="treatment_DEG_results_pvaluethres.csv")
#heatmap of the most significant genes
library("pheatmap")
#ddsMat <- DESeqDataSetFromMatrix(countData = cts_genes,
#                                 colData = coldata,
#                                 design = ~ cell + dex)
vsd <- vst(dds)
mat <- assay(vsd)[ head(order(res$pvalue),20), ]
mat <- mat - rowMeans(mat)
row.names(mat)[1] <- "Cdc5l"  
row.names(mat)[2] <- "Zfp937"  
row.names(mat)[3] <- "Fth1"  
row.names(mat)[4] <- "Nubp2" 
row.names(mat)[5] <- "3300002I08Rik"
row.names(mat)[6] <- "Zfp953"
row.names(mat)[7] <- "Selenow"
row.names(mat)[8] <- "Prdx2"
row.names(mat)[9] <- "Cpox"
row.names(mat)[10] <- "Gm42980"
row.names(mat)[11] <- "Gemin7"
row.names(mat)[12] <- "Cenps"
row.names(mat)[13] <- "Timm17a"
row.names(mat)[14] <- "Atp5j2"
row.names(mat)[15] <- "Gm3143"
row.names(mat)[16] <- "Gm30400"
row.names(mat)[17] <- "Phf5a"
row.names(mat)[18] <- "Mgl2"
row.names(mat)[19] <- "Cox5a"
row.names(mat)[20] <- "Gm13361"
temp_mat <- mat
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    temp_mat[i,j] = -1 * mat[i,j]
  }
}
df <- as.data.frame(rownames(temp_mat))
pheatmap(temp_mat, fontsize_row=15, fontsize_col=15, fontsize=10)
#plot PCA
BiocManager::install("pcaExplorer")
library(pcaExplorer)
vsd <- vst(dds, blind=FALSE)
options(repr.plot.width =9, repr.plot.height =10)
pcaData <- plotPCA(vsd, intgroup=c("treatment")) + theme(text = element_text(size = 20))
#Obtaining exon counts
BiocManager::install("DEXSeq")
library(DEXSeq)
BiocManager::install("HTSeqCounts")
exon_countfiles <- list.files("/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/", full.names=TRUE, pattern=".counts")
#exon_samplenames <- sapply(
#  strsplit(sapply(
#    strsplit(exon_countfiles, "gr"), 
#    "[[", 1), 
#    "\\/"), 
#  "[[", 1)
exon_conditions <- c("KO", "KO", "WT", "WT")
exon_type <- c("single", "single", "single", "single")
flattenedFile = list.files("/gpfs/data/myajima/ylee131/alternative_splicing/GSC/Splicing/mm10/", pattern="gtf$", full.names=TRUE)
flattenedFile <- flattenedFile[1]
#construct sample table
exon_des <- data.frame(row.names=c("KO1", "KO2", "WT1", "WT2"), condition=exon_conditions, libType=exon_type)
dxd = DEXSeqDataSetFromHTSeq(
  exon_countfiles,
  sampleData=exon_des,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
head( exonIDs(dxd) )
split( seq_len(ncol(dxd)), colData(dxd)$exon )
head( featureCounts(dxd), 5 )
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd )
plotDispEsts( dxd , mar=20)+theme(legend.text = element_text(size=10))
dxd <- testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 <- DEXSeqResults( dxd )
mcols(dxr1)$description
#how many exonic regions are significant with a false discovery rate of 10%
table ( dxr1$padj < 0.05 )
#how many genes are affected
table ( tapply( dxr1$padj < 0.05, dxr1$groupID, any ) )
rows_lowfdr_dxr <- which(dxr1$padj < 0.05)
sig_genes_dxr <- dxr1$groupID[rows_lowfdr_dxr]
plotMA( dxr1, cex=0.8 )
plotDEXSeq( dxr1, "Ddx4", legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
plotDEXSeq( dxr1, "ENSMUST00000122127.1+ENSMUST00000151646.1+ENSMUST00000051364.3+ENSMUST00000122055.1+ENSMUST00000208826.1+ENSMUST00000117222.1+ENSMUST00000119912.1", legend=TRUE, displayTranscripts=FALSE,expression=FALSE, splicing=TRUE, names=FALSE )
DEXSeqHTML( dxr1, FDR=0.05, color=c("#FF000080", "#0000FF80") )
dexseq_genes <- sig_genes_dxr
for (i in 1:length(dexseq_genes)) {
  dexseq_genes[i] <- gsub("\\..*","",dexseq_genes[i])
}
dexseq_genesymbols <- getBM(attributes = c('ensembl_transcript_id', 
                                           'ensembl_gene_id', 
                                           'mgi_symbol',
                                           'external_transcript_name',
                                           'external_gene_name'),
                            filters = 'ensembl_transcript_id', 
                            values = dexseq_genes,
                            mart=ensembl)
#pathway analysis of dexseq
enrgo_dexseq <- enrichGO(gene = dexseq_genesymbols,
                         universe = uni_gene_list,
                         OrgDb         = organism,
                         keyType = "SYMBOL",
                         ont           = "ALL")
dotplot(enrgo_dexseq, showCategory=30, font.size=10, label_format=100)
temp_dexseq <- list()
temp_dexseq_names <- list()
for (i in 1:length(uni_gene_names)) {
  for (j in 1:nrow(dexseq_genesymbols)) {
    if (uni_gene_names[i] == dexseq_genesymbols[j,5]) {
      temp_dexseq[length(temp_dexseq) + 1] <- uni_gene_list[i]
      temp_dexseq_names[length(temp_dexseq_names) + 1] <- uni_gene_names[i]
    }
  }
}
names(temp_dexseq) <- temp_dexseq_names
gsea_dexseq <- na.omit(temp_dexseq)
gsea_dexseq = sort(unlist(gsea_dexseq), decreasing = TRUE)
gse_dexseq <- gseGO(geneList=gsea_dexseq, 
                    ont ="ALL", 
                    keyType = "SYMBOL", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "none")
dotplot(gse_dexseq, showCategory=30, split=".sign", font.size=10, label_format=100) + facet_grid(.~.sign)
cnetplot(gse_dexseq, categorySize="pvalue", showCategory = 40, max.overlaps=50)
cnetplot(gse_dexseq, categorySize="pvalue", showCategory = c("RNA splicing", "ubiquitin protein ligase binding", "glutamatergic synapse", "cytoplasmic vesicle membrane", "carbohydrate homeostasis", "glucose homeostasis", "cytokine-mediated signaling pathway"))
install.packages("ggridges")
library(ggridges)
ridgeplot(gse_dexseq, decreasing=TRUE, showCategory = 30) + labs(x = "enrichment distribution")
#KEGG
dexseq_ids<-bitr(unlist(dexseq_genesymbols[5]), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb=organism)
dexseq_ids_cleaned = dexseq_ids[!duplicated(dexseq_ids[c("SYMBOL")]),]
temp_uni_genes_1 <- uni_genes
for (i in 1:nrow(temp_uni_genes_1)) {
  temp_uni_genes_1[i,1] <- gsub("\\..*","",temp_uni_genes_1[i,1])
}
df_dexseq = temp_uni_genes_1[temp_uni_genes_1$X %in% dexseq_ids_cleaned$ENSEMBL,]
df_dexseq_names <- bitr(unlist(df_dexseq[1]), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
df_dexseq_names = df_dexseq_names[!duplicated(df_dexseq_names[c("ENSEMBL")]),]
df_dexseq$Y = df_dexseq_names$ENTREZID
kegg_gene_list_dexseq <- df_dexseq$log2FoldChange
names(kegg_gene_list_dexseq) <- df_dexseq$Y
kegg_gene_list_dexseq<-na.omit(kegg_gene_list_dexseq)
kegg_gene_list_dexseq = sort(kegg_gene_list_dexseq, decreasing = TRUE)
kk2_dexseq <- gseKEGG(geneList     = kegg_gene_list_dexseq,
                      organism     = "mmu",
                      nPerm        = 10000,
                      minGSSize    = 3,
                      maxGSSize    = 800,
                      pvalueCutoff = 1,
                      pAdjustMethod = "none",
                      keyType       = "ncbi-geneid")
dotplot(kk2_dexseq, showCategory = 10, title = "Enriched Pathways" , split=".sign", font.size=15, label_format=70) + facet_grid(.~.sign)
cnetplot(kk2_dexseq, categorySize="pvalue", showCategory=30)
gseaplot(kk2_dexseq, by = "all", title = kk2_dexseq$Description[34], geneSetID = 34)

#see if there are matches
dexseq_miso_match <- list()
for (i in 1:nrow(dexseq_genesymbols)) {
  for (j in 1:nrow(miso_gene_symbols)) {
    if (dexseq_genesymbols[i, 3] == miso_gene_symbols[j, 2]) {
      dexseq_miso_match[length(dexseq_miso_match) + 1] <- miso_gene_symbols[j,2]
    }
  }
}
dexseq_deg_match <- list()
for (i in 1:nrow(dexseq_genesymbols)) {
  for (j in 1:nrow(final_pval0.001lfc1.5_symbols)) {
    if (dexseq_genesymbols[i, 3] == final_pval0.001lfc1.5_symbols[j, 2]) {
      dexseq_deg_match[length(dexseq_deg_match) + 1] <- final_pval0.001lfc1.5_symbols[j,2]
    }
  }
}
miso_deg_match <- list()
for (i in 1:nrow(miso_gene_symbols)) {
  for (j in 1:nrow(final_pval0.05_symbols)) {
    if (miso_gene_symbols[i, 2] == final_pval0.05_symbols[j, 2]) {
      miso_deg_match[length(miso_deg_match) + 1] <- final_pval0.05_symbols[j,2]
    }
  }
}
all_match <- list()
for (i in 1:length(miso_deg_match)) {
  for (j in 1:length(dexseq_miso_match)) {
    if (miso_deg_match[[i]] == dexseq_miso_match[[j]]) {
      all_match[length(all_match) + 1] <- dexseq_miso_match[[j]]
    }
  }
}

#DTU ANALYSIS
#Pre-filtering to remove low-expression genes
library(data.table)
library(GenomicFeatures)
library(tximport)
library(DEXSeq)
library(stageR)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
samps = data.frame(sample_id = sampleData$ENA_RUN, group = sampleData$paris_classification)
dxr<-na.omit(dxr)
qval = perGeneQValue(dxr)
dxr.g = data.frame(gene = names(qval), qval)
dxr.t = as.data.frame(dxr[, c("featureID","groupID","pvalue")])
#Set up a function to strip away the version numbers in the Ensembl gene and transcript IDs
strp <- function(x) substr(x,1,15)
pScreen = qval
names(pScreen) = strp(names(pScreen))
pConfirmation = matrix(dxr.t$pvalue, ncol=1)
dimnames(pConfirmation) = list(strp(dxr.t$featureID),"transcript")
tx2gene = data.frame(dxr.t[,c("featureID", "groupID")], 
                     dxr.t[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])
# qval used in pScreen, hence pScreenAdjusted=TRUE
library("stageR")
stageRObj = stageRTx(pScreen = pScreen, 
                     pConfirmation = pConfirmation, 
                     pScreenAdjusted = TRUE, 
                     tx2gene = tx2gene[1:2])
stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
dex.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = TRUE)
dex.padj = merge(tx2gene, dex.padj, by.x = c("groupID","featureID"), by.y = c("geneID","txID"))


#heatmap from filtered genes from the table (excel)
df_deg <- read.csv('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/DEG_filtered_0.001.csv', header=TRUE)
df_deg_new <- na.omit(df_deg)
df_deg_new[1]
vsd <- vst(dds)
mat_filtered <- assay(vsd)[head(order(df_deg_new$padj), 30),]
mat_filtered <- mat_filtered - rowMeans(mat_filtered)
#df_deg_new <- as.data.frame(rownames(mat_filtered))
pheatmap(mat_filtered)


#Pathway Analysis - ClusterProfiler & KEGG
library(clusterProfiler)
library('biomaRt')
#mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
#G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","mgi_symbol", "description"),values=geneList_KEGG,mart= mart)
mmu <- search_kegg_organism('Mus musculus', by='scientific_name')
#geneList_KEGG <- list(df_deg_new$'rownames(mat_filtered)')
#CONVERT ENSMBL TO GENE SYMBOL
geneList_KEGG <- rownames(mat_filtered)
data(geneList_KEGG, package="DOSE")
require(org.Mm.eg.db)
for (i in 1:length(geneList_KEGG)) {
  geneList_KEGG[i] <- gsub("\\..*","",geneList_KEGG[i])
}
require(biomaRt)
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
#converted ENSMBL to gene symbols
annot <- getBM(
  attributes = c(
    'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = geneList_KEGG,
  mart = ensembl)
annot <- c(annot)
#converted ENSMBL to ENTREZID ID
entr_genes <- bitr(geneList_KEGG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
#entr_genes <- list(entr_genes$ENTREZID)
#convert to KEGG ID
library(org.Mm.eg.db)
library(MetaboSignal)
#kegg_ids <- MS_convertGene(genes=geneList_KEGG, organism_code='mmu', organism_name='mouse', output = "matrix",
#               orthology = TRUE)
#kegg_ids <- list(kegg_ids[, KEGG_ID])
KEGG_IDs = mget(as.character(entr_genes),ifnotfound=NA)
kegg_ids = c("497097", "19888", "20671", "620009", "27395", "18777", "21399", "108664", "18387", "12421", "100418138", "240690", "319263", "71096", "59014", "76187", "17864", "70675", "170755", "108167613", "73824")

kk <- enrichKEGG(gene         = kegg_ids,
                 organism     = 'mmu',
                 keyType="kegg",
                 pvalueCutoff = 0.05)
mkk <- enrichMKEGG(gene = kegg_ids,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   keyType="kegg",
                   qvalueCutoff = 1)
head(kk)
genes_symbol <- c("Rp1", "Sox17", "Lypla1", "Oprk1","Rb1cc1","Sntg1", "Adhfe1","Mybl1","Sgk3","St18","Atp6v1h","Tcea1","Mrpl15","Vcpip1","Gm6195","Pcmtd1","Xkr4","Rrs1","Gm15452","Gm26901","Snhg6","Gm37381","Gm37225","Gm30414","Gm2147","Gm38008","Gm19026","Gm6123" )

#enrichGO THIS WORKS NOW
temp <- df_deg_new[1:30,c("X", "padj")]
#temp <- df_deg_new[1:56, 1]
for (i in 1:30) {
  temp[i,1] <- gsub("\\..*","",temp[i,1])
}
final_gsymbols <- df_deg_new$padj[1:30]
genesymbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp,
  mart = ensembl)
names(final_gsymbols) <- genesymbols$mgi_symbol
final_gsymbols = sort(final_gsymbols, decreasing = FALSE)
#universe
uni_genes <- read.csv('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/treatment_DEG_results.csv', header=TRUE)
uni_gene_list <- uni_genes$padj
#uni_gene_list <- uni_genes$log2FoldChange
temp_gene <- uni_genes$X
for (i in 1:length(temp_gene)) {
  temp_gene[i] <- gsub("\\..*","",temp_gene[i])
}
temp_genesymbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_gene,
  mart = ensembl)
names(uni_gene_list) <- temp_genesymbols$mgi_symbol
uni_gene_list<-na.omit(uni_gene_list)
uni_gene_list = sort(uni_gene_list, decreasing = TRUE)
#different filters
library(dplyr)
pvalue_0.05 <- filter(treatment_DEG_results, pvalue < 0.05)
lfc_1.5 <- filter(treatment_DEG_results, log2FoldChange > 1.5 | log2FoldChange < -1.5)
pval0.05_lfc1.5 <- filter(lfc_1.5, pvalue < 0.05)
pval0.001_lfc1.5 <- filter(lfc_1.5, pvalue < 0.001)
pvalue_0.001 <- filter(treatment_DEG_results, pvalue < 0.001)
fdr_0.1 <- filter(treatment_DEG_results, padj < 0.1)
fdr_0.05 <- filter(treatment_DEG_results, padj < 0.05)
pval0.001_lfc1.5_fdr0.05 <- filter(pval0.001_lfc1.5, padj < 0.05)
temp_pval_lfc <- pval0.001_lfc1.5_fdr0.05$...1
temp_pval_lfc_padj <- pval0.001_lfc1.5_fdr0.05[c("...1", "padj")]
for (i in 1:length(temp_pval_lfc)) {
  temp_pval_lfc[i] <- gsub("\\..*","",temp_pval_lfc[i])
}
final_pvallfc_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pval_lfc,
  mart = ensembl)
final_pvallfc_symbols_padj <- final_pvallfc_symbols$...1
names(final_pvallfc_symbols) <- genesymbols$mgi_symbol
final_gsymbols = sort(final_gsymbols, decreasing = FALSE)

temp_pval0.05_lfc1.5 <- pval0.05_lfc1.5$...1
for (i in 1:length(temp_pval0.05_lfc1.5)) {
  temp_pval0.05_lfc1.5[i] <- gsub("\\..*","",temp_pval0.05_lfc1.5[i])
}
final_pval0.05_lfc1.5_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pval0.05_lfc1.5,
  mart = ensembl)

temp_pval0.05 <- pvalue_0.05$...1
for (i in 1:length(temp_pval0.05)) {
  temp_pval0.05[i] <- gsub("\\..*","",temp_pval0.05[i])
}
final_pval0.05_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pval0.05,
  mart = ensembl)

temp_lfc1.5 <- lfc_1.5$...1
for (i in 1:length(temp_lfc1.5)) {
  temp_lfc1.5[i] <- gsub("\\..*","",temp_lfc1.5[i])
}
final_lfc1.5_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_lfc1.5,
  mart = ensembl)

temp_pvalue0.001 <- pvalue_0.001$...1
for (i in 1:length(temp_pvalue0.001)) {
  temp_pvalue0.001[i] <- gsub("\\..*","",temp_pvalue0.001[i])
}
final_pvalue0.001_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pvalue0.001,
  mart = ensembl)

temp_pval0.001lfc1.5 <- pval0.001_lfc1.5$...1
for (i in 1:length(temp_pval0.001lfc1.5)) {
  temp_pval0.001lfc1.5[i] <- gsub("\\..*","",temp_pval0.001lfc1.5[i])
}
final_pval0.001lfc1.5_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pval0.001lfc1.5,
  mart = ensembl)

temp_pval0.05lfc1.5 <- pval0.05_lfc1.5$...1
for (i in 1:length(temp_pval0.05lfc1.5)) {
  temp_pval0.05lfc1.5[i] <- gsub("\\..*","",temp_pval0.05lfc1.5[i])
}
final_pval0.05lfc1.5_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_pval0.05lfc1.5,
  mart = ensembl)

temp_fdr0.1 <- fdr_0.1$...1
for (i in 1:length(temp_fdr0.1)) {
  temp_fdr0.1[i] <- gsub("\\..*","",temp_fdr0.1[i])
}
final_fdr0.1_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_fdr0.1,
  mart = ensembl)

temp_fdr0.05 <- fdr_0.05$...1
for (i in 1:length(temp_fdr0.05)) {
  temp_fdr0.05[i] <- gsub("\\..*","",temp_fdr0.05[i])
}
final_fdr0.05_symbols <- getBM(
  attributes = c('ensembl_gene_id',
                 'mgi_symbol'),
  filters = 'ensembl_gene_id',
  values = temp_fdr0.05,
  mart = ensembl)

ego_pval <- enrichGO(gene = gene_list_cutoff,
                     universe = universe = uni_gene_list,
                     OrgDb         = organism,
                     keyType = "SYMBOL",
                     ont           = "ALL")
#enrgo
enrgo_filter <- enrichGO(gene          = genesymbols$mgi_symbol,
                         universe = uni_gene_list,
                         OrgDb         = organism,
                         keyType = "SYMBOL",
                         ont           = "ALL")
#enrgo lfc and padj
enrgo_lfcpadj <- enrichGO(gene          = final_pvallfc_symbols$mgi_symbol,
                          universe = uni_gene_list,
                          OrgDb         = organism,
                          keyType = "SYMBOL",
                          ont           = "ALL")
enrgo_lfcpval <- enrichGO(gene          = final_pval0.05_lfc1.5_symbols$mgi_symbol,
                          universe = uni_gene_list,
                          OrgDb         = organism,
                          keyType = "SYMBOL",
                          ont           = "ALL")
enrgo_pval <- enrichGO(gene          = final_pval0.05_symbols$mgi_symbol,
                       universe = uni_gene_list,
                       OrgDb         = organism,
                       keyType = "SYMBOL",
                       ont           = "ALL")
enrgo_lfc1.5 <- enrichGO(gene          = final_lfc1.5_symbols$mgi_symbol,
                         universe = uni_gene_list,
                         OrgDb         = organism,
                         keyType = "SYMBOL",
                         ont           = "ALL")
enrgo_pval0.001 <- enrichGO(gene          = final_pvalue0.001_symbols$mgi_symbol,
                            universe = uni_gene_list,
                            OrgDb         = organism,
                            keyType = "SYMBOL",
                            ont           = "ALL")
enrgo_pval0.001lfc1.5 <- enrichGO(gene          = final_pval0.001lfc1.5_symbols$mgi_symbol,
                                  universe = uni_gene_list,
                                  OrgDb         = organism,
                                  keyType = "SYMBOL",
                                  ont           = "ALL")
enrgo_pval0.05lfc1.5 <- enrichGO(gene          = final_pval0.05lfc1.5_symbols$mgi_symbol,
                                 universe = uni_gene_list,
                                 OrgDb         = organism,
                                 keyType = "SYMBOL",
                                 ont           = "ALL")
enrgo_fdr0.1 <- enrichGO(gene          = final_fdr0.1_symbols$mgi_symbol,
                         universe = uni_gene_list,
                         OrgDb         = organism,
                         keyType = "SYMBOL",
                         ont           = "ALL")
enrgo_fdr0.05 <- enrichGO(gene          = final_fdr0.05_symbols$mgi_symbol,
                          universe = uni_gene_list,
                          OrgDb         = organism,
                          keyType = "SYMBOL",
                          ont           = "ALL")

dotplot(enrgo_filter, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_lfcpadj, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_lfcpval, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_pval, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_lfc1.5, showCategory=30, font.size=10, label_format=100)
library(stringr)
dotplot(enrgo_pval0.001, showCategory=15, font.size=15, label_format=50)+theme(legend.text = element_text(size=10))
dotplot(enrgo_pval0.001lfc1.5, showCategory=15, font.size=10, label_format=100)
dotplot(enrgo_pval0.05lfc1.5, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_fdr0.1, showCategory=30, font.size=10, label_format=100)
dotplot(enrgo_fdr0.05, showCategory=30, font.size=10, label_format=100)
go_genes_pval0.05 <- enrgo_pval@result$geneID
names(go_genes_pval0.05) <- enrgo_pval@result$Description
#groupgo
ggo_pval0.05 <- groupGO(gene     = final_pval0.05_symbols$mgi_symbol,
                        OrgDb    = organism,
                        keyType = "SYMBOL",
                        level = 3,
                        ont      = "BP",
                        readable = TRUE)
ggo_pval0.05_lfc1.5 <- groupGO(gene     = final_pval0.05_lfc1.5_symbols$mgi_symbol,
                               OrgDb    = organism,
                               keyType = "SYMBOL",
                               level = 3,
                               ont      = "BP",
                               readable = TRUE)
#GSEA
gsea_pval0.001_genes <- -1 * pvalue_0.001$log2FoldChange
names(gsea_pval0.001_genes) <- temp_pvalue0.001
#names(gsea_pval0.05_genes) <- final_pval0.05_symbols$mgi_symbol
temp_gsea <- pvalue_0.001[c("...1", "log2FoldChange")]
for (i in 1:nrow(temp_gsea)) {
  temp_gsea[i,1] <- gsub("\\..*","",temp_gsea[i,1])
}
genesymbols_gsea <- bitr(temp_gsea$...1, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
temp_gsea__ <- temp_gsea[match(genesymbols_gsea[,1], temp_gsea$...1),]
gsea_pval0.001_genes <- merge(genesymbols_gsea,temp_gsea__, by.x="ENSEMBL", by.y="...1")
temp_again <- gsea_pval0.001_genes$log2FoldChange
names(temp_again) <- gsea_pval0.001_genes$SYMBOL
final_gsea <- na.omit(temp_again)
final_gsea = sort(final_gsea, decreasing = TRUE)

temp_1 <- pvalue_0.001$log2FoldChange
names(temp_1) <- temp_pvalue0.001
temp_1 = sort(temp_1, decreasing = TRUE)

ego_gsea_pval0.001 <- gseGO(geneList = temp_1,
                            OrgDb = organism,
                            ont ="ALL", 
                            keyType = "ENSEMBL", 
                            verbose = TRUE,
                            pAdjustMethod = "none")
dotplot(ego_gsea_pval0.001, showCategory=10, split=".sign", font.size=12, label_format=70) + facet_grid(.~.sign)+theme(legend.text = element_text(size=10))

gse_filter <- gseGO(geneList=final_gsymbols, 
                    ont ="ALL", 
                    keyType = "SYMBOL", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "none")
#KEGG for filtered genes
temp_pval0.001_kegg <- pvalue_0.001$...1
for (i in 1:length(temp_pval0.001_kegg)) {
  temp_pval0.001_kegg[i] <- gsub("\\..*","",temp_pval0.001_kegg[i])
}
ids_pval0.001<-bitr(temp_pval0.001_kegg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
dedup_ids_pval0.001 = ids_pval0.001[!duplicated(ids_pval0.001[c("ENSEMBL")]),]
temp_uni_pval0.001 <- uni_genes
for (i in 1:nrow(temp_uni_pval0.001)) {
  temp_uni_pval0.001[i,1] <- gsub("\\..*","",temp_uni_pval0.001[i,1])
}
df_pval0.001 = temp_uni_pval0.001[temp_uni_pval0.001$X %in% dedup_ids_pval0.001$ENSEMBL,]
df_pval0.001$Y = dedup_ids_pval0.001$ENTREZID
kegg_gene_list_pval0.001 <- -1 * df_pval0.001$log2FoldChange
names(kegg_gene_list_pval0.001) <- df_pval0.001$Y
kegg_gene_list_pval0.001<-na.omit(kegg_gene_list_pval0.001)
kegg_gene_list_pval0.001 = sort(kegg_gene_list_pval0.001, decreasing = TRUE)
kk2_pval0.001 <- gseKEGG(geneList     = kegg_gene_list_pval0.001,
                         organism     = "mmu",
                         nPerm        = 10000,
                         minGSSize    = 3,
                         maxGSSize    = 800,
                         pvalueCutoff = 1,
                         pAdjustMethod = "none",
                         keyType       = "ncbi-geneid")
dotplot(kk2_pval0.001, showCategory = 10, title = "Enriched Pathways" , split=".sign", font.size=15, label_format=70) + facet_grid(.~.sign)+theme(legend.text = element_text(size=10))
cnetplot(kk2_pval0.05, categorySize="pvalue", showCategory=kk2_pval0.05$Description[220])
gseaplot(kk2_pval0.001, by = "all", title = kk2_pval0.001$Description[47], geneSetID = 47)
dme_pval0.001 <- pathview(gene.data=kegg_gene_list_pval0.001, pathway.id=kk2_pval0.001$ID[47], species = "mmu")
dme_pval0.05 <- pathview(gene.data=kegg_gene_list_pval0.05, pathway.id=kk2_pval0.05$ID[220], species = "mmu")
#for all genes
temp <- treatment_DEG_results
temp <-na.omit(temp)
for (i in 1:nrow(temp)) {
  temp[i,1] <- gsub("\\..*","",temp[i,1])
}
temp_entrezid <- bitr(temp$...1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
temp__ <- temp[match(temp_entrezid[,1], temp$...1),]
final_entrezid_all_genes <- merge(temp_entrezid,temp__, by.x="ENSEMBL", by.y="...1")
final_entrezid_all_genes_lfc <- final_entrezid_all_genes$log2FoldChange
names(final_entrezid_all_genes_lfc) <- final_entrezid_all_genes$ENTREZID
final_entrezid_all_genes_lfc <-na.omit(final_entrezid_all_genes_lfc)
final_entrezid_all_genes_lfc = sort(final_entrezid_all_genes_lfc, decreasing = TRUE)

gse_filter_all <- gseGO(geneList=final_all_entrezid_padj, 
                        ont ="ALL", 
                        keyType = "ENTREZID", 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "none")

dotplot(gse_filter, showCategory=30, split=".sign", font.size=10, label_format=100) + facet_grid(.~.sign)
dotplot(gse_filter_30, showCategory=30, split=".sign", font.size=10, label_format=100) + facet_grid(.~.sign)

#MISO gene pathway analysis
miso_gene_list <- filtered_BF10$event_name
temp_sss <- uni_gene_list
for (i in 1:length(uni_gene_list)) {
  temp_sss[i] = -1 * uni_gene_list[i]
}
for (i in 1:length(miso_gene_list)) {
  miso_gene_list[i] <- gsub("\\..*","",miso_gene_list[i])
}
enrgo_miso <- enrichGO(gene = miso_gene_list,
                       universe = temp_sss,
                       OrgDb         = organism,
                       keyType = "ENSEMBL",
                       ont           = "ALL")
dotplot(enrgo_miso, showCategory=10, font.size=15, label_format=30)
miso_gene_symbols <- bitr(miso_gene_list, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
uni_gene_list_miso <- na.omit(temp_sss)
uni_gene_names <- uni_gene_list_miso$...1
uni_gene_names <- na.omit(uni_gene_names)
temp_miso <- list()
temp_miso_names <- list()
for (i in 1:length(miso_uni)) {
  for (j in 1:nrow(miso_gene_symbols)) {
    if (names(miso_uni)[i] == miso_gene_symbols[j,1]) {
      temp_miso[length(temp_miso) + 1] <- miso_uni[i]
      temp_miso_names[length(temp_miso_names) + 1] <- names(miso_uni)[i]
    }
  }
}
temp_names_miso <- names(gsea_miso)
temp_4 <- bitr(temp_names_miso, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
#names(temp_miso) <- temp_miso_names
gsea_miso <- na.omit(temp_miso)
gsea_miso = sort(unlist(gsea_miso), decreasing = TRUE)
gse_miso <- gseGO(geneList=gsea_miso,
                  OrgDb = organism,
                  ont ="ALL", 
                  keyType = "SYMBOL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  pAdjustMethod = "none")
dotplot(gse_miso, showCategory=10, split=".sign", font.size=15, label_format=100) + facet_grid(.~.sign)
gse_miso[35]$enrichmentScore = 0.625
gseaplot(gse_miso, by = "all", title = gse_miso$Description[12], geneSetID = 12)
cnetplot(gse_miso, categorySize="pvalue", showCategory = 10, max.overlaps=50)
cnetplot(gse_miso, categorySize="pvalue", showCategory = c("mRNA splicing, via spliceosome"))
install.packages("ggridges")
library(ggridges)
ridgeplot(gse_miso, decreasing=TRUE, showCategory = 30) + labs(x = "enrichment distribution")
#KEGG for MISO
ids<-bitr(temp_pval0.001_kegg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
temp_uni_genes <- uni_genes
for (i in 1:nrow(temp_uni_genes)) {
  temp_uni_genes[i,1] <- gsub("\\..*","",temp_uni_genes[i,1])
}
df2 = temp_uni_genes[temp_uni_genes$X %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$log2FoldChange
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
for (i in 1:length(kegg_gene_list)) {
  kegg_gene_list[i] = -1 * kegg_gene_list[i]
}
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "mmu",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 1,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign", font.size=15, label_format=70) + facet_grid(.~.sign)
cnetplot(kk2, categorySize="pvalue", showCategory=30)
gseaplot(kk2, by = "all", title = kk2$Description[37], geneSetID = 37)
library(pathview)
dme <- pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[37], species = "mmu")