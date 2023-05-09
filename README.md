# Uncovering the Role of RNA Helicase Vasa in Alternative Splicing in Mouse Oocytes Using Computational Methods
This project aims to address three objectives: 1) identify candidate genes that undergo and/or participate in differential alternative splicing in MVH-knockout samples, 2) identify downstream pathways in which VASA regulates, and 3) conduct a comparative analysis with the glioblastoma stem cells in humans based on genes and pathways. If VASA/DDX4 is shown to play a role in alternative splicing, it will facilitate a deeper understanding of the greater network of VASA/DDX4 in mice and potentially humans. 
\
\
The paper can be found here: [thesis paper](files/YLee_final.pdf)
\
\
This page lays out the general workflow of the project. The entire code can be found here: [code](files/oscar_script.R)
\
\
The full GitHub repository can be found here: [repo](https://github.com/ylee131/thesis)
\
\
***Outline:***
1. RNA-seq pre-processing and mapping
2. Differential expression gene analysis
3. Differential isoform expression analysis
4. Differential exon usage analysis
5. Pathway analysis

## RNA-seq pre-processing and mapping
RNA-seq data of two wild type (WT) and two MVH-knockout (KO) provided by the collaborator Dr. Azusa Inoue at RIKEN in Japan as fastq files
### FastQC: quality check
### Trimmomatic
pre-processed and trimmed 
```
TrimmomaticSE -threads 16 \
RNA_Mvh_WT_Rep2.fastq.gz \
RNA_Mvh_WT_Rep2_trimmed.fq \
ILLUMINACLIP:Nextera-PE-PE.fa:2:30:5:6:true \
HEADCROP:15 \
SLIDINGWINDOW:10:25 MINLEN:50
```
### STAR
mapped reads to mm10
```
STAR --genomeDir /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/Trimmed/ \
--runThreadN 16 \
--readFilesIn /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/Trimmed/RNA_Mvh_KO_Rep2_trimmed.fastq \
--outFileNamePrefix /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/Trimmed/KO_2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 
```
### Qualimap
quantified aligned and non-aligned reads
```
module load java/jdk-11
```
```
qualimap rnaseq \
-outdir qualimap/WT1 \
-bam WT_1Aligned.sortedByCoord.out.bam \
-gtf /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/mm10.gtf \
--java-mem-size=8G
```
### featureCounts
generated count matrix
```
featureCounts -T 4 \
-a /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/mm10.gtf \
-o /users/ylee131/data/ylee131/alternative_splicing/fastq_Mvh/Trimmed/sorted_bam/counts.tsv \
*.bam
```

## Differential expression gene analysis (DESeq2)
```
BiocManager::install("DESeq2")
library(DESeq2) 
cts <- read.table('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/counts.tsv', sep="\t", header=T,quote="", comment.char="", stringsAsFactors=F)
cts <- as.matrix(cts,sep="\t",row.names="GeneIDs")
cts_symbol <- cts
for (i in 1:length(cts[,1])) {
  cts_symbol[,1][i] <- gsub("\\..*","",cts[,1][i])
}
cts_genes <- read.table('/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/counts_modified.tsv', sep="\t", header=T,quote="", comment.char="", stringsAsFactors=FALSE)
coldata <- coldata_DEG
cts_genes_symbol <- cts_genes
for (i in 1:length(cts_genes$GeneIDs)) {
  cts_symbol$SYMBOL <- gsub("\\..*","",cts_genes$GeneIDs)
}
coldata
rownames(coldata) <- coldata$...1
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
```
### Filtering genes based on different thresholds
```
library(dplyr)
pvalue_0.05 <- filter(treatment_DEG_results, pvalue < 0.05)
lfc_1.5 <- filter(treatment_DEG_results, log2FoldChange > 1.5 | log2FoldChange < -1.5)
pval0.05_lfc1.5 <- filter(lfc_1.5, pvalue < 0.05)
pval0.001_lfc1.5 <- filter(lfc_1.5, pvalue < 0.001)
pvalue_0.001 <- filter(treatment_DEG_results, pvalue < 0.001)
fdr_0.1 <- filter(treatment_DEG_results, padj < 0.1)
fdr_0.05 <- filter(treatment_DEG_results, padj < 0.05)
pval0.001_lfc1.5_fdr0.05 <- filter(pval0.001_lfc1.5, padj < 0.05)
```

## Differential isoform expression analysis (MISO)
```
singularity exec -B /users/ylee131 MISO.sif bash
slurm-<jobid>.out
index_gff --index SE.mm10.gff indexed/
miso --run indexed/ KO_1Aligned.sortedByCoord.out.bam --output-dir miso_output_KO/ --read-len 61 --settings-filename settings.txt
```

## Differential exon usage analysis (DEXSeq)
```
BiocManager::install("DEXSeq")
library(DEXSeq)
BiocManager::install("HTSeqCounts")
exon_countfiles <- list.files("/gpfs/data/myajima/ylee131/alternative_splicing/bam_count/", full.names=TRUE, pattern=".counts")
exon_samplenames <- sapply(
  strsplit(sapply(
    strsplit(exon_countfiles, "gr"), 
    "[[", 1), 
    "\\/"), 
  "[[", 1)
exon_conditions <- c("KO", "KO", "WT", "WT")
exon_type <- c("single", "single", "single", "single")
flattenedFile = list.files("/gpfs/data/myajima/ylee131/alternative_splicing/GSC/Splicing/mm10/", pattern="gtf$", full.names=TRUE)
flattenedFile <- flattenedFile[1]
#construct sample table
exon_des <- data.frame(row.names=exon_samplenames, condition=exon_conditions, libType=exon_type)
dxd = DEXSeqDataSetFromHTSeq(
  exon_countfiles,
  sampleData=exon_des,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxr <- DEXSeqResults(dxd)
plotDEXSeq( dxr, "ENSMUST00000000001.4", legend=TRUE, displayTranscripts=TRUE )
```

## Pathway analysis (ClusterProfiler)
### ORA analysis using enrichGO
```
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
enrgo_pval0.001 <- enrichGO(gene          = final_pvalue0.001_symbols$mgi_symbol,
                            universe = uni_gene_list,
                            OrgDb         = organism,
                            keyType = "SYMBOL",
                            ont           = "ALL")
dotplot(enrgo_pval0.001, showCategory=15, font.size=15, label_format=50)+theme(legend.text = element_text(size=10))
```
### GSEA analysis
```
gsea_pval0.001_genes <- pvalue_0.001$log2FoldChange
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
```
### KEGG for filtered genes
```
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
kegg_gene_list_pval0.001 <- df_pval0.001$log2FoldChange
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
dotplot(kk2_pval0.05, showCategory = 10, title = "Enriched Pathways" , split=".sign", font.size=15, label_format=70) + facet_grid(.~.sign)+theme(legend.text = element_text(size=10))
cnetplot(kk2_pval0.05, categorySize="pvalue", showCategory=kk2_pval0.05$Description[220])
gseaplot(kk2_pval0.05, by = "all", title = kk2_pval0.05$Description[220], geneSetID = 220)
dme_pval0.05 <- pathview(gene.data=kegg_gene_list_pval0.05, pathway.id=kk2_pval0.05$ID[220], species = "mmu")
```
<sub>Made with :love_letter: by Yuna Lee</sub> \
<sub>Brown University Yajima Lab © [@ylee131](https://github.com/ylee131)</sub>
