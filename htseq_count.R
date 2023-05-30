# load libraries ----------------------------------------------------------
library(readr)
library(DESeq2)
library(pheatmap)
library(ashr)
library(EnhancedVolcano)
library(ggpubr)
library(ggplot2)

# create a metadata file---------------------------------------------------------------
# create and maintain order of sample names
samplenames <- c("SRR17007285","SRR17007286", "SRR17007287", "SRR17007288",
                "SRR17007289", "SRR17007290", "SRR17007291", "SRR17007292",
                "SRR17007293", "SRR17007294", "SRR17007295", "SRR17007296",
                "SRR17007297", "SRR17007298", "SRR17007299", "SRR17007300",
                "SRR17007301", "SRR17007302", "SRR17007303", "SRR17007304",
                "SRR17007305", "SRR17007306", "SRR17007307", "SRR17007308",
                "SRR17007309", "SRR17007310", "SRR17007311", "SRR17007312",
                "SRR17007313", "SRR17007314", "SRR17007315", "SRR17007316",
                "SRR17007317", "SRR17007318", "SRR17007319", "SRR17007320")

# create a column with technical replicate IDs
replicates <- c("SRX13197517","SRX13197517","SRX13197516","SRX13197516","SRX13197515","SRX13197515",
                "SRX13197514","SRX13197514","SRX13197513","SRX13197513","SRX13197512","SRX13197512",
                "SRX13197511","SRX13197511","SRX13197510","SRX13197510","SRX13197509","SRX13197509",
                "SRX13197508","SRX13197508","SRX13197507","SRX13197507","SRX13197506","SRX13197506",
                "SRX13197505","SRX13197505","SRX13197504","SRX13197504","SRX13197503","SRX13197503",
                "SRX13197502","SRX13197502","SRX13197501","SRX13197501","SRX13197500","SRX13197500")

# load sample files
samplefiles <- grep("tabular", list.files("./htseq-count on collection 307-2/"), value = TRUE)


# create and maintain order of sample conditionsset
samplecondition <- c("Kras_Off_PDA_media_3_days","Kras_Off_PDA_media_3_days",
                     "Kras_Off_PDA_media_3_days","Kras_Off_PDA_media_3_days",
                     "Kras_Off_PDA_media_3_days","Kras_Off_PDA_media_3_days",
                     "Kras_Off_PDA_media_5_days","Kras_Off_PDA_media_5_days",
                     "Kras_Off_PDA_media_5_days","Kras_Off_PDA_media_5_days",
                     "Kras_Off_PDA_media_5_days","Kras_Off_PDA_media_5_days",
                     "Kras_On_PDA_media","Kras_On_PDA_media","Kras_On_PDA_media",
                     "Kras_On_PDA_media","Kras_On_PDA_media","Kras_On_PDA_media",
                     "LPS","LPS","LPS","LPS","LPS","LPS","IL4","IL4","IL4",
                     "IL4","IL4","IL4","M_CSF","M_CSF","M_CSF","M_CSF","M_CSF",
                     "M_CSF")


# create sample table using samplenames, samplefiles, and samplecondition
sampletable <- data.frame(sampleName = samplenames,
                          fileName = samplefiles,
                          condition = samplecondition,
                          replicates_ID = replicates)  
  
# set factor levels
sampletable$condition <- factor(sampletable$condition)  
sampletable$replicates_ID <- factor(sampletable$replicates_ID)

#set M_CSF as reference level, to compare all other in vitro conditions to.
sampletable$condition <- relevel(sampletable$condition, ref = "M_CSF")

# generate count matrix  ---------------------------------------
# set path to sample files
# ensure path corresponds to the directory containing the htseq-count data
directory <-  "./htseq-count on collection 307-2/"

# read count data in deseq2 and create count matrix
ddsMat <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                      directory = directory,
                                      design= ~ condition)  

#collapse technical replicates
ddsMatColl <- collapseReplicates(ddsMat, ddsMat$replicates_ID, renameCols = TRUE)


# Filter lowly expressed genes
#examining number of genes with above or equal to 10 reads
dim(ddsMatColl)
dim(ddsMatColl[rowSums(counts(ddsMatColl)) >= 10, ])

#excluding lowly expressed genes
ddsMatCollFil <- ddsMatColl[rowSums(counts(ddsMatColl)) >= 10, ]

# differential gene expression --------------------------------------------
# perform differential gene expression
DiffGene <- DESeq(ddsMatCollFil)
results <- results(DiffGene)
summary(results)

#quality inspection: scatter plots (not very helpful)
pdf("scatterplots.pdf")
pairs(log2(counts(DiffGene)+1))
dev.off()

#quality inspection: PCA _ plotPCA function takes DESeqtransform object
rld <- rlog(DiffGene, blind=FALSE)
plotPCA(rld, intgroup=c("condition"))

# quality insepction: pair-wise correlation 
pdf("heatmap_rlog_correlation.pdf")
pheatmap(cor(assay(rld)))
dev.off()

#Pairwise contrast
#Krass comparison
res_krass <- results(DiffGene, contrast=c("condition", "Kras_Off_PDA_media_3_days", "Kras_On_PDA_media"))
res_krass_ordered <- res_krass[order(res_krass$pvalue),]
res_krass_ordered_shrunk <- lfcShrink(DiffGene, contrast=c("condition","Kras_Off_PDA_media_3_days","Kras_On_PDA_media"), res=res_krass_ordered, type = "ashr")

#macrophage in vitro polarization comparison
res_pol_M1 <- results(DiffGene, contrast=c("condition", "M_CSF", "LPS"))
res_pol_M2 <- results(DiffGene, contrast=c("condition", "M_CSF", "IL4"))

res_pol_M1_ordered <- res_pol_M1[order(res_pol_M1$pvalue),]
res_pol_M2_ordered <- res_pol_M2[order(res_pol_M2$pvalue),]

res_pol_M1_ordered_shrunk <- lfcShrink(DiffGene, contrast=c("condition", "M_CSF", "LPS"), res=res_pol_M1_ordered, type = "ashr")
res_pol_M2_ordered_shrunk <- lfcShrink(DiffGene, contrast=c("condition", "M_CSF", "IL4"), res=res_pol_M2_ordered, type = "ashr")

# explore results
summary(res_krass_ordered_shrunk, alpha=0.05)
summary(res_pol_M1_ordered_shrunk, alpha=0.05)
summary(res_pol_M2_ordered_shrunk, alpha=0.05)

head(res_krass_ordered_shrunk, 10)
head(res_pol_M1_ordered_shrunk, 10)     
head(res_pol_M2_ordered_shrunk, 10)  

##annotation
#Read Biomart annotation file
annotation <- read.table("mart_export.txt",header=TRUE,sep="\t",quote="\"")

#Annotate krass results
results_annot_krass <- merge(as.data.frame(res_krass_ordered_shrunk),annotation,by.x=0,by.y=1)

#Annotae polarization results
results_annot_polM1 <- merge(as.data.frame(res_pol_M1_ordered_shrunk),annotation,by.x=0,by.y=1)
results_annot_polM2 <- merge(as.data.frame(res_pol_M2_ordered_shrunk),annotation,by.x=0,by.y=1)

#MA plots
plotMA(res_krass_ordered_shrunk, ylim = c(-5, 5), main="krass wild-type vs krass mutant")
plotMA(res_pol_M1_ordered_shrunk, ylim = c(-10, 10), main="M0 vs M1 polarization")
plotMA(res_pol_M2_ordered_shrunk, ylim = c(-6, 6), main="M0 vs M2 polarization")

# advanced MAplot
#krass
ggmaplot(res_krass_ordered_shrunk, main = "krass-wild type vs krass mutant",
         fdr = 0.05, fc = 2, size = 0.6,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(results_annot_krass$Gene.name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold", 
         ggtheme = ggplot2::theme_minimal())

#M1
ggmaplot(res_pol_M1_ordered_shrunk, main = "M0 vs M1",
         fdr = 0.05, fc = 2, size = 0.6,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(results_annot_polM1$Gene.name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold", 
         ggtheme = ggplot2::theme_minimal())


#M2
ggmaplot(res_pol_M2_ordered_shrunk, main = "M0 vs M2",
         fdr = 0.05, fc = 2, size = 0.6,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(results_annot_polM2$Gene.name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold", 
         ggtheme = ggplot2::theme_minimal())

##valcano plots

#krass
EnhancedVolcano(res_krass_ordered_shrunk,
                lab = results_annot_krass$Gene.name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "krass-wild-type vs krass mutant")


#M1
EnhancedVolcano(res_pol_M1_ordered_shrunk,
                lab = results_annot_polM1$Gene.name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "M0 vs M1")

#M2
EnhancedVolcano(res_pol_M2_ordered_shrunk,
                lab = results_annot_polM2$Gene.name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "M0 vs M2")

##barplots
#krass
results_annot_krass <- results_annot_krass[order(results_annot_krass$padj),]
bar_krass <- results_annot_krass[1:50, c ( "Gene.name", "log2FoldChange", "padj")]

ggplot(bar_krass) +
  geom_col (aes(x=reorder(Gene.name, - log2FoldChange), y=log2FoldChange, fill=padj), color = "black", width = 1) +
  scale_fill_gradient(low="#2c4bc7", high="#ed0707") +
  geom_hline(yintercept = 0, linetype="solid") +
  geom_hline(yintercept = 0.5:-0.5, linetype="dashed", size=1) +
  labs(title = "Differential expression of the top 50 significant gene - log fold change", subtitle = "krass-wild-type vs krass mutant", x = "Genes") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), plot.background = element_rect(color = "white"),
        panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), 
        axis.text.x=element_text(size=10, angle = 90),
        plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(1,2,1,2, "cm")) + ylim(-5,7.5)  +
  scale_y_continuous(n.breaks = 8)



#M1
results_annot_polM1 <- results_annot_polM1[order(results_annot_polM1$padj),]
bar_pol1 <- results_annot_polM1[1:50, c ( "Gene.name", "log2FoldChange", "padj")]

ggplot(bar_pol1) +
  geom_col (aes(x=reorder(Gene.name, - log2FoldChange), y=log2FoldChange, fill=padj), color = "black", width = 1) +
  scale_fill_gradient(low="#2c4bc7", high="#ed0707") +
  geom_hline(yintercept = 0, linetype="solid") +
  geom_hline(yintercept = 0.5:-0.5, linetype="dashed", size=1) +
  labs(title = "Differential expression of the top 50 significant gene - log fold change", subtitle = "M0 vs M1", x = "Genes") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), plot.background = element_rect(color = "white"),
        panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), 
        axis.text.x=element_text(size=10, angle = 90),
        plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(1,2,1,2, "cm")) + ylim(-5,7.5)  +
  scale_y_continuous(n.breaks = 8)


#M2
results_annot_polM2 <- results_annot_polM2[order(results_annot_polM2$padj),]
bar_pol2 <- results_annot_polM2[1:50, c ( "Gene.name", "log2FoldChange", "padj")]

ggplot(bar_pol2) +
  geom_col (aes(x=reorder(Gene.name, - log2FoldChange), y=log2FoldChange, fill=padj), color = "black", width = 1) +
  scale_fill_gradient(low="#2c4bc7", high="#ed0707") +
  geom_hline(yintercept = 0, linetype="solid") +
  geom_hline(yintercept = 0.5:-0.5, linetype="dashed", size=1) +
  labs(title = "Differential expression of the top 50 significant gene - log fold change", subtitle = "M0 vs M2", x = "Genes") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), plot.background = element_rect(color = "white"),
        panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), 
        axis.text.x=element_text(size=10, angle = 90),
        plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(1,2,1,2, "cm")) + ylim(-5,7.5)  +
  scale_y_continuous(n.breaks = 8)


# store all results in text file
write_csv(results_annot_krass, "")
write_csv(results_annot_polM1, "")
write_csv(results_annot_polM2, "")

