###########################################################
### Main project #2: RNA-seq analysis after Salmon step ###
###########################################################

library(ensembldb)
library(AnnotationHub)
library(stringr)
library(dplyr)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tximeta)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(grid)
library(apeglm)
library(ashr)
library(Rtsne)
library(M3C)
library(reshape2)
library(ggpubr)
library(DESeq2)
library(edgeR)
library(clusterProfiler)
library(enrichR)
library(EnhancedVolcano)
library(tidyverse)

options(ggrepel.max.overlaps = Inf)


my.theme <- theme(axis.text = element_text(colour="black", size=15, face = "bold"),
                  text = element_text(size=16),
                  panel.background = element_rect(fill = 'gray99',
                                                  colour = "black",
                                                  linewidth=0.5),
                  axis.title.x=  element_text(vjust=-0.45, face = "bold"),
                  axis.title.y = element_text(vjust=1.2, face = "bold"),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line(),
                  panel.grid.major = element_line(colour = "lightgray", linetype="dotted"),
                  panel.grid.minor = element_line(colour = "lightgray", linetype="dashed"),
                  legend.title=element_blank(),
                  legend.text = element_text(size = 14))

#######################################################################################################
##### PART 1. - generate TPM dataset to be used for generating context-specific metabolic models. #####
#######################################################################################################


#connect ENSMUST to ENSMUSG IDs: (needed for ftINIT step)
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus")
edb = ah[["AH109655"]]
txs <- transcripts(edb, return.type = "DataFrame")
tx2gene = as.data.frame(cbind(txs@listData[["tx_id"]], txs@listData[["gene_id"]]))
names(tx2gene) = c("tx_id", "gene_id")

path = "/Users/rokosango/PhD/TRM_macs/TRM_RNAseq/salmon_results/"

paths4allFiles = paste0(path,
       list.files('/Users/rokosango/PhD/TRM_macs/TRM_RNAseq/salmon_results/'),
       "/quant.sf")

SampleIDs = str_extract(list.files('/Users/rokosango/PhD/TRM_macs/TRM_RNAseq/salmon_results/'),
            "[^_]*_[^_]*") #extract substring before a second underscore

sampleList = list()
#load all samples
# for(i in 1:length(paths4allFiles)) {                             
#   assign(SampleIDs[i],                                   
#          read.csv2(paths4allFiles[i],
#                    sep = "\t"))
# }

for(i in 1:length(paths4allFiles)) {                            
  sampleList[[i]] = read.csv2(paths4allFiles[i],
                   sep = "\t")
  names(sampleList)[i] = SampleIDs[i]
  
}

for(i in 1:length(sampleList)) {                            

  sampleList[[i]] = sampleList[[i]] %>%
    dplyr::mutate(TxName = sampleList[["AT_1"]]$Name) %>%
    mutate(TxName = sub("\\..*", "", TxName)) %>% #select everything before a dot
    mutate(GeneName = tx2gene$gene_id[match(TxName, tx2gene$tx_id)])
}

geneIDs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= sampleList[["AT_1"]]$GeneName, 
                             keytype = "GENEID", columns = c("SYMBOL","GENEID"))

for(i in 1:length(sampleList)) {
  
sampleList[[i]] = sampleList[[i]] %>%
  select(c(TxName, GeneName, TPM)) %>%
  group_by(GeneName) %>% #line 61-63: find identical ENSEMBL gene names, select the one with highest TPM and proceed to SYMBOL conversion
  slice_max(TPM) %>% 
  ungroup() %>% 
  distinct(GeneName, .keep_all = TRUE) %>%
  filter(GeneName %in% geneIDs$GENEID) %>%
  mutate(Symbol = ensembldb::select(EnsDb.Mmusculus.v79, keys = GeneName, 
                                   keytype = "GENEID", columns = c("SYMBOL","GENEID"))[, 1])

}

mutate(TxName = sub("\\..*", "", TxName))

test_match_order2 <- function(x,y) {
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

for (i in 1:length(sampleList)) {
  test_match_order2(sampleList[[i]]$GeneName, sampleList[[i]]$GeneName)
  
}

TPM4AllSamples = lapply(sampleList, `[`, "TPM") %>% #select TPM column from each list element
  as.data.frame() %>%
  `colnames<-`(SampleIDs) %>% #rename columns in one go
  mutate(GeneName = sampleList[["AT_1"]]$GeneName) %>%
  relocate(GeneName, .before = AT_1) %>%
  mutate(Symbol = sampleList[["AT_1"]]$Symbol) %>%
  relocate(Symbol, .after = GeneName) %>%
  mutate(Symbol = sub("\\..*", "", Symbol))

for (i in 1:length(sampleList)) {
  print(sampleList[[i]]$TPM[1])
  
}

#find best way to map as many genes: ensembldb, mapIds:

#ensembldb:
(dim(geneIDs)[1] / dim(test2)[1]) * 100 = 82.1 #percent
#mapIds:
table(duplicated(test2$SYMBOL))[1] / sum(table(duplicated(test2$SYMBOL))) * 100 = 73.3 #percent


# ensembldb:in 35627 individual entries for ENSEMBL GENE entries extracted from ENSEMBL TX entries,
# 29240 non-duplicated SYMBOL entries were recovered: 29240 / 35627 = 82%, 20 entries duplicates
# in mapIds function approach, 26144 / 35627 = 73% percent of gene SYMBOL names were recovered, with 9483 duplicates
# Hence, go with ensembldb approach.

#need to convert columns storing TPM values from "character" to "numeric":
TPM4AllSamples[, 3:dim(TPM4AllSamples)[2]] = sapply(TPM4AllSamples[, 3:dim(TPM4AllSamples)[2]], as.numeric)

write.csv(TPM4AllSamples, "TPM4AllSamples.csv", row.names = F)

colSums(TPM4AllSamples[, 3:dim(TPM4AllSamples)[2]])

#next step: feed TPM4AllSamples into ftINIT (MATLAB step).
#step after: import CompMat from MATLAB after generating ftINIT models for all 52 samples.

compMat = read.csv('/Users/rokosango/PhD/TRM_macs/Mouse-GEM/compMat.csv', header = F)
names(compMat) = SampleIDs

tsne(compMat, perplex = 10, labels = coldata$tissue, scale = 3, controlscale=TRUE, axistextsize = 1, seed = 1234)

NumRxns4AllModels = read.csv('/Users/rokosango/PhD/TRM_macs/Mouse-GEM/NumRxns4AllModels.csv', header = T)

rownames(NumRxns4AllModels) = SampleIDs
NumRxns4AllModels$Tissue = coldata$tissue

ggplot(NumRxns4AllModels, aes(x = rownames(NumRxns4AllModels), 
                              y = NumReactions, fill = Tissue)) +
  geom_bar(stat = "identity", color="black", width = 2) +
  labs(x = NULL, y = '# of Reactions in each model') +
  theme(axis.text.x = element_blank()) + my.theme
  

#######################################################################
##### PART 2. - differential gene expression with DESEq2 & edgeR. #####
#######################################################################

filt <- c('_5', '_6', '_7', '_8')

coldata = data.frame(
  names = SampleIDs,
  files = paths4allFiles,
  tissue = rep(c("AT", "BMDM", "Brain", "Colon", "Liver", "Lung", "Monocytes", "Peritoneal", "Spleen"),
               times = c(4, 4, 4, 8, 8, 4, 4, 8, 8)),
  sex = ifelse(grepl(paste(filt, collapse='|'), coldata$SampleName), "Female", "Male")
)

se <- tximeta(coldata)
gse <- summarizeToGene(se)

##############
### DESeq2 ###
##############

dds <- DESeqDataSet(gse, design = ~ tissue + sex)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
nrow(dds) #35627 initial, reduced to 20431 genes
vsd <- vst(dds, blind = FALSE)


# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix( sampleDists )
# rownames(sampleDistMatrix) <- paste( vsd$sex, vsd$tissue, sep = " - " )
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

options(ggrepel.max.overlaps = Inf)


pca = plotPCA(vsd, intgroup = c("sex", "tissue"), returnData = T)
percentVar <- round(100 * attributes(pca)[4]$percentVar)

ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = sex, colour = tissue), size = 2.5) + 
  stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = tissue)) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  my.theme

dds <- DESeq(dds)

uniq = unique(gse$tissue)
Comparison = c()

for (i in 1:length(uniq)) {
  Comparison[[i]] = print(paste(uniq[i], "vs", uniq[!uniq %in% uniq[i]]))
  names(Comparison)[i] = uniq[i]
}

DESeqRes4allComparisons = list()

for (i in names(Comparison)) {
  DESeqRes4allComparisons[[i]] = list()
}

for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    DESeqRes4allComparisons[[j]][[i]] = results(dds, contrast = c("tissue", word(Comparison[[j]][[i]], 1),
                                                                       word(Comparison[[j]][[i]], 3)),
                                      independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    
    DESeqRes4allComparisons[[j]][[i]] = lfcShrink(dds, contrast = c("tissue", word(Comparison[[j]][[i]], 1),
                                                                         word(Comparison[[j]][[i]], 3)),
                                        res = DESeqRes4allComparisons[[j]][[i]], type = "ashr")
    
    names(DESeqRes4allComparisons[[j]])[[i]] = Comparison[[j]][[i]]
  }
}



#############
### edgeR ###
#############

# directly from gse object after tximeta gene-level summarization:

cts <- assays(gse)[["counts"]]
normMat <- assays(gse)[["length"]]
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts, group = coldata$tissue)
y <- calcNormFactors(y)
y$samples

norm_counts <- cpm(y)
head(norm_counts)

y <- estimateDisp(y)


EdgeRes4allComparisons = list()

for (i in names(Comparison)) {
  EdgeRes4allComparisons[[i]] = list()
}

for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    EdgeRes4allComparisons[[j]][[i]] = exactTest(y, pair = c(word(Comparison[[j]][[i]], 3),
                                                             word(Comparison[[j]][[i]], 1)))
    
    EdgeRes4allComparisons[[j]][[i]] = topTags(EdgeRes4allComparisons[[j]][[i]], 
                                               n = dim(EdgeRes4allComparisons[[j]][[i]][["table"]])[1],
                                               sort.by = "PValue", p.value = 0.1)
    
    EdgeRes4allComparisons[[j]][[i]] = EdgeRes4allComparisons[[j]][[i]][["table"]]
    names(EdgeRes4allComparisons[[j]])[[i]] = Comparison[[j]][[i]]
  }
}

MetabolicGenes = read.csv("MouseGemMetabolicGenes.csv")
MetabolicGenes = ensembldb::select(EnsDb.Mmusculus.v79, keys = MetabolicGenes$Gene, 
                                   keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

CombinedRes4allComparisons = list()

for (i in names(Comparison)) {
  CombinedRes4allComparisons[[i]] = list()
}

for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    CombinedRes4allComparisons[[j]][[i]] = merge(as.data.frame(DESeqRes4allComparisons[[j]][[i]]), 
                                                 EdgeRes4allComparisons[[j]][[i]],
                                                 by = 'row.names', all = F)
    
    CombinedRes4allComparisons[[j]][[i]] = subset(CombinedRes4allComparisons[[j]][[i]],
                                               FDR < 0.1 & padj < 0.1)
    
    Genes = ensembldb::select(EnsDb.Mmusculus.v79, 
                                        keys = as.character(CombinedRes4allComparisons[[j]][[i]]$Row.names), 
                                        keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    
    CombinedRes4allComparisons[[j]][[i]] = CombinedRes4allComparisons[[j]][[i]] %>%
      filter(Row.names %in% Genes$GENEID) %>%
      mutate(Symbol = Genes$SYMBOL) %>%
      arrange(desc(log2FoldChange))
      #filter(abs(log2FoldChange) > 2)
      #filter(Row.names %in% MetabolicGenes$GENEID)
      
    
    names(CombinedRes4allComparisons[[j]])[[i]] = Comparison[[j]][[i]]
  }
}
  

#create bird's eye view of L2FC different for each tissue vs all other tissues:
  
l2fcMatrices = function(x, Heatmap=TRUE) {
  
  df_list = list()
  
  for(i in 1:(length(Comparison)-1)) {                                    
    new_element <- CombinedRes4allComparisons[[x]][[i]]      # Create new list element
    df_list[[length(df_list) + 1]] <- new_element               # Append new list element
    
  }
  
  new_names = sub(".*? ", "", names(CombinedRes4allComparisons[[x]])) #remove everything before first space
  
  test_df = df_list %>% 
    reduce(inner_join, by='Row.names') %>% #join all individual pairwise comparisons with a common ENSEMBL gene naames
    arrange(Row.names) %>%
    column_to_rownames(., var = "Row.names") %>%
    dplyr::select(starts_with("log2Fold")) %>%
    dplyr::select(1:8)
  
  test_df = test_df %>%
    dplyr::mutate(Symbol = ensembldb::select(EnsDb.Mmusculus.v79, 
                                             keys = as.character(rownames(test_df)), 
                                             keytype = "GENEID", columns = c("SYMBOL","GENEID"))[,1]) %>%
    dplyr::select(Symbol, everything())
  
  test_df$MeanL2FC = rowMeans(test_df[,2:9])
  
  test_df = test_df %>% arrange(desc(MeanL2FC))
  
  names(test_df) = c("Symbol", new_names, "MeanL2FC")
  
  if (Heatmap) {
    
    hm = pheatmap(test_df,
                  cluster_rows = T,
                  cluster_cols = F,
                  show_rownames = F,
                  fontsize_row = 5,
                  show_colnames = T,
                  main = x,
                  fontsize_col = 12,
                  cellwidth = 15,
                  angle_col = 45
    )
    
    return(hm)
    
  } else {
    
    return(test_df)
    
  }
  
}

for (i in unique(coldata$tissue)) {
  l2fcMatrices(i, Heatmap = T)
}

TissueL2FCMatList = list()

for (i in unique(coldata$tissue)) {
  TissueL2FCMatList[[i]] = l2fcMatrices(i, Heatmap = F)
}

#Volcano Plots for a few overlapping metabolic genes
  
VolcanoPlots = list()

for (i in names(Comparison)) {
  VolcanoPlots[[i]] = list()
}

for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    VolcanoPlots[[j]][[i]] = CombinedRes4allComparisons[[j]][[i]] %>%
      dplyr::filter(Row.names %in% rownames(TissueL2FCMatList[[j]]))
    
    names(VolcanoPlots[[j]])[[i]] = Comparison[[j]][[i]]

    VolcanoPlots[[j]][[i]] = EnhancedVolcano(VolcanoPlots[[j]][[i]],
                    lab =  VolcanoPlots[[j]][[i]]$Symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = names(VolcanoPlots[[j]])[[i]],
                    subtitle = NULL,
                    caption = paste(dim(VolcanoPlots[[j]][[i]])[1], "overlapping metabolic genes", sep = " "),
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 6.0)

  }
}


RadarPlots = function(tissue, n) {
  
  test = TissueL2FCMatList[[tissue]] %>%
    top_n(5) %>%
    select(-MeanL2FC) %>%
    remove_rownames() %>%
    column_to_rownames(., var = "Symbol")
  
  test["Max", ] = max(test) + 1
  test["Min", ] = 0
  
  test = rbind(test["Max", ],
               test["Min", ],
               head(test, -2))
  
  radarchart(test, axistype=1 , 
             plwd=3 , plty=1,
             cglcol="grey", 
             cglty=3, 
             axislabcol="black",
             caxislabels=seq(0,20,5), cglwd=2,
             vlcex=1,
             title = tissue)
  
  legend("bottomleft",
         legend = rownames(test)[3:length(rownames(test))],
         bty = "n", pch = 20, col = 1:5,
         text.col = "grey25", pt.cex = 2)
  
}

RadarPlots("AT", 5)


#Barplots for top n up and down metabolic genes

TopExprGenesBarplot = function(tissue, nGenes) {
  
  testOnlyPos = TissueL2FCMatList[[tissue]][rowSums(TissueL2FCMatList[[tissue]][2:9] > 0) == 8,]
  testOnlyNeg = TissueL2FCMatList[[tissue]][rowSums(TissueL2FCMatList[[tissue]][2:9] < 0) == 8,] %>% arrange(MeanL2FC)
  
  testOnlyPos = testOnlyPos %>%
    top_n(nGenes) %>%
    select(-MeanL2FC) %>%
    mutate(Status = "DiffExprUp") %>%
    melt() %>%
    rename(Tissue = variable, Log2FoldChange = value) %>%
    ggplot(., aes(x=reorder(Symbol, Log2FoldChange), y=Log2FoldChange, fill=Tissue)) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette="Dark2") +
    labs(x = NULL, title = "Expression Up") +
    coord_flip() +
    theme(axis.text = element_text(colour="black", size=10, face = "bold"),
          text = element_text(size=12),
          panel.background = element_rect(fill = 'gray99',
                                          colour = "black",
                                          linewidth=0.5),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.ticks = element_line(colour="black"),
          axis.line = element_line(),
          panel.grid.major = element_line(colour = "lightgray", linetype="dotted"),
          panel.grid.minor = element_line(colour = "lightgray", linetype="dashed"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          title = element_text(size = 10, face = "bold"))
  
  
  testOnlyNeg = testOnlyNeg %>%
    top_n(-nGenes) %>%
    select(-MeanL2FC) %>%
    mutate(Status = "DiffExprNeg") %>%
    melt() %>%
    rename(Tissue = variable, Log2FoldChange = value) %>%
    ggplot(., aes(x=reorder(Symbol, -Log2FoldChange), y=Log2FoldChange, fill=Tissue)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_brewer(palette="Dark2") +
    labs(x = NULL, title = "Expression Down") +
    theme(axis.text = element_text(colour="black", size=10, face = "bold"),
          text = element_text(size=12),
          panel.background = element_rect(fill = 'gray99',
                                          colour = "black",
                                          linewidth=0.5),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.ticks = element_line(colour="black"),
          axis.line = element_line(),
          panel.grid.major = element_line(colour = "lightgray", linetype="dotted"),
          panel.grid.minor = element_line(colour = "lightgray", linetype="dashed"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          title = element_text(size = 10, face = "bold"))
  
  fig = ggarrange(testOnlyPos, 
                  testOnlyNeg,
                  ncol = 2, nrow = 1,
                  common.legend = T,
                  legend = "bottom")
  fig = annotate_figure(fig, top = text_grob(paste(tissue, ": Top", nGenes, "Genes"), face = "bold"))
  fig
  
}

TopExprGenesBarplot("AT", 10)


#search for genes provided by the gene list from Peter Murray:
#look in each individual comparison of unique TRM vs BDMDM

PeterMurrayGeneList = read.csv('PeterMurrayGeneList.csv')
PeterMurrayGeneList$Symbol = str_to_title(PeterMurrayGeneList$Symbol)#convert all uppercase -> first letter upper everythine else low


BMDM_list = list()
BMDM_Volcano_Plots = list()
for (i in names(CombinedRes4allComparisons[["BMDM"]])) {
  BMDM_list[[i]] = CombinedRes4allComparisons[["BMDM"]][[i]] %>%
    dplyr::filter(CombinedRes4allComparisons[["BMDM"]][[i]]$Symbol %in% PeterMurrayGeneList$Symbol)
  
  BMDM_list[[i]] = BMDM_list[[i]] %>%
    dplyr::filter(abs(BMDM_list[[i]]$log2FoldChange) > 2) %>%
    arrange(log2FoldChange)
  
}

for (i in 1:length(BMDM_list)) {
  
  BMDM_Volcano_Plots[[i]] = EnhancedVolcano(BMDM_list[[i]],
                                            lab =  BMDM_list[[i]]$Symbol,
                                            x = 'log2FoldChange',
                                            y = 'padj',
                                            title = names(BMDM_list)[i],
                                            subtitle = NULL,
                                            caption = paste(dim(BMDM_list[[i]])[1], "genes", sep = " "),
                                            xlab = NULL,
                                            ylab = NULL,
                                            FCcutoff = 0.5,
                                            pCutoff = 0.05,
                                            pointSize = 3.0,
                                            labSize = 6.0,
                                            drawConnectors = TRUE)
                                            #widthConnectors = 0.75)
  
  names(BMDM_Volcano_Plots)[i] = names(BMDM_list)[i]
  
}

fig = ggarrange(BMDM_Volcano_Plots[[1]],
          BMDM_Volcano_Plots[[2]],
          BMDM_Volcano_Plots[[3]],
          BMDM_Volcano_Plots[[4]],
          ncol = 2, nrow = 2,
          common.legend = T,
          legend = "top")

annotate_figure(fig, left = textGrob(bquote(~-Log[10] ~ italic(P)), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob(bquote(~Log[2] ~ "fold change"), gp = gpar(cex = 1.3)))

fig = ggarrange(BMDM_Volcano_Plots[[5]],
          BMDM_Volcano_Plots[[6]],
          BMDM_Volcano_Plots[[7]],
          BMDM_Volcano_Plots[[8]],
          ncol = 2, nrow = 2,
          common.legend = T,
          legend = "top")

annotate_figure(fig, left = textGrob(bquote(~-Log[10] ~ italic(P)), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob(bquote(~Log[2] ~ "fold change"), gp = gpar(cex = 1.3)))

########################################################################
############ clustering of edgeR normalized counts: ####################
############  #1. change from ENSEMBL to SYMBOL ########################
############  #2. transpose so features (genes) are columns ############
############  #3. filter for those found significant in DiffExpr #######
############  #4. filter for high-variance genes #######################
############  #5. scale the dataset ####################################
########################################################################

Symbol = ensembldb::select(EnsDb.Mmusculus.v79, 
                                         keys = as.character(rownames(norm_counts)), 
                                         keytype = "GENEID", columns = c("SYMBOL","GENEID"))

test_norm_counts = norm_counts[rownames(norm_counts) %in% Symbol$GENEID,]
rownames(test_norm_counts) = Symbol$SYMBOL

test_norm_counts = t(test_norm_counts)

OuterJoinDegs = list()

for (i in unique(coldata$tissue)) {
  OuterJoinDegs[[i]] = CombinedRes4allComparisons[[i]] %>% 
    reduce(full_join, by='Symbol') %>%
    dplyr::select(Symbol)
}

OuterJoinDegs = unlist(OuterJoinDegs) %>% unique()

test_norm_counts = test_norm_counts[, colnames(test_norm_counts) %in% OuterJoinDegs]
test_norm_counts = test_norm_counts[, !duplicated(colnames(test_norm_counts))]

varVector = apply(test_norm_counts, 2, var)

summary(varVector)
index = varVector > summary(varVector)[3]
index = index[index == T]

test_norm_counts = test_norm_counts[,colnames(test_norm_counts) %in% names(index)]
dim(test_norm_counts) #8293 genes

scaled_test_norm_counts = scale(test_norm_counts)
#find NA columns:
colnames(scaled_test_norm_counts)[ apply(scaled_test_norm_counts, 2, anyNA) ]
# Elbow method
fviz_nbclust(scaled_test_norm_counts, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(scaled_test_norm_counts, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
set.seed(123)
fviz_nbclust(scaled_test_norm_counts, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#every method showed 9 clusters to be optimal
set.seed(123)
km.res <- kmeans(scaled_test_norm_counts, 9, nstart = 50)
km.res

###############################
######## cluster on genes: ####
###############################
fviz_nbclust(scale(t(test_norm_counts)), kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
fviz_nbclust(scale(t(test_norm_counts)), kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
set.seed(123)
fviz_nbclust(scaled_test_norm_counts, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

set.seed(123)
km.res.genes <- kmeans(scale(t(test_norm_counts)), 10, nstart = 50)
km.res.genes$size

t.test_norm_counts <- cbind(t(test_norm_counts), cluster = km.res.genes$cluster)
##############################################################
# clustering for genes done, continue with sample clustering #
##############################################################

scaled_test_norm_counts = as.data.frame(scaled_test_norm_counts)
scaled_test_norm_counts = t(scaled_test_norm_counts)

paletteLength <- 256
#setting 0-point for pheatmap in R:
myBreaks <- c(seq(min(scaled_test_norm_counts), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_test_norm_counts)/paletteLength, max(scaled_test_norm_counts), length.out=floor(paletteLength/2)))

annot = data.frame(
  Population = coldata$tissue,
  Cluster = km.res$cluster
)

annot = annot %>% arrange(Cluster)

annot$Cluster = factor(annot$Cluster, levels = unique(annot$Cluster))
annot$Population = factor(annot$Population, levels = unique(annot$Population))

ann_colors = list(
   Population = c(Monocytes = "#A6761D", Peritoneal = "#666666", AT = "#1B9E77", Colon = "#E7298A",
                  Brain = "#7570B3", Liver = "#66A61E", BMDM = "#D95F02",Lung = "#E6AB02", Spleen = "#D53E4F"),
   Cluster = c(`1` = "#8DD3C7", `2` = "#FFFFB3", `3` = "#BEBADA", `4` = "#FB8072", `5` = "#80B1D3",
               `6` = "#FDB462", `7` = "#B3DE69", `8` = "#FCCDE5", `9` = "#D9D9D9")
)

scaled_test_norm_counts = scaled_test_norm_counts[, rownames(annot)]
pheatmap(scaled_test_norm_counts,
         color = rev(colorRampPalette(c("red", "white", "blue"))(paletteLength)),
         breaks=myBreaks,
         cluster_rows = T, 
         cluster_cols = F,
         show_colnames = F, 
         show_rownames = F,
         angle_col = 45,
         annotation_col = annot,
         annotation_colors = ann_colors,
         treeheight_row = 0, 
         treeheight_col = 0,
         fontsize = 8,
         cellwidth = 5,
         cellheight = 0.04)

MatGenes = data.frame(Genes = rownames(scaled_test_norm_counts))


TopNMarkerGeneHeatmap = function(n) {
  
  MarkerGenes = list()
  for (i in 1:9) {
    
    MarkerGenes[[i]] = data.frame(Centroid = km.res$centers[i, ]) %>%
      arrange(desc(Centroid)) %>%
      top_n(n)
    
    MarkerGenes[[i]] = MarkerGenes[[i]] %>%
      dplyr::mutate(Gene = rownames(MarkerGenes[[i]])) %>%
      dplyr::select(Gene)
    
    names(MarkerGenes)[i] = paste("Cluster", i, sep = "_")
    
  }
  
  MarkerGeneVector = unlist(MarkerGenes) %>% unique()
  
  pheatmap(scaled_test_norm_counts[MarkerGeneVector,], 
           cluster_rows = T, 
           cluster_cols = F,
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           treeheight_row = 0, 
           treeheight_col = 0,
           annotation_col = annot,
           fontsize_row = 9,
           fontsize_col = 8,
           annotation_colors = ann_colors,
           cellheight = 10,
           cellwidth = 10,
           legend = T,
           main = paste("Top", n, "genes in each cluster"))
}

TopNMarkerGeneHeatmap(5)

#marker genes from Lavin et al. 2014 Cell:

MarkersToCheck = c("Sall1", "Clec4f", "Spic", "Car4", "Tgfb2", "Cd74", "Ccr2", "Irf7", "Vcam1",
                   "Siglech", "Pparg", "Cybb")

pheatmap(scaled_test_norm_counts[MarkersToCheck,], 
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T, 
         show_colnames = T,
         angle_col = 45,
         treeheight_row = 0, 
         treeheight_col = 0,
         annotation_col = annot,
         fontsize_row = 9,
         fontsize_col = 8,
         annotation_colors = ann_colors,
         cellheight = 10,
         cellwidth = 10,
         legend = T,
         main = "Genes from Lavin et al 2014 Cell")


# Subset the heatmap by the genes found from Metabolic Atlas
# and see how their expression varies across populations:

MetAtlasGlycolysis = read.csv("Spreadsheets/MetAtlasGlycolysis.tsv",
                              sep = "\t")
MetAtlasOxPhos = read.csv("Spreadsheets/MetAtlasOxPhos.tsv",
                              sep = "\t")
MetAtlasTCA = read.csv("Spreadsheets/MetAtlasTCA.tsv",
                          sep = "\t")
MetAtlasPPP = read.csv("Spreadsheets/MetAtlasPPP.tsv",
                       sep = "\t")
MetAtlasArgProMet = read.csv("Spreadsheets/MetAtlasArgProMet.tsv",
                       sep = "\t")
MetAtlasFattyAcidOxid = read.csv("Spreadsheets/MetAtlasFattyAcidOxid.tsv",
                             sep = "\t")


ExprPatternsByPathway = function(mat, Pathway, title) {
  
  Genes = Pathway$Genes
  Genes = unlist(str_split(Genes, ";")) #split each cell with delimiting semicolon, this creates a list, so need to unlist
  Genes = str_trim(Genes) #removes whitespace from start and end
  Genes = Genes[Genes != ""]
  Genes = sort(Genes)
  RowsInMat = rownames(mat)[grepl(paste(Genes, collapse = '|'), rownames(mat))] #grepl multiple patterns
  
  pheatmap(mat[RowsInMat,], cluster_rows = F, cluster_cols = F,
           show_rownames = T, 
           show_colnames = F, 
           annotation_col = annot,
           fontsize_row = 10,
           annotation_colors = ann_colors,
           main = title)
  
  
}

ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasGlycolysis, "Glycolysis")
ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasOxPhos, "OxPhos")
ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasTCA, "TCA cycle")
ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasPPP, "PPP")
ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasArgProMet, "Arginine/Proline Met")
ExprPatternsByPathway(scaled_test_norm_counts, MetAtlasFattyAcidOxid, "Fatty Acid Oxidation")


#summary stats: how many genes up and down for each tissue vs all tissues:

DegHeatmaps = function(x, Up=TRUE) {
  
  vec = list()
  
  for (i in names(Comparison)) {
    vec[[i]] = list()
  }
  
  for (i in 1:(length(Comparison)-1)) {
    
    for (j in 1:9) {
      
      if (Up) {
        vec[[j]][[i]] = length(which(CombinedRes4allComparisons[[j]][[i]]$log2FoldChange > 0))
        
      } 
        else {
          
        vec[[j]][[i]] = length(which(CombinedRes4allComparisons[[j]][[i]]$log2FoldChange < 0))
      }
      
      names(vec[[j]])[[i]] = Comparison[[j]][[i]]                                 
      
    }
  }
  
  vec = unlist(vec)
  
  vec = vec[grepl(paste0(x, "."), names(vec))]
  vec = as.matrix(vec)
  colnames(vec) = paste(x, "vs")
  
  if (colnames(vec) == "Monocytes vs") {
    colnames(vec) = sub("Monocytes vs", "Mono vs", colnames(vec))
  }
  
  rownames(vec) = sub(".*? ", "", rownames(vec))
  rownames(vec) = sub(".*? ", "", rownames(vec))
  rownames(vec) = sub("Monocytes", "Mono", rownames(vec))
  
  vec = as.matrix(vec)
  
   if (Up) {
     cols = colorRampPalette(c("white", "red"))(100)
   } else {
     cols = colorRampPalette(c("white", "blue"))(100)
   }

   p = pheatmap(t(vec),
            color = cols,
            cluster_rows = F, cluster_cols = F,
            show_colnames = T, show_rownames = T,
            fontsize_row = 12, fontsize_col = 10,
            display_numbers = T, number_color = "black",
            fontsize_number = 15, number_format = "%.0f",
            legend = F, angle_col = 0, cellwidth = 45, cellheight = 15)
  
   if (Up) {
     ggsave(paste0(x, "Up.pdf"), plot = p, path = "/Users/rokosango/PhD/TRM_macs/TRM_RNAseq/Plots/N_DEGs_Heatmaps",
           width = 6, height = 3, device = "pdf")
   } else {
     ggsave(paste0(x, "Down.pdf"), plot = p, path = "/Users/rokosango/PhD/TRM_macs/TRM_RNAseq/Plots/N_DEGs_Heatmaps",
           width = 6, height = 3, device = "pdf")
   }

}

vecUp = list()
vecDown = list()

for (i in unique(coldata$tissue)) {
  vecUp[[i]] = DegHeatmaps(i, Up = T)
}

for (i in unique(coldata$tissue)) {
  vecDown[[i]] = DegHeatmaps(i, Up = F)
}

vecUp = unlist(vecUp)
vecDown = unlist(vecDown)

MatUp = matrix(vecUp, nrow = 8, ncol = 9)
MatDown = matrix(vecDown, nrow = 8, ncol = 9)

BoxPlotDf = data.frame(rbind(MatUp, MatDown))
colnames(BoxPlotDf) = unique(coldata$tissue)
BoxPlotDf$Direction = rep(c("Up", "Down"), each = 8)
BoxPlotDf$Direction = as.factor(BoxPlotDf$Direction)
BoxPlotDf = melt(BoxPlotDf)

ggplot(BoxPlotDf, aes(x=variable, y=value, fill=Direction)) +
  geom_boxplot(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), dotsize = 0.8, binwidth = 100) +
  scale_fill_manual(values=c("#377EB8", "#E41A1C")) +
  labs(x = NULL, y = "Gene Count") + my.theme + coord_flip()


##########################################################################
##### PART 3. - GSEA on combined DEG results table, for each tissue. #####
##########################################################################

organism = "org.Mm.eg.db"

#######################
#### GENE ONTOLOGY ####
#######################


GOEnrichment = function(df) {
  
  original_gene_list <- df$log2FoldChange
  original_gene_list = as.numeric(original_gene_list)
  names(original_gene_list) <- df$Row.names
  gene_list<-na.omit(original_gene_list)
  
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               minGSSize = 3, 
               nPermSimple = 100000,
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE,
               eps = 0,
               OrgDb = organism, 
               pAdjustMethod = "fdr")
  
  gse.res = gse@result
  
  gse.res = gse.res %>% 
    mutate(GeneRatio = readr::parse_number(gse.res$leading_edge) / 100) %>% # parse_number great function for extracting numbers from strings 
    arrange(desc(NES)) %>%
    relocate(GeneRatio, .before = NES)
  
}

GOEnrichment2 = function (df) {
  
  GeneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys = as.character(df$Row.names),
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  dfUp = df %>%
    dplyr::filter(Row.names %in% GeneIDs$GENEID) %>%
    mutate(Symbol = GeneIDs$SYMBOL) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange > 0)
  
  dfDown = df %>%
    dplyr::filter(Row.names %in% GeneIDs$GENEID) %>%
    mutate(Symbol = GeneIDs$SYMBOL) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange < 0)
  
  setEnrichrSite("Enrichr")
  
  dbs <- listEnrichrDbs()
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
  
  enrichedUp = enrichr(dfUp$Symbol, dbs)
  
  enrichedUpBind = bind_rows(list(enrichedUp[["GO_Molecular_Function_2023"]],
                                  enrichedUp[["GO_Cellular_Component_2023"]],
                                  enrichedUp[["GO_Biological_Process_2023"]]), .id = "Term") #combining a list elements bind the same label
  
  enrichedUpBind$Direction = rep("Up", times = nrow(enrichedUpBind))
  
  enrichedUpBind$Term = as.character(c(enrichedUp[["GO_Molecular_Function_2023"]]$Term,
                                       enrichedUp[["GO_Cellular_Component_2023"]]$Term,
                                       enrichedUp[["GO_Biological_Process_2023"]]$Term))
  
  enrichedDown = enrichr(dfDown$Symbol, dbs)
  
  enrichedDownBind = bind_rows(list(enrichedDown[["GO_Molecular_Function_2023"]],
                                    enrichedDown[["GO_Cellular_Component_2023"]],
                                    enrichedDown[["GO_Biological_Process_2023"]]), .id = "Term")
  
  
  enrichedDownBind$Direction = rep("Down", times = nrow(enrichedDownBind))
  
  enrichedDownBind$Term = as.character(c(enrichedDown[["GO_Molecular_Function_2023"]]$Term,
                                         enrichedDown[["GO_Cellular_Component_2023"]]$Term,
                                         enrichedDown[["GO_Biological_Process_2023"]]$Term))
  
  
  enriched = as.data.frame(rbind(enrichedUpBind, enrichedDownBind))
  enriched  = enriched %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(desc(Combined.Score))
  
}

GO_list = list()
GO_list_2 = list()

for (i in names(Comparison)) {
  GO_list[[i]] = list()
  GO_list_2[[i]] = list()
}


for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    # GO_list[[j]][[i]] = GOEnrichment(CombinedRes4allComparisons[[j]][[i]])
    # names(GO_list[[j]])[[i]] = Comparison[[j]][[i]]
    
    GO_list_2[[j]][[i]] = GOEnrichment2(CombinedRes4allComparisons[[j]][[i]])
    names(GO_list_2[[j]])[[i]] = Comparison[[j]][[i]]
    
  }
}




#######################
#### KEGG PATHWAYS ####
#######################


KEGGEnrichment = function (df) {
  
  GeneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys = as.character(df$Row.names),
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  df = df %>%
    dplyr::filter(Row.names %in% GeneIDs$GENEID) %>%
    mutate(Symbol = GeneIDs$SYMBOL) %>%
    arrange(desc(log2FoldChange))
  
  original_gene_list <- df$log2FoldChange
  original_gene_list = as.numeric(original_gene_list)
  names(original_gene_list) = df$Symbol
  ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  names(original_gene_list) = ids$ENTREZID
  original_gene_list = original_gene_list[!duplicated(names(original_gene_list))]
  original_gene_list = sort(original_gene_list, decreasing = TRUE)
  
  kk2 <- gseKEGG(geneList     = original_gene_list,
                 organism     = "mmu",
                 #nPerm        = 10000,
                 nPermSimple = 100000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 eps = 0,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr", #"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                 keyType       = "ncbi-geneid")
  
  kegg.res = kk2@result
  
  kegg.res = kegg.res %>% 
    mutate(GeneRatio = readr::parse_number(kegg.res$leading_edge) / 100) %>% # parse_number great function for extracting numbers from strings 
    arrange(desc(NES)) %>%
    relocate(GeneRatio, .before = NES)
  
}

KEGGEnrichment2 = function (df) {
  
  GeneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys = as.character(df$Row.names),
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  dfUp = df %>%
    dplyr::filter(Row.names %in% GeneIDs$GENEID) %>%
    mutate(Symbol = GeneIDs$SYMBOL) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange > 0)
  
  dfDown = df %>%
    dplyr::filter(Row.names %in% GeneIDs$GENEID) %>%
    mutate(Symbol = GeneIDs$SYMBOL) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::filter(log2FoldChange < 0)
  
  setEnrichrSite("Enrichr")
  
  dbs <- listEnrichrDbs()
  dbs <- "KEGG_2019_Mouse"
  enrichedUp = enrichr(dfUp$Symbol, dbs)
  enrichedUp = as.data.frame(enrichedUp[["KEGG_2019_Mouse"]])
  enrichedUp$Direction = rep("Up", times = nrow(enrichedUp))
  
  enrichedDown = enrichr(dfDown$Symbol, dbs)
  enrichedDown = as.data.frame(enrichedDown[["KEGG_2019_Mouse"]])
  enrichedDown$Direction = rep("Down", times = nrow(enrichedDown))
  
  
  enriched = as.data.frame(rbind(enrichedUp, enrichedDown))
  enriched  = enriched %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(desc(Combined.Score))
  
}

KEGG_list = list()
KEGG_list_2 = list()

for (i in names(Comparison)) {
  #KEGG_list[[i]] = list()
  KEGG_list_2[[i]] = list()
  
}


for (i in 1:(length(Comparison)-1)) {
  for (j in 1:9) {
    
    KEGG_list[[j]][[i]] = KEGGEnrichment(CombinedRes4allComparisons[[j]][[i]])
    names(KEGG_list[[j]])[[i]] = Comparison[[j]][[i]]
    
    KEGG_list_2[[j]][[i]] = KEGGEnrichment2(CombinedRes4allComparisons[[j]][[i]])
    names(KEGG_list_2[[j]])[[i]] = Comparison[[j]][[i]]
    
    
  }
}

KEGG_list_targPathways = list()

for (i in names(Comparison)) {
  KEGG_list_targPathways[[i]] = list()
  
}

KEGG_list_targPathways = KEGG_list_targPathways[-2]
KEGG_list_2_no_BMDM = KEGG_list_2[-2]


for (i in 1:(length(Comparison)-1)) {
    
  ind = grepl("BMDM", names(KEGG_list_2_no_BMDM[[i]])) #logical: find "TRMi vs BMDM" in names
  
  KEGG_list_targPathways[[i]] = KEGG_list_2_no_BMDM[[i]][ind] #extract only "TRMi vs BMDM", disregard "TRMi vs TRMj"
  
  KEGG_list_targPathways[[i]] = KEGG_list_targPathways[[i]] %>%
    as.data.frame()
  
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], Term = ends_with("Term"))
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], Overlap = ends_with("Overlap"))
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], OldAdjustedPValue = ends_with("Old.Adjusted.P.value"))
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], AdjustedPValue = ends_with("Adjusted.P.value"))
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], CombinedScore = ends_with("Combined.Score"))
  KEGG_list_targPathways[[i]] = rename(KEGG_list_targPathways[[i]], Direction = ends_with("Direction"))
  
  KEGG_list_targPathways[[i]] = KEGG_list_targPathways[[i]] %>%
    dplyr::filter(Term %in% c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
                             "Pentose phosphate pathway", "Oxidative phosphorylation",
                             'Arginine and proline metabolism')) %>%
    dplyr::mutate(Population = rep(names(KEGG_list_targPathways)[i])) %>%
    dplyr::select(Term, Overlap, AdjustedPValue, CombinedScore, Direction, Population)

}

CombinedKeggFrame = as.data.frame(do.call(rbind, KEGG_list_targPathways))




###################################################
#### Common GO terms across tissue comparisons ####
###################################################

BarPlotsEnrichedTerms = function(x) {
  
  EnrichedTerms = list()
  
  for(i in 1:(length(Comparison)-1)) {                                    
    new_element <- GO_list_2[[x]][[i]]      # Create new list element
    EnrichedTerms[[length(EnrichedTerms) + 1]] <- new_element             # Append new list element
   
    EnrichedTerms[[i]] <- EnrichedTerms[[i]][order(EnrichedTerms[[i]]$Combined.Score, 
                                                         decreasing = TRUE), ]      
    EnrichedTerms[[i]] <- EnrichedTerms[[i]][!duplicated(EnrichedTerms[[i]]$Term), ] 
  }
  
  new_names = sub(".*? ", "", names(CombinedRes4allComparisons[[x]])) #remove everything before first space
  
  test_df = EnrichedTerms %>% 
    reduce(inner_join, by='Term') %>% #join all individual pairwise comparisons with a common ENSEMBL gene naames
    arrange(Term) %>%
    column_to_rownames(., var = "Term") %>%
    #as_tibble() %>%
    dplyr::select(starts_with("Direction")) %>%
    dplyr::select(1:8)
  
  names(test_df) = new_names
  
  test_df = test_df %>%
    mutate(UpReg = rowSums(test_df == "Up")) %>% #upregulated terms across comparisons
    mutate(DownReg = ncol(test_df) - UpReg) %>% #downregulated terms across comparisons
    mutate(Term = rownames(test_df)) %>%
    dplyr::select(UpReg, DownReg, Term) %>%
    mutate(Term = str_replace(Term, " \\s*\\([^\\)]+\\)", "")) %>%
    melt() %>%
    rename(Status = variable, Count = value) %>%
    ggplot(data=., aes(x=Term, y=Count, fill=Status)) +
    geom_bar(stat="identity")+
    coord_flip() +
    scale_fill_brewer(palette="Paired")+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
    labs(x = NULL, title = paste("GO Terms -", x, "vs all")) +
    theme_minimal() +
    theme(axis.text.y = element_text(colour="black", size=8, face = "bold"))
}

for (i in unique(coldata$tissue)) {
  print(BarPlotsEnrichedTerms(i))
}

# print(BarPlotsEnrichedTerms("AT"))
# print(BarPlotsEnrichedTerms("BMDM"))
# print(BarPlotsEnrichedTerms("Brain"))
# print(BarPlotsEnrichedTerms("Colon"))
# print(BarPlotsEnrichedTerms("Liver")) #not working, no common GO terms
# print(BarPlotsEnrichedTerms("Lung")) #not working, no common GO terms
# print(BarPlotsEnrichedTerms("Monocytes")) #not working, no common GO terms
# print(BarPlotsEnrichedTerms("Peritoneal")) #not working, no common GO terms
# print(BarPlotsEnrichedTerms("Spleen"))


########################################################
#### Common KEGG Pathways across tissue comparisons ####
########################################################

BarPlotsEnrichedPathways = function(x) {
  

EnrichedPathways = list()

for(i in 1:(length(Comparison)-1)) {                                    
  new_element <- KEGG_list_2[[x]][[i]]      # Create new list element
  EnrichedPathways[[length(EnrichedPathways) + 1]] <- new_element          # Append new list element
  
  EnrichedPathways[[i]] <- EnrichedPathways[[i]][order(EnrichedPathways[[i]]$Combined.Score, 
                                                       decreasing = TRUE), ]      # Order data
  EnrichedPathways[[i]] <- EnrichedPathways[[i]][!duplicated(EnrichedPathways[[i]]$Term), ]   # Unique rows of ordered data
  
}

new_names = sub(".*? ", "", names(CombinedRes4allComparisons[[x]])) #remove everything before first space

test_df = EnrichedPathways %>% 
  reduce(inner_join, by='Term') %>% #join all individual pairwise comparisons with a common ENSEMBL gene naames
  arrange(Term) %>%
  column_to_rownames(., var = "Term") %>%
  dplyr::select(starts_with("Direction")) %>%
  dplyr::select(1:8)

names(test_df) = new_names

test_df = test_df %>%
  mutate(UpReg = rowSums(test_df == "Up")) %>% #upregulated terms across comparisons
  mutate(DownReg = ncol(test_df) - UpReg) %>% #downregulated terms across comparisons
  mutate(Term = rownames(test_df)) %>%
  dplyr::select(UpReg, DownReg, Term) %>%
  mutate(Term = str_replace(Term, "AGE-RAGE signaling pathway in diabetic complications",
                            "AGE-RAGE signaling in diabetes")) %>%
  melt() %>%
  rename(Status = variable, Count = value) %>%
  ggplot(data=., aes(x=Term, y=Count, fill=Status)) +
  geom_bar(stat="identity")+
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 55)) +
  labs(x = NULL, title = paste("KEGG Pathways -", x, "vs all")) +
  theme_minimal() +1
  theme(axis.text.y = element_text(colour="black", size=8, face = "bold"))

}

for (i in unique(coldata$tissue)) {
  print(BarPlotsEnrichedPathways(i))
}

#############################################################
#### search genes for in ATM for Thomas Rauchenwald Graz ####
#############################################################

AtmGenesList = list()

for (i in 1:length(CombinedRes4allComparisons[["AT"]])) {
  
  AtmGenesList[[i]] = CombinedRes4allComparisons[["AT"]][[i]] %>%
    dplyr::filter(Symbol %in% c("Tnfrsf11a", "Tnfrsf11b", "Tnfsf11", "Pth1r")) %>%
    mutate(Comparison = rep(names(CombinedRes4allComparisons[["AT"]])[i]))
  
  names(AtmGenesList)[i] = names(CombinedRes4allComparisons[["AT"]])[i]
}


AtmGenesList = do.call( rbind, AtmGenesList)
write_xlsx(AtmGenesList, "AtmGenesList.xlsx")

