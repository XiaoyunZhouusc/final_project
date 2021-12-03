if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!require("ashr", quietly = TRUE))
  BiocManager::install("ashr")

if (!require("apeglm", quietly = TRUE))
  BiocManager::install("apeglm")

if (!require("IHW", quietly = TRUE))
  BiocManager::install("IHW")

if (!require("vsn", quietly = TRUE))
  BiocManager::install("vsn")

if (!require("ggdendro", quietly = TRUE))
  install.packages("ggdendro")

if (!require("gtable", quietly = TRUE))
  install.packages("gtable")

if (!require("ff", quietly = TRUE))
  install.packages("ff")

library(gtable)
library(grid)
library(ggdendro)
library(gridExtra)
library(png)
library(DESeq2)
library(tidyr)

library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

library(IHW)

library(vsn)
library(pheatmap)
library(RColorBrewer)

library(reshape2)

library(ff)

#library(BiocParallel)
#register(MulticoreParam(4))

differential_expression_analysis <- function(label, exposure, clinical, htseq_directory_of_covariate, list_to_include){

  dir.create("images", showWarnings = FALSE)

  unlink(file.path("images", label), recursive = TRUE)
  dir.create(file.path("images", label))

  non_label = paste0("non",".", label)

  clinical_exposure <- exposure(read.table("data/clinical/exposure.tsv", sep = '\t', header = TRUE,  quote = ""))

  png(filename = file.path("images", label, "PIE_CHART.png"))
  pie(c(length(clinical_exposure$condition[clinical_exposure$condition == FALSE]),
        length(clinical_exposure$condition[clinical_exposure$condition == TRUE])),
      c(non_label, label), main = "Overall"
  )
  dev.off()

  PIE_CHART <- readPNG(file.path("images", label, "PIE_CHART.png"))

  grid.arrange(rasterGrob(PIE_CHART),ncol=1)

  clinical_clinical <- read.table("data/clinical/clinical.tsv", sep = '\t', header = TRUE,  quote = "")

  odd <- seq(1,length(clinical_clinical$case_id),2)

  clinical_clinical <- clinical_clinical[odd,] # extract odd row. Somehow '' has duplicates

  clinical_clinical<-clinical(clinical_clinical)

  merged_clinical <- merge(clinical_exposure, clinical_clinical)

  by_race <- table(merged_clinical$condition, merged_clinical$race)

  png(filename = file.path("images", label, "Condition_vs_Non_Condition_by_Ethnicity.png"))
  barplot(by_race, main=paste0(label, " vs. ", non_label, " by Ethnicity"),
          xlab="Ethnicity",col=c("blue", "red"), ylab="Population",
          legend = paste0(label, ": ", rownames(by_race)))

  dev.off()

  bar_race <- readPNG(file.path("images", label, "Condition_vs_Non_Condition_by_Ethnicity.png"))

  by_gender <- table(merged_clinical$condition, merged_clinical$gender)

  png(filename = file.path("images", label, "Condition_vs_Non_Condition_by_Gender.png"))
  barplot(by_gender, main=paste0(label, " vs. ", non_label, " by Gender"),
          xlab="Gender",col=c("blue", "red"), ylab="Population",
          legend = paste0(label, ": ", rownames(by_gender)))

  dev.off()

  bar_gender <- readPNG(file.path("images", label, "Condition_vs_Non_Condition_by_Gender.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),ncol=3)


  directory <- "data/htseq.counts"

  htseq_table <- read.table(file.path(directory, "MANIFEST.txt"), header=TRUE)

  target_cases <- read.table(htseq_directory_of_covariate, header = TRUE)

  include_cases <- read.table(list_to_include, header = TRUE)

  htseq_table<-htseq_table[htseq_table$id != "\\N", ]

  htseq_table<-htseq_table[htseq_table$id %in% include_cases$id, ]

  htseq_table$condition = non_label

  htseq_table[htseq_table$id %in% target_cases$id, ]$condition <- label

  dds <- DESeqDataSetFromHTSeqCount(sampleTable = htseq_table,
                                         directory = directory,
                                         design= ~ condition)

  # Pre-filtering

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  # Differential expression analysis

  # set factor level, so no need to use "contrast" later
  dds$condition <- factor(dds$condition, levels = c(non_label, label))
  #dds$condition <- droplevels(dds$condition)
  dds <- DESeq(dds)
  #res <- results(dds)
  #res <- results(dds, name="condition_non.smoker_vs_smoker")
  res <- results(dds, contrast=c("condition",non_label,label))

  res

  # Log fold change shrinkage for visualization and ranking

  resultsNames(dds)

  resLFC <- lfcShrink(dds, coef=paste0("condition_", label, "_vs_", non_label), type="apeglm")
  resLFC

  # p-values and adjusted p-values

  resOrdered <- res[order(res$pvalue),]

  summary(res)

  sum(res$padj < 0.1, na.rm=TRUE)

  res05 <- results(dds, alpha=0.05)
  summary(res05)

  sum(res05$padj < 0.05, na.rm=TRUE)

  # Independent hypothesis weighting

  resIHW <- results(dds, filterFun=ihw)
  summary(resIHW)
  sum(resIHW$padj < 0.1, na.rm=TRUE)
  metadata(resIHW)$ihwResult

  # Exploring and exporting results

  png(filename = file.path("images", label, "Compare_between_Condition_and_non_Condition_LFC_Shrink.png"))
  plotMA(resLFC, ylim=c(-2,2), main=paste0("Compare between ", label, " and ", non_label, " (LFC Shrink)"))
  dev.off()

  png(filename = file.path("images", label, "Compare_between_Condition_and_non_Condition.png"))
  plotMA(res, ylim=c(-2,2), main=paste0("Compare between ", label, " and ", non_label))
  dev.off()

  img2 <- readPNG(file.path("images", label, "Compare_between_Condition_and_non_Condition_LFC_Shrink.png"))

  img3 <- readPNG(file.path("images", label, "Compare_between_Condition_and_non_Condition.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3), nrow = 2, ncol=4)

  #idx <- identify(res$baseMean, res$log2FoldChange)
  #rownames(res)[idx]

  # Alternative shrinkage estimators

  resultsNames(dds)

  resNorm <- lfcShrink(dds, coef=2, type="normal")
  resAsh <- lfcShrink(dds, coef=2, type="ashr")

  xlim <- c(1,1e5); ylim <- c(-3,3)

  png(filename = file.path("images", label, "apeglm.png"))
  plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
  dev.off()

  png(filename = file.path("images", label, "normal.png"))
  plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
  dev.off()

  png(filename = file.path("images", label, "ashr.png"))
  plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
  dev.off()

  img4 <- readPNG(file.path("images", label, "apeglm.png"))

  img5 <- readPNG(file.path("images", label, "normal.png"))

  img6 <- readPNG(file.path("images", label, "ashr.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6), nrow = 2, ncol=4)

  # Plot counts

  #plotCounts(dds, gene=which.min(res$padj), intgroup="condition", main = "Plot counts")

  d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
                  returnData=TRUE)

  ggplot(d, aes(x=condition, y=count)) +
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10(breaks=c(25,100,400))+ ggtitle("Plot counts")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave("PlotCounts.png", path = file.path("images", label), device = "png")

  img7 = readPNG(file.path("images", label, "PlotCounts.png"))

  # Coerce to a data frame
  deseq2ResDF <- as.data.frame(res)

  # Examine this data frame
  head(deseq2ResDF)

  # Set a boolean column for significance
  deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

  # Plot the results similar to DEseq2
  #ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

  # Let's add some more detail

  ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) +
    geom_point(size=1) +
    scale_y_continuous(limits=c(-3, 3), oob=squish) +
    scale_x_log10() +
    geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
    labs(x="mean of normalized counts", y="log fold change") +
    scale_colour_viridis(direction=-1, trans='sqrt') +
    theme_bw() + geom_density_2d(colour="black", size=2)+
    ggtitle(paste0("Significance in ",label,"/", non_label, " comparison"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave("Significance.png", path = file.path("images", label), device = "png")

  img8 <- readPNG(file.path("images", label, "Significance.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender), rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7), nrow = 3, ncol=4)

  # More information on results columns

  mcols(res)$description

  # Exporting results to CSV files

  write.csv(as.data.frame(resOrdered),
            file=file.path("images", label, paste0("condition_", label, "_vs_", non_label, ".csv")))

  resSig <- subset(resOrdered, padj < 0.1)
  resSig

  # Effects of transformations on the variance
  vsd <- vst(dds, blind=FALSE)
  #rld <- rlog(dds, blind=FALSE)
  deseq2VST <- assay(vsd)
  head(deseq2VST, 3)
  # this gives log2(n + 1)
  ntd <- normTransform(dds)

  png(filename = file.path("images", label, "meanSdPlot_btd.png"))
  meanSdPlot(assay(ntd))$gg+ggtitle("log2(n + 1) vs. mean")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  dev.off()

  img9 <- readPNG(file.path("images", label, "meanSdPlot_btd.png"))

  png(filename = file.path("images", label, "meanSdPlot_vst.png"))
  meanSdPlot(deseq2VST)$gg+ggtitle("Effects of transformations on the variance")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  dev.off()

  img10 <- readPNG(file.path("images", label, "meanSdPlot_vst.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),rasterGrob(img9),rasterGrob(img10), nrow = 3, ncol=4)


  # Convert the DESeq transformed object to a data frame

  deseq2VST <- as.data.frame(deseq2VST)
  deseq2VST$Gene <- rownames(deseq2VST)
  head(deseq2VST)

  # Keep only the significantly differentiated genes where the fold-change was at least 3
  sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
  deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

  # First compare wide vs long version
  deseq2VST_wide <- deseq2VST
  deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

  head(deseq2VST_wide)
  head(deseq2VST_long)

  # Now overwrite our original data frame with the long format
  deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

  # Make a heatmap
  heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #heatmap

  # Convert the significant genes back to a matrix for clustering
  deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
  rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
  deseq2VSTMatrix$Gene <- NULL

  # Compute a distance calculation on both dimensions of the matrix
  distanceGene <- dist(deseq2VSTMatrix)
  distanceSample <- dist(t(deseq2VSTMatrix))

  # Cluster based on the distance calculations
  clusterGene <- hclust(distanceGene, method="average")
  clusterSample <- hclust(distanceSample, method="average")

  sampleModel <- as.dendrogram(clusterSample)
  sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
  sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()+ggtitle("Heatmap of genes of each case")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Re-factor samples for ggplot2
  deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

  # Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
  heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #heatmap

  # Modify the ggplot objects
  sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
  heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

  # Convert both grid based objects to grobs
  sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
  heatmapGrob <- ggplotGrob(heatmap_1)

  # Check the widths of each grob
  sampleDendrogramGrob$widths
  heatmapGrob$widths

  # Add in the missing columns
  sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
  sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

  # Make sure every width between the two grobs is the same
  maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
  sampleDendrogramGrob$widths <- as.list(maxWidth)
  heatmapGrob$widths <- as.list(maxWidth)

  # Arrange the grobs into a plot
  finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

  # Draw the plot
  png(filename = file.path("images", label, "heatmap.png"))
  grid.draw(finalGrob)
  dev.off()

  img11 <- readPNG(file.path("images", label, "heatmap.png"))

  # Principal component plot of the samples

  #plotPCA(vsd, intgroup=c("condition"))

  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()+ggtitle("PCA")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave("PCA.png", path = file.path("images", label), device = "png")

  img12 <- readPNG(file.path("images", label, "PCA.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12), nrow = 4, ncol=4)

  # Wald test individual steps

  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)

  # Contrasts

  #results(dds, contrast=c("condition","C","B"))

  # Interactions

  dds$group <- factor(paste0(dds$genotype, dds$condition))
  design(dds) <- ~ group
  dds <- DESeq(dds)
  resultsNames(dds)
  results(dds, contrast=c("group", label, non_label))

  # Likelihood ratio test

  #dds <- DESeq(dds, test="LRT", reduced=~1)
  #res <- results(dds)

  #dds <- DESeq(dds, test="LRT", reduced=~batch)
  #res <- results(dds)

  # Extended section on shrinkage estimators

  resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)

  png(filename = file.path("images", label, "shrinkage_estimators.png"))
  plotMA(resApeT, ylim=c(-3,3), cex=.8, main="Extended section on shrinkage estimators")
  abline(h=c(-1,1), col="dodgerblue", lwd=2)
  dev.off()

  img13 <- readPNG(file.path("images", label, "shrinkage_estimators.png"))

  condition <- factor(rep(c("A","B","C"),each=2))
  model.matrix(~ condition)

  # to compare C vs B, make B the reference level,
  # and select the last coefficient
  condition <- relevel(condition, "B")
  model.matrix(~ condition)

  grp <- factor(rep(1:3,each=4))
  cnd <- factor(rep(rep(c("A","B"),each=2),3))
  model.matrix(~ grp + cnd + grp:cnd)

  # to compare condition effect in group 3 vs 2,
  # make group 2 the reference level,
  # and select the last coefficient
  grp <- relevel(grp, "2")
  model.matrix(~ grp + cnd + grp:cnd)

  grp <- factor(rep(1:2,each=4))
  ind <- factor(rep(rep(1:2,each=2),2))
  cnd <- factor(rep(c("A","B"),4))
  model.matrix(~grp + grp:ind + grp:cnd)

  # to compare condition effect across group,
  # add a main effect for 'cnd',
  # and select the last coefficient
  model.matrix(~grp + cnd + grp:ind + grp:cnd)

  # Approach to count outliers

  png(filename = file.path("images", label, "count_outliers.png"))
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main="count outliers")
  dev.off()

  img14 <- readPNG(file.path("images", label, "count_outliers.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12),rasterGrob(img13),rasterGrob(img14), nrow = 4, ncol=4)

  # Dispersion plot and fitting alternatives
  png(filename = file.path("images", label, "Dispersion_plot_and_fitting_alternatives.png"))
  plotDispEsts(dds, main="Dispersion plot and fitting alternatives")
  dev.off()

  img15 <- readPNG(file.path("images", label, "Dispersion_plot_and_fitting_alternatives.png"))

  # Supply a custom dispersion fit

  ddsCustom <- dds
  useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
  medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                       na.rm=TRUE)
  dispersionFunction(ddsCustom) <- function(mu) medianDisp
  ddsCustom <- estimateDispersionsMAP(ddsCustom)

  # Independent filtering of results

  metadata(res)$alpha
  metadata(res)$filterThreshold

  png(filename = file.path("images", label, "quantiles_of_filter.png"))
  plot(metadata(res)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter", main="Independent filtering of results")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()

  img16 <- readPNG(file.path("images", label, "quantiles_of_filter.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12),rasterGrob(img13),rasterGrob(img14),rasterGrob(img15),rasterGrob(img16), nrow = 5, ncol=4)


  resNoFilt <- results(dds, independentFiltering=FALSE)
  addmargins(table(filtering=(res$padj < .1),
                   noFiltering=(resNoFilt$padj < .1)))

  # Tests of log2 fold change above or below a threshold

  par(mfrow=c(2,2),mar=c(2,2,1,1))
  ylim <- c(-2.5,2.5)
  resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
  resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
  resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
  resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
  drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)

  png(filename = file.path("images", label, "GreaterAbs.png"))
  plotMA(resGA, ylim=ylim, main="Greater (Abs)"); drawLines()
  dev.off()

  png(filename = file.path("images", label, "GreaterLess.png"))
  plotMA(resLA, ylim=ylim, main="Less (Abs)"); drawLines()
  dev.off()

  png(filename = file.path("images", label, "Greater.png"))
  plotMA(resG, ylim=ylim, main="Greater"); drawLines()
  dev.off()

  png(filename = file.path("images", label, "Less.png"))
  plotMA(resL, ylim=ylim, main="Less"); drawLines()
  dev.off()

  img17 <- readPNG(file.path("images", label, "GreaterAbs.png"))
  img18 <- readPNG(file.path("images", label, "GreaterLess.png"))
  img19 <- readPNG(file.path("images", label, "Greater.png"))
  img20 <- readPNG(file.path("images", label, "Less.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),
               rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),
               rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),
               rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12),
               rasterGrob(img13),rasterGrob(img14),
               rasterGrob(img17),rasterGrob(img18),rasterGrob(img19),rasterGrob(img20),
               rasterGrob(img15),rasterGrob(img16), nrow = 6, ncol=4)


  # Access to all calculated values

  mcols(dds,use.names=TRUE)[1:4,1:4]
  substr(names(mcols(dds)),1,10)
  mcols(mcols(dds), use.names=TRUE)[1:4,]
  head(assays(dds)[["mu"]])
  head(assays(dds)[["cooks"]])
  head(dispersions(dds))
  head(mcols(dds)$dispersion)
  sizeFactors(dds)
  head(coef(dds))
  attr(dds, "betaPriorVar")
  #priorInfo(resLFC)
  priorInfo(resNorm)
  priorInfo(resAsh)
  dispersionFunction(dds)

  attr(dispersionFunction(dds), "dispPriorVar")
  metadata(dds)[["version"]]

  # Sample-/gene-dependent normalization factors

  coldata <- DataFrame(grp=factor(rep(c("X","Y"),each=6)),
                       ind=factor(rep(1:6,each=2)),
                       cnd=factor(rep(c("A","B"),6)))
  coldata

  as.data.frame(coldata)

  coldata$ind.n <- factor(rep(rep(1:3,each=2),2))
  as.data.frame(coldata)

  model.matrix(~ grp + grp:ind.n + grp:cnd, coldata)

  #results(dds, contrast=list("grpY.cndB","grpX.cndB"))

  # Levels without samples

  group <- factor(rep(1:3,each=6))
  condition <- factor(rep(rep(c("A","B","C"),each=2),3))
  d <- DataFrame(group, condition)[-c(17,18),]
  as.data.frame(d)

  m1 <- model.matrix(~ condition*group, d)
  colnames(m1)

  unname(m1)

  all.zero <- apply(m1, 2, function(x) all(x==0))
  all.zero

  idx <- which(all.zero)
  m1 <- m1[,-idx]
  unname(m1)

  # Count outlier detection

  W <- res$stat
  maxCooks <- apply(assays(dds)[["cooks"]],1,max)
  idx <- !is.na(W)
  #plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
  #     ylab="maximum Cook's distance per gene",
  #     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
  m <- ncol(dds)
  p <- 3

  # Filtering criteria

  png(filename = file.path("images", label, "Outlier_detection.png"))
  plot(res$baseMean+1, -log10(res$pvalue),
       log="x", xlab="mean of normalized counts",
       ylab=expression(-log[10](pvalue)),
       ylim=c(0,30),
       cex=.4, col=rgb(0,0,0,.3), title("Outlier detection (if any)"))

  abline(h=qf(.99, p, m - p))
  dev.off()

  img21 <- readPNG(file.path("images", label, "Outlier_detection.png"))

  #  p value histogram

  use <- res$baseMean > metadata(res)$filterThreshold
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")

  png(filename = file.path("images", label, "bar_plot.png"))
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "p value histogram", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste0(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  dev.off()

  img22 <- readPNG(file.path("images", label, "bar_plot.png"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),
               rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),
               rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),
               rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12),
               rasterGrob(img13),rasterGrob(img14),
               rasterGrob(img17),rasterGrob(img18),rasterGrob(img19),rasterGrob(img20),rasterGrob(img15),rasterGrob(img16),rasterGrob(img21),rasterGrob(img22), nrow = 6, ncol=4)

  pdf(file = file.path("images", label, "all.pdf"))

  grid.arrange(rasterGrob(PIE_CHART),rasterGrob(bar_race),rasterGrob(bar_gender),
               rasterGrob(img2), rasterGrob(img3),rasterGrob(img4),
               rasterGrob(img5),rasterGrob(img6),rasterGrob(img8),rasterGrob(img7),
               rasterGrob(img9),rasterGrob(img10),rasterGrob(img11),rasterGrob(img12),
               rasterGrob(img13),rasterGrob(img14),
               rasterGrob(img17),rasterGrob(img18),rasterGrob(img19),rasterGrob(img20),rasterGrob(img15),rasterGrob(img16),rasterGrob(img21),rasterGrob(img22), nrow = 6, ncol=4)

  dev.off()

  dir.create("output", showWarnings = FALSE)

  unlink(file.path("output", label), recursive = TRUE)

  file.move(file.path("images", label),file.path("output", label))
}

differential_expression_analysis("alcohol", function(clinical_exposure) {
  clinical_exposure <- clinical_exposure[clinical_exposure$alcohol_history != "Not Reported",]
  clinical_exposure$condition <- FALSE
  clinical_exposure[clinical_exposure$alcohol_history == "Yes",]$condition <- TRUE
  clinical_exposure
}, function(clinical) {
  clinical<-clinical[clinical$gender == "male",]
  clinical<-clinical[clinical$race == "white",]
  clinical
}, "static/alcohol_white_man.txt", "static/white_man.txt")
