# Final project of 510

### Description 

RNA-seq Analysis of alcohol as risk factor of pancreas cancer in white male. 

### Overview of project 

Pancreas cancer as known as the king of cancer, it showed poor 5-year survival rate after diagnosed. The cause of pancreas cancer remains unknown, but I believed alcohol using is one of the most significant risk factors to cause pancreas cancer.

To test whether if the alcohol using is one of the risk factors, I compared two groups’ stage II pancreas cancer patients with alcohol using history and non-alcohol users in the white male. 

The result showed that majority of pancreas cancer patients in white male has alcohol using history, provided that alcohol using could be a risk factor in pancreas cancer. Meanwhile, the alcohol using will cause the most AKT1S1, MARK2P8, OCIAD1-AS1, RP11-43F13.3, LINC00632 gene mutation that may lead to pancreas cancer.

### Data 

From the national cancer institute GDC data portal.
Cancer: Pancreas cancer 
Program: TCGA
Project: TCGA-PAAD
Race: White
Sex: Male

### Section 1: overall result 
From picture 1,2,3 we could know that majority of white male patients of pancreas cancer have alcohol using history.

pie(c(length(clinical_exposure$condition[clinical_exposure$condition == FALSE]), length(clinical_exposure$condition[clinical_exposure$condition == TRUE])), c(non_label, label), main = "Overall")

Used these codes to get picture 1.

> barplot(by_race, main=paste0(label, " vs. ", non_label, " by Ethnicity"),
>           xlab="Ethnicity",col=c("blue", "red"), ylab="Population",
>           legend = paste0(label, ": ", rownames(by_race)))

Used these codes to get picture 2.

> barplot(by_gender, main=paste0(label, " vs. ", non_label, " by Gender"),
>           xlab="Gender",col=c("blue", "red"), ylab="Population",
>           legend = paste0(label, ": ", rownames(by_gender)))

Used these codes to get picture 3. 

### Section 2: PCA and Plot Counts
From picture 4 (PCA) we know data of alcohol user and non-alcohol are mixed with no clear boundary. But in picture 5 (plot counts), they have obvious difference.

> Vsd <- vst(dds, blind=FALSE)
> plotPCA(vsd, intgroup=c(“condition”), returnData=TRUE)

Used these codes to get picture 4

> ggplot(d, aes(x=condition, y=count)) +
> geom_point(position=position_jitter(w=0.1,h=0)) +
> scale_y_log10(breaks=c(25,100,400))+ ggtitle(“Plot counts”)+
> theme(plot.title = element_text(hjust = 0.5, face = “bold”))

Used these codes to get picture 5

### Section 3: Exploring and Exporting result 
From picture 6 the raw graph and picture 7 the analyzed picture we could find out that the alcohol patients will change a lot of different genes, and there are no specific genes that would mutate due by the alcohol using. Both non-alcohol and alcohol patients of pancreas cancer will have high mutate on in different genes.
	
> Res <- results(dds, contrast=c(“condition”,non_label,label))
> plotMA(res, ylim=c(-2,2), main=paste0(“Compare between “, label, “ and “, non_label))

Used these codes to get picture 6. 

> resLFC <- lfcShrink(dds, coef=paste0(“condition_”, label, “_vs_”, non_label), type=”apeglm”)
> plotMA(resLFC, ylim=c(-2,2), main=paste0(“Compare between “, label, “ and “, non_label, “ (LFC Shrink)”))

Used these codes to get picture 7. 

### Section 4: Alternative shrinkage estimators
From picture 8, the gene mutation between non-alcohol and alcohol of pancreas cancer patients do not have much difference, while the picture 10 showed us that the specific genes of non-alcohol and alcohol patients. The specific genes of alcohol patients showed wild dispersed while the non-alcohol patients showed similar.

> resLFC <- lfcShrink(dds, coef=paste0(“condition_”, label, “_vs_”, non_label), type=”apeglm”)
> plotMA(resLFC, xlim=xlim, ylim=ylim, main=”apeglm”)

Used these codes to get picture 8.

> resNorm <- lfcShrink(dds, coef=2, type=”normal”)
> plotMA(resNorm, xlim=xlim, ylim=ylim, main=”normal”)

Used these codes to get picture 9.

> resAsh <- lfcShrink(dds, coef=2, type=”ashr”)
> plotMA(resAsh, xlim=xlim, ylim=ylim, main=”ashr”)

Used these codes to get picture 10

### Section 5: Significance in alcohol/non-alcohol comparison 
From picture 11 it showed us the probability of different gene mutation. The darker colors have higher probability occur in the patients, however, compare with picture 6 most of the high chance mutation came from the non-alcohol patients not alcohol patients. 

> Ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) +
> geom_point(size=1) +
> scale_y_continuous(limits=c(-3, 3), oob=squish) +
> scale_x_log10() +
> geom_hline(yintercept = 0, colour=”darkorchid4”, size=1, linetype=”longdash”) +
> labs(x=”mean of normalized counts”, y=”log fold change”) +
> scale_colour_viridis(direction=-1, trans=’sqrt’) +
> theme_bw() + geom_density_2d(colour=”black”, size=2)+
> ggtitle(paste0(“Significance in “,label,”/”, non_label, “ comparison”))+
> theme(plot.title = element_text(hjust = 0.5, face = “bold”))

Used these codes to get picture 11

### Section 6: Effects of transformations on the variance
From picture 12 and 13, we could notify that from rank 20,000 to 30,000 the gene mutation occurs more frequently due to the high sd. The overall of two pictures showed that the pancreas cancer patients have high genes mutations. 

> Ntd <- normTransform(dds)
> meanSdPlot(assay(ntd))$gg+ggtitle(“log2(n + 1) vs. mean”)+
> theme(plot.title = element_text(hjust = 0.5, face = “bold”))

Used these codes to get picture 12.

> Deseq2VST <- assay(vsd)
> meanSdPlot(deseq2VST)$gg+ggtitle(“Effects of transformations on the variance”)+
> theme(plot.title = element_text(hjust = 0.5, face = “bold”))

Used these codes to get picture 13.

### Section 7: Extended section on shrinkage estimators
From picture 14, it also showed that the alcohol patients have a lot of genes mutation, and the mutations are not specific.

> resApeT <- lfcShrink(dds, coef=2, type=”apeglm”, lfcThreshold=1)
> plotMA(resApeT, ylim=c(-3,3), cex=.8, main=”Extended section on shrinkage estimators”)

Used these codes to get picture 14

### Section 8: Tests of log2 fold change above or below a threshold
From picture 15 – 18, the four pictures of Wald tests, the specific genes from greater (abs), greater, less are limited. Most specific genes that are selected are from less (abs) range from 1e+00 – 1e+02.

> resGA <- results(dds, lfcThreshold=.5, altHypothesis=”greaterAbs”)
> plotMA(resGA, ylim=ylim, main=”Greater (Abs)”); drawLines()

Used these codes to get picture 15.

> resLA <- results(dds, lfcThreshold=.5, altHypothesis=”lessAbs”)
> plotMA(resLA, ylim=ylim, main=”Less (Abs)”); drawLines()

Used these codes to get picture 16.

> resG <- results(dds, lfcThreshold=.5, altHypothesis=”greater”)
> plotMA(resG, ylim=ylim, main=”Greater”); drawLines()

Used these codes to get picture 17.

> resL <- results(dds, lfcThreshold=.5, altHypothesis=”less”)
> plotMA(resL, ylim=ylim, main=”Less”); drawLines()

Used these codes to get picture 18.

### Section 9: Count outlier detection
From picture 19, we could know after using Cook’s distance, the distance per gene is small, and most of the gene do not have high mutation chance. From picture 20, we could notify the outliner from picture 19. The gene mutation frequency is almost equal except 0 group.

> Plot(res$baseMean+1, -log10(res$pvalue),
>        log=”x”, xlab=”mean of normalized counts”,
>        ylab=expression(-log[10](pvalue)),
>        ylim=c(0,30),
>        cex=.4, col=rgb(0,0,0,.3), title(“Outlier detection (if any)”))

Used these codes to get picture 19

> use <- res$baseMean >  metadata(res)$filterThreshold
> h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
> h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
> colori <- c(`do not pass`=”khaki”, `pass`=”powderblue”)
> 
> barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
>           col = colori, space = 0, main = “p value histogram”, ylab=”frequency”)
> text(x = c(0, length(h1$counts)), y = 0, label = paste0(c(0,1)),
>        adj = c(0.5,1.7), xpd=NA)
> legend(“topright”, fill=rev(colori), legend=rev(names(colori)))

Used these codes to get picture 20

### Section 10: Independent filtering of results
From picture 21, we could know this picture maximizes the number of rejections over the quantiles of a filter statistic.

> metadata(res)$alpha
> metadata(res)$filterThreshold
> 
> plot(metadata(res)$filterNumRej,
>        type="b", ylab="number of rejections",
>        xlab="quantiles of filter", main="Independent filtering of results")
> lines(metadata(res)$lo.fit, col="red")
> abline(v=metadata(res)$filterTheta)

Used these codes to get picture 21

### Section 11: Dispersion plot and fitting alternatives
From picture 22, we could find that a lot of final gene-est are not shrunk towards the fitted value, which means there are lot of outliers. 

> plotDispEsts(dds, main="Dispersion plot and fitting alternatives")

Used these codes to get picture 22

### Section 12: Approach to count outliers
From picture 23, the picture showed us after using Cook’s distances which gene has higher range than others. For all the genes, the minimum requirement is third times replications.

> boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main="count outliers")

Used these codes to get picture 23

### Section 13: Heatmap of genes of each case 
From picture 24, we could know from each case, which gene has higher mutation. The higher mutation rate, the color is brighter. 

> grid.draw(arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5)))

Used these codes to get picture 24

### Section 14: Mutated gene overlapping 

Using GSEA tested the top 5 mutation risk of gene, AKT1S1, MARK2P8, OCIAD1-AS1, RP11-43F13.3, LINC00632, and there is no overlapping between those five genes. Provided that alcohol using will cause random mutation that leads to pancreas cancer. 


# Conclusion

In the white male, the result showed that alcohol using patients has higher chance to get pancreas cancer than the non-alcohol users. Provided that the alcohol using could be one of the risk factor of pancreas cancer. 

Compare with alcohol using and non-alcohol patients, the most mutated gene including: AKT1S1, MARK2P8, OCIAD1-AS1, RP11-43F13.3, LINC00632 showed that those gene mutation could be a cause pf pancreas cancer and be a gene target for future treatment such as CRISPR-Cas9.  

However, for those non-alcohol patients, there are still a lot of unknowns that need further research. 
