#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#setwd(system("pwd"))

library(DESeq2)
library(tidyr)

directory <- "data/htseq.counts"

htseq_table <- read.table(file.path(directory, "MANIFEST.txt"), header=TRUE)

htseq_table<-htseq_table[htseq_table$id != "\\N", ]

#htseq_table$case_id <- str_extract(htseq_table$filename, "(?=.*)/.*(?=.htseq.counts.gz)")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = htseq_table,
                                       directory = directory,
                                       design= ~ 1)

clinical_exposure <- read.table("data/clinical/exposure.tsv", sep = '\t', header = TRUE,  quote = "")

clinical_exposure_years_smoked <- clinical_exposure$years_smoked

clinical_exposure_years_smoked[clinical_exposure_years_smoked <= 0] <- 0

clinical_exposure_years_smoked <- as.numeric(clinical_exposure_years_smoked)

pie(c(length(clinical_exposure_years_smoked[clinical_exposure_years_smoked == 0]),
      length(clinical_exposure_years_smoked[clinical_exposure_years_smoked > 0])),
    c("non-smoker", "smoker")
    )

#library(ggplot2)
#library(data.table)

#clinical<- read.table("data/clinical/clinical.tsv", sep = '\t', header = TRUE,  quote = "")

#cig_per_day_vs_days_to_live <- merge(clinical_exposure,
#                                     clinical, by="case_id")

#cig_per_day_vs_days_to_live <- merge(htseq_table,
#                                     cig_per_day_vs_days_to_live, by="case_id")





#cig_per_day_vs_days_to_live <- data.table(cig_per_day_vs_days_to_live$cigarettes_per_day, cig_per_day_vs_days_to_live$days_to_death)

#cig_per_day_vs_days_to_live$V1[cig_per_day_vs_days_to_live$V1 == "'--"] <- 0

#cig_per_day_vs_days_to_live$V1 <- as.numeric(cig_per_day_vs_days_to_live$V1)

#cig_per_day_vs_days_to_live[cig_per_day_vs_days_to_live$V2 <= 0] <- 0

#cig_per_day_vs_days_to_live$V2 <- as.numeric(cig_per_day_vs_days_to_live$V2)

#cig_per_day_vs_days_to_live[cig_per_day_vs_days_to_live$V1 > 0 ]


