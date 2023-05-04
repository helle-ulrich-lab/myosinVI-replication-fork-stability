#load packages
library(plyr)
library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(readr)
library(plotly)
library(gplots)
library(igraph)
library(RColorBrewer)
library(GGally)
library(viridisLite)
library(viridis)
library(pheatmap)

# check loaded libraries
(.packages())

#clear environment
rm(list = ls())

# set working directory, file and exp names
setwd("//insert/directory/address/here")

file_name <- "proteinGroups.txt"
hp_proteins <- "Ubdomain.txt"

data_dir <- file.path('.')

exp_names <- c("Exp_1", "Exp_2","Exp_5", "Exp_6")

#load data, filter
sites <- read.delim(file.path(data_dir, file_name),
                    stringsAsFactors=FALSE)

chosen_ones <- read.delim(file.path(data_dir, hp_proteins),
                          stringsAsFactors=FALSE)

no_reverse <- sites[,'Reverse'] != '+'
no_contaminant <- sites[,'Potential.contaminant'] != '+'
not_only_identified_bysite <- sites[,"Only.identified.by.site"] !="+"
unique_peptides<- sites[,"Unique.peptides"] >= 1
AScore <- sites[,'Score'] >=0

sites <-
  subset(sites,
         no_reverse & no_contaminant & not_only_identified_bysite & unique_peptides & AScore)


#set ratio groups
RatioCols_ML <- c("Ratio.M.L.Exp_1", "Ratio.M.L.Exp_2", "Ratio.M.L.Exp_5", "Ratio.M.L.Exp_6")
sites <-
  cbind(sites,
        log2=log2(sites[,RatioCols_ML]))

#define limma function
significance_limma = function(ratios, max_q = 0.05, min_count = 2, prefix = '') {
  result = data.frame(
    count = ncol(ratios) - rowSums(is.na(ratios)),
    mean = rowMeans(ratios, na.rm = T)
  )
  fit = as.data.frame(eBayes(lmFit(ratios[result$count >= min_count, ])))
  result[result$count >= min_count, 'p'] = fit$p.value
  result[result$count >= min_count, 'q'] = p.adjust(fit$p.value, method = 'fdr')
  results = mutate(
    result,
    up = count >= min_count & mean > 2 & q <= max_q,
    down = count >= min_count & mean < 1/2 & q <= max_q
  )
  names(results) = paste0(prefix, names(results))
  return(results)
}

#calculate limma statistics
ml_statistics = significance_limma(
  log2(sites[,RatioCols_ML]))
colnames(ml_statistics) = paste0('Ratio.M.L.normalized.limma.', colnames(ml_statistics))
sites = bind_cols(sites, ml_statistics)

#filter names column
sites$names <- gsub("^(.*?)\\;.*$", "\\1", sites$Gene.names)

#write a table with limma statistics
write.table(sites, file = "limma-calculation.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

##ggplot function (volcano plot)
#define fold change and FDR cutoffs
x1 = 2
y1 = 0.05

upreg_ML <- sites[,"Ratio.M.L.normalized.limma.mean"] > log2(x1) & 
  sites[,"Ratio.M.L.normalized.limma.q"] < y1

downreg_ML <- sites[,"Ratio.M.L.normalized.limma.mean"] < -log2(x1) & 
  sites[,"Ratio.M.L.normalized.limma.q"] < y1

notreg_ML <- !upreg_ML & !downreg_ML

#plot
b<- c(log2(x1),-log2(x1))
ggplot(sites, aes(x=Ratio.M.L.normalized.limma.mean,
                  y=-log10(Ratio.M.L.normalized.limma.q)))+
  geom_point(data = subset(sites, notreg_ML), color = "grey70", alpha =0.5, size=3)+
  geom_point(data = subset(sites, upreg_ML), color = "red3", alpha =0.1, size=3)+
  geom_point(data = subset(sites, downreg_ML), color = "royalblue4", alpha =0.3, size=3)+
  geom_point(data=filter(sites, (sites[,"names"] %in% chosen_ones[,1]) & upreg_ML == TRUE), shape=21, color = "black", fill="red3", alpha=0.99, size=3)+
  geom_hline(yintercept = -log10(y1), linetype="dashed", color="grey48", size=0.5)+
  geom_vline(xintercept = b, linetype="dashed", color="grey48", size=0.5)+
  scale_y_continuous(breaks = c(0,1,1.3,2,3))+
  scale_x_continuous(breaks = c(-5,-2.5,-1,0,1,2.5,5))+
  geom_text_repel(data=filter(sites, (sites[,"names"] %in% chosen_ones[,1]) & (upreg_ML == TRUE)), size=5,
                  max.iter = 1500, box.padding=0.5,  force=1.0, segment.alpha=0.8, max.overlaps = 10000,
                  segment.color = 'grey48',
                  aes(label=names))+
  xlab(expression("log"[2]*"(GST-MyUb/GST)")) + ylab(expression("-log"[10]*"(FDR)"))+
  theme_bw(base_size = 24) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

dev.copy2pdf(file="interactome-MyUb-limma.pdf",width=10,height=8)

#technical info
sessionInfo()