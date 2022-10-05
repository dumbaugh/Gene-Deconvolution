#### SpatialDecon
library(readxl)
library("SpatialDecon")
library(tidyverse)

#Objective: Use a previously created single cell expression matrix from APAP-treated mice to deconvolve bulk RNAseq data from mice
# treated with APAP300 or APAP300+MSC (cell based therapy)

### format properly for spatialDecon ++ read in single cell expression matrices ####
#rownames are genes, colnames are cell labels
#read in APAP 24h expression matrix
AvgExp_APAP24h <- read_excel("AvgExp_APAP24h.xlsx") #this is the single cell expression matrix that defines celltypes
AvgExp_APAP24h<- distinct(AvgExp_APAP24h) %>% na.omit() %>% 
  as.data.frame(AvgExp_APAP24h) 
rownames(AvgExp_APAP24h)<- make.names(AvgExp_APAP24h[,1], unique = TRUE)
AvgExp_APAP24h <- AvgExp_APAP24h[,-1]
AvgExp_APAP24h_mat <- as.matrix(AvgExp_APAP24h)


#read in 48h ####

MSC48_normalizedcounts <- read_excel("MSC48_normalizedcounts.xlsx")

#remove NAs, remove duplicated genes, ensure unique gene names
MSC48_normalizedcounts<- distinct(MSC48_normalizedcounts) %>% na.omit() %>% 
  as.data.frame(MSC48_normalizedcounts) 
rownames(MSC48_normalizedcounts)<- make.names(MSC48_normalizedcounts[,1], unique = TRUE)
MSC48_normalizedcounts <- MSC48_normalizedcounts[,-1]


## try with dplyr
#remove all rows (e.g. genes) that have a count<2
MSC48_normalizedcounts_filtered <- MSC48_normalizedcounts %>% 
  mutate(minGene=apply(MSC48_normalizedcounts,1,FUN=min)) %>% 
  filter(minGene>2)

#remove minGene column
MSC48_normalizedcounts_filtered <- MSC48_normalizedcounts_filtered[,-7]

MSC48_mat <- as.matrix(MSC48_normalizedcounts)
MSC48_mat_filtered <- as.matrix(MSC48_normalizedcounts_filtered)



#perform spatial deconvolution
msc48_decon_filtered <- spatialdecon(MSC48_mat_filtered,bg=1,X=AvgExp_APAP24h_mat,align_genes = T,maxit=1000)



layout(mat = (matrix(c(1, 2), 1)), widths = c(3, 3))
TIL_barplot(msc48_decon_filtered$prop_of_all, draw_legend = TRUE, 
            cex.names = 0.8)

## non filtered matrix
msc48_decon <- spatialdecon(MSC48_mat,bg=1,X=AvgExp_APAP24h_mat,align_genes = T,maxit=1000)

#visualize cell proportions
layout(mat = (matrix(c(1, 2), 1)), widths = c(5, 2))
TIL_barplot(msc48_decon$prop_of_all, draw_legend = TRUE, 
            cex.names = 1.2)

heatmap(msc48_decon$beta, cexCol = 2, cexRow = 2, margins = c(7,10),Colv=NA)
heatmap(msc48_decon_filtered$beta, cexCol = 2, cexRow = 2, margins = c(7,9),Colv = NA)

## reverse deconvolution
# can either add 1 to the matrix
msc48_decon_rev = reverseDecon(norm = MSC48_mat+1,
                               beta = msc48_decon$beta)
#or use the filtered matrix which will remove 0s
msc48_decon_filtered_rev = reverseDecon(norm = MSC48_mat_filtered,
                                        beta = msc48_decon_filtered$beta)

#evaluate 
plot(msc48_decon_rev$cors, msc48_decon_rev$resid.sd, col = 0)
showgenes = c("Cxcl14", "Lyz", "Nkg7")
text(msc48_decon_rev$cors[setdiff(names(msc48_decon_rev$cors), showgenes)], 
     msc48_decon_rev$resid.sd[setdiff(names(msc48_decon_rev$cors), showgenes)], 
     setdiff(names(msc48_decon_rev$cors), showgenes), cex = 0.5)
text(msc48_decon_rev$cors[showgenes], msc48_decon_rev$resid.sd[showgenes], 
     showgenes, cex = 0.75, col = 2)


plot(msc48_decon_filtered_rev$cors, msc48_decon_filtered_rev$resid.sd, col = 0)
showgenes = c("Cyp2e1", "Cyp2f2", "Ccr2","S100a8","Nkg7","Wnt2","Ntn4","Lyve1","Mmp8")
text(msc48_decon_filtered_rev$cors[setdiff(names(msc48_decon_filtered_rev$cors), showgenes)], 
     msc48_decon_filtered_rev$resid.sd[setdiff(names(msc48_decon_filtered_rev$cors), showgenes)], 
     setdiff(names(msc48_decon_filtered_rev$cors), showgenes), cex = 0.5)
text(msc48_decon_filtered_rev$cors[showgenes], msc48_decon_filtered_rev$resid.sd[showgenes], 
     showgenes, cex = 0.75, col = 2)

which(rownames(msc48_decon_filtered_rev[["coefs"]])=="Cyp2e1") #note no Ngp/Camp interesting!
msc48_decon_filtered_rev[["coefs"]][7622,]
rownames(msc48_decon_filtered_rev$coefs)[7622]

heatmap(pmax(pmin(msc48_decon_rev$resids, 2), -2))
heatmap(pmax(pmin(msc48_decon_filtered_rev$resids, 2), -2))
head(msc48_decon_filtered_rev$resids)
str(msc48_decon_filtered_rev)

top100 <- msc48_decon_filtered_rev$resids[1:200,1:6]
heatmap(pmax(pmin(top100, 2), -2))
heatmap(top100)
?pmax
#read in 72h ####
MSC72_normalizedcounts <- read_excel("MSC72_normalizedcounts.xlsx")

MSC72_normalizedcounts<- distinct(MSC72_normalizedcounts) %>% na.omit() %>% 
  as.data.frame(MSC72_normalizedcounts) 
rownames(MSC72_normalizedcounts)<- make.names(MSC72_normalizedcounts[,1], unique = TRUE)
MSC72_normalizedcounts <- MSC72_normalizedcounts[,-1]


## try with dplyr
#remove all rows (e.g. genes) that have a count<2
MSC72_normalizedcounts_filtered <- MSC72_normalizedcounts %>% 
  mutate(minGene=apply(MSC72_normalizedcounts,1,FUN=min)) %>% 
  filter(minGene>2)

#remove minGene column
MSC72_normalizedcounts_filtered <- MSC72_normalizedcounts_filtered[,-7]

MSC72_mat <- as.matrix(MSC72_normalizedcounts)
MSC72_mat_filtered <- as.matrix(MSC72_normalizedcounts_filtered)

msc72_decon_filtered <- spatialdecon(MSC72_mat_filtered,bg=1,X=AvgExp_APAP24h_mat,align_genes = T,maxit=1000)

layout(mat = (matrix(c(1, 2), 1)), widths = c(3, 3))
TIL_barplot(msc72_decon_filtered$prop_of_all, draw_legend = TRUE, 
            cex.names = 0.8)

msc72_decon <- spatialdecon(MSC72_mat,bg=1,X=AvgExp_APAP24h_mat,align_genes = T,maxit=1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(5, 2))
TIL_barplot(msc72_decon$prop_of_all, draw_legend = TRUE, 
            cex.names = 1.2)

heatmap(msc72_decon$beta, cexCol = 2, cexRow = 2, margins = c(7,10),Colv=NA)
heatmap(msc72_decon_filtered$beta, cexCol = 2, cexRow = 2, margins = c(7,9),Colv = NA)