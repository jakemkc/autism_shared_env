## March 4 2020
## Goal: heatmap of the MWAS sign in ASD

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
# cmd-shit-f10: restart R. reset loaded settings & library etc
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
#'
# ******** ---


# *****************************************************************************
## Notes
#' 1. mzid: v1 v2 v3..
#' 2. mz: 100.0023, 123.0039
# *****************************************************************************

library(tidyverse)

# ******** -----
# A. read-data ------------------------------------------------------------

load("results/hilic_plots.Rdata")

## Switch for saving
saveName = "hilic"
# saveName = "c18"


# ******** -----
# A. MWAS-Correlation ------------------------------------------------------------


autism_asd <- autism_ccEwas %>% filter(Role == "Proband")
autism_asdVar <- autism_asd %>% select(fdrsign$mzid, Age, Sex, Batch)

## recode sex and batch
autism_asdVar <- autism_asdVar %>% mutate(
    sexRe = case_when(
        Sex == "Female" ~ 0,
        Sex == "Male" ~ 1
    )
)

autism_asdVar$batchRe <- as.numeric(levels(autism_asdVar$Batch))[autism_asdVar$Batch]

## final dataframe ready
autism_asdVar <- autism_asdVar %>% select(-Sex, -Batch)

## corr
autismCorMat <- cor(autism_asdVar, use = "pairwise.complete.obs", method = "spearman")


## partial R
library(psych)

par.r <- partial.r(autismCorMat, c(1:nrow(fdrsign)), c("Age", "sexRe", "batchRe"))


## heatmap
library(gplots)

# Color function  # use cp group's common scheme for faster interpretation
heatmapColors <- function(numColors=16) {
    c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
    c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
    c3 <- c(c1,rev(c2));  # rev: reverse element
    return(c3)
}

## Rainbow explain
# rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
# n: no. of color
# s, v: satuartion and value to complete HSV colo
# start: hue in [0,1] to begin
# alpha: transpancy, in 0,1


quartz(height=6, width=8)

## No cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(par.r,
          main = "Spearman correlation of \n MWAS sign in ASD",
          trace="none", 
          margins=c(5,5), 
          col=heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,
          srtRow = 10,
          cexRow=0.2,  # row label size
          cexCol=0.25,  # column label size
          dendrogram= "none",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = F,
          Colv = F,
          symbreaks=T)   # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/-
# ColSideColors=varnamesToCorr$color,  # this is an extra color column best for showing group
# RowSideColors=varnamesToCorr$color)


outdirectory <- "results"
outfilename <- sprintf("heatmap_mwasSign_%s.pdf", saveName)
save(file=file.path(outdirectory,outfilename),
     

quartz.save(file.path(outdirectory,outfilename), type = "pdf", device = dev.cur(), dpi = 200)





# ******** -----
# B. MWAS-annotated-Corr ------------------------------------------------------------

readName <- sprintf("tablemz_sign_unique_%s.csv", saveName)
annoFDR <- read_csv(file.path("results", readName))


autism_asd <- autism_ccEwas %>% filter(Role == "Proband")
autism_asdVar <- autism_asd %>% select(annoFDR$mzid, Age, Sex, Batch)

## recode sex and batch
autism_asdVar <- autism_asdVar %>% mutate(
    sexRe = case_when(
        Sex == "Female" ~ 0,
        Sex == "Male" ~ 1
    )
)

autism_asdVar$batchRe <- as.numeric(levels(autism_asdVar$Batch))[autism_asdVar$Batch]

## final dataframe ready
autism_asdVar <- autism_asdVar %>% select(-Sex, -Batch)

names(autism_asdVar) <- c(annoFDR$Name, "Age", "sexRe", "batchRe")

## corr
autismCorMat <- cor(autism_asdVar, use = "pairwise.complete.obs", method = "spearman")


## partial R
library(psych)

par.r <- partial.r(autismCorMat, c(1:nrow(annoFDR)), c("Age", "sexRe", "batchRe"))


## heatmap
library(gplots)

# Color function  # use cp group's common scheme for faster interpretation
heatmapColors <- function(numColors=16) {
    c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
    c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
    c3 <- c(c1,rev(c2));  # rev: reverse element
    return(c3)
}

## Rainbow explain
# rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
# n: no. of color
# s, v: satuartion and value to complete HSV colo
# start: hue in [0,1] to begin
# alpha: transpancy, in 0,1


quartz(height=6, width=8)

## No cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(par.r,
          main = "Spearman correlation of \n MWAS sign in ASD",
          trace="none", 
          margins=c(5,7), 
          col=heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 20,
          srtRow = 60,
          cexRow=0.4,  # row label size
          cexCol=0.4,  # column label size
          dendrogram= "none",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = F,
          Colv = F,
          symbreaks=T)   # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/-
# ColSideColors=varnamesToCorr$color,  # this is an extra color column best for showing group
# RowSideColors=varnamesToCorr$color)


outdirectory <- "results"
outfilename <- sprintf("heatmap_mwasSign_anno_%s.pdf", saveName)


quartz.save(file.path(outdirectory,outfilename), type = "pdf", device = dev.cur(), dpi = 200)



# ******** -----
# Testing code ------------------------------------------------------------

# ## Testng code 030420
# library(psych)
# 
# jen <- make.hierarchical()    #make up a correlation matrix 
# lowerMat(jen[1:5,1:5])
# par.r <- partial.r(jen,c(1,3,5),c(2,4), method="spearman") # spearman or pearson won't affect because you provided corr matrix
# 
# par.r <- partial.r(jen, c(1,3,5), cs(V2,V4)) # same as above, cs a helper function in psy to convert list of text to char vector
# 
# 
# lowerMat(par.r)
# cp <- corr.p(par.r,n=98)  #assumes the jen data based upon n =100.
# print(cp,short=FALSE)  #show the confidence intervals as well
# #partial all from all correlations.
# lowerMat(partial.r(jen)) 
# 
# #Consider the Tal.Or data set.
# lowerCor(Tal.Or)
# #partial gender and age from these relations (they hardly change)
# partial.r(Tal.Or,1:4,cs(gender,age))
# #find the partial correlations between the first three variables and the DV (reaction)
# round(partial.r(Tal.Or[1:4])[4,1:3],2) #The partial correlations with the criterion