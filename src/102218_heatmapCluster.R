###' Mar 22 2018
###' Goal: heatmap

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
# cmd-shit-f10: restart R. reset loaded settings & library etc
cat("\014")   # same as ctrl-L
options(max.print = 3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# dev.new() # new plot window
# getwd()


# *****************************************************************************
###' Goal
# *****************************************************************************


# *****************************************************************************
## Notes
#' 1. 
# *****************************************************************************

# *****************************************************************************
## Checkpt cmts
# *****************************************************************************

# *****************************************************************************
## Hacks (Ref to previously)
#' 1) data dict in csv ready to view in number.app; "rvar" for quick view
#' 2) key input and output for script: Write down dim() and info; 
#'    Conceptual clarity esp for SQL
#' 3) orignal key data in CSV or STA (if big) for quick browsing
# *****************************************************************************

## load
library(tidyverse)

load("results/hilic_corr_sign.Rdata")


## Switch for saving
saveName = "hilic"
# saveName = "c18"


# ******** ----
# A. heatmap  ----

library(gplots)

heatmapColors <- function(numColors=16) {
    c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
    c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
    c3 <- c(c1,rev(c2));
    return(c3)
}


## notes
# the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
# default complete
# 
# he distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
# default euclidean

## asd
heatmap.2(cor_asd,
          main = "ASD correlation",
          trace = "none", 
          margins = c(5,5), 
          col = heatmapColors(4), 
          srtCol = 70,
          srtRow = 10,
          cexRow = 0.2,  
          cexCol = 0.25,  
          dendrogram = "both",  
          Rowv = T,
          Colv = T,
          # RowSideColors = heatmap_color_f,
          symbreaks = T, 
          symkey = T)   

# save
outdirectory <- "./results"
outfilename <- sprintf("heatmapAsd_%s.pdf", saveName)
quartz.save(file.path(outdirectory,outfilename), type = "pdf", device = dev.cur(), dpi = 300)


## parents
heatmap.2(cor_parents,
          main = "Parent correlation",
          trace = "none", 
          margins = c(5,5), 
          col = heatmapColors(4), 
          srtCol = 70,
          srtRow = 10,
          cexRow = 0.2,  
          cexCol = 0.25,  
          dendrogram = "both",  
          Rowv = T,
          Colv = T,
          # RowSideColors = heatmap_color_f,
          symbreaks = T, 
          symkey = T)   

# save
outdirectory <- "./results"
outfilename <- sprintf("heatmapParents_%s.pdf", saveName)
quartz.save(file.path(outdirectory,outfilename), type = "pdf", device = dev.cur(), dpi = 300)

## family
heatmap.2(cor_family,
          main = "Family correlation",
          trace = "none", 
          margins = c(5,5), 
          col = heatmapColors(4), 
          srtCol = 70,
          srtRow = 10,
          cexRow = 0.2,  
          cexCol = 0.25,  
          dendrogram = "both",  
          Rowv = T,
          Colv = T,
          # RowSideColors = heatmap_color_f,
          symbreaks = T, 
          symkey = T)   

# save
outdirectory <- "./results"
outfilename <- sprintf("heatmapFamily_%s.pdf", saveName)
quartz.save(file.path(outdirectory,outfilename), type = "pdf", device = dev.cur(), dpi = 300)



# ******** ----
# B. Cluster  ----
#' get group by cluster instead of RT
#' 1) run cluster analysis (use default) on ASD kids
#' 2) plot dendrogram, decide how many groups
#' 3) color plot the cut in dendrogram, make sure it is
#' 4) get dendrogram chemical name ordered vector
#' 5) get chemical name and group ID vector
#' 6) Cbind to create sorting and conversion table of chemical names and groups


## cluster analysis
res.hc <- cor_asd %>%
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "ward.D2")     # Compute hierachical clustering


## plot
plot(res.hc)
# about 8-10 groups

library("factoextra")

fviz_dend(res.hc, k = 9, # Cut in four groups
          cex = 0.6, # label size
          k_colors = c("#89C5DA", "#DA5724", "#689030", "#CE50CA"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

outdirectory <- "./results"
outfilename <- sprintf("clusterSign_%s.png", saveName)
dev.copy(png, filename=file.path(outdirectory,outfilename), width=1024, height=768)
dev.off()


#' use 9 groups

## Dendrogram chemical ordered vector 
res.hcD <- cor_asd %>%
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "ward.D2") %>% as.dendrogram()     # Compute hierachical clustering

chemOrder <- tibble(chemID = labels(res.hcD)) # tibble instead of df: string NOT factor, row names NULL

## chemical names and group ID
hclustercut <- cutree(res.hc, k = 9)

tmp1 <- tibble(group = hclustercut, chemID = names(hclustercut)) # tibble instead of df: string NOT factor, row names NULL


## create sorting and conversion dataframe

chemOrder <- left_join(chemOrder, tmp1, by = "chemID")

chemOrder$chemIDp <- str_replace(chemOrder$chemID, "_a$", "_p")
# str_view_all(chemOrder$chemID, "_a$")  # Testing. "_a" at the end of string. see Stringr cheatsheet

chemOrder$mzid <- str_replace(chemOrder$chemID, "_a$", "")

## save
outdirectory <- "./results"
outfilename <- sprintf("chemOrder_%s.rds", saveName)
saveRDS(chemOrder, file.path(outdirectory,outfilename))
# hello <- readRDS("results/chemOrder.rds")

