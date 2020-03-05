## Feb 21 2018
## Goal: FDR sign and annotation

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
## Goal: annotation of FDR sign
# Steps
#' 1. retframe from EWAS get mz columns
#' 1.1. remove duplicate mz id any
#'    special: get full retframe with mz and time, ready for mummichog & save to csv
#' 2. join 1.1) with hmdb
#'    ALL these emory annotation could have no unique mz matched name (e.g. mz 100; got 2 names)
#' 3. join 1.1 with kegg
#' 4. join 1.1 with lipid map
#' 5. merge the joins
#' 5.1. estimate unique hits etc
# ******** ---


# *****************************************************************************
## Notes
#' 1. mzid: v1 v2 v3..
#' 2. mz: 100.0023, 123.0039
# *****************************************************************************




# ******** -----
# A. read-data ------------------------------------------------------------

# load("results/hilic_plots.Rdata") 
load("results/hilicScale_plots.Rdata")
hilic_hmdb <- read.csv(file = "data/hil_new_Stage5_HMDB_pos.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
hilic_kegg<- read.csv(file = "data/hil_new_Stage5_KEGG_pos.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
hilic_lipid<- read.csv(file = "data/hil_new_Stage5_LM_pos.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
# MELTIN no batch d/l, leave it first


## Switch for saving
saveName = "hilic"
# saveName = "c18"



# ******** -----
# B.  join-annotation ------------------------------------------------------------
library(tidyverse)

rownames(fdrsign) <- fdrsign$mzid 
# Warning: Setting row names on a tibble is deprecated
# colnames(fdrsign) <- c("Estimate", "SE", "z", "p", "mzid", "OR", "fdr")

### add mzid column (code from explore.R)
hilic1 <- data_original
hilic1$mzid <- paste("v", c(1:nrow(hilic1)), sep="")
# hilic1$mzid <- paste("v", c(1:13098), sep="")
# hilic1 %>% t %>% .[, 1:3] %>% View
hilic1 <- hilic1[, c("mz","time","mzid")]


## left join to get mz and RT
fdrsign_1 <- left_join(fdrsign, hilic1, by = c("mzid" = "mzid")) # rownames gone
# colnames(fdrsign_1)
# head(fdrsign_1)
# str(fdrsign_1)
# dim(fdrsign_1) # 66


# remove duplicate mz if any
#' no duplicate in 101118 code
length(unique(fdrsign_1$mz)) # 66 sign mz; one is duplicate
which(duplicated(fdrsign_1$mz)) # index 34
fdrsign_1[33:34, ]
# Estimate        SE         z            p  mzid    OR        fdr       mz  time
# 33 -2.884716 0.7313469 -3.944388 0.0000800041 v1052 0.056 0.01273132 130.4919 200.3
# 34 -5.372138 1.5313192 -3.508176 0.0004511901 v1053 0.005 0.03566195 130.4919 273.8


# .\ prep-mummichog ----
    ###' For metaboanalyst mummichog run (retframe is not updated with the dropped var from EWAS.R)
    ###' mummichog needs full ewas loop results, with mz info (and RT?), and a p value cut off
    
# retFrame %>% colnames
# 
rownames(retFrame) <- retFrame$mzid
# Warning: Setting row names on a tibble is deprecated
# colnames(retFrame) <- c("Estimate", "SE", "z", "p", "mzid", "OR", "fdr", "FDR groups")

### add mzid column (code from explore.R)
# hilic1 <- data_original
# hilic1$mzid <- paste("v", c(1:13098), sep="")
# hilic1 <- hilic1[, c(1,2,270)]

## left join to get mz and RT
retFrame1 <- left_join(retFrame, hilic1, by = c("mzid" = "mzid"))
# colnames(retFrame1)
# head(retFrame1)
# str(retFrame1)
# dim(retFrame1)

# p value cut off for metaboanalyst
fdrsign$p.value %>% max # 0.0003942237
    ###' mummichog example uses tscore, but I am providing glm's zscore and p value instead

outdirectory <- "results"
outfilename <- sprintf("%s_ewas_retFrame.csv", saveName)
write.csv(retFrame1, file = file.path(outdirectory,outfilename))
# write.csv(retFrame1, file = "results/hilic_ewas_retFrame.csv")



# .\ join-hmdb ----

## cut hmdb size
# head(hilic_hmdb)
# str(hilic_hmdb)
# colnames(hilic_hmdb)
# colnames(hilic_hmdb)[1:11]
# [1] "MatchCategory"    "mz"               "MatchCategory.1"  "theoretical.mz"   "HMDBID"          
# [6] "Name"             "Formula"          "MonoisotopicMass" "Adduct"           "AdductMass"      
# [11] "time"  

# hilic_hmdb <- hilic_hmdb[, 1:11 ] # only get some useful category

## left join to get hmdb info
# row increases because hmdb non unique name match to mz
fdrsign_hmdb <- left_join(fdrsign_1, hilic_hmdb, by = c("mz", "mz"))

length(unique(fdrsign_hmdb$mz)) # 66 same as before joining

table(fdrsign_hmdb$MatchCategory)
# Multiple   Unique 
#   142       4
dim(fdrsign_hmdb)

fdrsign_hmdb_unique <- fdrsign_hmdb %>% filter(MatchCategory == "Unique")
fdrsign_hmdb_multi <- fdrsign_hmdb %>% filter(MatchCategory == "Multiple")

fdrsign_hmdb_unique$mz %>% unique %>% length  # 4 unique mz in unique HMDB ID
fdrsign_hmdb_multi$mz %>% unique %>% length  # 15 unique mz in mutiple HMDB ID

# no annoation count
tmp1 <- fdrsign_hmdb$MatchCategory %>% is.na(.)
tmp2 <- fdrsign_hmdb[tmp1, ] %>% .$mz %>% unique %>% length # 47
    # 47 + 4 + 15 = 66 unique mz


# .\ join-kegg ----
## cut kegg size
# head(hilic_kegg)
# str(hilic_kegg)
# colnames(hilic_kegg)
# hilic_kegg <- hilic_kegg[, 1:11 ]

## left join to get hmdb info
# row increases because source non unique name match to mz
fdrsign_kegg <- left_join(fdrsign_1, hilic_kegg, by = c("mz", "mz"))

length(unique(fdrsign_kegg$mz)) # 66 same as before joining

table(fdrsign_kegg$MatchCategory)
# Multiple   Unique 
# 87       7 
# dim(fdrsign_kegg)
fdrsign_kegg_unique <- fdrsign_kegg %>% filter(MatchCategory == "Unique")
fdrsign_kegg_multi <- fdrsign_kegg %>% filter(MatchCategory == "Multiple")

fdrsign_kegg_unique$mz %>% unique %>% length  # 7 unique mz in unique kegg ID
fdrsign_kegg_multi$mz %>% unique %>% length  # 11 unique mz in mutiple kegg ID

tmp1 <- fdrsign_kegg$MatchCategory %>% is.na(.)
tmp2 <- fdrsign_kegg[tmp1, ] %>% .$mz %>% unique %>% length # 48
# 48 + 7 + 11 = 66 unique mz


# .\ join-lipid ----
## cut kegg size
# head(hilic_lipid)
# str(hilic_lipid)
# colnames(hilic_lipid)
# hilic_lipid <- hilic_lipid[, 1:11 ]

## left join to get hmdb info
# row increases because source non unique name match to mz
fdrsign_lipid <- left_join(fdrsign_1, hilic_lipid, by = c("mz", "mz"))

length(unique(fdrsign_lipid$mz)) # 66 same as before joining

table(fdrsign_lipid$MatchCategory)
# Multiple   Unique 
# 325        1
# dim(fdrsign_lipid)
fdrsign_lipid_unique <- fdrsign_lipid %>% filter(MatchCategory == "Unique")
fdrsign_lipid_multi <- fdrsign_lipid %>% filter(MatchCategory == "Multiple")

fdrsign_lipid_unique$mz %>% unique %>% length  # 1 unique mz in unique kegg ID
fdrsign_lipid_multi$mz %>% unique %>% length  # 14 unique mz in mutiple HMDB ID

tmp1 <- fdrsign_lipid$MatchCategory %>% is.na(.)
tmp2 <- fdrsign_lipid[tmp1, ] %>% .$mz %>% unique %>% length # 51
# 51 + 1 + 14 = 66 unique mz



## Unique and multiple hit mz count
    ## set operations on 2 sets, I have 3 sets, with repeated potential. use them when I know more of the cmd.

tmph1 <- fdrsign_hmdb_unique$mz %>% unique  # 4 unique mz in unique HMDB ID
tmph2 <- fdrsign_hmdb_multi$mz %>% unique # 15 unique mz in mutiple HMDB ID

tmpk1 <- fdrsign_kegg_unique$mz %>% unique  # 7 unique mz in unique kegg ID
tmpk2 <- fdrsign_kegg_multi$mz %>% unique  # 11 unique mz in mutiple HMDB ID

tmpl1 <- fdrsign_lipid_unique$mz %>% unique  # 1 unique mz in unique kegg ID
tmpl2 <- fdrsign_lipid_multi$mz %>% unique  # 14 unique mz in mutiple HMDB ID

c(tmph1, tmpk1, tmpl1) %>% unique %>% length # 10 uniquely named
c(tmph2, tmpk2, tmpl2) %>% unique %>% length # 22 mutiple named


## make unique match table
colnames(fdrsign_hmdb_unique) %>% dput # just mod the dataseID to a general ID

# need same colnames for rbind
colnames(fdrsign_hmdb_unique) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", 
                                   "conf.high", "mzid", "OR", "ORCIlow", "ORCIhigh", "fdr", "by", 
                                   "mz", "time.x", "chemical_ID", "Confidence", "score", "Module_RTclust", 
                                   "time.y", "MatchCategory", "theoretical.mz", "delta_ppm", "Name", 
                                   "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "mean_int_vec", 
                                   "MD")

fdrsign_hmdb_both <- bind_rows(fdrsign_hmdb_unique, fdrsign_hmdb_multi)

colnames(fdrsign_kegg_unique) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", 
                                   "conf.high", "mzid", "OR", "ORCIlow", "ORCIhigh", "fdr", "by", 
                                   "mz", "time.x", "chemical_ID", "Confidence", "score", "Module_RTclust", 
                                   "time.y", "MatchCategory", "theoretical.mz", "delta_ppm", "Name", 
                                   "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "mean_int_vec", 
                                   "MD")

fdrsign_kegg_both <- bind_rows(fdrsign_kegg_unique, fdrsign_kegg_multi)

colnames(fdrsign_lipid_unique) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", 
                                    "conf.high", "mzid", "OR", "ORCIlow", "ORCIhigh", "fdr", "by", 
                                    "mz", "time.x", "chemical_ID", "Confidence", "score", "Module_RTclust", 
                                    "time.y", "MatchCategory", "theoretical.mz", "delta_ppm", "Name", 
                                    "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "mean_int_vec", 
                                    "MD")

fdrsign_lipid_both <- bind_rows(fdrsign_lipid_unique, fdrsign_lipid_multi)



# rbind 3 unique match tables
tmp3 <- rbind(fdrsign_hmdb_both, fdrsign_kegg_both, fdrsign_lipid_both) # 12 obs mz

# extract the unique among the rbinded above
tmp4 <- c(tmph1, tmpk1, tmpl1) %>% unique #10 obs mz
tmp_logic <- duplicated(tmp3$mz)

tablemz_sign_unique <- tmp3[which(!tmp_logic), ]

# both
tmp4 <- tmp3 %>% filter(Adduct == "M+H")


## special: revising table 052019


# ## datapasta::df_paste()
# tmp1 <- data.frame(stringsAsFactors=FALSE,
#         dbID = c("C13673", "HMDB60005", "HMDB40901", "HMDB29622", "C10061",
#                  "C09191", "HMDB38748", "HMDB01932", "C20973", "LMFA03020066",
#                  "C17770", "C16367", "HMDB37539", "C08566"),
#         name = c("(S)-ACPA;(S)-2-Amino-3-(3-carboxy-5-methyl-4-isoxazolyl)porionic acid", "Methylgallic acid-O-sulphate",
#                  "13-Hydroxy-9-methoxy-10-oxo-11-octadecenoic acid",
#                  "Isopetasoside", "Euxanthone;17-Dihydroxyxanthone", "Fruticosonine",
#                  "Osmanthuside A", "Metoprolol",
#                  "Angiotensin (5-7);Angiotensin-(5-7)", "pentahydroxyicosa-681014-tetraenoic acid", "Sarcostin",
#                  "Delphinidin 3-O-(6-caffeoyl-beta-D-glucoside)", "6-Caffeoylhyperin",
#                  "Portulacaxanthin III"),
#    chemclass = c("Amino acid", "Benzenoids", "Lipids and lipid-like molecules",
#                  "Lipids and lipid-like molecules", "Xanthones", "Tryptamines",
#                  "Glycoside", "Benzenoids", "Peptide",
#                  "Lipids and lipid-like molecules", "Hormone", "Glycoside", "Glycoside", "Amino acid")
# )
# 
# 
# 
# 
# tmp1 <- tmp1 %>% select(dbID, chemclass)
# tablemz_sign_unique <- tablemz_sign_unique %>% left_join(tmp1, by = "dbID")
# tmp2 <- tablemz_sign_unique %>% select(MonoisotopicMass, OR, ORCIlow, ORCIhigh, fdr, dbID, Name, chemclass)
# tmp2 <- tmp2 %>% arrange(OR)
# tmp2 %>% copydat::Copy()





outdirectory <- "results"
outfilename <- sprintf("tablemz_sign_both_supp_%s.csv", saveName)
write.csv(tmp4, file = file.path(outdirectory,outfilename))


outdirectory <- "results"
outfilename <- sprintf("tablemz_sign_both_supp_%s.rds", saveName)
saveRDS(tmp4, file = file.path(outdirectory,outfilename))
    # hello <- readRDS("results/tablemz_sign_unique.rds")


# ******** -----
## S. save-rdata ----

# outdirectory <- "results"
# outfilename <- "1_1_50_subset.Rdata"
# # outfilename <- sprintf("%s_reg_7.Rdata", depVariable)
# save(file=file.path(outdirectory,outfilename), 
#      cleanup1, cleanup2, cleanup3, drop0, dhs_original)

# Load data
# load("results/1_1_50_subset.Rdata") 
# 
# outdirectory <- "results"
# outfilename <- sprintf("fdrsign_hmdb_unique_%s.csv", saveName)
# write.csv(fdrsign_hmdb_unique, file = file.path(outdirectory,outfilename))
# 
# outfilename <- sprintf("fdrsign_kegg_unique_%s.csv", saveName)
# write.csv(fdrsign_kegg_unique, file = file.path(outdirectory,outfilename))
# 
# outfilename <- sprintf("fdrsign_lipid_unique_%s.csv", saveName)
# write.csv(fdrsign_lipid_unique, file = file.path(outdirectory,outfilename))




# ******** -----
# TODO -----------------------------------------------------------------




# ******** -----
# Supplement -----------------------------------------------------------



# .\ sub-tasks ----
# .\ sub-tasks \ tasks ----

# 
# 
# # ******** -----
# # Playground, Learning: Send to Quiver
# 
# spearM$row2 <- gsub("_m$", "", spearM$row, ignore.case = T)
# 
# spearM$col_row_pair_2 <- sprintf("%s-%s", spearM$col2, spearM$row2)
# 
# 
# 
# 
# # ******** -----
# # P1.  try cluster with all mz ------------------------------------------------------------
# 
# ### extract sign mz
# 
# 
# fdrsign_hclust <- autism_n3
# 
# 
# colnames(fdrsign_hclust) %>% length # 193 + 9 = 202
# 
# ### subset with proband
# fdrsign_hclust1 <- fdrsign_hclust %>% filter(. , grepl("Proband", Role))
# 
# fdrsign_hclust1$Diag
# 
# fdrsign_hclust1$Diag %>% table(useNA = "ifany")
# 
# # 1-autism          2-asd      3-nonspectrum   could-not-determine-cut-off-value 
# # 58                7          6               4 
# # no-diagnosis      <NA> 
# # 1                 6 
# 
# fdrsign_hclust1$subid <- as.character(NA)
# 
# tmp_ind1 <- which(fdrsign_hclust1$Diag %in% "1-autism")
# fdrsign_hclust1$subid[tmp_ind1] <- sprintf("%02d", 1:58)
# 
# tmp_ind2 <- which(fdrsign_hclust1$Diag %in% "2-asd")
# fdrsign_hclust1$subid[tmp_ind2] <- sprintf("%02d", 1:7)
# 
# tmp_ind3 <- which(fdrsign_hclust1$Diag %in% "3-nonspectrum")
# fdrsign_hclust1$subid[tmp_ind3] <- sprintf("%02d", 1:6)
# 
# tmp_ind4 <- which(fdrsign_hclust1$Diag %in% "could-not-determine-cut-off-value")
# fdrsign_hclust1$subid[tmp_ind4] <- sprintf("%02d", 1:4)
# 
# tmp_ind5 <- which(fdrsign_hclust1$Diag %in% "no-diagnosis")
# fdrsign_hclust1$subid[tmp_ind5] <- sprintf("%02d", 1)
# 
# ### check
# fdrsign_hclust1$subid %>% table(useNA = "ifany")
# 
# ### create unique row id
# fdrsign_hclust1$rowdiag <- sprintf("%s_%s", fdrsign_hclust1$Diag, fdrsign_hclust1$subid)
# fdrsign_hclust2 <- na.omit(fdrsign_hclust1)
# rownames(fdrsign_hclust2) <- fdrsign_hclust2$rowdiag
# 
# colnames(fdrsign_hclust2)
# 
# fdrsign_hclust3 <- fdrsign_hclust2[, 1:11935]
# 
# 
# ### Compute hierarchical clustering
# library("factoextra")
# 
# res.hc <- fdrsign_hclust3 %>%
#     +1 %>% 
#     log10 %>% 
#     scale %>%                    # Scale the data
#     dist(method = "euclidean") %>% # Compute dissimilarity matrix
#     hclust(method = "ward.D2")     # Compute hierachical clustering
# 
# ### Visualize using factoextra
# ### Cut in 4 groups and color by groups
# fviz_dend(res.hc, k = 4, # Cut in four groups
#           cex = 0.8, # label size
#           k_colors = c("#89C5DA", "#DA5724", "#689030", "#CE50CA"),
#           color_labels_by_k = TRUE, # color labels by groups
#           rect = TRUE # Add rectangle around groups
# )
# 
# # dev.copy(png, filename="results/hcluster_ASD_82_all.png", width=1024, height=768)
# # dev.off()
# 
