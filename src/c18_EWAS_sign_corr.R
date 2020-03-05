## Feb 14 2019
## Goal: correlation of EWAS sign. meta to 1) family 2) sex 3) age

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
# cmd-shit-f10: restart R. reset loaded settings & library etc
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
#' see: "analysis plan on correlation 021219.txt" for full
#' 
#' A) prep Data for correlation
#' 1. From autism df (row = 197 obs, col: mz + meta info), select only EWAS sign. mz
#' 2. filter to get asd and parents
#' 3. create familyID (exposome globe lines match in family lv); for both asd ad parents
#' 3.1. remove dup ID
#' 3.2. average values per family ID, for both asd and parents (familyID as rows)
#' 4. intersect the row between asd and parents for both asd and parent df (match rows)
#' 5 remove mz var = 0 in both
#' 6. rename mz with asd or parent tag in both dataframe
#' 
#' Special: value 0
#' 1) now I adjust for covarite, and get residue, and I basically have no "0" to distort correlation
#' 
#' 
#' C) correlation
#' 1. I have 2 row matched dataframe with row = family, col = sign mz; run correlation; 
#'    * no need familyID for corr() but know rows match between asd and parents
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

load("results/c18_plots.Rdata")
# load("results/_= c18/_ 20 percent/c18_plots.Rdata")

## Switch for saving
# saveName = "hilic"
saveName = "c18"


# ******** -----
# B1. family 050719 ------------------------------------------------------------

## *********** 
## Data

## ASD and parent mzID
# for correlation variable labels (parents and autism)
mzmatch_table <- fdrsign
mzmatch_table$mzid_a <- sprintf("%s_a", mzmatch_table$mzid)
mzmatch_table$mzid_p <- sprintf("%s_p", mzmatch_table$mzid)


## autism df pick only EWAS sign mz
# name of meta in autism df
tmpMetaName <- autism %>% select(-starts_with("v")) %>% colnames # 11

# autism df with sign mz and meta col only
autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))


# check. from corr.test warning. 050819
#' Y needs at need 2 points to get the regression and predicted/residual
# autism_sign %>% filter(Role == "Mother") %>% select(v6856)
# autism_sign %>% filter(Role == "Mother") %>% select(v6870)
# autism_sign %>% filter(Role == "Mother") %>% select(v7103)
# autism_sign %>% filter(Role == "Father") %>% select(v8255)
# autism_sign %>% filter(Role == "Mother") %>% select(v8824)


## *********** 
# B1. adjust and residual prep

# adjust age and sex; switch the whole section on/off

l <- list() # set tmp

## impute missing age, by role 
# (they are diff age group, not a single distribution)
l[["a"]] <- autism_sign %>% filter(. , grepl("Father|Mother", Role)) %>% select(Age)
l[["parentsub"]] <- Hmisc::impute(l[["a"]][[1]]) %>% as.numeric #median impute
autism_sign$Age[which(autism_sign$Role %in% c("Father", "Mother"))] <- l[["parentsub"]]

l[["b"]] <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
l[["asdsub"]] <- Hmisc::impute(l[["b"]][[1]]) %>% as.numeric #median impute
autism_sign$Age[which(autism_sign$Role %in% c("Proband"))] <- l[["asdsub"]]


## reference EST corr, create dataframe for each 
autism_sign_case <- autism_sign %>% filter(Role %in% c("Proband"))
autism_sign_mother <- autism_sign %>% filter(Role %in% c("Mother"))
autism_sign_father <- autism_sign %>% filter(Role %in% c("Father"))


## filter X with min 20% data in each df

# function: percentage of "0"
percentzero <- function (x) {
    100*mean(!(x))
}

filterDF <- function(autism_sign_case) {
    loop_mz <- autism_sign_case %>% select(starts_with("v")) %>% colnames
    meta <- autism_sign_case %>% select(-starts_with("v"))
    
    tmp1 <- autism_sign_case %>% select(loop_mz)
    tmp2 <- sapply(tmp1, percentzero)   # sum(autism_n$v1 %in% 0) / nrow(autism_n) * 100
    tmp3 <- names(tmp2)[which(tmp2 < 80)] # find names < 80% "0", i.e., keep mzid names
    tmp4 <- tmp1[, tmp3] # df with needed variables
    tmp5 <- bind_cols(tmp4, meta)
}

autism_sign_case1 <- filterDF(autism_sign_case) 
dim(autism_sign_case1) # 82,119 with meta
autism_sign_mother1 <- filterDF(autism_sign_mother)
dim(autism_sign_mother1) # 42, 96
autism_sign_father1 <- filterDF(autism_sign_father)
dim(autism_sign_father1) # 44, 98


## get mutual Xs across pair
# Reduce(intersect, list(names(autism_sign_case1), names(autism_sign_mother1), names(autism_sign_father1)))

# asd and mother
tmp1 <- intersect(names(autism_sign_case1), names(autism_sign_mother1))
autism_cmX <- autism_sign_case1 %>% select(tmp1)
autism_mX <- autism_sign_mother1 %>% select(tmp1)


# asd and father
tmp1 <- intersect(names(autism_sign_case1), names(autism_sign_father1))
autism_cfX <- autism_sign_case1 %>% select(tmp1)
autism_fX <- autism_sign_father1 %>% select(tmp1)


# get residual
# stratify instead of get residual from a single model 050719

library(broom)

adjFunASDMother <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_cmX, Role == "Proband")) %>% augment()
    ret[[".resid"]]
}


adjFunASDFather <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_cfX, Role == "Proband")) %>% augment()
    ret[[".resid"]]
}


adjFunMother <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_mX, Role == "Mother")) %>% augment()
    ret[[".resid"]]
}

adjFunFather <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_fX, Role == "Father")) %>% augment()
    ret[[".resid"]]
}


# asd mohter
loop_mz <- autism_cmX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunASDMother(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_cmX %>% select(-starts_with("v"))

autism_cmX_res <- cbind(residual, nonmz)


# mother
loop_mz <- autism_mX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunMother(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_mX %>% select(-starts_with("v"))

autism_mX_res <- cbind(residual, nonmz)


# asd father
loop_mz <- autism_cfX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunASDFather(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_cfX %>% select(-starts_with("v"))

autism_cfX_res <- cbind(residual, nonmz)


# fahter
loop_mz <- autism_fX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunFather(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_fX %>% select(-starts_with("v"))

autism_fX_res <- cbind(residual, nonmz)



## create familyID from BCH "ID"

autism_cmX_res$ID2 <- autism_cmX_res$ID
autism_cmX_res <- autism_cmX_res %>%
    separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")

autism_cfX_res$ID2 <- autism_cfX_res$ID
autism_cfX_res <- autism_cfX_res %>%
    separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")

autism_mX_res$ID2 <- autism_mX_res$ID
autism_mX_res <- autism_mX_res %>%
    separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")

autism_fX_res$ID2 <- autism_fX_res$ID
autism_fX_res <- autism_fX_res %>%
    separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")


## remove duplicate n

# remove duplicate, take first case
autism_cmX_res <- autism_cmX_res[!duplicated(autism_cmX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz and familyid variables only
autism_cmX_res <- autism_cmX_res %>% select(starts_with("v"), "familyid") # fdrsign = 113


# remove duplicate, take first case
autism_cfX_res <- autism_cfX_res[!duplicated(autism_cfX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz and familyid variables only
autism_cfX_res <- autism_cfX_res %>% select(starts_with("v"), "familyid") # fdrsign = 113


# remove duplicate, take first case
autism_mX_res <- autism_mX_res[!duplicated(autism_mX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz and familyid variables only
autism_mX_res <- autism_mX_res %>% select(starts_with("v"), "familyid") # fdrsign = 113


# remove duplicate, take first case
autism_fX_res <- autism_fX_res[!duplicated(autism_fX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz and familyid variables only
autism_fX_res <- autism_fX_res %>% select(starts_with("v"), "familyid") # fdrsign = 113



################ 
## Simplify interpretatin of chord. Get a set X for asd, mother, and father 052019 (this is only done in family cor, not in concordance below)
## corr matrix of ASD in mother and ASD in father won't be 100% same because the residual estimated before this intersection()

## check dim
autism_cmX_res %>% dim
autism_mX_res %>% dim

autism_cfX_res %>% dim
autism_fX_res %>% dim

## get intersect
tmp1 <- Reduce(intersect, list(names(autism_cmX_res), names(autism_mX_res), names(autism_cfX_res), names(autism_fX_res)))

autism_cmX_res <- autism_cmX_res %>% select(tmp1)
autism_mX_res <- autism_mX_res %>% select(tmp1)

autism_cfX_res <- autism_cfX_res %>% select(tmp1)
autism_fX_res <- autism_fX_res %>% select(tmp1)



# ## new 050719
# # they have 2 ASD in a family
# # table(autism_case$familyid)
# # autism_case$familyid %>% unique %>% length()
# # ASD n = 75, family n = 69
# #
# # join instead of average ASD in each family. think as having 2n in a family
# ## 
# 
# 
# 
# ### creaet mother table
# autism_mother <- autism_sign2 %>% filter(. , grepl("Mother", Role))
# 
# # remove duplicate, take first case
# autism_mother <- autism_mother[!duplicated(autism_mother$ID), ] # from 86 rows to 79
# 
# # pick mz and familyid variables only
# autism_mother <- autism_mother %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
# 
# 
# # table(autism_mother$familyid)
# # autism_mother$familyid %>% unique %>% length
# # mother n = 40, family n = 40
# 
# 
# ### check: some familyid only has parents wo ASD. e.g. 0050, 0061, 0062
# # autism_case2$familyid %>% table
# # autism_parents2$familyid %>% table()



# left join to autism-mother
# .x is asd, .y is mother
tmp1 <- left_join(autism_cmX_res, autism_mX_res, by = "familyid")
tmp2 <- tmp1[complete.cases(tmp1), ]
tmp2$familyid %>% unique %>% length()
# 35 asd matched with mother, from 31 family

# get matrix ready for corr (no family id)
autism_casemotherCor <- tmp2 %>% select(matches(".x"))
autism_motherCor <- tmp2 %>% select(matches(".y"), -"familyid")




# 
# ### creaet fahter table
# autism_father <- autism_sign2 %>% filter(. , grepl("Father", Role))
# 
# # remove duplicate, take first case
# autism_father <- autism_father[!duplicated(autism_father$ID), ] # from 86 rows to 79
# 
# # pick mz and familyid variables only
# autism_father <- autism_father %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
# 
# 
# # table(autism_father$familyid)
# # autism_father$familyid %>% unique %>% length
# # father n = 39, family n = 39
# 
# 
# ### check: some familyid only has parents wo ASD. e.g. 0050, 0061, 0062
# # autism_case2$familyid %>% table
# # autism_parents2$familyid %>% table()



# left join to autism-father
# .x is asd, .y is mother
tmp1 <- left_join(autism_cfX_res, autism_fX_res, by = "familyid")
tmp2 <- tmp1[complete.cases(tmp1), ]
tmp2$familyid %>% unique %>% length()
# 35 asd matched with mother, from 31 family

# get matrix ready for corr (no family id)
autism_casefatherCor <- tmp2 %>% select(matches(".x")) # should be same as autism_casemother2 because we have both parents matched
autism_fatherCor <- tmp2 %>% select(matches(".y"), -"familyid")


## *********** 
# B2. correlation 
library(psych)

## mother
cor_mother <- as.matrix(corr.test(autism_casemotherCor, autism_motherCor, method = "spearman")$r)
cor_motherIntra <- as.matrix(corr.test(autism_motherCor, method = "spearman")$r)
cor_mothercaseIntra <- as.matrix(corr.test(autism_casemotherCor, method = "spearman")$r)

# scale() won't change the spearman matrix
# log() has very minor effects (rounding mostly) to the spearman matrix

mother_r <- corr.test(autism_casemotherCor, autism_motherCor, method = "spearman")$r %>% diag
names(mother_r) <- names(autism_motherCor)


# DONT default p.adjust the whole matrix as I only need diagnoal then fdr
mother_p <- corr.test(autism_casemotherCor, autism_motherCor, method = "spearman", adjust = "none")$p %>% diag

names(mother_p) <- names(autism_motherCor)
mother_fdr <- p.adjust(mother_p, method='fdr')
filter(as.tibble(mother_fdr), value <=0.05) 
# SO? 2 mz from EWAS, after filter (x = 82) have sth to do with family


## ggplot
hist(mother_r)
p <- ggplot(as.tibble(mother_r), aes(x = mother_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
# positive more

outdirectory <- "results"
outfilename <- sprintf("%s_hist_ewas_sign_corr_mother_r.png", saveName)

# ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)





## father

cor_father <- as.matrix(corr.test(autism_casefatherCor, autism_fatherCor, method = "spearman")$r)
cor_fatherIntra <- as.matrix(corr.test(autism_fatherCor, method = "spearman")$r)
cor_fathercaseIntra <- as.matrix(corr.test(autism_casefatherCor, method = "spearman")$r)
# scale() won't change the spearman matrix
# log() has very minor effects (rounding mostly) to the spearman matrix

father_r <- corr.test(autism_casefatherCor, autism_fatherCor, method = "spearman")$r %>% diag
names(father_r) <- names(autism_fatherCor)


# DONT default p.adjust the whole matrix as I only need diagnoal then fdr
father_p <- corr.test(autism_casefatherCor, autism_fatherCor, method = "spearman", adjust = "none")$p %>% diag

names(father_p) <-  names(autism_fatherCor)
father_fdr <- p.adjust(father_p, method='fdr')
filter(as.tibble(father_fdr), value <=0.05) 
# SO? 1 mz from EWAS, after filter (x = 85) have sth to do with family


## ggplot
hist(father_r)
p <- ggplot(as.tibble(father_r), aes(x = father_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
# positive more

outdirectory <- "results"
outfilename <- sprintf("%s_hist_ewas_sign_corr_father_r.png", saveName)

# ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)



## violin plot

tmp1 <- mother_r %>% as.tibble
tmp1$Role <- "ASD-Mother"

tmp2 <- father_r %>% as.tibble
tmp2$Role <- "ASD-Father"

tmp3 <- bind_rows(tmp1, tmp2)

# plot (expert in one package is better)
p <- ggplot(data = tmp3, aes(x = Role, y = value)) 
p <- p + geom_violin(aes(fill = Role)) + geom_boxplot(width = 0.15)
p <- p + geom_jitter(shape=16, alpha = 0.8, size = 1.3, height = 0, width = 0.1)
p <- p + ylim(-0.8, 0.8)


outdirectory <- "results"
outfilename <- sprintf("%s_violin_asd_parents.pdf", saveName)

# ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)


# ******** -----
## S. save-rdata ----

outdirectory <- "results"
outfilename <- sprintf("%s_ewas_sign_corr.Rdata", saveName)
# outfilename <- "hilic_corr_sign.Rdata"
save(file=file.path(outdirectory,outfilename), 
     data_original, data_pheno,
     autism,
     autism_ccEwas, retFrame,
     fdrsign,
     mzmatch_table,
     cor_mother,
     cor_motherIntra,
     cor_mothercaseIntra,
     cor_father,
     cor_fatherIntra,
     cor_fathercaseIntra
)

# Load data
# load("results/hilic_ewas_sign_corr.Rdata") 

# ## for adhoc_cp calculation (back in early 2018)
# saveRDS(cor_asd, "results/cor_asd_163sign.rds")
# saveRDS(cor_parents, "results/cor_parents_163sign.rds")
# saveRDS(cor_family, "results/cor_family_163sign.rds")
# # hello <- readRDS("results/tablemz_sign_unique.rds")


# ******** -----
# B2. concordance 050819 ------------------------------------------------------------

## *********** 
## Data

## ASD and parent mzID
# for correlation variable labels (parents and autism)
mzmatch_table <- fdrsign
mzmatch_table$mzid_a <- sprintf("%s_a", mzmatch_table$mzid)
mzmatch_table$mzid_p <- sprintf("%s_p", mzmatch_table$mzid)


## autism df pick only EWAS sign mz
# name of meta in autism df
tmpMetaName <- autism %>% select(-starts_with("v")) %>% colnames # 11

# autism df with sign mz and meta col only
autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))


# check. from corr.test warning. 050819
#' Y needs at need 2 points to get the regression and predicted/residual
# autism_sign %>% filter(Role == "Mother") %>% select(v6856)
# autism_sign %>% filter(Role == "Mother") %>% select(v6870)
# autism_sign %>% filter(Role == "Mother") %>% select(v7103)
# autism_sign %>% filter(Role == "Father") %>% select(v8255)
# autism_sign %>% filter(Role == "Mother") %>% select(v8824)

l <- list() # set tmp

## impute missing age, by role 
# (they are diff age group, not a single distribution)
l[["a"]] <- autism_sign %>% filter(. , grepl("Father|Mother", Role)) %>% select(Age)
l[["parentsub"]] <- Hmisc::impute(l[["a"]][[1]]) %>% as.numeric #median impute
autism_sign$Age[which(autism_sign$Role %in% c("Father", "Mother"))] <- l[["parentsub"]]


l[["b"]] <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
l[["asdsub"]] <- Hmisc::impute(l[["b"]][[1]]) %>% as.numeric #median impute
autism_sign$Age[which(autism_sign$Role %in% c("Proband"))] <- l[["asdsub"]]

l[["b"]] <- autism_sign %>% filter(. , grepl("Control", Role)) %>% select(Age)
l[["ctsub"]] <- Hmisc::impute(l[["b"]][[1]]) %>% as.numeric #median impute
autism_sign$Age[which(autism_sign$Role %in% c("Control"))] <- l[["ctsub"]]


## reference EST corr, create dataframe for each 
autism_sign_case <- autism_sign %>% filter(Role %in% c("Proband"))
autism_sign_ct <- autism_sign %>% filter(Role %in% c("Control"))

autism_sign_mother <- autism_sign %>% filter(Role %in% c("Mother"))
autism_sign_father <- autism_sign %>% filter(Role %in% c("Father"))

## filter X with min 20% data in each df

# function: percentage of "0"
percentzero <- function (x) {
    100*mean(!(x))
}

filterDF <- function(autism_sign_case) {
    loop_mz <- autism_sign_case %>% select(starts_with("v")) %>% colnames
    meta <- autism_sign_case %>% select(-starts_with("v"))
    
    tmp1 <- autism_sign_case %>% select(loop_mz)
    tmp2 <- sapply(tmp1, percentzero)   # sum(autism_n$v1 %in% 0) / nrow(autism_n) * 100
    tmp3 <- names(tmp2)[which(tmp2 < 80)] # find names < 80% "0", i.e., keep mzid names
    tmp4 <- tmp1[, tmp3] # df with needed variables
    tmp5 <- bind_cols(tmp4, meta)
}

autism_sign_case1 <- filterDF(autism_sign_case) 
dim(autism_sign_case1) # 82,119 with meta
autism_sign_ct1 <- filterDF(autism_sign_ct)
dim(autism_sign_ct1) # 29, 118

autism_sign_mother1 <- filterDF(autism_sign_mother)
dim(autism_sign_mother1) # 42, 96
autism_sign_father1 <- filterDF(autism_sign_father)
dim(autism_sign_father1) # 44, 98



## get mutual Xs across pair
# asd and ct
tmp1 <- intersect(names(autism_sign_case1), names(autism_sign_ct1))
autism_asdX <- autism_sign_case1 %>% select(tmp1)
autism_ctX <- autism_sign_ct1 %>% select(tmp1)

# asd and mother
tmp1 <- intersect(names(autism_sign_case1), names(autism_sign_mother1))
autism_cmX <- autism_sign_case1 %>% select(tmp1)
autism_mX <- autism_sign_mother1 %>% select(tmp1)


# asd and father
tmp1 <- intersect(names(autism_sign_case1), names(autism_sign_father1))
autism_cfX <- autism_sign_case1 %>% select(tmp1)
autism_fX <- autism_sign_father1 %>% select(tmp1)


# get residual
# stratify instead of get residual from a single model 050719

library(broom)

adjFunASD <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_asdX, Role == "Proband")) %>% augment()
    ret[[".resid"]]
}


adjFunCt <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_ctX, Role == "Control")) %>% augment()
    ret[[".resid"]]
}



adjFunASDMother <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_cmX, Role == "Proband")) %>% augment()
    ret[[".resid"]]
}


adjFunASDFather <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_cfX, Role == "Proband")) %>% augment()
    ret[[".resid"]]
}



adjFunMother <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_mX, Role == "Mother")) %>% augment()
    ret[[".resid"]]
}

adjFunFather <- function (mz) {
    setform <- sprintf("log2(I(%s+1)) ~ Age + Batch", mz) #nature ln = log
    ret <- lm(setform, filter(autism_fX, Role == "Father")) %>% augment()
    ret[[".resid"]]
}



# asd
loop_mz <- autism_asdX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunASD(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_asdX %>% select(-starts_with("v"))

autism_asdX_res <- cbind(residual, nonmz)



# ct
loop_mz <- autism_ctX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunCt(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_ctX %>% select(-starts_with("v"))

autism_ctX_res <- cbind(residual, nonmz)




# asd mohter
loop_mz <- autism_cmX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunASDMother(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_cmX %>% select(-starts_with("v"))

autism_cmX_res <- cbind(residual, nonmz)


# mother
loop_mz <- autism_mX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunMother(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_mX %>% select(-starts_with("v"))

autism_mX_res <- cbind(residual, nonmz)


# asd father
loop_mz <- autism_cfX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunASDFather(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_cfX %>% select(-starts_with("v"))

autism_cfX_res <- cbind(residual, nonmz)


# fahter
loop_mz <- autism_fX %>% select(starts_with("v")) %>% colnames

residual <- loop_mz %>% map(safely(~adjFunFather(.))) %>% setNames(loop_mz)

residual <- residual %>% transpose() %>% .[["result"]]

residual <- as.data.frame(residual)

nonmz <- autism_fX %>% select(-starts_with("v"))

autism_fX_res <- cbind(residual, nonmz)





## ASD remove duplicate, take first case
autism_asdX_res <- autism_asdX_res[!duplicated(autism_asdX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_asdX_res <- autism_asdX_res %>% select(starts_with("v")) # fdrsign = 113


## CT remove duplicate, take first case
autism_ctX_res <- autism_ctX_res[!duplicated(autism_ctX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_ctX_res <- autism_ctX_res %>% select(starts_with("v")) # fdrsign = 113



## asd mother
# remove duplicate, take first case
autism_cmX_res <- autism_cmX_res[!duplicated(autism_cmX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_cmX_res <- autism_cmX_res %>% select(starts_with("v")) # fdrsign = 113


## asd father
# remove duplicate, take first case
autism_cfX_res <- autism_cfX_res[!duplicated(autism_cfX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_cfX_res <- autism_cfX_res %>% select(starts_with("v")) # fdrsign = 113


## mohter
# remove duplicate, take first case
autism_mX_res <- autism_mX_res[!duplicated(autism_mX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_mX_res <- autism_mX_res %>% select(starts_with("v")) # fdrsign = 113


## father
# remove duplicate, take first case
autism_fX_res <- autism_fX_res[!duplicated(autism_fX_res$ID), ] #29 controls + 75 cases = 104 total
# pick mz
autism_fX_res <- autism_fX_res %>% select(starts_with("v")) # fdrsign = 113






## correlation
library(psych)
cor_asd <- as.matrix(corr.test(autism_asdX_res, method = "spearman")$r)

cor_ct <- as.matrix(corr.test(autism_ctX_res, method = "spearman")$r)

# asdmotther
cor_cmX <- as.matrix(corr.test(autism_cmX_res, method = "spearman")$r)

cor_mX <- as.matrix(corr.test(autism_mX_res, method = "spearman")$r)

# asdfather
cor_cfX <- as.matrix(corr.test(autism_cfX_res, method = "spearman")$r)

cor_fX <- as.matrix(corr.test(autism_fX_res, method = "spearman")$r)




## correlation matrx to dataframe




prefConcordDF <- function(cor_asd, cor_ct) {
  ## ASD
  
  indASD <- which(upper.tri(cor_asd, diag = F) , arr.ind = T)  # arr.ind: return array indices
  
  spearASD <- data.frame(col = dimnames(cor_asd)[[2]][indASD[,2]],  # picking the selected (with repeat) from the list, NOT extract
                       row = dimnames(cor_asd)[[1]][indASD[,1]],
                       val = cor_asd[indASD])  # Extracting the df with the indF matrix's coordinated pair into a vector
  
  # Create pair ID
  spearASD$col_row_pair <- sprintf("%s-%s", spearASD$col, spearASD$row)
  
  
  ## CT
  
  indCT <- which(upper.tri(cor_ct, diag = F) , arr.ind = T)  # arr.ind: return array indices
  
  spearCT <- data.frame(col = dimnames(cor_ct)[[2]][indCT[,2]],  # picking the selected (with repeat) from the list, NOT extract
                         row = dimnames(cor_ct)[[1]][indCT[,1]],
                         val = cor_ct[indCT])  # Extracting the df with the indF matrix's coordinated pair into a vector
  
  # Create pair ID
  spearCT$col_row_pair <- sprintf("%s-%s", spearCT$col, spearCT$row)
  
  
  
  # Left join
  spearASD <- left_join(spearASD, spearCT, by = "col_row_pair") #
  spearASD <- spearASD %>% filter(!is.na(val.y)) # remove NA # val.x is asd, val.y is ct
  return(spearASD)
}

asdct <- prefConcordDF(cor_asd, cor_ct)
asdmother <- prefConcordDF(cor_cmX, cor_mX)
asdfather <- prefConcordDF(cor_cfX, cor_fX)



## concordance

# pearson corrlelation
corr.test(asdct$val.x, asdct$val.y, method = "pearson")$r
# 0.2468878

corr.test(asdmother$val.x, asdmother$val.y, method = "pearson")$r
#  0.3576137

corr.test(asdfather$val.x, asdfather$val.y, method = "pearson")$r
#  0.2880821


## plot

# asd vs ct
library(ggplot2)
quartz(height=6, width=8)

p <- ggplot(data = asdct, aes(x = val.x, y = val.y))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("ASD, spearman corr coef") + 
    ylab("CT, spearman corr coef") +
    ggtitle("Concordance (pearson)")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/concordance_r_asd_ct.png", scale=1, dpi=400)


# asd vs mother
library(ggplot2)
quartz(height=6, width=8)

p <- ggplot(data = asdmother, aes(x = val.x, y = val.y))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("ASD, spearman corr coef") + 
    ylab("Mother, spearman corr coef") +
    ggtitle("Concordance (pearson)")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/concordance_r_asd_mother.png", scale=1, dpi=400)



# asd vs father
library(ggplot2)
quartz(height=6, width=8)

p <- ggplot(data = asdfather, aes(x = val.x, y = val.y))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("ASD, spearman corr coef") + 
    ylab("Father, spearman corr coef") +
    ggtitle("Concordance (pearson)")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/concordance_r_asd_father.png", scale=1, dpi=400)




## -----
## ***************
## Testing if family is using the family matched and asd-ct unmatch as reference
## why: did not match family ASD-mother is useless to say the "family" effect. it's similar to permutation if no family match.
## 061419 
## ***************

## rename corr matrxi
tmp1 <- colnames(cor_motherIntra)
colnames(cor_motherIntra) <- str_replace(tmp1, ".y$", "")
tmp1 <- rownames(cor_motherIntra)
rownames(cor_motherIntra) <- str_replace(tmp1, ".y$", "")


tmp1 <- colnames(cor_mothercaseIntra)
colnames(cor_mothercaseIntra) <- str_replace(tmp1, ".x$", "")
tmp1 <- rownames(cor_mothercaseIntra)
rownames(cor_mothercaseIntra) <- str_replace(tmp1, ".x$", "")


tmp1 <- colnames(cor_fatherIntra)
colnames(cor_fatherIntra) <- str_replace(tmp1, ".y$", "")
tmp1 <- rownames(cor_fatherIntra)
rownames(cor_fatherIntra) <- str_replace(tmp1, ".y$", "")


tmp1 <- colnames(cor_fathercaseIntra)
colnames(cor_fathercaseIntra) <- str_replace(tmp1, ".x$", "")
tmp1 <- rownames(cor_fathercaseIntra)
rownames(cor_fathercaseIntra) <- str_replace(tmp1, ".x$", "")



## create pairs
asdct <- prefConcordDF(cor_asd, cor_ct)
# 1711 obs
asdmother <- prefConcordDF(cor_motherIntra, cor_mothercaseIntra)
# 946 obs
asdfather <- prefConcordDF(cor_fatherIntra, cor_fathercaseIntra)
# 946 obs



## concordance

# pearson corrlelation, 5050 obs
corr.test(asdct$val.x, asdct$val.y, method = "pearson")$r
corr.test(asdct$val.x, asdct$val.y, method = "pearson")$ci
#  0.246887

corr.test(asdmother$val.x, asdmother$val.y, method = "pearson")$r
corr.test(asdmother$val.x, asdmother$val.y, method = "pearson")$ci
#  0.4284781

corr.test(asdfather$val.x, asdfather$val.y, method = "pearson")$r
corr.test(asdfather$val.x, asdfather$val.y, method = "pearson")$ci
#  0.3098309

## ----------



# ******** -----

#' 
#' # ******** -----
#' # B2. family old 050719 ------------------------------------------------------------
#' # 
#' # ## *********** 
#' # ## Data
#' # 
#' # ## ASD and parent mzID
#' # # for correlation variable labels (parents and autism)
#' # mzmatch_table <- fdrsign
#' # mzmatch_table$mzid_a <- sprintf("%s_a", mzmatch_table$mzid)
#' # mzmatch_table$mzid_p <- sprintf("%s_p", mzmatch_table$mzid)
#' # 
#' # 
#' # ## autism df pick only EWAS sign mz
#' # # name of meta in autism df
#' # tmpMetaName <- autism %>% select(-starts_with("v")) %>% colnames # 11
#' # 
#' # # autism df with sign mz and meta col only
#' # autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))
#' # 
#' # 
#' # 
#' # 
#' # 
#' # ## *********** 
#' # # B1. adjust and residual prep
#' # 
#' # # adjust age and sex; switch the whole section on/off
#' # 
#' # l <- list() # set tmp
#' # 
#' # # impute missing age, by role
#' # l[["a"]] <- autism_sign %>% filter(. , grepl("Father|Mother", Role)) %>% select(Age)
#' # l[["parentsub"]] <- Hmisc::impute(l[["a"]][[1]]) %>% as.numeric #median impute
#' # autism_sign$Age[which(autism_sign$Role %in% c("Father", "Mother"))] <- l[["parentsub"]]
#' # 
#' # l[["b"]] <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
#' # l[["asdsub"]] <- Hmisc::impute(l[["b"]][[1]]) %>% as.numeric #median impute
#' # autism_sign$Age[which(autism_sign$Role %in% c("Proband"))] <- l[["asdsub"]]
#' # 
#' # 
#' # # get residual
#' # ## for both autism kid and parents
#' # # in a single model (father mother proband as one pop) instead of stratify 050719
#' # 
#' # 
#' # # take only ASD and parents
#' # autism_sign2 <- autism_sign %>%  filter(Role %in% c("Father", "Mother", "Proband"))
#' # 
#' # library(broom)
#' # 
#' # loop_mz <- autism_sign2 %>% select(starts_with("v")) %>% colnames
#' # 
#' # adjFun <- function (mz) {
#' #     setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
#' #     ret <- lm(setform, autism_sign2) %>% augment()
#' #     ret[[".resid"]]
#' # }
#' # 
#' # residual <- loop_mz %>% map(safely(~adjFun(.))) %>% setNames(loop_mz)
#' # 
#' # residual <- residual %>% transpose() %>% .[["result"]]
#' # 
#' # residual <- as.data.frame(residual)
#' # 
#' # nonmz <- autism_sign2 %>% select(-starts_with("v"))
#' # 
#' # autism_sign2 <- cbind(residual, nonmz)
#' # 
#' # 
#' # ## *********** 
#' # # B2. prep table for correlation after fixing lots of "0" issue
#' # 
#' # # filter for rows (family)
#' # autism_sign2 <- autism_sign2 %>% filter(. , grepl("Proband|Father|Mother", Role))
#' # 
#' # 
#' # # create familyID from BCH "ID"
#' # autism_sign2$ID2 <- autism_sign2$ID
#' # 
#' # autism_sign2 <- autism_sign2 %>%
#' #     separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")
#' # 
#' # 
#' # ## case and parents tables, summerize row by family ID
#' # 
#' # ### creaet case table
#' # autism_case <- autism_sign2 %>% filter(. , grepl("Proband", Role))
#' # 
#' # # remove duplicate, take first case
#' # autism_case <- autism_case[!duplicated(autism_case$ID), ] #29 controls + 75 cases = 104 total
#' # 
#' # # pick mz and familyid variables only
#' # autism_case <- autism_case %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
#' # 
#' # # tmp1 <- autism_case[, 100:126]
#' # 
#' # 
#' # ## new 050719
#' # # they have 2 ASD in a family
#' # # table(autism_case$familyid)
#' # # autism_case$familyid %>% unique %>% length()
#' # # ASD n = 75, family n = 69
#' # #
#' # # join instead of average ASD in each family. think as having 2n in a family
#' # ## 
#' # 
#' # 
#' # 
#' # ### creaet mother table
#' # autism_mother <- autism_sign2 %>% filter(. , grepl("Mother", Role))
#' # 
#' # # remove duplicate, take first case
#' # autism_mother <- autism_mother[!duplicated(autism_mother$ID), ] # from 86 rows to 79
#' # 
#' # # pick mz and familyid variables only
#' # autism_mother <- autism_mother %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
#' # 
#' # 
#' # # table(autism_mother$familyid)
#' # # autism_mother$familyid %>% unique %>% length
#' # # mother n = 40, family n = 40
#' # 
#' # 
#' # ### check: some familyid only has parents wo ASD. e.g. 0050, 0061, 0062
#' # # autism_case2$familyid %>% table
#' # # autism_parents2$familyid %>% table()
#' # 
#' # 
#' # 
#' # # left join to autism
#' # # .x is asd, .y is mother
#' # tmp1 <- left_join(autism_case, autism_mother, by = "familyid")
#' # tmp2 <- tmp1[complete.cases(tmp1), ]
#' # tmp2$familyid %>% unique %>% length()
#' # # 35 asd matched with mother, from 31 family
#' # 
#' # # get matrix ready for corr (no family id)
#' # autism_casemother2 <- tmp2 %>% select(matches(".x"))
#' # autism_mother2 <- tmp2 %>% select(matches(".y"), -"familyid")
#' # 
#' # 
#' # 
#' # 
#' # 
#' # ### creaet fahter table
#' # autism_father <- autism_sign2 %>% filter(. , grepl("Father", Role))
#' # 
#' # # remove duplicate, take first case
#' # autism_father <- autism_father[!duplicated(autism_father$ID), ] # from 86 rows to 79
#' # 
#' # # pick mz and familyid variables only
#' # autism_father <- autism_father %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
#' # 
#' # 
#' # # table(autism_father$familyid)
#' # # autism_father$familyid %>% unique %>% length
#' # # father n = 39, family n = 39
#' # 
#' # 
#' # ### check: some familyid only has parents wo ASD. e.g. 0050, 0061, 0062
#' # # autism_case2$familyid %>% table
#' # # autism_parents2$familyid %>% table()
#' # 
#' # 
#' # 
#' # # left join to autism
#' # # .x is asd, .y is mother
#' # tmp1 <- left_join(autism_case, autism_father, by = "familyid")
#' # tmp2 <- tmp1[complete.cases(tmp1), ]
#' # tmp2$familyid %>% unique %>% length()
#' # # 35 asd matched with mother, from 31 family
#' # 
#' # # get matrix ready for corr (no family id)
#' # autism_casefather2 <- tmp2 %>% select(matches(".x")) # should be same as autism_casemother2 because we have both parents matched
#' # autism_father2 <- tmp2 %>% select(matches(".y"), -"familyid")
#' # 
#' # 
#' # 
#' # 
#' # ## *********** 
#' # # B2. correlation 
#' # library(psych)
#' # 
#' # ## mother
#' # 
#' # cor_mother <- as.matrix(corr.test(autism_casemother2, autism_mother2, method = "spearman")$r)
#' # # scale() won't change the spearman matrix
#' # # log() has very minor effects (rounding mostly) to the spearman matrix
#' # 
#' # mother_r <- corr.test(autism_casemother2, autism_mother2, method = "spearman")$r %>% diag
#' # names(mother_r) <- mzmatch_table$mzid
#' # 
#' # 
#' # # DONT default p.adjust the whole matrix as I only need diagnoal then fdr
#' # mother_p <- corr.test(autism_casemother2, autism_mother2, method = "spearman", adjust = "none")$p %>% diag
#' # 
#' # names(mother_p) <- mzmatch_table$mzid
#' # mother_fdr <- p.adjust(mother_p, method='fdr')
#' # filter(as.tibble(mother_fdr), value <=0.05) 
#' # # SO? 33 mz from EWAS have sth to do with family
#' # 
#' # 
#' # ## ggplot
#' # hist(mother_r)
#' # p <- ggplot(as.tibble(mother_r), aes(x = mother_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
#' # 
#' # outdirectory <- "results"
#' # outfilename <- sprintf("%s_hist_ewas_sign_corr_mother_r.png", saveName)
#' # 
#' # ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)
#' # 
#' # 
#' # 
#' # 
#' # 
#' # ## father
#' # 
#' # cor_father <- as.matrix(corr.test(autism_casefather2, autism_father2, method = "spearman")$r)
#' # # scale() won't change the spearman matrix
#' # # log() has very minor effects (rounding mostly) to the spearman matrix
#' # 
#' # father_r <- corr.test(autism_casefather2, autism_father2, method = "spearman")$r %>% diag
#' # names(father_r) <- mzmatch_table$mzid
#' # 
#' # 
#' # # DONT default p.adjust the whole matrix as I only need diagnoal then fdr
#' # father_p <- corr.test(autism_casefather2, autism_father2, method = "spearman", adjust = "none")$p %>% diag
#' # 
#' # names(father_p) <- mzmatch_table$mzid
#' # father_fdr <- p.adjust(father_p, method='fdr')
#' # filter(as.tibble(father_fdr), value <=0.05) 
#' # # SO? 40 mz from EWAS have sth to do with family
#' # 
#' # 
#' # ## ggplot
#' # hist(father_r)
#' # p <- ggplot(as.tibble(father_r), aes(x = father_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
#' # 
#' # outdirectory <- "results"
#' # outfilename <- sprintf("%s_hist_ewas_sign_corr_father_r.png", saveName)
#' # 
#' # ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)
#' 
#' 
#' 
#' # ******** -----
#' # B3. family old ------------------------------------------------------------
#' 
#' ## *********** 
#' ## Data
#' 
#' ## ASD and parent mzID
#' # for correlation variable labels (parents and autism)
#' mzmatch_table <- fdrsign
#' mzmatch_table$mzid_a <- sprintf("%s_a", mzmatch_table$mzid)
#' mzmatch_table$mzid_p <- sprintf("%s_p", mzmatch_table$mzid)
#' 
#' 
#' ## autism df pick only EWAS sign mz
#' # name of meta in autism df
#' tmpMetaName <- autism %>% select(-starts_with("v")) %>% colnames # 11
#' 
#' # autism df with sign mz and meta col only
#' autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))
#' 
#' 
#' 
#' 
#' 
#' ## *********** 
#' # B1. adjust and residual prep
#' 
#' # adjust age and sex; switch the whole section on/off
#' 
#' l <- list() # set tmp
#' 
#' # impute missing age, by role
#' l[["a"]] <- autism_sign %>% filter(. , grepl("Father|Mother", Role)) %>% select(Age)
#' l[["parentsub"]] <- Hmisc::impute(l[["a"]][[1]]) %>% as.numeric #median impute
#' autism_sign$Age[which(autism_sign$Role %in% c("Father", "Mother"))] <- l[["parentsub"]]
#' 
#' l[["b"]] <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
#' l[["asdsub"]] <- Hmisc::impute(l[["b"]][[1]]) %>% as.numeric #median impute
#' autism_sign$Age[which(autism_sign$Role %in% c("Proband"))] <- l[["asdsub"]]
#' 
#' 
#' # get residual
#' ## for both autism kid and parents
#' 
#' library(broom)
#' 
#' loop_mz <- autism_sign %>% select(starts_with("v")) %>% colnames
#' 
#' adjFun <- function (mz) {
#'     setform <- sprintf("log2(I(%s+1)) ~ Age + Sex + Batch", mz) #nature ln = log
#'     ret <- lm(setform, autism_sign) %>% augment()
#'     ret[[".resid"]]
#' }
#' 
#' residual <- loop_mz %>% map(safely(~adjFun(.))) %>% setNames(loop_mz)
#' 
#' residual <- residual %>% transpose() %>% .[["result"]]
#' 
#' residual <- as.data.frame(residual)
#' 
#' nonmz <- autism_sign %>% select(-starts_with("v"))
#' 
#' autism_sign <- cbind(residual, nonmz)
#' 
#' 
#' ## *********** 
#' # B2. prep table for correlation after fixing lots of "0" issue
#' 
#' # filter for rows (family)
#' autism_sign <- autism_sign %>% filter(. , grepl("Proband|Father|Mother", Role))
#' 
#' # create familyID from BCH "ID"
#' autism_sign$ID2 <- autism_sign$ID
#' 
#' autism_sign <- autism_sign %>%
#'     separate(ID2, into = c("dump1", "familyid", "dump2"), sep = "-")
#' 
#' 
#' ## case and parent tables, summerize row by family ID
#' 
#' ### creaet case table
#' autism_case <- autism_sign %>% filter(. , grepl("Proband", Role))
#' 
#' # remove duplicate, take first case
#' autism_case <- autism_case[!duplicated(autism_case$ID), ] #29 controls + 75 cases = 104 total
#' 
#' # pick mz and familyid variables only
#' autism_case <- autism_case %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
#' 
#' ### get unique rows by familyid
#' autism_case2 <- autism_case %>% group_by(familyid) %>% summarise_all(funs(mean)) # they may have 2 asd in a family, I average them for now
#'     # row = 69
#' 
#' ### check:OK 
#'     # > View(autism_case2[, c(1, 194)])
#'     # > View(autism_case2[, c(1,2)])
#' 
#' 
#' 
#' ### creaet parent table
#' autism_parents <- autism_sign %>% filter(. , grepl("Father|Mother", Role))
#' 
#' # remove duplicate, take first case
#' autism_parents <- autism_parents[!duplicated(autism_parents$ID), ] # from 86 rows to 79
#' 
#' # pick mz and familyid variables only
#' autism_parents <- autism_parents %>% select(c(mzmatch_table$mzid, "familyid")) # fdrsign = 113
#' 
#' ### get unique rows by familyid
#' autism_parents2 <- autism_parents %>% group_by(familyid) %>% summarise_all(funs(mean)) # ☐ they may have 2 asd in a family, I average them for now
#'     # row = 43
#' 
#' ### check: some familyid only has parents wo ASD. e.g. 0050, 0061, 0062
#' # autism_case2$familyid %>% table
#' # autism_parents2$familyid %>% table()
#' 
#' 
#' 
#' ## row (familyid) match between asd and parent dfs
#' #' find the mutual set row
#' tmp_familyid <- intersect(autism_case2$familyid, autism_parents2$familyid)
#' tmp_familyid %>% length # 33
#' 
#' # extract mutual row match dataset
#' # case
#' autism_case2 <- autism_case2 %>% filter(grepl(paste(tmp_familyid, collapse="|"), familyid))
#' 
#' # parent
#' autism_parents2 <- autism_parents2 %>% filter(grepl(paste(tmp_familyid, collapse="|"), familyid))
#' 
#' # same order?
#' all(autism_case2$familyid == autism_parents2$familyid) # OK, sorted
#' # all(c(1,3,2) == c(1,2,3))
#' 
#' ## final dataset ready for corr (drop family ID variable)
#' autism_case2 <- autism_case2 %>% select(-familyid)
#' autism_parents2 <- autism_parents2 %>% select(-familyid)
#' 
#' ### rename mz # so corr matrix label is clear from autism or parents
#' # case
#' # autism_case3 <-  autism_case3 %>% select(starts_with("v")) # corr take numeric only
#' tmp_index_mzcase <- match(colnames(autism_case2), mzmatch_table$mzid)
#' colnames(autism_case2) <- mzmatch_table$mzid_a[tmp_index_mzcase]
#' 
#' 
#' # parents
#' # autism_parents3 <-  autism_parents3 %>% select(starts_with("v")) # corr take numeric only
#' tmp_index_mzparents <- match(colnames(autism_parents2), mzmatch_table$mzid)
#' colnames(autism_parents2) <- mzmatch_table$mzid_p[tmp_index_mzparents]
#' 
#' 
#' # # add familyID meta if needed
#' # # autism_case3 <- bind_cols(familyid = autism_case2$familyid, autism_case3) #cbind default char = factor
#' # # autism_parents3 <- bind_cols(familyid = autism_parents2$familyid, autism_parents3)
#' # 
#' # # same order?
#' # # all(autism_case3$familyid == autism_parents3$familyid) # OK, sorted
#' # # all(c(1,3,2) == c(1,2,3))
#' 
#' 
#' outdirectory <- "results"
#' outfilename <- sprintf("corr_family_case2_%s.csv", saveName)
#' write.csv(autism_case2, file = file.path(outdirectory,outfilename))
#' 
#' outfilename <- sprintf("corr_family_parents2_%s.csv", saveName)
#' write.csv(autism_parents2, file = file.path(outdirectory,outfilename))
#' 
#' 
#' 
#' ## *********** 
#' # B2. correlation 
#' library(psych)
#' 
#' cor_asd <- as.matrix(corr.test(autism_case2, method = "spearman")$r)
#' 
#' cor_parents <- as.matrix(corr.test(autism_parents2, method = "spearman")$r)
#' 
#' cor_family <- as.matrix(corr.test(autism_case2, autism_parents2, method = "spearman")$r)
#'     # scale() won't change the spearman matrix
#'     # log() has very minor effects (rounding mostly) to the spearman matrix
#' 
#' family_r <- corr.test(autism_case2, autism_parents2, method = "spearman")$r %>% diag
#' names(family_r) <- mzmatch_table$mzid
#' 
#' 
#' # DONT default p.adjust as I only need diagnoal then fdr
#' family_p <- corr.test(autism_case2, autism_parents2, method = "spearman", adjust = "none")$p %>% diag
#' 
#' names(family_p) <- mzmatch_table$mzid
#' family_fdr <- p.adjust(family_p, method='fdr')
#' filter(as.tibble(family_fdr), value <=0.05) 
#' # SO? 26 mz from EWAS have sth to do with family
#' 
#' ## ggplot
#' hist(family_r)
#' p <- ggplot(as.tibble(family_r), aes(x = family_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
#' 
#' outdirectory <- "results"
#' outfilename <- sprintf("%s_hist_ewas_sign_corr_family_r.png", saveName)
#' 
#' ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)
#' 
#' 
#' # ******** -----
#' # C. age ------------------------------------------------------------
#' 
#' ## prep basic df
#' # autism df with sign mz and meta col only
#' autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))
#' 
#' # select autism only
#' autism_sign <- autism_sign %>% filter(Role == "Proband")
#' 
#' # remove duplicate, take first autism case
#' autism_sign <- autism_sign[!duplicated(autism_sign$ID), ] #29 controls + 75 cases = 104 total
#' 
#' # ******** -----
#' ## adjust & get residual
#' 
#' # impute missing age, by role
#' tmp1 <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
#' tmp2 <- Hmisc::impute(tmp1$Age) %>% as.numeric #median impute
#' autism_sign$Age <- tmp2
#' 
#' # adjust
#' library(broom)
#' 
#' loop_mz <- autism_sign %>% select(starts_with("v")) %>% colnames
#' 
#' adjFun <- function (mz) {
#'     setform <- sprintf("log2(I(%s+1)) ~ Sex + Batch", mz) #nature ln = log
#'     ret <- lm(setform, autism_sign) %>% augment()
#'     ret[[".resid"]]
#' }
#' 
#' residual <- loop_mz %>% map(safely(~adjFun(.))) %>% setNames(loop_mz)
#' 
#' residual <- residual %>% transpose() %>% .[["result"]]
#' 
#' residual <- as.data.frame(residual)
#' 
#' nonmz <- autism_sign %>% select(-starts_with("v"))
#' 
#' autism_sign <- cbind(residual, nonmz)
#' 
#' 
#' ## extract only mz for correlation
#' autism_sign_mz <- autism_sign %>% select(starts_with("v"))
#' autism_sign_age <- autism_sign %>% select(Age)
#' 
#' ## *********** 
#' # B2. correlation 
#' library(psych)
#' 
#' age_r <- corr.test(autism_sign_mz, autism_sign_age, , method = "spearman")$r
#' 
#' # use adjust = "none" for raw p
#' age_ci <- corr.test(autism_sign_mz, autism_sign_age, , method = "spearman", adjust = "fdr")$ci
#' 
#' age_p <- as.tibble(corr.test(autism_sign_mz, autism_sign_age, , method = "spearman", adjust = "fdr")$p)
#' age_p$mzid <- mzmatch_table$mzid
#' filter(age_p, Age <=0.05) # non pass fdr sign
#' 
#' 
#' ## ggplot
#' hist(age_r)
#' p <- ggplot(as.tibble(age_r), aes(x = age_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
#' 
#' outdirectory <- "results"
#' outfilename <- sprintf("%s_hist_ewas_sign_corr_age_r.png", saveName)
#' 
#' ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)
#' 
#' 
#' # ******** -----
#' # D. sex ------------------------------------------------------------
#' 
#' ## prep basic df
#' # autism df with sign mz and meta col only
#' autism_sign <- autism %>% select(one_of(mzmatch_table$mzid, tmpMetaName))
#' 
#' # select autism only
#' autism_sign <- autism_sign %>% filter(Role == "Proband")
#' 
#' # remove duplicate, take first autism case
#' autism_sign <- autism_sign[!duplicated(autism_sign$ID), ] #29 controls + 75 cases = 104 total
#' 
#' # ******** -----
#' ## adjust & get residual
#' 
#' # impute missing age, by role
#' tmp1 <- autism_sign %>% filter(. , grepl("Proband", Role)) %>% select(Age)
#' tmp2 <- Hmisc::impute(tmp1$Age) %>% as.numeric #median impute
#' autism_sign$Age <- tmp2
#' 
#' # adjust
#' library(broom)
#' 
#' loop_mz <- autism_sign %>% select(starts_with("v")) %>% colnames
#' 
#' adjFun <- function (mz) {
#'     setform <- sprintf("log2(I(%s+1)) ~ Age + Batch", mz) #nature ln = log
#'     ret <- lm(setform, autism_sign) %>% augment()
#'     ret[[".resid"]]
#' }
#' 
#' residual <- loop_mz %>% map(safely(~adjFun(.))) %>% setNames(loop_mz)
#' 
#' residual <- residual %>% transpose() %>% .[["result"]]
#' 
#' residual <- as.data.frame(residual)
#' 
#' nonmz <- autism_sign %>% select(-starts_with("v"))
#' 
#' autism_sign <- cbind(residual, nonmz)
#' 
#' 
#' ## extract only mz for correlation
#' autism_sign_mz <- autism_sign %>% select(starts_with("v"))
#' autism_sign_sex <- autism_sign %>% select(Sex)
#' tmp1 <- c("Female" = 0, "Male" = 1)
#' autism_sign_sex <- tmp1[autism_sign_sex$Sex]
#' 
#' 
#' ## *********** 
#' # B2. correlation 
#' 
#' 
#' ## *********** 
#' # Now use pearson to approximate point biserial
#' #
#' #' ref: https://stackoverflow.com/questions/35880910/point-biserial-and-p-value
#' # # library(psych)
#' # # biserial(autism_sign_mz, autism_sign_sex) #no p values
#' # 
#' # # ref: https://stackoverflow.com/questions/16281667/p-value-for-polyserial-correlation
#' # # ref: https://stackoverflow.com/questions/35880910/point-biserial-and-p-value
#' # library(polycor)
#' # polyserial(autism_sign_mz$v21, autism_sign_sex$Sex, ML = TRUE, std.err = TRUE, maxcor=.9999, bins=4)
#' # 
#' # biSerialFunc <- function(autism_sign_mz) {
#' #     polyserial(autism_sign_mz, autism_sign_sex$Sex, ML = TRUE, std.err = TRUE, maxcor=.9999, bins=4)
#' # }
#' # 
#' # 
#' # tmp1 <- map(autism_sign_mz, biSerialFunc)
#' # map(tmp1, 2)
#' # map(tmp1, `[`, c("rho", "n"))
#' # map_dfr(tmp1, `[`, c("rho", "n"))
#' ## *********** 
#' 
#' 
#' 
#' #  pearson to approximate point biserial (Y = continous, X is binary normial)
#' # cor.test(autism_sign_mz$v21, autism_sign_sex) # test
#' 
#' sex_r <- corr.test(autism_sign_mz, autism_sign_sex, , method = "pearson")$r
#' 
#' sex_ci <- corr.test(autism_sign_mz, autism_sign_sex, , method = "pearson", adjust = "fdr")$ci
#' 
#' sex_p <- as.tibble(corr.test(autism_sign_mz, autism_sign_sex, , method = "pearson", adjust = "fdr")$p)
#' sex_p$mzid <- mzmatch_table$mzid
#' filter(sex_p, V1 <=0.05) # one sign
#' 
#' 
#' ## ggplot
#' hist(sex_r)
#' p <- ggplot(as.tibble(sex_r), aes(x = sex_r)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha=.2, fill = "#E69F00")
#' 
#' outdirectory <- "results"
#' outfilename <- sprintf("%s_hist_ewas_sign_corr_sex_r.png", saveName)
#' ggsave(filename = file.path(outdirectory,outfilename), plot = p, scale = 1, dpi = 400)






# ******** -----
# TODO -----------------------------------------------------------------




# ******** -----
# Supplement -----------------------------------------------------------



# .\ sub-tasks ----
# .\ sub-tasks \ tasks ----



# ******** -----
# Playground, Learning: Send to Quiver


