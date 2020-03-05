## Feb 14 2018
## Goal: Test autism data

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
## Goals
# 1. data_pheno: 1) missing at NA; 2)set as.integer
# 2. hilicdata: 1) create long format (cols: mz + emoryVT); 2) mzid as v1 v2..of mz; 3) emoryVT: emory queue ID
#    note: hilic data is already averaged, according to sekwon
# 
# role: data collector
# 3. select (clean) data columns: 1) remove mz > 80% missing; 2) remove no variance mz; 3) remove top 10% high CV
# 4. check QCs: 1) PCA; 2) CV
# 5. do the same cleaning on autism data. ready for next analysis
# ******** ---


# ******** -----
# A. read-data ------------------------------------------------------------

data_original <- read.csv(file = "data/hilic_ComBat_mzcalibrated_untargeted_averaged_featuretable.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
data_pheno <- read.csv(file = "data/hilic_pheno.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


## Switch for saving
saveName = "hilic"
# saveName = "c18"


# ******** -----
# B.  clean-up ------------------------------------------------------------

library(tidyverse)
# library(DT)


## check werid values
# head(data_original)
    # ✔ mz: continous with 0
    # ✔ excel browsed. no NA or werid, just number and 0.
    # datatable(data_original)
        # ✔ manual sampled. sorted without weird symbol
    # str(data_original)
    # colnames(data_original) # one of them has a werid VT for c18 data
    
    
    # head(data_pheno, n = 20)
        # ✖︎ empty cells, NA, <NA>
    # datatable(data_pheno)
        # ✔ empty
    # table(data_pheno$Diag)
        # ✔ NA
# colnames(data_pheno)


## replace "" (blk empty cell) & " " (space) as NA
# str(data_pheno)
    # 'data.frame':	777 obs. of  11 variables:
    # $ FileName: chr  "VT_170719_M213_001" "VT_170719_M213_003" "VT_170719_M213_005" "VT_170719_M213_007" ...
    # $ SampleID: chr  "nist_1" "nist_3" "nist_5" "q3June2014_1a_1" ...
    # $ Batch   : int  1 1 1 1 1 1 1 1 1 1 ...
    # $ BioID   : chr  "nist_1" "nist_3" "nist_5" "q3June20" ...
    # $ ID      : chr  "" "" "" "" ...
    # $ Sex     : chr  "" "" "" "" ...
    # $ Age     : int  NA NA NA NA NA NA NA NA NA NA ...
    # $ Role    : chr  "" "" "" "" ...
    # $ Orig    : chr  "" "" "" "" ...
    # $ Diag    : chr  "" "" "" "" ...
    # $ c18name : chr  "VT_170719_M213_002" "" "" "VT_170719_M213_008" ...

# sub blank etc as NA
data_pheno <- as.data.frame(apply(data_pheno, 2, function(x) gsub("^$|^ $", NA, x)), stringsAsFactors = FALSE) 
    # ^ should be front, $ should be end, | should be OR

# class
data_pheno$Batch <- as.integer(data_pheno$Batch) # 1-5
data_pheno$Age <- as.integer(data_pheno$Age) # in month
# str(data_pheno)


# ******** -----
# C. processing ------------------------------------------------------------

## Transpose data_original

### add mzid column
hilic1 <- data_original
hilic1$mzid <- paste("v", c(1:nrow(hilic1)), sep="")

# # 061919
# tmp1 <- hilic1 %>% select(mz, time, mzid)
# tmp2 <- c("v1051", "v10677", "v1391", "v1663", "v1827", "v24", "v2926", "v3231", "v3455", "v4564", "v4608", "v4887", "v5936", "v6000", "v6326", "v6519", "v6650", "v7108", "v8559", "v8589", "v8983")
# tmp3 <- tmp1 %>% filter(mzid %in% tmp2)


### name rows, columns, extract data
    # colnames(hilic1)
rownames(hilic1) <- hilic1$mzid # t() on colnames & rownames
    # head(hilic1[1:10, 1:10]) # non sample id such as time, mz.min
hilic1 <- hilic1 %>% select(starts_with("VT_")) # only take sample columns (VT_xxx), forget "mz", "time", "mzid" ...


## transpose
hilic2 <- as.data.frame(t(hilic1)) # samples = row
    class(hilic2)
#    colnames(hilic2)[1:10]
#    hilic2[1:10, 1:10]
hilic2$emoryVT <- rownames(hilic2) # ensure matching order from t()
rownames(hilic2) <- NULL
    #colnames(hilic2)[13000:13099]
    # head(hilic2)[,1:10]

# check data
# hilic2 %>% tibble::as_tibble(.) %>% dplyr::select(1:6)
# hilic2 %>% tibble::as_tibble(.) %>% dplyr::select((ncol(.)-6):ncol(.))
    
## skewness check
# par(cex.axis = 1.2, cex = 0.4, las = 2, fg="black", col="black", mar = c(15, 4, 4, 2) + 0.1)
# boxplot(scale(hilic2[, sample(1:13099, 200)]), pars = list(boxwex = 0.6, outcol="black"))
# dev.copy(png, filename="results/boxplot_scale_random.png", width=1024, height=768)
# dev.off()

# boxplot(scale(log((hilic2[, sample(1:13099, 200)]))+1), pars = list(boxwex = 0.6, outcol="black"))
# dev.copy(png, filename="results/boxplot_scale_log_random.png", width=1024, height=768)
# dev.off()
    # ✔ log and scale: OK distribution.

# ******** -----
# D. transform-explore ------------------------------------------------------------

## QC data

### join LC data & BCH pheno df
str(data_pheno)
str(hilic2$emoryVT)
join1 <- left_join(hilic2, data_pheno, by = c("emoryVT" = "FileName")) # VT_170719_M213_005 etc

# check what's at the end (last 15 var)
head(join1)[,(ncol(join1)-15):ncol(join1)]
str(join1[,(ncol(join1)-15):ncol(join1)])

### extract QC: chearqc to a df
table(join1$SampleID)
tmp1 <- grep("chearplasma", join1$SampleID, value = FALSE) # value = F give me index match
chearqc <- join1[tmp1, ]
# str(chearqc)
# chearqc$Batch

### extract QC: q3qc to a df
tmp2 <- grep("q3June2014", join1$SampleID, value = FALSE) # value = F give me index match
q3qc <- join1[tmp2, ]
# str(q3qc)
# q3qc$Batch

### get autism df (remove QC only samples)
# only QC samples in emory queue has no "ID" in data_pheno (BCH ID)
is.na(join1$ID) %>% which %>% length # 62 NA; 
autism_n <- join1[!is.na(join1$ID), ] # 197 obs (case, control, parents); some maybe repeated dummy by BCD

# tmp3 <- autism_n[, (ncol(autism_n)-10):ncol(autism_n)]


## name vector
# create non mzid (vxxx) in autism_n (length : 11)
nonMZnames <- autism_n %>% select(-matches("^v[0-9].*$")) %>% names

# [1] "emoryVT"      "Sample.ID" "Batch"     "BioID"     "ID"        "Sex"       "Age"       "Role"     
# [9] "Orig"      "Diag"      "c18.name" 


# ******** -----
# F1. chearqc ------------------------------------------------------------

### This section select mz columns for cleaning
#' 1. cut into mz only (for cleaning columns) and meta columns
#' 2. clean x3 on mz columns
#' 3. merge final mz columns to meta columns 


# create mz and meta datasets
clean1 <- chearqc %>% select(starts_with("v")) # select mz only, remove 11 meta var

chearqcMeta <- chearqc %>% select(-starts_with("v")) # select 11 meta only



### 1. remove columns > 80% "0" (No ND in the data but as "0"; i.e., keep mz >20% detected)

# sapply(clean1, FUN = function(x) typeof(x)) %>%  table # make sure all variable's typeof is compatible to function
    # table(sapply(clean1, FUN = function(x) typeof(x)))  # same


# get the "0" % vector
# function: percentage of "0"
percentzero <- function (x) {
    100*mean(!(x))
}

zero1 <- sapply(clean1, percentzero)   # sum(chearqc$v1 %in% 0) / nrow(chearqc) * 100
zero2 <- names(zero1)[which(zero1 < 80)] # find names < 80% "0", i.e., keep mzid names
clean2 <- clean1[, zero2] # df with needed variables #[] select columns from the shorter-length named vector



### 2. Remove mz columns with no variance
# vapply each column, find a unique value list (char or number) and count how many (into numeric), then create a T/F vector from a "test =1"
# F = no variant column


noVar1 <- vapply(clean2, function(x) length(unique(x)) > 1, logical(1L)) 
#1L just to make it less memory taking tha 1
# table(noVar1)
clean3 <- clean2[, noVar1] # df with needed variables; [] select columns from the matching-length logic vector



### 3. remove mz columns top 10% high CV (00818 new)

#' **** reminder
#' Intermediate vector to get CV is LOGGED to have normal
#' If good data, QC is repeated and should be normal; not for this omics data
#' FINAL data from this step is NOT logged as I selected by named variables
#' ****

CV <- function (x) {
    sd(x)/mean(x)*100
}

# tmpx1 <- chearqc3 %>% select(matches("^v[0-9].*$")) 

# sapply(tmpx1, FUN = function(x) typeof(x)) %>%  table # make sure all variable's typeof is compatible to function
clean4 <- log(clean3+1) # normal dist
highCV1 <- map(clean4, CV) # get CV
highCV1 <- highCV1 %>% as_tibble
highCV1 <- as_tibble(cbind(mzid = names(highCV1), t(highCV1))) # transpose and create tibble with mzid and CV cols
names(highCV1) <- c("mzid", "V1")
highCV1$V1 <- highCV1$V1 %>% as.numeric # for order
highCV2 <- highCV1 %>% arrange(V1) %>% slice(1:(round(0.80*nrow(highCV1)))) # keep lowest x x% CV
chearqchighCV2 <- highCV2

# remove mzid in chearQC
badmzChear <- highCV1 %>% arrange(desc(V1)) %>% slice(1:(round(0.20*nrow(highCV1)))) 
# sort from high to low CV, pick the highest 1-x% bad mz out

x <- highCV1 %>% arrange(desc(V1))

# Histogram of CV
# hist(highCV2$V1, breaks = 40, freq = FALSE)
# hist(highCV2$V1, breaks = 60)

# remove the high CV cols in qc data; 
# NOT in log scale
clean4 <- clean3 %>% select(highCV2$mzid)

summarytools::dfSummary(clean4[, 1:100]) %>% summarytools::view()

## chearqc: PCA

# log and scale
pca1 <- scale(log(clean4+1))

# add meta sample ID 
#' can only use row names with factorextra, I guess

rownames(pca1) <- chearqcMeta$SampleID # pca plot uses row names
rownames(pca1)[17] <- "chearplasma_3e_1" # damx they have typo. duplicated not allow in fviz
rownames(pca1)[18] <- "chearplasma_3f_1"

### pca
library(factoextra)
pca2 <- prcomp(pca1, retx = TRUE, center = FALSE, scale = FALSE)

fviz_pca_ind(pca2,
             col.ind = "cos2", # Color by the quality of representation; cos2 and contrib are auto controlled, see ?
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


## SO
#' 1. 2e1 and 2f1 different, 
#' but I can't use an interval to remove samples. 
#' because They were consecutive measurement, sandwiched by other QC that seems to be okay


outfilename <- sprintf("chearqc_score_plot_%s.png", saveName)
dev.copy(png, filename = file.path("./results", outfilename), width=1024, height=768)
dev.off()

# dev.copy(png, filename="results/chearqc_score_plot_low_res.png", width=800, height=600)
# dev.off()



## chearqc: CV dist

### use log2 for normal dist thus CV means sth; not scale() as it is mean zero, unit variance
qccv <- log(clean4+1)
qccv2 <- sapply(qccv, CV)


# qccv2 <- sapply(clean4, CV) # non log value

# # Histogram of CV
# hist(qccv2, breaks = 40, freq = FALSE)
# hist(qccv2, breaks = 60)

# dev.copy(png, filename="results/chearqc_cv_hist.png", width=1024, height=768)
# dev.off()



# cuulative dist (all 30 QCs)
qccv3 <- ecdf(qccv2)

plot(qccv3, xlab="percentage CV", ylab="mz fraction",main="")

outfilename <- sprintf("chearqc_cv_cdf_%s.png", saveName)
dev.copy(png, filename = file.path("./results", outfilename), width=1024, height=768)
dev.off()

library(ggplot2)
p <- ggplot(as.data.frame(qccv2), aes(x=qccv2)) + 
    geom_density(adjust = 1.5) # adj density curvey
p

outdirectory <- "./results"
outfilename <- sprintf("chearqc_cv_density_%s.png", saveName)
ggsave(file.path(outdirectory,outfilename), scale=1, dpi=400)


## cummulative CV dist (by batch, QCs)

nestCV <- function (x) {
    map(x, CV) %>% unlist
}


# get meta into data
qccv <- log(clean4+1)
qccv <- cbind(qccv, chearqcMeta["Batch"])

# too slow with summareize_all
# tmp1 <- qccv %>% group_by(Batch) %>% select(Batch, matches("^v[0-9].*$")) %>% summarize_all(CV)

# CV of all mz, by group
tmp1 <- qccv %>% group_by(Batch) %>% 
    nest %>%
    mutate(data2 = map(data, nestCV)) %>% select(-data) %>% unnest



# plot
library(ggplot2)

p <- ggplot(tmp1, aes(data2, colour = as.factor(Batch))) + stat_ecdf()
p <- p + ggtitle("Cumulative Density plot of CVs") +
    xlab("percentage CV") + 
    ylab("cumulative density")
p


outdirectory <- "./results"
outfilename <- sprintf("chearqc_cv_cdf_batch_%s.png", saveName)
ggsave(file.path(outdirectory,outfilename), scale=1, dpi=400)



# ******** -----
# F2. q3qc ------------------------------------------------------------

### This section select mz columns for cleaning
#' 1. cut into mz only (for cleaning columns) and meta columns
#' 2. clean x3 on mz columns
#' 3. merge final mz columns to meta columns 


# create mz and meta datasets
clean1 <- q3qc %>% select(starts_with("v"))

q3qcMeta <- q3qc %>% select(-starts_with("v"))



### 1. remove columns > 80% "0" (No ND in the data but as "0"; i.e., keep mz >20% detected)

# sapply(clean1, FUN = function(x) typeof(x)) %>%  table # make sure all variable's typeof is compatible to function
# table(sapply(clean1, FUN = function(x) typeof(x)))  # same


# get the "0" % vector
# function: percentage of "0"
percentzero <- function (x) {
    100*mean(!(x))
}

zero1 <- sapply(clean1, percentzero)   # sum(q3qc$v1 %in% 0) / nrow(q3qc) * 100
zero2 <- names(zero1)[which(zero1 < 80)] # find names < 80% "0", i.e., keep mzid names
clean2 <- clean1[, zero2] # df with needed variables



### 2. Remove mz columns with no variance
# vapply each column, find a unique value list (char or number) and count how many (into numeric), then create a T/F vector from a "test =1"
# F = no variant column


noVar1 <- vapply(clean2, function(x) length(unique(x)) > 1, logical(1L)) 
#1L just to make it less memory taking tha 1
# table(noVar1)
clean3 <- clean2[, noVar1] # df with needed variables



### 3. remove mz columns top 10% high CV (00818 new)

#' **** reminder
#' Intermediate vector to get CV is LOGGED to have normal
#' If good data, QC is repeated and should be normal; not for this omics data
#' FINAL data from this step is NOT logged as I selected by named variables
#' ****

CV <- function (x) {
    sd(x)/mean(x)*100
}

# tmpx1 <- q3qc3 %>% select(matches("^v[0-9].*$")) 

# sapply(tmpx1, FUN = function(x) typeof(x)) %>%  table # make sure all variable's typeof is compatible to function
clean4 <- log(clean3+1) # normal dist
highCV1 <- map(clean4, CV) # get CV
# highCV1 <- highCV1 %>% as.tibble
highCV1 <- highCV1 %>% as_tibble
highCV1 <- as_tibble(cbind(mzid = names(highCV1), t(highCV1))) # transpose and create tibble with mzid and CV cols
names(highCV1) <- c("mzid", "V1")
highCV1$V1 <- highCV1$V1 %>% as.numeric # for order
highCV2 <- highCV1 %>% arrange(V1) %>% slice(1:(round(0.80*nrow(highCV1)))) # remove largest 10% CV
q3qchighCV2 <- highCV2


# remove mzid in chearQC
badmzq3 <- highCV1 %>% arrange(desc(V1)) %>% slice(1:(round(0.20*nrow(highCV1)))) # keep low %, 


# Histogram of CV
# hist(highCV2$V1, breaks = 40, freq = FALSE)
# hist(highCV2$V1, breaks = 60)

# remove the high CV cols in qc data
clean4 <- clean3 %>% select(highCV2$mzid)


## q3qc: PCA

# log and scale
pca1 <- scale(log(clean4+1))

# add meta sample ID 
#' can only use row names with factorextra, I guess

rownames(pca1) <- q3qcMeta$SampleID # pca plot uses row names

### pca
library(factoextra)
pca2 <- prcomp(pca1, retx = TRUE, center = FALSE, scale = FALSE)

fviz_pca_ind(pca2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


## SO
#' 1. 2e1 and 2f1 different, 
#' but I can't use an interval to remove samples. 
#' because They were consecutive measurement, sandwiched by other QC that seems to be okay

outfilename <- sprintf("q3qc_score_plot_%s.png", saveName)
dev.copy(png, filename = file.path("./results", outfilename), width=1024, height=768)
dev.off()



## q3qc: CV dist

### use log2 for normal dist thus CV means sth; not scale() as it is mean zero, unit variance
qccv <- log(clean4+1)
qccv2 <- sapply(qccv, CV)

# # Histogram of CV
# hist(qccv2, breaks = 40, freq = FALSE)
# hist(qccv2, breaks = 60)

# dev.copy(png, filename="results/q3qc_cv_hist.png", width=1024, height=768)
# dev.off()

# cuulative dist
qccv3 <- ecdf(qccv2)

plot(qccv3, xlab="percentage CV", ylab="mz fraction",main="")

outfilename <- sprintf("q3qc_cv_cdf_%s.png", saveName)
dev.copy(png, filename = file.path("./results", outfilename), width=1024, height=768)
dev.off()



## cummulative CV dist (by batch, QCs)

nestCV <- function (x) {
    map(x, CV) %>% unlist
}


# get meta into data
qccv <- log(clean4+1)
qccv <- cbind(qccv, chearqcMeta["Batch"])

# too slow with summareize_all
# tmp1 <- qccv %>% group_by(Batch) %>% select(Batch, matches("^v[0-9].*$")) %>% summarize_all(CV)

# CV of all mz, by group
tmp1 <- qccv %>% group_by(Batch) %>% 
    nest %>%
    mutate(data2 = map(data, nestCV)) %>% select(-data) %>% unnest



# plot
library(ggplot2)

p <- ggplot(tmp1, aes(data2, colour = as.factor(Batch))) + stat_ecdf()
p <- p + ggtitle("Cumulative Density plot of CVs") +
    xlab("percentage CV") + 
    ylab("cumulative density")
p


outdirectory <- "./results"
outfilename <- sprintf("q3qc_cv_cdf_batch_%s.png", saveName)
ggsave(file.path(outdirectory,outfilename), scale=1, dpi=400)




# ******** -----
# F3. autism ------------------------------------------------------------
#' Final prep for autism data, role as measurer, ready for next analytical


### This section select mz columns for cleaning
#' 1. cut into mz only (for cleaning columns) and meta columns
#' 2. clean x3 on mz columns
#' 3. merge final mz columns to meta columns 


# create mz and meta datasets
clean1 <- autism_n %>% select(starts_with("v"))

autism_nMeta <- autism_n %>% select(-starts_with("v"))



### 1. remove columns > 80% "0" (No ND in the data but as "0"; i.e., keep mz >20% detected)

# sapply(clean1, FUN = function(x) typeof(x)) %>%  table # make sure all variable's typeof is compatible to function
# table(sapply(clean1, FUN = function(x) typeof(x)))  # same


# get the "0" % vector
# function: percentage of "0"
percentzero <- function (x) {
    100*mean(!(x))
}

zero1 <- sapply(clean1, percentzero)   # sum(autism_n$v1 %in% 0) / nrow(autism_n) * 100
zero2 <- names(zero1)[which(zero1 < 80)] # find names < 80% "0", i.e., keep mzid names
clean2 <- clean1[, zero2] # df with needed variables



### 2. Remove mz columns with no variance
# vapply each column, find a unique value list (char or number) and count how many (into numeric), then create a T/F vector from a "test =1"
# F = no variant column


noVar1 <- vapply(clean2, function(x) length(unique(x)) > 1, logical(1L)) 
#1L just to make it less memory taking tha 1
# table(noVar1)
clean3 <- clean2[, noVar1] # df with needed variables



### 3. remove mz columns top 10% high CV (based on chearQC)
# chearQC is plasma and close to real sample matrix

# remove the high CV cols (based on chearQC)

# autism <- clean3 %>% select(one_of(chearqchighCV2$mzid)) 
# could be a different of 9905 mz left vs 12059 mz left
autism <- clean3 %>% select(-one_of(badmzChear$mzid)) 
autism <- cbind(autism, autism_nMeta)


# ******** -----
## S. save-rdata ------------------------------------------------------------

outdirectory <- "results"
# outfilename <- "hilic_QC.Rdata"
outfilename <- sprintf("QC_%s.Rdata", saveName)
save(file=file.path(outdirectory,outfilename),
     data_original, data_pheno,
     chearqchighCV2, q3qchighCV2
)

# Load data
# load(file.path(outdirectory,outfilename)) 



outdirectory <- "results"
# outfilename <- "hilic.Rdata"
outfilename <- sprintf("%s.Rdata", saveName)
save(file=file.path(outdirectory,outfilename),
     data_original, data_pheno,
     autism
)


# Load data
# load(file.path(outdirectory,outfilename)) 


#
# haven::write_dta(autism, "./tmp/autism.dta")




# ******** -----
# TODO -----------------------------------------------------------------






# ******** -----
# Supplement -----------------------------------------------------------
#' 
#' # tidyverse version
#' 
#' c18_original <- read_csv(file = "tmp/c18_ComBat_mzcalibrated_untargeted_averaged_featuretable.csv", col_names = TRUE)
#' c18_original %>% as.tibble(.) %>% print(n=20)
#' c18_original %>% as.tibble(.) %>% glimpse(.)
#' 
#' 
#' glimpse(as.tibble(c18_original))
#' 
#' # view data
#' library(DT)
#' datatable(data_original)
#' 
#' library(listviewer)
#' jsonedit(c18_original)
#' jsonedit(data_original)
#' 
#' 
#' # .\ sub-tasks ----
#' # .\ sub-tasks \ tasks ----
#' 
#' 
#' 
#' 
#' 
#' 
#'