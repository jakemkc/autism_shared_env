## Feb 20 2018
## Goal: autism EWAS

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
## Goal
# 1. assign variable class
# 2. filter case control
# 3. 
# ******** ---


# ******** -----
# A. read-data ------------------------------------------------------------

# Load data
load("results/c18.Rdata") 
# load("results/_= c18/_ 20 percent/c18.Rdata")
    # removed mz with > 80% "0"
    # removed duplicate mz
    # removed top 10% high CV mz based on chearqc
    # rows: case, control, and parents
    # columns: mz, pheno meta

## Switch for saving
# saveName = "hilic"
saveName = "c18"


# ******** -----
# B.  transform-explore ------------------------------------------------------------

library(tidyverse)

## No need to do row-wise normalization
# (apply(autism[, 1:8742], 1, sum) %>% sd) / (apply(autism[, 1:8742], 1, sum) %>% sum) * 100
# CV: 0.08%

## check class()
# autism %>% sapply(class) %>% table
# capture.output(str(autism, list.len = ncol(autism))) %>% listviewer::jsonedit(mode = "view")

## what classes do you have, and what are variables in that class
tmp1 <- split(names(autism),sapply(autism, function(x) paste(class(x), collapse=" ")))
# names(tmp1)
    # [1] "character" "integer"   "numeric"
# tmp1$character
    # [1] "emoryVT"  "SampleID" "BioID"    "ID"       "Sex"      "Role"     "Orig"    
    # [8] "Diag"     "c18name" 



## assign class()
# ******** ---
## pheno class:
# fct: sex, role, batch
#     - role as factor after extracting case-control
# num: age
# chr: let the rest be string first
#     - orig, diag could be factors
# ******** ---

fctnames <- c("Sex", "Batch") # FACTOR() Role done in later step
numnames <- c("Age")

# lapply(autism[,fctnames], as.factor)
autism[, fctnames] <- lapply(autism[, fctnames], as.factor)

autism[, numnames] <- autism[, numnames] %>% as.numeric


## check class()
# tmp1 <- split(names(autism),sapply(autism, function(x) paste(class(x), collapse=" ")))
# names(tmp1)
# [1] "character" "factor"    "numeric" 

# sapply(autism, is.factor) %>% autism[ ,.] %>% sapply(., levels)
    # $Batch
    # [1] "1" "2" "3" "4" "5"
    # 
    # $Sex
    # [1] "Female" "Male"  


# ******** -----
# C.  EWAS ------------------------------------------------------------

# "Role" colunm has case and control
    # xs code: use proband and control in role:  ae.mean.cc<-ae.mean.clean[,which(phenoClean$Role =="Proband" | phenoClean$Role =="Control")] # use role proband, not diag?
    # xs code: using log10, then log10. or I use log2: beta is 2 fold increase in X rather than 10 fold increase in X


## 1) Extract case and control
autism %>% filter(. , Role == "Proband") %>% dim # count is dplyr() and slow (because without the groupby?)
    # 82 cases, with dup

# unique case ("ID" is unique from BCH)
tmp2 <- autism %>% filter(. , Role == "Proband") %>% .$ID %>% table %>%  as.data.frame()
tmp2$Freq %>% length # 75 unique case by ID
sort(tmp2[,2]) # some duplicated 2 or 3 times

autism %>% filter(. , Role == "Control") %>% dim
    # 29 controls, no dup

# unique control
tmp3 <- autism %>% filter(. , Role == "Control") %>% .$ID %>% table %>%  as.data.frame()
tmp3$Freq %>% length # 29 unique case by ID
sort(tmp3[,2])

# Extract both case and control
autism_cc <- autism %>% filter(. , grepl("Proband|Control", Role))  #29 control, 82 proband
    # autism %>% filter(. , Role == "Proband" | Role == "Control") %>% dim

# 1.1) remove duplicate rows (by ID; pick first one)
autism_cc <- autism_cc[!duplicated(autism_cc$ID), ] #29 controls + 75 cases = 104 total



## 2) Factorize outcome 
# this can be done in the previous steps, even the factor level has case, control, father, mother; data only case and control
autism_cc$Role <- as.factor(autism_cc$Role)

# autism_cc$Role %>% str
# Factor w/ 2 levels "Control","Proband": 1 1 2 2 2 2 2 2 2 1 ...

# which is 0 
autism_cc$Role %>% contrasts
#           Proband
# Control       0
# Proband       1
#' it means it will create 1 dummy column named "proband" and it's equal to proband as it has a "1"


## 3) split case control
autism_case <- autism_cc %>% filter(Role == "Proband")
autism_control <- autism_cc %>% filter(Role == "Control")

metaNameVec <- autism_cc %>% select(-starts_with("v")) %>% colnames

# ## 4.1) removed mz with > 80% "0" per group
# #' Not doing this. explore.R done this already and is optional
# 
# percentzero <- function (x) {
#     100*mean(!(x))
# }
# 
# # case
# tmp4 <- autism_case %>% select(starts_with("v")) %>% sapply(., percentzero)
# tmp41 <- names(tmp4)[which(tmp4 <= 80)] # find names < 80% "0", i.e., keep mzid names
# autism_case <- autism_case %>% select(tmp41, metaNameVec)
# 
# # control
# tmp4 <- autism_control %>% select(starts_with("v")) %>% sapply(., percentzero)
# tmp41 <- names(tmp4)[which(tmp4 <= 80)] # find names < 80% "0", i.e., keep mzid names
# autism_control <- autism_control %>% select(tmp41, metaNameVec)


## 4.2) remove mz with no variance, per group (jake: this is not correlation but EWAS, if larger N and have time, dont do it per group to keep the rare exposure in case or control)# F = no variant column

# case
tmp5 <- autism_case %>% select(starts_with("v")) %>% vapply(., function(x) length(unique(x)) > 1, logical(1L))
tmp51 <- names(which(tmp5 == TRUE))
autism_case <- autism_case %>% select(tmp51, metaNameVec)

# control
tmp5 <- autism_control %>% select(starts_with("v")) %>% vapply(., function(x) length(unique(x)) > 1, logical(1L))
tmp51 <- names(which(tmp5 == TRUE))
autism_control <- autism_control %>% select(tmp51, metaNameVec)


## 5) bind the case and control back

autism_case <- autism_case[intersect(colnames(autism_case), colnames(autism_control))]
autism_control <- autism_control[intersect(colnames(autism_control), colnames(autism_case))]

autism_ccEwas <- rbind(autism_case, autism_control)



## glm function

autism_glm <- function(depvar = "Role", indvar, dat) {
    setform <- sprintf("%s ~ log2(%s+1) + Age + Sex + Batch", depvar, indvar) #nature ln = log
    with(data = dat, glm(as.formula(setform), family = binomial))
}

autism_glm_scale <- function(depvar = "Role", indvar, dat) {
    setform <- sprintf("%s ~ scale(log2(%s+1)) + Age + Sex + Batch", depvar, indvar) #nature ln = log
    with(data = dat, glm(as.formula(setform), family = binomial))
}

## special preparation
### error in glm, due to removing duplicate then autism_cc$v2961 is all "0"
# autism_cc <- subset(autism_cc, select=-c(v2961))

## glm loop
library(broom)
retFrame <- data.frame()

loop_mz <- autism_ccEwas %>% select(starts_with("v")) %>% colnames

for(i in 1:length(loop_mz)) {
    tmp1_var <- loop_mz[i]
    #mod <- autism_glm(indvar = tmp1_var, dat = autism_ccEwas) %>% tidy
    modRaw <- autism_glm(indvar = tmp1_var, dat = autism_ccEwas)
    tmp1 <- modRaw %>% tidy()
    tmp2 <- modRaw %>% confint_tidy() # it calls confint(), it's 2.5% LW and 97.5% HI (run modRaw %>% confint)
    mod <- cbind(tmp1, tmp2)
    mod$mzid <- tmp1_var
    retFrame <- rbind(retFrame, mod[2,])
    print(tmp1_var) # shown "after" the warnings
}


###' warning
###' Jake: It has lots of same type of warning as below:
###' [1] "v6457"
###' Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
###' [1] "v6458"
###' Run:  autism_cc$v6458 %>% table(., autism_cc$Role)
###' It's basically as what I search online, the perfect separation in the model, close to 0 or 1, likelihood is close to infinity and SE is huge


## Estimates from log odds to odds ratio (nothing to do with log2 X)
retFrame <- retFrame %>% mutate(OR = exp(estimate), ORCIlow = exp(conf.low), ORCIhigh = exp(conf.high)) %>% mutate_at(c("OR", "ORCIlow", "ORCIhigh"), round, 3)

#' it's a tibble bug
#' I can use "View()" to see this df, but the colnmae name is "or.estimate" and I can't access the object from terminal
# retFrame$OR <- 2^(retFrame[,2]) # cant read from console
# retFrame$OR <- 2^(retFrame[,2]) %>% round(digits = 3) # cant read from console
# 
# retFrame <- retFrame %>% as.data.frame() # this works
# retFrame$OR <- 2^(retFrame[,2])


## FDR
retFrame$fdr <- p.adjust(retFrame$p.value, method='fdr')  # min(fdr)
retFrame$by <- p.adjust(retFrame$p.value, method='BY')  # min(fdr)

# retFrame %>% filter(fdr < 0.05) %>% nrow # 113

# retFrame$fdr <- p.adjust(retFrame$"Pr(>|z|)", method='fdr')  # min(fdr)


### Done these for autism meeting poster; no need in Oct2018 version
### Drop useless var in plot
###' Try do it in plot.R after fdr, as the total variable affects fdr estimate. 
###' Unless I know they are bad and should be excluded from EWAS, it's just for "plot"
#     # from volcano plot there are some with high coef but p value close to 1
#     # Use the code in playground to discover the issue: control values are all "0"
# tmp_dropmz <- retFrame %>% filter(p.value >0.9, estimate > 0.5)
# 
# tmp_dropmz_index <- which(retFrame$variable %in% tmp_dropmz$variable)
# 
# retFrame <- retFrame[-tmp_dropmz_index, ]

## glm loop scaled 011519
library(broom)
retFrameScale <- data.frame()

loop_mz <- autism_ccEwas %>% select(starts_with("v")) %>% colnames

for(i in 1:length(loop_mz)) {
    tmp1_var <- loop_mz[i]
    #mod <- autism_glm_scale(indvar = tmp1_var, dat = autism_ccEwas) %>% tidy
    modRaw <- autism_glm_scale(indvar = tmp1_var, dat = autism_ccEwas)
    tmp1 <- modRaw %>% tidy()
    tmp2 <- modRaw %>% confint_tidy() # it calls confint(), it's 2.5% LW and 97.5% HI (run modRaw %>% confint)
    mod <- cbind(tmp1, tmp2)
    mod$mzid <- tmp1_var
    retFrameScale <- rbind(retFrameScale, mod[2,])
    print(tmp1_var) # shown "after" the warnings
}

## Estimates from log odds to odds ratio (nothing to do with log2 X)
retFrameScale <- retFrameScale %>% mutate(OR = exp(estimate), ORCIlow = exp(conf.low), ORCIhigh = exp(conf.high)) %>% mutate_at(c("OR", "ORCIlow", "ORCIhigh"), round, 3)

## FDR
retFrameScale$fdr <- p.adjust(retFrameScale$p.value, method='fdr')  # min(fdr)
retFrameScale$by <- p.adjust(retFrameScale$p.value, method='BY')  # min(fdr)



# ******** -----
# D.  extra \ fold-change ------------------------------------------------------------

## calculated fold change for poster
# need annotation results

# ### case
# autism_case <- autism_ccEwas %>% filter(. , grepl("Proband", Role))  #29 control, 75 proband
# 
# # All read file path needs to be changed manually
# tablemz_sign_unique <- readRDS("results/tablemz_sign_unique_c18.rds") # from annotation.R
# 
# tmp_index1 <- which(colnames(autism_case) %in% tablemz_sign_unique$mzid)
# 
# autism_fc_case <- autism_case[, tmp_index1]
# 
# autism_fc_case <- log10(autism_fc_case+1) # make it normal for reasonable mean ratio (fold change)
# 
# case_log_mean <- sapply(autism_fc_case, mean)
# 
# ### control
# autism_control <- autism_ccEwas %>% filter(. , grepl("Control", Role))  #29 control, 75 proband
# 
# tablemz_sign_unique <- readRDS("results/tablemz_sign_unique_c18.rds")
# 
# tmp_index1 <- which(colnames(autism_control) %in% tablemz_sign_unique$mzid)
# 
# autism_fc_control <- autism_control[, tmp_index1]
# 
# autism_fc_control <- log10(autism_fc_control+1)
# 
# control_log_mean <- sapply(autism_fc_control, mean)
# 
# ### fold change
# foldchange_table <- case_log_mean/control_log_mean
# 
# # save
# outdirectory <- "results"
# outfilename <- sprintf("foldchange_table_%s.csv", saveName)
# write.csv(foldchange_table, file = file.path(outdirectory,outfilename))




# ******** -----
## S. save-rdata ----

outdirectory <- "results"
outfilename <- sprintf("%s_ewas.Rdata", saveName)
save(file=file.path(outdirectory,outfilename),
     data_original, data_pheno,
     autism,
     autism_ccEwas, retFrame, retFrameScale
)

# Load data
# load(file.path(outdirectory,outfilename)) 







# ******** -----
# TODO -----------------------------------------------------------------




# ******** -----
# Supplement -----------------------------------------------------------



# .\ sub-tasks ----
# .\ sub-tasks \ tasks ----



# ******** -----
# Playground, Learning: Send to Quiver


## EWAS volcano plot shown some are high estimate but large P, why?

## 1. checked the QC CV, if there are found in QC, the CV seems OK

# ## 2. check case and control
# CV <- function (x) {
#     sd(x)/mean(x)*100
# }
# 
# autism_casecheck <- autism %>% filter(. , grepl("Proband", Role))  #29 control, 82 proband
# 
# tmp1 <- select(autism_casecheck, "v2664", "v3247", "v8715", "v8356")
# tmp1 <- log10(tmp1+1)
# tmp2 <- sapply(tmp1, CV)
# 
# # v2664     v3247     v8715     v8356 
# # 85.01448 133.02722 108.40533 132.61363 
# 
# autism_controlcheck <- autism %>% filter(. , grepl("Control", Role))  #29 control, 82 proband
# 
# tmp1 <- select(autism_controlcheck, "v2664", "v3247", "v8715", "v8356")
# tmp1 <- log10(tmp1+1)
# tmp2 <- sapply(tmp1, CV)
# 
# # v2664 v3247 v8715 v8356 
# # NaN   NaN   NaN   NaN 
#     # why? because they are all 0, no variance of these in control.  We have variance in 
# 

## parallel

# library(foreach) 
# library(doParallel)
# registerDoParallel(cores = detectCores() - 1)
# 
# 
# tmpx <- 
#     foreach(i = 1:length(loop_mz))  %dopar% {
#         tmpvar <- loop_mz[i]
#         mod <- autism_glm(indvar = tmp1_var, dat = autism_ccEwas)
#         frm <- as.data.frame(coef(summary(mod)))
#     }
# 
# 
# stopImplicitCluster()
# 
# 



