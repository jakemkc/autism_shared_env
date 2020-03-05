## Feb 21 2018
## Goal: plot EWAS results

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings
# quartz(height=6, width=8)
# dev.off() # reset par
# getwd()


# ******** ---
## Goal
# 1. create volcano plot
# 2. create qq plot
# 
# ******** ---


# ******** -----
# read-data ------------------------------------------------------------

load("results/c18_ewas.Rdata") 

## Switch for saving

## no scale data for plot
# saveName = "hilic"
# saveName = "c18"

## scaled data for plot
# ggrepelReadName = "hilic"
# saveName = "hilicScale"
# retFrame <- retFrameScale
ggrepelReadName = "c18"
saveName = "c18Scale"
retFrame <- retFrameScale

# ******** -----
# A. Prep plots  -------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggrepel)

# colnames(retFrame) <- c("Estimate", "SE", "z", "p", "mz", "OR")

# ******** -----
# B. exploratory plots  -------------------------------------------

# ## exploratory
# hist(retFrame$p.value,
#      main="mv_8_semen, p",
#      xlab="p",
#      border="blue",
#      col="green",
#      ylim = c(0, 2),
#      # xlim=c(100,700),
#      # las=1,  # labels are perpendicular to axis, 0 = parallel, 1 = perpendicular
#      # breaks=5,
#      prob = TRUE  #probability density expressed through the y-axis instead of the regular frequency.
# )


# ******** -----
# C. Manahatten plot ------------------------------------------------------
# 
# title <- "autism risk of mz"
# 
# # def threshold for horizontal line
# p1_e3 <- 1e-3       # 0.001
# p5_e3 <- 5e-3       # 0.005
# pbon <- 0.05/nrow(retFrame)   # 4.200269e-06
# 
# # def what to label later 
# pbonlab <- subset(retFrame, p < pbon)
# p1_e3lab <- subset(retFrame, p < 1e-3)
# p5_e3lab <- subset(retFrame, p < 5e-3)

# 
# 
# ### Drop useless var in plot
# ###' Try do it in plot.R after fdr, as the total variable affects fdr estimate. 
# ###' Unless I know they are bad and should be excluded from EWAS, it's just for "plot"
# #     # from volcano plot there are some with high coef but p value close to 1
# #     # Use the code in playground to discover the issue: control values are all "0"
# # 
# tmp_dropmz <- retFrame %>% filter(p >0.95, Estimate > 0.5) # 30
#  
# tmp_dropmz_index <- which(retFrame$mz %in% tmp_dropmz$mz)
#  
# retFrame <- retFrame[-tmp_dropmz_index, ]
# 
# 
# 
# # plot
# library(ggplot2)
# p <- ggplot(retFrame, aes(mz, -log10(p)))
# p <- p + geom_point(size = 0.5) 
# 
# # Threshold lines
# # Bonferroni correction, divide the critical P value (α) by the number of comparisons being made.
# p <- p + geom_hline(yintercept=-log10(pbon), col="red") 
# p <- p + geom_hline(yintercept=-log10(pfdr), col="orange")
# # p <- p + geom_hline(yintercept=-log10(p1_e3), col="green")
# # p <- p + geom_hline(yintercept=-log10(0.01), col="yellow")
# 
# # label axis
# p <- p + xlab("mz") + ggtitle(sprintf("%s. GLM age sex adjusted", title))
# 
# # label dot by threshod
# # p <- p + geom_text(data = pbonlab, aes(mz, -log10(p), label=as.character(mz)), size=3, vjust=0)
# 
#     # p <- p + geom_text(data=p1_e3lab, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), size=3, vjust=0)
#     # p <- p + geom_text(data=p5_e3lab, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), size=3, vjust= -0.8)
#     # vjust: (default: 0.5) position of the anchor (0=bottom edge, 1=top edge), 
#     # can go below 0 or above 1; similar for hjust
# 
# # label dot by fdr threshold (red)
# 
# # p <- p + geom_text(data=fdrsign, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), hjust=0, vjust = -0.8, color = "red")
# 
# p <- p + annotate("text", x = 40, y = 5.5, label = "bonf correction", size = 3)
# p <- p + annotate("text", x = 20, y = 3.22, label = "fdr", size = 3)
#     # p <- p + annotate("text", x = 6, y = 2.4, label = "p = 0.005", size = 3)
#     # p <- p + theme_bw() # use a b/w color theme
# p <- p + theme(plot.title = element_text(hjust = 0.5), axis.text.x  = element_text(angle=75, vjust=0.7, hjust= 0.7, size=2)) # move to the left a bit (using vjust since the labels are rotated)
# p
# 
# # quartz(height=6, width=8)
# 
# 
# 
# ## for ggplot only: you can specify or it will take default including to save the last ggplot
# savingpath <- "results/manha_hilic_glm.png"
# ggsave(savingpath, scale=1, dpi=400)
# # ggsave(savingpath, scale=1, dpi=200)



# ******** -----
# E. Volcano plot ---------------------------------------------

# 042919
# retFrame <- retFrame %>% select(-by)
# write_csv(retFrame, path = file.path("results", "c18EWASplot.csv"))

# FDR label
fdrsign <- retFrame %>% filter(fdr <= 0.05)  #row = 113
# fdrlab <- subset(retFrame, fdr < 0.05)

# FDR equ p value (to draw FDR line)
fdrsign %>% arrange(p.value) %>% tail
# fdrsign$p.value %>% order %>% fdrsign[., ] # corresponding p value of largest fdr: 4.966986e-04
# manual modifiy
pfdr <- 0.000394

## create color vector for volcano plot
retFrame$vgrp <- 5

tmpindex1 <- which(retFrame$fdr <= 0.05 & retFrame$fdr > 0.025 & retFrame$estimate >0)
tmpindex2 <- which(retFrame$fdr < 0.025 & retFrame$estimate >0)
tmpindex3 <- which(retFrame$fdr <= 0.05 & retFrame$fdr > 0.025 & retFrame$estimate <0)
tmpindex4 <- which(retFrame$fdr < 0.025 & retFrame$estimate <0)

retFrame$vgrp[tmpindex1] <- 1
retFrame$vgrp[tmpindex2] <- 2
retFrame$vgrp[tmpindex3] <- 3
retFrame$vgrp[tmpindex4] <- 4
retFrame$vgrp %>% as.factor %>% levels
retFrame$vgrp %>% table

## rename level using tidyverse function
retFrame$vgrp <- retFrame$vgrp %>% as.factor
retFrame$vgrp <- fct_recode(retFrame$vgrp, "FDR<0.05, +ve" = "1", "FDR<0.025, +ve" = "2", "FDR<0.05, -ve" = "3", "FDR<0.025, -ve" = "4", "Not Sign" = "5")

## rename factor for volcano plot
colnames(retFrame)[which(colnames(retFrame) == "vgrp")] <- "FDR groups"

## color of the group in volcano plot
# cols <- c("1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange", "5" = "black")
cols <- c("FDR<0.05, +ve" = "#CE50CA", "FDR<0.025, +ve" = "red", "FDR<0.05, -ve" = "#74D944", "FDR<0.025, -ve" = "blue", "Not Sign" = "black")


# FDR range
quantile(retFrame$fdr)
hist(retFrame$fdr)


## plot
title <- "autism risk of mz"

p <- ggplot(retFrame, aes(estimate, -log10(p.value)))
# p <- p + geom_point(size = 0.7) + xlim(-20, 20)
p <- p + geom_point(aes(colour = `FDR groups`), size = 0.7) + xlim(-15, 15) # mind it's a plot trick, not full data
p <- p + scale_colour_manual(values = cols)
# p <- p + geom_point() + xlim(-.8, .8) # special setting for logi_catstcri


# Threshold lines
# Bonferroni correction, divide the critical P value (α) by the number of comparisons being made.
# p <- p + geom_hline(yintercept=-log10(pbon), col="red") 
p <- p + geom_hline(yintercept=-log10(pfdr), col="orange")
# p <- p + geom_hline(yintercept=-log10(p1_e3), col="green")
# p <- p + geom_hline(yintercept=-log10(0.01), col="yellow")

# label axis
p <- p + xlab("Log odds") + ggtitle(sprintf("%s. GLM age sex adjusted", title))
p <- p + ylab(expression(paste(-log[10], (italic(p))))) #?plotmath, from LIFE plots_mv_EWAS code

# annotate those p < some value with ggrel

# ******** -----
#' #' off this section if you don't have annotation.R results yet
tablemz_sign_unique <- readRDS(file.path("results", sprintf("tablemz_sign_unique_%s.rds", ggrepelReadName))) #from annotatin.R
retFrame_1 <- retFrame[which(retFrame$mzid %in% tablemz_sign_unique$mzid), ]
retFrame_1 <- retFrame_1 %>% filter(p.value <=3.951633e-05) # top 4
retFrame_1 <- retFrame_1 %>% left_join(tablemz_sign_unique[, c("mzid", "Name")], by = "mzid")
# tmp_index1 <- which(tablemz_sign_unique$mzid %in% retFrame_1$mzid)
# retFrame_1 <- bind_cols(retFrame_1, as.data.frame(tablemz_sign_unique[tmp_index1, "Name"]))

p <- p + geom_text_repel(data = retFrame_1, aes(label = Name), nudge_x = 0.5,
                         size = 3,
                         box.padding = 1,
                         point.padding = 0.4)


# label dot by threshod
# p <- p + geom_text(data = pbonlab, aes(Estimate, -log10(p), label=as.character(mz)), size=3, vjust=0)

    # p <- p + geom_text(data=p1_e3lab, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), size=3, vjust=0)
    # p <- p + geom_text(data=p5_e3lab, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), size=3, vjust= -0.8)
    # vjust: (default: 0.5) position of the anchor (0=bottom edge, 1=top edge), 
    # can go below 0 or above 1; similar for hjust

# label dot by fdr threshold (red)
# p <- p + geom_text(data=fdrsign, aes(indvar, -log10(two_tail_p), label=as.character(indvar)), hjust=0, vjust = -0.8, color = "red")

# p <- p + annotate("text", x = -18, y = 5.5, label = "bonf correction", size = 3)
p <- p + annotate("text", x = -19, y = 3.22, label = "FDR = 0.05", size = 3)
# p <- p + theme_bw() # use a b/w color theme
p

# quartz(height=6, width=8)

outdirectory <- "results"
outfilename <- sprintf("volca_glm_%s.pdf", saveName)
savingpath <- file.path(outdirectory,outfilename)
ggsave(savingpath, scale = 1, dpi = 400)

outfilename <- sprintf("volca_glm_%s.svg", saveName)
ggsave(file.path(outdirectory,outfilename), scale=1, dpi=400)


# ******** -----
# D. QQ plot --------------------------------------------------------------------------

# the code is using chi square sth, as what I found online
# ii <- order(-log10(retData7$pvalue))
# this has no use, the formula included the sort
# this is like Excel select table and sort by 1 column, other also changed
# this one return the index, so you can use df [, ii] to rearrange

# Using plot function
#plot(-log10((nrow(retData7):1)/nrow(retData7)), sort(-log10(retData7$pvalue)), xlab='expected -log10(pvalue) under null', ylab='actual -log10(pvalue)', main='', ylim=c(0,5))
#abline(0,1)

# use ggplot
actualp <- -log10((nrow(retFrame):1)/nrow(retFrame))
expectedp <- sort(-log10(retFrame$p.value))
qqplot <- cbind(actualp,expectedp)
qqplot <- data.frame(qqplot)
str(qqplot)

p <- ggplot(qqplot, aes(actualp, expectedp))
p <- p + geom_point()
p <- p + xlab("expected -log10(pvalue) under null") + ylab("actual -log10(pvalue)") + ggtitle(sprintf("%s. QQ plot. GLM", title))
p <- p + scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,7))
p <- p + geom_abline(intercept = 0, slope = 1, color = "red")
p

# quartz(height=6, width=8)

outdirectory <- "results"
outfilename <- sprintf("qqplot_glm_%s.pdf", saveName)
savingpath <- file.path(outdirectory,outfilename)
ggsave(savingpath, scale=1, dpi=400)



# ******** -----
## S. save-rdata ----

outdirectory <- "results"
# outfilename <- "hilic_plots.Rdata"
outfilename <- sprintf("%s_plots.Rdata", saveName)
save(file=file.path(outdirectory,outfilename), 
     data_original, data_pheno,
     autism,
     autism_ccEwas, retFrame,
     fdrsign
)


# Load data
# load(file.path(outdirectory,outfilename)) 

outdirectory <- "results"
outfilename <- sprintf("fdrsign_%s.csv", saveName)
write.csv(fdrsign, file = file.path(outdirectory,outfilename))




# ******** -----
# E. Special analysis on OR --------------------------------------------------------------------------
# 
# 
# ## no of sign mz
# dim(fdrsign) # 119 significant mz
# 
# ## OR and fdr
# hist(fdrsign$fdr)
# hist(fdrsign$OR) # one guy seems get a huge OR
# summary(fdrsign$OR)
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # 0.000     0.083     0.242    59.214     0.371 11262.768 
# 
# ## which guy
# which(fdrsign$OR >500) %>% fdrsign[.,]
# # Estimate       SE        z            p    mz       OR         fdr
# # I(scale(log10(v2364 + 1))) 9.329258 2.178871 4.281694 1.854755e-05 v2364 11262.77 0.005825397
# 
# ## check raw data (or run the glm)
# table(autism_cc$v2364, autism_cc$Role)
# 
# ## It's mz v2364 sorted. all controls in the lower range
# #           Control Proband
# # 0                3       0
# # 1130339.2        1       0
# # 1221940.2        1       0
# # 1399056.4        0       1
# # 1503086.6        1       0
# # 1537649.8        1       0
# # 1572249.6        1       0
# # 1761526.8        1       0
# # 2014695.3        1       0
# # 2069901          1       0
# # 2081031.2        1       0
# # 2241261.6        1       0
# # 2310939.6        1       0
# # 2433647.4        0       1
# # 2438876.5        1       0
# # 2528779          1       0
# # 2581418.5        1       0
# # 2592753.7        0       1
# # 2721619.9        1       0
# # 2832264.3        1       0
# # 2999083          0       1
# # 3120139.3        0       1
# # 3192687.4        0       1
# # 3278231.2        1       0
# # 3361557.6        0       1
# # 3375415.5        0       1
# # 3404593.9        1       0
# # 3599358.7        1       0
# # 3624305.4        0       1
# # 3794532          0       1
# # 3935893.2        1       0
# # 3999511.4        0       1
# # 4504336.5        1       0
# # 4533482.6        1       0
# # 4895482.5        1       0
# # 5003313.9        0       1
# # 5351645.8        0       1
# # 5442153          1       0
# # 5893925.9        0       1
# # 6578684.5        1       0
# # 6750587.7        0       1
# # 6875145          0       1
# # 6960248.5        0       1
# # 7092352.2        1       0
# # 7191947.9        0       1
# # 7318306          0       1
# # 7465010.1        0       1
# # 8917918.5        0       1
# # 8969402.3        0       1
# # 9513209.5        0       1
# # 9530319.1        0       1
# # 9668939.3        0       1
# # 9731388.9        0       1
# # 10002855.8       0       1
# # 10281792.3       0       1
# # 10331178.8       0       1
# # 10349930.5       0       1
# # 10433157.4       0       1
# # 10596589.6       0       1
# # 10661356.2       0       1
# # 10812745.6       0       1
# # 10894717.5       0       1
# # 10966616.4       0       1
# # 11409658.4       0       1
# # 11446378.6       0       1
# # 11658098.3       0       1
# # 11852743.1       0       1
# # 11890394         0       1
# # 12003080.2       0       1
# # 12251533.6       0       1
# # 12325130.4       0       1
# # 12377250.8       0       1
# # 12389668.9       0       1
# # 12806918.9       0       1
# # 12812123.2       0       1
# # 13078045.1       0       1
# # 13345794.9       0       1
# # 13422102.4       0       1
# # 13509499.8       0       1
# # 13616570.1       0       1
# # 13686650.2       0       1
# # 13748142.4       0       1
# # 13835439.2       0       1
# # 13948587.1       0       1
# # 13984294.6       0       1
# # 13990947.7       0       1
# # 14377321.3       0       1
# # 14737235.4       0       1
# # 15100618         0       1
# # 15151426.6       0       1
# # 15721364.1       0       1
# # 16135941.2       0       1
# # 16139459.5       0       1
# # 16253337.8       0       1
# # 16269912         0       1
# # 16696469.7       0       1
# # 16718184.9       0       1
# # 16851945.3       0       1
# # 17445773.5       0       1
# # 18496538.4       0       1
# # 18529723.8       0       1
# # 18761099.5       0       1
# # 18819951         0       1
# # 19078217.9       0       1
# # 19278537.9       0       1
# # 20358235.7       0       1
# # 21022239.8       0       1
# # 24994544.8       0       1
# # 27696691.7       0       1
# 
# ## as a plot
# p <- ggplot(retFrame, aes(OR, -log10(p)))
# p <- p + geom_point() + xlim(-50, 100)
# 
# # ******** -----
# # TODO -----------------------------------------------------------------
# 
# 
# 
# 
# # ******** -----
# # Supplement -----------------------------------------------------------
# 
# 
# 
# # .\ sub-tasks ----
# # .\ sub-tasks \ tasks ----
# 
# 
# 
# # ******** -----
# # Playground, Learning: Send to Quiver
# 
# 
# 
# p <- ggplot(mtcars, aes(mpg, wt)) +
#     geom_point(aes(colour = factor(cyl)))
# p + scale_colour_manual(values = c("red", "blue", "green"))
# 
# 
# summary(fdrsign$p)
# c(-log10(min(fdrsign$p)), -log10(max(fdrsign$p))) %>% mean # 4.761408
# 10^-4.761408 # p value cut: 1.732176e-05
# 
# 
# 
# 
# fdrsign$vgrp <- 5
# 
# tmpindex1 <- which(fdrsign$fdr > 0.005 & fdrsign$Estimate >0)
# tmpindex2 <- which(fdrsign$fdr <= 0.005 & fdrsign$Estimate >0)
# tmpindex3 <- which(fdrsign$fdr > 0.005 & fdrsign$Estimate <0)
# tmpindex4 <- which(fdrsign$fdr <= 0.005 & fdrsign$Estimate <0)
# 
# fdrsign$vgrp[tmpindex1] <- 1
# fdrsign$vgrp[tmpindex2] <- 2
# fdrsign$vgrp[tmpindex3] <- 3
# fdrsign$vgrp[tmpindex4] <- 4
# fdrsign$vgrp %>% table

