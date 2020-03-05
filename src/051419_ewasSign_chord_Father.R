## Feb 26 4 2018
## Goal: create exposome globe


## Initialize the workspace
rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
# dev.off()  # reset graphic par
# quartz(height=6, width=6)


# *****************************************************************************
## Notes
#' 1. EWAS sign = 113 mz; correlation matrix = 104 (filter bad m/z or not enough obs)
#' 2. mzmatch_table (row = 113): conversion table with ewas results + ID, class, color etc for chord; from corr.R
#'    Section A: create list (RT class x10)
#'    a) conversion table: row = 104 and create time class for parents
#'    b) create the chemical-time-class list for ASD from a) 
#'    c) create the chemical-time-class list for parent from a)
#'    
#'    Switch: 
#'    Section A: create list (cluster class x9)
#' 3. next time modify to have map() to run array of setting for best color/width setting
# *****************************************************************************





# ******** -----
# A. Create-lists  -------------------------------------------

library(tidyverse)

# load("results/hilic_corr_sign.Rdata") 
load("results/hilic_ewas_sign_corr.Rdata") 

## Switch for saving
saveName = "hilic"
# saveName = "c18"



# # *********
# ##' By RT class
# 
# ## 1 mzmatch table matching the new mz list in corr (i.e. lesser mz)
# tmp_index_match <- match(mzmatch_table$mzid_a, colnames(autism_case3))
# 
# mzmatch_table_2 <- mzmatch_table[!is.na(tmp_index_match), ]
# 
# mzmatch_table_2$timeclassp <- sprintf("%s.", mzmatch_table_2$timeclass) # timeclass label for parents
# 
# ## 2 chem-RT list for ASD
# 
# ### spread into short format (timeclass and mz partial dummy)
# mzmatch_long <- spread(mzmatch_table_2, timeclass, mzid_a) # on the whole df because spread will say error if 2 rows are the same. You can remove columns later
# mzmatch_long <- mzmatch_long %>% select("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
# 
# ### list
# tmp_mzmatch_list <- as.list(mzmatch_long)
# autism_case3_list <- lapply(tmp_mzmatch_list, function(x) x[!is.na(x)])
# 
# 
# ## 3 chem-RT list for parents
# 
# ### spread into short format (timeclass and mz partial dummy)
# mzmatch_long <- spread(mzmatch_table_2, timeclassp, mzid_p) # on the whole df because spread will say error if 2 rows are the same. You can remove columns later
# mzmatch_long <- mzmatch_long %>% select("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.", "J.")
# 
# ### list
# tmp_mzmatch_list <- as.list(mzmatch_long)
# tmp_parents_list <- lapply(tmp_mzmatch_list, function(x) x[!is.na(x)])
# 
# ## chem-RT parents-asd list
# 
# chord_list <- c(tmp_parents_list, autism_case3_list)


# *********
##' By cluster class 

## 1 create conversion dataframe

# read

outdirectory <- "./results"
# outfilename <- sprintf("chemOrder_fatherCaseIntra_%s.rds", saveName)
outfilename <- sprintf("chemOrder_motherCaseIntra_%s.rds", saveName) #using mother for consistent
chemOrder <- readRDS(file.path(outdirectory, outfilename))


# replace with correct mz name
chemOrder$chemIDp <- str_replace(chemOrder$chemIDp, ".x$", ".y") # .x is asd, .y is parent (mother or father)
chemOrder$mzid <- str_replace(chemOrder$mzid, ".x$", "") # .x is asd, .y is parent (mother or father)


chemOrder$group <- as.character(chemOrder$group) # numeric grp to char for asd

chemOrder$groupp <- paste(chemOrder$group, ".", sep = "") # char grp for parents

# create color for groups

colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                "#673770")

names(colorCodes) <- table(chemOrder$group) %>% names %>% dput()

which(names(colorCodes) %in% chemOrder$group)

chemOrder$grpcolor <- "n____"

for(i in 1:length(colorCodes)) {
    mtchC <- which(chemOrder$group %in% names(colorCodes)[i])
    chemOrder$grpcolor[mtchC] <- colorCodes[i]
}

## save
outdirectory <- "./results"
outfilename <- sprintf("chemOrderClusTab_Father_%s.rds", saveName)
saveRDS(chemOrder, file.path(outdirectory,outfilename))
# hello <- readRDS("results/chemOrder.rds")


## 2 chem-cluster list for ASD

# spread into short format (timeclass and mz partial dummy)
mzmatch_long <- spread(chemOrder, group, chemID) # on the whole df because spread will say error if 2 rows are the same. You can remove columns later
mzmatch_long <- mzmatch_long %>% select(table(chemOrder$group) %>% names %>% dput())

# list
tmp_mzmatch_list <- as.list(mzmatch_long)
autism_case3_list <- lapply(tmp_mzmatch_list, function(x) x[!is.na(x)])

# rearrange the cluster group (after initial plotting)
autism_case3_list <- autism_case3_list[c("9", "8", "5", "4", "3", "2", "6", "1", "7")]

## save
outdirectory <- "./results"
outfilename <- sprintf("chemOrderListASD_Father_%s.rds", saveName)
saveRDS(autism_case3_list, file.path(outdirectory,outfilename))
# hello <- readRDS("results/chemOrder.rds")


## 3 chem-cluster list for parents

# spread into short format (timeclass and mz partial dummy)
mzmatch_long <- spread(chemOrder, groupp, chemIDp) # on the whole df because spread will say error if 2 rows are the same. You can remove columns later
mzmatch_long <- mzmatch_long %>% select(table(chemOrder$groupp) %>% names %>% dput())

# list
tmp_mzmatch_list <- as.list(mzmatch_long)
tmp_parents_list <- lapply(tmp_mzmatch_list, function(x) x[!is.na(x)])

# rearrange the cluster group (after initial plotting)
tmp_parents_list <- tmp_parents_list[c("9.", "8.", "5.", "4.", "3.", "2.", "6.", "1.", "7.")]

## 4 chem-cluster parents-asd list
chord_list <- c(tmp_parents_list, autism_case3_list)



# ******** -----
# B. chord-prep  -------------------------------------------
# LIFE chord's code

# Make parents (topright) and asd (topleft) symmetric in the plot: 
# reverse asd chem cateogry order

ordertemp <- names(chord_list)
ordertemp[(length(autism_case3_list)+1):length(ordertemp)] <- rev(ordertemp[(length(autism_case3_list)+1):length(ordertemp)])
chord_list <- chord_list[ordertemp]



# revser asd chemical order within each cateogry

# For cl male category has _m suffix
# for(i in grep("m$", names(cl), ignore.case = T)) {
#     cl[[i]] <- rev(cl[[i]])
# }

for(i in 11:length(ordertemp)) {
    chord_list[[i]] <- rev(chord_list[[i]])
}



## create a chemical and chemical category directory vector. unlist(cl, use.names = F) can avoid extract funny name from cd
cd = structure(rep(names(chord_list), times = sapply(chord_list, length)),
               names = unlist(chord_list))

# structure(1:6, dim = 2:3); structure(1:9, dim = c(3,3))

# number of chem groups
n_group = length(chord_list)


# color x10 for RT
# colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
#                 "#673770", "#D3D93E") # "#38333E", "#508578", "#D7C1B1"

# color x9 for cluster
colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                "#673770")

# Chord category label (for making line break)
# label in 2 rows test with name "pcb \n mcd" etc
# linebrk <- c("PCBs", "OCPs", "Polybrominated_\ncpds", "PFASs", "Blood_metals", 
#   "Cotinine", "Phytoestrogens", "Phthalates", "Phenols", "Anti_microbial_cpds", 
#   "Paracetamols", "Urine_metals", "Urine_metalloids", "Urine_metalloids.m", 
#   "Urine_metals.m", "Paracetamols.m", "Anti_microbial_cpds.m", 
#   "Phenols.m", "Phthalates.m", "Phytoestrogens.m", "Cotinine.m", 
#   "Blood_metals.m", "PFASs.m", "Polybrominated_cpds.m", "OCPs.m", 
#   "PCBs.m")

# C. Chord diagram  -------------------------------------------

library(circlize)
circos.clear()
circos.par(
    gap.degree = c(rep(2, n_group/2-1), 10, rep(2, n_group/2-1), 10), 
    cell.padding = c(0, 0, 0, 0), 
    start.degree = 84,
    canvas.xlim = c(-1.5, 1.5),  # 1.2
    canvas.ylim = c(-1.5, 1.5)
)
# Why using without cell.padding, I see unequal cell width among sector?
# the concept: dawing mechanism for circos rectangle. You have the whole things defined, then to draw what, its has "border" concept
# You adj. boarder then it will auto draw the line to reveal rectangle.
# thus, if I don't set cell.padding to 0 (i.e. board gap of cells within sector), only the whole rectangle sector is in sub equal width
# and those with only 1 cell vs 30 cells, the 1 cell will be narrow
# see figure 7 in the vingette (1. circlize intro.pdf)
# after you set "cell.padding" to 0, doesn't matter you have "bg.border" or not as they will be overlapped
# this "limit", "padding" and "your object" x3 concept should apply to canvas area and ploting area in figure!
# You have have many invisible settings/ boundary/ padding(gap)/ concept in R graphic, implicityly needed for graphic
# Table 1 in vingette got the default parameters

circos.initialize(factors = names(chord_list), xlim = cbind(c(rep(0, n_group)), sapply(chord_list, length)))

circos.trackPlotRegion(
    ylim = c(0, 1), 
    bg.border = NA,
    bg.col = c(colorCodes, rev(colorCodes)),
    track.height = 0.04,
    panel.fun = function(x, y) {
        nm = get.cell.meta.data("sector.index")
        # ind = get.cell.meta.data("sector.numeric.index")  For label line break
        r = chord_list[[nm]]
        n = length(r)
        circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(1, n), lwd = 0.5)
        # circos.text(1:n - 0.5, rep(0.5, n), r, niceFacing = TRUE, cex = 0.7)  #cell text
        circos.text(n/2, 3.7, nm, adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = 0.65) 
        # circos.text(n/2, 3.5, linebrk[ind], adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = 0.6) # for label line break
    }
)




## color
# col_fun = colorRamp2(c(-0.75, 0, 0.75), c("darkgreen", "white", "red"))  # 0.75 range

# col_fun = colorRamp2(c(-1.4, -0.25, 0.05, 0.9), c("darkgreen", "white", "white", "red"))  
# col_fun = colorRamp2(c(-0.9, -0.45, 0.45, 0.9), c("darkgreen", "white", "white", "red"))  

# max range best double of what I want to show at with top 20% sth
col_fun = colorRamp2(c(-1.3, 0, 1.3), c("darkgreen", "white", "red"))  # 0.75 range
# col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))  # 0.75 range
# https://jsfiddle.net/teddyrised/g02s07n4/embedded/result/  site to read 8 digits hex code to rgb and color

# C. links for within family  -------------------------------------------
mat_fm = as.matrix(cor_father)

# hist(mat_fm)
###' start trial at rougly with 20% as cutoff

mat_fm[mat_fm < 0.5 & mat_fm > -0.5] <- NA  # set NA threshold



n = nrow(mat_fm)
rn = rownames(mat_fm)
cn = colnames(mat_fm)


for(i in 1:n) {
    for(j in 1:n) {
        g1 = cd[rn[i]]
        g2 = cd[cn[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rn[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == cn[j])
        if (is.na(mat_fm[i, j])) {
            # circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_mf[i, j]*4,
            #            col = "#00000000") # set transparency
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_fm[i, j]/1.5, h = 0.5, # h = 0.6,
                        col = col_fun(mat_fm[i, j]) )#, rou1 = 0.8, rou2 = 0.8) # rou1 = 0.5, rou2 = 0.5 are optional to shrink
            # lwd: use log to max the scale diff vs simple: mat_m[i, j]; log base: pick a value output 0-1 range, for the given input r range.
            # h's (abs(j-i)/(nm-1): vary height by distance; h's *value: how tight the line band together; h's + value: how far the band away track
        }
    }
}

# circos.clear()

# circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_fm[i, j]*2, h = 0.6, # h = 0.6,


# D. links for within parents  -------------------------------------------

mat_f = as.matrix(cor_fatherIntra)

mat_f[mat_f < 0.5 & mat_f > -0.5] <- NA  # set NA threshold

nf = nrow(mat_f)
rnf = rownames(mat_f)

for(i in 1:(nf-1)) {
    for(j in seq(i+1, nf)) {  # rev(seq(i+1, nf))
        g1 = cd[rnf[i]]
        g2 = cd[rnf[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rnf[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == rnf[j])
        if (is.na(mat_f[i, j])) {
            # circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*2, h = (abs(j-i)/(nf-1)*0.7),
            #            col = "#00000000") # set transparet
        } else if (g1 == g2) { # within class; k1-0.5 should be element cell in sector. width 1 by default in my layout, -0.5 = line from middle of the element cell
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]/1.5, h = (abs(j-i)/(nf-1)*-0.10), # h = (abs(j-i)/(nf-1)*0.68 + 0.02)
                        col = col_fun(mat_f[i, j]), rou1 = 1.01, rou2 = 1.01)
        } else if (mat_f[i, j] > 0) {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = log(abs(mat_f[i, j]), 2)+1, h = (abs(j-i)/(nf-1)*0.10 + 0.1),
                        col = col_fun(mat_f[i, j])) 
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = log(abs(mat_f[i, j]), 2)+1, h = (abs(j-i)/(nf-1)*0.03 + 0.26),
                        col = col_fun(mat_f[i, j])) 
            # lwd: use log to max the scale diff vs simple: mat_m[i, j]; log base: pick a value output 0-1 range, for the given input r range.
            # h's (abs(j-i)/(nm-1): vary height by distance; h's *value: how tight the line band together; h's + value: how far the band away track
        }
    }
    
}
# circos.clear()

# # E. links for within asd  -------------------------------------------
# 
# mat_m = as.matrix(cor_asd)
# 
# mat_m[mat_m < 0.6 & mat_m > -0.6] <- NA  # set NA threshold
# 
# nm = nrow(mat_m)
# rnm = rownames(mat_m)
# 
# for(i in 1:(nm-1)) {
#     for(j in seq(i+1, nm)) {
#         g1 = cd[rnm[i]]
#         g2 = cd[rnm[j]]
#         r1 = cd[cd == g1]
#         k1 = which(names(r1) == rnm[i])
#         r2 = cd[cd == g2]
#         k2 = which(names(r2) == rnm[j])
#         if (is.na(mat_m[i, j])) {
#             # circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.7),
#             #            col = "#00000000") # set transparet
#         } else if (g1 == g2) {
#             circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]/1.5, h = (abs(j-i)/(nm-1)*-0.10), # h = (abs(j-i)/(nm-1)*0.68 + 0.02)
#                         col = col_fun(mat_m[i, j]), rou1 = 1.01, rou2 = 1.01)
#         } else {
#             circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j], h = (abs(j-i)/(nm-1)*0.20 + 0.04),
#                         col = col_fun(mat_m[i, j]))
#             # h's (abs(j-i)/(nm-1): vary height by distance; h's *value: how tight the line band together; h's + value: how far the band away track
#         }
#     }
# }
# 
# circos.clear()


# E. links for within asd  -------------------------------------------

mat_m = as.matrix(cor_fathercaseIntra) 

mat_m[mat_m < 0.5 & mat_m > -0.5] <- NA  # set NA threshold

nm = nrow(mat_m)
rnm = rownames(mat_m)

for(i in 1:(nm-1)) {
    for(j in seq(i+1, nm)) {
        g1 = cd[rnm[i]]
        g2 = cd[rnm[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rnm[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == rnm[j])
        if (is.na(mat_m[i, j])) {
            # circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.7),
            #            col = "#00000000") # set transparet
        } else if (g1 == g2) {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]/1.5, h = (abs(j-i)/(nm-1)*-0.10), # h = (abs(j-i)/(nm-1)*0.68 + 0.02)
                        col = col_fun(mat_m[i, j]), rou1 = 1.01, rou2 = 1.01)
        } else if (mat_m[i, j] > 0) {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = log(abs(mat_m[i, j]), 2)+1, h = (abs(j-i)/(nm-1)*0.10 + 0.1),
                        col = col_fun(mat_m[i, j])) 
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = log(abs(mat_m[i, j]), 2)+1, h = (abs(j-i)/(nm-1)*0.03 + 0.26),
                        col = col_fun(mat_m[i, j]))  # lwd = # mat_m[i, j]
            # lwd: use log to max the scale diff vs simple: mat_m[i, j]; log base: pick a value output 0-1 range, for the given input r range.
            # h's (abs(j-i)/(nm-1): vary height by distance; h's *value: how tight the line band together; h's + value: how far the band away track
        }
    }
}

circos.clear()


# quartz.save("results/hilic_chord.png", type = "png", device = dev.cur(), dpi = 400)
outdirectory <- "./results"
outfilename <- sprintf("chord_father_%s.pdf", saveName)
quartz.save(file.path(outdirectory,outfilename), type = "pdf")


# lwd test:
# log(abs(0.65), 1.7)+1
# log(abs(0.68), 1.7)+1
# log(abs(0.8), 1.7)+1
# log(abs(0.9), 1.7)+1
# log(abs(1), 1.7)+1


# ******** -----
# x. TODO -----------------------------------------------------------
# link draw seq to have intra the uppest or as rou setting per group to show what I want as discussed in the paper
# For this goal on links:

# 1. done
# color
# thickness lwd
# rou
# h
# w

# 2. done 
# # the PCB to POP cluster, h seems got problem
# non linear h for male
# 
# > g1
# pcb177amt_M 
# "PCBsM" 
# > g2
# pophcbamt_M 
# "POPsM" 
# > i
# [1] 25
# > j
# [1] 48
# > (abs(j-i)/(nm-1)*0.7)
# [1] 0.1219697
# > 
#     
# so,
# it;s because the count of cell is in order, and female, it's left to right
# for male: categ order reversed, but not the within class cell order. it's PCBsM left to right, categ class from right to left reversed


# 3. 
# adj visual to you get the msg ASAP
# longer link first, then shorter on top



# done longer links draw first
# male
# links_m <- data.frame()
# for(i in 1:(nm-1)) {
#     for(j in seq(i+1, nm)) {
#         g1 = cd[rnm[i]]
#         # print(i)
#         g2 = cd[rnm[j]]
#         # print(j)
#         r1 = cd[cd == g1]
#         k1 = which(names(r1) == rnm[i])
#         r2 = cd[cd == g2]
#         k2 = which(names(r2) == rnm[j])
#         links_m_tmp = data.frame(i, g1, k1, j, g2, k2) # vector can only hold 1 class of information
#         links_m <- rbind(links_m, links_m_tmp)
#     }
# }
# 
# links_m$dist <- abs(links_m$j-links_m$i)
# links_m <- links_m[order(links_m$dist, decreasing = T),]
# links_m$g1 <- as.character(links_m$g1)
# links_m$g2 <- as.character(links_m$g2)
# 
# for(i in 1:nrow(links_m)) {
#     m <- links_m[i,]
#     if (is.na(mat_m[m$i, m$j])) {
#     } else {
#         circos.link(m$g1, m$k1 - 0.5, m$g2, m$k2 - 0.5, lwd = mat_m[m$i, m$j]*2, h = (abs(m$j-m$i)/(nm-1)*0.7),
#                     col = col_fun(mat_m[m$i, m$j]))
#     }
# }


# 4. done
# update the plot region. after setting the region in initalized, create andother class list (e.g. no PCBsF)
# and have it to draw the circle withtout PCBsF, then I can update PCBsF with thicker track


# circos.trackPlotRegion(
#     ylim = c(0, 1), 
#     bg.border = NA,
#     bg.col = c(colorCodes, rev(colorCodes)),
#     track.height = 0.04,
#     panel.fun = function(x, y) {
#         nm = get.cell.meta.data("sector.index")
#         r = cl[[nm]]
#         n = length(r)
#         if (nm == "PCBsF") {
#         } else {
#             circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(1, n), lwd = 0.5)
#             # circos.text(1:n - 0.5, rep(0.5, n), r, niceFacing = TRUE, cex = 0.7)  #cell text
#             circos.text(n/2, 1.5, nm, adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = 0.8)
#         }
#     }
# )
# 
# 
# circos.updatePlotRegion(sector.index = "PCBsF", track.index = 1, bg.col = NA, bg.border = NA)
# r = cl[["PCBsF"]]
# n = length(r)
# circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(2, n), lwd = 0.5, col = "#89C5DA")
# circos.text(n/2, 3, "PCBsF", adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = 0.8)
# 



# ******** -----
# X. Errors -----------------------------------------------------------------

# Note: 1 point is out of plotting region in sector 'PCBsF', track '1'.
# Can ignore it, it will still draw. It's just saying it's outside the region

# ******** -----
# W. Learning -----------------------------------------------------------------

######################################################### #
## \\ W1. Learning.[Color_space_gradient_comparison]  ----
######################################################### #
# library(circlize)
# 
# space = c("RGB", "HSV", "LAB", "XYZ", "sRGB", "LUV")
# par(xpd = NA)
# plot(NULL, xlim = c(-1, 1), ylim = c(0, length(space)+1), type = "n")
# for(i in seq_along(space)) {
#     f = colorRamp2(c(-1, 0, 1), c("green", "black", "red"), space = space[i])
#     x = seq(-1, 1, length = 200)
#     rect(x-1/200, i-0.5, x+1/200, i+0.5, col = f(x), border = NA)
#     text(1, i, space[i], adj = c(-0.2, 0.5))
# }
# par(xpd = FALSE)
# 
# ######################################################### #
# ## \\ W1. Learning.[Gradient_ramp_with_transparency]  ----
# ######################################################### #
# 
# myColours = c(1, "steelblue", "#FFBB00", rgb(0.4, 0.2, 0.3))
# adjustcolor(myColours, 0.4)

# also see my firefox bookmark addalpha and colorramppalettealpha package

# most color space has no transpancy (alpha channel). it's addon
# HEX code with 2 more number in front to indicate
# !! HEX code is *NOT* a color space as c("RGB", "HSV", "LAB", "XYZ", "sRGB", "LUV"), but a notation
# RGBA color space has alpha
# R can read HEX and thus with alpha; but colorramp etc use the color space without alpha component


## rect example
# 
# rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
# cat("\014")   # same as ctrl-L
# require(grDevices)
# ## set up the plot region:
# op <- par(bg = "thistle")
# plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",
#      main = "2 x 11 rectangles; 'rect(100+i,300+i,  150+i,380+i)'")
# the above just set the title, nth to do with "i"

# i <- 4*(0:10)
# ## draw rectangles with bottom left (100, 300)+i
# ## and top right (150, 380)+i
# rect(100+i, 300+i, 150+i, 380+i, col = rainbow(11, start = 0.7, end = 0.1))
# rect(240-i, 320+i, 250-i, 410+i, col = heat.colors(11), lwd = i/5)




# ******** -----
# X. Playground ----------------------------------------------------------------

# X1 rou  -------------------------------------------
# 
# mat_m = as.matrix(corrubinm) 
# 
# mat_m[mat_m < 0.25 & mat_m > -0.25] <- NA  # set NA threshold
# 
# nm = nrow(mat_m)
# rnm = rownames(mat_m)
# 
# for(i in 1:(nm-1)) {
#     for(j in seq(i+1, nm)) {
#         g1 = cd[rnm[i]]
#         g2 = cd[rnm[j]]
#         r1 = cd[cd == g1]
#         k1 = which(names(r1) == rnm[i])
#         r2 = cd[cd == g2]
#         k2 = which(names(r2) == rnm[j])
#         if (is.na(mat_m[i, j])) {
#             # circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.7),
#             #            col = "#00000000") # set transparet
#         } 
#         else {
#             if (all(c("PCBsM", "POPsM") %in% c(g1, g2))){
#                 circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.7),
#                             rou = 0.76,
#                             col = col_fun(mat_m[i, j]))    
#             }
#             else {
#                 circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.7),
#                             col = col_fun(mat_m[i, j]))
#             }
#         }
#     }
# }
# circos.clear()