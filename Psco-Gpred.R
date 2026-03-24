
-------------------------------------------------------------------
  # title: "Pseudophryne corroboree Genomic Predictions"
  # author: "Mikaeylah Davidson"
  -------------------------------------------------------------------
  
  # All data is stored within the Psco-Gpred project
  
  # Packages ----
if (!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, OpenMx, gplots, dplyr, corrplot, ggplot2, data.table, jsonlite, 
               asreml, ASRgenomics, BiocManager, ggrepel, ASRgwas, ASRtools, cowplot,
               MASS, qqman, gaston, reshape2, scico, pcadapt, randomcoloR, RColorBrewer,
               netview, htmlwidgets,  pheatmap
)
devtools::install_github("esteinig/netview")


#-------------------------------------------------------------------
#   PHENOTYPE FILES ----
#-------------------------------------------------------------------
# Read in phenotype file
pheno <- read.table("data/phenotypes-895.txt", header=TRUE)

# Mutate data 
pheno <- pheno %>%
  mutate(across(c(ID, zoo, experiment), as.factor))

# create new variable, with ZOO as numeric
# MZ = 1
# TZ = 2
pheno <- pheno %>%
  mutate(zoo_numeric = ifelse(zoo == "MZ", 1, ifelse(zoo == "TZ", 2, NA))) 

#-------------------------------------------------------------------
#   GRM ----
#-------------------------------------------------------------------

#--- Gaston GRM ----
# Load bed matrix 
bedM <- read.bed.matrix("data/filteredPsco.37768") 
# Create GRM
gaston_GRM <-GRM(bedM) 

# save gaston GRM as .csv to look at in excel 
write.table(gaston_GRM, file = "output/gaston-grm.csv", sep = ",")

#- Inverse GRM 
gaston_gInv <- ginv(gaston_GRM)

# Add IDs to G inverse
# Need to load pheno file before this
attr(gaston_gInv,"rowNames") <- as.character(pheno$ID)
attr(gaston_gInv,"colNames") <- as.character(pheno$ID)
dimnames(gaston_gInv)[[1]] <- as.character(pheno$ID)
dimnames(gaston_gInv)[[2]] <- as.character(pheno$ID)
attr(gaston_gInv,"INVERSE") <- TRUE

#- Dominance matrix ----
gaston_DOM <-DM(bedM) 
dim(gaston_DOM)

#- Inverse GRM 
gaston_dInv <- ginv(gaston_DOM)

# Add IDs to G inverse
# Need to load pheno file before this
attr(gaston_dInv,"rowNames") <- as.character(pheno$ID)
attr(gaston_dInv,"colNames") <- as.character(pheno$ID)
dimnames(gaston_dInv)[[1]] <- as.character(pheno$ID)
dimnames(gaston_dInv)[[2]] <- as.character(pheno$ID)
attr(gaston_dInv,"INVERSE") <- TRUE


#--- Plot GRM as heat map ---- 
heatmap.2(gaston_GRM, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

# Customise heatmap & save 
custom_col <- colorRampPalette(c("white", "#b6d9f6", "#4b6b85" ))(10)  

min_val <- min(gaston_GRM, na.rm = TRUE)
max_val <- max(gaston_GRM, na.rm = TRUE)

GRM_heatmap <- pheatmap(gaston_GRM,
                        col = custom_col,            
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_rownames = FALSE,
                        show_colnames = FALSE,
                        legend_breaks = c(min_val, (min_val + max_val)/2, max_val),
                        legend_labels = c("0", "0.5", "1")
)


pdf("figures/formatted/FigureS5-GRM.pdf", width = 8, height = 8)
GRM_heatmap
dev.off()


#--- GRM Variance ----
# For GCTA-GREML power calculator, need the SNP derived genetic relationship variance
# This is the variance of the off diagonal elements 

# Extract upper triangular elements
off_diagonal_elements <- gaston_GRM[upper.tri(gaston_GRM, diag = FALSE)]  

# Calculate the variance of the off-diagonal elements
var(off_diagonal_elements)
# 0.009229
# Variance for power caluclator 



#----------------------------------------------------------------
#    HERITABILITY - ASREML ----
#----------------------------------------------------------------
# ASREML is a licensed product available at: https://asreml.kb.vsni.co.uk 
# the script is provided to conduct the analysis, but a license is required and not provided here. 
#----------------------
#   SURVIVAL ----
#----------------------
# DATA ----
# this is the slope of the log transformed Bd loads 
hist(pheno$days_survived)
summary(pheno$days_survived)
#--- NUL----
daysSurv_nul <- asreml(fixed=days_survived~1,
                       random=~vm(ID, gaston_gInv),
                       residual=~idv(units),
                       data=pheno,
                       ai.sing=TRUE,
                       workspace="2gb", 
                       threads=16)

plot(daysSurv_nul)
# variance components
summary(daysSurv_nul)$varcomp
# fixed effects
summary(daysSurv_nul, coef=TRUE)$coef.fixed

# heritabilty
vpredict(daysSurv_nul,h2~V1/(V1+V2))$Estimate 
# SE of heritability
round(vpredict(daysSurv_nul,h2~(V1)/(V1+V2))$SE,3)

# calculate GEBV ----
gebv_dsurv_nul <- as.data.frame(summary(daysSurv_nul,coef=T)$coef.random[,1])
gebv_dsurv_nul
summary(gebv_dsurv_nul)
# plot
plot(summary(daysSurv_nul,coef=T)$coef.random[,1])

# write GEBVs out
write.csv(gebv_dsurv_nul, "output/ASRemlgebv_slope_nul.csv", row.names = FALSE)




#--- ZOO ----
daysSurv_zoo <- asreml(fixed=days_survived~1 + zoo,
                       random=~vm(ID,gaston_gInv),
                       residual=~idv(units),
                       data=pheno,
                       ai.sing=TRUE,
                       workspace="2gb", 
                       threads=16)
plot(daysSurv_zoo)
# variance components
summary(daysSurv_zoo)$varcomp
# fixed effects
summary(daysSurv_zoo, coef=TRUE)$coef.fixed

# habitability 
vpredict(daysSurv_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(daysSurv_zoo,h2~(V1)/(V1+V2))$SE,3)

# Phenotypic variation  = genetic+residual variance 
35.23027 + 301.06051
# Get phenotypic SE = sqrt of genetic SE ^2 + sqrt residual SE ^2
sqrt(14.19^2 + 17.21^2)


#--- AGE ----
daysSurv_age <- asreml(fixed=days_survived~1 + age_start,
                       random=~vm(ID,gaston_gInv),
                       residual=~idv(units),
                       data=pheno,
                       ai.sing=TRUE,
                       workspace="2gb", 
                       threads=16)
plot(daysSurv_age)
# variance components
summary(daysSurv_age)$varcomp
# fixed effects
summary(daysSurv_age, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(daysSurv_age,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(daysSurv_age,h2~(V1)/(V1+V2))$SE,3)









#--- AGE + ZOO ----
daysSurv_age_zoo <- asreml(fixed=days_survived~ age_start + zoo,
                           random=~vm(ID,gaston_gInv),
                           residual=~idv(units),
                           data=pheno,
                           ai.sing=TRUE,
                           workspace="2gb", 
                           threads=16)
plot(daysSurv_age_zoo)
# variance components
summary(daysSurv_age_zoo)$varcomp
# fixed effects
summary(daysSurv_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(daysSurv_age_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(daysSurv_age_zoo,h2~(V1)/(V1+V2))$SE,3)
# GEBV ----
gebv_daysSurv_age_zoo <- as.data.frame(summary(daysSurv_age_zoo,coef=T)$coef.random[,1])
gebv_daysSurv_age_zoo
summary(gebv_daysSurv_age_zoo)
# plot
plot(summary(summary(daysSurv_age_zoo,coef=T)$coef.random[,1]))

# write GEBVs out
write.csv(gebv_daysSurv_age_zoo, "output/ASReml-gebv_dSurv-AgeZoo.csv", row.names = FALSE)




#--- DOMINANCE ----
#--- Age + Zoo ----
# Add both the additive and the dominance inverse GRMs
DdaysSurv_age_zoo <- asreml(fixed=days_survived~ age_start + zoo,
                            random=~vm(ID, gaston_gInv) + vm(ID, gaston_dInv),
                            residual=~idv(units),
                            data=pheno,
                            ai.sing=TRUE,
                            workspace="2gb", 
                            threads=16)
plot(DdaysSurv_age_zoo)
# variance components
summary(DdaysSurv_age_zoo)$varcomp
# fixed effects
summary(DdaysSurv_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
# additive 
vpredict(DdaysSurv_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$Estimate
# dominance 
vpredict(DdaysSurv_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$Estimate

# SE of heritability 
# Additive 
round(vpredict(DdaysSurv_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$SE, 3)
# Dominance 
round(vpredict(DdaysSurv_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$SE, 3)




#----------------------
#  PATHOGEN LOAD ----
#----------------------
# CHECK DATA ----
# this is the slope of the log transformed bd loads 
hist(pheno$slope_bd)


#--- NUL ----
slope_nul <- asreml(fixed=slope_bd~1,
                    random=~vm(ID,gaston_gInv),
                    residual=~idv(units),
                    data=pheno,
                    ai.sing=TRUE,
                    workspace="2gb", 
                    threads=16)
plot(slope_nul)
# variance components
summary(slope_nul)$varcomp
# fixed effects
summary(slope_nul, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(slope_nul,h2~V1/(V1+V2))$Estimate # = Va 
# SE of heritability 
round(vpredict(slope_nul,h2~(V1)/(V1+V2))$SE,3)

# random effects = for each individual as only random effect is ID
summary(slope_nul, coef=TRUE)$coef.random

# Incremental wald tests
wald.asreml(slope_nul, denDF='algebraic')$wald



# GEBV ----
gebv_slope_nul <- as.data.frame(summary(slope_nul,coef=T)$coef.random[,1])
gebv_slope_nul
summary(gebv_slope_nul)
# plot
plot(summary(slope_nul,coef=T)$coef.random[,1])

# write GEBVs out
write.csv(gebv_slope_nul, "output/ASReml-gebv_slope_null.csv", row.names = FALSE)







#--- ZOO ----
slope_zoo <- asreml(fixed=slope_bd~1 + zoo,
                    random=~vm(ID,gaston_gInv),
                    residual=~idv(units),
                    data=pheno,
                    ai.sing=TRUE,
                    workspace="2gb", 
                    threads=16)
plot(slope_zoo)
# variance components
summary(slope_zoo)$varcomp
# fixed effects
summary(slope_zoo, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(slope_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(slope_zoo,h2~(V1)/(V1+V2))$SE,3)








#--- AGE ----
slope_age <- asreml(fixed=slope_bd~1 + age_start,
                    random=~vm(ID,gaston_gInv),
                    residual=~idv(units),
                    data=pheno,
                    ai.sing=TRUE,
                    workspace="2gb", 
                    threads=16)
plot(slope_age)
# variance components
summary(slope_age)$varcomp
# fixed effects
summary(slope_age, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(slope_age,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(slope_age,h2~(V1)/(V1+V2))$SE,3)









#--- AGE + ZOO ----
slope_age_zoo <- asreml(fixed=slope_bd~ age_start + zoo,
                        random=~vm(ID,gaston_gInv),
                        residual=~idv(units),
                        data=pheno,
                        ai.sing=TRUE,
                        workspace="2gb", 
                        threads=16)
plot(slope_age_zoo)
# variance components
summary(slope_age_zoo)$varcomp
# fixed effects
summary(slope_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(slope_age_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(slope_age_zoo,h2~(V1)/(V1+V2))$SE,3)

# GEBV ----
gebv_slope_age_zoo <- as.data.frame(summary(slope_age_zoo,coef=T)$coef.random[,1])
gebv_slope_age_zoo
summary(gebv_slope_age_zoo)
# plot
plot(summary(slope_age_zoo,coef=T)$coef.random[,1])

# write GEBVs out
write.csv(gebv_slope_age_zoo, "output/ASReml-gebv_slope-AgeZoo.csv", row.names = FALSE)



#--- DOMINANCE ----
#--- Age + Zoo ----
# Add both the additive and the dominance inverse GRMs
Dslope_age_zoo <- asreml(fixed=slope_bd~ age_start + zoo,
                         random=~vm(ID, gaston_gInv) + vm(ID, gaston_dInv),
                         residual=~idv(units),
                         data=pheno,
                         ai.sing=TRUE,
                         workspace="2gb", 
                         threads=16)
plot(Dslope_age_zoo)
# variance components
summary(Dslope_age_zoo)$varcomp
# fixed effects
summary(Dslope_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
# additive 
vpredict(Dslope_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$Estimate
# dominance 
vpredict(Dslope_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$Estimate

# SE of heritability 
# Additive 
round(vpredict(Dslope_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$SE, 3)
# Dominance 
round(vpredict(Dslope_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$SE, 3)


#----------------------
#  BODY CONDITION ----
#----------------------
# DATA ----
hist(pheno$bc_start)

#--- NUL ----
bc_nul <- asreml(fixed=bc_start~1,
                 random=~vm(ID,gaston_gInv),
                 residual=~idv(units),
                 data=pheno,
                 ai.sing=TRUE,
                 workspace="2gb", 
                 threads=16)
plot(bc_nul)
# variance components
summary(bc_nul)$varcomp
# fixed effects
summary(bc_nul, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(bc_nul,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(bc_nul,h2~(V1)/(V1+V2))$SE,3)


#--- ZOO ----
bc_zoo <- asreml(fixed=bc_start~1 + zoo,
                 random=~vm(ID,gaston_gInv),
                 residual=~idv(units),
                 data=pheno,
                 ai.sing=TRUE,
                 workspace="2gb", 
                 threads=16)
plot(bc_zoo)
# variance components
summary(bc_zoo)$varcomp
# fixed effects
summary(bc_zoo, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(bc_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(bc_zoo,h2~(V1)/(V1+V2))$SE,3)




#--- AGE ----
bc_age <- asreml(fixed=bc_start~1 + age_start,
                 random=~vm(ID,gaston_gInv),
                 residual=~idv(units),
                 data=pheno,
                 ai.sing=TRUE,
                 workspace="2gb", 
                 threads=16)
plot(bc_age)
# variance components
summary(bc_age)$varcomp
# fixed effects
summary(bc_age, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(bc_age,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(bc_age,h2~(V1)/(V1+V2))$SE,3)



#--- AGE + ZOO ----
bc_age_zoo <- asreml(fixed=bc_start~ age_start + zoo,
                     random=~vm(ID,gaston_gInv),
                     residual=~idv(units),
                     data=pheno,
                     ai.sing=TRUE,
                     workspace="2gb", 
                     threads=16)
plot(bc_age_zoo)
# variance components
summary(bc_age_zoo)$varcomp
# fixed effects
summary(bc_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
vpredict(bc_age_zoo,h2~V1/(V1+V2))$Estimate 
# SE of heritability 
round(vpredict(bc_age_zoo,h2~(V1)/(V1+V2))$SE,3)


# GEBV ----
gebv_bc_age_zoo <- as.data.frame(summary(bc_age_zoo,coef=T)$coef.random[,1])
gebv_bc_age_zoo
summary(gebv_bc_age_zoo)
# plot
plot(summary(bc_age_zoo,coef=T)$coef.random[,1])

# write GEBVs out
write.csv(gebv_bc_age_zoo, "output/ASReml-gebv_BC-AgeZoo.csv", row.names = FALSE)




#--- DOMINANCE ----
#--- Age + Zoo ----
# Add both the additive and the dominance inverse GRMs
Dbc_age_zoo <- asreml(fixed=bc_start~ age_start + zoo,
                      random=~vm(ID, gaston_gInv) + vm(ID, gaston_dInv),
                      residual=~idv(units),
                      data=pheno,
                      ai.sing=TRUE,
                      workspace="2gb", 
                      threads=16)
plot(Dbc_age_zoo)
# variance components
summary(Dbc_age_zoo)$varcomp
# fixed effects
summary(Dbc_age_zoo, coef=TRUE)$coef.fixed

# heritabilty 
# additive 
vpredict(Dbc_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$Estimate
# dominance 
vpredict(Dbc_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$Estimate

# SE of heritability 
# Additive 
round(vpredict(Dbc_age_zoo, h2_additive ~ V1/(V1 + V2 + V3))$SE, 3)
# Dominance 
round(vpredict(Dbc_age_zoo, h2_dominance ~ V2/(V1 + V2 + V3))$SE, 3)

#----------------------------------------------------------
#   TRAIT CORRELATIONS ---- 
#----------------------------------------------------------
# Phenotypic and genetic correlations for
  # Days survived
  # Slope of infection
  # Infection status 
  # Body condition

# pathogen load and survival ----
plot(pheno$slope_bd, pheno$days_survived)

slope_survival <- asreml(fixed=cbind(slope_bd, days_survived)~trait,
                         random=~us(trait, init = c(1, 0.1, 1)):vm(ID, gaston_gInv), # genetic variance
                         coruh = ~ id(units):us(trait, init = c(1, 0.01, 1)), # residual variance 
                         ai.sing = TRUE,
                         data=pheno,
                         workspace="2gb", 
                         threads=16)

summary(slope_survival)$varcomp
vpredict(slope_survival,rg~V2/sqrt(V1*V3)) # genetic correlation 
vpredict(slope_survival,re~V6/sqrt(V5*V7)) # environmental correlation
vpredict(slope_survival,rp~(V2 + V6)/sqrt((V1 + V5) * (V3 + V7))) # phenotypic correlation. 

# pathogen load and bc ----
plot(pheno$slope_bd, pheno$bc_start)

slope_bc <- asreml(fixed=cbind(slope_bd, bc_start)~trait,
                   random=~us(trait, init = c(1, 0.1, 1)):vm(ID,gaston_gInv), # genetic variance
                   coruh = ~ id(units):us(trait, init = c(1, 0.01, 1)), # residual variance 
                   ai.sing = TRUE,
                   data=pheno,
                   maxit=200)

summary(slope_bc)$varcomp
vpredict(slope_bc,rg~V2/sqrt(V1*V3)) # genetic correlation 
vpredict(slope_bc,re~V6/sqrt(V5*V7)) # environmental correlation
vpredict(slope_bc,rp~(V2 + V6)/sqrt((V1 + V5) * (V3 + V7))) # phenotypic correlation. 



# survival and body condition ----
plot(pheno$days_survived, pheno$bc_start)

survival_bc <- asreml(fixed=cbind(days_survived, bc_start)~trait,
                      random=~us(trait, init = c(1, 0.1, 1)):vm(ID,gaston_gInv), # genetic variance
                      coruh = ~ id(units):us(trait, init = c(1, 0.01, 1)), # residual variance 
                      ai.sing = TRUE,
                      data=pheno)

summary(survival_bc)$varcomp
vpredict(survival_bc,rg~V2/sqrt(V1*V3)) # genetic correlation 
vpredict(survival_bc,re~V6/sqrt(V5*V7)) # environmental correlation
vpredict(survival_bc,rp~(V2 + V6)/sqrt((V1 + V5) * (V3 + V7))) # phenotypic correlation. 


#--- Plot correlations ----
# Summarised the correlation in excell for all traits 
traitCor <- read_xlsx("output/Trait-correlations.xlsx", sheet = "trait-cor")

#--- Genetic correlations ----
rG <- read_xlsx("output/Trait-correlations.xlsx", sheet = "rG")
# Remove the first column 
rG <- rG[, -1] 
# Set row names
rownames(rG) <- colnames(rG)
# Convert to matrix 
rG <- as.matrix(rG)

# Genetic Correlation Plot
corrplot(rG, 
         method="color", 
         col = rev(scico(100, palette = "vik")),
         type = "lower",   
         title="Genetic Correlation (rG)", 
         tl.col="black", 
         cl.pos="b",
         sig.level = 0.05,                     
         insig = "blank",                      
         tl.srt = 45,                          
         addCoef.col = "black",    
         mar=c(0,0,5,0), 
         diag = FALSE)



#--- Phenotypic correlations ----
rP <- read_xlsx("output/Trait-correlations.xlsx", sheet = "rP")
# Remove the first column 
rP <- rP[, -1] 
# Set row names
rownames(rP) <- colnames(rP)
# Convert to matrix 
rP <- as.matrix(rP)
# Phenotypic Correlation Plot
corrplot(rP, 
         method="color", 
         col = rev(scico(100, palette = "vik")),
         type = "lower",  
         title="Phenotypic Correlation (rP)", 
         tl.col="black", 
         cl.pos="b",
         sig.level = 0.05,                     
         insig = "blank",                      
         tl.srt = 45,                          
         addCoef.col = "black",
         mar=c(0,0,5,0), 
         diag = FALSE)


#--- Environmental correlations ----
rE <- read_xlsx("output/Trait-correlations.xlsx", sheet = "rE")
# Remove the first column 
rE <- rE[, -1] 
# Set row names
rownames(rE) <- colnames(rE)
# Convert to matrix 
rE <- as.matrix(rE)                   
# Environmental Correlation Plot
corrplot(rE, 
         method="color", 
         col = rev(scico(100, palette = "vik")),
         type = "lower",   
         title="Environmental Correlation (rE)", 
         tl.col="black", 
         cl.pos="b",
         sig.level = 0.05,                     
         insig = "blank",                      
         tl.srt = 45,                          
         addCoef.col = "black",   
         mar=c(0,0,5,0), 
         diag = FALSE)


#--- Combine phenotype and genotype correlations ----
# Environmental correlations 
rPrG <- read_xlsx("output/Trait-correlations.xlsx", sheet = "rPtop-rGlower")
# Remove the first column 
rPrG <- rPrG[, -1] 
# Set row names
rownames(rPrG) <- colnames(rPrG)
# Convert to matrix 
rPrG <- as.matrix(rPrG)                   
# combine plot
corrplot(rPrG, 
         method = 'pie',
         col = rev(scico(100, palette = "vik")),
         title="rG lower half & rP upper half", 
         tl.col="black", 
         cl.pos="b",           # position of key
         sig.level = 0.05,                     
         tl.pos = 'd',       # position of labels                   
         mar=c(0,0,5,0))

# Set custom labels
# splitting labels over two lines
custom_labels <- c("Survival", "Pathogen\nload","Body\ncondition")  

# Update column and row names 
colnames(rPrG) <- custom_labels
rownames(rPrG) <- custom_labels

# Test for significance 
testRes = cor.mtest(rPrG, conf.level = 0.95)
# no sig correlations, so can remove from plot

# Plot 
corrplot(rPrG, 
         method = 'square',    
         col = rev(scico(100, palette = "vik")),
         tl.pos = 'd',
         tl.col="black", 
         cl.pos="b",
         cl.ratio = 0.2,                     
         p.mat = testRes$p,
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 1.2,
         insig = 'label_sig',                    
         pch.col = 'black',
         mar=c(0,0,0,0))

#--- Figure 2. Combine correlations w/ heritabilities ----
# Modify labs to include h2
custom_labels2 <- c(
  "Survival\ntime\nh² = 0.16 ± 0.047",
  "Pathogen\nload\n\nh² = 0.17 ± 0.049",
  "Body\ncondition\n\nh² = 0.41 ± 0.056"
)

# Update column and row names 
colnames(rPrG) <- custom_labels2
rownames(rPrG) <- custom_labels2

pdf("figures/Figure2-correlations.pdf", width = 8, height = 8)
corrplot(rPrG, 
         method = 'color', 
         addCoef.col = 'black',
         col = COL1(sequential = c('Blues'), n=1000),
         tl.pos = 'd',
         tl.col="black", 
         cl.pos="b",
         cl.ratio = 0.2, 
         addgrid.col = 'white',
         mar=c(0,0,0,0))
dev.off()



#-------------------------------------------------------------------
#   PCA ----
#-------------------------------------------------------------------
#~~~ All 895 animals ----
# BED file of all individuals 
psco_bedAll <- "data/filteredPsco.37768.bed"
psco_Allpcadapt <- read.pcadapt(psco_bedAll, type = "bed")
#~~~ Scree ------------------- 
# Proportion of variance 
psco_pcadapt_kplotAll <- pcadapt(input = psco_Allpcadapt, K = 30) 
plot(psco_pcadapt_kplotAll, option = "screeplot")

# Pull variance to plot manually 
singular_values <- psco_pcadapt_kplotAll$singular.values
# Square the singular values
squared_singular_values <- singular_values^2
# Calculate total variance
total_variance <- sum(squared_singular_values)
# Calculate the proportion of variance explained
variance_explained <- squared_singular_values / total_variance
# create df 
pc_number <- 1:length(variance_explained)  
scree_data <- data.frame(PC = pc_number, Variance_Explained = variance_explained)

# Plot 
pScree <- ggplot(scree_data, aes(x = PC, y = Variance_Explained * 100)) +  # Multiply by 100 for percentage
  labs( x = "Principal Component", y = "Variance Explained (%)") +
  geom_line(size = 0.5, colour = "lightblue") +
  geom_point(size = 1.5, colour = '#6baed6'  ) + 
  theme_cowplot() 

pScree




#~~~ Scores ------------------- 
plot(psco_pcadapt_kplotAll, option = "scores", i = 1 , j = 2)
# Pull scores so can plot manually
pca_scoresAll <- as.data.frame(psco_pcadapt_kplotAll$scores)  
# Rename cols 
# Only need to change the first two cols as using only using PC 1 and 2
colnames(pca_scoresAll) <- c("PC1", "PC2")  
# Plot
ggplot(pca_scoresAll, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.5 ) + 
  labs( x = "PC1", y = "PC2") +      
  theme_cowplot()




# Pull variance to plot manually 
singular_values <- psco_pcadapt_kplotAll$singular.values
# Square the singular values
squared_singular_values <- singular_values^2
# Calculate total variance
total_variance <- sum(squared_singular_values)
# Calculate the proportion of variance explained
variance_explained <- squared_singular_values / total_variance
print(variance_explained)
# PC 1 = 12.88%
# PC 2 = 8.77%



# ~~~ Highlight PCA from breeding simulations--------
# Read .fam file to get individual IDs
fam_file <- sub("\\.bed$", ".fam", psco_bedAll)
fam <- read.table(fam_file, stringsAsFactors = FALSE)
individual_ids <- fam$V2  # Column 2 is the individual ID

# Add IDs to PCA scores
pca_scoresAll$ID <- individual_ids

# Subset PCA to only have PC1 and PC2 and ID info
pca_12 <- pca_scoresAll[, c(1,2,31)]
## 200 Matings ----
# ~~~ OCS 6  --------
# Highlight file
highlight_OCS <- "data/AlphaMate/200matings-OCS/ContributorsModeOptTarget6.txt"
# Read the file
OCS_data <- read.table(highlight_OCS, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_OCSids <- OCS_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightOCS <- ifelse(pca_12$ID %in% highlight_OCSids, "highlight", "other")

#plot 
pScoreOCS6 <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightOCS)) +
  geom_point(
    data = subset(pca_12, highlightOCS == "other"),
    color = "black", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightOCS == "highlight"),
    color = "blue", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreOCS6 

# check numbers highlighted 
table(pca_12$highlightOCS) #highlight = 323   other = 572 

#### MULTIPLE Matings
# Amend colours, so males who mate multuple times are in different colours. 
# col 1 = ID, col 7 = number of matings
mate_info <- OCS_data[, c(1, 7)]
colnames(mate_info) <- c("ID", "MateCount")

# Merge mating info into PCA data
pca_12 <- merge(pca_12, mate_info, by = "ID", all.x = TRUE)

# Create a new column for color group — NA if not in highlight_data
pca_12$MateGroupPlot <- as.character(pca_12$MateCount)
pca_12$MateGroupPlot[!pca_12$MateGroupPlot %in% c("1", "2", "3")] <- "Other"
pca_12$MateGroupPlot <- factor(
  pca_12$MateGroupPlot,
  levels = c("1", "2", "3", "Other")
)

# Define colours for number of matings 
mate_colors <- c(
  "1" = "#FEC44F",
  "2" = "#FE9929",
  "3" = "#CC4C02",
  "Other" = "black"
)

# Plot
ggplot(pca_12, aes(x = PC1, y = PC2, color = MateGroupPlot)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = mate_colors, name = "No. Matings") +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()


# Check number highlighted
table(pca_12$MateGroupPlot)
# 1      2       3    Other 
# 255    59      9    572 



# ~~~ Maximum criterion/gains  --------
# Maximum Criterion file
highlight_MaxCri <- "data/AlphaMate/200matings-maxCri/ContributorsModeMaxCriterion.txt"
# Read the file
MaxCri_data <- read.table(highlight_MaxCri, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_Maxids <- MaxCri_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMaxCri <- ifelse(pca_12$ID %in% highlight_Maxids, "highlight", "other")

#plot 
pScoreMaxCri <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightMaxCri)) +
  geom_point(
    data = subset(pca_12, highlightMaxCri == "other"),
    color = "black", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightMaxCri == "highlight"),
    color = "#FEC44F", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMaxCri

# check numbers highlighted 
table(pca_12$highlightMaxCri) #highlight = 267   other = 628 


# ~~~ Minimum inbreeding  --------
# Min inbreed file
highlight_MinInb <- "data/AlphaMate/200matings-minInbreed/ContributorsModeMinInbreeding.txt"
# Read the file
MinInb_data <- read.table(highlight_MinInb, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_MinInbids <- MinInb_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMinInb <- ifelse(pca_12$ID %in% highlight_MinInbids, "highlight", "other")

#plot 
pScoreMinInb <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightMinInb)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(
    values = c("highlight" = "#FEC44F", "other" = "black"),
    guide = "none"  # Hide legend
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMinInb

# check numbers highlighted 
table(pca_12$highlightMinInb) #highlight = 393   other = 502 


# ~~~ Minimum coancestry  --------
# Min inbreed file
highlight_MinCoA <- "data/AlphaMate/200matings-OCS/ContributorsModeMinCoancestry.txt"
# Read the file
MinCoA_data <- read.table(highlight_MinCoA, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_MinCoAids <- MinCoA_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMinCoA <- ifelse(pca_12$ID %in% highlight_MinCoAids, "highlight", "other")

#plot 
pScoreMinCoA <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightMinCoA)) +
  geom_point(
    data = subset(pca_12, highlightMinCoA == "other"),
    color = "black", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightMinCoA == "highlight"),
    color = "#FEC44F", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMinCoA

# check numbers highlighted 
table(pca_12$highlightMinCoA) #highlight = 392   other = 503 


# ~~~ EBVs - top 10%  --------
# EBV file
ebv_data <- read_excel("output/EBVs.xlsx", sheet = "EBVs")

# Calculate the 80th percentile threshold for EBV_days_survived
threshold <- quantile(ebv_data$EBV_days_survived, probs = 0.8, na.rm = TRUE)
# 90th percentile
#threshold <- quantile(ebv_data$EBV_days_survived, probs = 0.9, na.rm = TRUE)

# Filter to get the top % IDs
top_ids <- ebv_data %>%
  filter(EBV_days_survived >= threshold) %>%
  pull(ID)  

# Add highlight column based on whether the ID is in the top 10% list
pca_12$highlightEBVs <- ifelse(pca_12$ID %in% top_ids, "highlight", "other")

# Plot
pScoreEBV20 <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightEBVs)) +
  geom_point(
    data = subset(pca_12, highlightEBVs == "other"),
    color = "black", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightEBVs == "highlight"),
    color = "#FEC44F", alpha = 0.8, size = 1
  ) +  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreEBV20

# check numbers highlighted 
table(pca_12$highlightEBVs) #highlight = 179   other = 716 





# ~~~ Combine Scores plots ----
#pPCAcombined <- plot_grid(pScoreMaxCri, pScoreMinInb, pScoreOCS6, pScoreEBV10,
 #                         labels = "AUTO", rel_widths = c(1, 1.2))

pPCAcombined <- plot_grid(
  pScoreMaxCri  + ggtitle("MaxGains"),
  pScoreMinInb  + ggtitle("MinInbreeding"),
  pScoreOCS6    + ggtitle("OCS-6"),
  pScoreEBV20   + ggtitle("EBVs 10%"),
  labels = "AUTO",
  rel_widths = c(1, 1.2)
)

pPCAcombined


## 100 Matings ----
# ~~~ OCS 4 - 100 matings --------
# Highlight file
highlight_OCS100 <- "data/AlphaMate/100matings-OCS/ContributorsModeOptTarget4.txt"
# Read the file
OCS100_data <- read.table(highlight_OCS100, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_OCSids100 <- OCS100_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightOCS100 <- ifelse(pca_12$ID %in% highlight_OCSids100, "highlight", "other")

#plot 
pScoreOCS4_100 <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightOCS100)) +
  geom_point(
    data = subset(pca_12, highlightOCS100 == "other"),
    color = "#D3D3D3", alpha = 1, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightOCS100 == "highlight"),
    color = "#407EC9", alpha = 1, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreOCS4_100 

# check numbers highlighted 
table(pca_12$highlightOCS100) #highlight = 155   other = 740 



# ~~~ Maximum criterion/gains --------
# Maximum Criterion file
highlight_MaxCri100 <- "data/AlphaMate/100matings-maxCri/ContributorsModeMaxCriterion.txt"
# Read the file
MaxCri100_data <- read.table(highlight_MaxCri100, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_Maxids100 <- MaxCri100_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMaxCri100 <- ifelse(pca_12$ID %in% highlight_Maxids100, "highlight", "other")

#plot 
pScoreMaxCri100 <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightMaxCri100)) +
  geom_point(
    data = subset(pca_12, highlightMaxCri100 == "other"),
    color = "#D3D3D3", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightMaxCri100 == "highlight"),
    color = "#407EC9", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMaxCri100

# check numbers highlighted 
table(pca_12$highlightMaxCri100) #highlight = 134   other = 761 






# ~~~ Minimum inbreeding --------
# Min inbreed file
highlight_MinInb100 <- "data/AlphaMate/100matings-minInbreed/ContributorsModeMinInbreeding.txt"
# Read the file
MinInb100_data <- read.table(highlight_MinInb100, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_MinInbids100 <- MinInb100_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMinInb100 <- ifelse(pca_12$ID %in% highlight_MinInbids100, "highlight", "other")

#plot 
pScoreMinInb100 <- ggplot(pca_12, aes(x = PC1, y = PC2, color = highlightMinInb100)) +
  geom_point(
    data = subset(pca_12, highlightMinInb100 == "other"),
    color = "#D3D3D3", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightMinInb100 == "highlight"),
    color = "#407EC9", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMinInb100

# check numbers highlighted 
table(pca_12$highlightMinInb100) #highlight = 189   other = 706 


# ~~~ Minimum coansestry --------
# Min inbreed file
highlight_MinCoA100 <- "data/AlphaMate/100matings-OCS/ContributorsModeMinCoancestry.txt"
# Read the file
MinCoA100_data <- read.table(highlight_MinCoA100, header = TRUE, stringsAsFactors = FALSE)
# Extract vector of IDs
highlight_MinCoA100ids <- MinCoA100_data[[1]]  

# Add highlight column: TRUE if in highlight list, FALSE otherwise
pca_12$highlightMinCoA100 <- ifelse(pca_12$ID %in% highlight_MinCoA100ids, "highlight", "other")

#plot 
pScoreMinCoA100 <- ggplot(pca_12, aes(x = PC1, y = PC2)) +
  geom_point(
    data = subset(pca_12, highlightMinCoA100 == "other"),
    color = "black", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = subset(pca_12, highlightMinCoA100 == "highlight"),
    color = "#FEC44F", alpha = 0.8, size = 1
  ) +
  labs(x = "PC1", y = "PC2") +
  theme_cowplot()

pScoreMinCoA100

# check numbers highlighted 
table(pca_12$highlightMinCoA100) #highlight = 199   other = 696 


# ~~~ Combine Scores plots 100 matings ----

pPCAcombined_100 <- plot_grid(
  pScoreMinInb100    + ggtitle("MinInbreed100"),
  pScoreMaxCri100  + ggtitle("MaxGains100"),
  pScoreOCS4_100    + ggtitle("OCS4-100"),
  labels = "AUTO",
  nrow = 1,  
  ncol = 3
)

pPCAcombined_100

# ggsave("figures/PCA100.pdf", plot = pPCAcombined_100, width = 12, height = 4)

# No titles:
pPCAcombined_100_notitle <- plot_grid(
  pScoreMinInb100,
  pScoreMaxCri100,
  pScoreOCS4_100,
  labels = "AUTO",
  nrow = 1,  
  ncol = 3
)

pPCAcombined_100_notitle


# Save
ggsave("figures/Figure3-PCA.pdf", plot = pPCAcombined_100_notitle, width = 12, height = 4)



# ~~~ Combine Scores plots 100 + 200 matings ----
pPCAcombined_12 <- plot_grid(
  pScoreMinCoA100    + ggtitle("MinCoAncestry 100"),
  pScoreMaxCri100  + ggtitle("MaxGains 100"),
  pScoreOCS4_100    + ggtitle("OCS5 100"),
  pScoreMinCoA    + ggtitle("MinCoAncestry 200"),
  pScoreMaxCri  + ggtitle("MaxGains 200"),
  pScoreOCS6    + ggtitle("OCS6 200"),
  labels = "AUTO"
)

pPCAcombined_12


# Save
ggsave("figures/PCA100+200.pdf", plot = pPCAcombined_12, width = 12, height = 6)


#----------------------------------------------------------------
#    EBVs ----
#----------------------------------------------------------------
# All EBVs caculated from ASReml ----
EBVs<- read_excel("output/EBVs.xlsx", sheet = "EBVs")

#--- Survival ----

EBV_daysurv <- read_excel("output/EBVs.xlsx", sheet = "ASReml_daysSurv")

# sort by EBVs
EBV_daysurv <- EBVs %>% 
  arrange(desc(EBV_days_survived))

# create rank to have the top X % of population
EBV_daysurv <- EBV_daysurv %>%
  mutate(
    rank = row_number(),                              
    total = n(),                                      
    top_percent = ceiling(rank / total * 100)        
  )


# Save as CSV to get number of animals per % for selection 
# write.table(EBV_daysurv, file = "output/EBV_daysSurvived-rankings.csv", sep = ",")


# Calculate the % of animals in the top X %, from each zoo 
# list of top % want to calculate 
top_x_percent_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50) 

top_percent_dsurv <- data.frame()

#---% of animals by zoo ----
# Loop through each top X% level
for (top_x_percent in top_x_percent_list) {
  # Filter data for the current top X% level
  top_animals <- EBV_daysurv %>%
    filter(top_percent <= top_x_percent)
  
  # Calculate the count and percentage for each zoo in the current top X% group
  result <- top_animals %>%
    group_by(zoo) %>%
    summarise(
      count_in_top = n(),                             
      total_in_top = nrow(top_animals),               
      percent_in_top = (count_in_top / total_in_top) * 100  
    ) %>%
    mutate(top_percent_group = paste0("Top ", top_x_percent, "%")) 
  
  # Combine the result with top_percent_dsurv data frame
  top_percent_dsurv <- bind_rows(top_percent_dsurv, result)
}

# Display the combined results
print(top_percent_dsurv)

#--- % of families ---- 
# Initialize data frame to store results
uniqueFams_dsurv <- data.frame()

# Loop through each top X% level
for (top_x_percent in top_x_percent_list) {
  # Filter data for the current top X% level
  top_animals <- EBV_daysurv %>%
    filter(top_percent <= top_x_percent)
  
  # Calculate the total unique families across all zoos
  total_unique_families <- top_animals %>%
    summarise(unique_families = n_distinct(fam)) %>%
    pull(unique_families)  # Extract the value as a scalar
  
  # Add the result to the output data frame
  result <- data.frame(
    top_percent_group = paste0("Top ", top_x_percent, "%"),
    total_unique_families = total_unique_families
  )
  
  # Combine the result with top_percent_dsurv data frame
  uniqueFams_dsurv <- bind_rows(top_percent_dsurv, result)
}

# Display the combined results
print(uniqueFams_dsurv)


#----------------------------------------------------------------
#    GENOME WIDE ASSOCIATION ANALYSIS ----
#----------------------------------------------------------------


#------------------------------
#   SIGNIFICANCE THRESHOLD ----
#------------------------------
# Calculate genome wide significance and suggested thresholds
# Number of SNPS = 37768 ----
# Number of unlinked SNPS = 5795 ----
# SNPs with a r2 < 0.2 

# Use unlinked SNPs

# Genomewide line ----
# Bonferonni correction (0.05) / number snps 
# 0.05/5795
# 8.628128e-06
# log transform 
genoSigLine <- -log10(8.628128e-06)
# 5.064083 == genome wide significance threshold 

# Suggestive line ----
# one false positive / number snps 
# 1/5795
# 0.0001725626
# log transform 
suggestSigLine <- -log10(0.0001725626)
# 3.763053 == Suggestive significance threshold 

#------------------------------
#   POWER ANALYSIS ----
#------------------------------
# https://github.com/kaustubhad/gwas-power
# download source script and upload to spartan
source("../../programs/power_calc_functions.R")

# n = Sample size (has to be >= 0).
# qsq = Fraction of trait variance explained by the SNP (has to be between 0 and 1). Denoted as q^2 or q-squared; other studies use h^2 too, to indicate it's similarity with heritability.
# pval = P-value threshold for significance. 

power_n_hsq(n = (1:9)*100, qsq = (1:10)/100, pval=0.000009)


#----------------------------------------------------------------
#    Genome-Wide Association Study ----
#----------------------------------------------------------------
# This analysis uses ASREML which requires a license - see above
# The code is provided, but the license to run the analysis is not. 

# GENOTYPES ---- 
geno <- read.bed.matrix("data/filteredPsco.37768.bed", verbose = TRUE)
# make map
map <- geno@snps[,c(2,1,4)]
# Rename 
colnames(map) <- c("marker","chrom","pos")
# Convert genotype file to matrix 
X <- as.matrix(geno)
X[1:5,1:5]
# Check individuals genotyped X no snps
dim(geno) # 895 37768

#-----------------------------
#   SURVIVAL ----
#-----------------------------
#--- Model ----
# Preparation for gwas 
pre.dSurvEBV <- pre.gwas(pheno.data=EBV_daysurv, 
                         indiv="ID",
                         resp="EBV_daysSurv",
                         fixedf = "zoo",
                         cov = "age_start",
                         geno.data=X, 
                         map.data = map,
                         Q.method="K",
                         method="VanRaden",
                         maf=0.01, 
                         marker.callrate=0.8,
                         ind.callrate=0.20, 
                         heterozygosity=0.8, 
                         Fis=0.98, 
                         impute=TRUE)

# Check values match, and what data cleaning has removed
dim(pre.dSurv$geno.data) # 895 37768

# Scree plot 
pre.dSurv$plot.scree

# Matches original data, with very low variation, =  7 looks okay... 
npc <-  7

# getting error with invalid numberic version, ASReml Tech team gave me this to input to stop version errors. 
Sys.setenv('_R_CHECK_STOP_ON_INVALID_NUMERIC_VERSION_INPUTS_'='false')

#--- GWAS ----
# Prepare data for gwas = checking, cleaning, imputing 
# Set p val threshold 
# Threads need to be used if on the HPC / on local not needed
gwasDaysSurv <- gwas.asreml(pheno.data=pre.dSurv$pheno.data, 
                            resp="days_survived",
                            fixedf = "zoo",
                            cov = "age_start",
                            gen="ID", 
                            Kinv=pre.dSurv$Kinv, 
                            Q=pre.dSurv$Q,
                            npc=7, 
                            geno.data=pre.dSurv$geno.data, 
                            map.data=pre.dSurv$map, 
                            pvalue.thr= 0.0001725626,
                            bonferroni=FALSE,
                            P3D=TRUE, 
                            #threads = 1
)
# A total of 3 markers identified based on a p-value threshold of 0.0001725626
# The observed FDR for the set of significant markers is: 217.2448%.

# Marker effects 
head(gwasDaysSurv$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasDaysSurv$gwas.all, file = "output/ASReml-gwas_DaysSurvived.csv", sep = ",")

# List of significant SNPs at the pval set in mod
gwasDaysSurv$gwas.sel 

#--- QQ plot ----
pdf("figures/ASReml/qq_daysSurvived.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasDaysSurv$gwas.all) +
  labs(title="Days survived") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


#--- Manhattan plot ----
manhattan.plot(gwas.table = gwasDaysSurv$gwas.all, 
               point.size = 1.2) + 
  geom_hline(yintercept =  genoSigLine, color = "darkgrey", linetype = "solid") + 
  geom_hline(yintercept =  suggestSigLine,  color = "darkgrey", linetype = "dashed") +
  ylim(0, 6)



pdf("figures/ASReml/manhattan_dSurv.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasDaysSurv$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  -log10(7.640587e-06) , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  -log10(0.0001528117) , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Days survived") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()





#--- USING EBVS ----
# EBV data is already loaded, and has zoo and age in file 
pre.dSurvEBV <- pre.gwas(pheno.data=EBV_daysurv, 
                         indiv="ID",
                         resp="EBV_daysSurv",
                         geno.data=X, 
                         map.data = map,
                         Q.method="K",
                         method="VanRaden",
                         maf=0.01, 
                         marker.callrate=0.8,
                         ind.callrate=0.20, 
                         heterozygosity=0.8, 
                         Fis=0.98, 
                         impute=TRUE)

# Check values match, and what data cleaning has removed
dim(pre.dSurvEBV$geno.data) # 895 37768

# Scree plot 
pre.dSurvEBV$plot.scree

# Matches original data, with very low variation, =  7 looks okay... 
npc <-  7

# Prepare data for gwas = checking, cleaning, imputing 
# Set p val threshold 
# Threads need to be used if on the HPC / on local not needed
gwasDaysSurvEBV <- gwas.asreml(pheno.data=pre.dSurvEBV$pheno.data, 
                               resp="EBV_daysSurv",
                               gen="ID", 
                               Kinv=pre.dSurvEBV$Kinv, 
                               Q=pre.dSurvEBV$Q,
                               npc=7, 
                               geno.data=pre.dSurvEBV$geno.data, 
                               map.data=pre.dSurvEBV$map, 
                               pvalue.thr= 8.628128e-06,         #0.0001725626,
                               bonferroni=FALSE,
                               P3D=TRUE, 
                               #threads = 1
)

#A total of 11 markers identified based on a p-value threshold of 8.628128e-06.
#The observed FDR for the set of significant markers is: 2.9624%.

# 122 markers identified based on a p-value threshold of 0.0001725626.
# The observed FDR for the set of significant markers is: 5.3421%. 
# Marker effects 
head(gwasDaysSurvEBV$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasDaysSurvEBV$gwas.all, file = "output/ASReml-gwas_DaysSurvived-EBV.csv", sep = ",")

# List of significant SNPs at the pval set in mod
gwasDaysSurvEBV$gwas.sel 

pdf("figures/ASReml/qq_daysSurvivedEBV.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasDaysSurvEBV$gwas.all) +
  labs(title="Days survived EBV") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


pdf("figures/ASReml/manhattan_dSurvEBV.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasDaysSurvEBV$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  -log10(7.640587e-06) , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  -log10(0.0001528117) , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Days survived using EBV") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()




#--- Look at chromosomes separately ---- 
# Need to pull the data out of the gwas df, and into own df & rename cols 
gwas_days_renamed <- gwasDaysSurvEBV$gwas.all %>%
  rename(CHR = chrom, 
         BP = pos, 
         P = p.value, 
         SNP = marker)


qqman::manhattan(subset(gwas_days_renamed, CHR == 6), 
                 ylim = c(0, 6), 
                 main = "Chr 6", 
                 suggestiveline = suggestSigLine, 
                 genomewideline = genoSigLine)

# format plot 
qqman::manhattan(subset(gwas_days_renamed, CHR == 6), 
                 main = "Days Survived EBV CHR 6", ylim = c(0, 6),
                 col = c("#648fb0", "#b6d9f6"), 
                 suggestiveline = FALSE,
                 genomewideline =  FALSE, 
                 cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
abline(h = genoSigLine, col = "darkgrey", lty = 1, lwd=2)  
abline(h = suggestSigLine, col = "darkgrey", lty = 2, lwd=2)  


# Look into regions ----
# Zoom into interesting region, and annotate SNPs
qqman::manhattan(subset(gwas_days_renamed, CHR == 7), 
                 ylim=c(0,6), main = "Chr 7", 
                 xlim = c(240000000,270000000), 
                 suggestiveline = suggestSigLine,
                 genomewideline = genoSigLine, 
                 annotatePval = 0.001, annotateTop = FALSE)

# Highlight markers with significant P value
qqman::manhattan(subset(gwas_days_renamed, CHR == 7), 
                 ylim=c(0,6), main = "Chr 7", 
                 annotatePval = 0.01, 
                 suggestiveline = suggestSigLine,
                 genomewideline = genoSigLine)

# Annotate SNPs that meet a threshold:
qqman::manhattan(gwas_days_renamed, 
                 annotatePval = 0.0001528117, 
                 annotateTop = TRUE, 
                 ylim=c(0,6),
                 suggestiveline = suggestSigLine,
                 genomewideline = genoSigLine)

#---  Format figure CHR7 ----
# A ----
# Format into DF to be able to use GGplot to plot for more control 

# Convert into DF not ASReml list
gwas_DS <- gwasDaysSurvEBV$gwas.all
# Transform p vals
gwas_DS$logP <- -log10(gwas_DS$p.value)  

# make chromosome sizes, to be able to plot as manhattan
chromosome_sizes <- gwas_DS %>%
  group_by(chrom) %>%
  summarize(chr_len = max(pos)) %>%
  mutate(chr_start = lag(cumsum(as.numeric(chr_len)), default = 0))

# Join chromsome sizes to df and calculate cumulative position
gwas_DS <- gwas_DS %>%
  left_join(chromosome_sizes, by = "chrom") %>%
  mutate(cum_pos = pos + chr_start)

# Calculate the midpoint for each chrom for labeling 
axis_df <- gwas_DS %>%
  group_by(chrom) %>%
  summarize(center = mean(cum_pos))

# Plot Manhattan 
Fig3A <- ggplot(gwas_DS, aes(x = cum_pos, y = -log10(p.value), color = as.factor(chrom))) +
  geom_hline(yintercept = genoSigLine, color = "grey40", linetype = "solid", linewidth = 0.6) +
  geom_hline(yintercept = suggestSigLine, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_point(alpha = 1, size = 1) +  
  scale_color_manual(values = rep(c("#648fb0", "#b6d9f6"), 
                                  length.out = length(unique(gwas_DS$chrom)))) +  
  labs(title = "", 
       x = "Chromosome", 
       y = expression(-log[10] (italic(p) * "-value"))
  ) +
  ylim(0, 6) + 
  scale_x_continuous(labels = axis_df$chrom, breaks = axis_df$center, expand = c(0, 0)) +  
  theme_cowplot() +
  theme(
    legend.position = "none" 
  )

Fig3A 


# B ----
# Plot QQ plot 
# Calculate expected log(p) 
expected_logP <- -log10(ppoints(length(gwas_DS$p.value)))  

# Sort the observed p-values
observed_logP <- -log10(sort(gwas_DS$p.value))  

Fig3B <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "#b6d9f6", linetype = "solid", size = 1) +  
  geom_point(aes(x = expected_logP, y = observed_logP), alpha = 0.6, color = "#648fb0") +  
  labs(
    title = "",
    x = expression(Expected ~ log[10] ~ "(" * italic(p) * ")"),  
    y = expression(Observed ~ log[10] ~ "(" * italic(p) * ")")
  ) +
  scale_x_continuous(breaks = seq(1, 5, by = 1)) + 
  scale_y_continuous(breaks = seq(1, 5, by = 1)) +  
  theme_cowplot()


# C ----
# Chromosome 7 
gwas_DS_chr7 <- subset(gwas_DS, chrom == 7)  # Subset for chromosome 7
gwas_DS_chr7 <- gwas_DS_chr7 %>%
  mutate(cum_pos = cum_pos - min(cum_pos) + 1)  # Adjust the positions to start at 1

# Identify the most significant SNP 
most_significant_snp <- gwas_DS_chr7 %>%
  filter(p.value == min(p.value))

# Add a new column to mark the most sign SNP
gwas_DS_chr7 <- gwas_DS_chr7 %>%
  mutate(is_significant = ifelse(p.value == min(p.value), "yes", "no"))

# Plot only CHR 7
Fig3C <- 
  ggplot(gwas_DS_chr7, aes(x = cum_pos / 1e6, y = logP)) +
  geom_hline(yintercept = genoSigLine, color = "grey40", linetype = "solid", linewidth = 0.6) +
  geom_hline(yintercept = suggestSigLine, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(color = is_significant, size = is_significant), alpha = 1) + 
  scale_color_manual(values = c("yes" = "#b6d9f6", "no" = "#648fb0")) +  
  scale_size_manual(values = c("yes" = 2, "no" = 1)) +  
  labs(
    title = "",
    x = "Chromosome 7 (Mb)", 
    y = expression(-log[10] (italic(p) * "-value"))
  ) +
  geom_text(
    data = most_significant_snp,
    aes(
      x = cum_pos / 1e6, 
      y = logP, 
      label = marker
    ),
    hjust = -0.1,  
    size = 3,  
    color = "black"  
  ) +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  ylim(0, 6) + 
  theme_cowplot() +
  theme(
    legend.position = "none" 
  )


Fig3C


# D ----

# Plot Days survived significant SNP 
Fig3D <- ggplot(snp_AX.691207752, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey40", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey40") +
  geom_jitter(shape = 19, color = "#b6d9f6", alpha = 1, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "GG", "1" = "AG", "2" = "AA")) +
  labs( 
    x = "SNP AX_691207752 Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")

Fig3D 


# Combine figure ----
Fig3 <- ggarrange(
  ggarrange(Fig3A, Fig3B, ncol = 2, labels = c("A", "B"), 
            widths = c(2, 1)), 
  ggarrange(Fig3C, Fig3D, ncol = 2, labels = c("C", "D"),
            widths = c(2, 1)), 
  nrow = 2
) 

Fig3

ggsave("figures/FigureS1-GWAS-survival.pdf", plot = Fig3, 
       width = 20, height = 10)


















#--- Format figure of Chr11 ----
# A ----
# Chromosome 11
gwas_DS_chr11 <- subset(gwas_DS, chrom == 11)
gwas_DS_chr11 <- gwas_DS_chr11 %>%
  mutate(cum_pos = cum_pos - min(cum_pos) + 1)

# Identify the most significant SNP 
most_significant_snp <- gwas_DS_chr11 %>%
  filter(p.value == min(p.value))

# Add a new column to mark the most sign SNP
gwas_DS_chr11 <- gwas_DS_chr11 %>%
  mutate(is_significant = ifelse(p.value == min(p.value), "yes", "no"))

# Plot only CHR 11
Fig9A <- 
  ggplot(gwas_DS_chr11, aes(x = cum_pos / 1e6, y = logP)) +
  geom_hline(yintercept = genoSigLine, color = "grey40", linetype = "solid", linewidth = 0.6) +
  geom_hline(yintercept = suggestSigLine, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(color = is_significant, size = is_significant), alpha = 1) + 
  scale_color_manual(values = c("yes" = "#b6d9f6", "no" = "#648fb0")) +  
  scale_size_manual(values = c("yes" = 2, "no" = 1)) +  
  labs(
    title = "",
    x = "Chromosome 11 (Mb)", 
    y = expression(-log[10] (italic(p) * "-value"))
  ) +
  geom_text(
    data = most_significant_snp,
    aes(
      x = cum_pos / 1e6, 
      y = logP, 
      label = marker
    ),
    hjust = -0.1,  
    size = 3,  
    color = "black"  
  ) +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  ylim(0, 6) + 
  theme_cowplot() +
  theme(
    legend.position = "none" 
  )


Fig9A


# B ----
Fig9B <- ggplot(snp_AX.691449917, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey40", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey40") +
  geom_jitter(shape = 19, color = "#b6d9f6", alpha = 1, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "GG", "1" = "AG", "2" = "AA")) +
  labs( 
    x = "SNP AX_691449917 Genotype",
    y = "GEBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")

Fig9B




# Combine ----
Fig9 <- ggarrange(
  ggarrange(Fig9A, Fig9B, ncol = 2, labels = c("A", "B"), 
            widths = c(2, 1))
) 

Fig9

ggsave("figures/FigureS3-GWAS-survivalCHR11.pdf", plot = FigS9, 
       width = 20, height = 5)

















#----------------------------
#   PATHOGEN LOAD ----
#----------------------------
#--- Model ----
# Preparation for gwas 
pre.slope <- pre.gwas(pheno.data=pheno, 
                      indiv="ID",
                      resp="slope_bd",
                      fixedf = "zoo",
                      cov = "age_start",
                      geno.data=X, 
                      map.data = map,
                      Q.method="K",
                      method="VanRaden",
                      maf=0.01, 
                      marker.callrate=0.8,
                      ind.callrate=0.20, 
                      heterozygosity=0.8, 
                      Fis=0.98, 
                      impute=TRUE)

#--- GWAS ----
gwasSlope <- gwas.asreml(pheno.data=pre.slope$pheno.data, 
                         resp="slope_bd",
                         fixedf = "zoo",
                         cov = "age_start",
                         gen="ID", 
                         Kinv=pre.slope$Kinv, 
                         Q=pre.slope$Q,
                         npc=7, 
                         geno.data=pre.slope$geno.data, 
                         map.data=pre.slope$map, 
                         pvalue.thr= 0.0001725626,
                         bonferroni=FALSE,
                         P3D=TRUE, 
                         #threads = 1
)

# A total of 3 markers identified based on a p-value threshold of 0.0001725626.
# The observed FDR for the set of significant markers is: 217.2448%.

# Marker effects 
head(gwasSlope$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasSlope$gwas.all, file = "output/ASReml-gwas_Slope.csv", sep = ",")

# List of significant SNPs at the pval set in model 
gwasSlope$gwas.sel 


#--- QQ plot ----
pdf("figures/ASReml/qq_Slope.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasSlope$gwas.all) +
  labs(title="Slope of Infection") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


#--- Manhattan plot ----
manhattan.plot(gwas.table = gwasSlope$gwas.all, 
               point.size = 1.2) + 
  geom_hline(yintercept =  genoSigLine, color = "darkgrey", linetype = "solid") + 
  geom_hline(yintercept =  suggestSigLine,  color = "darkgrey", linetype = "dashed") +
  ylim(0, 6)

pdf("figures/ASReml/manhattan_slope.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasSlope$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  genoSigLine , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  suggestSigLine , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Slope of Infection") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()






#--- USING EBVS ----
# EBV data is already loaded, and has zoo and age in file 
pre.slopeEBV <- pre.gwas(pheno.data=EBV, 
                         indiv="ID",
                         resp="EBV_slope",
                         geno.data=X, 
                         map.data = map,
                         Q.method="K",
                         method="VanRaden",
                         maf=0.01, 
                         marker.callrate=0.8,
                         ind.callrate=0.20, 
                         heterozygosity=0.8, 
                         Fis=0.98, 
                         impute=TRUE)

# Check values match, and what data cleaning has removed
dim(pre.slopeEBV$geno.data) # 895 37768

# Scree plot 
pre.slopeEBV$plot.scree

# Matches original data, with very low variation, =  7 looks okay... 
npc <-  7

# Prepare data for gwas = checking, cleaning, imputing 
# Set p val threshold 
# Threads need to be used if on the HPC / on local not needed
gwasSlopeEBV <- gwas.asreml(pheno.data=pre.slopeEBV$pheno.data, 
                            resp="EBV_slope",
                            gen="ID", 
                            Kinv=pre.slopeEBV$Kinv, 
                            Q=pre.slopeEBV$Q,
                            npc=7, 
                            geno.data=pre.slopeEBV$geno.data, 
                            map.data=pre.slopeEBV$map, 
                            pvalue.thr=        0.0001725626, #8.628128e-06,  
                            bonferroni=FALSE,
                            P3D=TRUE, 
                            #threads = 1
)

#A total of 24 markers identified based on a p-value threshold of 8.628128e-06.
#The observed FDR for the set of significant markers is: 1.3578%.

# A total of 60 markers identified based on a p-value threshold of 0.0001725626.
# The observed FDR for the set of significant markers is: 10.8622%.

# Marker effects 
head(gwasSlopeEBV$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasSlopeEBV$gwas.all, file = "output/ASReml-gwas_slope-EBV.csv", sep = ",")

# List of significant SNPs at the pval set in mod
gwasSlopeEBV$gwas.sel 

pdf("figures/ASReml/qq_slopeEBV.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasSlopeEBV$gwas.all) +
  labs(title="Slope EBV") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


pdf("figures/ASReml/manhattan_slopeEBV.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasSlopeEBV$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  -log10(7.640587e-06) , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  -log10(0.0001528117) , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Slope using EBV") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()


#--- Look at chromosomes separately ---- 
# Need to pull the data out of the gwas df, and into own df & rename cols 
gwas_slope_renamed <- gwasSlopeEBV$gwas.all %>%
  rename(CHR = chrom, 
         BP = pos, 
         P = p.value, 
         SNP = marker)


qqman::manhattan(subset(gwas_slope_renamed, CHR == 7), 
                 ylim = c(0, 7), 
                 main = "Chr 7", 
                 suggestiveline = suggestSigLine, 
                 genomewideline = genoSigLine)

# format plot 
qqman::manhattan(subset(gwas_slope_renamed, CHR == 7), 
                 main = "Slope of infection EBV CHR 7", ylim = c(0, 7),
                 col = c("#648fb0", "#b6d9f6"), 
                 suggestiveline = FALSE,
                 genomewideline =  FALSE, 
                 cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
abline(h = genoSigLine, col = "darkgrey", lty = 1, lwd=2)  
abline(h = suggestSigLine, col = "darkgrey", lty = 2, lwd=2)  


# Look into regions 
# Zoom into interesting region, and annotate SNPs
qqman::manhattan(subset(gwas_slope_renamed, CHR == 7), 
                 ylim=c(0,7), main = "Chr 7", 
                 xlim = c(240000000,270000000), 
                 suggestiveline = suggestSigLine,
                 genomewideline = genoSigLine, 
                 annotatePval = 0.001, annotateTop = FALSE)






# Fortmat figure  ----
# EBV GWAS for Slope
# A ----
# Format into DF to be able to use GGplot to plot for more control 

# Convert into DF not ASReml list
gwas_SL <- gwasSlopeEBV$gwas.all
# Transform p vals
gwas_SL$logP <- -log10(gwas_SL$p.value)  

# make chromosome sizes, to be able to plot as manhattan
chromosome_sizes <- gwas_SL %>%
  group_by(chrom) %>%
  summarize(chr_len = max(pos)) %>%
  mutate(chr_start = lag(cumsum(as.numeric(chr_len)), default = 0))

# Join chromsome sizes to df and calculate cumulative position
gwas_SL <- gwas_SL %>%
  left_join(chromosome_sizes, by = "chrom") %>%
  mutate(cum_pos = pos + chr_start)

# Calculate the midpoint for each chrom for labeling 
axis_df <- gwas_SL %>%
  group_by(chrom) %>%
  summarize(center = mean(cum_pos))

# Plot Manhattan 
Fig7A <- ggplot(gwas_SL, aes(x = cum_pos, y = logP, color = as.factor(chrom))) +
  geom_hline(yintercept = genoSigLine, color = "grey40", linetype = "solid", linewidth = 0.6) +
  geom_hline(yintercept = suggestSigLine, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_point(alpha = 1, size = 1) +  
  scale_color_manual(values = rep(c("#648fb0", "#b6d9f6"), 
                                  length.out = length(unique(gwas_DS$chrom)))) +  
  labs(title = "", 
       x = "Chromosome", 
       y = expression(-log[10] (italic(p) * "-value"))
  ) +
  ylim(0, 6.5) + 
  scale_x_continuous(labels = axis_df$chrom, breaks = axis_df$center, expand = c(0, 0)) +  
  theme_cowplot() +
  theme(
    legend.position = "none" 
  )

Fig7A


# B ----
# Plot QQ plot 
# Calculate expected log(p) 
expected_logP <- -log10(ppoints(length(gwas_SL$p.value)))  

# Sort the observed p-values
observed_logP <- -log10(sort(gwas_SL$p.value))  

Fig7B <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "#b6d9f6", linetype = "solid", size = 1) +  
  geom_point(aes(x = expected_logP, y = observed_logP), alpha = 0.6, color = "#648fb0") +  
  labs(
    title = "",
    x = expression(Expected ~ log[10] ~ "(" * italic(p) * ")"),  
    y = expression(Observed ~ log[10] ~ "(" * italic(p) * ")")
  ) +
  scale_x_continuous(breaks = seq(1, 6, by = 1)) + 
  scale_y_continuous(breaks = seq(1, 6, by = 1)) +  
  theme_cowplot()

Fig7B


# C ----
# Chromosome 7 
gwas_SL_chr7 <- subset(gwas_SL, chrom == 7)  # Subset for chromosome 7
gwas_SL_chr7 <- gwas_SL_chr7 %>%
  mutate(cum_pos = cum_pos - min(cum_pos) + 1)  # Adjust the positions to start at 1

# Identify the most significant SNP 
most_significant_snp <- gwas_SL_chr7 %>%
  filter(p.value == min(p.value))

# Add a new column to mark the most sign SNP
gwas_SL_chr7 <- gwas_SL_chr7 %>%
  mutate(is_significant = ifelse(p.value == min(p.value), "yes", "no"))

# Plot only CHR 7
Fig7C <- 
  ggplot(gwas_SL_chr7, aes(x = cum_pos / 1e6, y = logP)) +
  geom_hline(yintercept = genoSigLine, color = "grey40", linetype = "solid", linewidth = 0.6) +
  geom_hline(yintercept = suggestSigLine, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(color = is_significant, size = is_significant), alpha = 1) + 
  scale_color_manual(values = c("yes" = "#b6d9f6", "no" = "#648fb0")) +  
  scale_size_manual(values = c("yes" = 2, "no" = 1)) +  
  labs(
    title = "",
    x = "Chromosome 7 (Mb)", 
    y = expression(-log[10] (italic(p) * "-value"))
  ) +
  geom_text(
    data = most_significant_snp,
    aes(
      x = cum_pos / 1e6, 
      y = logP, 
      label = marker
    ),
    hjust = -0.1,  
    size = 3,  
    color = "black"  
  ) +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  ylim(0, 6.5) + 
  theme_cowplot() +
  theme(
    legend.position = "none" 
  )


Fig7C



# D ----

# Plot Days survived significant SNP 
Fig7D <- ggplot(snp_AX.691200909, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey40", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey40") +
  geom_jitter(shape = 19, color = "#b6d9f6", alpha = 1, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "TT", "1" = "CT", "2" = "CC")) +
  labs( 
    x = "SNP AX_691200909 Genotype",
    y = "GEBV Slope of Infection") +
  theme_cowplot() +
  theme(legend.position = "none")

Fig7D




# Combine ----
Fig7 <- ggarrange(
  ggarrange(Fig7A, Fig7B, ncol = 2, labels = c("A", "B"), 
            widths = c(2, 1)), 
  ggarrange(Fig7C, Fig7D, ncol = 2, labels = c("C", "D"),
            widths = c(2, 1)), 
  nrow = 2
) 

Fig7

ggsave("figures/FigureS2GWAS-pathLoad.pdf", plot = FigS7, 
       width = 20, height = 10)
















#----------------------------
#    BODY CONDITION ----
#----------------------------
#--- Model ----
# Preparation for gwas 
pre.bc <- pre.gwas(pheno.data=pheno, 
                   indiv="ID",
                   resp="bc_start",
                   fixedf = "zoo",
                   cov = "age_start",
                   geno.data=X, 
                   map.data = map,
                   Q.method="K",
                   method="VanRaden",
                   maf=0.01, 
                   marker.callrate=0.8,
                   ind.callrate=0.20, 
                   heterozygosity=0.8, 
                   Fis=0.98, 
                   impute=TRUE)

#--- GWAS ----
gwasBC <- gwas.asreml(pheno.data=pre.bc$pheno.data, 
                      resp="bc_start",
                      fixedf = "zoo",
                      cov = "age_start",
                      gen="ID", 
                      Kinv=pre.bc$Kinv, 
                      Q=pre.bc$Q,
                      npc=7, 
                      geno.data=pre.bc$geno.data, 
                      map.data=pre.bc$map, 
                      pvalue.thr= 0.0001725626,
                      bonferroni=FALSE,
                      P3D=TRUE)

# A total of 8 markers identified based on a p-value threshold of 0.0001725626.
# The observed FDR for the set of significant markers is: 81.4668%.

# Marker effects 
head(gwasBC$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasBC$gwas.all, file = "output/ASReml-gwas_BC.csv", sep = ",")

# List of significant SNPs at the pval set in model 
gwasBC$gwas.sel 


#--- QQ plot ----
pdf("figures/ASReml/qq_BC.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasBC$gwas.all) +
  labs(title="Body Condition") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


#--- Manhattan plot ----
manhattan.plot(gwas.table = gwasBC$gwas.all, 
               point.size = 1.2) + 
  geom_hline(yintercept =  genoSigLine, color = "darkgrey", linetype = "solid") + 
  geom_hline(yintercept =  suggestSigLine,  color = "darkgrey", linetype = "dashed") +
  ylim(0, 6)

pdf("figures/ASReml/manhattan_BC.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasBC$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  genoSigLine , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  suggestSigLine , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Body Condition") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()




#--- USING EBVS ----
# EBV data is already loaded, and has zoo and age in file 
pre.bcEBV <- pre.gwas(pheno.data=EBV, 
                      indiv="ID",
                      resp="EBV_bc",
                      geno.data=X, 
                      map.data = map,
                      Q.method="K",
                      method="VanRaden",
                      maf=0.01, 
                      marker.callrate=0.8,
                      ind.callrate=0.20, 
                      heterozygosity=0.8, 
                      Fis=0.98, 
                      impute=TRUE)

# Check values match, and what data cleaning has removed
dim(pre.bcEBV$geno.data) # 895 37768

# Scree plot 
pre.bcEBV$plot.scree

# Matches original data, with very low variation, =  7 looks okay... 
npc <-  7

# Prepare data for gwas = checking, cleaning, imputing 
# Set p val threshold 
# Threads need to be used if on the HPC / on local not needed
gwasbcEBV <- gwas.asreml(pheno.data=pre.bcEBV$pheno.data, 
                         resp="EBV_bc",
                         gen="ID", 
                         Kinv=pre.bcEBV$Kinv, 
                         Q=pre.bcEBV$Q,
                         npc=7, 
                         geno.data=pre.bcEBV$geno.data, 
                         map.data=pre.bcEBV$map, 
                         pvalue.thr=       0.0001725626, #8.628128e-06,   #
                         bonferroni=FALSE,
                         P3D=TRUE, 
                         #threads = 1
)

#A total of 0 markers identified based on a p-value threshold of 8.628128e-06.

# A total of 11 markers identified based on a p-value threshold of 0.0001725626.
# The observed FDR for the set of significant markers is: 59.2486%.

# Marker effects 
head(gwasbcEBV$gwas.all)

# Write out P vals of all SNPS, their effect size, and % of explained variance 
write.table(gwasbcEBV$gwas.all, file = "output/ASReml-gwas_bc-EBV.csv", sep = ",")

# List of significant SNPs at the pval set in mod
gwasbcEBV$gwas.sel 

pdf("figures/ASReml/qq_bcEBV.pdf", width = 6, height = 5)
qq.plot(gwas.table=gwasbcEBV$gwas.all) +
  labs(title="Body condition EBV") +
  geom_point(color = "#648fb0", shape = 19) + 
  geom_abline(intercept = 0, slope = 1, color = "#b6d9f6")  
dev.off()


pdf("figures/ASReml/manhattan_bcEBV.pdf", width = 14, height = 5)
manhattan.plot(gwas.table = gwasbcEBV$gwas.all, 
               point.size = 1.2) + 
  geom_point(shape = 16, size = 1.2) +
  geom_hline(yintercept =  -log10(7.640587e-06) , color = "grey20", linetype = "solid") + 
  geom_hline(yintercept =  -log10(0.0001528117) , color = "grey20", linetype = "dashed") +
  ylim(0, 6) +
  labs(title="Body condition using EBV") +
  theme_cowplot() +
  theme(legend.position="none")+
  scale_colour_manual(values = rep(c("#648fb0", "#b6d9f6"), 6))
dev.off()



#--------------------------------------------------------
#    SIGNIFICANT SNPS ----
#--------------------------------------------------------
# Load genotype data in 0/1/2 format converted in plink
geno012 <- read.table("data/significant-SNPs.raw", header = TRUE)
# Rename ID col 
colnames(geno012)[2] <- "ID"


# Merge geno012 with EBVs based on the common AnimalID
geno012 <- merge(geno012, EBVs[, c("ID", "EBV_days_survived", "EBV_slope", "EBV_maxLoad")], 
                 by = "ID")

# Convert to long formar
geno012_long <- na.omit(melt(
  geno012,
  id.vars = c("ID", "EBV_days_survived", "EBV_slope", "EBV_maxLoad"), 
  measure.vars = grep("^AX\\.", names(geno012), value = TRUE),         
  variable.name = "SNP",                                              
  value.name = "Genotype"                                              
))



#-----------------------------
#   SURVIVAL ----
#-----------------------------
# Plot all SNPs against days survived 
ggplot(geno012_long, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  geom_boxplot(aes(fill = as.factor(Genotype))) +
  facet_wrap(~ SNP, scales = "free_x") +
  labs(
    x = "Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")

# Only want Sig SNPs for days survived
# Make list of sig snps 

Dsurv_snps_to_keep <- c("AX.691207752_G", "AX.691200698_T", "AX.691200909_T", "AX.691200836_C", 
                        "AX.691201489_C", "AX.691201632_C", "AX.691201709_A", "AX.691219601_G", 
                        "AX.691219874_A", "AX.691219652_G", "AX.691449917_G")

# Subset DF for only these 11 SNPs
sigSNPs_dsurv_All <- geno012_long %>%
  filter(SNP %in% Dsurv_snps_to_keep)

# Plot against days survived 
ggplot(sigSNPs_dsurv_All, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  geom_boxplot(aes(fill = as.factor(Genotype))) +
  facet_wrap(~ SNP, scales = "free_x") +
  labs(
    x = "Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")

ggplot(sigSNPs_dsurv_All, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey") +
  facet_wrap(~ SNP, scales = "free_x") +
  labs(
    x = "Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")+
  
  geom_jitter(shape = 1, color = "grey50", alpha = 0.7, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) 



# Subset the SNP of highest significance = AX-691207752
snp_AX.691207752 <- subset(geno012_long, SNP == "AX.691207752_G")

# Count number of animals with each genotype 
snp_AX.691207752 %>%
  group_by(Genotype) %>%   
  summarise(count = n())  
# 0   479
# 1   336
# 2    50

# Plot Days survived significant SNP 
sigSNP_dsurv <- ggplot(snp_AX.691207752, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey") +
  geom_jitter(shape = 1, color = "grey50", alpha = 0.7, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "GG", "1" = "AG", "2" = "AA")) +
  labs( 
    x = "SNP AX-691207752 Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")


ggsave("figures/GWAS/sigSNP_dSurv_AX-691207752.pdf", plot = sigSNP_dsurv, 
       width = 7, height = 5)

# SIG SNP on CHR 11 
# Subset the SNP of highest significance = AX.691449917
snp_AX.691449917 <- subset(geno012_long, SNP == "AX.691449917_G")

# Count number of animals with each genotype 
snp_AX.691449917 %>%
  group_by(Genotype) %>%   
  summarise(count = n())  
# 0   479
# 1   354
# 2    41

# Plot Days survived significant SNP 
ggplot(snp_AX.691449917, aes(x = as.factor(Genotype), y = EBV_days_survived)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey") +
  geom_jitter(shape = 1, color = "grey50", alpha = 0.7, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "GG", "1" = "AG", "2" = "AA")) +
  labs( 
    x = "SNP AX_691449917 Genotype",
    y = "EBV Days Survived") +
  theme_cowplot() +
  theme(legend.position = "none")


#-----------------------------
#   PATHOGEN LOAD ----
#-----------------------------
# Subset the SNP of highest significace = AX_691200909
snp_AX.691200909 <- subset(geno012_long, SNP == "AX.691200909_T")

# Count number of animals with each genotype 
snp_AX.691200909 %>%
  group_by(Genotype) %>%   
  summarise(count = n())  
# 0   552
# 1   294
# 2    49

# Plot Days survived significant SNP 
sigSNP_slope <- ggplot(snp_AX.691200909, aes(x = as.factor(Genotype), y = EBV_slope)) +
  stat_boxplot(geom ='errorbar', 
               colour = "grey", 
               width = 0.2) + 
  geom_boxplot(aes(fill = as.factor(Genotype)),
               outlier.shape = NA,
               width = 0.5, 
               color = "grey") +
  geom_jitter(shape = 1, color = "grey50", alpha = 0.7, width = 0.1, size =1.5) + 
  scale_fill_manual(values = c("#9db4c5", "#9db4c5", "#9db4c5")) + 
  scale_x_discrete(labels = c("0" = "TT", "1" = "CT", "2" = "CC")) +
  labs( 
    x = "SNP AX-691200909 Genotype",
    y = "EBV Slope of Infection") +
  theme_cowplot() +
  theme(legend.position = "none")

sigSNP_slope 

ggsave("figures/GWAS/sigSNP_slope_AX-691200909.pdf", plot = sigSNP_slope, 
       width = 7, height = 5)


#----------------------------------------------------------------
#    SELECTION ----
#----------------------------------------------------------------
#-----------------------------
#   ALPHAMATE BY ZOO ----
#-----------------------------
# Subset by each zoo, and make GRMs to run in AlphaMate
#--- MZ ----
# Load bed matrix 
bedMZ <- read.bed.matrix("data/selection/MZ") 
# Create GRM
MZ_GRM <-GRM(bedMZ) 

# save gaston GRM as .csv to look at in excel 
write.table(MZ_GRM, file = "output/MZ-grm.csv", sep = ",")

#--- TZ ----
# Load bed matrix 
bedTZ <- read.bed.matrix("data/selection/TZ") 
# Create GRM
TZ_GRM <-GRM(bedTZ) 

# save gaston GRM as .csv to look at in excel 
write.table(TZ_GRM, file = "output/TZ-grm.csv", sep = ",")

