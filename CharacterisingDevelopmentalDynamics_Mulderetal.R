##################################
##### What makes clocks tick #####
##################################
### Script by R H Mulder, 2024
### Manuscript: 'Characterising developmental dynamics of adult epigenetic clock sites'



#########################
### Set working directory
#########################
setwd("")   # set to local directory


#################
### Load packages
#################
library(openxlsx)
library(meffil)


#############
### Load data
#############
#load epidelta results (see also http://epidelta.mrcieu.ac.uk/)
load(file="/input/epidelta/epidelta_results_pvals_derived_annot_201013.R")


#load clock sites
frst_horv             <- read.csv("/input/clocks/AdditionalFile3.csv")
frst_hannum           <- read.xlsx("/input/clocks/Hannum 2013 Suppl Table 1.xlsx")
frst_weidner          <- read.xlsx("/input/clocks/13059_2013_3218_MOESM2_ESM.xlsx", 3)
frst_zhang            <- read.table("/input/clocks/fromgithub_en.coef", header=TRUE)

scnd_phenoage_levine  <- read.csv("/input/clocks/Levine_2018_aging-10-101414-s002.csv")
scnd_telomereclock    <- read.xlsx("/input/clocks/aging-11-102173-s003.xlsx", startRow=6)
load("/input/clocks/fromgithub_methylclock_coefDunedinPACE.rda")             # coefDunedinPACE   


#load correlations over time in epidelta table
cors <- readRDS("/input/epidelta/correlations_0to17_ALSPAC_M3_240625.RDS")  # (0 to 17)
load(file="/input/epidelta/epidelta_randsds_M1_cors_230910.RData")   #m1_cors  (random effects)


#load prenatal exposures
smok        <- read.table("/input/PACE_prenexps_EWAScat/27040690_maternal_smoking_in_pregnancy_sustained_maternal_smoking_in_pregnancy_effect_on_newborns_adjusted_for_cell_counts.tsv", sep="\t", header=T)
bmi         <- read.table("/input/PACE_prenexps_EWAScat/29016858_sharp-gc_maternal_body_mass_index_meta-analysis_adjusted_for_cell_composition.tsv", sep="\t", header=T)
overwob     <- read.table("/input/PACE_prenexps_EWAScat/29016858_sharp-gc_maternal_overweight_obesity_meta-analysis_adjusted_for_cell_composition.tsv", sep="\t", header=T)
hypertens   <- read.table("/input/PACE_prenexps_EWAScat/31230546_kazmi-n_hypertensive_disorders_of_pregnancy_meta_analysis_adjusted_model.tsv", sep="\t", header=T)
hemoglob    <- read.table("/input/PACE_prenexps_EWAScat/33331245_ronkainen-j_maternal_haemoglobin_levels_during_pregnancy_meta-analysis_at_birth.tsv", sep="\t", header=T)
poll_no2    <- read.table("/input/PACE_prenexps_EWAScat/27448387_prenatal_no2_exposure_adjusted_for_cell_counts.tsv", sep="\t", header=T)
poll_part   <- read.table("/input/PACE_prenexps_EWAScat/31148503_gruzieva-h._particulate_matter__2_5um__pm2_5__discovery_and_replication.tsv", sep="\t", header=T)
anx_gener   <- read.table("/input/PACE_prenexps_EWAScat/sammallahti-s_maternal_anxiety_during_pregnancy_meta-analysis,_using_450k-only_sites.tsv", sep="\t", header=T)
anx_pregrel <- read.table("/input/PACE_prenexps_EWAScat/sammallahti-s_maternal_pregnancy-related_anxiety_during_pregnancy_meta-analysis_using_450k-only_sites.tsv", sep="\t", header=T)


#load WBC info Luo et al
luo_periph  <- read.xlsx("/input/WBCs/Luo 2023 Suppl.xlsx", sheet=2, startRow = 2)



#############################
### Prep correlations 0 to 17
#############################
dim(cors)

cors.df <- as.data.frame(cors)
cors.df.t <- tFrame(cors.df)

cors.df.t[1:3,]

cors.df.t$cor_0to17_sign     <- ifelse(cors.df.t$p.value<1e-7, "yes", "no")
cors.df.t$cor_0to17_estimate <- cors.df.t$estimate
cors.df.t$cor_0to17_p        <- cors.df.t$p.value

table(cors.df.t$cor_0to17_sign, useNA="always")
prop.table(table(cors.df.t$cor_0to17_sign, useNA="always"))*100



############################################
### Merge 0 to 17 correlations with epidelta
############################################
#pre check
all(rownames(epidelta)%in%rownames(cors.df.t))
#TRUE

#merge
epidelta <- merge(epidelta, cors.df.t[,c("cor_0to17_estimate","cor_0to17_p","cor_0to17_sign")], by="row.names", all.x=TRUE)

rownames(epidelta) <- epidelta$Row.names
epidelta$Row.names <- NULL


#######################################################
### Prep random SD correlations and merge with epidelta
#######################################################
m1_cors_df  <- as.data.frame(m1_cors[1:length(m1_cors)], names_m1_cors)
m1_cors_df  <- tFrame(m1_cors_df)

epidelta2 <- merge(epidelta, m1_cors_df, by="row.names", all.x=TRUE)
rownames(epidelta2) <- epidelta2$Row.names

epidelta2$Row.names <- NULL

epidelta <- epidelta2   

epidelta$M1.SD.int.age.pears.sign        <- ifelse(epidelta$M1.SD.int.age.pears.p < 1e-07,"yes","no")
epidelta$M1.SD.int.age.pears.dir         <- ifelse(epidelta$M1.SD.int.age.pears.p < 1e-07 & epidelta$M1.SD.int.age.pears.cor>0,"pos",
                                                   ifelse(epidelta$M1.SD.int.age.pears.p < 1e-07 & epidelta$M1.SD.int.age.pears.cor<0,"neg",NA))


#######################################################
### Prep prenatal exposure data and merge with epidelta
#######################################################
summary(smok$p)
summary(bmi$p)
summary(overwob$p)
summary(hypertens$p)
summary(hemoglob$p)
summary(poll_no2$p)
summary(poll_part$p)
summary(anx_gener$p)
summary(anx_pregrel$p)

#recheck smallest p-values
min(smok$p)
min(hypertens$p)


envexp.list <- list(smok, bmi, overwob, hypertens, hemoglob, poll_no2, poll_part, anx_gener, anx_pregrel)
envexp.names <- c("smok", "bmi", "overwob", "hypertens", "hemoglob", "poll_no2", "poll_part", "anx_gener", "anx_pregrel")

for(i in 1:length(envexp.list)){
  
  envexp.df   <- envexp.list[[i]]
  envexp.name <- envexp.names[i]
  
  print(paste0(envexp.name, ": ", nrow(envexp.df) ))
  
}

envexps <- rbind(smok, bmi, overwob, hypertens, hemoglob, poll_no2, poll_part, anx_gener, anx_pregrel)


#combine into unique file
envexps.uniq           <- data.frame(matrix(ncol = 7, nrow = length(unique(envexps$cpg))))
colnames(envexps.uniq) <- c("cpg","pace_expos","pace_author","pace_pmid","pace_beta","pace_se","pace_p")
envexps.uniq$cpg       <- unique(envexps$cpg)

dim(envexps)
length(unique(envexps$cpg))


###Creat file with unique CpGs, for which multiple exposures are collapsed in the same cells
for(cg in unique(envexps$cpg)){
  
  envexps.sub <- envexps[envexps$cpg %in% cg,]
  
  envexps.uniq[envexps.uniq$cpg==cg,"pace_exps"]   <- paste(envexps.sub$exposure, collapse=" / ")
  envexps.uniq[envexps.uniq$cpg==cg,"pace_author"] <- paste(envexps.sub$author, collapse=" / ")
  envexps.uniq[envexps.uniq$cpg==cg,"pace_pmid"]   <- paste(envexps.sub$pmid, collapse=" / ")
  envexps.uniq[envexps.uniq$cpg==cg,"pace_beta"]   <- paste(envexps.sub$beta, collapse=" / ")
  envexps.uniq[envexps.uniq$cpg==cg,"pace_se"]     <- paste(envexps.sub$se, collapse=" / ")
  envexps.uniq[envexps.uniq$cpg==cg,"pace_p"]      <- paste(envexps.sub$p, collapse=" / ")
  
}

dim(envexps.uniq)

###Merge with epidelta and clock files
epidelta$cpg       <- rownames(epidelta)
epidelta           <- merge(epidelta, envexps.uniq, by="cpg", all.x=TRUE) 
rownames(epidelta) <- epidelta$cpg



####################################################
### Prep and add WBC info for supplementary analyses
####################################################
dim(luo_periph)

epidelta$in_luo <- ifelse(rownames(epidelta)%in%luo_periph$probeID,"yes","no")

table(epidelta$in_luo, useNA="always")

salas_periph = meffil:::get.cell.type.reference("blood gse167998")$beta
gervin_cord  = meffil:::get.cell.type.reference("combined cord blood")$beta

dim(salas_periph)
colnames(salas_periph)

dim(gervin_cord)
colnames(gervin_cord)

gervin_cord      <- as.data.frame(gervin_cord)

epidelta$in_salas       <- ifelse(rownames(epidelta)%in%rownames(salas_periph),"yes","no")
epidelta$in_gervin      <- ifelse(rownames(epidelta)%in%rownames(gervin_cord),"yes","no")

table(epidelta$in_salas, useNA="always")
table(epidelta$in_gervin, useNA="always")



#######################
### Reformat clock info
#######################
scnd_dunedinpace <- coefDunedinPACE
rm(coefDunedinPACE)

frst_horv$intercept <- frst_horv[frst_horv$CpGmarker=="(Intercept)","CoefficientTraining"]
frst_horv           <- frst_horv[2:354,c(1,2,24)]
colnames(frst_horv) <- c("cpg","coef","intercept")
frst_horv$clock     <- "frst_horv"

frst_hannum$intercept <- 0                              # dit voor allemaal nog checken
frst_hannum           <- frst_hannum[,c(1,6,7)]
colnames(frst_hannum) <- c("cpg","coef","intercept")
frst_hannum$clock     <- "frst_hannum"

frst_weidner$intercept <- 0
frst_weidner           <- frst_weidner[6:107,c(1,4,580)]
colnames(frst_weidner) <- c("cpg","coef","intercept")
frst_weidner$clock     <- "frst_weidner"
frst_weidner$coef      <- as.numeric(frst_weidner$coef)

frst_zhang$intercept <- 0
frst_zhang           <- frst_zhang[2:515,]
colnames(frst_zhang) <- c("cpg","coef","intercept")
frst_zhang$clock     <- "frst_zhang"


scnd_phenoage_levine$intercept <- scnd_phenoage_levine[scnd_phenoage_levine$CpG=="Intercept","Weight"]
scnd_phenoage_levine           <- scnd_phenoage_levine[2:514,c(1,6,10)]
colnames(scnd_phenoage_levine) <- c("cpg","coef","intercept")
scnd_phenoage_levine$clock     <- "scnd_phenoage_levine"

scnd_telomereclock$intercept <- scnd_telomereclock[scnd_telomereclock$Variable=="Intercept","Coefficient"]
scnd_telomereclock           <- scnd_telomereclock[2:141,c(1,2,23)]
colnames(scnd_telomereclock) <- c("cpg","coef","intercept")
scnd_telomereclock$clock     <- "scnd_telomereclock"

scnd_dunedinpace$intercept <- scnd_dunedinpace[scnd_dunedinpace$CpGmarker=="(Intercept)","CoefficientTraining"]
scnd_dunedinpace$intercept <- 0      #nog checken
scnd_dunedinpace           <- scnd_dunedinpace[2:174,c(1,2,4)]
colnames(scnd_dunedinpace) <- c("cpg","coef","intercept")
scnd_dunedinpace$clock     <- "scnd_dunedinpace"


frst_horv$coef_dir            <- ifelse(frst_horv$coef>0,"positive","negative")
frst_hannum$coef_dir          <- ifelse(frst_hannum$coef>0,"positive","negative")
frst_weidner$coef_dir         <- ifelse(frst_weidner$coef>0,"positive","negative")
frst_zhang$coef_dir           <- ifelse(frst_zhang$coef>0,"positive","negative")

scnd_phenoage_levine$coef_dir <- ifelse(scnd_phenoage_levine$coef>0,"positive","negative")
scnd_telomereclock$coef_dir   <- ifelse(scnd_telomereclock$coef<0,"positive","negative")   #note REVERSED, because you want longer telomeres
scnd_dunedinpace$coef_dir     <- ifelse(scnd_dunedinpace$coef>0,"positive","negative")



###################
### Add to EpiDelta
###################
epidelta$M2.nonlinear.any     <- ifelse(epidelta$M2.agemin6.sign=="yes" | epidelta$M2.agemin9.sign=="yes", "yes", "no")
epidelta$M2.neut.change6      <- ifelse(epidelta$M2.nonlinear.short=="pos_neut" | epidelta$M2.nonlinear.short=="neg_neut", "yes", "no")
epidelta$meqtl_bios_cistrans  <- ifelse(epidelta$meqtl_bios_cis=="yes" | epidelta$meqtl_bios_trans=="yes", "yes", "no")



#######################
### Merge with EpiDelta
#######################
frst_horv_epi            <- merge(frst_horv, epidelta, by.x="cpg", by.y="row.names")
frst_hannum_epi          <- merge(frst_hannum, epidelta, by.x="cpg", by.y="row.names")
frst_weidner_epi         <- merge(frst_weidner, epidelta, by.x="cpg", by.y="row.names")
frst_zhang_epi           <- merge(frst_zhang, epidelta, by.x="cpg", by.y="row.names")

scnd_phenoage_levine_epi <- merge(scnd_phenoage_levine, epidelta, by.x="cpg", by.y="row.names")
scnd_telomereclock_epi   <- merge(scnd_telomereclock, epidelta, by.x="cpg", by.y="row.names")
scnd_dunedinpace_epi     <- merge(scnd_dunedinpace, epidelta, by.x="cpg", by.y="row.names")



##########
### Checks
##########
dim(frst_horv)
dim(frst_hannum)
dim(frst_weidner)
dim(frst_zhang)

dim(scnd_phenoage_levine)
dim(scnd_telomereclock)
dim(scnd_dunedinpace)



##################
### Combine clocks
##################
frst_gens <- rbind(frst_horv, frst_hannum, frst_weidner, frst_zhang)  # , frst_wu
scnd_gens <- rbind(scnd_phenoage_levine, scnd_telomereclock, scnd_dunedinpace)
clocks    <- rbind(frst_gens, scnd_gens)

#check duplicated sites
clocks$dups <- duplicated(clocks$cpg)
table(clocks$dups)
dups <- clocks[clocks$dups==TRUE,"cpg"]
table(clocks[clocks$cpg %in% dups,"clock"])

#combine wiht epidelta info
frst_gens_epi <- rbind(frst_horv_epi, frst_hannum_epi, frst_weidner_epi, frst_zhang_epi)    # , frst_wu_epi
scnd_gens_epi <- rbind(scnd_phenoage_levine_epi, scnd_telomereclock_epi, scnd_dunedinpace_epi)

clocks_epi <- rbind(frst_gens_epi, scnd_gens_epi)

#list lcoks
frst_gens_epi_list <- list(frst_horv_epi, frst_hannum_epi, frst_weidner_epi, frst_zhang_epi)
scnd_gens_epi_list <- list(scnd_phenoage_levine_epi, scnd_telomereclock_epi, scnd_dunedinpace_epi)


#for enrichment
epidelta$clock_frst_gen <- ifelse(rownames(epidelta)%in%frst_gens_epi$cpg, "yes", "no")
epidelta$clock_frst_gen <- ifelse(rownames(epidelta)%in%scnd_gens_epi$cpg & epidelta$clock_frst_gen=="no", NA, epidelta$clock_frst_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks
 
epidelta$clock_scnd_gen <- ifelse(rownames(epidelta)%in%scnd_gens_epi$cpg, "yes", "no")
epidelta$clock_scnd_gen <- ifelse(rownames(epidelta)%in%frst_gens_epi$cpg & epidelta$clock_scnd_gen=="no", NA, epidelta$clock_scnd_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks

table(epidelta$clock_frst_gen, useNA="always")
table(epidelta$clock_scnd_gen, useNA="always")

#check number of sites
length(frst_gens_epi$cpg)
length(scnd_gens_epi$cpg)
length(unique(frst_gens_epi$cpg))
length(unique(scnd_gens_epi$cpg))

###select only uniques
clocks_epi$dups    <- duplicated(clocks_epi$cpg)
frst_gens_epi$dups <- duplicated(frst_gens_epi$cpg)
scnd_gens_epi$dups <- duplicated(scnd_gens_epi$cpg)

clocks_epi_uniq    <- clocks_epi[clocks_epi$dups==FALSE,]
frst_gens_epi_uniq <- frst_gens_epi[frst_gens_epi$dups==FALSE,]
scnd_gens_epi_uniq <- scnd_gens_epi[scnd_gens_epi$dups==FALSE,]

#check
nrow(frst_gens_epi_uniq)
nrow(scnd_gens_epi_uniq)



#######################################
### Change and direction of change - M1
#######################################
###epidelta
prop.table(table(epidelta$M1.age.dir))*100

###first generation
prop.table(table(frst_gens_epi_uniq$M1.age.dir))*100

for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M1.age.dir"]))*100)
}

fisher.test(epidelta$clock_frst_gen, epidelta$M1.age.sign)
fisher.test(epidelta$clock_frst_gen, epidelta$M1.age.sign)[[1]]

###second generation
prop.table(table(scnd_gens_epi_uniq$M1.age.dir))*100

for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M1.age.dir"]))*100)
}

fisher.test(epidelta$clock_scnd_gen, epidelta$M1.age.sign)
fisher.test(epidelta$clock_scnd_gen, epidelta$M1.age.sign)[[1]]


###Direction of change consistent
#first generation 
fisher.test(frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","M1.age.dir"], frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","coef_dir"])  #not taking uniq cpgs bc directions might not be consistent between clocks
fisher.test(frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","M1.age.dir"], frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","coef_dir"])[[1]]

#second generation
fisher.test(scnd_gens_epi[scnd_gens_epi$M1.age.sign=="yes","M1.age.dir"], scnd_gens_epi[scnd_gens_epi$M1.age.sign=="yes","coef_dir"])



#########################
### Nonlinear change - M2
#########################
###epidelta
prop.table(table(epidelta$M2.nonlinear.any))*100
prop.table(table(epidelta$M2.agemin9.sign))*100

prop.table(table(epidelta$M2.nonlinear.short))*100
prop.table(table(epidelta$M2.neut.change6))*100


###first generation - nonlinear which
prop.table(table(frst_gens_epi_uniq$M2.nonlinear.short))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.nonlinear.short"]))*100)
}


###first generation - slope change at 6, neutral after
prop.table(table(frst_gens_epi_uniq$M2.neut.change6))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.neut.change6"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M2.neut.change6)
fisher.test(epidelta$clock_frst_gen, epidelta$M2.neut.change6)[[1]]

###first generation - slope change at 9
prop.table(table(frst_gens_epi_uniq$M2.agemin9.sign))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin9.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M2.agemin9.sign)
fisher.test(epidelta$clock_frst_gen, epidelta$M2.agemin9.sign)[[1]]


###second generation - nonlinear any
prop.table(table(scnd_gens_epi_uniq$M2.nonlinear.any))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.nonlinear.any"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.nonlinear.any)
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.nonlinear.any)[[1]]

###second generation - nonlinear which
prop.table(table(scnd_gens_epi_uniq$M2.nonlinear.short))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.nonlinear.short"]))*100)
}

###second generation - slope change at 6, neutral after
prop.table(table(scnd_gens_epi_uniq$M2.neut.change6))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.neut.change6"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.neut.change6)
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.neut.change6)[[1]]

###second generation - slope change at 9
prop.table(table(scnd_gens_epi_uniq$M2.agemin9.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin9.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.agemin9.sign)



######################################################
### Inter-individual variation at birth - intercept M1
######################################################
###epidelta
prop.table(table(epidelta$M1.intercept.rand.sign))*100

###first generation
prop.table(table(frst_gens_epi_uniq$M1.intercept.rand.sign))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M1.intercept.rand.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M1.intercept.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M1.intercept.rand.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M1.intercept.rand.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M1.intercept.rand.sign)
fisher.test(epidelta$clock_scnd_gen, epidelta$M1.intercept.rand.sign)[[1]]



#######################################################
### Inter-individual variation in slope from birth - M2
#######################################################
###epidelta
prop.table(table(epidelta$M2.age.rand.sign))*100

###first generation
prop.table(table(frst_gens_epi_uniq$M2.age.rand.sign))*100

for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.age.rand.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M2.age.rand.sign)
fisher.test(epidelta$clock_frst_gen, epidelta$M2.age.rand.sign)[[1]]

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.age.rand.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.age.rand.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.age.rand.sign)



###################################################
### Inter-individual variation in slope from 6 - M2
###################################################
###epidelta
prop.table(table(epidelta$M2.agemin6.rand.sign))*100

###first generation
prop.table(table(frst_gens_epi_uniq$M2.agemin6.rand.sign))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin6.rand.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M2.agemin6.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin6.rand.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin6.rand.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.agemin6.rand.sign)



###################################################
### Inter-individual variation in slope from 9 - M2
###################################################
###epidelta
prop.table(table(epidelta$M2.agemin9.rand.sign))*100

###first generation
prop.table(table(frst_gens_epi_uniq$M2.agemin9.rand.sign))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin9.rand.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M2.agemin9.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin9.rand.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M2.agemin9.rand.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M2.agemin9.rand.sign)



########################
### Correlations 0 to 17
########################
###epidelta
prop.table(table(epidelta$cor_0to17_sign))*100

###first generation
prop.table(table(frst_gens_epi_uniq$cor_0to17_sign))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"cor_0to17_sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$cor_0to17_sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$cor_0to17_sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"cor_0to17_sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$cor_0to17_sign)
fisher.test(epidelta$clock_scnd_gen, epidelta$cor_0to17_sign)[[1]]



################################################
### Correlations random SDs intercept and change
################################################
###epidelta
prop.table(table(epidelta$M1.SD.int.age.pears.sign))*100


###first generation
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"M1.SD.int.age.pears.sign"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$M1.SD.int.age.pears.sign)


###second generation
prop.table(table(scnd_gens_epi_uniq$M1.SD.int.age.pears.sign))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"M1.SD.int.age.pears.sign"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$M1.SD.int.age.pears.sign)





##########
### meQTLs 
##########
##### Main analyses: ALSPAC meQTLS in cordblood
###epidelta
prop.table(table(epidelta$meqtl_als_cord))*100


###first generation
prop.table(table(frst_gens_epi_uniq$meqtl_als_cord))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"meqtl_als_cord"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$meqtl_als_cord)

###second generation
prop.table(table(scnd_gens_epi_uniq$meqtl_als_cord))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"meqtl_als_cord"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$meqtl_als_cord)
fisher.test(epidelta$clock_scnd_gen, epidelta$meqtl_als_cord)[[1]]


### Sensitivity analyses: same in GoDMC
###epidelta
prop.table(table(epidelta$meqtl_bios_cistrans))*100


###first generation
prop.table(table(frst_gens_epi_uniq$meqtl_bios_cistrans))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"meqtl_bios_cistrans"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$meqtl_bios_cistrans)
fisher.test(epidelta$clock_frst_gen, epidelta$meqtl_bios_cistrans)[[1]]


###second generation
prop.table(table(scnd_gens_epi_uniq$meqtl_bios_cistrans))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"meqtl_bios_cistrans"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$meqtl_bios_cistrans)
fisher.test(epidelta$clock_scnd_gen, epidelta$meqtl_bios_cistrans)[[1]]



######################
### Prenatal exposures
######################
###epidelta
prop.table(table(epidelta$pace_exps_any))*100


###first generation
prop.table(table(frst_gens_epi_uniq$pace_exps_any))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"pace_exps_any"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$pace_exps_any)


###second generation
prop.table(table(scnd_gens_epi_uniq$pace_exps_any))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"pace_exps_any"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$pace_exps_any)


##which exposures
sort(table(frst_gens_epi_uniq$pace_exps, useNA="always"), decreasing=TRUE)
sort(table(scnd_gens_epi_uniq$pace_exps, useNA="always"), decreasing=TRUE)



#####################
### Figure: Bar plots
#####################
library(ggplot2)
library(wesanderson)
library(ggsignif)


comparisons <- c("M1.age.sign","M2.agemin9.sign","M2.intercept.rand.sign","M2.age.rand.sign","M2.agemin6.rand.sign","M2.agemin9.rand.sign","cor_0to17_sign","M1.SD.int.age.pears.sign","meqtl_als_cord","pace_exps_any")   
titles      <- c("change","nonlinear at 9y","variation at birth","variation from birth","variation from 6y","variation from 9y","birth and 17y","birth and change","genetic associations","prenatal environment") 

for(i in 1:length(comparisons)){
  
  comparison <- comparisons[i]
  title      <- titles[i]
  
  print(comparison)
  print(title)
  
  barplot.df           <- data.frame(matrix(NA, 3, 5))
  colnames(barplot.df) <- c("sites","prop.yes","prop.yes.min","prop.yes.max","sign")
  barplot.df$sites     <- c("all 450K sites", "first-generation clock sites", "second/third-generation clock sites")
  
  frsts.prop <- c()
  for(i in frst_gens_epi_list){
    frsts.prop <- c(frsts.prop, (prop.table(table(i[,comparison]))*100)["yes"][[1]])
  }
  
  scnds.prop <- c()
  for(i in scnd_gens_epi_list){
    scnds.prop <- c(scnds.prop, (prop.table(table(i[,comparison]))*100)["yes"][[1]])
  }
  
  barplot.df[barplot.df$sites=="all 450K sites","prop.yes"]                      <- (prop.table(table(epidelta[,comparison]))*100)[2][[1]]
  barplot.df[barplot.df$sites=="first-generation clock sites","prop.yes"]        <- (prop.table(table(epidelta[epidelta$clock_frst_gen=="yes",comparison]))*100)[2][[1]]
  barplot.df[barplot.df$sites=="second/third-generation clock sites","prop.yes"] <- (prop.table(table(epidelta[epidelta$clock_scnd_gen=="yes",comparison]))*100)[2][[1]]
  
  barplot.df[barplot.df$sites=="all 450K sites","prop.yes.min"]                      <- NA
  barplot.df[barplot.df$sites=="first-generation clock sites","prop.yes.min"]        <- min(frsts.prop, na.rm=TRUE)
  barplot.df[barplot.df$sites=="second/third-generation clock sites","prop.yes.min"] <- min(scnds.prop, na.rm=TRUE)
  
  barplot.df[barplot.df$sites=="all 450K sites","prop.yes.max"]                      <- NA
  barplot.df[barplot.df$sites=="first-generation clock sites","prop.yes.max"]        <- max(frsts.prop, na.rm=TRUE)
  barplot.df[barplot.df$sites=="second/third-generation clock sites","prop.yes.max"] <- max(scnds.prop, na.rm=TRUE)
  
  barplot.df[barplot.df$sites=="all 450K sites","sign"]                      <- ""
  barplot.df[barplot.df$sites=="first-generation clock sites","sign"]        <- ifelse(fisher.test(epidelta[,comparison], epidelta[,"clock_frst_gen"])[[1]]<0.05,"*","")
  barplot.df[barplot.df$sites=="second/third-generation clock sites","sign"] <- ifelse(fisher.test(epidelta[,comparison], epidelta[,"clock_scnd_gen"])[[1]]<0.05,"*","")
  
  barplot.df$label.loc <- barplot.df$prop.yes.max+1
  
  savename <- paste0("/output/barplots/barplot_",gsub("\\.","_",comparison),"_240730.png")
  
  plotthis <- ggplot(barplot.df, aes(x=sites, y=prop.yes, fill=sites)) +
    scale_fill_manual(values = rev(wes_palette("Zissou1", n=3)) ) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_errorbar(aes(ymin=prop.yes.min, ymax=prop.yes.max), width=.2, size=3) +
    geom_text(aes(label = sign, x = sites, y = label.loc), position = position_dodge(width = 0.8), vjust = -0.6, size=25, fontface="bold")  + 
    ylim(-5, 105) + 
    ggtitle(title) +
    xlab("") + 
    ylab("% of sites") +
    theme(plot.title = element_text(size=100),
          axis.title.y = element_text(size=90),
          axis.text.y = element_text(size=70),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    theme(aspect.ratio=1)
  #dev.off()
  ggsave(savename,plot=plotthis, height=50, width=50, units="cm", dpi=300)
  
}



#######################################
### Supplementary: Directions per clock
#######################################
###First generation
table(frst_horv_epi$M1.age.dir, frst_horv_epi$coef_dir)
fisher.test(frst_horv_epi[frst_horv_epi$M1.age.sign=="yes","M1.age.dir"], frst_horv_epi[frst_horv_epi$M1.age.sign=="yes","coef_dir"])

table(frst_hannum_epi$M1.age.dir, frst_hannum_epi$coef_dir)
fisher.test(frst_hannum_epi[frst_hannum_epi$M1.age.sign=="yes","M1.age.dir"], frst_hannum_epi[frst_hannum_epi$M1.age.sign=="yes","coef_dir"])

table(frst_weidner_epi$M1.age.dir, frst_weidner_epi$coef_dir)
fisher.test(frst_weidner_epi[frst_weidner_epi$M1.age.sign=="yes","M1.age.dir"], frst_weidner_epi[frst_weidner_epi$M1.age.sign=="yes","coef_dir"])
fisher.test(frst_weidner_epi[frst_weidner_epi$M1.age.sign=="yes","M1.age.dir"], frst_weidner_epi[frst_weidner_epi$M1.age.sign=="yes","coef_dir"])[[1]]

table(frst_zhang_epi$M1.age.dir, frst_zhang_epi$coef_dir)
fisher.test(frst_zhang_epi[frst_zhang_epi$M1.age.sign=="yes","M1.age.dir"], frst_zhang_epi[frst_zhang_epi$M1.age.sign=="yes","coef_dir"])
fisher.test(frst_zhang_epi[frst_zhang_epi$M1.age.sign=="yes","M1.age.dir"], frst_zhang_epi[frst_zhang_epi$M1.age.sign=="yes","coef_dir"])[[1]]

table(frst_gens_epi$M1.age.dir, frst_gens_epi$coef_dir)
fisher.test(frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","M1.age.dir"], frst_gens_epi[frst_gens_epi$M1.age.sign=="yes","coef_dir"])


###Second generation
table(scnd_phenoage_levine_epi$M1.age.dir, scnd_phenoage_levine_epi$coef_dir)
fisher.test(scnd_phenoage_levine_epi[scnd_phenoage_levine_epi$M1.age.sign=="yes","M1.age.dir"], scnd_phenoage_levine_epi[scnd_phenoage_levine_epi$M1.age.sign=="yes","coef_dir"])

table(scnd_telomereclock_epi$M1.age.dir, scnd_telomereclock_epi$coef_dir)
fisher.test(scnd_telomereclock_epi[scnd_telomereclock_epi$M1.age.sign=="yes","M1.age.dir"], scnd_telomereclock_epi[scnd_telomereclock_epi$M1.age.sign=="yes","coef_dir"])

table(scnd_dunedinpace_epi$M1.age.dir, scnd_dunedinpace_epi$coef_dir)
fisher.test(scnd_dunedinpace_epi[scnd_dunedinpace_epi$M1.age.sign=="yes","M1.age.dir"], scnd_dunedinpace_epi[scnd_dunedinpace_epi$M1.age.sign=="yes","coef_dir"])



########################################
### Supplementary: Enrichments per clock
########################################
library(stringr)

epidelta$horv          <- ifelse(rownames(epidelta) %in% frst_horv$cpg, "yes", "no")
epidelta$hannum        <- ifelse(rownames(epidelta) %in% frst_hannum$cpg, "yes", "no")
epidelta$weidner       <- ifelse(rownames(epidelta) %in% frst_weidner$cpg, "yes", "no")
epidelta$zhang         <- ifelse(rownames(epidelta) %in% frst_zhang$cpg, "yes", "no")
epidelta$phenoage      <- ifelse(rownames(epidelta) %in% scnd_phenoage_levine$cpg, "yes", "no")
epidelta$telomereclock <- ifelse(rownames(epidelta) %in% scnd_telomereclock$cpg, "yes", "no")
epidelta$dunedinpace   <- ifelse(rownames(epidelta) %in% scnd_dunedinpace$cpg, "yes", "no")


enrich.df <- as.data.frame(matrix(data=NA, nrow=11, ncol=22))

colnames(enrich.df) <- c("epidelta_PERC",
                         "frst_horv_epi_PERC","frst_horv_epi_OR","frst_horv_epi_PVAL", 
                         "frst_hannum_epi_PERC","frst_hannum_epi_OR","frst_hannum_epi_PVAL", 
                         "frst_weidner_epi_PERC","frst_weidner_epi_OR","frst_weidner_epi_PVAL", 
                         "frst_zhang_epi_PERC","frst_zhang_epi_OR","frst_zhang_epi_PVAL", 
                         "scnd_phenoage_levine_epi_PERC","scnd_phenoage_levine_epi_OR","scnd_phenoage_levine_epi_PVAL", 
                         "scnd_telomereclock_epi_PERC","scnd_telomereclock_epi_OR","scnd_telomereclock_epi_PVAL", 
                         "scnd_dunedinpace_epi_PERC","scnd_dunedinpace_epi_OR","scnd_dunedinpace_epi_PVAL")

rownames(enrich.df) <- c("M1.age.sign","M2.neut.change6","M2.agemin9.sign","M2.intercept.rand.sign","M2.age.rand.sign","M2.agemin6.rand.sign","M2.agemin9.rand.sign","cor_0to17_sign","M1.SD.int.age.pears.sign","meqtl_als_cord","pace_exps_any")

clocks.list  <- list(frst_horv_epi, frst_hannum_epi, frst_weidner_epi, frst_zhang_epi, scnd_phenoage_levine_epi, scnd_telomereclock_epi, scnd_dunedinpace_epi)
clocks.names <- c("horv", "hannum", "weidner", "zhang", "phenoage", "telomereclock", "dunedinpace")


for(i in rownames(enrich.df)){
  
  props.epidelta <- prop.table(table(epidelta[,i]))*100
  enrich.df[i,"epidelta_PERC"] <- sum(c(0, props.epidelta["yes"][[1]]), na.rm=TRUE)
  
  for(j in 1:length(clocks.list)){
    
    clock.which <- clocks.list[[j]]
    clock.name  <- clocks.names[j]
    cnames      <- colnames(enrich.df)[grep(clock.name, colnames(enrich.df))]
    cname_PERC  <- cnames[grep("PERC", cnames)]
    cname_OR    <- cnames[grep("OR", cnames)]
    cname_PVAL  <- cnames[grep("PVAL", cnames)]
    
    props.clock <- prop.table(table(clock.which[,i]))*100
    enrich.df[i,cname_PERC] <- sum(c(0,props.clock["yes"][[1]]), na.rm=TRUE)
    
    fish.clock <- fisher.test(epidelta[,i], epidelta[,clock.name])
    enrich.df[i, cname_OR]   <- fish.clock[[3]]
    enrich.df[i, cname_PVAL] <- fish.clock[[1]]
    
  }
  
}

write.table(enrich.df, file="/output/enrichment_perclock_240731.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

enrich.df





#############################################################
##### Supplementary analysis: GWAScat enrichment meQTLs #####
#############################################################

######################################################
### Load GWAScat results for ageing related phenotypes
######################################################
aging           <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0005422-withChildTraits.tsv")  # this one has on GWAScatalog EFO_0022597
skin_age        <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0005422.tsv")
skin_age_meas   <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0008006.tsv")
hippo_scler     <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0005678.tsv")
age_death       <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0005056.tsv")
longevity       <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0007796-withChildTraits.tsv")  # this one has on GWAScatalog EFO_0004300
healthspan      <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0009762.tsv")
age_meno        <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0004704.tsv")
parental_longev <- read_tsv("/input/GWAScat/gwas-association-downloaded_2024-06-11-EFO_0007796.tsv")


###########################
### Check number of studies
###########################
#number of publications
length(table(aging$PUBMEDID))
length(table(skin_age$PUBMEDID))
length(table(skin_age_meas$PUBMEDID))
length(table(hippo_scler$PUBMEDID))
length(table(age_death$PUBMEDID))
length(table(longevity$PUBMEDID))
length(table(healthspan$PUBMEDID))
length(table(age_meno$PUBMEDID))
length(table(parental_longev$PUBMEDID))


#number of studies
length(table(aging$"DISEASE/TRAIT"))
length(table(skin_age$"DISEASE/TRAIT"))
length(table(skin_age_meas$"DISEASE/TRAIT"))
length(table(hippo_scler$"DISEASE/TRAIT"))
length(table(age_death$"DISEASE/TRAIT"))
length(table(longevity$"DISEASE/TRAIT"))
length(table(healthspan$"DISEASE/TRAIT"))
length(table(age_meno$"DISEASE/TRAIT"))
length(table(parental_longev$"DISEASE/TRAIT"))



########
## Merge
########
gwascat <- do.call("rbind", list(aging, skin_age, skin_age_meas, hippo_scler, age_death, longevity, healthspan, age_meno, parental_longev))

dim(gwascat)

length(table(gwascat$PUBMEDID))
length(table(gwascat$"DISEASE/TRAIT"))



#############################
### Select epigenetic studies
#############################
epigen <- gwascat[grepl("methylation|epigenetic|Epigenetic|Biological age", gwascat$"DISEASE/TRAIT"),]

dim(epigen)

table(epigen$"DISEASE/TRAIT")

length(table(epigen$"DISEASE/TRAIT"))
length(table(epigen$PUBMEDID))
length(unique(epigen$SNPS))
length(unique(epigen$MAPPED_GENE))


gwascat_noepi <- gwascat[!grepl("methylation|epigenetic|Epigenetic|Biological age", gwascat$"DISEASE/TRAIT"),]

table(gwascat_noepi$"DISEASE/TRAIT")

length(table(gwascat_noepi$"DISEASE/TRAIT"))
length(table(gwascat_noepi$PUBMEDID))

length(unique(gwascat_noepi$SNPS))
length(unique(gwascat_noepi$MAPPED_GENE))



##############
### Enrichment
##############
######ARIES
frst_clock_meqtl_als_cord_snps <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),"meqtl_als_cord_snps"], collapse=","),",")))
scnd_clock_meqtl_als_cord_snps <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),"meqtl_als_cord_snps"], collapse=","),",")))

nonclock_meqtl_als_cord_snps   <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_frst_gen=="no" & epidelta$clock_scnd_gen=="no","meqtl_als_cord_snps"], collapse=","),",")))

#per clock
frst_horv_meqtl_als_cord_snps            <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_horv_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
frst_hannum_meqtl_als_cord_snps          <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_hannum_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
frst_weidner_meqtl_als_cord_snps         <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_weidner_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
frst_zhang_meqtl_als_cord_snps           <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_zhang_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
scnd_phenoage_levine_meqtl_als_cord_snps <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_phenoage_levine_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
scnd_telomereclock_meqtl_als_cord_snps   <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_telomereclock_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))
scnd_dunedinpace_meqtl_als_cord_snps     <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_dunedinpace_epi$cpg,"meqtl_als_cord_snps"], collapse=","),",")))


#take out "" (empty)
frst_clock_meqtl_als_cord_snps <- setdiff(frst_clock_meqtl_als_cord_snps, "")
scnd_clock_meqtl_als_cord_snps <- setdiff(scnd_clock_meqtl_als_cord_snps, "")
nonclock_meqtl_als_cord_snps   <- setdiff(nonclock_meqtl_als_cord_snps, "")

frst_horv_meqtl_als_cord_snps            <- setdiff(frst_horv_meqtl_als_cord_snps, "")
frst_hannum_meqtl_als_cord_snps          <- setdiff(frst_hannum_meqtl_als_cord_snps, "")
frst_weidner_meqtl_als_cord_snps         <- setdiff(frst_weidner_meqtl_als_cord_snps, "")
frst_zhang_meqtl_als_cord_snps           <- setdiff(frst_zhang_meqtl_als_cord_snps, "")
scnd_phenoage_levine_meqtl_als_cord_snps <- setdiff(scnd_phenoage_levine_meqtl_als_cord_snps, "")
scnd_telomereclock_als_cord_snps         <- setdiff(scnd_telomereclock_meqtl_als_cord_snps, "")
scnd_dunedinpace_als_cord_snps           <- setdiff(scnd_dunedinpace_meqtl_als_cord_snps, "")


length(frst_clock_meqtl_als_cord_snps)
length(scnd_clock_meqtl_als_cord_snps)
length(nonclock_meqtl_als_cord_snps)

length(frst_horv_meqtl_als_cord_snps)
length(frst_hannum_meqtl_als_cord_snps)
length(frst_weidner_meqtl_als_cord_snps)
length(frst_zhang_meqtl_als_cord_snps)
length(scnd_phenoage_levine_meqtl_als_cord_snps)
length(scnd_telomereclock_meqtl_als_cord_snps)
length(scnd_dunedinpace_meqtl_als_cord_snps)


#take out mqtls that are in first or second clock
nonclock_meqtl_als_cord_snps <- setdiff(nonclock_meqtl_als_cord_snps, c(frst_clock_meqtl_als_cord_snps, scnd_clock_meqtl_als_cord_snps))
length(nonclock_meqtl_als_cord_snps)


meqtl_als_cord_set <- as.data.frame(unique(c(frst_clock_meqtl_als_cord_snps, scnd_clock_meqtl_als_cord_snps, nonclock_meqtl_als_cord_snps)))
colnames(meqtl_als_cord_set) <- "meqtl"


meqtl_als_cord_set$clock_frst_gen <- ifelse(meqtl_als_cord_set$meqtl %in% frst_clock_meqtl_als_cord_snps, "yes", "no")
meqtl_als_cord_set$clock_frst_gen <- ifelse(meqtl_als_cord_set$meqtl %in% scnd_clock_meqtl_als_cord_snps & meqtl_als_cord_set$clock_frst_gen=="no", NA, meqtl_als_cord_set$clock_frst_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks

meqtl_als_cord_set$clock_scnd_gen <- ifelse(meqtl_als_cord_set$meqtl %in% scnd_clock_meqtl_als_cord_snps, "yes", "no")
meqtl_als_cord_set$clock_scnd_gen <- ifelse(meqtl_als_cord_set$meqtl %in% frst_clock_meqtl_als_cord_snps & meqtl_als_cord_set$clock_scnd_gen=="no", NA, meqtl_als_cord_set$clock_scnd_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks

meqtl_als_cord_set$in_age_gwas     <- ifelse(meqtl_als_cord_set$meqtl %in% gwascat_noepi$SNPS, "yes", "no")
meqtl_als_cord_set$in_epiage_gwas  <- ifelse(meqtl_als_cord_set$meqtl %in% epigen$SNPS, "yes", "no")
meqtl_als_cord_set$in_age_gwas_any <- ifelse(meqtl_als_cord_set$meqtl %in% c(gwascat_noepi$SNPS,epigen$SNPS), "yes", "no")


frst_horv_meqtl_als_cord            <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% frst_horv_meqtl_als_cord_snps, ]
frst_hannum_meqtl_als_cord          <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% frst_hannum_meqtl_als_cord_snps, ]
frst_weidner_meqtl_als_cord         <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% frst_weidner_meqtl_als_cord_snps, ]
frst_zhang_meqtl_als_cord           <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% frst_zhang_meqtl_als_cord_snps, ]
scnd_phenoage_levine_meqtl_als_cord <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% scnd_phenoage_levine_meqtl_als_cord_snps, ]
scnd_telomereclock_meqtl_als_cord   <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% scnd_telomereclock_meqtl_als_cord_snps, ]
scnd_dunedinpace_meqtl_als_cord     <- meqtl_als_cord_set[meqtl_als_cord_set$meqtl %in% scnd_dunedinpace_meqtl_als_cord_snps, ]


#list clocks
frst_gens_meqtl_als_cord_list <- list(frst_horv_meqtl_als_cord, frst_hannum_meqtl_als_cord, frst_weidner_meqtl_als_cord, frst_zhang_meqtl_als_cord)
scnd_gens_meqtl_als_cord_list <- list(scnd_phenoage_levine_meqtl_als_cord, scnd_telomereclock_meqtl_als_cord, scnd_dunedinpace_meqtl_als_cord)


#proportions
table(meqtl_als_cord_set$clock_frst_gen, useNA="always")
table(meqtl_als_cord_set$clock_scnd_gen, useNA="always")

table(meqtl_als_cord_set$in_age_gwas, useNA="always")
prop.table(table(meqtl_als_cord_set$in_age_gwas, useNA="always"))*100


table( meqtl_als_cord_set$in_epiage_gwas, useNA="always")
prop.table(table( meqtl_als_cord_set$in_epiage_gwas, useNA="always"))*100


###Enrichment ageing GWASs
table(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_age_gwas)
prop.table(table(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_age_gwas))*100
prop.table(table(meqtl_als_cord_set[meqtl_als_cord_set$clock_frst_gen=="yes", "in_age_gwas"]))*100
for(i in frst_gens_meqtl_als_cord_list){
  print(prop.table(table(i[,"in_age_gwas"]))*100)
}
fisher.test(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_age_gwas)



table(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_age_gwas)
prop.table(table(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_age_gwas))*100
prop.table(table(meqtl_als_cord_set[meqtl_als_cord_set$clock_scnd_gen=="yes", "in_age_gwas"]))*100
for(i in scnd_gens_meqtl_als_cord_list){
  print(prop.table(table(i[,"in_age_gwas"]))*100)
}
fisher.test(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_age_gwas)



###Enrichment epigenetic ageing GWASs
table(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_epiage_gwas)
prop.table(table(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_epiage_gwas))*100
prop.table(table(meqtl_als_cord_set[meqtl_als_cord_set$clock_frst_gen=="yes", "in_epiage_gwas"]))*100
for(i in frst_gens_meqtl_als_cord_list){
  print(prop.table(table(i[,"in_epiage_gwas"]))*100)
}
fisher.test(meqtl_als_cord_set$clock_frst_gen, meqtl_als_cord_set$in_epiage_gwas)



table(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_epiage_gwas)
prop.table(table(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_epiage_gwas))*100
prop.table(table(meqtl_als_cord_set[meqtl_als_cord_set$clock_scnd_gen=="yes", "in_epiage_gwas"]))*100
for(i in scnd_gens_meqtl_als_cord_list){
  print(prop.table(table(i[,"in_epiage_gwas"]))*100)
}
fisher.test(meqtl_als_cord_set$clock_scnd_gen, meqtl_als_cord_set$in_epiage_gwas)



######BIOS
frst_clock_meqtl_bios_snps <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),"meqtl_bios_snps"], collapse=","),",")))
scnd_clock_meqtl_bios_snps <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),"meqtl_bios_snps"], collapse=","),",")))

nonclock_meqtl_bios_snps   <- unique(unlist(strsplit(paste(epidelta[epidelta$clock_frst_gen=="no" & epidelta$clock_scnd_gen=="no","meqtl_bios_snps"], collapse=","),",")))

#per clock
frst_horv_meqtl_bios_snps            <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_horv_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
frst_hannum_meqtl_bios_snps          <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_hannum_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
frst_weidner_meqtl_bios_snps         <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_weidner_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
frst_zhang_meqtl_bios_snps           <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%frst_zhang_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
scnd_phenoage_levine_meqtl_bios_snps <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_phenoage_levine_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
scnd_telomereclock_meqtl_bios_snps   <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_telomereclock_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))
scnd_dunedinpace_meqtl_bios_snps     <- unique(unlist(strsplit(paste(epidelta[rownames(epidelta)%in%scnd_dunedinpace_epi$cpg,"meqtl_bios_snps"], collapse=","),",")))


#take out "" (empty)
frst_clock_meqtl_bios_snps <- setdiff(frst_clock_meqtl_bios_snps, "")
scnd_clock_meqtl_bios_snps <- setdiff(scnd_clock_meqtl_bios_snps, "")
nonclock_meqtl_bios_snps   <- setdiff(nonclock_meqtl_bios_snps, "")

frst_horv_meqtl_als_cord_snps            <- setdiff(frst_horv_meqtl_als_cord_snps, "")
frst_hannum_meqtl_als_cord_snps          <- setdiff(frst_hannum_meqtl_als_cord_snps, "")
frst_weidner_meqtl_als_cord_snps         <- setdiff(frst_weidner_meqtl_als_cord_snps, "")
frst_zhang_meqtl_als_cord_snps           <- setdiff(frst_zhang_meqtl_als_cord_snps, "")
scnd_phenoage_levine_meqtl_als_cord_snps <- setdiff(scnd_phenoage_levine_meqtl_als_cord_snps, "")
scnd_telomereclock_als_cord_snps         <- setdiff(scnd_telomereclock_meqtl_als_cord_snps, "")
scnd_dunedinpace_als_cord_snps           <- setdiff(scnd_dunedinpace_meqtl_als_cord_snps, "")


length(frst_clock_meqtl_bios_snps)
length(scnd_clock_meqtl_bios_snps)
length(nonclock_meqtl_bios_snps)

length(frst_horv_meqtl_bios_snps)
length(frst_hannum_meqtl_bios_snps)
length(frst_weidner_meqtl_bios_snps)
length(frst_zhang_meqtl_bios_snps)
length(scnd_phenoage_levine_meqtl_bios_snps)
length(scnd_telomereclock_meqtl_bios_snps)
length(scnd_dunedinpace_meqtl_bios_snps)


#take out mqtls that are in first or second clock
nonclock_meqtl_bios_snps <- setdiff(nonclock_meqtl_bios_snps, c(frst_clock_meqtl_bios_snps, scnd_clock_meqtl_bios_snps))
length(nonclock_meqtl_bios_snps)


meqtl_bios_set <- as.data.frame(unique(c(frst_clock_meqtl_bios_snps, scnd_clock_meqtl_bios_snps, nonclock_meqtl_bios_snps)))
colnames(meqtl_bios_set) <- "meqtl"


meqtl_bios_set$clock_frst_gen <- ifelse(meqtl_bios_set$meqtl %in% frst_clock_meqtl_bios_snps, "yes", "no")
meqtl_bios_set$clock_frst_gen <- ifelse(meqtl_bios_set$meqtl %in% scnd_clock_meqtl_bios_snps & meqtl_bios_set$clock_frst_gen=="no", NA, meqtl_bios_set$clock_frst_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks

meqtl_bios_set$clock_scnd_gen <- ifelse(meqtl_bios_set$meqtl %in% scnd_clock_meqtl_bios_snps, "yes", "no")
meqtl_bios_set$clock_scnd_gen <- ifelse(meqtl_bios_set$meqtl %in% frst_clock_meqtl_bios_snps & meqtl_bios_set$clock_scnd_gen=="no", NA, meqtl_bios_set$clock_scnd_gen)  # set to NA if in second-gen clocks, to compare first gen clocks with non-clocks

meqtl_bios_set$in_age_gwas     <- ifelse(meqtl_bios_set$meqtl %in% gwascat_noepi$SNPS, "yes", "no")
meqtl_bios_set$in_epiage_gwas  <- ifelse(meqtl_bios_set$meqtl %in% epigen$SNPS, "yes", "no")
meqtl_bios_set$in_age_gwas_any <- ifelse(meqtl_bios_set$meqtl %in% c(gwascat_noepi$SNPS,epigen$SNPS), "yes", "no")


frst_horv_meqtl_bios            <- meqtl_bios_set[meqtl_bios_set$meqtl %in% frst_horv_meqtl_bios_snps, ]
frst_hannum_meqtl_bios          <- meqtl_bios_set[meqtl_bios_set$meqtl %in% frst_hannum_meqtl_bios_snps, ]
frst_weidner_meqtl_bios         <- meqtl_bios_set[meqtl_bios_set$meqtl %in% frst_weidner_meqtl_bios_snps, ]
frst_zhang_meqtl_bios           <- meqtl_bios_set[meqtl_bios_set$meqtl %in% frst_zhang_meqtl_bios_snps, ]
scnd_phenoage_levine_meqtl_bios <- meqtl_bios_set[meqtl_bios_set$meqtl %in% scnd_phenoage_levine_meqtl_bios_snps, ]
scnd_telomereclock_meqtl_bios   <- meqtl_bios_set[meqtl_bios_set$meqtl %in% scnd_telomereclock_meqtl_bios_snps, ]
scnd_dunedinpace_meqtl_bios     <- meqtl_bios_set[meqtl_bios_set$meqtl %in% scnd_dunedinpace_meqtl_bios_snps, ]


#list clocks
frst_gens_meqtl_bios_list <- list(frst_horv_meqtl_bios, frst_hannum_meqtl_bios, frst_weidner_meqtl_bios, frst_zhang_meqtl_bios)
scnd_gens_meqtl_bios_list <- list(scnd_phenoage_levine_meqtl_bios, scnd_telomereclock_meqtl_bios, scnd_dunedinpace_meqtl_bios)


#proportions
table(meqtl_bios_set$clock_frst_gen, useNA="always")
table(meqtl_bios_set$clock_scnd_gen, useNA="always")

table(meqtl_bios_set$in_age_gwas, useNA="always")
prop.table(table(meqtl_bios_set$in_age_gwas, useNA="always"))*100
table( meqtl_bios_set$in_epiage_gwas, useNA="always")
prop.table(table( meqtl_bios_set$in_epiage_gwas, useNA="always"))*100


###Enrichment ageing GWASs
table(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_age_gwas)
prop.table(table(meqtl_bios_set[meqtl_bios_set$clock_frst_gen=="yes", "in_age_gwas"]))*100
prop.table(table(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_age_gwas))*100
for(i in frst_gens_meqtl_bios_list){
  print(prop.table(table(i[,"in_age_gwas"]))*100)
}
fisher.test(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_age_gwas)


table(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_age_gwas)
prop.table(table(meqtl_bios_set[meqtl_bios_set$clock_scnd_gen=="yes", "in_age_gwas"]))*100
prop.table(table(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_age_gwas))*100
for(i in scnd_gens_meqtl_bios_list){
  print(prop.table(table(i[,"in_age_gwas"]))*100)
}
fisher.test(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_age_gwas)



###Enrichment epigenetic ageing GWASs
table(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_epiage_gwas)
prop.table(table(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_epiage_gwas))*100
prop.table(table(meqtl_bios_set[meqtl_bios_set$clock_frst_gen=="yes", "in_epiage_gwas"]))*100
for(i in frst_gens_meqtl_bios_list){
  print(prop.table(table(i[,"in_epiage_gwas"]))*100)
}
fisher.test(meqtl_bios_set$clock_frst_gen, meqtl_bios_set$in_epiage_gwas)



table(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_epiage_gwas)
prop.table(table(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_epiage_gwas))*100
prop.table(table(meqtl_bios_set[meqtl_bios_set$clock_scnd_gen=="yes", "in_epiage_gwas"]))*100
for(i in scnd_gens_meqtl_bios_list){
  print(prop.table(table(i[,"in_epiage_gwas"]))*100)
}
fisher.test(meqtl_bios_set$clock_scnd_gen, meqtl_bios_set$in_epiage_gwas)



#########################################################
### Supplementary table: Enrichments per clock for meQTLs
#########################################################
library(stringr)

meqtl_als_cord_set$horv          <- ifelse(meqtl_als_cord_set$meqtl %in% frst_horv_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$hannum        <- ifelse(meqtl_als_cord_set$meqtl %in% frst_hannum_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$weidner       <- ifelse(meqtl_als_cord_set$meqtl %in% frst_weidner_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$zhang         <- ifelse(meqtl_als_cord_set$meqtl %in% frst_zhang_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$phenoage      <- ifelse(meqtl_als_cord_set$meqtl %in% scnd_phenoage_levine_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$telomereclock <- ifelse(meqtl_als_cord_set$meqtl %in% scnd_telomereclock_meqtl_als_cord$meqtl, "yes", "no")
meqtl_als_cord_set$dunedinpace   <- ifelse(meqtl_als_cord_set$meqtl %in% scnd_dunedinpace_meqtl_als_cord$meqtl, "yes", "no")


enrich.df <- as.data.frame(matrix(data=NA, nrow=2, ncol=22))

colnames(enrich.df) <- c("epidelta_PERC",
                         "frst_horv_epi_PERC","frst_horv_epi_OR","frst_horv_epi_PVAL", 
                         "frst_hannum_epi_PERC","frst_hannum_epi_OR","frst_hannum_epi_PVAL", 
                         "frst_weidner_epi_PERC","frst_weidner_epi_OR","frst_weidner_epi_PVAL", 
                         "frst_zhang_epi_PERC","frst_zhang_epi_OR","frst_zhang_epi_PVAL", 
                         "scnd_phenoage_levine_epi_PERC","scnd_phenoage_levine_epi_OR","scnd_phenoage_levine_epi_PVAL", 
                         "scnd_telomereclock_epi_PERC","scnd_telomereclock_epi_OR","scnd_telomereclock_epi_PVAL", 
                         "scnd_dunedinpace_epi_PERC","scnd_dunedinpace_epi_OR","scnd_dunedinpace_epi_PVAL")

rownames(enrich.df) <- c("in_age_gwas","in_epiage_gwas")

clocks.list  <- list(frst_horv_meqtl_als_cord, frst_hannum_meqtl_als_cord, frst_weidner_meqtl_als_cord, frst_zhang_meqtl_als_cord,
                     scnd_phenoage_levine_meqtl_als_cord, scnd_telomereclock_meqtl_als_cord, scnd_dunedinpace_meqtl_als_cord)
clocks.names <- c("horv", "hannum", "weidner", "zhang", "phenoage", "telomereclock", "dunedinpace")


for(i in rownames(enrich.df)){
  
  props.meqtl_als_cord_set <- prop.table(table(meqtl_als_cord_set[,i]))*100
  enrich.df[i,"epidelta_PERC"] <- sum(c(0, props.meqtl_als_cord_set["yes"][[1]]), na.rm=TRUE)
  
  for(j in 1:length(clocks.list)){
    
    clock.which <- clocks.list[[j]]
    clock.name  <- clocks.names[j]
    cnames      <- colnames(enrich.df)[grep(clock.name, colnames(enrich.df))]
    cname_PERC  <- cnames[grep("PERC", cnames)]
    cname_OR    <- cnames[grep("OR", cnames)]
    cname_PVAL  <- cnames[grep("PVAL", cnames)]
    
    props.clock <- prop.table(table(clock.which[,i]))*100
    enrich.df[i,cname_PERC] <- sum(c(0,props.clock["yes"][[1]]), na.rm=TRUE)
    
    fish.clock <- fisher.test(meqtl_als_cord_set[,i], meqtl_als_cord_set[,clock.name])
    enrich.df[i, cname_OR]   <- fish.clock[[3]]
    enrich.df[i, cname_PVAL] <- fish.clock[[1]]
    
  }
  
}

write.table(enrich.df, file="/output/enrichment_perclock_gwascatalog_240731.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

enrich.df




###################################################
##### Supplementary analysis: enrichment WBCs #####
###################################################

##########
### Gervin
##########
###epidelta
prop.table(table(epidelta$in_gervin))*100

###first generation
prop.table(table(frst_gens_epi_uniq$in_gervin))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"in_gervin"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$in_gervin)

###second generation
prop.table(table(scnd_gens_epi_uniq$in_gervin))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"in_gervin"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$in_gervin)



#########
### Salas
#########
###epidelta
prop.table(table(epidelta$in_salas))*100

###first generation
prop.table(table(frst_gens_epi_uniq$in_salas))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"in_salas"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$in_salas)

###second generation
prop.table(table(scnd_gens_epi_uniq$in_salas))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"in_salas"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$in_salas)



#########################
###Check enrichment Luo
#########################
###epidelta
prop.table(table(epidelta$in_luo))*100

###first generation
prop.table(table(frst_gens_epi_uniq$in_luo))*100
for(i in frst_gens_epi_list){
  print(prop.table(table(i[,"in_luo"]))*100)
}
fisher.test(epidelta$clock_frst_gen, epidelta$in_luo)

###second generation
prop.table(table(scnd_gens_epi_uniq$in_luo))*100
for(i in scnd_gens_epi_list){
  print(prop.table(table(i[,"in_luo"]))*100)
}
fisher.test(epidelta$clock_scnd_gen, epidelta$in_luo)




################################################################
##### Supplementary analysis: enrichment in matched probes #####
################################################################

########################################
### Match probes and check distributions
########################################
###load packages
library(MatchIt)


###Sample for 1st gen clock sites - random intercept
epidelta_frsts <- epidelta[!is.na(epidelta$clock_frst_gen),]

matched_frst_randints <- vector()

for(i in 1:100){
  epidelta_frsts            <- epidelta_frsts[!rownames(epidelta_frsts)%in%matched_frst_randints,]
  matched_frst_randint      <- matchit(clock_frst_gen ~ M1.intercept.rand.sd, data=epidelta_frsts,
                                       method="nearest", distance="euclidean")
  matched_frst_randints_new <- setdiff(rownames(match.data(matched_frst_randint)),unique(frst_gens_epi$cpg))
  matched_frst_randints     <- c(matched_frst_randints, matched_frst_randints_new)
  print(Sys.time())
}

length(matched_frst_randints)

#save
saveRDS(matched_frst_randints, file="/output/matched_frst_randints_240730.RDS")

summary(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),c("M1.intercept.rand.sd")])
summary(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),c("M1.intercept.rand.sd")])[[2]]
summary(epidelta[rownames(epidelta)%in%matched_frst_randints,c("M1.intercept.rand.sd")])    



###Sample for 2nd gen clock sites - random intercept
epidelta_scnds <- epidelta[!is.na(epidelta$clock_scnd_gen),]

matched_scnd_randints <- vector()

for(i in 1:100){
  epidelta_scnds            <- epidelta_scnds[!rownames(epidelta_scnds)%in%matched_scnd_randints,]
  matched_scnd_randint      <- matchit(clock_scnd_gen ~ M1.intercept.rand.sd, data=epidelta_scnds,
                                       method="nearest", distance="euclidean")
  matched_scnd_randints_new <- setdiff(rownames(match.data(matched_scnd_randint)),unique(scnd_gens_epi$cpg))
  matched_scnd_randints     <- c(matched_scnd_randints, matched_scnd_randints_new)
  print(Sys.time())
}

length(matched_scnd_randints)

#save
saveRDS(matched_scnd_randints, file="/output/matched_scnd_randints_240730.RDS")

summary(epidelta[epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),c("M1.intercept.rand.sd")])
summary(epidelta[rownames(epidelta)%in%matched_scnd_randints,c("M1.intercept.rand.sd")])


###Sample for 1st gen clock sites - age
epidelta_frsts <- epidelta[!is.na(epidelta$clock_frst_gen),]

matched_frst_fixedages <- vector()

for(i in 1:100){
  epidelta_frsts            <- epidelta_frsts[!rownames(epidelta_frsts)%in%matched_frst_fixedages,]
  matched_frst_randage      <- matchit(clock_frst_gen ~ M1.age.estimate + M1.age.se, data=epidelta_frsts,
                                       method="nearest", distance="euclidean")
  matched_frst_fixedages_new <- setdiff(rownames(match.data(matched_frst_randage)),unique(frst_gens_epi$cpg)) 
  matched_frst_fixedages     <- c(matched_frst_fixedages, matched_frst_fixedages_new)
  print(Sys.time())
}

length(matched_frst_fixedages)

#save
saveRDS(matched_frst_fixedages, file="/output/matched_frst_fixedages_240730.RDS")


summary(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),c("M1.age.estimate")])
summary(epidelta[rownames(epidelta)%in%matched_frst_fixedages,c("M1.age.estimate")])

summary(epidelta[epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),c("M1.age.se")])
summary(epidelta[rownames(epidelta)%in%matched_frst_fixedages,c("M1.age.se")])




###Sample for 2nd gen clock sites - age
epidelta_scnds <- epidelta[!is.na(epidelta$clock_scnd_gen),]

matched_scnd_fixedages <- vector()

for(i in 1:100){
  epidelta_scnds            <- epidelta_scnds[!rownames(epidelta_scnds)%in%matched_scnd_fixedages,]
  matched_scnd_randage      <- matchit(clock_scnd_gen ~ M1.age.estimate + M1.age.se, data=epidelta_scnds,
                                       method="nearest", distance="euclidean")
  matched_scnd_fixedages_new <- setdiff(rownames(match.data(matched_scnd_randage)),unique(scnd_gens_epi$cpg)) 
  matched_scnd_fixedages     <- c(matched_scnd_fixedages, matched_scnd_fixedages_new)
  print(Sys.time())
}

length(matched_scnd_fixedages)

#save
saveRDS(matched_scnd_fixedages, file="/output/matched_scnd_fixedages_240730.RDS")


summary(epidelta[epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),c("M1.age.estimate")])
summary(epidelta[rownames(epidelta)%in%matched_scnd_fixedages,c("M1.age.estimate")])

summary(epidelta[epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),c("M1.age.se")])
summary(epidelta[rownames(epidelta)%in%matched_scnd_fixedages,c("M1.age.se")])



###Add to epidelta
epidelta$clock_frst_gen_matched_ints <- ifelse(epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),"yes",
                                               ifelse(rownames(epidelta)%in%matched_frst_randints,"no",NA))
epidelta$clock_scnd_gen_matched_ints <- ifelse(epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),"yes",
                                               ifelse(rownames(epidelta)%in%matched_scnd_randints,"no",NA))

epidelta$clock_frst_gen_matched_ages <- ifelse(epidelta$clock_frst_gen=="yes" & !is.na(epidelta$clock_frst_gen),"yes",
                                               ifelse(rownames(epidelta)%in%matched_frst_fixedages,"no",NA))
epidelta$clock_scnd_gen_matched_ages <- ifelse(epidelta$clock_scnd_gen=="yes" & !is.na(epidelta$clock_scnd_gen),"yes",
                                               ifelse(rownames(epidelta)%in%matched_scnd_fixedages,"no",NA))


#check
table(epidelta$clock_frst_gen_matched_ints, useNA="always")
table(epidelta$clock_scnd_gen_matched_ints, useNA="always")

table(epidelta$clock_frst_gen_matched_ages, useNA="always")
table(epidelta$clock_scnd_gen_matched_ages, useNA="always")




#######################################
### Change and direction of change - M1 - matched on random intercept SD
#######################################
###first generation
prop.table(table(frst_gens_epi_uniq$M1.age.dir))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M1.age.dir"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.age.sign)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.age.sign)[[1]]

###second generation
prop.table(table(scnd_gens_epi_uniq$M1.age.dir))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M1.age.dir"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M1.age.sign)


###Increasing sites
#first gen
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.positive)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.positive)[[1]]

#second gen
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M1.positive)
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M1.positive)[[1]]


###Decreasing sites
#first gen
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.negative)

#second gen
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M1.negative)



#########################
### Nonlinear change - M2 - matched on random intercept SD
#########################
###first generation - slope change at 6, neutral after
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M2.neut.change6"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.neut.change6)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.neut.change6)[[1]]

###first generation - slope change at 9
prop.table(table(frst_gens_epi_uniq$M2.agemin9.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M2.agemin9.sign"]))*100

fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.agemin9.sign)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.agemin9.sign)[[1]]


###second generation - slope change at 6, neutral after
prop.table(table(scnd_gens_epi_uniq$M2.neut.change6))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M2.neut.change6"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M2.neut.change6)

###second generation - slope change at 9
prop.table(table(scnd_gens_epi_uniq$M2.agemin9.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M2.agemin9.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M2.agemin9.sign)



########################
### Correlations 0 to 17 - matched on random intercept SD
########################
###first generation 
prop.table(table(frst_gens_epi_uniq$cor_0to17_sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","cor_0to17_sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$cor_0to17_sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$cor_0to17_sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","cor_0to17_sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$cor_0to17_sign)



##################################################
### Correlations random intercept and random slope - matched on random intercept SD
##################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M1.SD.int.age.pears.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M1.SD.int.age.pears.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M1.SD.int.age.pears.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M1.SD.int.age.pears.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M1.SD.int.age.pears.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M1.SD.int.age.pears.sign)



#######################################################
### Inter-individual variation in slope from birth - M2 - matched on random intercept SD
#######################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.age.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M2.age.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.age.rand.sign)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.age.rand.sign)[[1]]

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.age.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M2.age.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M2.age.rand.sign)



###################################################
### Inter-individual variation in slope from 6 - M2 - matched on random intercept SD
###################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.agemin6.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M2.agemin6.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.agemin6.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin6.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M2.agemin6.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M2.agemin6.rand.sign)



###################################################
### Inter-individual variation in slope from 9 - M2 - matched on random intercept SD
###################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.agemin9.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","M2.agemin9.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$M2.agemin9.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin9.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","M2.agemin9.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$M2.agemin9.rand.sign)



##########
### meQTLs - matched on random intercept SD
##########
##### Main analyses: ALSPAC meQTLS in cordblood
###first generation 
prop.table(table(frst_gens_epi_uniq$meqtl_als_cord))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","meqtl_als_cord"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$meqtl_als_cord)

###second generation
prop.table(table(scnd_gens_epi_uniq$meqtl_als_cord))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","meqtl_als_cord"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$meqtl_als_cord)


####### Supplemental analyses: GoDMC
###first generation 
prop.table(table(frst_gens_epi_uniq$meqtl_bios_cistrans))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","meqtl_bios_cistrans"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$meqtl_bios_cistrans)
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$meqtl_bios_cistrans)[[1]]

###second generation
prop.table(table(scnd_gens_epi_uniq$meqtl_bios_cistrans))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","meqtl_bios_cistrans"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$meqtl_bios_cistrans)



######################
### Prenatal exposures - matched on random intercept SD
######################
###first generation 
prop.table(table(frst_gens_epi_uniq$pace_exps_any))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ints=="no","pace_exps_any"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ints, epidelta$pace_exps_any)

###second generation
prop.table(table(scnd_gens_epi_uniq$pace_exps_any))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ints=="no","pace_exps_any"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ints, epidelta$pace_exps_any)



#######################################
### Inter-individual variation at birth - matched on change
#######################################
###first generation
prop.table(table(frst_gens_epi_uniq$M1.intercept.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ages=="no","M1.intercept.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ages, epidelta$M1.intercept.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M1.intercept.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ages=="no","M1.intercept.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ages, epidelta$M1.intercept.rand.sign)
fisher.test(epidelta$clock_scnd_gen_matched_ages, epidelta$M1.intercept.rand.sign)[[1]]



#######################################################
### Inter-individual variation in slope from birth - M2 - matched on change
#######################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.age.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ages=="no","M2.age.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ages, epidelta$M2.age.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.age.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ages=="no","M2.age.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ages, epidelta$M2.age.rand.sign)



###################################################
### Inter-individual variation in slope from 6 - M2 - matched on change
###################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.agemin6.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ages=="no","M2.agemin6.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ages, epidelta$M2.agemin6.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin6.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ages=="no","M2.agemin6.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ages, epidelta$M2.agemin6.rand.sign)



###################################################
### Inter-individual variation in slope from 9 - M2 - matched on change
###################################################
###first generation 
prop.table(table(frst_gens_epi_uniq$M2.agemin9.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_frst_gen_matched_ages=="no","M2.agemin9.rand.sign"]))*100
fisher.test(epidelta$clock_frst_gen_matched_ages, epidelta$M2.agemin9.rand.sign)

###second generation
prop.table(table(scnd_gens_epi_uniq$M2.agemin9.rand.sign))*100    
prop.table(table(epidelta[epidelta$clock_scnd_gen_matched_ages=="no","M2.agemin9.rand.sign"]))*100
fisher.test(epidelta$clock_scnd_gen_matched_ages, epidelta$M2.agemin9.rand.sign)



