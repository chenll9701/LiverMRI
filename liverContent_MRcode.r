##risk factors: BMI, HC, WC, WHR, BFP, LBM, HDL, LDL, TC, TG, 
##HbA1c, FI, FG, 2hGlu, leptin, ferritin, iron, transferrin saturation
##SBP, DBP, T2D, NAFLD, CHD
options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
{###All these above are for risk factor data preparation
HbA1c <- as.data.frame(fread("C:/2345Downloads/MAGIC1000G_HbA1c_EUR.tsv.gz"))
FI <- as.data.frame(fread("C:/2345Downloads/MAGIC1000G_FI_EUR.tsv.gz"))
FG <- as.data.frame(fread("C:/2345Downloads/MAGIC1000G_FG_EUR.tsv.gz"))
Glu2h <- as.data.frame(fread("C:/2345Downloads/MAGIC1000G_2hGlu_EUR.tsv.gz"))
format_glycemicData <- function(mydata, snps = NULL, pheno) {
	if(!is.null(snps)) {mydata <- mydata[which(mydata$rsid%in%snps),]}
	if(is.null(snps)) {mydata <- mydata[which(mydata$effect_allele_frequency > 0.01&mydata$p_value<5e-8),]}
	mydata <- TwoSampleMR::format_data(
				mydata,
				type = "exposure",
				snp_col = "variant",
				beta_col = "beta",
				se_col = "standard_error",
				eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value",
				phenotype_col = pheno
				)
	if(is.null(snps)) {mydata <- TwoSampleMR::clump_data(mydata, clump_r2 = 0.01, clump_kb = 10000)}
	return(mydata)
	}

myHbA1c <- format_glycemicData(HbA1c, snps = NULL, pheno = 'HbA1c')
myFI <- format_glycemicData(FI, snps = NULL, pheno = 'FI')
myFG <- format_glycemicData(FG, snps = NULL, pheno = 'FG')
my2hGlu <- format_glycemicData(Glu2h, snps = NULL, pheno = '2hGlu')
fwrite(myFI, "C:/Users/86151/Documents/MyMaster/IV_for_MR/FI_IV.csv")
fwrite(myFG, "C:/Users/86151/Documents/MyMaster/IV_for_MR/FG_IV.csv")
fwrite(myHbA1c, "C:/Users/86151/Documents/MyMaster/IV_for_MR/HbA1c_IV.csv")
fwrite(my2hGlu, "C:/Users/86151/Documents/MyMaster/IV_for_MR/2hGlu_IV.csv")


##Liver volume, fat and iron data
ironData <- fread("liver_iron.tsv.gz", header=T)
fatData <- fread("liver_fat.tsv.gz", header=T)
volumeData <- fread("liver_volume.tsv.gz", header=T)
format_LiverData <- function(mydata, snps = NULL, pheno) {
	if(!is.null(snps)) {mydata <- mydata[which(mydata$rsid%in%snps),]}
	if(is.null(snps)) {mydata <- mydata[which(mydata$effect_allele_frequency > 0.01&mydata$p_value<5e-8),]}
	mydata <- TwoSampleMR::format_data(
				mydata,
				type = "exposure",
				snp_col = "variant_id",
				beta_col = "beta",
				se_col = "standard_error",
				eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value",
				phenotype_col = pheno
				)
	if(is.null(snps)) {mydata <- TwoSampleMR::clump_data(mydata, clump_r2 = 0.01, clump_kb = 10000)}
	return(mydata)
	}
myIron <- format_LiverData(ironData, snps = NULL, pheno = 'iron')
myFat <- format_LiverData(fatData, snps = NULL, pheno = 'fat')
myVolume <- format_LiverData(volumeData, snps = NULL, pheno = 'volume')
fwrite(myIron, "C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_iron_IV.csv")
fwrite(myFat, "C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_fat_IV.csv")
fwrite(myVolume, "C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_volume_IV.csv")

format_lipidData <- function(mydata, snps = NULL) {
	if(!is.null(snps)) {mydata <- mydata[which(mydata$rsid%in%snps),]}
	mydata$pval <- as.numeric(mydata$`P-value`)
	if(is.null(snps)) {mydata <- mydata[which(mydata$Freq.A1.1000G.EUR > 0.01&mydata$pval<5e-8),]}
	mydata <- TwoSampleMR::format_data(
				mydata,
				type = "exposure",
				snp_col = "rsid",
				beta_col = "beta",
				se_col = "se",
				eaf_col = "Freq.A1.1000G.EUR",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "pval",
				phenotype_col = "pheno"
				)
	if(is.null(snps)) {mydata <- TwoSampleMR::clump_data(mydata, clump_r2 = 0.01, clump_kb = 10000)}
	return(mydata)
	}
hdl <- as.data.frame(fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/lipid_exposure/jointGwasMc_HDL.txt.gz"), header=T)
hdl$pheno <- "HDL"
myHDL <- format_lipidData(hdl)
fwrite(myHDL, "C:/Users/86151/Documents/MyMaster/IV_for_MR/HDL_IV.csv")
ldl <- as.data.frame(fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/lipid_exposure/jointGwasMc_LDL.txt.gz"), header=T)
ldl$pheno <- "LDL"
myLDL <- format_lipidData(ldl)
fwrite(myLDL, "C:/Users/86151/Documents/MyMaster/IV_for_MR/LDL_IV.csv")
tc <- as.data.frame(fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/lipid_exposure/jointGwasMc_TC.txt.gz"), header=T)
tc$pheno <- "TC"
myTC <- format_lipidData(tc)
fwrite(myTC, "C:/Users/86151/Documents/MyMaster/IV_for_MR/TC_IV.csv")
tg <- as.data.frame(fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/lipid_exposure/jointGwasMc_TG.txt.gz"), header=T)
tg$pheno <- "TG"
myTG <- format_lipidData(tg)
fwrite(myTG, "C:/Users/86151/Documents/MyMaster/IV_for_MR/TG_IV.csv")

##smoking and drinking
drinking <- extract_instruments("ieu-b-73", r2 = 0.01, access_token = NULL) #drinking
fwrite(drinking, "C:/Users/86151/Documents/MyMaster/IV_for_MR/drinking_IV.csv")
smoking <- extract_instruments("ieu-b-142", r2 = 0.01, access_token = NULL) #smoking PMID:30643251
fwrite(smoking, "C:/Users/86151/Documents/MyMaster/IV_for_MR/smoking_IV.csv")
#body fat percentage
body_fat_percent <- extract_instruments("ebi-a-GCST003435", r2 = 0.01, access_token = NULL) #PMID:26833246
fwrite(body_fat_percent,"C:/Users/86151/Documents/MyMaster/IV_for_MR/body_fat_percent_IV.csv")
lean_body_mass <- extract_instruments("ebi-a-GCST004770", p1 = 1e-5, r2 = 0.01, access_token = NULL) #lean body mass #PMID:28743860
fwrite(lean_body_mass,"C:/Users/86151/Documents/MyMaster/IV_for_MR/lean_body_mass_1e-5_IV.csv")
##leptin （PMID：26833098)
leptin <- extract_instruments("ebi-a-GCST003367", p1 = 1e-5, r2 = 0.01, access_token = NULL) #  circulating leptin levels
fwrite(leptin,"C:/Users/86151/Documents/MyMaster/IV_for_MR/leptin_1e-5_IV.csv")
leptinadjBMI <- extract_instruments("ebi-a-GCST003368", r2 = 0.01, access_token = NULL) #  circulating leptin levels adjusted for BMI
fwrite(leptinadjBMI, "C:/Users/86151/Documents/MyMaster/IV_for_MR/leptinadjBMI_IV.csv")

##adiponectin
adiponectin <- extract_instruments("ieu-a-1", r2 = 0.01, access_token = NULL) #adiponectin PMID: 22479202
fwrite(adiponectin, "C:/Users/86151/Documents/MyMaster/IV_for_MR/adiponectin_IV.csv")

##GIS consortium (PMID；25352340）
iron <- extract_instruments("ieu-a-1049", r2 = 0.01, access_token = NULL) #iron
fwrite(iron, "C:/Users/86151/Documents/MyMaster/IV_for_MR/iron_IV.csv")
ferritin <- extract_instruments("ieu-a-1050", r2 = 0.01, access_token = NULL) #ferritin
fwrite(ferritin, "C:/Users/86151/Documents/MyMaster/IV_for_MR/ferritin_IV.csv")
transferrin <- extract_instruments("ieu-a-1052", r2 = 0.01, access_token = NULL) #Transferrin
fwrite(transferrin, "C:/Users/86151/Documents/MyMaster/IV_for_MR/transferrin_IV.csv")
transferrin_saturation <- extract_instruments("ieu-a-1051", r2 = 0.01, access_token = NULL) #Transferrin Saturation
fwrite(transferrin_saturation, "C:/Users/86151/Documents/MyMaster/IV_for_MR/transferrin_saturation_IV.csv")

##T2D
T2D <- extract_instruments("ieu-a-23", r2 = 0.01, access_token = NULL) #T2D PMID: 24509480
fwrite(T2D, "C:/Users/86151/Documents/MyMaster/IV_for_MR/T2D_IV.csv")
##CHD
CHD <- extract_instruments("ieu-a-9", r2 = 0.01, access_token = NULL) #CHD PMID: 23202125
fwrite(CHD, "C:/Users/86151/Documents/MyMaster/IV_for_MR/CHD_IV.csv")
###END
}



options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
mainMR_risk_factors <- function(exp_dat, out_dat, expName, outName, samplesize.outcome) {
	out_dat <- out_dat[which(out_dat$variant_id%in%exp_dat$SNP),]
	if(!is.null(out_dat)){
		out_dat <- TwoSampleMR::format_data(
					out_dat,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
		mydata <- harmonise_data(exp_dat, out_dat)
		mydata$samplesize.outcome <- samplesize.outcome
		if(dim(mydata)[1] >=1) {
			res <- mr(mydata, method_list = c("mr_wald_ratio","mr_ivw", "mr_weighted_median", "mr_egger_regression"))
			steiger <- directionality_test(mydata)
			res$exposure_r2 <- steiger[1,5]
			res$outcome_r2 <- steiger[1,6]
			res$direction <- steiger[1,7]
			res$PRESSO_beta <- NA
			res$PRESSO_pval <- NA
			res$IVW_Qpval <- NA
			res$egger_pleiotropy_pval <- NA
			if(dim(mydata)[1] > 3) {
				presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
										SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
										data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
				res$PRESSO_beta <- presso$`Main MR results`$`Causal Estimate`[2]
				res$PRESSO_pval <- presso$`Main MR results`$`P-value`[2]
				het <- mr_heterogeneity(mydata)
				pleio <- mr_pleiotropy_test(mydata)
				res$IVW_Qpval <- het$Q_pval[2]
				res$egger_pleiotropy_pval <- pleio$pval
				}
			res$expName <- expName
			res$outName <- outName
			}
		return(res)
	}
}

IV_dir <- "C:/Users/86151/Documents/MyMaster/IV_for_MR"
IV_list <- list.files(IV_dir)
#test <- read.csv(paste0(IV_dir,"/",IV_list[1]), header = T)
ironData <- fread("liver_iron.tsv.gz", header=T)
risk_factor_LiverIron <- {}
for (i in 9:length(IV_list)) {
	expName <- gsub(IV_list[i], pattern="_IV.csv", replacement = "")
	outName <- "Liver_iron"
	exp_dat <- read.csv(paste0(IV_dir,"/",IV_list[i]), header = T)
	myres <- mainMR_risk_factors(exp_dat, out_dat = ironData, expName, outName, samplesize.outcome = 32858)
	risk_factor_LiverIron <- rbind(myres, risk_factor_LiverIron)
}
fwrite(risk_factor_LiverIron, "risk_factor_LiverIron_results.csv")

IV_dir <- "C:/Users/86151/Documents/MyMaster/IV_for_MR"
IV_list <- list.files(IV_dir)
fatData <- fread("liver_fat.tsv.gz", header=T)
risk_factor_LiverFat <- {}
for (i in 1:length(IV_list)) {
	expName <- gsub(IV_list[i], pattern="_IV.csv", replacement = "")
	outName <- "Liver_fat"
	exp_dat <- read.csv(paste0(IV_dir,"/",IV_list[i]), header = T)
	myres <- mainMR_risk_factors(exp_dat, out_dat = fatData, expName, outName, samplesize.outcome = 32858)
	risk_factor_LiverFat <- rbind(myres, risk_factor_LiverFat)
}
fwrite(risk_factor_LiverFat, "risk_factor_LiverFat_results.csv")

IV_dir <- "C:/Users/86151/Documents/MyMaster/IV_for_MR"
IV_list <- list.files(IV_dir)
volumeData <- fread("liver_volume.tsv.gz", header=T)
risk_factor_LiverVolume <- {}
for (i in 1:length(IV_list)) {
	expName <- gsub(IV_list[i], pattern="_IV.csv", replacement = "")
	outName <- "Liver_volume"
	exp_dat <- read.csv(paste0(IV_dir,"/",IV_list[i]), header = T)
	myres <- mainMR_risk_factors(exp_dat, out_dat = volumeData, expName, outName, samplesize.outcome = 32860)
	risk_factor_LiverVolume <- rbind(myres, risk_factor_LiverVolume)
}
fwrite(risk_factor_LiverVolume, "risk_factor_LiverVolume_results.csv")

options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
mainMR_diseases <- function(exp_dat, out_dat, expName, outName) {
	mydata <- harmonise_data(exp_dat, out_dat)
	if(length(mydata$SNP[which(mydata$mr_keep == TRUE)]) >0 ) {
		res <- mr(mydata, method_list = c("mr_wald_ratio","mr_ivw", "mr_weighted_median", "mr_egger_regression"))
		res$PRESSO_beta <- NA
		res$PRESSO_pval <- NA
		res$IVW_Qpval <- NA
		res$egger_pleiotropy_pval <- NA
		if(dim(mydata)[1] > 3) {
			presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
									SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
									data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
			res$PRESSO_beta <- presso$`Main MR results`$`Causal Estimate`[2]
			res$PRESSO_pval <- presso$`Main MR results`$`P-value`[2]
			het <- mr_heterogeneity(mydata)
			pleio <- mr_pleiotropy_test(mydata)
			res$IVW_Qpval <- het$Q_pval[2]
			res$egger_pleiotropy_pval <- pleio$pval
		}
		res$expName <- expName
		res$outName <- outName
		return(res)
	}
}

disease_id <- as.vector(read.table("C:/Users/86151/Documents/MyMaster/ieu-diseases.txt", header=F)$V1)

liver_iron_exp <- read.csv("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_iron_IV.csv", header=T)
iron_diseases = {}
for (i in 135:length(disease_id)) {
	out_dat <- extract_outcome_data(snps = liver_iron_exp$SNP, outcomes = disease_id[i], maf_threshold = 0.01)
	if(!is.null(out_dat)) {
		myres <- mainMR_diseases(liver_iron_exp, out_dat, expName = "Liver_iron", outName = disease_id[i])
		iron_diseases = rbind(myres, iron_diseases)
	}
}
fwrite(iron_diseases, "liver_iron_diseases_results.csv")

liver_fat_exp <- read.csv("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_fat_IV.csv", header=T)
fat_diseases = {}
for (i in 1:length(disease_id)) {
	out_dat <- extract_outcome_data(snps = liver_fat_exp$SNP, outcomes = disease_id[i], maf_threshold = 0.01)
	if(!is.null(out_dat)) {
		myres <- mainMR_diseases(liver_fat_exp, out_dat, expName = "Liver_fat", outName = disease_id[i])
		fat_diseases = rbind(myres, fat_diseases)
	}
}
fwrite(fat_diseases, "liver_fat_diseases_results.csv")

liver_volume_exp <- read.csv("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_volume_IV.csv", header=T)
volume_diseases = {}
for (i in 1:length(disease_id)) {
	out_dat <- extract_outcome_data(snps = liver_volume_exp$SNP, outcomes = disease_id[i], maf_threshold = 0.01)
	if(!is.null(out_dat)) {
		myres <- mainMR_diseases(liver_volume_exp, out_dat, expName = "Liver_volume", outName = disease_id[i])
		volume_diseases = rbind(myres, volume_diseases)
	}
}
fwrite(volume_diseases, "liver_volume_diseases_results.csv")



setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent/Clean_Results")
library(data.table)
fat <- read.csv("liverFat_IVW.csv", header=T)
fat$PRESSO_se <- abs(fat$PRESSO_beta/qnorm(fat$PRESSO_pval/2, lower.tail = F))
colnames(fat)[1] <- "Exposure"
iron <- read.csv("liverIron_IVW.csv", header=T)
iron$PRESSO_se <- abs(iron$PRESSO_beta/qnorm(iron$PRESSO_pval/2, lower.tail = F))
colnames(iron)[1] <- "Exposure"
volume <- read.csv("liverVolume_IVW.csv", header=T)
volume$PRESSO_se <- abs(volume$PRESSO_beta/qnorm(volume$PRESSO_pval/2, lower.tail = F))
colnames(volume)[1] <- "Exposure"
fwrite(fat, "liverFat_IVW.csv")
fwrite(iron, "liverIron_IVW.csv")
fwrite(volume, "liverVolume_IVW.csv")

data1 <- read.csv("liverFat_IVW.csv", header = T)
data1 <- data1[order(data1$pval, decreasing = F),]
data1$fdr <- p.adjust(data1$pval, "fdr")
fwrite(data1, "liverFat_IVW.csv")

data2 <- read.csv("liverIron_IVW.csv", header = T)
data2 <- data2[order(data2$pval, decreasing = F),]
data2$fdr <- p.adjust(data2$pval, "fdr")
fwrite(data2, "liverIron_IVW.csv")

data3 <- read.csv("liverVolume_IVW.csv", header = T)
data3 <- data3[order(data3$pval, decreasing = F),]
data3$fdr <- p.adjust(data3$pval, "fdr")
fwrite(data3, "liverVolume_IVW.csv")


##coloc code
options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
library(coloc)
IV_dir <- "C:/Users/86151/Documents/MyMaster/IV_for_MR"

###Liver fat( WHR 33 SNP)
whr <- fread(paste0(IV_dir, "/WHR_IV.csv"), header=T)
fat <- fread("liver_fat.tsv.gz", header=T)

mainMR_risk_factors <- function(exp_dat, out_dat, expName, outName, samplesize.outcome) {
	out_dat <- out_dat[which(out_dat$variant_id%in%exp_dat$SNP),]
	if(!is.null(out_dat)){
		out_dat <- TwoSampleMR::format_data(
					out_dat,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
		mydata <- harmonise_data(exp_dat, out_dat)
		mydata$samplesize.outcome <- samplesize.outcome
		if(dim(mydata)[1] >=1) {
			res <- mr(mydata, method_list = c("mr_wald_ratio","mr_ivw", "mr_weighted_median", "mr_egger_regression"))
			index = {}
			for (j in 1:dim(mydata)[1]) {
				steiger <- directionality_test(mydata[j,])
				if(steiger[1,7] == FALSE) {
					index <- c(index, j)
					}
				}
			mydata = mydata[-index,]
			steiger <- directionality_test(mydata)
			res$exposure_r2 <- steiger[1,5]
			res$outcome_r2 <- steiger[1,6]
			res$direction <- steiger[1,7]
			res$PRESSO_beta <- NA
			res$PRESSO_pval <- NA
			res$IVW_Qpval <- NA
			res$egger_pleiotropy_pval <- NA
			if(dim(mydata)[1] > 3) {
				presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
										SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
										data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
				res$PRESSO_beta <- presso$`Main MR results`$`Causal Estimate`[2]
				res$PRESSO_pval <- presso$`Main MR results`$`P-value`[2]
				het <- mr_heterogeneity(mydata)
				pleio <- mr_pleiotropy_test(mydata)
				res$IVW_Qpval <- het$Q_pval[2]
				res$egger_pleiotropy_pval <- pleio$pval
				}
			res$expName <- expName
			res$outName <- outName
			}
		return(res)
	}
}
expName <- "WHR"
outName <- "Liver_fat"
myres1 <- mainMR_risk_factors(whr, out_dat = fat, expName, outName, samplesize.outcome = 32858)

###Liver iron (iron HbA1c ferritin)
32858
iron = fread("liver_iron.tsv.gz", header = T)
IV_list = paste0(c("iron", "HbA1c", "ferritin"), "_IV.csv")
risk_factor_LiverIron = {}
for (i in 1:length(IV_list)) {
	expName <- gsub(IV_list[i], pattern="_IV.csv", replacement = "")
	outName <- "Liver_iron"
	exp_dat <- read.csv(paste0(IV_dir,"/",IV_list[i]), header = T)
	myres <- mainMR_risk_factors(exp_dat, out_dat = iron, expName, outName, samplesize.outcome = 32858)
	risk_factor_LiverIron <- rbind(myres, risk_factor_LiverIron)
}
###Liver (WC HIP FI 2hGlu)
32860
volume = fread("liver_volume.tsv.gz", header = T)
IV_list = paste0(c("WC", "HIP", "FI", "2hGlu"), "_IV.csv")
risk_factor_LiverVolume = {}
for (i in 1:length(IV_list)) {
	expName <- gsub(IV_list[i], pattern="_IV.csv", replacement = "")
	outName <- "Liver_volume"
	exp_dat <- read.csv(paste0(IV_dir,"/",IV_list[i]), header = T)
	myres <- mainMR_risk_factors(exp_dat, out_dat = volume, expName, outName, samplesize.outcome = 32860)
	risk_factor_LiverVolume <- rbind(myres, risk_factor_LiverVolume)
}

result <- rbind(myres1, risk_factor_LiverIron, risk_factor_LiverVolume)
result$PRESSO_se <- abs(result$PRESSO_beta/qnorm(result$PRESSO_pval/2, lower.tail = F))

options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
library(coloc)
library(ieugwasr)
###for liver fat
#significant traits: FI, WC, BMI, HDL, T2D, body_fat_percent

###for liver iron
#significant traits: transferrin_saturation; TC; HbA1c; HIP; ferritin; drinking; BMI

##for liver volume
#significnat traits: WC; BMI; HIP; SBP; body_fat_percent; LDL; HDL; NAFLD

bmi_gwas <-fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/SNP_gwas_mc_merge_nogc.tbl.uniq", header = T)
#fat = fread("liver_fat.tsv.gz", header=T)
#fat = fread("liver_iron.tsv.gz", header=T)
#fat = fread("liver_volume.tsv.gz", header=T)
bmi_sigSNP <- bmi_gwas[which(bmi_gwas$p< 5e-8),c("SNP", "p")]
fat_sigSNP <- fat[which(fat$p_value< 5e-8),c("variant_id", "p_value")]
colnames(bmi_sigSNP) <- c("rsid", "pval")
colnames(fat_sigSNP) <- c("rsid", "pval")
mySNP <- rbind(bmi_sigSNP, fat_sigSNP)
mySNP <- ld_clump(mySNP, clump_kb = 500, clump_r2 = 0.1)
mySNP <- mySNP[which(mySNP$rsid%in%bmi_gwas$SNP),]
mySNP <- mySNP[which(mySNP$rsid%in%fat$variant_id),]
#fwrite(mySNP, "BMI_Liver_coloc_SNP.csv")
SNP_info <- fat[which(fat$variant_id%in%mySNP$rsid), c("variant_id", "chromosome", "base_pair_location")]
coloc_calc <- function(i) {
	data1 <- with(fat,fat[which(chromosome == SNP_info$chromosome[i]&base_pair_location > SNP_info$base_pair_location[i]-200000&base_pair_location < SNP_info$base_pair_location[i]+200000),])
	data1 <- data1[!duplicated(data1$variant_id),]
	data2 <- bmi_gwas[which(bmi_gwas$SNP%in%data1$variant_id),]
	data2 <- data2[!duplicated(data2$SNP),]
	data1 <- data1[which(data1$variant_id%in%data2$SNP),]
	data1 <- data1[order(data1$variant_id),]
	data2 <- data2[order(data2$SNP),]
	data2$chr <- data1$chromosome
	data2$pos <- data1$base_pair_location
	data1 <- TwoSampleMR::format_data(
					data1,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
	data2 <- TwoSampleMR::format_data(
					data2,
					type = "exposure",
					snp_col = "SNP",
					beta_col = "b",
					se_col = "se",
					eaf_col = "Freq1.Hapmap",
					effect_allele_col = "A1",
					other_allele_col = "A2",
					pval_col = "p"
					#phenotype_col = outName
					)
	mydata <- harmonise_data(data2, data1)
	mydata <- mydata[which(!is.na(mydata$eaf.exposure)&!is.na(mydata$eaf.outcome)&mydata$mr_keep == TRUE),]
	data1 <- mydata[,c("SNP", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome")]
	data2 <- mydata[,c("SNP", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
	data1$maf <- ifelse(data1$eaf.outcome > 0.5, 1-data1$eaf.outcome, data1$eaf.outcome)
	data2$maf <- ifelse(data2$eaf.exposure > 0.5, 1-data2$eaf.exposure, data2$eaf.exposure)
	myres <- coloc.abf(dataset1 = list(snp = data1$SNP, beta = data1$beta.outcome, MAF = data1$maf, N = 344369, pvalues = data1$pval.outcome,type = "quant"),
				dataset2 = list(snp = data2$SNP, beta = data2$beta.exposure, MAF = data2$maf, N = 32860, pvalues = data2$pval.exposure, type = "quant"))
	return(myres)
}

results = {}
for (i in 1:dim(SNP_info)[1]) {
	res <- coloc_calc(i)
	results <- rbind(results, res$summary)
}
results <- as.data.frame(results)
rownames(results) <- SNP_info$variant_id
results$SNP <- rownames(results)
#fwrite(results, "fat_bmi_coloc.csv")
#fwrite(results, "iron_bmi_coloc.csv")
fwrite(results, "volume_bmi_coloc.csv")
snp <- results$SNP[which(results$PP.H4.abf>0.9)]


options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
library(MVMR)
fat = fread("liver_fat.tsv.gz", header=T)
iron = fread("liver_iron.tsv.gz", header=T)
volume = fread("liver_volume.tsv.gz", header=T)
FatMVMR_adjBMI <- function(expID, outName) {
	mvexp <- mv_extract_exposures(c(expID,"ieu-a-2"),clump_r2 = 0.01,clump_kb = 10000, access_token = NULL)
	mvout <- fat[which(fat$variant_id%in%mvexp$SNP),]
	mvout <- format_data(
					mvout,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
	mydata <- mv_harmonise_data(mvexp, mvout)
	nexp <- length(mydata$expname$exposure)
	#res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
	#res1$type <- "mv_lasso"
	res2 <- mv_ivw(mydata)$result
	res2$type <- "mv_ivw"
	F.data <- format_mvmr(BXGs = mydata$exposure_beta,
						  BYG = mydata$outcome_beta,
						  seBXGs = mydata$exposure_se,
						  seBYG = mydata$outcome_se,
						  RSID = rownames(mydata$exposure_beta)
						  )

	sres <- strength_mvmr(r_input = F.data, gencov = 0)
	print(sres)
	pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
	print(pres)
	#res <- ivw_mvmr(r_input = F.data)
	#res$egger_pleiotropy_pval <- pleio$pval
	resA <- res2
	resA$outcome <- outName
	return(resA)
}

IronMVMR_adjBMI <- function(expID, outName) {
	mvexp <- mv_extract_exposures(c(expID,"ieu-a-2"),clump_r2 = 0.01,clump_kb = 10000, access_token = NULL)
	mvout <- iron[which(iron$variant_id%in%mvexp$SNP),]
	mvout <- format_data(
					mvout,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
	mydata <- mv_harmonise_data(mvexp, mvout)
	nexp <- length(mydata$expname$exposure)
	#res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
	#res1$type <- "mv_lasso"
	res2 <- mv_ivw(mydata)$result
	res2$type <- "mv_ivw"
	F.data <- format_mvmr(BXGs = mydata$exposure_beta,
						  BYG = mydata$outcome_beta,
						  seBXGs = mydata$exposure_se,
						  seBYG = mydata$outcome_se,
						  RSID = rownames(mydata$exposure_beta)
						  )

	sres <- strength_mvmr(r_input = F.data, gencov = 0)
	print(sres)
	pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
	print(pres)
	#res <- ivw_mvmr(r_input = F.data)
	#res$egger_pleiotropy_pval <- pleio$pval
	resA <- res2
	resA$outcome <- outName
	return(resA)
}

VolumeMVMR_adjBMI <- function(expID, outName) {
	mvexp <- mv_extract_exposures(c(expID,"ieu-a-2"),clump_r2 = 0.01,clump_kb = 10000, access_token = NULL)
	mvout <- iron[which(volume$variant_id%in%mvexp$SNP),]
	mvout <- format_data(
					mvout,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
	mydata <- mv_harmonise_data(mvexp, mvout)
	nexp <- length(mydata$expname$exposure)
	#res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
	#res1$type <- "mv_lasso"
	res2 <- mv_ivw(mydata)$result
	res2$type <- "mv_ivw"
	F.data <- format_mvmr(BXGs = mydata$exposure_beta,
						  BYG = mydata$outcome_beta,
						  seBXGs = mydata$exposure_se,
						  seBYG = mydata$outcome_se,
						  RSID = rownames(mydata$exposure_beta)
						  )

	sres <- strength_mvmr(r_input = F.data, gencov = 0)
	print(sres)
	pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
	print(pres)
	#res <- ivw_mvmr(r_input = F.data)
	#res$egger_pleiotropy_pval <- pleio$pval
	resA <- res2
	resA$outcome <- outName
	return(resA)
}

#FI, WC, HDL, T2D and BFP
mvmr_fat_id <- c("ieu-b-115", "ieu-a-66", "ieu-a-299", "ieu-a-23", "ebi-a-GCST003435") #FI, WC, HDL, T2D and BFP
fat_mvmr_res <- {}
for (id in mvmr_fat_id) {
	res <- FatMVMR_adjBMI(id, "liver_fat")
	fat_mvmr_res <- rbind(fat_mvmr_res, res)
}
#saturation; TC; HbAa1c;Hip;ferritin; drinking 
mvmr_iron_id <- c("ieu-a-1051", "ieu-a-301", "ieu-b-103", "ieu-a-54", "ieu-a-1050", "ieu-b-73")
iron_mvmr_res <- {}
for (id in mvmr_iron_id) {
	res <- IronMVMR_adjBMI(id, "liver_iron")
	iron_mvmr_res <- rbind(iron_mvmr_res, res)
}
#wc; (hip); sbp; bfp; LDL; HDL-C (NAFLD).
mvmr_volume_id <- c("ieu-a-66", "ieu-b-38","ieu-a-54", "ebi-a-GCST003435", "ieu-a-300", "ieu-a-299")
volume_mvmr_res <- {}
for (id in mvmr_volume_id) {
	res <- VolumeMVMR_adjBMI(id, "liver_volume")
	volume_mvmr_res <- rbind(volume_mvmr_res, res)
}
data = rbind(fat_mvmr_res, iron_mvmr_res, volume_mvmr_res)


bmi_gwas <-fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/SNP_gwas_mc_merge_nogc.tbl.uniq", header = T)
#nafld_gwas <- fread("C:/Users/86151/Documents/DingXX/GCST90011885_buildGRCh37.tsv.gz", header=T)
#nafld_gwas <- fread("C:/2345Downloads/MAGIC1000G_FI_EUR.tsv.gz", header=T)
#fat = fread("liver_fat.tsv.gz", header=T)
#fat = fread("liver_iron.tsv.gz", header=T)
#fat = fread("liver_volume.tsv.gz", header=T)
bmi_sigSNP <- bmi_gwas[which(bmi_gwas$p< 5e-8),c("SNP", "p")]
nafld_sigSNP <- nafld_gwas[which(nafld_gwas$p_value< 5e-8),c("variant", "p_value")]
colnames(bmi_sigSNP) <- c("rsid", "pval")
colnames(nafld_sigSNP) <- c("rsid", "pval")
mySNP <- rbind(bmi_sigSNP, nafld_sigSNP)
mySNP <- ieugwasr::ld_clump(mySNP, clump_kb = 500, clump_r2 = 0.1)
mySNP <- mySNP[which(mySNP$rsid%in%bmi_gwas$SNP),]
mySNP <- mySNP[which(mySNP$rsid%in%nafld_gwas$variant),]
bmi_dat <- bmi_gwas[which(bmi_gwas$SNP%in%mySNP$rsid),]
bmi_dat$pheno <- "BMI"
nafld_dat <- nafld_gwas[which(nafld_gwas$variant%in%mySNP$rsid),]
nafld_dat$pheno <- "FI"
bmi_dat <- TwoSampleMR::format_data(
					bmi_dat,
					type = "exposure",
					snp_col = "SNP",
					beta_col = "b",
					se_col = "se",
					eaf_col = "Freq1.Hapmap",
					effect_allele_col = "A1",
					other_allele_col = "A2",
					pval_col = "p",
					phenotype_col = "pheno"
					)
nafld_dat <- TwoSampleMR::format_data(
				nafld_dat,
				type = "exposure",
				snp_col = "variant",
				beta_col = "beta",
				se_col = "standard_error",
				eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value",
				phenotype_col = "pheno"
				)
#nafld_dat <- nafld_dat[,colnames(bmi_dat)]
mvexp <- rbind(bmi_dat, nafld_dat)
#mvout <- volume[which(volume$variant_id%in%mvexp$SNP),]
mvout <- fat[which(fat$variant_id%in%mvexp$SNP),]
mvout <- format_data(
					mvout,
					type = "outcome",
					snp_col = "variant_id",
					beta_col = "beta",
					se_col = "standard_error",
					eaf_col = "effect_allele_frequency",
					effect_allele_col = "effect_allele",
					other_allele_col = "other_allele",
					pval_col = "p_value"
					#phenotype_col = outName
					)
mydata <- mv_harmonise_data(mvexp, mvout)
nexp <- length(mydata$expname$exposure)
res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
res1$type <- "mv_lasso"
res2 <- mv_ivw(mydata)$result
res2$type <- "mv_ivw"
data = rbind(fat_mvmr_res, iron_mvmr_res, volume_mvmr_res, res2)

options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
library(MVMR)
mainMR_diseases <- function(exp_dat, out_dat, expName, outName) {
	out_dat <- out_dat[which(out_dat$variant_id%in%exp_dat$SNP),]
	out_dat <- TwoSampleMR::format_data(
				out_dat,
				type = "outcome",
				snp_col = "variant_id",
				beta_col = "lnOR",
				se_col = "standard_error",
				#eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value",
				phenotype_col = "pheno"
				)
	mydata <- harmonise_data(exp_dat, out_dat)
	if(length(mydata$SNP[which(mydata$mr_keep == TRUE)]) >0 ) {
		res <- mr(mydata, method_list = c("mr_wald_ratio","mr_ivw", "mr_weighted_median", "mr_egger_regression"))
		res$PRESSO_beta <- NA
		res$PRESSO_pval <- NA
		res$IVW_Qpval <- NA
		res$egger_pleiotropy_pval <- NA
		if(dim(mydata)[1] > 3) {
			presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
									SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
									data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
			res$PRESSO_beta <- presso$`Main MR results`$`Causal Estimate`[2]
			res$PRESSO_pval <- presso$`Main MR results`$`P-value`[2]
			het <- mr_heterogeneity(mydata)
			pleio <- mr_pleiotropy_test(mydata)
			res$IVW_Qpval <- het$Q_pval[2]
			res$egger_pleiotropy_pval <- pleio$pval
		}
		res$expName <- expName
		res$outName <- outName
		return(res)
	}
}

fat_IV <- fread("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_fat_IV.csv")
iron_IV <- fread("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_iron_IV.csv")
volume_IV <- fread("C:/Users/86151/Documents/MyMaster/IV_for_MR/Liver_volume_IV.csv")
nafld_gwas <- fread("C:/Users/86151/Documents/DingXX/GCST90011885_buildGRCh37.tsv.gz", header=T)
res1 <- mainMR_diseases(fat_IV, nafld_gwas, "Liver_fat", "NAFLD") 
res2 <- mainMR_diseases(iron_IV, nafld_gwas, "Liver_iron", "NAFLD") 
res3 <- mainMR_diseases(volume_IV, nafld_gwas, "Liver_volume", "NAFLD") 

ieuID <- fread("disease_gwasInfo.csv", header = T)
data1 = fread("liver_fat_diseases_results.csv", header = T)
data2 = fread("liver_iron_diseases_results.csv", header = T)
data3 = fread("liver_volume_diseases_results.csv", header = T)
data1 = data1[which(data1$outName%in%c(ieuID$id,"NAFLD")),]
data2 = data2[which(data2$outName%in%c(ieuID$id,"NAFLD")),]
data3 = data3[which(data3$outName%in%c(ieuID$id,"NAFLD")),]
data1$PRESSO_se <- abs(data1$PRESSO_beta/qnorm(data1$PRESSO_pval/2, lower.tail = F))
data2$PRESSO_se <- abs(data2$PRESSO_beta/qnorm(data2$PRESSO_pval/2, lower.tail = F))
data3$PRESSO_se <- abs(data3$PRESSO_beta/qnorm(data3$PRESSO_pval/2, lower.tail = F))
fwrite(data1, "liver_fat_diseases_results.csv")
fwrite(data2, "liver_iron_diseases_results.csv")
fwrite(data3, "liver_volume_diseases_results.csv")


setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent/Clean_Results")
library(data.table)
data1 <- read.csv("liverFat_disease_IVW.csv", header = T)
data1 <- data1[order(data1$pval, decreasing = F),]
data1$fdr <- p.adjust(data1$pval, "fdr")
fwrite(data1, "liverFat_disease_IVW.csv")

data2 <- read.csv("liverIron_disease_IVW.csv", header = T)
data2 <- data2[order(data2$pval, decreasing = F),]
data2$fdr <- p.adjust(data2$pval, "fdr")
fwrite(data2, "liverIron_disease_IVW.csv")

data3 <- read.csv("liverVolume_disease_IVW.csv", header = T)
data3 <- data3[order(data3$pval, decreasing = F),]
data3$fdr <- p.adjust(data3$pval, "fdr")
fwrite(data3, "liverVolume_disease_IVW.csv")


options(stringsAsFactors = F)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent")
library(data.table)
library(TwoSampleMR)
library(MVMR)
LiverMVMR_adjBMI <- function(mvexp, outName) {
	mvout <- extract_outcome_data(mvexp$SNP, outName, access_token = NULL)
	mydata <- mv_harmonise_data(mvexp, mvout)
	nexp <- length(mydata$expname$exposure)
	#res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
	#res1$type <- "mv_lasso"
	res2 <- mv_ivw(mydata)$result
	res2$type <- "mv_ivw"
	F.data <- format_mvmr(BXGs = mydata$exposure_beta,
						  BYG = mydata$outcome_beta,
						  seBXGs = mydata$exposure_se,
						  seBYG = mydata$outcome_se,
						  RSID = rownames(mydata$exposure_beta)
						  )

	sres <- strength_mvmr(r_input = F.data, gencov = 0)
	print(sres)
	pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
	print(pres)
	#res <- ivw_mvmr(r_input = F.data)
	#res$egger_pleiotropy_pval <- pleio$pval
	resA <- res2
	#resA$outcome <- outName
	return(resA)
}

LiverMVMR_adjBMI_NAFLD <- function(mvexp, outName) {
	mvout <- nafld_gwas[which(nafld_gwas$variant_id%in%mvexp$SNP),]
	mvout <- TwoSampleMR::format_data(
				mvout,
				type = "outcome",
				snp_col = "variant_id",
				beta_col = "lnOR",
				se_col = "standard_error",
				#eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value"
				#phenotype_col = "pheno"
				)
	mydata <- mv_harmonise_data(mvexp, mvout)
	nexp <- length(mydata$expname$exposure)
	#res1 <- mv_subset(mydata, features = mv_lasso_feature_selection(mydata))$result
	#res1$type <- "mv_lasso"
	res2 <- mv_ivw(mydata)$result
	res2$type <- "mv_ivw"
	F.data <- format_mvmr(BXGs = mydata$exposure_beta,
						  BYG = mydata$outcome_beta,
						  seBXGs = mydata$exposure_se,
						  seBYG = mydata$outcome_se,
						  RSID = rownames(mydata$exposure_beta)
						  )

	sres <- strength_mvmr(r_input = F.data, gencov = 0)
	print(sres)
	pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
	print(pres)
	#res <- ivw_mvmr(r_input = F.data)
	#res$egger_pleiotropy_pval <- pleio$pval
	resA <- res2
	resA$outcome <- outName
	return(resA)
}

##MVMR in liver content adjusting BMI
bmi_gwas <-fread("C:/Users/86151/Documents/DingJ/colorectalCa_MR/SNP_gwas_mc_merge_nogc.tbl.uniq", header = T)
fat_gwas = fread("liver_fat.tsv.gz", header=T)
iron_gwas = fread("liver_iron.tsv.gz", header=T)
volume_gwas = fread("liver_volume.tsv.gz", header=T)

MVMR_IV <- function(liver_gwas, name) {
	bmi_sigSNP <- bmi_gwas[which(bmi_gwas$p< 5e-8),c("SNP", "p")]
	liver_sigSNP <- liver_gwas[which(liver_gwas$p_value< 5e-8),c("variant_id", "p_value")]
	colnames(bmi_sigSNP) <- c("rsid", "pval")
	colnames(liver_sigSNP) <- c("rsid", "pval")
	mySNP <- rbind(bmi_sigSNP, liver_sigSNP)
	mySNP <- ieugwasr::ld_clump(mySNP, clump_kb = 500, clump_r2 = 0.1)
	mySNP <- mySNP[which(mySNP$rsid%in%bmi_gwas$SNP),]
	mySNP <- mySNP[which(mySNP$rsid%in%liver_gwas$variant_id),]
	bmi_dat <- bmi_gwas[which(bmi_gwas$SNP%in%mySNP$rsid),]
	bmi_dat$pheno <- "BMI"
	liver_dat <- liver_gwas[which(liver_gwas$variant_id%in%mySNP$rsid),]
	liver_dat$pheno <- name
	bmi_dat <- TwoSampleMR::format_data(
					bmi_dat,
					type = "exposure",
					snp_col = "SNP",
					beta_col = "b",
					se_col = "se",
					eaf_col = "Freq1.Hapmap",
					effect_allele_col = "A1",
					other_allele_col = "A2",
					pval_col = "p",
					phenotype_col = "pheno"
					)
	liver_dat <- TwoSampleMR::format_data(
				liver_dat,
				type = "exposure",
				snp_col = "variant_id",
				beta_col = "beta",
				se_col = "standard_error",
				eaf_col = "effect_allele_frequency",
				effect_allele_col = "effect_allele",
				other_allele_col = "other_allele",
				pval_col = "p_value",
				phenotype_col = "pheno"
				)
	bmi_fat <- rbind(bmi_dat, liver_dat)
	return(bmi_fat)
}

bmi_fat <- MVMR_IV(fat_gwas, "Liver_fat")
bmi_iron <- MVMR_IV(iron_gwas, "Liver_iron")
bmi_volume <- MVMR_IV(volume_gwas, "Liver_volume")
ieuID <- as.vector(fread("disease_gwasInfo.csv", header = T)$id)
fat_mvmr_res = {}
for (i in 58:length(ieuID)) {
	res <- LiverMVMR_adjBMI(bmi_fat, ieuID[i])
	fat_mvmr_res <- rbind(fat_mvmr_res, res)
}
fwrite(fat_mvmr_res, "fatadjBMI_disease.csv")

iron_mvmr_res = {}
for (i in 40:length(ieuID)) {
	res <- LiverMVMR_adjBMI(bmi_iron, ieuID[i])
	iron_mvmr_res <- rbind(iron_mvmr_res, res)
}
fwrite(iron_mvmr_res, "ironadjBMI_disease.csv")

volume_mvmr_res = {}
for (i in 30:length(ieuID)) {
	res <- LiverMVMR_adjBMI(bmi_volume, ieuID[i])
	volume_mvmr_res <- rbind(volume_mvmr_res, res)
}
fwrite(volume_mvmr_res, "volumeadjBMI_disease.csv")

nafld_gwas <- fread("C:/Users/86151/Documents/DingXX/GCST90011885_buildGRCh37.tsv.gz", header=T)
res1 <- LiverMVMR_adjBMI_NAFLD(bmi_fat, "NAFLD")
res2 <- LiverMVMR_adjBMI_NAFLD(bmi_iron, "NAFLD")
res3 <- LiverMVMR_adjBMI_NAFLD(bmi_volume, "NAFLD")

###forest plot （remove internalizing problem)
setwd("C:/Users/86151/Documents/MyMaster/study1_liverContent/Clean_Results")
library(data.table)
library(forestplot)
mydata1 <- read.csv("liverFat_IVW.csv", header=T, stringsAsFactors=F)
mydata1 <- mydata1[,c(1:5)]
colnames(mydata1) <- c('Risk factor', 'NSNP', 'BETA', 'SE','P')
mydata1$`95%LCI` <- with(mydata1, BETA-1.96*SE)
mydata1$`95%UCI` <- with(mydata1, BETA+1.96*SE)

mydata2 <- read.csv("liverFat_disease_plot.csv", header=T, stringsAsFactors=F)
mydata2 <- mydata2[,c(1:5)]
colnames(mydata2) <- c('Disease', 'NSNP', 'OR', '95%LCI','95%UCI')

{## code for forestplot of cholelithiasis
tab1 <- mydata1[,c('Risk factor','NSNP','BETA','SE','P')]
lab1 <- c('Risk factor','NSNP','BETA','SE','P')
outcome1 <- c(tab1$`Risk factor`)
nsnp1 = mydata1$NSNP
mean1 = c(sprintf("%.3f",mydata1[,"BETA"]))
se1 = c(sprintf("%.3f",mydata1[,"SE"]))
p1 = c(sprintf("%.3f",mydata1[,"P"]))
newdata1 = cbind(outcome1,nsnp1,mean1,se1,p1)
tabletext1 <- rbind(lab1,newdata1)

tab2 <- mydata2[,c('Disease','NSNP','OR','95%LCI','95%UCI')]
lab2 <- c('Disease','NSNP','OR','95%LCI','95%UCI')
outcome2 <- c(tab2$Disease)
nsnp2 = mydata2$NSNP
mean2 = c(sprintf("%.3f",mydata2[,"OR"]))
lower2 = c(sprintf("%.3f",mydata2[,"95%LCI"]))
upper2 = c(sprintf("%.3f",mydata2[,"95%UCI"]))
newdata2 = cbind(outcome2,nsnp2,mean2,lower2,upper2)
tabletext2 <- rbind(lab2,newdata2)

#par(mfrow=c(2,2))
pdf("Figure2.pdf",width = 10, height = 14)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
pushViewport(viewport(layout.pos.row = 1))
forestplot(tabletext1, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "A. Liver fat content as the outcome",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
           lwd.ci = 3, lty.ci = "solid",
		   mean = c(NA,as.numeric(sprintf("%.3f",mydata1[,"BETA"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata1[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata1[,"95%UCI"]))),
           clip =c(-0.2,1.3),
		   xlog=F,
           col=fpColors(box=c("black")),
           grid = structure(c(-0.2,1.4), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(-0.2,-0.1,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4),
           xlab="BETA",
		   new_page = FALSE)
popViewport()
pushViewport(viewport(layout.pos.row = 2))
forestplot(tabletext2, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "B. Liver fat content as the exposure",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
		   lwd.ci = 3, lty.ci = "solid",
           mean = c(NA,as.numeric(sprintf("%.3f",mydata2[,"OR"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata2[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata2[,"95%UCI"]))),
           clip =c(0.6,12.2),
		   xlog=T,
           col=fpColors(box=c("black")),
           grid = structure(c(0.6,12.2), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(0.6,0.7, 0.8, 0.9, 1.0,1.1, 1.2,2.4,6,12.2),
           xlab="Odds ratio",
		   new_page = FALSE)
popViewport(2)
dev.off()
}

mydata3 <- read.csv("liverIron_IVW.csv", header=T, stringsAsFactors=F)
mydata3 <- mydata3[,c(1:5)]
colnames(mydata3) <- c('Risk factor', 'NSNP', 'BETA', 'SE','P')
mydata3$`95%LCI` <- with(mydata3, BETA-1.96*SE)
mydata3$`95%UCI` <- with(mydata3, BETA+1.96*SE)

mydata4 <- read.csv("liverIron_disease_plot.csv", header=T, stringsAsFactors=F)
mydata4 <- mydata4[,c(1:5)]
colnames(mydata4) <- c('Disease', 'NSNP', 'OR', '95%LCI','95%UCI')

{## code for forestplot of cholelithiasis
tab3 <- mydata3[,c('Risk factor','NSNP','BETA','SE','P')]
lab3 <- c('Risk factor','NSNP','BETA','SE','P')
outcome3 <- c(tab3$`Risk factor`)
nsnp3 = mydata3$NSNP
mean3 = c(sprintf("%.3f",mydata3[,"BETA"]))
se3 = c(sprintf("%.3f",mydata3[,"SE"]))
p3 = c(sprintf("%.3f",mydata3[,"P"]))
newdata3 = cbind(outcome3,nsnp3,mean3,se3,p3)
tabletext3 <- rbind(lab3,newdata3)

tab4 <- mydata4[,c('Disease','NSNP','OR','95%LCI','95%UCI')]
lab4 <- c('Disease','NSNP','OR','95%LCI','95%UCI')
outcome4 <- c(tab4$Disease)
nsnp4 = mydata4$NSNP
mean4 = c(sprintf("%.3f",mydata4[,"OR"]))
lower4 = c(sprintf("%.3f",mydata4[,"95%LCI"]))
upper4 = c(sprintf("%.3f",mydata4[,"95%UCI"]))
newdata4 = cbind(outcome4,nsnp4,mean4,lower4,upper4)
tabletext4 <- rbind(lab4,newdata4)

#par(mfrow=c(2,2))
pdf("Figure3.pdf",width = 10, height = 14)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
pushViewport(viewport(layout.pos.row = 1))
forestplot(tabletext3, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "A. Liver iron content as the outcome",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
           lwd.ci = 3, lty.ci = "solid",
		   mean = c(NA,as.numeric(sprintf("%.3f",mydata3[,"BETA"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata3[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata3[,"95%UCI"]))),
           clip =c(-0.3,3.0),
		   xlog=F,
           col=fpColors(box=c("black")),
           grid = structure(c(-0.3,3.0), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(-0.3,-0.2, -0.1, 0.0,0.2,0.4,0.6,0.8,1.0,2.0,3.0),
           xlab="BETA",
		   new_page = FALSE)
popViewport()
pushViewport(viewport(layout.pos.row = 2))
forestplot(tabletext4, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "B. Liver iron content as the exposure",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
		   lwd.ci = 3, lty.ci = "solid",
           mean = c(NA,as.numeric(sprintf("%.3f",mydata4[,"OR"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata4[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata4[,"95%UCI"]))),
           clip =c(0.5,1.2),
		   xlog=T,
           col=fpColors(box=c("black")),
           grid = structure(c(0.5,1.2), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2),
           xlab="Odds ratio",
		   new_page = FALSE)
popViewport(2)
dev.off()
}


mydata5 <- read.csv("liverVolume_IVW.csv", header=T, stringsAsFactors=F)
mydata5 <- mydata5[,c(1:5)]
colnames(mydata5) <- c('Risk factor', 'NSNP', 'BETA', 'SE','P')
mydata5$`95%LCI` <- with(mydata5, BETA-1.96*SE)
mydata5$`95%UCI` <- with(mydata5, BETA+1.96*SE)

mydata6 <- read.csv("liverVolume_disease_plot.csv", header=T, stringsAsFactors=F)
mydata6 <- mydata6[,c(1:5)]
colnames(mydata6) <- c('Disease', 'NSNP', 'OR', '95%LCI','95%UCI')

{## code for forestplot of cholelithiasis
tab5 <- mydata5[,c('Risk factor','NSNP','BETA','SE','P')]
lab5 <- c('Risk factor','NSNP','BETA','SE','P')
outcome5 <- c(tab5$`Risk factor`)
nsnp5 = mydata5$NSNP
mean5 = c(sprintf("%.3f",mydata5[,"BETA"]))
se5 = c(sprintf("%.3f",mydata5[,"SE"]))
p5 = c(sprintf("%.3f",mydata5[,"P"]))
newdata5 = cbind(outcome5,nsnp5,mean5,se5,p5)
tabletext5 <- rbind(lab5,newdata5)

tab6 <- mydata6[,c('Disease','NSNP','OR','95%LCI','95%UCI')]
lab6 <- c('Disease','NSNP','OR','95%LCI','95%UCI')
outcome6 <- c(tab6$Disease)
nsnp6 = mydata6$NSNP
mean6 = c(sprintf("%.3f",mydata6[,"OR"]))
lower6 = c(sprintf("%.3f",mydata6[,"95%LCI"]))
upper6 = c(sprintf("%.3f",mydata6[,"95%UCI"]))
newdata6 = cbind(outcome6,nsnp6,mean6,lower6,upper6)
tabletext6 <- rbind(lab6,newdata6)

#par(mfrow=c(2,2))
pdf("Figure4.pdf",width = 10, height = 14)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
pushViewport(viewport(layout.pos.row = 1))
forestplot(tabletext5, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "A. Liver volume as the outcome",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
           lwd.ci = 3, lty.ci = "solid",
		   mean = c(NA,as.numeric(sprintf("%.3f",mydata5[,"BETA"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata5[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata5[,"95%UCI"]))),
           clip =c(-0.4,1.2),
		   xlog=F,
           col=fpColors(box=c("black")),
           grid = structure(c(-0.4,1.2), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2),
           xlab="BETA",
		   new_page = FALSE)
popViewport()
pushViewport(viewport(layout.pos.row = 2, width =10, height = 4 ))
forestplot(tabletext6, 
           #legend = c("Discovery", "Replication"),
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
		    txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                          gpar(fontfamily = "",col = "black")
										  ),
                            ticks = gpar(fontfamily = "", cex=0.8),
                            xlab  = gpar(fontfamily = "", cex = 1.1)),
		   title = "B. Liver fat content as the exposure",
		   #title = "B. Randomization Stage for Target Diseases in T2D Patients",
		   is.summary=c(TRUE,rep(FALSE,26)),
		   #is.summary=c(TRUE,rep(FALSE,10)),
           boxsize = .2, # We set the box size to better visualize the type
           line.margin = 0.1, # We need to add this to avoid crowding
		   lwd.ci = 3, lty.ci = "solid",
           mean = c(NA,as.numeric(sprintf("%.3f",mydata6[,"OR"]))),
           lower = c(NA,as.numeric(sprintf("%.3f",mydata6[,"95%LCI"]))),
           upper = c(NA,as.numeric(sprintf("%.3f",mydata6[,"95%UCI"]))),
           clip =c(0.7,6.6),
		   xlog=T,
           col=fpColors(box=c("black")),
           grid = structure(c(0.7,6.6), 
                            gp = gpar(lty = 2, col = "#CCCCFF")),
		   xticks = c(0.7, 0.8, 0.9, 1.0,1.1,1.2,2.0,2.5,5,6.6),
           xlab="Odds ratio",
		   new_page = FALSE)
popViewport(2)
dev.off()
}

##coloc code
options(stringsAsFact
The sample size of glycemic traits is 281416.
The sample size of body fat percentage is 65831. The sample size of body lean mass is 8327.
The sample size of circulating leptin level is 33987. The sample of iron metabolism is 23986. 
The sample size of CHD is ~100,000. The sample size of adiponectin is 51874. 
The sample size of NAFLD is 4804 adults (PMID: 21423719). 
The sample size of drinking is ~500000 and that of smoking > 200000. 
The sample size of T2D is ~ 100000.





