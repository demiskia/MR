### Hello! start by giving the name of your exposure trait and set the working directory ####

exposure <- "exposure_name"
directory <- "directory_name"
scatter_y_label <- "Exposure name that will be presented in the plot"

setwd("working_directory_name")
dir.create(sprintf("%s", directory))
dir.create(sprintf("%s/results", directory))
dir.create(sprintf("%s/figures", directory))
dir.create(sprintf("%s/qc", directory))

### Read in your files...
exposure_snps <- read.delim("exposure_file", header=TRUE, sep="\t", stringsAsFactors=FALSE)
pd_snps <-read.table("pd_gwas_file", header=TRUE, sep="\t", stringsAsFactors=FALSE)

### Load in required libraries (plyr tends to make data handling easier)
# library(plyr)
# install.packages("ggplot2")
# install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

### Let's sort out the exposure dataset first ####
head(exposure_snps)

exposure_snps$locus <- NULL
exposure_snps$Gene <- NULL
exposure_snps$P <- NULL
# exposure_snps$MAF <- NULL
exposure_snps$Trait <- NULL

head(exposure_snps)

### And, again, change the column names
names(exposure_snps) <- c('MarkerName','exp_EA','exp_OA','eaf_exp','beta_exp','se_exp')
# names(urate_snps) <- c('MarkerName','EA_pd','OA_pd','beta_pd','se_pd')
head(exposure_snps)
dim(exposure_snps)

### Now let's sort out the PD dataset...
head(pd_snps)

### 1. Read in the PD file summary file. Here, for each SNP, we need:
# i) a column for SNP id (rsid),
# ii) a column with the effect allele (effect_allele_pd)
# iii) a column with the other allele (other_allele_pd)
# iv) a column for the beta (beta_pd)
# v) a column for the standard error (standard_error_pd)
### Include header=TRUE, if the input file has a header, and replace sep="\t" if tab-delimited



### Check that the file has been read in correctly
### If unnecessary columns, can get rid of them with pd_snps$columnname <- NULL
pd_snps$FreqSE <- NULL
pd_snps$MinFreq <- NULL
pd_snps$MaxFreq <- NULL
pd_snps$P.value <- NULL
pd_snps$Direction <- NULL
pd_snps$HetChiSq <- NULL
pd_snps$df <- NULL
pd_snps$P.value.1 <- NULL
pd_snps$Chr <- NULL
pd_snps$Bp <- NULL
head(pd_snps)


### 2. Change column names for consistency
names(pd_snps) <- c('MarkerName','EA_pd','OA_pd',"eaf_pd",'beta_pd','se_pd')
# names(pd_snps) <- c('MarkerName','exp_EA','exp_OA','se_exp','beta_exp')


### 5. Now, let's merge our two tables into one complete table, using the join() function from the plyr library. If the naming of the columns has been done consistently, it should match them based on the only matchin column (rsid).
# R also has the merge() function for this, but I tend to avoid it, as it does not preserve the order of the entries, which in some situations can become problematic.
mr_dataset <- merge(exposure_snps,pd_snps, by="MarkerName")

# Let's check our table with the head() and tail() function. This should look nice and clean now. If there are weird extra columns that have appeared, use mr_dataset$weirdcolumn <- NULL to get rid of it.
head(mr_dataset)
tail(mr_dataset)
dim(mr_dataset)

# #Change upper case and lower case
mr_dataset$exp_EA[mr_dataset$exp_EA=="C"]<-"c"
mr_dataset$exp_EA[mr_dataset$exp_EA=="T"]<-"t"
mr_dataset$exp_EA[mr_dataset$exp_EA=="A"]<-"a"
mr_dataset$exp_EA[mr_dataset$exp_EA=="G"]<-"g"
mr_dataset$exp_OA[mr_dataset$exp_OA=="C"]<-"c"
mr_dataset$exp_OA[mr_dataset$exp_OA=="T"]<-"t"
mr_dataset$exp_OA[mr_dataset$exp_OA=="A"]<-"a"
mr_dataset$exp_OA[mr_dataset$exp_OA=="G"]<-"g"


mr_dataset$EA_pd[mr_dataset$EA_pd=="C"]<-"c"
mr_dataset$EA_pd[mr_dataset$EA_pd=="T"]<-"t"
mr_dataset$EA_pd[mr_dataset$EA_pd=="A"]<-"a"
mr_dataset$EA_pd[mr_dataset$EA_pd=="G"]<-"g"
mr_dataset$OA_pd[mr_dataset$OA_pd=="C"]<-"c"
mr_dataset$OA_pd[mr_dataset$OA_pd=="T"]<-"t"
mr_dataset$OA_pd[mr_dataset$OA_pd=="A"]<-"a"
mr_dataset$OA_pd[mr_dataset$OA_pd=="G"]<-"g"

head(mr_dataset)
tail(mr_dataset)
dim(mr_dataset)

### 6. Now, let's check the dimensions of our complete dataset. This should, naturally, have the number of rows equivalent to our number of SNPs, and 9 columns.
dim(mr_dataset)
# We'll also now create some variables that store values that we'll be using often. First, we'll make one for the number of our alleles in our instrument (called n_snps)
n_snps <- nrow(mr_dataset)
# Next, we'll make a vector of our SNP IDs:
snp_ids <- mr_dataset$MarkerName

### 7. Now, let's check that we don't have any funny SNPs where the SNPs don't match at all. This can be the case if one GWAS reports the results using the alleles on the forward strand, and the other GWAS using the alleles on the reverse strand. This shouldn't be an issue, but just to be extra safe, the for-loop below checks that the alleles match, and returns a logical vector (bad_alleles) that says TRUE if the alleles match, and FALSE if they don't. (I haven't tested this, so fingers crossed it doesn't have a typo and works fine...)
bad_alleles <- vector(length=n_snps)
for (i in 1:n_snps){
  bad_alleles[i] <- ((as.character(mr_dataset$EA_pd[i]) == as.character(mr_dataset$exp_EA[i])) & (as.character(mr_dataset$OA_pd[i]) == as.character(mr_dataset$exp_OA[i]))) | ((as.character(mr_dataset$EA_pd[i]) == as.character(mr_dataset$exp_OA[i])) & (as.character(mr_dataset$OA_pd[i]) == as.character(mr_dataset$exp_EA[i])))
}

# Now, using the nifty all() function, we can check if all the values in our logical vector are TRUE. If the result if this function is TRUE, that means we're safe. If FALSE, it means we have a problem somewhere.
all(bad_alleles)

# To check where the potential problem is, this command returns the dodgy rows, where the alleles in PD GWAS don't match the alleles in the intermediate phenotype dataset. I would remove these for now, and then we can manually check them from ensembl (or other equivalent source) to reinclude them correctly.
mr_dataset[which(!bad_alleles),]

### 8. Now, let's get rid of other strand ambiguities and palindromic SNPs. I would initially check the SNPs manually for palindromic SNPs (A/T, C/G), and remove them from initial analysis. These can then later be checked from ensembl, made sure they're included correctly, then readded to the dataset. We can look at the full thing using the View() function (remember capital V).
# View(mr_dataset)

#another method for flagging the palindromic SNPs
palindromic_snps <- mr_dataset[(mr_dataset$exp_EA %in% c("g","c") & mr_dataset$exp_OA %in% c("g","c")) | (mr_dataset$exp_EA %in% c("t","a") & mr_dataset$exp_OA %in% c("t","a")), ]

write.table(palindromic_snps, file=sprintf("%s/qc/%s_palindromic_snps.txt", directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)

# And to remove palindromic SNPs, make a note of the the rows for these, then do:
# mr_dataset_no_palindromic <- mr_dataset[-c(row,numbers,separated,by,commas),]

### 9. Note, at this stage, we have to be sure that the beta_pd column is in fact betas and not odds ratios. If odds ratios, this can easily be corrected by mr_dataset$beta_pd <- log(mr_dataset$beta_pd) The log() function takes the natural logarithm by default. Similarly, if we only have confidence intervals but no standard errors, the standard error can be derived (to a rough approximation) from the confidence intervals, but I'm trusting that the dataset has the original standard errors included, so there's probably no need.

### 10. Forgot to mention, same thing obviously applies to the intermediate phenotype odds ratios / betas, if a binary trait (smoking initiation etc), so repeat as above if necessary.

### 11. Now that we have our nice dataset with betas, we make one final check to ensure the directions of effect are consistent between our pd and intermediate phenotype data. This is easy, as the betas can simply be multiplied by -1 to get the effect estimate with respect to the other allele, if the alleles have been reversed. We can achieve this with this little for-loop. As we have already checked that the allele pairs match, a non-matching pair here always means that the alleles have flipped:
head(mr_dataset)

for(i in 1:n_snps){
    if (as.character(mr_dataset$EA_pd[i]) != as.character(mr_dataset$exp_EA[i])){
        mr_dataset$beta_exp[i] <- -mr_dataset$beta_exp[i]
    }
}
head(mr_dataset)

# Finally, before getting on to actual analysis, we'll ensure that the exposures are positive, while maintaining the relationship between the directions of effect for the exposure and the outcome
mr_dataset$beta_pd <- mr_dataset$beta_pd*sign(mr_dataset$beta_exp)
mr_dataset$beta_exp <- abs(mr_dataset$beta_exp)

head(mr_dataset)

### 12. We are now ready to roll! The dataset at this stage should be nice and clean. Probably worth checking everything manually just in case, but there should be no inconsistencies. Hopefully, no SNPs have been removed because of unmatching alleles or strand ambiguities.

# Let's use similar notation to Bowden/Burgess and friends:

BetaYG <- mr_dataset$beta_pd
BetaXG <- mr_dataset$beta_exp
seBetaYG <- mr_dataset$se_pd
seBetaXG <- mr_dataset$se_exp

# Let's first create a vector of associations for individual snps

individual_betas <- vector(length=n_snps)
for (i in 1:n_snps){
  individual_betas[i] <- BetaYG[i]/BetaXG[i]
}
individual_ORs <- exp(individual_betas)

# Now, we need the variances, using the delta method:
individual_standard_errors <- vector(length=n_snps)
for (i in 1:n_snps){
  individual_standard_errors[i] <- sqrt(((seBetaYG[i]^2)/(BetaXG[i]^2))+(((BetaYG[i]^2)*(seBetaXG[i]^2))/(BetaXG[i]^4)))
}

CI_lower <- exp(individual_betas-1.96*individual_standard_errors)
CI_upper <- exp(individual_betas+1.96*individual_standard_errors)

individual_results_matrix <- cbind(snp_ids, individual_ORs, CI_lower, CI_upper)
individual_results_dataframe <- as.data.frame(individual_results_matrix, stringsAsFactors = FALSE)
individual_results_dataframe[,-1] <- lapply(individual_results_dataframe[-1], function(x) as.numeric(x))

names(individual_results_dataframe) <- c("snp", "odds_ratio", "CI_lower", "CI_upper")

# That's it, we can view and save the estimates based on each individual SNP
individual_results_dataframe$present <- paste(format(round(individual_results_dataframe$odds_ratio, 2),nsmall=2, trim=TRUE),format(round(individual_results_dataframe$CI_lower, 2), nsmall = 2, trim=TRUE), sep=" (")
individual_results_dataframe$presentable <- paste(individual_results_dataframe$present,format(round(individual_results_dataframe$CI_upper, 2),nsmall = 2, trim=TRUE), sep="-")
individual_results_dataframe$presentably <- paste(individual_results_dataframe$presentable, ")", sep="")
individual_results_dataframe$present <- NULL
individual_results_dataframe$presentable <- NULL
# View(individual_results_dataframe)
write.table(individual_results_dataframe, file=sprintf("%s/results/%s_individual_betas_result.txt", directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
# Now moving on to combining these individual effect estimates using IVW, Egger, and median

    
### Now do Inverse variance weighted method

BetaYG <- mr_dataset$beta_pd
BetaXG <- mr_dataset$beta_exp
seBetaYG <- mr_dataset$se_pd
seBetaXG <- mr_dataset$se_exp

IVWfit <- summary(lm(BetaYG ~ -1+BetaXG,weights=1/seBetaYG^2))

# That's it! This is our point estimate (the standard error will be incorrect):
IVWfit$coef

# Let's correct our standard errors and derive the correct p-value and confidence intervals:
DF <- length(BetaYG)-1
IVWBeta <- IVWfit$coef[1,1]
SE <- IVWfit$coef[1,2]/min(1,IVWfit$sigma)
IVW_p <- 2*(1-pt(abs(IVWBeta/SE),DF))
IVW_CI <- IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE

# And let's merged all of the elements above into one "IVWResults" object:    
IVWResults <- c(IVWBeta,IVW_CI,SE,IVWBeta/SE,IVW_p)

# IVWResults now includes the point estimate, 95% Confidence interval, corrected standard error, t-statistic, p-value

# Let's now exponentiate back to obtain the odds ratio:
round(exp(IVWResults[1]),3)

# Some for 95% upper bound:
round(exp(IVWResults[2]),2)

# And lower bound:
round(exp(IVWResults[3]),2)

#Let's export & save our IVW result. Replace the word exposure here with the actual exposure the get correct file name
IVWResults_OR <- c(exp(IVWBeta),exp(IVW_CI),SE,IVWBeta/SE,IVW_p)
IVWResults_OR_matrix <- matrix(IVWResults_OR, nrow=1, ncol=length(IVWResults))
IVWResults_OR_dataframe <- as.data.frame(IVWResults_OR_matrix)
names(IVWResults_OR_dataframe) <- c("odds_ratio", "CI_lower", "CI_upper", "standard_error", "test_statistic", "p")

IVWResults_OR_dataframe$present <- paste(format(round(IVWResults_OR_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(IVWResults_OR_dataframe$CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
IVWResults_OR_dataframe$presentable <- paste(IVWResults_OR_dataframe$present,format(round(IVWResults_OR_dataframe$CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
IVWResults_OR_dataframe$presentably <- paste(IVWResults_OR_dataframe$presentable, ")", sep="")
IVWResults_OR_dataframe$present <- NULL
IVWResults_OR_dataframe$presentable <- NULL
# View(IVWResults_OR_dataframe)

write.table(IVWResults_OR_dataframe, file=sprintf("%s/results/%s_IVW_results_exposure.txt", directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
### And that's it for the Toby-Johnson method. MR Egger coming soon...


### MR EGGER is here! ###
MREggerFit <- summary(lm(BetaYG ~ BetaXG,weights=1/seBetaYG^2))

# That's it! This is our point estimate (the standard error will be incorrect):
MREggerFit$coef

# Let's correct our standard errors and derive the correct p-value and confidence intervals:
MREggerBeta0 <- MREggerFit$coef[1,1]
MREggerBeta1 <- MREggerFit$coef[2,1]
SE0 <- MREggerFit$coef[1,2]/min(1,MREggerFit$sigma)
SE1 <- MREggerFit$coef[2,2]/min(1,MREggerFit$sigma)
DF <- length(BetaYG)-2
MRBeta0_p <- 2*(1-pt(abs(MREggerBeta0/SE0),DF))
MRBeta1_p <- 2*(1-pt(abs(MREggerBeta1/SE1),DF))
MRBeta0_CI <- MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
MRBeta1_CI <- MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1


# Now let's create a matrix "MREggerResults" with point estimate, 95% Confidence interval, corrected standard error, t-statistic, and p-value for the intercept (AKA pleiotropy) in row 1
MREggerResults <- matrix(nrow = 2,ncol = 6)
MREggerResults[1,] <- c(MREggerBeta0,MRBeta0_CI,SE0,MREggerBeta0/SE0,MRBeta0_p)

# And we'll do the same for slope (AKA effect estimate) in row 2
MREggerResults[2,] <- c(MREggerBeta1,MRBeta1_CI,SE1,MREggerBeta1/SE1,MRBeta1_p)

# Let's now exponentiate back to obtain the odds ratio:
round(exp(MREggerResults[2,][1]),3)

# Same for 95% upper bound:
round(exp(MREggerResults[2,][2]),3)

# And lower bound:
round(exp(MREggerResults[2,][3]),3)

### That's it!

### Alternatively, we can use bootstrapping to obtain the confidence intervals for our effect estimate:
boot <- NULL; straps = 10000
for (i in 1:straps) {
  BYG_boot <- rnorm(length(BetaYG), mean=BetaYG, sd=seBetaYG)
  BXG_boot <- rnorm(length(BetaXG), mean=BetaXG, sd=seBetaXG)
  BYG_boot <- BYG_boot*sign(BXG_boot)
  BXG_boot <- abs(BXG_boot)
  boot[i] <- summary(lm(BYG_boot~BXG_boot,weights=seBetaYG^-2))$coef[2,1]
 }
boot_upper <- sort(boot)[9751]
boot_lower <- sort(boot)[250]
boot_se <- sd(boot)

# And, again, let's merge all of these into one neat object with the beta, and 95% confidence interval, and standard error:
MREggerBoot <- c(MREggerBeta1,boot_lower,boot_upper,boot_se)

# Let's now exponentiate back to obtain the odds ratio:
round(exp(MREggerBoot[1]),3)

# Some for 95% upper bound:
round(exp(MREggerBoot[2]),3)

# And lower bound:
round(exp(MREggerBoot[3]),3)

### All good - add the exponentiation for the MR-Egger analysis


# Let's now do the new and improved MR Egger, which is supposed to behave better (in the recent unpublished paper)
alt_betaEGGER <- summary(lm(BetaYG~BetaXG, weights=seBetaYG^-2))

alt_MREggerBeta0 <- alt_betaEGGER$coef[1,1]
alt_MREggerBeta1 <- alt_betaEGGER$coef[2,1]
alt_SE0 <- alt_betaEGGER$coef[1,2]/min(1,alt_betaEGGER$sigma)
alt_SE1 <- alt_betaEGGER$coef[2,2]/min(1,alt_betaEGGER$sigma)
alt_DF <- length(BetaYG)-2
alt_MRBeta0_p <- 2*(1-pt(abs(MREggerBeta0/SE0),DF))
alt_MRBeta1_p <- 2*(1-pt(abs(MREggerBeta1/SE1),DF))
alt_MRBeta0_CI <- alt_MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
alt_MRBeta1_CI <- alt_MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1

alt_MREggerResults <- matrix(nrow = 2,ncol = 6)
alt_MREggerResults[1,] <- c(alt_MREggerBeta0,alt_MRBeta0_CI,alt_SE0,alt_MREggerBeta0/SE0,alt_MRBeta0_p)
alt_MREggerResults[2,] <- c(alt_MREggerBeta1,alt_MRBeta1_CI,alt_SE1,alt_MREggerBeta1/SE1,alt_MRBeta1_p)

alt_boot <- NULL; straps = 10000
for (i in 1:straps) {
  alt_BYG_boot <- rnorm(length(BetaYG), mean=BetaYG, sd=seBetaYG)
  alt_BXG_boot <- rnorm(length(BetaXG), mean=BetaXG, sd=seBetaXG)
  alt_BYG_boot <- alt_BYG_boot*sign(alt_BXG_boot)
  alt_BXG_boot <- abs(alt_BXG_boot)
  alt_boot[i] <- summary(lm(alt_BYG_boot~alt_BXG_boot,weights=seBetaYG^-2))$coef[2,1]
}
alt_boot_upper <- sort(boot)[9751]
alt_boot_lower <- sort(boot)[250]
alt_boot_se <- sd(boot)

alt_MREggerBoot <- c(alt_MREggerBeta1,alt_boot_lower,alt_boot_upper,alt_boot_se)


#Let's export & save our MRE result. Replace the word exposure here with the actual exposure the get correct file name
old_MREggerResults_OR <- matrix(nrow = 2,ncol = 6)
old_MREggerResults_OR[1,] <- c(MREggerBeta0,MRBeta0_CI,SE0,MREggerBeta0/SE0,MRBeta0_p)
old_MREggerResults_OR[2,] <- c(exp(MREggerBeta1),exp(MRBeta1_CI),SE1,MREggerBeta1/SE1,MRBeta1_p)
old_MREggerBoot_OR <- c(exp(MREggerBeta1),exp(boot_lower),exp(boot_upper),boot_se)

old_MREggerBoot_OR_matrix <- matrix(c(rep(NA,times=3),old_MREggerBoot_OR[-1]), nrow=2, ncol=length(old_MREggerBoot_OR[-1]), byrow=TRUE)
old_MREgger_OR_combined_matrix <- cbind(old_MREggerResults_OR, old_MREggerBoot_OR_matrix)

old_MREggerResults_OR_dataframe <- as.data.frame(old_MREgger_OR_combined_matrix)
names(old_MREggerResults_OR_dataframe) <- c("odds_ratio", "CI_lower", "CI_upper", "standard_error", "test_statistic", "p", "boot_CI_lower", "boot_CI_upper", "boot_se")

alt_MREggerResults_OR <- matrix(nrow = 2,ncol = 6)
alt_MREggerResults_OR[1,] <- c(alt_MREggerBeta0,alt_MRBeta0_CI,alt_SE0,alt_MREggerBeta0/alt_SE0,alt_MRBeta0_p)
alt_MREggerResults_OR[2,] <- c(exp(alt_MREggerBeta1),exp(alt_MRBeta1_CI),alt_SE1,alt_MREggerBeta1/alt_SE1,alt_MRBeta1_p)
alt_MREggerBoot_OR <- c(exp(alt_MREggerBeta1),exp(alt_boot_lower),exp(alt_boot_upper),alt_boot_se)

alt_MREggerBoot_OR_matrix <- matrix(c(rep(NA,times=3),alt_MREggerBoot_OR[-1]), nrow=2, ncol=length(alt_MREggerBoot_OR[-1]), byrow=TRUE)
alt_MREgger_OR_combined_matrix <- cbind(alt_MREggerResults_OR, alt_MREggerBoot_OR_matrix)

alt_MREggerResults_OR_dataframe <- as.data.frame(alt_MREgger_OR_combined_matrix)
names(alt_MREggerResults_OR_dataframe) <- c("new_odds_ratio", "new_CI_lower", "new_CI_upper", "new_standard_error", "new_test_statistic", "new_p", "new_boot_CI_lower", "new_boot_CI_upper", "new_boot_se")

rowtitles <- matrix(c("intercept", "slope"))
rowtitles_dataframe <- as.data.frame(rowtitles)
names(rowtitles_dataframe) <- "parameter"

MREggerResults_OR_dataframe <- cbind(rowtitles_dataframe, old_MREggerResults_OR_dataframe, alt_MREggerResults_OR_dataframe)
MREggerResults_OR_dataframe[,-1] <- lapply(MREggerResults_OR_dataframe[-1], function(x) as.numeric(x))

MREggerResults_OR_dataframe$present <- paste(format(round(MREggerResults_OR_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerResults_OR_dataframe$CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerResults_OR_dataframe$presentable <- paste(MREggerResults_OR_dataframe$present,format(round(MREggerResults_OR_dataframe$CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerResults_OR_dataframe$presentably <- paste(MREggerResults_OR_dataframe$presentable, ")", sep="")
MREggerResults_OR_dataframe$present <- NULL
MREggerResults_OR_dataframe$presentable <- NULL

MREggerResults_OR_dataframe$boot_present <- paste(format(round(MREggerResults_OR_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerResults_OR_dataframe$boot_CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerResults_OR_dataframe$boot_presentable <- paste(MREggerResults_OR_dataframe$boot_present,format(round(MREggerResults_OR_dataframe$boot_CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerResults_OR_dataframe$boot_presentably <- paste(MREggerResults_OR_dataframe$boot_presentable, ")", sep="")
MREggerResults_OR_dataframe$boot_present <- NULL
MREggerResults_OR_dataframe$boot_presentable <- NULL

MREggerResults_OR_dataframe$new_present <- paste(format(round(MREggerResults_OR_dataframe$new_odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerResults_OR_dataframe$new_CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerResults_OR_dataframe$new_presentable <- paste(MREggerResults_OR_dataframe$new_present,format(round(MREggerResults_OR_dataframe$new_CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerResults_OR_dataframe$new_presentably <- paste(MREggerResults_OR_dataframe$new_presentable, ")", sep="")
MREggerResults_OR_dataframe$new_present <- NULL
MREggerResults_OR_dataframe$new_presentable <- NULL

MREggerResults_OR_dataframe$new_boot_present <- paste(format(round(MREggerResults_OR_dataframe$new_odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerResults_OR_dataframe$new_boot_CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerResults_OR_dataframe$new_boot_presentable <- paste(MREggerResults_OR_dataframe$new_boot_present,format(round(MREggerResults_OR_dataframe$new_boot_CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerResults_OR_dataframe$new_boot_presentably <- paste(MREggerResults_OR_dataframe$new_boot_presentable, ")", sep="")
MREggerResults_OR_dataframe$new_boot_present <- NULL
MREggerResults_OR_dataframe$new_boot_presentable <- NULL

# View(MREggerResults_OR_dataframe)

write.table(MREggerResults_OR_dataframe, file=sprintf("%s/results/%s_MRE_results_exposure.txt", directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)


### MEDIAN-BASED METHOD ###

betaYG <- mr_dataset$beta_pd
betaXG <- mr_dataset$beta_exp
sebetaYG <- mr_dataset$se_pd
sebetaXG <- mr_dataset$se_exp

weighted_median <- function(betaIV_in, weights_in) {
  betaIV_order <- betaIV_in[order(betaIV_in)]
  weights_order <- weights_in[order(betaIV_in)]
  weights_sum <- cumsum(weights_order)-0.5*weights_order
  weights_sum <- weights_sum/sum(weights_order)
  below <- max(which(weights_sum<0.5))
  weighted_est <- betaIV_order[below] + (betaIV_order[below+1]-betaIV_order[below])*
    (0.5-weights_sum[below])/(weights_sum[below+1]-weights_sum[below])
  return(weighted_est) }

# Bootstrap
weighted_median_boot <- function(betaXG_in, betaYG_in, sebetaXG_in, sebetaYG_in, weights_in){
  med <- NULL
  for(i in 1:1000){
    betaXG_boot <- rnorm(length(betaXG_in), mean=betaXG_in, sd=sebetaXG_in)
    betaYG_boot <- rnorm(length(betaYG_in), mean=betaYG_in, sd=sebetaYG_in)
    betaIV_boot <- betaYG_boot/betaXG_boot
    med[i] <- weighted_median(betaIV_boot, weights_in)
  }
  return(sd(med)) }


betaIV <- betaYG/betaXG # ratio estimates
weights <- (sebetaYG/betaXG)^-2 # inverse-variance weights
betaIVW <- sum(betaYG*betaXG*sebetaYG^-2)/sum(betaXG^2*sebetaYG^-2) # IVW estimate

# penalized weights:
penalty <- pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
pen_weights <- weights*pmin(1, penalty*20)

# weighted median estimate:
betaWM <- weighted_median(betaIV, weights)

# weighted median standard error:
sebetaWM <- weighted_median_boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)

# penalized weighted median estimate:
betaPWM <- weighted_median(betaIV, pen_weights)

# penalized weighted median standard error:
sebetaPWM <- weighted_median_boot(betaXG, betaYG, sebetaXG, sebetaYG, pen_weights)

# Alternative weighted median method
weighted_median_empirical <- function(betaIV_in, weights_in) {
  betaIV_order <- betaIV_in[order(betaIV_in)]
  weights_order <- weights_in[order(betaIV_in)]
  weights_sum <- cumsum(weights_order)
  weights_sum <- weights_sum/sum(weights_order)
  which_below <- max(which(weights_sum<0.5))
  return(betaIV_order[which_below+1]) }

# alternative weighted median estimate
betaWME <- weighted_median_empirical(betaIV, weights)

# Can use inverse standard errors instead of inverse variances as weights with the following:
# weights <- (sebetaYG/betaXG)^-1

# Let's collect these results together:

weighted_median_odds_ratio <- exp(betaWM)
weighted_median_CI_lower <- exp(betaWM-1.96*sebetaWM)
weighted_median_CI_upper <- exp(betaWM+1.96*sebetaWM)

penalised_weighted_median_odds_ratio <- exp(betaPWM)
penalised_weighted_median_CI_lower <- exp(betaPWM-1.96*sebetaPWM)
penalised_weighted_median_CI_upper <- exp(betaPWM+1.96*sebetaPWM)

empirical_weighted_median_odds_ratio <- exp(betaWME)


#Let's export & save our MRE result. Replace the word exposure here with the actual exposure the get correct file name
weighted_median_OR_results_matrix_nonames <- matrix(c(weighted_median_odds_ratio, weighted_median_CI_lower, weighted_median_CI_upper, penalised_weighted_median_odds_ratio, penalised_weighted_median_CI_lower, penalised_weighted_median_CI_upper, empirical_weighted_median_odds_ratio, NA, NA), ncol = 3, byrow = TRUE)
weighted_median_OR_results_matrix <- cbind(c("weighted_median", "penalised_weighted_median", "empirical_weighted_median"), weighted_median_OR_results_matrix_nonames)

weighted_median_OR_results_dataframe <- as.data.frame(weighted_median_OR_results_matrix, stringsAsFactors = FALSE)
names(weighted_median_OR_results_dataframe) <- c("method", "odds_ratio", "CI_lower", "CI_upper")
weighted_median_OR_results_dataframe[,-1] <- lapply(weighted_median_OR_results_dataframe[-1], function(x) as.numeric(x))

# View(weighted_median_OR_results_dataframe)

weighted_median_OR_results_dataframe$present <- paste(format(round(weighted_median_OR_results_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(weighted_median_OR_results_dataframe$CI_lower, 2), nsmall=2,trim=TRUE), sep=" (")
weighted_median_OR_results_dataframe$presentable <- paste(weighted_median_OR_results_dataframe$present,format(round(weighted_median_OR_results_dataframe$CI_upper, 2),nsmall = 2,trim=TRUE), sep="-")
weighted_median_OR_results_dataframe$presentably <- paste(weighted_median_OR_results_dataframe$presentable, ")", sep="")
weighted_median_OR_results_dataframe$present <- NULL
weighted_median_OR_results_dataframe$presentable <- NULL

# View(weighted_median_OR_results_dataframe)

write.table(weighted_median_OR_results_dataframe, file=sprintf("%s/results/%s_weighted_median_results_exposure.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)

### Now let's do some tests for heterogeneity...
# First, we'll obtain individual effect estimates for each SNP (again, using urate for demonstration purposes):
BetaYG <- mr_dataset$beta_pd
BetaXG <- mr_dataset$beta_exp
seBetaYG <- mr_dataset$se_pd
seBetaXG <- mr_dataset$se_exp

individual_betas <- vector(length=n_snps)
for (i in 1:n_snps){
  individual_betas[i] <- BetaYG[i]/BetaXG[i]
}

# Now, we need the variances, using the delta method:
individual_variances <- vector(length=n_snps)

for (i in 1:n_snps){
  individual_variances[i] <- ((seBetaYG[i]^2)/(BetaXG[i]^2))+(((BetaYG[i]^2)*(seBetaXG[i]^2))/(BetaXG[i]^4))
}

# Now we need the inverse variances to use as weights for Cochran's Q:
inverse_variances <- 1/individual_variances

# Now we can make the calculation for Cochran's Q
weighted_mean <- sum(inverse_variances*individual_betas)/sum(inverse_variances)
Q <- sum(inverse_variances*((individual_betas-weighted_mean)^2))

# And to get the p-value for Cochran's Q (if significant, evidence that there is pleiotropy):
pvalue_pleiotropy <- pchisq(Q,n_snps,lower.tail=FALSE)
print(pvalue_pleiotropy)

# And the I2 statistic:
if (Q >= n_snps-1){
  I2 <- ((Q-(n_snps-1))/Q)*100
} else {
  I2 <- 0
}

#I2 95% confidence interval upper and bounds:
H <- sqrt(Q/(n_snps-1))

if (Q > n_snps){
  se_ln_H <- (1/2)*((log(Q-log(n_snps-1)))/(sqrt(2*Q)-sqrt(2*n_snps-3)))
} else {
  se_ln_H <- sqrt((1/(2*(n_snps-2)))*(1-(1/3*(n_snps-2)^2)))
}

I2_CI_upper <- 100*(1-(1/(exp(log(H)+1.96*se_ln_H))^2))
I2_CI_lower <- 100*(1-(1/(exp(log(H)-1.96*se_ln_H))^2))

print(I2)
print(I2_CI_lower)
print(I2_CI_upper)

#Let's now save and export these results
heterogeneity_results <- c(Q, pvalue_pleiotropy, I2, I2_CI_lower, I2_CI_upper)
heterogeneity_results_matrix <- matrix(heterogeneity_results, nrow=1, ncol=length(heterogeneity_results))
heterogeneity_results_dataframe <- as.data.frame(heterogeneity_results_matrix)
names(heterogeneity_results_dataframe) <- c("cochran's_q", "p", "I2", "I2_CI_lower", "I2_CI_upper")

# View(heterogeneity_results_dataframe)
write.table(heterogeneity_results_dataframe, file=sprintf("%s/qc/%s_heterogeneity_results_exposure.txt", directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)

### NOW WE CAN PLOT OUR RESULTS ###

# Let's first create a forest plot of our estimates
### NOW REPEAT WITH X-AXIS IN LOGARITHMIC SCALE ###
library(scales)

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)


plotting_data <- data.frame(snps = factor(c(individual_results_dataframe$snp, "IVW", "MR_Egger", "MR_Egger_boot", "new_MR_Egger", "new_MR_Egger_boot", "weighted_median", "penalised_weighted_median", "empirical_Weighted_median", "G"), levels=c(rev(c(individual_results_dataframe$snp, "IVW", "MR_Egger", "MR_Egger_boot", "new_MR_Egger", "new_MR_Egger_boot", "weighted_median", "penalised_weighted_median", "empirical_Weighted_median")), "G")),
                            odds_ratios = c(individual_results_dataframe$odds_ratio, IVWResults_OR_dataframe$odds_ratio, rep(MREggerResults_OR_dataframe$odds_ratio[2], times=2), rep(MREggerResults_OR_dataframe$new_odds_ratio[2], times=2), weighted_median_OR_results_dataframe$odds_ratio, NA),
                            lower = c(individual_results_dataframe$CI_lower, IVWResults_OR_dataframe$CI_lower, MREggerResults_OR_dataframe$CI_lower[2], MREggerResults_OR_dataframe$boot_CI_lower[2], MREggerResults_OR_dataframe$new_CI_lower[2], MREggerResults_OR_dataframe$new_boot_CI_lower[2], weighted_median_OR_results_dataframe$CI_lower, NA),
                            upper = c(individual_results_dataframe$CI_upper, IVWResults_OR_dataframe$CI_upper, MREggerResults_OR_dataframe$CI_upper[2], MREggerResults_OR_dataframe$boot_CI_upper[2], MREggerResults_OR_dataframe$new_CI_upper[2], MREggerResults_OR_dataframe$new_boot_CI_upper[2], weighted_median_OR_results_dataframe$CI_upper, NA))

maxrange <- round(max(plotting_data$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange <- round(min(plotting_data$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals <- c(rev(seq(1, minrange, -round(((maxrange)-(minrange))/10, 1))), seq(1, maxrange, round(((maxrange)-(minrange))/10, 1)))

p <- ggplot(plotting_data, aes(odds_ratios, snps))
p <- p + geom_point(size=1, shape=19)
p <- p + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p <- p + geom_vline(xintercept = 1, linetype = "longdash")
p <- p + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma)
p <- p + labs(x="Odds ratio (95% CI)", y="")
p <- p + labs(y=NULL)
lab <- data.frame(V0 = factor(c(individual_results_dataframe$snp, "IVW", "MR_Egger", "MR_Egger_boot", "new_MR_Egger", "new_MR_Egger_boot", "weighted_median", "penalised_weighted_median", "empirical_Weighted_median", "G", individual_results_dataframe$snp, "IVW", "MR_Egger", "MR_Egger_boot", "new_MR_Egger", "new_MR_Egger_boot", "weighted_median", "penalised_weighted_median", "empirical_Weighted_median", "G"),levels = rev(c(individual_results_dataframe$snp, "IVW", "MR_Egger", "MR_Egger_boot", "new_MR_Egger", "new_MR_Egger_boot", "weighted_median", "penalised_weighted_median", "empirical_Weighted_median", "G"))),
                  V05 = rep(c(1,2), each=nrow(plotting_data)),
                  V1 = c("Instrumental variable", individual_results_dataframe$snp, "IVW", "MR Egger", "MR Egger Bootstrapping", "New MR Egger", "New MR Egger Bootstrapping", "Weighted Median", "Penalised Weighted Median", "Empirical Weighted Median", "OR (95% CI)", individual_results_dataframe$presentably, IVWResults_OR_dataframe$presentably, MREggerResults_OR_dataframe$presentably[2], MREggerResults_OR_dataframe$boot_presentably[2], MREggerResults_OR_dataframe$new_presentably[2], MREggerResults_OR_dataframe$new_boot_presentably[2], weighted_median_OR_results_dataframe$presentably[c(1,2)],format(round(weighted_median_OR_results_dataframe$odds_ratio[3],2),nsmall=2,trim=TRUE)))

data_table <- ggplot(lab, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table <- data_table + geom_text(size = 3.5, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table <- data_table + geom_hline(aes(yintercept=nrow(plotting_data)-0.5))
data_table <- data_table + geom_hline(aes(yintercept=nrow(plotting_data)+0.5))
data_table <- data_table + theme(panel.grid.major = element_blank(),
                                 legend.position = "none",
                                 panel.border = element_blank(), 
                                 axis.text.x = element_text(colour="white"),#element_blank(),
                                 axis.text.y = element_blank(), 
                                 axis.ticks = element_line(colour="white"),#element_blank(),
                                 plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))
data_table
grid.arrange(data_table, p, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_results <- arrangeGrob(data_table, p, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_results_exposure_full_logscale.pdf", directory, exposure), forest_plot_results)


# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_results <- arrangeGrob(data_table, p, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_results_exposure_full_logscale_modified.pdf", directory, exposure), forest_plot_results)


### THAT'S IT (ggplot2 gives us pretty much unparalleled power to tweak and modify our figures until we're happy with the result)!

### And now we can plot a scatter plot of the betas with lines for our estimates:
BetaYG <- mr_dataset$beta_pd
BetaXG <- mr_dataset$beta_exp
seBetaYG <- mr_dataset$se_pd
seBetaXG <- mr_dataset$se_exp

theme_set(theme_bw())
theme_update(legend.position = "bottom")
plot <- ggplot(mr_dataset, aes(x=BetaXG, y=BetaYG))
plot <- plot + geom_point(size=3, shape=19)
plot <- plot + geom_errorbar(aes(ymin = BetaYG-1.96*seBetaYG,ymax = BetaYG+1.96*seBetaYG, width = 0.001))
plot <- plot + geom_errorbarh(aes(xmin = BetaXG-1.96*seBetaXG ,xmax = BetaXG+1.96*seBetaXG, height = 0.001))
plot <- plot + geom_abline(aes(intercept=0, slope=IVWResults[1], colour="line1", show.legend = TRUE))
plot <- plot + geom_abline(aes(intercept=MREggerResults[1,1], slope=MREggerResults[2,1], colour="line2", show.legend = TRUE))
plot <- plot + geom_abline(aes(intercept=0, slope=betaPWM, colour="line3", show.legend = TRUE))
plot <- plot + scale_colour_manual(values=c(line1="red", line2="blue", line3="green"), labels = c("Inverse variance weighted", "Egger", "Penalised weighted median"))
plot <- plot + xlab("Parkinson's disease beta")
plot <- plot + ylab(sprintf("%s beta",scatter_y_label))

plot

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
ggsave(file=sprintf("%s/figures/%s_scatter_plot_results_exposure.pdf", directory, exposure), plot)

### Next step will be to run some simulations to screen for pleiotropy, possibly multivariable MR, and some plots. I'll add some scripts for this tomorrow.

### Let's repeat the analysis now, each time leaving one SNP out (may have bugs, need to go thourgh this to make sure works). This script loops around, running both IVW and Egger method on the dataset, each time leaving one SNP out. The end result is three tables:
# i) IVWResults_oneout, which has 7 columns: "missing snp", "beta", "standard error", "CI lower bound", "CI upper bound", "test statistic", "p-value"
# ii) MREggerResults_oneout, which has 13 columns: "missing snp", "intercept point estimate", "standard error of intercept", "CI lower bound (intercept)", "CI upper bound (intercept)", "test statistic (intercept)", "p-value (intecept)","effect point estimate", "standard error of effect", "CI lower bound (effect)", "CI upper bound (effect)", "test statistic (effect)", "p-value (effect)"
# iii) MREggerBoot_oneout, which has 5 columns: "missing snp","beta","bootstrap standard error","bootstrap CI lower bound","bootstrap CI upper bound"

### I'll extend this into a simulation where each time a random sample of three SNPs is removed (this method works particularly well when the instrument is made of many SNPs, such as BMI)
IVWResults_oneout <- matrix(nrow=n_snps, ncol=7)
MREggerResults_oneout <- matrix(nrow=n_snps ,ncol = 13)
MREggerBoot_oneout <- matrix(nrow=n_snps, ncol=5)

IVWResults_OR_oneout <- matrix(nrow=n_snps, ncol=7)
MREggerResults_OR_oneout <- matrix(nrow=n_snps ,ncol = 13)
MREggerBoot_OR_oneout <- matrix(nrow=n_snps, ncol=5)
weighted_median_OR_oneout <- matrix(nrow=n_snps, ncol=4)
penalised_weighted_median_OR_oneout <- matrix(nrow=n_snps, ncol=4)
empirical_weighted_median_OR_oneout <- matrix(nrow=n_snps, ncol=2)

for (i in 1:n_snps){
  mr_dataset_oneout <- mr_dataset[-i,]
  
  BetaYG_oneout <- mr_dataset_oneout$beta_pd
  BetaXG_oneout <- mr_dataset_oneout$beta_exp
  seBetaYG_oneout <- mr_dataset_oneout$se_pd
  seBetaXG_oneout <- mr_dataset_oneout$se_exp
  
  IVWfit_oneout <- summary(lm(BetaYG_oneout ~ -1+BetaXG_oneout,weights=1/seBetaYG_oneout^2))
  
  DF_oneout <- length(BetaYG_oneout)-1
  IVWBeta_oneout <- IVWfit_oneout$coef[1,1]
  SE_oneout <- IVWfit_oneout$coef[1,2]/min(1,IVWfit_oneout$sigma)
  IVW_p_oneout <- 2*(1-pt(abs(IVWBeta_oneout/SE_oneout),DF_oneout))
  IVW_CI_oneout <- IVWBeta_oneout + c(-1,1)*qt(df=DF_oneout, 0.975)*SE_oneout
  
  IVWResults_oneout[i,] <- c(snp_ids[i],IVWBeta_oneout,IVW_CI_oneout,SE,IVWBeta_oneout/SE_oneout,IVW_p_oneout)
  IVWResults_OR_oneout[i,] <- c(snp_ids[i],exp(IVWBeta_oneout),exp(IVW_CI_oneout),SE,IVWBeta_oneout/SE_oneout,IVW_p_oneout)  
  
  ### And same for Egger...  
  MREggerFit_oneout <- summary(lm(BetaYG_oneout ~ BetaXG_oneout,weights=1/seBetaYG_oneout^2))
  
  MREggerBeta0_oneout <- MREggerFit_oneout$coef[1,1]
  MREggerBeta1_oneout <- MREggerFit_oneout$coef[2,1]
  SE0_oneout <- MREggerFit_oneout$coef[1,2]/min(1,MREggerFit_oneout$sigma)
  SE1_oneout <- MREggerFit_oneout$coef[2,2]/min(1,MREggerFit_oneout$sigma)
  DF_oneout <- length(BetaYG_oneout)-2
  MRBeta0_p_oneout <- 2*(1-pt(abs(MREggerBeta0_oneout/SE0_oneout),DF_oneout))
  MRBeta1_p_oneout <- 2*(1-pt(abs(MREggerBeta1_oneout/SE1_oneout),DF_oneout))
  MRBeta0_CI_oneout <- MREggerBeta0_oneout + c(-1,1)*qt(df=DF_oneout, 0.975)*SE0_oneout
  MRBeta1_CI_oneout <- MREggerBeta1_oneout + c(-1,1)*qt(df=DF_oneout, 0.975)*SE1_oneout
  
  MREggerResults_oneout[i,] <- c(snp_ids[i],MREggerBeta0_oneout,SE0_oneout,MRBeta0_CI_oneout,MREggerBeta0_oneout/SE0_oneout,MRBeta0_p_oneout,MREggerBeta1_oneout,SE1_oneout,MRBeta1_CI_oneout,MREggerBeta1_oneout/SE1_oneout,MRBeta1_p_oneout)
  MREggerResults_OR_oneout[i,] <- c(snp_ids[i],MREggerBeta0_oneout,MRBeta0_CI_oneout,SE0_oneout,MREggerBeta0_oneout/SE0_oneout,MRBeta0_p_oneout,exp(MREggerBeta1_oneout),exp(MRBeta1_CI_oneout),SE1_oneout,MREggerBeta1_oneout/SE1_oneout,MRBeta1_p_oneout)
  
  boot_oneout <- NULL; straps = 10000
  for (j in 1:straps) {
    BYG_boot_oneout <- rnorm(length(BetaYG_oneout), mean=BetaYG_oneout, sd=seBetaYG_oneout)
    BXG_boot_oneout <- rnorm(length(BetaXG_oneout), mean=BetaXG_oneout, sd=seBetaXG_oneout)
    BYG_boot_oneout <- BYG_boot_oneout*sign(BXG_boot_oneout)
    BXG_boot_oneout <- abs(BXG_boot_oneout)
    boot_oneout[j] <- summary(lm(BYG_boot_oneout~BXG_boot_oneout,weights=seBetaYG_oneout^-2))$coef[2,1]
  }
  boot_upper_oneout <- sort(boot_oneout)[9751]
  boot_lower_oneout <- sort(boot_oneout)[250]
  boot_se_oneout <- sd(boot_oneout)
  
  MREggerBoot_oneout[i,] <- c(snp_ids[i],MREggerBeta1_oneout,boot_se_oneout,boot_lower_oneout,boot_upper_oneout)
  MREggerBoot_OR_oneout[i,] <- c(snp_ids[i],exp(MREggerBeta1_oneout),exp(boot_lower_oneout),exp(boot_upper_oneout),boot_se_oneout)

  ### MEDIAN-BASED METHOD ###
  
  BYG_median_oneout <- mr_dataset_oneout$beta_pd
  BXG_median_oneout <- mr_dataset_oneout$beta_exp
  seBYG_median_oneout <- mr_dataset_oneout$se_pd
  seBXG_median_oneout <- mr_dataset_oneout$se_exp
  
  weighted_median <- function(betaIV_in, weights_in) {
    betaIV_order <- betaIV_in[order(betaIV_in)]
    weights_order <- weights_in[order(betaIV_in)]
    weights_sum <- cumsum(weights_order)-0.5*weights_order
    weights_sum <- weights_sum/sum(weights_order)
    below <- max(which(weights_sum<0.5))
    weighted_est <- betaIV_order[below] + (betaIV_order[below+1]-betaIV_order[below])*
      (0.5-weights_sum[below])/(weights_sum[below+1]-weights_sum[below])
    return(weighted_est) }
  
  # Bootstrap
  weighted_median_boot <- function(betaXG_in, betaYG_in, sebetaXG_in, sebetaYG_in, weights_in){
    med <- NULL
    for(i in 1:1000){
      betaXG_boot <- rnorm(length(betaXG_in), mean=betaXG_in, sd=sebetaXG_in)
      betaYG_boot <- rnorm(length(betaYG_in), mean=betaYG_in, sd=sebetaYG_in)
      betaIV_boot <- betaYG_boot/betaXG_boot
      med[i] <- weighted_median(betaIV_boot, weights_in)
    }
    return(sd(med)) }
  
  
  betaIV_oneout <- BYG_median_oneout/BXG_median_oneout # ratio estimates
  weights_oneout <- (seBYG_median_oneout/BXG_median_oneout)^-2 # inverse-variance weights
  betaIVW_oneout <- sum(BYG_median_oneout*BXG_median_oneout*seBYG_median_oneout^-2)/sum(BXG_median_oneout^2*seBYG_median_oneout^-2) # IVW estimate
  
  # penalized weights:
  penalty_oneout <- pchisq(weights*(betaIV_oneout-betaIVW_oneout)^2, df=1, lower.tail=FALSE)
  pen_weights_oneout <- weights_oneout*pmin(1, penalty*20)
  
  # weighted median estimate:
  betaWM_oneout <- weighted_median(betaIV_oneout, weights_oneout)
  
  # weighted median standard error:
  sebetaWM_oneout <- weighted_median_boot(BXG_median_oneout, BYG_median_oneout, seBXG_median_oneout, seBYG_median_oneout, weights_oneout)
  
  # penalized weighted median estimate:
  betaPWM_oneout <- weighted_median(betaIV_oneout, pen_weights_oneout)
  
  # penalized weighted median standard error:
  sebetaPWM_oneout <- weighted_median_boot(BXG_median_oneout, BYG_median_oneout, seBXG_median_oneout, seBYG_median_oneout, pen_weights)
  
  # Alternative weighted median method
  weighted_median_empirical <- function(betaIV_in, weights_in) {
    betaIV_order <- betaIV_in[order(betaIV_in)]
    weights_order <- weights_in[order(betaIV_in)]
    weights_sum <- cumsum(weights_order)
    weights_sum <- weights_sum/sum(weights_order)
    which_below <- max(which(weights_sum<0.5))
    return(betaIV_order[which_below+1]) }
  
  # alternative weighted median estimate
  betaWME_oneout <- weighted_median_empirical(betaIV_oneout, weights_oneout)
  
  # Can use inverse standard errors instead of inverse variances as weights with the following:
  # weights <- (sebetaYG/betaXG)^-1
  
  # Let's collect these results together:
  
  weighted_median_odds_ratio_oneout <- exp(betaWM_oneout)
  weighted_median_CI_lower_oneout <- exp(betaWM_oneout-1.96*sebetaWM_oneout)
  weighted_median_CI_upper_oneout <- exp(betaWM_oneout+1.96*sebetaWM_oneout)
  
  penalised_weighted_median_odds_ratio_oneout <- exp(betaPWM_oneout)
  penalised_weighted_median_CI_lower_oneout <- exp(betaPWM_oneout-1.96*sebetaPWM_oneout)
  penalised_weighted_median_CI_upper_oneout <- exp(betaPWM_oneout+1.96*sebetaPWM_oneout)
  
  empirical_weighted_median_odds_ratio_oneout <- exp(betaWME_oneout)
  
  
  weighted_median_OR_oneout[i,] <- c(snp_ids[i], weighted_median_odds_ratio_oneout, weighted_median_CI_lower_oneout, weighted_median_CI_upper_oneout)
  penalised_weighted_median_OR_oneout[i,] <- c(snp_ids[i], penalised_weighted_median_odds_ratio_oneout, penalised_weighted_median_CI_lower_oneout, penalised_weighted_median_CI_upper_oneout)
  empirical_weighted_median_OR_oneout[i,] <- c(snp_ids[i], empirical_weighted_median_odds_ratio_oneout)
}
  
### Now, let's convert these matrices to data frames, to make it easier to view and plot the results etc
IVWResults_OR_oneout_dataframe <- as.data.frame(IVWResults_OR_oneout, stringsAsFactors=FALSE)
MREggerResults_OR_oneout_dataframe <- as.data.frame(MREggerResults_OR_oneout, stringsAsFactors=FALSE)
MREggerBoot_OR_oneout_dataframe <- as.data.frame(MREggerBoot_OR_oneout, stringsAsFactors=FALSE)
weighted_median_OR_oneout_dataframe <- as.data.frame(weighted_median_OR_oneout, stringsAsFactors=FALSE)
penalised_weighted_median_OR_oneout_dataframe <- as.data.frame(penalised_weighted_median_OR_oneout, stringsAsFactors=FALSE)
empirical_weighted_median_OR_oneout_dataframe <- as.data.frame(empirical_weighted_median_OR_oneout, stringsAsFactors=FALSE)

# Next we coerce all columns except the snp name column into type numeric
IVWResults_OR_oneout_dataframe[,-1] <- lapply(IVWResults_OR_oneout_dataframe[-1], function(x) as.numeric(x))
MREggerResults_OR_oneout_dataframe[,-1] <- lapply(MREggerResults_OR_oneout_dataframe[-1], function(x) as.numeric(x))
MREggerBoot_OR_oneout_dataframe[,-1] <- lapply(MREggerBoot_OR_oneout_dataframe[-1], function(x) as.numeric(x))
weighted_median_OR_oneout_dataframe[,-1] <- lapply(weighted_median_OR_oneout_dataframe[-1], function(x) as.numeric(x))
penalised_weighted_median_OR_oneout_dataframe[,-1] <- lapply(penalised_weighted_median_OR_oneout_dataframe[-1], function(x) as.numeric(x))
empirical_weighted_median_OR_oneout_dataframe[,-1] <- lapply(empirical_weighted_median_OR_oneout_dataframe[-1], function(x) as.numeric(x))

# We add column names
names(IVWResults_OR_oneout_dataframe) <- c("missing_snp", "odds_ratio", "CI_lower_bound", "CI_upper_bound", "standard_error", "test_statistic", "p-value") 
names(MREggerResults_OR_oneout_dataframe) <- c("missing_snp", "intercept_point_estimate", "CI_lower_bound_(intercept)", "CI_upper_bound_(intercept)", "standard_error_of_intercept", "test_statistic_(intercept)", "p-value_(intecept)","effect_odds_ratio", "CI_lower_bound_(effect)", "CI_upper_bound_(effect)", "standard_error_of_effect", "test_statistic_(effect)", "p-value_(effect)")
names(MREggerBoot_OR_oneout_dataframe) <- c("missing_snp","odds_ratio","bootstrap_CI_lower_bound","bootstrap_CI_upper_bound","bootstrap_standard_error")
names(weighted_median_OR_oneout_dataframe) <- c("missing_snp", "odds_ratio", "CI_lower_bound", "CI_upper_bound") 
names(penalised_weighted_median_OR_oneout_dataframe) <- c("missing_snp", "odds_ratio", "CI_lower_bound", "CI_upper_bound") 
names(empirical_weighted_median_OR_oneout_dataframe) <- c("missing_snp", "odds_ratio") 


# Let's add a column with ORs and CIs

IVWResults_OR_oneout_dataframe$present <- paste(format(round(IVWResults_OR_oneout_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(IVWResults_OR_oneout_dataframe$CI_lower_bound, 2), nsmall=2,trim=TRUE), sep=" (")
IVWResults_OR_oneout_dataframe$presentable <- paste(IVWResults_OR_oneout_dataframe$present,format(round(IVWResults_OR_oneout_dataframe$CI_upper_bound, 2),nsmall = 2,trim=TRUE), sep="-")
IVWResults_OR_oneout_dataframe$presentably <- paste(IVWResults_OR_oneout_dataframe$presentable, ")", sep="")
IVWResults_OR_oneout_dataframe$present <- NULL
IVWResults_OR_oneout_dataframe$presentable <- NULL

MREggerResults_OR_oneout_dataframe$present <- paste(format(round(MREggerResults_OR_oneout_dataframe$effect_odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerResults_OR_oneout_dataframe$`CI_lower_bound_(effect)`, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerResults_OR_oneout_dataframe$presentable <- paste(MREggerResults_OR_oneout_dataframe$present,format(round(MREggerResults_OR_oneout_dataframe$`CI_upper_bound_(effect)`, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerResults_OR_oneout_dataframe$presentably <- paste(MREggerResults_OR_oneout_dataframe$presentable, ")", sep="")
MREggerResults_OR_oneout_dataframe$present <- NULL
MREggerResults_OR_oneout_dataframe$presentable <- NULL

MREggerBoot_OR_oneout_dataframe$present <- paste(format(round(MREggerBoot_OR_oneout_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(MREggerBoot_OR_oneout_dataframe$bootstrap_CI_lower_bound, 2), nsmall=2,trim=TRUE), sep=" (")
MREggerBoot_OR_oneout_dataframe$presentable <- paste(MREggerBoot_OR_oneout_dataframe$present,format(round(MREggerBoot_OR_oneout_dataframe$bootstrap_CI_upper_bound, 2),nsmall = 2,trim=TRUE), sep="-")
MREggerBoot_OR_oneout_dataframe$presentably <- paste(MREggerBoot_OR_oneout_dataframe$presentable, ")", sep="")
MREggerBoot_OR_oneout_dataframe$present <- NULL
MREggerBoot_OR_oneout_dataframe$presentable <- NULL

weighted_median_OR_oneout_dataframe$present <- paste(format(round(weighted_median_OR_oneout_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(weighted_median_OR_oneout_dataframe$CI_lower_bound, 2), nsmall=2,trim=TRUE), sep=" (")
weighted_median_OR_oneout_dataframe$presentable <- paste(weighted_median_OR_oneout_dataframe$present,format(round(weighted_median_OR_oneout_dataframe$CI_upper_bound, 2),nsmall = 2,trim=TRUE), sep="-")
weighted_median_OR_oneout_dataframe$presentably <- paste(weighted_median_OR_oneout_dataframe$presentable, ")", sep="")
weighted_median_OR_oneout_dataframe$present <- NULL
weighted_median_OR_oneout_dataframe$presentable <- NULL

penalised_weighted_median_OR_oneout_dataframe$present <- paste(format(round(penalised_weighted_median_OR_oneout_dataframe$odds_ratio, 2),nsmall=2,trim=TRUE),format(round(penalised_weighted_median_OR_oneout_dataframe$CI_lower_bound, 2), nsmall=2,trim=TRUE), sep=" (")
penalised_weighted_median_OR_oneout_dataframe$presentable <- paste(penalised_weighted_median_OR_oneout_dataframe$present,format(round(penalised_weighted_median_OR_oneout_dataframe$CI_upper_bound, 2),nsmall = 2,trim=TRUE), sep="-")
penalised_weighted_median_OR_oneout_dataframe$presentably <- paste(penalised_weighted_median_OR_oneout_dataframe$presentable, ")", sep="")
penalised_weighted_median_OR_oneout_dataframe$present <- NULL
penalised_weighted_median_OR_oneout_dataframe$presentable <- NULL

### Now we can view our results nicely
# View(IVWResults_OR_oneout_dataframe)
# View(MREggerResults_OR_oneout_dataframe)
# View(MREggerBoot_OR_oneout_dataframe)
# View(weighted_median_OR_oneout_dataframe)
# View(penalised_weighted_median_OR_oneout_dataframe)
# View(empirical_weighted_median_OR_oneout_dataframe)

# Let's export and save these tables, so they can be viewed in excel etc
write.table(IVWResults_OR_oneout_dataframe, file=sprintf("%s/results/%s_IVW_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
write.table(MREggerResults_OR_oneout_dataframe, file=sprintf("%s/results/%s_EGGER_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
write.table(MREggerBoot_OR_oneout_dataframe, file=sprintf("%s/results/%s_BootEGGER_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
write.table(weighted_median_OR_oneout_dataframe, file=sprintf("%s/results/%s_weighted_median_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
write.table(penalised_weighted_median_OR_oneout_dataframe, file=sprintf("%s/results/%s_penalised_weighted_median_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)
write.table(empirical_weighted_median_OR_oneout_dataframe, file=sprintf("%s/results/%s_empirical_weighted_median_oneout_results.txt",directory, exposure), quote=FALSE, sep="\t", row.names=FALSE)



### NOW WE CAN DRAW A FOREST PLOT OF THE "LEAVE-ONE-OUT" RESULTS"

oneout_snp_names <- IVWResults_OR_oneout_dataframe$missing_snp
oneout_odds_ratios <- IVWResults_OR_oneout_dataframe$odds_ratio
oneout_lower_bounds <- IVWResults_OR_oneout_dataframe$CI_lower_bound
oneout_upper_bounds <- IVWResults_OR_oneout_dataframe$CI_upper_bound
oneout_presentable_odds_ratios <- IVWResults_OR_oneout_dataframe$presentably

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
                                   odds_ratios = c(oneout_odds_ratios, NA),
                                   lower = c(oneout_lower_bounds, NA),
                                   upper = c(oneout_upper_bounds, NA))
number_of_rows <- nrow(plotting_data_oneout)

maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))

p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
p_oneout <- p_oneout + geom_point(size=1, shape=19)
p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
p_oneout <- p_oneout + labs(y=NULL)

lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
                  V05 = rep(c(1,2), each=number_of_rows),
                  V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))

data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
                                 legend.position = "none",
                                 panel.border = element_blank(), 
                                 axis.text.x = element_text(colour="white"),#element_blank(),
                                 axis.text.y = element_blank(), 
                                 axis.ticks = element_line(colour="white"),#element_blank(),
                                 plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))

grid.arrange(data_table_oneout, p_oneout, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_oneout_results_exposure.pdf", directory, exposure), forest_plot_oneout_results)


#### SAME WITH EGGER ###

oneout_snp_names <- MREggerResults_OR_oneout_dataframe$missing_snp
oneout_odds_ratios <- MREggerResults_OR_oneout_dataframe$effect_odds_ratio
oneout_lower_bounds <- MREggerResults_OR_oneout_dataframe$`CI_lower_bound_(effect)`
oneout_upper_bounds <- MREggerResults_OR_oneout_dataframe$`CI_upper_bound_(effect)`
oneout_presentable_odds_ratios <- MREggerResults_OR_oneout_dataframe$presentably

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
                                   odds_ratios = c(oneout_odds_ratios, NA),
                                   lower = c(oneout_lower_bounds, NA),
                                   upper = c(oneout_upper_bounds, NA))
number_of_rows <- nrow(plotting_data_oneout)

maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))

p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
p_oneout <- p_oneout + geom_point(size=1, shape=19)
p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
p_oneout <- p_oneout + labs(y=NULL)

lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
                         V05 = rep(c(1,2), each=number_of_rows),
                         V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))

data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
                                               legend.position = "none",
                                               panel.border = element_blank(), 
                                               axis.text.x = element_text(colour="white"),#element_blank(),
                                               axis.text.y = element_blank(), 
                                               axis.ticks = element_line(colour="white"),#element_blank(),
                                               plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))

grid.arrange(data_table_oneout, p_oneout, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_oneout_Egger_results_exposure.pdf", directory, exposure), forest_plot_oneout_results)


### SAME WITH BOOT EGGER ###

oneout_snp_names <- MREggerBoot_OR_oneout_dataframe$missing_snp
oneout_odds_ratios <- MREggerBoot_OR_oneout_dataframe$odds_ratio
oneout_lower_bounds <- MREggerBoot_OR_oneout_dataframe$bootstrap_CI_lower_bound
oneout_upper_bounds <- MREggerBoot_OR_oneout_dataframe$bootstrap_CI_upper_bound
oneout_presentable_odds_ratios <- MREggerBoot_OR_oneout_dataframe$presentably

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
                                   odds_ratios = c(oneout_odds_ratios, NA),
                                   lower = c(oneout_lower_bounds, NA),
                                   upper = c(oneout_upper_bounds, NA))
number_of_rows <- nrow(plotting_data_oneout)

maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))

p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
p_oneout <- p_oneout + geom_point(size=1, shape=19)
p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
p_oneout <- p_oneout + labs(y=NULL)

lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
                         V05 = rep(c(1,2), each=number_of_rows),
                         V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))

data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
                                               legend.position = "none",
                                               panel.border = element_blank(), 
                                               axis.text.x = element_text(colour="white"),#element_blank(),
                                               axis.text.y = element_blank(), 
                                               axis.ticks = element_line(colour="white"),#element_blank(),
                                               plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))

grid.arrange(data_table_oneout, p_oneout, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_oneout_Boot_Egger_results_exposure.pdf", directory, exposure), forest_plot_oneout_results)


### SAME WITH WEIGHTED MEDIAN

oneout_snp_names <- weighted_median_OR_oneout_dataframe$missing_snp
oneout_odds_ratios <- weighted_median_OR_oneout_dataframe$odds_ratio
oneout_lower_bounds <- weighted_median_OR_oneout_dataframe$CI_lower_bound
oneout_upper_bounds <- weighted_median_OR_oneout_dataframe$CI_upper_bound
oneout_presentable_odds_ratios <- weighted_median_OR_oneout_dataframe$presentably

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
                                   odds_ratios = c(oneout_odds_ratios, NA),
                                   lower = c(oneout_lower_bounds, NA),
                                   upper = c(oneout_upper_bounds, NA))
number_of_rows <- nrow(plotting_data_oneout)

maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))

p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
p_oneout <- p_oneout + geom_point(size=1, shape=19)
p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
p_oneout <- p_oneout + labs(y=NULL)

lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
                         V05 = rep(c(1,2), each=number_of_rows),
                         V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))

data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
                                               legend.position = "none",
                                               panel.border = element_blank(), 
                                               axis.text.x = element_text(colour="white"),#element_blank(),
                                               axis.text.y = element_blank(), 
                                               axis.ticks = element_line(colour="white"),#element_blank(),
                                               plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))

grid.arrange(data_table_oneout, p_oneout, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_weighted_median_oneout_results.pdf", directory, exposure), forest_plot_oneout_results)


### SAME WITH PENALISED WEIGHTED MEDIAN

oneout_snp_names <- penalised_weighted_median_OR_oneout_dataframe$missing_snp
oneout_odds_ratios <- penalised_weighted_median_OR_oneout_dataframe$odds_ratio
oneout_lower_bounds <- penalised_weighted_median_OR_oneout_dataframe$CI_lower_bound
oneout_upper_bounds <- penalised_weighted_median_OR_oneout_dataframe$CI_upper_bound
oneout_presentable_odds_ratios <- penalised_weighted_median_OR_oneout_dataframe$presentably

theme_set(theme_bw())
theme_update(
  axis.line.x = element_line(colour = "red", size=1),
  axis.line.y = element_line(colour = "red", size=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
                                   odds_ratios = c(oneout_odds_ratios, NA),
                                   lower = c(oneout_lower_bounds, NA),
                                   upper = c(oneout_upper_bounds, NA))
number_of_rows <- nrow(plotting_data_oneout)

maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))

p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
p_oneout <- p_oneout + geom_point(size=1, shape=19)
p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
p_oneout <- p_oneout + labs(y=NULL)

lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
                         V05 = rep(c(1,2), each=number_of_rows),
                         V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))

data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
                                               legend.position = "none",
                                               panel.border = element_blank(), 
                                               axis.text.x = element_text(colour="white"),#element_blank(),
                                               axis.text.y = element_blank(), 
                                               axis.ticks = element_line(colour="white"),#element_blank(),
                                               plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))

grid.arrange(data_table_oneout, p_oneout, ncol=2)

# Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
ggsave(file=sprintf("%s/figures/%s_forest_plot_penalised_weighted_median_oneout_results.pdf", directory, exposure), forest_plot_oneout_results)


# ### SAME WITH EMPIRICAL WEIGHTED MEDIAN
# 
# oneout_snp_names <- empirical_weighted_median_OR_oneout_dataframe$missing_snp
# oneout_odds_ratios <- empirical_weighted_median_OR_oneout_dataframe$odds_ratio
# oneout_lower_bounds <- empirical_weighted_median_OR_oneout_dataframe$CI_lower_bound
# oneout_upper_bounds <- empirical_weighted_median_OR_oneout_dataframe$CI_upper_bound
# oneout_presentable_odds_ratios <- empirical_weighted_median_OR_oneout_dataframe$presentably
# 
# theme_set(theme_bw())
# theme_update(
#   axis.line.x = element_line(colour = "red", size=1),
#   axis.line.y = element_line(colour = "red", size=1),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(),
#   panel.border = element_blank(),
#   panel.background = element_blank(),
#   axis.text.y = element_blank(),
#   axis.ticks.y = element_blank(),
#   plot.margin = unit(c(0,0,0,0), "lines")
# )
# 
# plotting_data_oneout <- data.frame(snps = factor(c(oneout_snp_names, "G"), levels=c(rev(oneout_snp_names), "G")),
#                                    odds_ratios = c(oneout_odds_ratios, NA),
#                                    lower = c(oneout_lower_bounds, NA),
#                                    upper = c(oneout_upper_bounds, NA))
# number_of_rows <- nrow(plotting_data_oneout)
# 
# maxrange_oneout <- round(max(plotting_data_oneout$upper, na.rm=TRUE)/0.2)*0.2+0.2
# minrange_oneout <- round(min(plotting_data_oneout$lower, na.rm=TRUE)/0.2)*0.2-0.2
# intervals_oneout <- c(rev(seq(1, minrange_oneout, -round(((maxrange_oneout)-(minrange_oneout))/10, 1))), seq(1, maxrange_oneout, round(((maxrange_oneout)-(minrange_oneout))/10, 1)))
# 
# p_oneout <- ggplot(plotting_data_oneout, aes(odds_ratios, snps))
# p_oneout <- p_oneout + geom_point(size=1, shape=19)
# p_oneout <- p_oneout + geom_errorbarh(aes(xmax = upper, xmin = lower, height = 0.15))
# p_oneout <- p_oneout + geom_vline(xintercept = 1, linetype = "longdash")
# p_oneout <- p_oneout + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = comma_format(digits=2))
# p_oneout <- p_oneout + labs(x="Odds Ratio (95% CI)", y="")
# p_oneout <- p_oneout + labs(y=NULL)
# 
# lab_oneout <- data.frame(V0 = factor(c(oneout_snp_names, "G", oneout_snp_names, "G"),levels = rev(c(oneout_snp_names,"G"))),
#                          V05 = rep(c(1,2), each=number_of_rows),
#                          V1 = c("Missing SNP", oneout_snp_names, "OR (95% CI)", oneout_presentable_odds_ratios))
# 
# data_table_oneout <- ggplot(lab_oneout, aes(x = V05, y = V0, label = format(V1, nsmall = 1)))
# data_table_oneout <- data_table_oneout + geom_text(size = 4, hjust=0, vjust=0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
# data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows-0.5))
# data_table_oneout <- data_table_oneout + geom_hline(aes(yintercept=number_of_rows+0.5))
# data_table_oneout <- data_table_oneout + theme(panel.grid.major = element_blank(),
#                                                legend.position = "none",
#                                                panel.border = element_blank(), 
#                                                axis.text.x = element_text(colour="white"),#element_blank(),
#                                                axis.text.y = element_blank(), 
#                                                axis.ticks = element_line(colour="white"),#element_blank(),
#                                                plot.margin = unit(c(0,-0.5,0,0), "lines")) + labs(x="",y="") + coord_cartesian(xlim=c(1,2.5))
# 
# grid.arrange(data_table_oneout, p_oneout, ncol=2)
# 
# # Now we save this image as a pdf (remember to replace the word "exposure" in the filename with the actual exposure)
# forest_plot_oneout_results <- arrangeGrob(data_table_oneout, p_oneout, ncol=2)
# ggsave(file=sprintf("%s/figures/%s_forest_plot_empirical_weighted_median_oneout_results.pdf", directory, exposure), forest_plot_oneout_results)
# 
# 


# ### NOW, LET'S MOVE ON TO MULTIVARIABLE MR ###
# # First, we need to construct our mvmr_dataset (as we did in the beginning of this script). However, this time add in all the SNPs that are associated with any one of the exposures (a total of n SNPs)
# # In addition to the "main" risk factor column (urate in our example), we need additional identical columns for each additional risk factor. Note, that this requires that we have access to the betas and standard errors of all the SNPs in our dataset for all of the risk factors (this may prove problematic, if data other than the "top hits" are not available publicly)
# # First, we need to create the following vectors:
# # i) MVBetaYG (as before), is a length n vector of gene-outcome betas for n SNPs
# # ii) MVseBetaYG (as before) is a length n vector of standard errors for the above gene-outcome betas of n SNPs
# # iii) MVBetaXG1, MVBetaXG1, MVBetaXG2, MVBetaXG3 ... MVBetaXGm are m vectors of length n, with gene-exposure betas for exposures XG1, XG2, XG3...XGm, for all n SNPs
# # iv) MVseBetaXG1, MVBetaXG1, MVBetaXG2, MVBetaXG3 ... MVBetaXGm are m vectors of length n, with standard errors of the above gene-exposure betas for exposures XG1, XG2, XG3...XGm, for all n SNPs
# 
# MVBetaYG <- mvmr_dataset$beta_pd
# MVBetaX1 <- mvmr_dataset$beta_riskfactor1
# MVBetaX2 <- mvmr_dataset$beta_riskfactor2
# MVBetaX3 <- mvmr_dataset$beta_riskfactor3
# # Add for every additional risk factor
# 
# MVseBetaYG <- mvmr_dataset$standard_error_pd
# MVseBetaX1 <- mvmr_dataset$standard_error_riskfactor1
# MVseBetaX2 <- mvmr_dataset$standard_error_riskfactor2
# MVseBetaX3 <- mvmr_dataset$standard_error_riskfactor3
# # Add for every additional risk factor
# 
# # Now, run regress the gene-outcome betas on the "interfering" risk-factors, take the residuals of this model, and regress the residuals on the risk-factor of interest. This gives us the "corrected" effect estimate for that risk-factor
# effect_estimate_riskfactor1 <- lm(lm(MVBetaYG~MVBetaX2+MVBetaX3)$res~MVBetaX1)$coef[2]
# effect_estimate_riskfactor2 <- lm(lm(MVBetaYG~MVBetaX1+MVBetaX3)$res~MVBetaX2)$coef[2]
# effect_estimate_riskfactor3 <- lm(lm(MVBetaYG~MVBetaX1+MVBetaX2)$res~MVBetaX3)$coef[2]
# 
# 
# # And we pull out the standard errors for these effect estimates:
# standard_error_riskfactor1 <- summary(lm(lm(MVBetaYG~MVBetaX2+MVBetaX3)$res~MVBetaX1))$coef[2,2]
# standard_error_riskfactor2 <- summary(lm(lm(MVBetaYG~MVBetaX1+MVBetaX3)$res~MVBetaX2))$coef[2,2]
# standard_error_riskfactor3 <- summary(lm(lm(MVBetaYG~MVBetaX1+MVBetaX2)$res~MVBetaX3))$coef[2,2]

### NOTE: While this approach works, it does not take into consideration the standard errors of our gene-outcome or gene-exposure betas. Indeed, as described in the original paper, "it is an ad hoc approach which has no clear theoretical basis and which ignores the uncertainty in the beta coefficients". Therefore, we can use it and it gives us a rough idea of the effects, but if we decide that we want to do multivariable MR, we should verify these results with the Bayesian likelihood-based method, implemented in something like WinBUGS (I'll prepare the script for this later)... ###