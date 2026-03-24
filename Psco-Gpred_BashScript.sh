# Mikaeylah Davidson, 2025
# Scripts used for genomic predictions and evaluations in Pseudophryne corroboree 
# all analysis was either run using
    # 1. This script, on The Univeristy of Melbournes HPC 'SPARTAN' or 
    # 2. R script, using R on local drive
# All data and scripts is stored within the R project associated with this manuscript "Psco-Gpred.R"


# SNP data was originally called and filtered, in the manuscript: 
    # Design and application of a genome-wide SNP array to improve conservation outcomes in the Critically Endangered southern corroboree frog, Davidson et al. 
    # Data is available at: DataDryad https://doi.org/10.5061/dryad.ngf1vhj4b
    # This step takes the previously filtered data, and filters further for the purposes of this analysis.


###--------------- SNP FILE ---------------###
# Original SNP file from PSCO array work, includes controls and non-F1 samples 
# This was subsetted in that work to remove 
            # Any samples not from Bd exposure expeirment (972)
            # Controls (54)
            # 23 aniamls which failed QC
# Resulting in
    # 895 out of 1087 Individuals
    # 39681 out of a possible 39681 Sites  

    # Copy plink files from that work for use
cp plink-infectionPsco . 


#--- UPDATE .FAM FILE
    # add population data to file 
        # created file from phenotypes file which is 
        # clutch, id 
        #        **NOTE** when creating IDs CANT be only alpha, must be alphanumeric or numeric
        #           + make sure pheno file is in same order as geno file 
    nano update-plink-fam.txt
    # need to add  0, 0, 0, -9 for formatting
    #          sireID, motherid, sexCode, phenotype
    awk '{print $1, $2, 0, 0, 0, -9}' update-plink-fam.txt > plink-infectionPsco.fam


###--------------- HWE & MISSINGNESS FILTERING ---------------###
    # Get HWE pvalues & filter 
    mkdir filtering
    # get HWE pvalues 
    plink --bfile ../plink-infectionPsco --hardy 

    # filter by p = 10^-6
    # missingness = 95% CR 
    plink --bfile ../plink-infectionPsco --hwe 0.0000001 --geno 0.05 --recode --out hwe-0.0000001
        # --hwe: 1597 variants removed 
        # --geno: 316 variants removed
        # 37768 variants and 895 

###--- REFORMAT FILTERED FILE
    # copy and move HWE filtered data into main dir 
    cp hwe-0.0000001.map ../../filteredPsco.37768.map
    cp hwe-0.0000001.ped ../../filteredPsco.37768.ped

    # generate bed file 
    plink --file filteredPsco.37768 --make-bed --noweb --out filteredPsco.37768
    # generate vcf 
        # need to use vcf-iid to keep ID as ONLY psco_#, if use just vcf then adds pop to start of ID
    plink --bfile filteredPsco.37768   --recode vcf-iid --out filteredPsco.37768


###--------------- PHENOTYPE FILE ---------------###
# Upload phenotype file "phenotypes-895.txt" from local drive


###--- FILTERED & FORMATTED FILE
# Genotype data 
filteredPsco.37768
    # 37768 SNPs
    # 895 animals 
# Phenotype data 
phenotypes-895.txt


###---------------   GWAS   ---------------###
###--------------- SIGNIFICANCE THRESHOLD ---------------###
    # Calculate significance thresholds for GWAS on unlinked SNPs
    # Unlinked = r2 < 0.2 
    # Prune SNP file based on r2 of 0.2
    plink --bfile ../filteredPsco.37768 --indep-pairwise 50 5 0.2 --out LD_0.2
    #  31973 of 37768 variants removed.
    #  5795 unlinked SNPs


###--------------- GWAS SIGNIFICANT SNPS---------------###
    mkdir gwas-asreml
    #  snp  files, so they are in 012 format
    # convert to plink file 
    vcftools --vcf ../filteredPsco.37768.vcf --plink --out reml
    # convert to bed 
    plink --file reml --make-bed --noweb --keep-allele-order --out reml
    # reformat into 012 format 

# Get ALT and Ref SNPs for GWAS significant SNP output 
    # Make list of sgniifcant SNPs 
    snp-list 
    # pull the SNPs with their alt allele FIRST, then reference allele 
    # note == this will make them into alphabetical order 
    awk 'NR==FNR {snp[$1]; next} $3 in snp {print $3, $5, $4}' snp-list ../filteredPsco.37768.vcf

    # NOTE -- if SNP id uses a "-" not "_" in snp list file then use this
    awk 'NR==FNR {gsub("_", "-", $1); snp[$1]; next} $3 in snp {print $3, $5, $4}' snp-list ../filteredPsco.37768.vcf

# Plot phenotypes against genotypes for significant snps 
    # make list of significant snps 
    nano snp-significant
    # Subset the bed file with only the significant SNPs
    # convert bed to 0/1/2 format 
    plink --bfile ../filteredPsco.37768 --extract snp-list --recodeA --out significant-SNPs


###---------------  EBVS ---------------###
###--------------- GBLUP ---------------###
    # Make GBLUP to preduct the genetic value (Breeding values, GEBVs)
    # make BLUP of the 
    #   1. individuals for each trait 
    #       col $4 is the total genetic variance = EBV
    #   2. the SNPs (SNP effects) for each trait

    #  Download all to local drive
###--- Survival
    # For individuals 
    gcta64 --reml --grm $GRM/gcta --pheno ../heritabilities/days-survived/pheno_daysSurvived --reml-pred-rand --out daysSurvived
        # Columns = 
        # family ID, individual ID, an intermediate variable, the total genetic value, another intermediate variable and the residual. 
    # Zoo + Age
    gcta64 --reml --grm $GRM/gcta --pheno ../heritabilities/days-survived/pheno_daysSurvived --covar ../heritabilities/covar-zoo --qcovar ../heritabilities/covar-age  --reml-pred-rand --out daysSurvived+AgeZoo

    # SNP effects
    gcta64 --bfile ../filteredPsco.37768 --blup-snp daysSurvived.indi.blp --out daysSurvived
        # Columns =
        # SNP ID, reference allele, BLUP of SNP effect, and residual effect

###--- Pathogen load
    # For individuals 
    gcta64 --reml --grm $GRM/gcta --pheno ../heritabilities/bd-slope/pheno_bdSlope --reml-pred-rand --out bdSlope
      # SNP effects
    gcta64 --bfile ../filteredPsco.37768 --blup-snp bdSlope.indi.blp --out bdSlope

###--- Body Condition
    # For individuals 
    gcta64 --reml --grm $GRM/gcta --pheno ../heritabilities/body-condition/pheno_bc --reml-pred-rand --out bc
      # SNP effects
    gcta64 --bfile ../filteredPsco.37768 --blup-snp bc.indi.blp --out bc
 


###--------------- SELECTIVE BREEDING MODELS ---------------###

###--------------- ALPHAMATE ---------------###
    # Windows application == need to run on vitual machine 
    # Need to create:  
    # 1. CoancestryMatrixFile
        # GRM without the top row of IDs in txt format
    # 2. GenderFile
        # we dont know genders, so used dummy 1,2 replicated for all animals 
    # 3. SelCriterionFile 
        # ID followed by the EBV for the trait slecting for 

    # Create AlphaMateSpec.txt 
        # this is the specification file to change selection parameters 

###--------------- BY ZOO ---------------###
    # split the animals up by zoo, and run selection seperatly 
    
    # Created keep lists for each zoo in excel 
    nano keep-MZ
    nano keep-TZ
    # Subset vcf to make new bed file for each zoo
    vcftools --vcf ../filteredPsco.37768.vcf --keep keep-MZ --recode-INFO-all --recode --out MZ
    vcftools --vcf ../filteredPsco.37768.vcf --keep keep-TZ --recode-INFO-all --recode --out TZ
    # convert to .bed to be able to make GRM in R
    vcftools --vcf MZ.recode.vcf --plink --out MZ
    plink --file MZ --make-bed --noweb --out MZ
    vcftools --vcf TZ.recode.vcf --plink --out TZ
    plink --file TZ --make-bed --noweb --out TZ
    # download files into R





