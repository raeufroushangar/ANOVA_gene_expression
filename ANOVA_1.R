library(data.table)



# Analysis of Variance:
#----------------------


# read ANOVA archive for annotation
setwd( '/Volumes/RaMain/GEOData/AML_Analyses/')

# this data archive has 2761 rows (arrays/subjects) and 5 columns in the following order: 
# ID_REF,	Age,	Sex,	Sample_source,	Disease_state

# use file to annotate the 44754 probesets from the batch effects corrected gene expression data
array_annotation =read.csv("AML_and_Healthy_archive_with_predected_sex_and_sample_source_for_ANOVA_without_Study_geo_accession.csv", header = TRUE, sep = ",",as.is=TRUE)

# age binning used to group subjects/patients age into 7 different age-groups
# 0 = 0 to 19 
# 1 = 20 to 29 
# 2 = 30 to 39
# 3 = 40 to 49
# 4 = 50 to 59 
# 5 = 60 to 69 
# 6 = 70 to 79 
# 7 = 80 and above
array_annotation$Age=as.factor(findInterval(array_annotation$Age, c(20, 30, 40, 50, 60, 70, 80)))
array_annotation$Age


setwd('/Volumes/RaMain/GEOData/AML_Analyses/Each_DataSet_Corrected_Datasetwise/Correct_for_All_Factors_SampleSource_DiseaseState_Batch/')


# using data.table package since file is too large
eset=as.data.frame(fread(file="All_2761_Corrected_for_All_Factors_SampleSource_DiseaseState_Batch_Datasetwise_2213_AML_1st_548_Healthy_2nd_and_removed_613_dummy_with_44754_probsets_RMA_Normalized_Log2Trans_Zscore_Standardized_Transposed_Data.csv", header=TRUE, sep=",",colClasses = "numeric"))


# assign probset-id to row name
dim(eset)
rownames(eset) = eset[,1]
eset[,1] = NULL
rownames(eset)


# EXPR is dependent variable (gene expression probeset)
# array_annotation$EXPR is used to clear before next iteration
array_annotation$EXPR
geneNames=row.names(eset[1])
for (EXPR in 1:nrow(eset)) {  
  array_annotation$EXPR
  array_annotation$EXPR = as.numeric(c(eset[EXPR,]))
  # linear model
  aovres = aov(EXPR~ Age + Sex + Sample_source + Disease_state + 
                     Age:Sex + Age:Sample_source + Age:Disease_state + 
                     Sex:Sample_source + Sex:Disease_state + 
                     Sample_source:Disease_state, 
                     data = array_annotation )
  AVz = anova(aovres)
  # write p-values and probesets names into csv file
  write(c(geneNames[EXPR], AVz$`Pr(>F)`), "ANOVA_P_Values_For_All_2761_Corrected_for_All_Factors_SampleSource_DiseaseState_Batch_Datasetwise_2213_AML_1st_548_Healthy_2nd_and_removed_613_dummy_with_44754_probsets_RMA_Normalized_Log2Trans_Zscore_Standardized_Transposed_Data.csv",sep=",", append=TRUE,ncolumns=11)
  }




# Adjusted P-Values of Variance:
#------------------------------

# add headers to the output file from ANOVA: 
# headers: ProbeSet, Age, Sex, Sample_source, Disease_state, Age:Sex, Age:Sample_source, Age:Disease_state, Sex:Sample_source, Sex:Disease_state, Residuals

# set directory of ANOVA p-value file  
setwd('/Volumes/RaMain/GEOData/AML_Analyses/Each_DataSet_Corrected_Datasetwise/Correct_for_All_Factors_SampleSource_DiseaseState_Batch/')


read_pvalues_file1 =read.csv("ANOVA_P_Values_For_All_2761_Corrected_for_All_Factors_SampleSource_DiseaseState_Batch_Datasetwise_2213_AML_1st_548_Healthy_2nd_and_removed_613_dummy_with_44754_probsets_RMA_Normalized_Log2Trans_Zscore_Standardized_Transposed_Data.csv", header = TRUE, sep = ",",as.is=TRUE)

BH_p_adjus_for_disease_state1 = p.adjust(read_pvalues_file1$Disease_state, "BH")<0.01
length((read_pvalues_file1$ProbeSet[BH_p_adjus_for_disease_state1]))

bonferroni_p_adjus_for_disease_state1 = p.adjust(read_pvalues_file1$Disease_state, "bonferroni")<0.01
length((read_pvalues_file1$ProbeSet[bonferroni_p_adjus_for_disease_state1]))
bonferroni_p_adjus_for_disease_state_probsets1 = read_pvalues_file1$ProbeSet[bonferroni_p_adjus_for_disease_state1]


write.csv(bonferroni_p_adjus_for_disease_state_probsets1, file="Bonferroni_P_Adjus_Disease_state_15694_probsets_from_ANOVA_P_Values_For_All_2761_Corrected_for_All_Factors_SampleSource_DiseaseState_Batch_Datasetwise.csv")




#-------------------#-------------------#-------------------#-------------------#-------------------

# set directory of ANOVA p-value file  
setwd('/Volumes/RaMain/GEOData/AML_Analyses/All_DataSets_Corrected_Together/')



read_pvalues_file2 =read.csv("ANOVA_P_Values_For_All_2761_Corrected_for_All_Factors_2213_AML_1st_548_Healthy_2nd_and_removed_613_dummy_with_44754_probsets_RMA_Normalized_Log2Trans_Zscore_Standardized_Transposed_Data.csv", header = TRUE, sep = ",",as.is=TRUE)

# after ANOVA, we filtered output for signinfigant genes using benjamini<0.01 on disease state.
BH_p_adjus_for_disease_state2 = p.adjust(read_pvalues_file2$Disease_state, "BH")<0.01
length((read_pvalues_file2$ProbeSet[BH_p_adjus_for_disease_state2]))

# after ANOVA, we filtered output for signinfigant genes using bonferroni<0.01 on disease state.
bonferroni_p_adjus_for_disease_state2= p.adjust(read_pvalues_file2$Disease_state, "bonferroni")<0.01
length((read_pvalues_file2$ProbeSet[bonferroni_p_adjus_for_disease_state2]))
bonferroni_p_adjus_for_disease_state_probsets2 = read_pvalues_file2$ProbeSet[bonferroni_p_adjus_for_disease_state2]


write.csv(bonferroni_p_adjus_for_disease_state_probsets2, file="Bonferroni_P_Adjus_for_Disease_state_16292_probsets_from_ANOVA_P_Values_For_All_2761_Corrected_for_All_Factors_not_datasetwise.csv")






