
library(tictoc)
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


eset[1,1]
array_annotation$EXPR
array_annotation$EXPR = as.numeric(c(expr_data[1,]))
array_annotation$EXPR = as.numeric(c(eset[4,]))
names(array_annotation$EXPR)

tic() # Starts timer 
aovres = aov(EXPR~ Study_geo_accession + Age + Sex + Sample_source + Disease_state + Study_geo_accession:Age + Study_geo_accession:Sex + Study_geo_accession:Sample_source + Study_geo_accession:Disease_state + Age:Sex + Age:Sample_source + Age:Disease_state + Sex:Sample_source + Sex:Disease_state + Sample_source:Disease_state, data = array_annotation )

aovres = aov(EXPR~ Age + Sex + Sample_source + Disease_state + Age:Sex + Age:Sample_source + Age:Disease_state + Sex:Sample_source + Sex:Disease_state + Sample_source:Disease_state, data = array_annotation )

AVz = anova(aovres)
toc() # End timer 
AVz
AVz$`Pr(>F)`
final_result = c(row.names(eset)[[1]],AVz$`Pr(>F)`)
write.table(t(final_result), "ANOVA_P_Values copy.csv",sep=",", append=TRUE, row.names = FALSE,col.names=FALSE) 

row.names(AVz$`Pr(>F)`)
array_annotation
write.table(AVz$`Pr(>F)` , file="44754_ProbSets_ANOVA_P_Values_After_BECorrection.csv", sep=",", row.names = FALSE)


AVz$`Pr(>F)`
tk = TukeyHSD(aovres)
(tk$`Age:Disease_state`)


1:nrow(expr_data) #find number of columns in file  
AVz = rep(NA, nrow(expr_data)) # Create a table, AVz with the same number of rows as in expr_data
AVz
for (i in 2:nrow(expr_data)) {
  each_row = row.names(expr_data[i])
  array_annotation$EXPR=expr_data[1,]
aovres=aov(EXPR~ Study_geo_accession + Age + Sex + Sample_source + Disease_state + Study_geo_accession:Age + Study_geo_accession:Sex + Study_geo_accession:Sample_source + Study_geo_accession:Disease_state + Age:Sex + Age:Sample_source + Age:Disease_state + Sex:Sample_source + Sex:Disease_state + Sample_source:Disease_state, data = array_annotation )
AVz = anova(aovres)
tk = TukeyHSD(aovres)



expr_data[1,1:2] # 1st row and 1st & 2ns columns
expr_data[1,] # 1st row across all columns

my_table = read.csv(file="AML_and_Healthy_archive_with_predected_sex_and_sample_source_for_ANOVA.csv", header=TRUE, sep=",")
my_table[1:6,][["ID_REF"]]

summary (aml_model) #this will give you your p values
TukeyHSD(aml_model, conf.level = 0.99)




