#########################################################################################################################
# This script normalizes gene expression data from micrroarrays
# IMPORTANT!: - when dataset is 'StudioFormat' = 'yes', means it is in BeadStudio or GenomeStudio output file format
#              (illuminia proprietary softwares for analyzing data output by the scanning system), and we need to:
#               * manually gunzip the non-normalized data folder and extract the 'non_normalized' data file (name varies)
#               * change file name to <Id>_non-normalized_data.txt and save outside of untared folder, remove the zipped file
#               * manually change file header to be 'ProbeID', sample.AVG_Signal', 'sample.Detection Pval'
#                 when not in correct format
#########################################################################################################################

### Load Packages
library(GEOquery)
library(oligo)
library(limma)
source('<home>/<usr>/bRNAsPipe/scr/R/MicroarrayNormalize_Functions.R')

### User defined settings
dataDir <- '<home>/<usr>/bRNAsPipe/data' # directory with data for a project
modelFile <- '<home>/<usr>/bRNAsPipe/MetModels/iHuman_processedInfoEMEM.xlsx' # path to xlsx file with model info
studySettingsTableName <- 'Studies_Microarrays.xlsx'
studyId <- 'M8.2.3'

### Read study info
StudyInfo = ReadStudyInfo(studyId, dataDir, studySettingsTableName)

### Get pheno (sample) and feature (gene) metadata:
PhenoFeature = getPhenoFeature(StudyInfo$studyDir, StudyInfo$FeatFileName, StudyInfo$GPLId, modelFile, StudyInfo$Id)

### Preprocess raw expression data provided at GEO:
rawData = PreProcessRawExpr(StudyInfo$studyDir, StudyInfo$Id, PhenoFeature$microarrayType, StudyInfo$sampleCorr, StudyInfo$illumnStudioFormat)

### Normalization of expression data
NormExpr = NormaliseExpr(PhenoFeature$microarrayType, rawData, StudyInfo$Id)

### Join normalized expression data, sample (pheno) and feature data:
JoinedExpPhenoFeat = JoinExpPhenoFeat(StudyInfo$samplesList, PhenoFeature$pheno, NormExpr$normData, PhenoFeature$feature, StudyInfo$treatments, PhenoFeature$InitialFinalIds, StudyInfo$Id, NormExpr$ToKeep, PhenoFeature$microarrayType, modelFile)

### Do a summarization to gene level of transcript expression values - for sample and sample groups:
SummNormExpr = summarizationFunct(JoinedExpPhenoFeat$FinalIdsExprVal, StudyInfo$organism) # for each sample
GroupNormExpr = GroupNormExp(studySettingsTableName, StudyInfo$studyId, SummNormExpr) # for each sample group
 # with R version 4 gives a warning - it is normal

### save gene normalized transcript expression data, phenoData (sample data) and featureData (probe data):
SaveNormExpr(StudyInfo$studyDir, SummNormExpr, GroupNormExpr, JoinedExpPhenoFeat$eset, StudyInfo$Id) # saves normalized transcript expression for each sample


