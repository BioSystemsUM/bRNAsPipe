#################
# Load Packages #
#################

library(GEOquery)
library(oligo)
library(limma)
library(beadarray)
library(openxlsx)
library(ArrayExpress)
source("D:/Tania_2/PhD/CSCs/scr/R/ConvertIds.R")

#############
# Functions #
#############

ReadStudyInfo <- function(studyId, dataDir, studySettingsTableName){
  studyDir <- paste0(dataDir,'/',studyId) # path of directory of specific study
  dir.create(studyDir, suppressWarnings('already exists'))
  studySettingsTable <- openxlsx::read.xlsx(paste0(dataDir,'/', studySettingsTableName), sheet=1, colNames=T)
  study = subset(studySettingsTable, studySettingsTable$Study==studyId)
  organism <- unique(study$Organism)
  Id <- unique(study$DatabaseStudyId) # Id of the study/series
  GPLId <- unique(study$GPL_Id) # Id of microarray platform
  FeatFileName = study$ArrayExpressFeatureData[1] # when feature data is a supplfile (ArrayExpress), this is name of file containing feature data - has to be manually downloaded
  #get list in which each element is a treatment:
  treatments <- unique(study$Conditions)
  # get list in which each elemnt is sample ids of a certain treatment:
  samplesList <- list()
  for (trt in treatments) {samplesList[[trt]] <- study[study$Conditions==trt, 'SampleId'] }
  ## sometimes samples are not identified by their ids
  ## and we need to get correspondence between original sample name and sample id:
  if (!is.na(study$Sample_Name_Correspondence[1])) { # i there is info on correspondence between sample name and sample id
    sampleCorr <- list() # list which elements will have sample id, and initial sample names as names
    corresp <- unlist(strsplit(study$Sample_Name_Correspondence, ', '))
    for (el in corresp) {
      el_split <- unlist(strsplit(el, '='))
      Initial <- el_split[1] # initial sample name
      Final <- el_split[2] # final sample id
      sampleCorr[[Initial]] <- Final # give initial sample name as name and sample id as value of each element of list
    }
  } else {
    sampleCorr <- list()
  }
  illumnStudioFormat <- unique(study$StudioFormat) # If dataset is from illumina, is it in BeadStudio or GenomeStudio output file format (illuminia proprietary softwares for analyzing data output by the scanning system)
  return(list(studyId = studyId, studyDir = studyDir, organism = organism, Id = Id, GPLId = GPLId, treatments = treatments, samplesList = samplesList, sampleCorr = sampleCorr, illumnStudioFormat = illumnStudioFormat, FeatFileName = FeatFileName))
}


# getPhenoFeature - gets pheno (sample) and feature (gene) metadata from GEO
#                 - adds ensembl and entrez id and symbol to each probe set id on feature data
#                 - When there is more than one entrez id for same probe set id, the entrez id that is also used for cell metabolic model is chosen.
#                 - If any of those optional entrez ids are not in model it chooses one arbitrarily
#                 - Note! use this function when metabolic model we'll work with has entrez gene ids - cause this uses functions that translate better from microarray to entrez
#                         Also this is optimized for recon3d01, as that model has info on gene splice variant instead of gene, which is processed here
#         - Arguments: * Id: string with id of study/series (for e.g. GEO id, eg.'GSE22247')
#                      * GPLId: string with GEO microarray platform identifier, eg.'GPL6947'
#                      * modelFile: path to file with info on base metabolic model.
#                                   has to have a tab 'GENES' with mRNA variant in entrez id format
#         - Values: * pheno: AnnotatedDataframe with pheno data (sample)
#                   * feature: AnnotatedDataframe with feature data (probe set)
#                   * InitialFinalIds: dataframe with initial probe set ids in 1st column and entrez gene ids in 2nd column
#                   * microarrayType: sentence describing microarray type

getPhenoFeature <- function (studyDir, FeatFileName, GPLId, modelFile, Id) {
  if (grepl('E-', Id)) { # if study from ArrayExpress
    ## Get feature (genes) data:
    rawDataDir <- paste0(studyDir, '/', Id) # path to directory where raw data will be downloaded
    dir.create(rawDataDir, suppressWarnings('already exists')) # creates directory where raw data will be downloaded to
    feature <- read.table(paste0(rawDataDir, '/', FeatFileName), sep = '\t', stringsAsFactors = F, header = T, fill=T)
    ids <- as.character(feature[,1]) # probe ids
    platformData <- getGEO(GPLId, GSEMatrix = T) # get microarray platform info from GEO using GPLId
    microarrayType <- Meta(platformData)$title # get microarray type from GEO
    ## Add ensembl and entrez id and symbol to probe ids on feature data:
    lst = AddIdSymbol(feature, microarrayType, ids, modelFile)
    lst[["microarrayType"]] <- microarrayType
    lst[["pheno"]] <- c('already in raw data')
  } else { # if study from GEO
    ## Get pheno (sample) and feature (probeset) metadata from GEO:
    seriesMatrixFilesList <- getGEO(Id, GSEMatrix = T) # GSEMatrix = T retrieves GSE Series Matrix files from GEO and output them as a list of ExpressionSets
    for (i in 1:length(seriesMatrixFilesList)) { # seriesMatrixFilesList has sometimes different expressionsets for different platforms used in same experiment. here appropriate platform is selected
      if (as.character(unique(seriesMatrixFilesList[[i]]$platform_id)) == GPLId) {
        seriesMatrixFile <- seriesMatrixFilesList[[i]]
      }
    }
    pheno <- phenoData(seriesMatrixFile) # retrieve pheno data
    feature <- featureData(seriesMatrixFile) # retrieve feature data
    ids <- pData(feature)$ID # get platform probe ids 
    gse <- getGEO(Id, GSEMatrix = F) # object of class GSE
    microarrayType <- Meta(GPLList(gse)[[GPLId]])$title # get microarray type
    ## Add ensembl and entrez id and symbol to probe ids on feature data:
      # some studies Affy 'HuGene' type don't have correct correspondence between feature data and raw data
      # so we exclude feature data:
    if(grepl('HuGene', microarrayType) | grepl('HuEx', microarrayType) | grepl('HTA', microarrayType)){ 
      lst = list(feature = vector(), InitialFinalIds = vector())
    } else {
      lst = AddIdSymbol(feature, microarrayType, ids, modelFile)
    }
    lst[["microarrayType"]] <- microarrayType
    lst[["pheno"]] <- pheno
  }
    return(lst)
}

# PreProcessRawExprILMN - read Illumina raw expression data from .txt file
#         - Arguments: * studyDir: string with path to directory of specific study
#                      * GEOId: string with GEO microarray platform identifier, eg.'GPL6947'
#         - Values: a dataframe with raw expression values
PreProcessRawExprILMN <- function (studyDir, Id) {
  # path to raw data txt.gz file:
  txtGzpath <- paste0(studyDir, '/', Id,'/', Id, '_non-normalized_data.txt.gz')
  # gunzip txt file:
  gunzip(txtGzpath, paste0(studyDir, '/', Id,'/', Id, '_non-normalized_data.txt'), overwrite = T, remove = F)
  # path to txt file:
  txtpath <- paste0(studyDir, '/', Id,'/', Id, '_non-normalized_data.txt')
  # read txt file:
  rawData <- read.table(paste0(txtpath), strip.white = T, sep='\t', blank.lines.skip = T, comment.char = '#' , header = T, row.names = 1, stringsAsFactors = F)
  # sometimes the file is sent from GEO in incorrect format (with lines containing ' "# ') so we need to overwrite it and open again:
  if (colnames(rawData)[1] == 'X') {
    rawData = rawData[!apply(rawData == "", 1, all),] # remove empty rows
    write.table(rawData, txtpath, quote = F, sep =  '\t', col.names = F) # first overwrite
    rawData <- read.table(paste0(txtpath), strip.white = T, sep='\t', blank.lines.skip = TRUE , header = T, row.names = 1, stringsAsFactors = F)
    
  }
  return (rawData)
}

# PreProcessRawExprILMNBeadStudioFormat - read Illumina raw expression data from file in BeadStudio/GenomeStudio's output file format (illuminia proprietary softwares for analyzing data
#         - Arguments: * studyDir: string with path to directory of specific study
#                      * Id: string with GEO microarray platform identifier, eg.'GPL6947'
#         - Values: a dataframe with raw expression values
PreProcessRawExprILMNBeadStudioFormat <- function (studyDir, Id) {
  filesInFolder <- list.files(paste0(studyDir, '/', Id))
  File <- filesInFolder[grepl('non-normalized', filesInFolder) & !grepl('.gz', filesInFolder)]
  # path to txt file:
  txtpath <- paste0(studyDir, '/', Id,'/', File)
  # read txt file:
  rawData <- readBeadSummaryData(dataFile = txtpath, controlID = "ProbeID", skip = 0, qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"))
  rawData <- data.frame(exprs(rawData), check.names = F)
  return (rawData)
}

# PreProcessRawExprAffy - read Affymetrix raw expression data from .CEL files downloaded from GEO
#         - Arguments: * studyDir: string with path to directory of specific study
#                      * GEOId: string with microarray study identifier, eg.'GSE...'
#         - Values: a expressionset containing raw expression values
PreProcessRawExprAffy <- function(studyDir, Id) {
  # path to tarfile:
  tarfile <- paste0(studyDir, '/', Id,'/',Id,'_RAW.tar')
  # list file names inside tarfile:
  FilesInside <- untar(tarfile, list = T)
  # select names of CEL files only;
  CELfiles <- grep('cel', FilesInside, value = T, ignore.case = T)
  # untar CEL files into subdirectory /CEL:
  untar(tarfile, files = CELfiles, exdir = paste0(studyDir,'/', Id, '/CEL'))
  # set path of CEL files:
  pathCELfiles <- paste0(studyDir, '/', Id,'/CEL/', CELfiles)
  # read CEL files:
  rawData <- read.celfiles(pathCELfiles)
  return (rawData)
}


# PreProcessRawExprAgil - read Agilent raw expression data from files downloaded from GEO
#         - Arguments: * studyDir: string with path to directory of specific study
#                      * GEOId: string with microarray study identifier, eg.'GSE...'
#         - Values: EListRaw containing raw expression values
PreProcessRawExprAgil <- function(studyDir, Id) {
  # path to tarfile:
  tarfile <- paste0(studyDir, '/', Id,'/',Id,'_RAW.tar')
  # list file names inside tarfile:
  FilesInside <- untar(tarfile, list = T)
  # untar files into subdirectory /CEL:
  untar(tarfile, files = FilesInside, exdir = paste0(studyDir,'/', Id, '/CEL'))
  # set path of files:
  pathCELfiles <- paste0(studyDir, '/', Id,'/CEL/', FilesInside)
  # read the gene expression intensity data:
  rawData <- read.maimages(pathCELfiles, source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")
  # extra column gIsWellAboveBG records whether intensity of each spot is above the background level, latter used for probe filtering
  return(rawData)
}

PreProcessRawExpr <- function(studyDir, Id, microarrayType, sampleCorr, illumnStudioFormat){
  if (grepl('E-', Id)) { # if it is from ArrayExpress
    if (grepl('Affy', microarrayType)) { # if it is affymetrix
      rawData = ArrayExpress(Id, path = paste0(studyDir, '/', Id), save = T) # reads raw expression data
    }
  } else { # if it is from GEO
    if (grepl('Illumina', microarrayType)) { # if it is illumina
      if (!is.na(illumnStudioFormat)) { # if dataset is in GenomeStudio or BeadStudio file output format (proprietary analysis softwares from illumina)
        print("Array with illumina Studio Format: manually download and process data as described on top of 'MicroarrayNormalize.R'")
        rawData = PreProcessRawExprILMNBeadStudioFormat(studyDir, Id)
      } else {
        getGEOSuppFiles(Id, baseDir = studyDir) # gets files with study raw data (it's a supplementary file at Id page)
        rawData = PreProcessRawExprILMN(studyDir, Id)
      }
    } else if (grepl('Affy', microarrayType)) { # if it is affymetrix
      getGEOSuppFiles(Id, baseDir = studyDir) # gets files with study raw data (it's a supplementary file at Id page)
      rawData = PreProcessRawExprAffy(studyDir, Id) 
    } else if (grepl('Agilent', microarrayType)) {
      getGEOSuppFiles(Id, baseDir = studyDir) # gets files with study raw data (it's a supplementary file at Id page)
      rawData = PreProcessRawExprAgil(studyDir, Id)
    }
  }
  ## preprocess sample names (to turn them into sample ids) on raw data:
  # if info on correspondence between initial sample id and Id sample id is provided,
  # otherwise use file names to get sample geo ids:
  if (length(sampleCorr) != 0) { # files have other alternative names
    if (class(rawData) == 'data.frame' | class(rawData) == 'ExonFeatureSet') {
      colnames(rawData) <- unlist(sampleCorr[colnames(rawData)])
    } else if (class(rawData) == 'ExpressionFeatureSet' | class(rawData) == 'ExpressionSetIllumina') {
      colnames(rawData) <- unlist(sampleCorr[colnames(exprs(rawData))])
    }
  } else if (length(sampleCorr) == 0) { # files just have the sample name (GEO sample Id)
    if (class(rawData) == 'data.frame' | class(rawData) == 'ExonFeatureSet' | class(rawData) == 'HTAFeatureSet') {
      colnames(rawData) <- unlist(lapply(strsplit(colnames(rawData), '(\\.)|(_)'), function(x) x[1])) # replace sample file names by sample names on matrix of expression raw values
    } else if (class(rawData) == 'ExpressionFeatureSet' | class(rawData) == 'GeneFeatureSet') { 
      colnames(rawData) <- unlist(lapply(strsplit(colnames(exprs(rawData)), '(\\.)|(_)'), function(x) x[1])) # replace sample file names by sample names on matrix of expression raw values
      # '|' means 'or'(an alternative separator) and '\\' is used for '.' to have literal meaning (instead of special meaning)
    } else if (class(rawData) == 'EListRaw'){
      colnames(rawData) <-  unlist(lapply(strsplit(colnames(rawData), '(\\.)|(/CEL/)'), function(x) unlist(strsplit(x[4], '_'))[1]))
    }
  }
  return(rawData)
}

# NormaliseExprILMN - normalises raw Illumina expression values with log2 transformation and quantile normalisation
#         - Arguments: * rawData: dataframe with Illumina raw expression values
#         - Values: matrix with Illumina expression values normalised with log2 and quantile transformations
NormaliseExprILMN <- function(rawData) {
  # apply log2 and transform into a matrix:
  if (class(rawData) == 'ExpressionSetIllumina') {
    rawDataCorr <- log2(backgroundCorrect(as.matrix(exprs(rawData)), method="normexp"))
      # normalizeBetweenArrays (below) requires log2-transformed input values when input is a matrix
  } else {
    rawDataCorr <- log2(backgroundCorrect(as.matrix(rawData), method="normexp"))
  }
  # apply quantile normalization: 
  normData <- normalizeBetweenArrays(rawDataCorr, method = "quantile") # appropriate for single-channel - 'withinarray' normalization is only for 2-channel microarray
  return(normData)
}

NormaliseExprAgil <- function(rawData) {
  # apply background correction ('normexp') and quantile normalization:
  rawDataBC <- backgroundCorrect(rawData, method="normexp")
  normData <- normalizeBetweenArrays(rawDataBC, method="quantile")
  ToKeep <- normData$genes[(rowSums(normData$other$gIsWellAboveBG > 0) >= 1) & (normData$genes$ControlType!=1L),"ProbeName"] # identify probes where signall is above background in at least one sample and that are not control probes
  normData = as.data.frame.EList(avereps(normData, ID=normData$genes$ProbeName))
  return(list(normData = normData, ToKeep = ToKeep))
}

# NormaliseExprAffyOld - normalises raw Affymetrix expression values with RMA normalisation for old affymetrix platforms (like U133 platforms)
#                      - it has 3 steps:
#                        * background correction
#                        * quantile normalization
#                        * calculation of gene expression: takes the different probes that exist for same probe set and outputs a single value
#                          With this method we run quantum normalization on the probe level & the probe level was summarized into the probe set level.
#                          if we had a normalization on probe set level. these distributions would have been identical
#                      - Arguments: * rawData: expression set with Affymetrix raw expression values
#                      - Values: matrix with Affymetrix expression values normalised with RMA
NormaliseExprAffyOld <- function(rawData, Id) {
  normData <- oligo::rma(rawData, background=TRUE, normalize=TRUE) # outputs an ExpressionSet.
  if (!grepl('E-', Id)) { # if study is NOT from ArrayExpress
    # convert into matrix:
    normData <- exprs(normData)
  }
  return(normData)
}

# NormaliseExprAffyNew - normalises raw Affymetrix expression values with RMA normalisation for new 'GENE' or 'EXON' affymetrix platforms (like ST platforms)
#                      - it has 3 steps:
#                        * background correction
#                        * quantile normalization
#                        * calculation of gene expression: takes the different probes that exist for same probe set and outputs a single value
#                          With this method we run quantum normalization on the probe level & the probe level was summarized into the probe set level.
#                          if we had a normalization on probe set level. these distributions would have been identical
#                      - Arguments: * rawData: expression set with Affymetrix raw expression values
#                      - Values: matrix with Affymetrix expression values normalised with RMA
NormaliseExprAffyNew <- function(rawData, Id) {
  normData <- oligo::rma(rawData, background=TRUE, normalize=TRUE, target="core") # outputs an ExpressionSet.
                                                                                      # target: character vector describing the summarization target.
                                                                                      #         Valid values are ’probeset’, ’core’ (Gene/Exon), ’full’ (Exon) and ’extended’ (Exon)
  if (!grepl('E-', Id)) { # if study is NOT from ArrayExpress
    # convert into matrix:
    normData <- exprs(normData)
  }
  return(normData)
}

NormaliseExpr <- function(microarrayType, rawData, Id){
  # apply RMA normalisation if microarray platform is affymetrix
  # log2 and quantile transformation if microarray platform is illumina
  # and background correction with quantile normalization for agilent
  if (grepl('Illumina', microarrayType)) {
    normData = NormaliseExprILMN(rawData)
    ToKeep = vector()
  } else if (grepl('Affy', microarrayType)) {
    if (grepl('HuGene', microarrayType) | grepl('HuEx', microarrayType)) { # for when Affy array platform is new
      normData = NormaliseExprAffyNew(rawData, Id)
      ToKeep = vector()
    } else {
      normData = NormaliseExprAffyOld(rawData, Id) # for when Affy array platform is old
      ToKeep = vector()
    }
  } else if (grepl('Agilent', microarrayType)) {
    normRes = NormaliseExprAgil(rawData)
    normData = normRes$normData
    ToKeep = normRes$ToKeep
  }
  return(list(normData = normData, ToKeep = ToKeep))
}


JoinExpPhenoFeat <- function(samplesList, pheno, normData, feature, treatments, InitialFinalIds, Id, ToKeep, microarrayType, modelFile, ...){
  if (grepl('E-', Id)) { # if it is from ArrayExpress
    ## Join feature data into expression set:
    # order feature data rows according to normData rows:
    rownames(feature) <- feature[,1] # give 1st column of feature as row names of feature
    featureFinal <- feature[rownames(normData),] # order rows
    # create a class AnnotatedDataFrame for feature data:
    featureAnn = AnnotatedDataFrame(featureFinal)
    # add feature data to expression set:
    featureData(normData) <- featureAnn
    ## Subset expressionset to have only samples that are relevant for comparisons we want to make:
    sampleNames <- unlist(samplesList) # Ids of samples we want to use
    eset <- normData[,sampleNames]
  } else {
    ## Make normalized expression matrix correspond to pheno and feature data:
    sampleNames <- unlist(samplesList) # Ids of samples we want to use
    # order sample names in columns of normalized expression matrix by order of sample names in pheno data:
    normData <- normData[, sampleNames(pheno)]
    # Affy 'HuGene' studies sometimes don' have feature data corresponding to normData,
    # for those tudies feature data is not included, otherwise it is included:
    if((grepl('HuGene', microarrayType) | grepl('HuEx', microarrayType) | grepl('HTA', microarrayType)) & grepl('Affy', microarrayType)){
      # Join normalized expression data and pheno data into a expressionset:
      eset <- new("ExpressionSet", exprs = normData, phenoData = pheno)
      ids = rownames(eset)
      #sum(!is.na(ConversionTable$ensembl_gene_id)) # 376741
      ConversionTableEnsembl = Microarray2EnsemblTableBiomart(microarrayType, ids)
      ConversionTableEnsembl = ConversionTableEnsembl[!is.na(ConversionTableEnsembl$ensembl_gene_id),]
      # among initial probeids corresponding to more than one final gene id, select the one that is also used in the metabolic model we will work with:
      modelGenestab <-  openxlsx::read.xlsx(modelFile, sheet = 'GENES') # all gene ids of recon3d01 metabolic model
      modelGenes <- modelGenestab[,'ID']
      finalGeneIdType <- colnames(ConversionTableEnsembl)[2]
      InitialFinalIds = chooseId(ConversionTableEnsembl, modelGenes, finalGeneIdType)
    } else {
      # order feature names in rows of normalized expression matrix by order of feature names in feature data:
      normData <- normData[featureNames(feature), ]
      # Join normalized expression data, pheno and feature data into expressionset:
      eset <- new("ExpressionSet", exprs = normData, phenoData = pheno, featureData = feature)
      #eset <- new("ExpressionSet", exprs = normData, phenoData = pheno)
    }
    ## Subset expressionset to have only samples that are relevant for comparisons we want to make:
    eset <- eset[,sampleNames]
    ## when microarray is Agilent only keep probes where signall is above background in at least one sample and that are not control probes:
    if (length(ToKeep) != 0) {
      eset <- eset[ToKeep,]
    }
  }
  # create a new column in pheno data (called 'group') that has type of treatment given to each sample:
  for (trt in treatments) {
    phenoData(eset)$group[sampleNames(eset)%in%samplesList[[trt]]] <- trt
  }
  ## Add ensembl and entrez id and symbol to each gene on normalized expression data:
  # Merge matrix of normalized expression data and initial-final gene ids dataframe:
  FinalIdsExprVal = mergeFinalIdsExprVal(eset, InitialFinalIds, "ensembl_gene_id")
  return(list(FinalIdsExprVal = FinalIdsExprVal, eset = eset))
}

summarizationFunct <- function(FinalIdsExprVal, organism){
  # this means we average expression of probe sets for same gene
  # NAs on Ensembl column are automatically removed
  ExpFinalGene <- aggregate(list(FinalIdsExprVal[colnames(FinalIdsExprVal)[3:ncol(FinalIdsExprVal)]]), by = list(ENSEMBL = FinalIdsExprVal$ensembl_gene_id), mean)
  colnames(ExpFinalGene) <- c('ENSEMBL', colnames(FinalIdsExprVal)[3:ncol(FinalIdsExprVal)])
  Ensembl2Entrez <- Ensembl2EntrezVectorAnnotationDb(organism, as.character(ExpFinalGene$ENSEMBL))
  Ensembl2Symbol <- Ensembl2SymbolVectorAnnotationDb(organism, as.character(ExpFinalGene$ENSEMBL))
  ExpFinalGene_all <- data.frame(Ensembl = ExpFinalGene$ENSEMBL, Entrez = Ensembl2Entrez, Symbol = Ensembl2Symbol, ExpFinalGene[,2:ncol(ExpFinalGene)], stringsAsFactors = F, check.names = F)
  return(ExpFinalGene_all)
}

SaveNormExpr <- function(studyDir, SummNormExpr, GroupNormExpr, eset, Id){
  normExpDir <- paste0(studyDir, '/NormData') # path to output directory
  dir.create(file.path(normExpDir), suppressWarnings('already exists')) # creates the new directory for normalized expression data
  sampleExpPath <- paste0(normExpDir,'/NormGeneExp_Samples.tab') # path to normalized expression file - samples
  groupExpPath <- paste0(normExpDir,'/NormGeneExp_Groups.tab') # path to normalized expression file - sample groups
  phenoPath <- paste0(normExpDir,'/sample_data.tab') # path to pheno (sample) data
  featurePath <- paste0(normExpDir,'/probe_set_data.csv') # path to feature (probe set id) data
  
  write.table(SummNormExpr, sampleExpPath, sep ='\t', quote=F, row.names = F) # write normalized expression gene values by sample
  write.table(GroupNormExpr, groupExpPath, sep ='\t', quote=F, row.names = F) # write normalized expression gene values by sample groups
  
  write.table(pData(phenoData(eset)), phenoPath, sep =  '\t', quote = F, row.names = F) # write pheno (sample) data
  write.table(pData(featureData(eset)), featurePath, sep = ',', quote = F, row.names = T) # write feature (probe set id) data
  # save expression set:
  saveRDS(eset, paste0(normExpDir, '/', Id, 'ExpressionSet.Rds'))
}

# AddIdSymbol -
## Add ensembl and entrez id and symbol to probe ids on feature data:
AddIdSymbol <- function(feature, microarrayType, ids, modelFile){
  ConversionEntrez = Microarray2EntrezVectorAnnotationDb(microarrayType, ids) # convert probe ids to gene entrez ids
  ConversionSymbol = Microarray2SymbolVectorAnnotationDb(microarrayType, ids) # convert probe ids to gene symbol
  feature$Entrez = ConversionEntrez
  feature$SYMBOL = ConversionSymbol
  ConversionTableEnsembl = Microarray2EnsemblTableBiomart(microarrayType, ids) # convert probe ids to entrez ids. gives table with one to many correspondence
  # among initial probeids corresponding to more than one final gene id, select the one that is also used in the metabolic model we will work with:
  modelGenestab <-  openxlsx::read.xlsx(modelFile, sheet = 'GENES') # all gene ids of recon3d01 metabolic model
  modelGenes <- modelGenestab[,'ID']
  finalGeneIdType <- colnames(ConversionTableEnsembl)[2]
  InitialFinalIds = chooseId(ConversionTableEnsembl, modelGenes, finalGeneIdType)
  # order ensembl gene id to correspond to probe id:
  dtf <- as.data.frame(rownames(feature))
  colnames(dtf) <- 'Microarray_transcript_ID'
  dtf2 <- merge(dtf, InitialFinalIds, by = 'Microarray_transcript_ID', all.x = T) # all.x = T allows to order rows by rownames(feature)
  feature$Ensembl = as.character(dtf2$ensembl_gene_id)
  return(list(feature = feature, InitialFinalIds = InitialFinalIds))
}

GroupNormExp = function(studySettingsTableName, studyId, SummNormExpr) {
  fToApply = function(x) {
    if(is.na(x['Patient']) & is.na(x['Subtype'])) {
      paste(x['Conditions'])
    } else if (is.na(x['Patient'])) {
      paste(x['Conditions'], x['Subtype'], sep = '_')
    } else if (is.na(x['Subtype'])) {
      paste(x['Conditions'], x['Patient'], sep = '_')
    } else {
      paste(x['Conditions'], x['Patient'], x['Subtype'], sep = '_')
    }
  }
  SamplesInfo <- openxlsx::read.xlsx(paste0(dataDir,'/',studySettingsTableName), sheet=1, colNames=T)
  studyInfo <- SamplesInfo[SamplesInfo$Study == studyId,]
  studyInfo$SampleGroup = apply(studyInfo, 1, function(x) fToApply(x))
  SampleGroupLst = lapply(split(studyInfo$SampleId, studyInfo$SampleGroup), unique)
  groupDataLst = lapply(SampleGroupLst, function(x) sapply(SummNormExpr[,x], as.numeric))
  groupMeanLst = lapply(groupDataLst, function(x) if(class(x) != 'numeric') {rowMeans(x)} else {x})
  GroupNormExpr = cbind(SummNormExpr[,1:3], do.call(cbind, groupMeanLst))
  return(GroupNormExpr)
}

