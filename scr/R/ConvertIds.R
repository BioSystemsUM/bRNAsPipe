#################
# Load Packages #
#################

library(biomaRt)
library(openxlsx)

#############
# Functions #
#############

# Microarray2EnsemblTableBiomart - * converts microarray gene id to ensembl gene id USING BIOMART
#                                  * gives one to many correspondence
#                                  * when there is no correspondent ensembl id, the original microarray gene id is REMOVED from result,
#                                    so if we remove duplicated entries, number of obtained ensembl ids is SMALLER that number of microarray ids
#                                  * order of microarrays ids is NOT SAME in output as in input (ids)   
#                                  * sometimes BIOMART is capable of converting more ids sometimes is ANNOTATIONDBI 
#                                - Arguments: * microarrayType: string of format "Illumina HumanHT-12 V3.0 expression beadchip"
#                                                               and got from 'Meta(GPLList(gse)[[GPLId]])$title' where gse is of class GSE (from GEOQuery) and GPLId is of format 'GPL6947'
#                                                               ** if from affymetrix - string of format 'pd.hg.u133.plus.2' got from 'annotation(seriesMatrixFile)' where seriesMatrixFile is an ExpressionSet
#                                             * ids: vector strings of gene ids from a microarray type eg."1007_s_at" "1053_at"   "117_at"    "121_at"    "1255_g_at" "1294_at" 
#                                - Values: * it is a dataframe of microarray specific gene ids and ensembl gene ids
#                                          * may have more than one ensembl gene id for a microarray specifc id
#                                          * has NAs as ensembl gene id when no id was found for a array gene id
Microarray2EnsemblTableBiomart <- function (microarrayType, ids) { 
  if(grepl('HuGene', microarrayType) | grepl('HuEx', microarrayType) | grepl('HTA', microarrayType)) { # if platform is affymetrix 'HuGene'
    st <- unlist(strsplit(microarrayType,']'))[1]
    microarrayPackage <- paste0(gsub('-|_','',tolower(substr(st,2,nchar(st)))), 'transcriptcluster.db')
    library(microarrayPackage, character.only=T)
    ConversionTable <- select(eval(parse(text = microarrayPackage)), keys = as.character(ids), column = c('PROBEID','ENSEMBL'), keytype='PROBEID')
    colnames(ConversionTable)[2] = 'ensembl_gene_id'
    } else if (startsWith(microarrayType, 'Illumina')) {
    org <- substring(microarrayType, 10, 14)
    st <- gsub(' expression beadchip', '', microarrayType)
    micrTypeFilter = gsub('.0', '', tolower(gsub(' |-', '_', st)))
    dataset = "hsapiens_gene_ensembl"
    ensembl = useEnsembl("ensembl",dataset)
    ConversionTable <- getBM(attributes=c(micrTypeFilter,'ensembl_gene_id'),filters = micrTypeFilter, values = ids ,mart = ensembl)
  } else if (grepl('Affy', microarrayType) & !grepl('st', microarrayType)) { # if affymetrix NOT 'ST'
    org <- substring(microarrayType,2,3)
    st <- unlist(strsplit(microarrayType,']'))[1]
    micrTypeFilter = paste0('affy_', tolower(gsub('-', '_', substr(st,2,nchar(st)))))
    #if (grepl('st', micrTypeFilter)) {micrTypeFilter = paste0(micrTypeFilter, '_v1')}
    dataset = "hsapiens_gene_ensembl"
    ensembl = useEnsembl("ensembl",dataset)
    ConversionTable <- getBM(attributes=c(micrTypeFilter,'ensembl_gene_id'),filters = micrTypeFilter, values = ids ,mart = ensembl)
  } else if (grepl('Agilent', microarrayType)) {
    if (grepl('SurePrint.*G3.*GE.*v2', microarrayType)) {
      micrTypeFilter = 'agilent_sureprint_g3_ge_8x60k_v2' # checks if microarrayType has 'SurePrint', 'G3', 'GE' and 'v2' with variable number of characters in between
      dataset = "hsapiens_gene_ensembl"
      ensembl = useEnsembl("ensembl",dataset)
      ConversionTable <- getBM(attributes=c(micrTypeFilter,'ensembl_gene_id'),filters = micrTypeFilter, values = ids ,mart = ensembl)
    } else if (grepl('SurePrint.*G3.*GE', microarrayType)) {
      micrTypeFilter = 'agilent_sureprint_g3_ge_8x60k'
      dataset = "hsapiens_gene_ensembl"
      ensembl = useEnsembl("ensembl",dataset)
      ConversionTable <- getBM(attributes=c(micrTypeFilter,'ensembl_gene_id'),filters = micrTypeFilter, values = ids ,mart = ensembl)
    } else if (grepl('Agilent-014850', microarrayType)) {
      micrTypeFilter = 'agilent_wholegenome_4x44k_v1'
      dataset = "hsapiens_gene_ensembl"
      ensembl = useEnsembl("ensembl",dataset)
      ConversionTable <- getBM(attributes=c(micrTypeFilter,'ensembl_gene_id'),filters = micrTypeFilter, values = ids ,mart = ensembl)
    } else if (grepl('Agilent-027114', microarrayType)) {
      microarrayPackage <- paste0(tolower(gsub('-', '', unlist(strsplit(microarrayType,' '))[1])), '.db')
      library(microarrayPackage, character.only=T)
      ConversionTable <- select(eval(parse(text = microarrayPackage)), keys = as.character(ids), column = c('PROBEID','ENSEMBL'), keytype='PROBEID')
      colnames(ConversionTable)[2] = 'ensembl_gene_id'
    }
  }
  # when previous doesn't work use other archived versions:
  # ensembl = useEnsembl("ensembl",dataset, mirror = "useast")
  # ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL", dataset = dataset)
  # ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2012.archive.ensembl.org")
  # ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org")
  return (ConversionTable)
}

# Microarray2EntrezVectorAnnotationDb - * converts microarray probe id to entrez gene id USING ANNOTATIONDBI
#                                        * gives one to one correspondence (arbitrarily selects id)
#                                     - Arguments: * microarrayType: string with this format 'pd.hg.u133.plus.2'
#                                                  * ids: vector strings of probe ids from a microarray type eg."1007_s_at" "1053_at"  "117_at"
#                                     - Values: * it is a vector of microarray specific pobe ids and entrez gene ids
#                                               * only has one gene id for one microarray specific probe id
Microarray2EntrezVectorAnnotationDb <- function(microarrayType, ids) {
  if(startsWith(microarrayType, 'Illumina')) {
    st <- gsub(' expression beadchip', '', microarrayType)
    st2 <- gsub('.0', '.db',gsub('\\d\\d','',gsub(' |-|HT','', st)))
    microarrayPackage <- gsub('V','v',paste0(tolower(substr(st2,1,1)), substr(st2,2,nchar(st2))))
  } else if (grepl('Affy', microarrayType)) {
    if (grepl('HuGene', microarrayType)) { # if affy microarray platform is of recent type
      st <- unlist(strsplit(microarrayType,']'))[1]
      microarrayPackage <- paste0(gsub('-|_','',tolower(substr(st,2,nchar(st)))), 'transcriptcluster.db')
    } else { # if affy microarray platform is of old type
      st <- unlist(strsplit(microarrayType,']'))[1]
      microarrayPackage <- paste0(gsub('-|_','',tolower(substr(st,2,nchar(st)))), '.db')
    }
  } else if (grepl('Agilent', microarrayType)) {
    microarrayPackage <- paste0(tolower(gsub('-', '', unlist(strsplit(microarrayType,' '))[1])), '.db')
  }
  library(microarrayPackage, character.only=T)
  ConversionEntrez <- mapIds(eval(parse(text = microarrayPackage)), keys = as.character(ids), column = 'ENTREZID', keytype='PROBEID')
  return(ConversionEntrez)
}

# Microarray2SymbolVectorAnnotationDb - * converts microarray probe id to gene symbol USING ANNOTATIONDBI
#                                       * gives one to one correspondence (arbitrarily selects symbol)
#                                     - Arguments: * microarrayType: string with this format 'pd.hg.u133.plus.2'
#                                                  * ids: vector strings of probe ids from a microarray type eg."1007_s_at" "1053_at" "117_at" 
#                                     - Values: * it is a vector of microarray specific probe ids and a vector of gene symbols
#                                               * only has one gene symbol for one microarray specifc gene id
Microarray2SymbolVectorAnnotationDb <- function(microarrayType, ids) {
  if(startsWith(microarrayType, 'Illumina')) {
    st <- gsub(' expression beadchip', '', microarrayType)
    st2 <- gsub('.0', '.db',gsub('\\d\\d','',gsub(' |-|HT','', st)))
    microarrayPackage <- gsub('V','v',paste0(tolower(substr(st2,1,1)), substr(st2,2,nchar(st2))))
  } else if(grepl('Affy', microarrayType)) {
    if (grepl('HuGene', microarrayType)) { # if affy microarray platform is of recent type
      st <- unlist(strsplit(microarrayType,']'))[1]
      microarrayPackage <- paste0(gsub('-|_','',tolower(substr(st,2,nchar(st)))), 'transcriptluster.db')
    } else { # if affy microarray platform is of old type
      st <- unlist(strsplit(microarrayType,']'))[1]
      microarrayPackage <- paste0(gsub('-|_','',tolower(substr(st,2,nchar(st)))), '.db')
    }
  } else if (grepl('Agilent', microarrayType)) {
    microarrayPackage <- paste0(tolower(gsub('-', '', unlist(strsplit(microarrayType,' '))[1])), '.db')
  }
  library(microarrayPackage, character.only=T)
  ConversionSymbol <- mapIds(eval(parse(text = microarrayPackage)), keys = as.character(ids), column = 'SYMBOL', keytype='PROBEID')
  return(ConversionSymbol)
}

# Ensembl2EntrezVectorAnnotationDb - * converts enseml gene id to entrez gene ID USING ANNOTATIONDBI
#                                       * gives one to one correspondence (arbitrarily selects entrez Id)
#                                       * sometimes ANNOTATIONDBI is capable of converting more ids sometimes BIOMART
#                                     - Arguments: * organism: string with this format 'homo_sapiens'
#                                                  * ids: vector strings of ensembl gene ids eg."ENSG00000000003" "ENSG00000001036" 
#                                     - Values: * it is a vector of ensembl gene ids and a vector of entrez gene ids
Ensembl2EntrezVectorAnnotationDb <- function(organism, ids) {
  if (organism == 'homo_sapiens') {organismPackage <- 'org.Hs.eg.db'}
  library(organismPackage, character.only=T) # character.only=T cause organismPackage is in string format
  ConversionEntrez <- mapIds(eval(parse(text = organismPackage)), keys = ids, column = 'ENTREZID', keytype='ENSEMBL') # we use mapIds() instead of select() although select() gives not only final gene id but also initial gene id
  # cause otherwise we would have to manualy choose one alternative final id when one to many correspondence exist (note: select() with multiVals = 'first' doesn't work)
  return(ConversionEntrez)
}

# Ensembl2SymbolVectorAnnotationDb - * converts enseml gene id to gene symbol USING ANNOTATIONDBI
#                                       * gives one to one correspondence (arbitrarily selects symbol)
#                                       * sometimes ANNOTATIONDBI is capable of converting more ids sometimes BIOMART
#                                     - Arguments: * organism: string with this format 'homo_sapiens'
#                                                  * ids: vector strings of ensembl gene ids eg."ENSG00000000003" "ENSG00000001036" 
#                                     - Values: * it is a vector of ensembl gene ids and a vector of gene symbols
# 
Ensembl2SymbolVectorAnnotationDb <- function(organism, ids) {
  if (organism == 'homo_sapiens') {organismPackage <- 'org.Hs.eg.db'}
  library(organismPackage, character.only=T) # character.only=T cause organismPackage is in string format
  ConversionEntrez <- mapIds(eval(parse(text = organismPackage)), keys = ids, column = 'SYMBOL', keytype='ENSEMBL') # we use mapIds() instead of select() although select() gives not only final gene id but also initial gene id
  # cause otherwise we would have to manualy choose one alternative final id when one to many correspondence exist (note: select() with multiVals = 'first' doesn't work)
  return(ConversionEntrez)
}

# chooseSameId - it is called inside chooseId function
#              - it is applied to each element of the list (to each initial id) 
#              - it selects among the optional final ids corresponding to the same initial id, the final id that is also in a list we provide
chooseSameId <- function(x, providedGenes) {
  optionalFinalIds <- x[,2][x[,2]%in%providedGenes] # final gene ids that are also gene ids in list provided
  if (length(x[,2]) == 1) { # if for one initial gene id there is only one final gene id
    chosenid <- x[,2] # final id is kept
  } else  if (length(optionalFinalIds) == 1) { # if for one initial gene id there are several final gene ids, but only one final gene id is also present in provided list
    chosenid <- optionalFinalIds # final id is kept
  } else if (length(optionalFinalIds) > 1) { # if for one initial gene id there are several final gene ids, and more than one of those is also present in provided list
    chosenid <- optionalFinalIds[1] # the first final id to be in provided list is selected
  } else if (length(optionalFinalIds) == 0) { # if for one initial gene id there are several final gene ids, but none of them is present in provided list
    chosenid <- x[,2][1] # the first final id is selected
  }
  InitialFinal <- c(unique(x[,1]), chosenid)
  return(InitialFinal)
}

# chooseId - Often, during conversion of gene ids, some initial gene ids correspond to more than one final gene id (one to many correspondence)
#            Among those this function selects the one that is also used in a list of ids the user provides
#         - Arguments: * ConversionTable: a table with initial ids in 1st column (may contain repeated ids) and final ids in the 2nd
#                      * providedGenes: character vector with ids that should be kept if are also final ids in ConversionTable
#                      * finalGeneIdType: string to identify the type of final gene ids, eg. 'ensembl_gene_id
#         - Values: a table with initial ids (not repeated) in 1st column and final ids in 2nd column.
#                   selected final ids are the ones also present in a provided list
chooseId <- function(ConversionTable, providedGenes, finalGeneIdType){
  # get a list where each element is a dataframe corresponding to one initial gene id and has all corresponding final gene ids:
  groupedDuplicates <- split(ConversionTable, ConversionTable[,1]) # initial gene ids
  InitialFinalIdsList <- lapply(groupedDuplicates, chooseSameId, providedGenes) # to each group of initial gene ids apply a function "chooseSameId"
  InitialFinalIds <- unique(Reduce(rbind, InitialFinalIdsList)) # joins elements of list into a dataframe
  colnames(InitialFinalIds) <- c('Microarray_transcript_ID', finalGeneIdType)
  return(InitialFinalIds)
}

# mergeFinalIdsExprVal - merges initial-final gene ids dataframe with matrix of initial gene ids and normalized expression values for each sample
#                      - Arguments: * eset: matrix with normalized gene expression values
#                                   * InitialFinalIds: dataframe with initial (1st column) and final gene ids (2nd column).
#                                                      one initial id corresponds to only one final id and vice-versa
#                                   * finalGeneIdType: string to identify the type of final gene ids, eg. 'ensembl_gene_id
mergeFinalIdsExprVal <- function(eset, InitialFinalIds, finalGeneIdType) {
  # convert normalized data from matrix into dataframe and add rownames as new colum:
  ExprValFrame <- data.frame(rownames(eset),exprs(eset), stringsAsFactors =  F, check.names = F)
  colnames(ExprValFrame)[1] <- 'Microarray_transcript_ID'
  # merge dataframes:
  FinalIdsExprVal <- merge(ExprValFrame, InitialFinalIds, by = 'Microarray_transcript_ID', all = T)
  FinalIdsExprVal <- FinalIdsExprVal[c('Microarray_transcript_ID', finalGeneIdType, colnames(eset))]
  return(FinalIdsExprVal)
}

