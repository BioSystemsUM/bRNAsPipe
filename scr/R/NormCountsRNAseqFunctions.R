# chooseId - Often, during conversion of gene ids, some initial gene ids correspond to more than one final gene id (one to many correspondence)
#           so among those we try to select the one that is also used in a list of ids the user provides
#         - Arguments: * ConversionTable: a table with initial ids in 1st column (also with repeated ids) and final ids in the 2nd
#                      * providedGenes: character vector with ids that should be kept if are also final ids in ConConversionTable
#                      * finalGeneIdType: string to identify the type of final gene ids, eg. 'ensembl_gene_id
#         - Values: a table with initial ids (not repeated) in 1st column and final ids in 2nd column.
#                   selected final ids (from optional several final ids of same initial id) are the ones also present in a provided list
# chooseSameId - it is called inside chooseId function
#              - it is applied to each element of the list (to each initial id) 
#              - it selects among the optional final ids corresponding to the same initial id, the final id that is also in a list we provide 
chooseSameId <- function(x, providedGenes) {
  optionalFinalIds <- x[,2][x[,2]%in%providedGenes] # get a vector with final gene ids that are also gene ids in base model (eg.HMR2)
  if (length(x[,2]) == 1) { # if for one initial gene id there is only one final gene id
    chosenid <- x[,2] # then the final id is kept
  } else  if (length(optionalFinalIds) == 1) { # if for one initial gene id there are several final gene ids, but only one final gene id is also present in the model
    chosenid <- optionalFinalIds # then that final id is kept
  } else if (length(optionalFinalIds) > 1) { # if for one initial gene id there are several final gene ids, and more than one of those is also present in the cell model
    chosenid <- optionalFinalIds[1] # then one of the final ids that are in the cell model is abitrarily selected (the 1st in the list)
  } else if (length(optionalFinalIds) == 0) { # if for one initial gene id there are several final gene ids, but none of them is present in the cell model
    chosenid <- x[,2][1] # then we arbitrarily select one of the final gene ids (select the 1st in the list)
  }
  InitialFinal <- c(unique(x[,1]), chosenid)
  return(InitialFinal)
}

chooseId <- function(ConversionTable, providedGenes, initialGeneIdType, finalGeneIdType){
  # dataframe with all final and initial gene ids is converted to a list.
  # each element of the list is a dataframe corresponding to one initial gene id and has all final gene ids corresponding to that initial gene id.
  groupedDuplicates <- split(ConversionTable, ConversionTable[,1]) # column 1 of ConversionTable has initial gene ids
  InitialFinalIdsList <- lapply(groupedDuplicates, chooseSameId, providedGenes) # to each group of initial gene ids apply a function "chooseSameId"
  # in lapply, extra arguments (here 'modelGenes') of the function to apply to each list element must be provided after
  # like 'apply(list, function1 ,extrargument)' where function1 is defined like function1 <- function(x, extrargument){}
  InitialFinalIds <- unique(Reduce(rbind, InitialFinalIdsList)) # joins elements of list into a dataframe
  colnames(InitialFinalIds) <- c(initialGeneIdType, finalGeneIdType) # rename columns of initial final id convertion ConversionTable: initial ID is Microarray_Id and final ID is Gene_ID
  return(InitialFinalIds)
}


# JoinRawCounts - joins samples' raw counts into a table (with samples of same condition side by side)
#               - adds entrez and ensembl gene ids
#               - Arguments: * CountsDir: path to directory where each file is raw counts of rnaseq alignments of one sample
#                            * SamplesInfoFile: path to file containing study samples info. This file has: Study, StudyId, Conditions, Patient, CellSubtype, GEOSamples, ENAExperimentAccession, ENARunAcession, ArrayExpressExperiment, ArrayExpressSample, Link
#               - Values: * Table: a table with gene 
JoinRawCounts = function(CountsDir, SamplesInfoFile) {
  # Join datasets' raw counts into a table:
  CountsFilesPath = list.files(paste(CountsDir, sep = '/') , pattern = '.*\\.RawCounts', full.names = T)
  Table = read.table(CountsFilesPath[1], stringsAsFactors = F, header=F)[,1]
  for (Path in CountsFilesPath) {
    Counts = read.table(Path, stringsAsFactors = F, header=F)[,2]
    Table=cbind(Table, Counts)
  }
  Table = data.frame(Table, row.names = 1, stringsAsFactors = F)
  colnames(Table) = sub('.RawCounts', '', basename(CountsFilesPath))
  # order samples so that samples of same condition are side by side:
  SamplesInfo <- read.table(SamplesInfoFile, stringsAsFactors = F, header = T, fill=T)
  SampleNamesOrdered <- unique(SamplesInfo[SamplesInfo$Study == Study, 'SampleId'])
  Table = Table[,SampleNamesOrdered]
  return(Table)
}

AddEntrezEnsemblSymbol = function(CountsTable) {
  if (substr(rownames(CountsTable)[1], 1,4) == 'ENSG') { # if primary ID is ensembl, get entrez id and symbol:
    # convert ensembl id to entrez id (get only one entrez id from those corresponding to same ensembl id):
    Entrez = mapIds(org.Hs.eg.db, keys = rownames(CountsTable), column = 'ENTREZID', keytype ='ENSEMBL')
    # remove genes where conversion to entrez id was not found:
    Entrez = Entrez[!is.na(Entrez)]
    # get symbol from ensembl id:
    Symbol = mapIds(org.Hs.eg.db, keys = Entrez, column = 'SYMBOL', keytype ='ENTREZID')
    # merge ids into the table with counts:
    Int = cbind(Entrez, Symbol)
    CountsTableOrd = CountsTable[row.names(Int),] # order table with counts
    CountsTableOrd = sapply(CountsTableOrd, as.numeric)
    FinalTable = data.frame('Ensembl' = as.character(row.names(Int)), 'Entrez' = as.character(Int[,'Entrez']), 'Symbol' = as.character(Int[,'Symbol']) , CountsTableOrd, stringsAsFactors = F)
  } else { # else, primary ID is symbol and we need to get ensembl and entrez ids:
    # convert gene symbol to entrez id (get only one entrez id from those corresponding to same symbol):
    Entrez = mapIds(org.Hs.eg.db, keys = rownames(CountsTable), column = 'ENTREZID', keytype ='SYMBOL')
    # remove genes where conversion to entrez id was not found:
    Entrez = Entrez[!is.na(Entrez)]
    # get ensembl ids from entrez id:
    Ensembl = mapIds(org.Hs.eg.db, keys = Entrez, column = 'ENSEMBL', keytype ='ENTREZID')
    # merge ids into the table with counts:
    Int = cbind(Entrez, Ensembl)
    CountsTableOrd = CountsTable[row.names(Int),] # order table with counts
    CountsTableOrd = sapply(CountsTableOrd, as.numeric)
    FinalTable = data.frame('Entrez' = as.character(Int[,'Entrez']), 'Ensembl' = as.character(Int[,'Ensembl']), 'Symbol' = as.character(row.names(Int)), CountsTableOrd, stringsAsFactors = F)
  }
  return(FinalTable)
}

normalise_GeTMM = function(RawCounts, GenLenPath){
  GenLen = read.table(GenLenPath, row.names=1, header=T, stringsAsFactors=F)
  if (substr(rownames(GenLen)[1], 1,4) == 'ENSG') { # if primary ID is ensembl:
    geneMetaData = RawCounts[,1:3]
    geneCountData = data.matrix(RawCounts[,4:ncol(RawCounts)])
    # Determine reads per kilobase (RPK):
    geneEnsembl = geneMetaData$Ensembl
    for(samp in colnames(geneCountData)){
      sampleValues = geneCountData[,samp]
      CountsPerLen = sampleValues/(GenLen[geneEnsembl,'length'])
      RpkSample = CountsPerLen * (10^3)
      geneCountData[,samp] = RpkSample
    }
  } else { # else, primary ID is symbol:
    # Note: Gtftools is not able to calculate gene length for all genes in gtf file,
    #       pseudogenes of NCBI gtf are excluded.
    #       So, we have to exclude those pseudogenes from analysis
    RawCounts = RawCounts[RawCounts$Symbol%in%row.names(GenLen),]
    geneMetaData = RawCounts[,1:3]
    geneCountData = data.matrix(RawCounts[,4:ncol(RawCounts)])
    # Determine reads per kilobase (RPK):
    geneSymbol = geneMetaData$Symbol
    for(samp in colnames(geneCountData)){
      sampleValues = geneCountData[,samp]
      CountsPerLen = sampleValues/(GenLen[geneSymbol,'length'])
      RpkSample = CountsPerLen * (10^3)
      geneCountData[,samp] = RpkSample
    }
  }
  # Apply TMM normalization from edgeR package:
  edgeRdata = edgeR::DGEList(geneCountData)
  edgeRdata = edgeR::calcNormFactors(edgeRdata)
  edgeRdata = edgeR::cpm(edgeRdata, normalized.lib.sizes=T)
  final = data.frame(geneMetaData, edgeRdata, stringsAsFactors = F)
  return(final)
}

GroupNormExp = function(SamplesInfoFile, Study, SampleNormCounts) {
  fToApply = function(x) {
    if(is.na(x['Patient']) & is.na(x['CellSubtype'])) {
      x['Conditions']
    } else if (is.na(x['Patient'])) {
      paste(x['Conditions'], x['CellSubtype'], sep = '.')
    } else if (is.na(x['CellSubtype'])) {
      paste(x['Conditions'], x['Patient'], sep = '.')
    } else {
      paste(x['Conditions'], x['Patient'], x['CellSubtype'], sep = '.')
    }
  }
  SamplesInfo <- read.table(SamplesInfoFile, stringsAsFactors = F, header = T, fill=T)
  studyInfo <- SamplesInfo[SamplesInfo$Study == Study,]
  studyInfo$SampleGroup = apply(studyInfo, 1, function(x) fToApply(x))
  SampleGroupLst = lapply(split(studyInfo$SampleId, studyInfo$SampleGroup), unique)
  groupDataLst = lapply(SampleGroupLst, function(x) sapply(SampleNormCounts[,x], as.numeric))
  groupMeanLst = lapply(groupDataLst, function(x) ifelse(class(x) != 'numeric', rowMeans(x), x))
  GroupNormCounts = cbind(SampleNormCounts[,1:3], do.call(cbind, groupMeanLst))
  return(GroupNormCounts)
}


