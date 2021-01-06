# Load Packages and scripts:
library('org.Hs.eg.db')
library('openxlsx')
library('edgeR')

args = commandArgs(T)
Script2Load = args[1]
source(Script2Load)


# Variables:
Study = args[2]
BaseDir= args[3]
SamplesInfoFile = paste(BaseDir, args[4], sep='/')
CountsDir = paste(BaseDir, args[5], sep = '/')
AnnotationDir= args[6]
NormDir= paste(BaseDir, args[7], sep = '/')

### Join sample raw counts into a table (with samples of same condition side by side):
CountsTable = JoinRawCounts(CountsDir, SamplesInfoFile)

### Add entrez, ensembl and symbol gene ids to table with raw counts:
RawCounts = AddEntrezEnsemblSymbol(CountsTable)

### Normalise counts with GeTMM:
GenLenPath = paste(BaseDir, AnnotationDir, 'GeneLength.bed', sep='/')
SampleNormCounts = normalise_GeTMM(RawCounts, GenLenPath)
# save normalised counts for each sample:
write.table(SampleNormCounts, paste(NormDir, 'NormGeneExp_Samples.tab', sep = '/'), quote=F, row.names = F, col.names = T, sep ='\t')
GroupNormCounts = GroupNormExp(SamplesInfoFile, Study, SampleNormCounts)
# save normalised counts for each group of samples:
write.table(GroupNormCounts, paste(NormDir, 'NormGeneExp_Groups.tab', sep = '/'), quote=F, row.names = F, col.names = T, sep ='\t')

