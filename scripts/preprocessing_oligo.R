source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")
biocLite("oligo")
biocLite("Biobase")
library(oligo)
library(Biobase)


library(affy)
library(ArrayExpress)

setwd('/home/rstudio/r/article-microarrays')
getwd()

rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc/', sep='/')


# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
studies$processingSoft


# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# install cdf annotation files for all listed microarray platforms
for (array in levels(studies$platformAbbr)){
  install.brainarray(array)
}

i = 7
source('~/r/article-microarrays/scripts/install.brainarray.R')
install.brainarray(studies$platformAbbr[[i]])



# i = 6 E-GEOD-36083


current_path = paste(rawspath, '/', studies$accession[[i]], sep='')
if (! dir.exists(current_path)){
  dir.create(current_path)
}

aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  # local = TRUE,
  type = 'raw')

z = Biobase::read.AnnotatedDataFrame(filename = paste(current_path, aeData$sdrf, sep='/'))

# z@data$Source.Name = paste(z@data$Assay.Name, " 1", sep = '')  # for E-GEOD-73374


# merge ArrayExpress phenodata with IGEA phenodata
pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')


pd$Experiment
rownames(pd) = pd$Array.Data.File

pd$Experiment

rownames(pd)


path = '/home/rstudio/r/article-microarrays/raws/affymetrix'

gseString = "E-GEOD-73374"

# gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste(path, gseString, sep = "/")

setwd(celFilesPath)





filePaths = paste(celFilesPath, rownames(pd), sep = "/")
oligoData = oligo::read.celfiles(filenames = filePaths, phenoData = pd)

normalizedData = oligo::rma(oligoData)

arrayQualityMetrics::arrayQualityMetrics(expressionset = normalizedData,
                                         outdir = "QC_oligo",
                                         force = TRUE,
                                         intgroup = "Target")

