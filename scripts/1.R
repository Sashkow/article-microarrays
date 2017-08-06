library(sva)
path = '/home/rstudio/r/article-microarrays/preprocessed/affymetrix/E-GEOD-14722_preprocessed_affymetrix.tsv'
pathA = '/home/rstudio/r/article-microarrays/preprocessed/affymetrix/E-GEOD-14722_preprocessed_affymetrixA.tsv'
pathB = '/home/rstudio/r/article-microarrays/preprocessed/affymetrix/E-GEOD-14722_preprocessed_affymetrixB.tsv'

exprs = read.table(path, header = TRUE, sep = '\t')
exprsA = read.table(pathA, header = TRUE, sep = '\t')
exprsB = read.table(pathB, header = TRUE, sep = '\t')

exprsA$group = 'A'
exprsB$group = 'B'

ab = intersect(rownames(exprsA), rownames(exprsB))
ab

length(ab)

a.remove = exprsA[!rownames(exprsA) %in% ab, ]
nrow(a.remove)

intersect(rownames(a.remove), rownames(exprsB))

colnames(a.remove) = colnames(exprsB)
a.remove
exprsB
aandb = rbind(a.remove, exprsB)
nrow(aandb)

pdata = data.frame(aandb$group)
rownames(pdata) = rownames(aandb)

aandb <- subset(aandb, select = -c(group))



batch = as.factor(pdata$aandb.group)

aandb.t = t(aandb)

View(aandb)



mod = model.matrix(~1, data=pdata)
exprs.nobatch = ComBat(dat=aandb.t, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs.nobatch

write.table(t(exprs.nobatch), path, sep="\t")
t = read.table(path, sep="\t")



exprs = cbind(t(exprs.nobatch),pdata)
View(exprs)
expr
exprs = within(exprs, rm(aandb.group))

library(ggplot2)
library(reshape)

ncol(exprs)
melted = melt(exprs)
melted$group = 'A'
melted[1:nrow(melted)/2,]$group = 'B'



ggplot(melted, aes(x=value)) +
  geom_density(aes(colour=group, group=variable))




indA <- sample.int(12322, size=1232)
indB <- sample.int(9142, size=914)
meltedA <- melt(exprsA)
meltedA$group <- "A"
meltedB <- melt(exprsB)
meltedB$group <- "B"
melted <- rbind(meltedA, meltedB)

meltedA <- melt(exprsA[indA,])
meltedA$group <- "A"
meltedB <- melt(exprsB[indB,])
meltedB$group <- "B"

melted <- rbind(meltedA, meltedB)

ggplot(melted, aes(x=value)) +
  geom_density(aes(group=variable, colour=group))

### Remove batch caused by two different studies


