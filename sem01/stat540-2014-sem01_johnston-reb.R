# STAT540 Seminar 01: Exploration of small gene expression dataset
# Gene expression data from photo receptor cells in mice
# Rebecca Johnston 11/01/2014

# DATA DETAILS
# Gene expression data from various developmental stages, two genotypes
# Here we will work with gene expression data for only 3 genes
# Gene names have been replaced with 3 Pokemon attacks
# crabHammer, eggBomb, and poisonFang



# Load GSE4051_MINI.txt "prDat" -------------------------------------------

# Save GSE4051_MINI.txt to wd from STAT540 github page
# https://github.com/jennybc/stat540_2014/tree/master/examples/photoRec/data
prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)


# Basic exploration of prDat ----------------------------------------------

str(prDat)
summary(prDat)

# How many rows are there?
nrow(prDat)

# How many columns or variables are there?
ncol(prDat)

# Inspect the first few observations or the last few or a random sample
head(prDat)
tail(prDat)
prDat[5:16, ]
prDat[sample(nrow(prDat), size = 6), ] # From Jenny's code! Cool!

# What does row correspond to -- different genes or different mice?
# Answer = Rows correspont to different mice (samples)

# What are the variable names?
names(prDat)

# What "flavor" is each variable, i.e. numeric, character, factor?
str(prDat)

# For "sample", does each integer between 1 and nrows occur exactly once?
nrow(prDat) # length(prDat$sample)
table(prDat$sample)
# Sort sample numbers in ascending order
sort(prDat$sample)
# Create sequence of numbers from 1 to number of rows in dataset
seq(1, nrow(prDat))
# Determine if these two variables are identical
sort(prDat$sample) == seq(1, nrow(prDat))
all.equal(sort(prDat$sample), seq(1, nrow(prDat)))
identical(sort(prDat$sample), seq(1, nrow(prDat)))

# For each factor variable, what are the levels? 
levels(prDat$devStage)
levels(prDat$gType)

# How many observations do we have for each level of devStage? For gType?
summary(prDat$devStage)
summary(prDat$gType)

# Perform a cross-tabulation of devStage and gType
table(prDat$devStage, prDat$gType)

# If you had to take a wild guess, what do you think the intended experimental 
# design was? What actually happened in real life?
# Answer = From the cross-tabulation above of genotype (gType) and developmental 
# stage (devStage), I would assume that the intended experimental design was to
# compare two different genotypes in mice (wild type and knock out) across 5
# different developmental stages. There were meant to be 4 mice for each 
# of the 5 developmental stages, and 4 mice for each of the two genotypes, but
# it appears in "real life", one of the knock out mice (NrlKO) at developmental
# stage E16 was not able to be obtained (e.g. KO coukd have been lethal?)

# For each quantitative variable, what are the extremes? 
# How about average or median?
# Since there are multiple calculations we want to perform on multiple variables
# Can we turn this into a function?
# Yes but not sure how to do this when arg returns multiple values (e.g. range)

statSumm <- function(x){
  data.frame(min = min(x), max = max(x), range = range(x))
}
statSumm(prDat$eggBomb)

# Could simply use summary for the variable names we want
summary(prDat[ , c("crabHammer", "eggBomb", "poisonFang")])
# Other arguments
range(prDat$crabHammer) # returns min and max
fivenum (prDat$crabHammer) # returns Tukey's five number summary


# Indexing and subsetting prDat -------------------------------------------

# Create a new data.frame called weeDat only containing observations for which
# expression of poisonFang is above 7.5
prDat$poisonFang > 7.5
# use "which" argument and indexing
weeDat <- prDat$poisonFang[which(prDat$poisonFang > 7.5)]
str(weeDat)

# For how many observations poisonFang > 7.5? 
# How do they break down by genotype and developmental stage?
length(weeDat)
# prDat[prDat$poisonFang %in% weeDat, ] 
prDat[prDat$poisonFang %in% weeDat, c("poisonFang", "gType", "devStage")]
# prDat[prDat$poisonFang > 7.5, ]

# Print the observations with row names "Sample_16" and "Sample_38" to screen,
# showing only the 3 gene expression variables.
prDat[c("Sample_16", "Sample_38"), c("crabHammer", "eggBomb", "poisonFang")]

# Which samples have expression of eggBomb less than the 0.10 quantile?
sort(prDat$eggBomb)
quantile(prDat$eggBomb, probs = 0.1) # Answer = 6.1844
# quantile does not return a number. Can you code this instead of quoting it?
prDat[prDat$eggBomb < 6.1844, c("eggBomb", "sample")]
# Yes you can! Just add the two lines together and voila! Jenny's code:
rownames(prDat[prDat$eggBomb < quantile(prDat$eggBomb, 0.1), ])
