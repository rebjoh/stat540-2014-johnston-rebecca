# STAT540 Sem05 Fitting and interpreting linear models (low volume)
# Rebecca Johnston 05/02/2014


# Load required packages and photoRec data --------------------------------
library("ggplot2") # for graphing
library("reshape2") # for converting data from wide to tall format
library("plyr") # for data aggregation tasks

# Load gene expression dataset
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)

# Load design matrix
prDes <- readRDS("GSE4051_design.rds")
str(prDes)


# Method A: Function to prepare a dataset for given genes  ----------------

# Write a function that takes any given affy probe ID and turns their associated
# values into a data frame.

# # METHOD A: constructing data frame manually
# # Based on Jenny's code from sem04 creating miniDat
# 
# # Index prDat with chosen genes
# luckyGenes <- c("1419655_at","1438815_at")
# luckyDf <- prDat[luckyGenes, ]
# 
# # Design data frame manually by column, starting with data in luckyDf, "gExp"
# # and "gene" then binding with prDes.
# 
# # For "gExp":
# # Convert luckyDf into a matrix, transpose, then convert into a vector. This is
# # done by row, so first 39 numbers are expression levels across 39 samples for
# # probe 1, then next 39 numbers are expression levels across 39 samples for
# # probe 2.
# 
# # For "gene": make rownames (probe ID) of luckyDf a factor, and repeat them each
# # 39 times (the number of columns in luckyDf)
# 
# luckyDf <- data.frame(gExp = as.vector(t(as.matrix(luckyDf))),
#                       gene = factor(rep(rownames(luckyDf), 
#                                     each = ncol(luckyDf)), levels = luckyGenes))
# 
# # Check samples are in the same order before column binding
# identical(rownames(t(luckyDf)), prDes$sidChar)
# luckyDf <- suppressWarnings(data.frame(prDes, luckyDf))
# 
# # We can now write the function where variable input is "luckyGenes"
# luckyGenes <- c("1419655_at","1438815_at")
# prepData <- function(x) {
#   luckyDf <- prDat[x, ]
#   luckyDf <- data.frame(gExp = as.vector(t(as.matrix(luckyDf))),
#                         gene = factor(rep(rownames(luckyDf), 
#                                           each = ncol(luckyDf)), 
#                                       levels = luckyGenes))
#   luckyDf <- suppressWarnings(data.frame(prDes, luckyDf))  
# }
# jDat <- prepData(luckyGenes)
# str(jDat)
# head(jDat)


# Method B: Function to prepare a dataset for given genes  ----------------

# Write a function that takes any given affy probe ID and turns their associated
# values into a data frame.

# METHOD B: melting and merging
# Index prDat with chosen genes
luckyGenes <- c("1419655_at","1438815_at")
luckyDf <- prDat[luckyGenes, ]
head(luckyDf)

# Reshape and bind with design matrix
# Add rownames (probe ID) as a column for melt function first
luckyDf <- data.frame(gene = rownames(luckyDf), luckyDf) 
luckyDf <- melt(luckyDf, id.vars = "gene", variable.name = "sidChar", 
                value.name = "gExp")
head(luckyDf)
luckyDf <- merge(luckyDf, prDes, by = "sidChar")
# Reorder using plyr function arrange
luckyDf <- arrange(luckyDf, devStage, gType)

# We can now write the function where variable input is "luckyGenes"
luckyGenes <- c("1419655_at","1438815_at")
prepData <- function(x) {
  luckyDf <- prDat[x, ]
  luckyDf <- data.frame(gene = rownames(luckyDf), luckyDf) 
  luckyDf <- melt(luckyDf, id.vars = "gene", variable.name = "sidChar", 
                  value.name = "gExp")
  luckyDf <- merge(luckyDf, prDes, by = "sidChar")
  luckyDf <- arrange(luckyDf, devStage, gType)
}
jDat <- prepData(luckyGenes)
str(jDat)
head(jDat)

# Plot data using ggplot
ggplot(jDat, aes(gExp, gType, colour = gType)) + 
  geom_point(alpha = 0.5) + facet_wrap(~ gene)

ggplot(jDat, aes(devStage, gExp, colour = gType, group = gType)) + 
  geom_point() + facet_wrap(~ gene) +
  stat_summary(fun.y = mean, geom = "line")


# Write a function to stripplot a mini-dataset ----------------------------

# You will probably make lots of these plots. Why not write a function for this?
makeStripplot <- function(x) {
  ggplot(x, aes(devStage, gExp, colour = gType, group = gType)) + 
    geom_point() + facet_wrap(~ gene) +
    stat_summary(fun.y = mean, geom = "line")
}
makeStripplot(jDat)

# You can use both of your functions together and create a minidatset and plot
# it all at once:
makeStripplot(newDat <- prepData("1456341_a_at"))
str(newDat)
head(newDat)


# Do a two-sample t-test --------------------------------------------------

# Let's test for a difference in expected gene expression for probeset
# "1456341_a_at" at developmental stage P2 vs. 4 weeks post-natal (ignoring
# genotype, i.e. lump the wild types and knockouts together). Let's assume a
# common variance in the two groups.

levels(newDat$devStage)
t.test(newDat$gExp[newDat$devStage == "P2"], 
       newDat$gExp[newDat$devStage == "4_weeks"])

# OR subset data!
t.test(gExp ~ devStage, newDat, subset = devStage == "P2" | 
         devStage == "4_weeks")


# Fit a linear model with a categorical covariate -------------------------

# In other words, do "one-way ANOVA". Focus on probeset "1438786_a_at".

# Use function we already prepared to subset data and draw stripplot:
makeStripplot(mDat <- prepData("1438786_a_at"))
str(mDat)

# Now call linear model using "lm" function
mDatlm <- lm(gExp ~ devStage, mDat, subset = gType == "wt")
summary(mDatlm)

# This result looks plausible: the intercept of 8.52 is close to the mean of E16
# (the reference). And the values for the other devStages change accordingly
# (-/+ slope wrt devStageE16)


# Perform inference for a contrast ----------------------------------------

# The "W" shape of the expression profile for "1438786_a_at" means that the
# expression values for developmental stages P2 and P10 are quite similar. We
# could formally test whether the P2 and P10 effects are equal or, equivalently,
# whether their difference is equal to zero.

# First extract the parameter estimates from the linear model you fit above.
# Hint: the coef() function will pull parameter estimates out of a wide array of
# fitted model objects in R.
coef(mDatlm)

# Now you need to construct the contrast matrix to form the difference between
# the P2 and P10 effects. I called mine contMat. Hint: it's OK for a matrix to
# have just one row.

# Check order of devStage levels
levels(mDat$devStage)
# Contrast matrix for observed difference in sample mean for wt mice at P2 - P10
contMat <- matrix(c(0, 1, 0, -1, 0), nrow = 1)
(obsDiff <- contMat %*% coef(mDatlm)) # %*% means matrix multiplication!

# Let's check that this really is the observed difference in sample mean for the
# wild type mice, P2 vs. P10.
(sampMeans <- aggregate(gExp ~ devStage, mDat, FUN = mean,
                        subset = gType == "wt"))

with(sampMeans, gExp[devStage == "P2"] - gExp[devStage == "P10"])

# Yes! Agrees with the observed difference we computed by multiplying our
# contrast matrix and the estimated parameters. If you don't get agreement, you
# have a problem ... probably with your contrast matrix.

# Now we need the (estimated) standard error for our contrast. The
# variance-covariance matrix of the parameters estimated in the original model
# can be obtained with vcov() and is equal to (X^TX)^−1σ^^2.
vcov(mDatlm)

# Let's check that this is really true. If we take the diagonal elements and
# take their square root, they should be exactly equal to the standard errors
# reported for out original model. Are they?
summary(mDatlm)$coefficients[ , "Std. Error"]
sqrt(diag(vcov(mDatlm)))

# Yes! Note for the future that you can get the typical matrix of inferential
# results from most fitted model objects for further computing like so:
summary(mDatlm)$coefficients

# Returning to our test of the P2 vs. P10 contrast, recall that the
# variance-covariance matrix of a contrast obtained as Cα^ is C(X^TX)^−1C^Tσ^^2.
(estSe <- contMat %*% vcov(mDatlm) %*% t(contMat))

# Now we form a test statistic as an observed effect divided by its estimated
# standard error:
(testStat <- obsDiff/estSe)

# Under the null hypothesis that the contrast equals zero, i.e. there is no true
# difference in mean for expression at P2 and P10 in wild type mice for this
# gene, the test statistic has a t distribution with n−p=20−5=15 degrees of
# freedom. We compute a two-sided p-value and we're done.
2 * pt(abs(testStat), df = df.residual(mDatlm), lower.tail = FALSE)

# Not surprisingly, this p-value is rather large and we conclude there is no
# difference.


# Fit a linear model with two categorical covariates ----------------------

# Let's focus on probeset "1448690_at". Use your functions to prepare the data
# and plot it.
makeStripplot(oDat <- prepData("1448690_at"))
str(oDat)

# Fit a linear model with covariates gType and devStage and include their
# interactions.

# Include interactions = multiply categorical covariates?
oFitBig <- lm(gExp ~ gType * devStage, oDat)
round(summary(oFitBig)$coefficients, digits = 5)

# Vet the results. Is the intercept plausible? How about the various effects? Do
# the ones with small p-values, e.g. meeting a conventional cut-off of 0.05,
# look 'real' to you?

# P-values for wtdevStageP10 and wtdevStageP6 look real, significant p-value for
# wtdevStage4_weeks also seems real. As for the intercepts for the interaction 
# of gType and devStage, not sure about this? Why are the estimate values 
# positive? Now I understand. This is the interaction effect, the observed
# average - expected average (intercept + KO effect + devStageP2 effect)

# Omit interactions = sum categorical covariates?
oFitSml <- lm(gExp ~ gType + devStage, oDat)
# No interaction means the lines must be parallel
summary(oFitSml)$coefficients

# Now let's determine if the interaction terms are truly necessary. From the
# plot, the case for interaction seems very weak. This can be assessed with an F
# test that essentially looks at the reduction in the sum of squared residuals
# due to using a larger, more complicated model and determines if it is "big
# enough" given the number of additional parameters used. Recall the anova()
# function can take two fitted models, one nested within the other, and perform
# this test. (anova() can also be used on a single model to assess significance
# of terms, but remember the problem with the standard anova() function and
# unbalanced data. See references given in lecture for remedies.)
anova(oFitSml, oFitBig)

# With a p-value awfully close to one, we confirm that, no, there is no evidence
# for interaction in this particular case.

# If you'd like to get a more exciting result, take a look at probeset
# "1429225_at". Here are my plots, excerpts from the fitted model reports, and
# the F test for interaction

makeStripplot(pDat <- prepData("1429225_at"))
str(pDat)
# Big model, with interaction = 10 parameters: 
# gType * devStage 
# = 2 * 5 = 10
pFitBig <- lm(gExp ~ gType * devStage, pDat)
round(summary(pFitBig)$coefficient, digits = 5)
# Sml model, no interaction = 6 parameters: 
# (gType - 1) + (devStage - 1) + 1 
# = 1 + 4 + 1 = 6
pFitSml <- lm(gExp ~ gType + devStage, pDat)
round(summary(pFitSml)$coefficient, digits = 5)
anova(pFitSml, pFitBig)

# Not surprisingly, the interaction here is highly statistically significant.


# Ideas for further work --------------------------------------------------

# We wrote functions to prepare and plot data for more than 1 gene. But when we
# started fitting models and conducting tests, we only worked with 1 gene at a
# time. Can you use data aggregation strategies from last week to do some of the
# same work for small sets of genes?

# 

