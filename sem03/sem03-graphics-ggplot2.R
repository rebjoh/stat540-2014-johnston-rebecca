# STAT540 sem03
# Introduction to R graphics: Exploration of ggplot2
# Rebecca Johnston 22/01/13


# Load required packages --------------------------------------------------
library("ggplot2") # for graphing
library("hexbin") # for making hexbin objects
library("reshape2") # for reshaping data from wide to tall format (melt)


# Basic ggplot2 concepts and args -----------------------------------------

# Layers: ggplot2 graphics are built using different layers
# Aesthetics: "aes()"
# Other arguments:
apropos("^geom_") # Geometries "geom_"
apropos("^stat_") # Statistics "stat_"
apropos("^scale_") # Scale "scale_"
apropos("^facet_") # Facetting "facet_"


# Read photoRec data ------------------------------------------------------

# Save rds file from STAT540 github
# https://github.com/jennybc/stat540_2014/blob/master/examples/photoRec/data
# Note this is a different format compared to sem01
# There are now 7 variables including sidChar and sidNum
kDat <- readRDS("GSE4051_MINI.rds")

# Sanity checks to ensure the data loaded successfully:
str(kDat)
head(kDat)
table(kDat$devStage, kDat$gType)


# Scatter plots -----------------------------------------------------------

# Scatter plots "geom_point"
(p <- ggplot(kDat, aes(x = crabHammer, y = eggBomb)) + 
  geom_point())

# Add smoothing line "stat_smooth"
(p <- p + stat_smooth())

# Change background and add labels
(p <- p + theme_bw() + # white background
   xlab("Expression of crabHammer") + # x axis label 
   ylab("Expression of eggBomb") + # y axis label 
   ggtitle("Scatterplot for expression levels")) # plot title


# RESAHPE DATA
# Convert into long format so there is a single column of gene expression values
# Only keep genes "eggBomb" and "poisonFang"
# Not sure why we keep crabHammer in a separate column? Use it as a reference? Weird.
nDat <-
  with(kDat,
       data.frame(sidChar, sidNum, devStage, gType, crabHammer,
                  probeset = factor(rep(c("eggBomb", "poisonFang"), each = nrow(kDat))),
                  geneExp = c(eggBomb, poisonFang)))
str(nDat)

# ADD COLOUR AND SMOOTHING LINES
# Assign colour by probeset column under aes arg
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point())

# Add a smoothing line
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point() + 
   stat_smooth(se = FALSE)) # se parameter: display standard error around smooth

# Overrule groupings in aes by specifiying new group so only one smoothing line is plotted
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point() + 
   stat_smooth(se = F, aes(group = 1)))


# FACETTING
(p <- ggplot(nDat, aes(crabHammer, geneExp)) + 
   geom_point() + 
   facet_wrap(~ probeset))

# Add colour by gType
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = gType)) + 
   geom_point() + 
   facet_wrap(~ probeset))

# Add colour by devStage
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = devStage)) +
   geom_point() +
   facet_wrap(~ probeset))


# Strip plots -------------------------------------------------------------

# Also use geom_point for strip plots, but coordinates are mapped to a factor
# RESHAPE DATA: Long format so all genes in "probeset" this time!
oDat <-
  with(kDat,
       data.frame(sidChar, sidNum, devStage, gType,
                  probeset = factor(rep(c("crabHammer", "eggBomb",
                                          "poisonFang"), each = nrow(kDat))),
                  geneExp = c(crabHammer, eggBomb, poisonFang)))
str(oDat)

# Plot geneExp by probe using geom_point
(p <- ggplot(oDat, aes(geneExp, probeset)) + 
   geom_point())

# Add jitter
(p <- ggplot(oDat, aes(geneExp, probeset)) + 
   geom_point(position = position_jitter(height = 0.1)))

# Plot geneExp by devStage
# Note order in aes matters! Specify x and y!
(p <- ggplot(oDat, aes(devStage, geneExp)) + 
   geom_point())

# Facetting by probeset
(p <- p + facet_wrap(~ probeset))

# Colour by gType
(p <- p + aes(color = gType))

# Add some stats: mean for each group
(p <- p + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4))


# Density plots -----------------------------------------------------------

# Gene expression across entire dataset
(p <- ggplot(oDat, aes(geneExp)) + 
   geom_density())

# Remove border lines
(p <- ggplot(oDat, aes(geneExp)) + 
   stat_density(geom = "line", position = "identity"))

# Add strip plot with jitter
(p <- ggplot(oDat, aes(geneExp)) + 
   stat_density(geom = "line", position = "identity") + 
   geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

# Change bandwidth with adjust arg in stat_density
(p <- ggplot(oDat, aes(geneExp)) + 
   stat_density(geom = "line", position = "identity", adjust = 0.5) + 
   geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

# Define color by gType
(p <- ggplot(oDat, aes(geneExp, color = gType)) + 
   stat_density(geom = "line", position = "identity") + 
   geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

# Define colour by devStage
(p <- ggplot(oDat, aes(geneExp, color = devStage)) + 
   stat_density(geom = "line", position = "identity", adjust = 0.75) + 
   geom_point(aes(y = 0.05), position = position_jitter(height = 0.01)))

# Add facet wrap by devStage
(p <- p + facet_wrap(~ devStage))


# Boxplots  ---------------------------------------------------------------

# Plot geneExp by devStage using geom_boxplot
(p <- ggplot(oDat, aes(devStage, geneExp)) + 
   geom_boxplot())

# Add facet wrap by gType
(p <- p + facet_wrap(~ gType))

# As a violin plot instead
(p <- ggplot(oDat, aes(devStage, geneExp)) + 
   geom_violin())


# Overplotting and plot matrix --------------------------------------------

# Load entire photorec dataset into R (save file from STAT540 github)
# Gene expression data
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)
head(prDat)
tail(prDat)

# Load design matrix
prDes <- readRDS("GSE4051_design.rds")
str(prDes)

# Choose 2 samples (columns) at random
set.seed(2)
(yo <- sample(1:ncol(prDat), size = 2))

# Create new dataframe with only these 2 random samples
bDat <- data.frame(y = prDat[[yo[1]]], z = prDat[[yo[2]]])
str(bDat)

# Plot samples against each other using geom_point
(p <- ggplot(bDat, aes(z, y)) + 
   geom_point())

# Reduce transparency using alpha arg in geom_point
(p <- ggplot(bDat, aes(z, y)) + 
   geom_point(alpha = 0.1))

# Use stat_density2d function
(p <- ggplot(bDat, aes(z, y)) + 
   stat_density2d())

# Fill with colour
(p <- ggplot(bDat, aes(z, y)) + 
   stat_density2d(geom = "tile", contour = FALSE, aes(fill = ..density..)) + 
   scale_fill_gradient(low = "white", high = "blue"))
# Some functions (especially stat) return their own calculated values
# Use these values by calling ..[value name].., e.g. the use of fill = ..density.. here

# Use stat_binhex
(p <- ggplot(bDat, aes(z, y)) + 
   stat_binhex())

# Use plotmatrix for pairwise scatterplots
set.seed(3)
(yo <- sample(1:ncol(prDat), size = 4))
pairDat <- subset(prDat, select = yo)
str(pairDat)
(p <- plotmatrix(pairDat) + 
   stat_binhex())


# Take home problem -------------------------------------------------------

# The full photoRec dataset has 39 samples and 29,949 probesets. Choose 2 or 20
# or 200 random probesets/genes and look for gene expression differences between
# the two genotypes, wild type versus knockout. Make use of the graphing 
# techniques discussed this week such as scatter plots, box plot, etc.

# Choose 20 random genes to start
set.seed(50)
(someGene <- sample(nrow(prDat), 20))

# Create new dataframe with only these 20 random genes
tDat <- prDat[someGene,]
str(tDat, max.level = 0)
head(tDat)

# What structure do we need the data to be in? To use ggplot2, we will need
# separate columns for: probeID, sampleID, gExp, gType.

# # What about creating a dataframe manually by writing a function? 
# # First "stack" the gene expression values:
# stDat <- stack(tDat)
# head(stDat) # Lose probeID
# 
# # Write function
# makeDf <- function(a, b, c){
#   # a = tDat
#   # b = stDat
#   # c = prDes
#   probeID <- as.factor(row.names(a))
#   gExp <- b$values # I don't like this code :S better way to get gExp values?
#   gType <- c$gType
#   sID <- c$sidChar 
#   data.frame(probeID, gExp, gType)
# }
# 
# # Call function
# newDf <- makeDf(tDat, stDat, prDes)
# head(newDf)

# I don't know if I trust the above result... particularly whether the gene 
# expression values match the correct sample... Is there a simple way I could
# check this? I will use the reshape package instead.

# First need to make row names a column
tDat <- data.frame(probeID = row.names(tDat), tDat)
head(tDat, n = 3)

# Reshape to tall format using melt, where ID variable is "probeID"
tDat <- melt(tDat, id.vars = "probeID")
head(tDat, n = 3)
# Rename columns
names(tDat) <- c("probeID", "sidChar", "gExp")
# Merge prDes to obtain gType using "sidChar" as common column
tDat <- merge(tDat, prDes, by = "sidChar")
# For simplicity, rename probes to smaller names
levels(tDat$probeID) <- c(LETTERS[1:length(levels(tDat$probeID))])

# ggplot strip plot
# Gene expression per gene
(strPlot <- ggplot(tDat, aes(gExp, probeID, colour = gType)) +
  geom_point())

# This is a bit messy... Let's reorder probeID based on gExp
tDat <- within(tDat, probeID <- reorder(probeID, gExp))
# Replot
(strPlot <- ggplot(tDat, aes(gExp, probeID, colour = gType)) +
   geom_point())

# ggplot boxplot
# I couldn't figure out how to put box plots horizontally?

# Gene expression across all probes comparing gType
(bxPlot <- ggplot(tDat, aes(gType, gExp)) +
   geom_boxplot())

# Gene expression per probeID comparing gType
(bxPlot <- ggplot(tDat, aes(probeID, gExp, 
                            fill = interaction(gType))) +
   geom_boxplot() + 
   theme()) 

# ggplot density plot
# Gene expression across all 20 genes comparing gType
(densPlot <- ggplot(tDat, aes(gExp, colour = gType)) +
   stat_density(geom = "line", position = "identity"))

# Gene expression per probeID comparing gType
(densPlot <- ggplot(tDat, aes(gExp, colour = gType)) +
  stat_density(geom = "line", position = "identity") +
  facet_wrap(~ probeID))

# Using ggplot visualisation tools alone, it appears that for the 20 probes that
# have been randomly chosen here, there is no significant difference between
# expression for wt and KO mice.