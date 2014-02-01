# STAT540 sem02: Playing with probability distributions and simulations
# Rebecca Johnston 15/01/2014


# Generate random numbers from a normal distribution ----------------------

# Set seed to regenerate the same random numbers
set.seed(540) # choose any positive integer
rnorm(10) # random numbers from a normal distribution, default mean = 0, sd = 1


# Generate lots of data, the R way ----------------------------------------

# Generate a matrix
# Specify 2 dimensions, each named under a variable
# Use these 2 variables within the rnorm arg together with matrix arg
n <- 10
B <- 4
(x <- matrix(rnorm(n * B), nrow = n))
str(x)
head(x)
tail(x)
summary(x)


# Generate lots of data, the awkward ways ---------------------------------

# For loop
x <- matrix(0, nrow = n, ncol = B)
for(j in 1:B) {
  x[ , j] <- rnorm(n)
}

# Gluing stuff together
sample1 <- rnorm(n)
sample2 <- rnorm(n)
sample3 <- rnorm(n)
sample4 <- rnorm(n)
x <- cbind(sample1, sample2, sample3, sample4)
# Repetitive code is a sure sign you're attacking from the wrong angle
# But the lack of scalability is the main problem here

# Gluing stuff together inside a for loop
x <- rnorm(n)
for(j in 1:(B - 1)) {
  x <- cbind(x, rnorm(n))
}


# Compute sample means and explore ----------------------------------------

# Create a matrix
x <- matrix(rnorm(n * B), nrow = n)

# Create matrix row and column names using sprintf arg
# sprintf = C-style command, returns a character vector containing a formatted
# combination of text and variable values
# must start with % and end with one of the letters in the set aAdifeEgGosxX%
# %02d = format number with up to 2 leading zeroes
rownames(x) <- sprintf("obs%02d", 1:n)
colnames(x) <- sprintf("samp%02d", 1:B)
x
dim(x)

# Compute sample means
mean(x[ , 2]) # sample mean for 2nd sample
colMeans(x)
apply(x, 2, mean) # same answer as colMeans, as 2 in "apply" denotes by column

# Exercise: Recall the claim that the expected value of the sample mean is the
# true mean. Compute the average of the 4 sample means we have. Is it (sort of)
# close to the true mean? Feel free to change n or B at any point.
mean(colMeans(x))
# Yes, the average of the sample means is close to the true mean (since
# generated from rnorm with defaults, the true mean is 0).

# Let's try other n and b values:
n <- 100
B <- 100
x <- matrix(rnorm(n * B), nrow = n)
mean(colMeans(x))

n <- 500
B <- 500
x <- matrix(rnorm(n * B), nrow = n)
mean(colMeans(x))

n <- 1000
B <- 1000
x <- matrix(rnorm(n * B), nrow = n)
mean(colMeans(x))
# Takes a while to compute, but the larger the sample size, the closer the
# sample mean to the true mean

# Exercise: Recall the Weak Law of Large Numbers said that, as the sample size
# gets bigger, the distribution of the sample means gets more concentrated
# around the true mean. Recall also that the variance of the sample mean is
# equal to the true data-generating variance divided by the sample size n.
# Explore these probability facts empirically. Don't go crazy -- just pick a few
# different sample sizes, compute sample means, and explore the variability of
# the sample means as a function of sample size.

# Sample size defined by number of rows
# Sample means calculated by column - keep this value constant
B <- 1000
x10 <- matrix(rnorm(10 * B), nrow = 10)
x100 <- matrix(rnorm(100 * B), nrow = 100)
x1000 <- matrix(rnorm(1000 * B), nrow = 1000)
x10000 <- matrix(rnorm(10000 * B), nrow = 10000)
x100000 <- matrix(rnorm(100000 * B), nrow = 100000)

# Calculate sample means using colMeans
xBar10 <- colMeans(x10)
xBar100 <- colMeans(x100)
xBar1000 <- colMeans(x1000)
xBar10000 <- colMeans(x10000)
xBar100000 <- colMeans(x100000)

# Calculate SD of colMeans using sd(colMeans)
xBarSd10 <- sd(colMeans(x10))
xBarSd100 <- sd(colMeans(x100))
xBarSd1000 <- sd(colMeans(x1000))
xBarSd10000 <- sd(colMeans(x10000))
xBarSd100000 <- sd(colMeans(x100000))

# Generate dataframe to summarise findings
cbind(sampSize = c(10, 100, 1000, 10000, 100000),
      trueSEM = 1 / sqrt(c(10, 100, 1000, 10000, 100000)),
      obsSEM = c(xBarSd10, xBarSd100, xBarSd1000, xBarSd10000, xBarSd100000))

# Jenny's more sophisticated solution, using more advanced R!!!
B <- 1000
n <- 10^(1:4)
names(n) <- paste0("n", n)
getSampleMeans <- function(n, B) colMeans(matrix(rnorm(n * B), nrow = n))
x <- sapply(n, getSampleMeans, B, simplify = FALSE)
cbind(sampSize = n,
      trueSEM = 1 / sqrt(n),
      obsSEM = sapply(x, sd),
      sampMeanIQR = sapply(x, IQR), # IQR = inter-quartile range
      sampMeanMad = sapply(x, mad)) # MAD = median absolute deviation

# Exercise: Repeat the above for a different distribution
# I have chosen the exponential distribution
# Jenny's code is super clean and generalisable: 
# only need to change "rnorm" to "rexp"
B <- 1000
n <- 10^(1:4)
names(n) <- paste0("n", n)
getSampleMeans <- function(n, B) colMeans(matrix(rexp(n * B), nrow = n))
x <- sapply(n, getSampleMeans, B, simplify = FALSE)
cbind(sampSize = n,
      trueSEM = 1 / sqrt(n),
      obsSEM = sapply(x, sd),
      sampMeanIQR = sapply(x, IQR), # IQR = inter-quartile range
      sampMeanMad = sapply(x, mad)) # MAD = median absolute deviation


# Compare probabilities and observed relative frequencies -----------------

# Exercise: Generate a reasonably large sample from some normal distribution (it
# need not be standard normal!). Pick a threshhold. What is the CDF at that
# threshhold, i.e. what's the true probability of seeing an observation less
# than or equal to the threshhold? Use your large sample to compute the observed
# proportion of observations that are less than threshhold. Are the two numbers
# sort of close? Hint: If x is a numeric vector, then mean(x <= threshhold)
# computes the proportion of values less than or equal to threshhold

# Define parameters
distSize <- 1000
distMean <- 100
distSd <- 50

# Generate large sample
set.seed(10)
x <- rnorm(distSize, distMean, distSd)
head(x)

# Choose a threshhold
thresh <- 138

# Expected CDF: use pnorm
(expCdf <- pnorm(thresh, distMean, distSd))

# Observed CDF: use less than or equal to symbol
(obsCdf <- mean(x <= thresh))

obsCdf - expCdf
# Yes the observed proportion and true probability are sort of close!


# Do the same for a different distribution
# I've chosen the Binomial distribution
# The parameters for rbinom and pbinom were confusing to figure out :S
# I think I got there in the end!
nObs <- 1000
sizeTrial <- 100
probTrial <- 0.8
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
thresh <- 75
(expCdf <- pbinom(thresh, sizeTrial, probTrial))
(obsCdf <- mean(x <= thresh))
obsCdf - expCdf

# Do the same for a variety of sample sizes. Do the two numbers tend to be
# closer for larger samples?
# n = 100
nObs <- 100
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
thresh <- 75
(expCdf <- pbinom(thresh, sizeTrial, probTrial))
(obsCdf <- mean(x <= thresh))
obsCdf - expCdf

# n = 10000
nObs <- 10000
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
thresh <- 75
(expCdf <- pbinom(thresh, sizeTrial, probTrial))
(obsCdf <- mean(x <= thresh))
obsCdf - expCdf

# n = 1000000
nObs <- 1000000
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
thresh <- 75
(expCdf <- pbinom(thresh, sizeTrial, probTrial))
(obsCdf <- mean(x <= thresh))
obsCdf - expCdf
# Yes, for larger sample sizes, the observed CDF approaches the expected CDF


# Instead of focusing on values less than the threshhold, focus on values
# greater than the threshhold
# Within rbinom and pbinom args use lower.tail = FALSE
nObs <- 10000
sizeTrial <- 100
probTrial <- 0.8
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
thresh <- 75
(expCdf <- pbinom(thresh, sizeTrial, probTrial, lower.tail = FALSE))
(obsCdf <- mean(x >= thresh))
obsCdf - expCdf


# Instead of focusing on tail probabilities, focus on the probability of the
# observed values falling in an interval
nObs <- 10000
sizeTrial <- 100
probTrial <- 0.8
x <- rbinom(n = nObs, size = sizeTrial, prob = probTrial)
lowThresh <- 73
uppThresh <- 79
mean(x >= lowThresh & x <= uppThresh)
pbinom(uppThresh, sizeTrial, probTrial) - 
  pbinom(lowThresh, sizeTrial, probTrial)
 

# Explore the distribution of sample means and the CLT --------------------
# CLT: central limit theorem

# Let's use ggplot2
install.packages("ggplot2")
library("ggplot2")

# Empirical distribution of sample means for various sample sizes
B <- 1000
n <- round(10^(seq(from = 1, to = 2.5, length = 4)), 0)
names(n) <- paste0("n", n)
getSampleMeans <- function(n, B) colMeans(matrix(rnorm(n * B), nrow = n))
x <- data.frame(sapply(n, getSampleMeans, B))

# Put data in tall format for ggplot2
xTallSkinny <- stack(x)
names(xTallSkinny) <- c("x","n")
xTallSkinny$n <- factor(xTallSkinny$n, levels = colnames(x))
head(xTallSkinny)
tail(xTallSkinny)

# Add histogram layer
p <- ggplot(xTallSkinny, aes(x, group = n, colour = n)) + 
  geom_density()

# Add strip plot layer
p <- p + geom_point(aes(x, y = -0.5, group = n, colour = n), 
                    position = position_jitter(height = 0.1),
                    alpha = 0.3)

# Add labels
p <- p + xlab("Sample means") +
  ylab("Density") +
  ggtitle("Exploring the Central Limit Theorem") +
  theme(plot.title = element_text(face = "bold"))
p
