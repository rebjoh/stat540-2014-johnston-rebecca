STAT540 Seminar 10: Supervised learning, classification, cross validation, variable selection
---
_Rebecca Johnston 19/03/2014_

* <dim id="1a">[Introduction](#1b)
    * <dim id="2a">[Load required packages](#2b)
    * <dim id="3b">[Load dataset and perform data preparation](#3b)
* <dim id="4a">[Classification](#4b)
    * <dim id="5a">[Feature and Model Selection](#5b)
        * <dim id="6a">[Cross validation splits](#6b)
        * <dim id="7a">[Loop for feature selection and modeling](#7b)
        * <dim id="8a">[Error Rates](#8b)
        * <dim id="9a">[Results of the CV](#9b)
    * <dim id="10a">[Testing the selected model](#10b)
    * <dim id="11a">[CMA](#11b)


### <dim id="1b">[Introduction](#1b) 

In this Seminar we go over packages and codes for performing supervised learning and evaluation in R. In supervised learning, one is given a training data set with known response and covariates and our goal is to predict the response in a new test set for which we only know the covariates. A supervised learning process can be decomposed into the following steps:

*Step 1*: Select Features. Before training a model, in many applications, it is usually important to perform a pre-filtering step in which one retains only the most informative features (e.g., genes) as candidate "biomarkers". The amount of features retained to select and train a model is up to the analyst and the methods used in the next steps. For example, some methods may be unfeasible or too slow to run with a large number of features.

*Step 2*: Select and train a classifier. Once the set of candidate markers have been selected, the next step is to select and train a model to predict the labels of a test data.

*Step 3*: Test. Finally, a model is chosen and used to predict labels of a test data.


**Cross-validation**  
Although it makes sense to proceed as described above, many methods are available and many constants within methods need to be selected in these steps. Thus, a *cross-validation* is usually required to *evaluate* how well different trained models work and select the *best* model to proceed. Note that although you may only want to select among different choices available in Step 2, the cross-validation needs to start in Step 1. Why? The results of the cross-validation will be over-optimistic and biased if the samples in the test sets of the cross-validation (i.e., left-out folds) were used to *select* the most promising features in Step 1!! For example, if the performance of a complex model is (artificially) good, you may not penalize regression coefficients enough in Step 2, and may yield to a poor performance in Step 3.

In many studies, in the absence of a test set, cross-validation is used to estimate performance. In those cases, a *nested cross-validation* is required! The inner cross-validation will be used to select features and tune parameters, the outer cross-validation will be used to a selected test model. Further readings in this topic are available in [1][Ref1] and [2][Ref2].

In this seminar, we will work with a dataset that has both a training and a test set. Thus, we will not do a nested cross-validation. However, keep it in mind for your project or future work!


#### <dim id="2b">[Load required packages](#2b)

There are many packages for performing supervised learning in R, each of which may implement one or more algorithms. There have also been at least two major efforts to unify these libraries under a common framework to make them easier to use: `MLInterfaces` and `CMA`. Although these may be useful and save you a lot of time in your analysis, it is important that you understand what these packages are doing and what they are *not* doing. Thus, I will not use these packages in this Seminar but I encourage you to reproduce the analysis using at least one of them! (I recommend `CMA`).

Install the following packages from Bioconductor: `CMA` and `GEOquery`, and from CRAN: `ROCR`, `car`, `e1071` (for SVM), and `glmnet`.





```r
library(MASS)
library(reshape)
library(car)
library(limma)
library(e1071)
library(glmnet)
library(ROCR)
library(CMA)
library(GEOquery)
library(lattice)
library(ggplot2)
library(class)
```



#### <dim id="3b">[Load dataset and perform data preparation](#3b)

This seminar is based on a dataset that comes from a paper by [Smeets et al. 2010](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23177), who studied Affymetrix expression profiles from primary breast tumors. Smeets group was interested in whether whether tumors which had spread to lymph nodes (LN positive, generally a bad sign) have different gene expression profiles than LN negative tumors. If so, a gene expression signature can be use to predict tumor class.

Their data set contains 24236 genes on 116 samples. The status of the lymph node is known for each sample, with 59 LN positive and 57 LN negative. Samples were divided into two parts: 96 samples (48 LN positive and 48 LN negative) were used as a "training" set and 20 samples (11 LN positive and 9 LN negative) were used as a "test" set. There is also a quantitative measure, "LnRatio", the fraction of affected lymph nodes, presumably reflecting "how bad" the LnStatus is. Thus, we can use this dataset to illustrate classification and regularization methods! This seminar will focus on the first task, i.e., classification. In the past, Paul selected this dataset to illustrate the challenges of supervised learning tasks!

In the paper, the authors trained a support vector machine classifier to distinguish between LN positive and LN negative tumors (i.e., classification), and evaluated the results using ROC curves. After some optimization, they got an area under the ROC curve (AUC) of 0.66 on the training set and 0.65 on the test data. This is better than chance, but still not very convincing results about the relevance of the derived molecular signature (random chance would give an AUC of 0.5; perfect classification would give an AUC of 1.0).

**Data Preparation**  
First, let's retrieve our datasets from GEO with `getGEO` from `GEOquery` package. Warning: this may take several minutes! So to avoid re-downloading in the future, save the data once you get it into a good shape.


```r
if (file.exists("class_LNstatus.Rdata")) {
    # if previously downloaded
    load("class_LNstatus.Rdata")
} else {
    # if downloading for the first time takes a several mins! returns a list
    datgeo <- getGEO("GSE23177", GSEMatrix = TRUE)
    dat <- datgeo[[1]]  # Note that dat is an ExpressionSets
    str(pData(dat), max.level = 0)
    
    # extract only those variables of interest
    pData(dat) <- subset(pData(dat), select = c("characteristics_ch1.2", "characteristics_ch1.3", 
        "characteristics_ch1"))
    names(pData(dat)) <- c("LnStatus", "LnRatio", "Set")
    
    # Note: LNRatio will not be used in this Seminar. However, you can use it to
    # try some of the regularization techniques learned in class
    
    # split the ExpressionSet into training and test sets.
    train.es <- dat[, dat$Set == "patient type: training set"]
    test.es <- dat[, dat$Set != "patient type: training set"]
    
    # Re-label factor
    pData(train.es)$LnStatus <- recode(pData(train.es)$LnStatus, "levels(pData(train.es)$LnStatus)[1] =' neg'; else = 'pos'", 
        levels = c("neg", "pos"))
    
    pData(test.es)$LnStatus <- recode(pData(test.es)$LnStatus, "levels(pData(test.es)$LnStatus)[1] = 'neg'; else = 'pos'", 
        levels = c("neg", "pos"))
    
    # create data matrices with expression values (probesets in rows). Some of
    # the functions we will use do not take ExpressionSets as objects
    trainDat <- exprs(train.es)
    testDat <- exprs(test.es)
    
    # Redefine the quantitative variable LnRatio to make it a numeric variable.
    ntrain <- dim(pData(train.es))[1]
    ntest <- dim(pData(test.es))[1]
    
    pData(train.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(train.es)$LnRatio)), 
        ":", fixed = TRUE))[(1:ntrain) * 2])
    pData(test.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(test.es)$LnRatio)), 
        ":", fixed = TRUE))[(1:ntest) * 2])
    
    # save the data to avoid future re-downloading
    save(dat, trainDat, testDat, train.es, test.es, file = "class_LNstatus.Rdata")
}
```


Now, we can do some exploratory analysis of the data before trying some classification methods:

```r
# undestand your data for classification
table(pData(train.es)$LnStatus)
```

```
## 
## neg pos 
##  48  48
```

```r
table(pData(test.es)$LnStatus)
```

```
## 
## neg pos 
##   9  11
```

```r
# understand the continuous response
tapply(pData(train.es)$LnRatio, pData(train.es)$LnStatus, summary)
```

```
## $neg
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0       0       0       0       0       0 
## 
## $pos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.040   0.070   0.110   0.194   0.228   0.960
```

```r
tapply(pData(test.es)$LnRatio, pData(test.es)$LnStatus, summary)
```

```
## $neg
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0       0       0       0       0       0 
## 
## $pos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.050   0.180   0.500   0.457   0.610   0.940
```

```r
# look at the expression of 3 randomly picked genes in both training and
# test sets
set.seed(1234)
(getMe <- sample(1:nrow(train.es), size = 3))  ## [1] 2756 15082 14766
```

```
## [1]  2756 15082 14766
```

```r
# training data
trDat <- trainDat[getMe, ]
str(trDat)
```

```
##  num [1:3, 1:96] 7.19 9.17 8.38 7.13 9.38 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:3] "201513_at" "223139_s_at" "222698_s_at"
##   ..$ : chr [1:96] "GSM570518" "GSM570519" "GSM570520" "GSM570521" ...
```

```r
trDat <- data.frame(LnStatus = pData(train.es)$LnStatus, Set = rep("train", 
    nrow(pData(train.es))), t(trDat))
str(trDat)
```

```
## 'data.frame':	96 obs. of  5 variables:
##  $ LnStatus    : Factor w/ 2 levels "neg","pos": 1 1 2 1 1 2 1 1 2 1 ...
##  $ Set         : Factor w/ 1 level "train": 1 1 1 1 1 1 1 1 1 1 ...
##  $ X201513_at  : num  7.19 7.13 7.39 6.86 6.96 ...
##  $ X223139_s_at: num  9.17 9.38 9.03 9.55 9.5 ...
##  $ X222698_s_at: num  8.38 8.24 7.23 7.87 8.45 ...
```

```r
plotDat.train <- melt(trDat, id = c("LnStatus", "Set"), variable_name = "gene")
colnames(plotDat.train)[colnames(plotDat.train) == "value"] = "gExp"

# test data
tDat <- testDat[getMe, ]
str(tDat)
```

```
##  num [1:3, 1:20] 6.05 9.15 7.55 6.87 8.95 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:3] "201513_at" "223139_s_at" "222698_s_at"
##   ..$ : chr [1:20] "GSM570498" "GSM570499" "GSM570500" "GSM570501" ...
```

```r
tDat <- data.frame(LnStatus = pData(test.es)$LnStatus, Set = rep("test", nrow(pData(test.es))), 
    t(tDat))
str(tDat)
```

```
## 'data.frame':	20 obs. of  5 variables:
##  $ LnStatus    : Factor w/ 2 levels "neg","pos": 1 1 1 1 1 1 1 1 1 2 ...
##  $ Set         : Factor w/ 1 level "test": 1 1 1 1 1 1 1 1 1 1 ...
##  $ X201513_at  : num  6.05 6.87 6.71 8 6.54 ...
##  $ X223139_s_at: num  9.15 8.95 9.09 9.81 9.2 ...
##  $ X222698_s_at: num  7.55 8.34 8.32 7.33 8.14 ...
```

```r
plotDat.test <- melt(tDat, id = c("LnStatus", "Set"), variable_name = "gene")
colnames(plotDat.test)[colnames(plotDat.test) == "value"] = "gExp"

plotDat <- rbind(plotDat.train, plotDat.test)
head(plotDat)
```

```
##   LnStatus   Set       gene  gExp
## 1      neg train X201513_at 7.192
## 2      neg train X201513_at 7.131
## 3      pos train X201513_at 7.395
## 4      neg train X201513_at 6.858
## 5      neg train X201513_at 6.955
## 6      pos train X201513_at 7.464
```

```r
# plot 3 randomly picked genes in both training and test sets using ggplot2
ggplot(plotDat, aes(LnStatus, gExp, colour = LnStatus)) + geom_point(alpha = 0.5) + 
    facet_wrap(~gene + Set)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 






### <dim id="4b">[Classification](#4b)

The prediction of a discrete response is usually refer to as *classification*. A response taking values over a finite set of labels is essentially the same thing as a factor. We will use the dataset from Smeets et al. to find the *best-trained* classifier and use it to predict the `LnStatus` of the 20 samples in the test set, i.e., classify those as "lymph node positive" or "negative".

#### <dim id="5b">[Feature and Model Selection](#5b)

We will first identify the best set of features that we will use to train the model using a cross-validation. Thus, I will divide the training set into 6 folds (the authors used 10 folds). We also want the proportion of positive and negative examples in each split to be approximately the same as for the full data set (i.e., stratified 6-fold CV with 8 positive and 8 negative samples within each fold). For each round of cross-validation, we use one fold as the test data and the rest of the data as training to select features and train different classifier.


#### <dim id="6b">[Cross validation splits](#6b)

This is not the only way to create splits of the training data to run a cross-validation. Note that if the samples can not be evenly divided into the nfolds you specified, then you need to complete the matrices below with NAs and call for entries different from NA at those folds.


```r
nfold <- 6
tabTrain <- table(train.es$LnStatus)
indlist <- sapply(names(tabTrain), function(z) which(train.es$LnStatus == z), 
    simplify = FALSE)
set.seed(1234)
# Each row contains 8 pos and 8 negative samples.
fold.pos <- matrix(sample(indlist[["pos"]]), nrow = nfold)
fold.neg <- matrix(sample(indlist[["neg"]]), nrow = nfold)
```


*Note*: with `CMA` you can use the command `GenerateLearningsets` to split the training data into folds. However, it does not show you how the data was split. Thus, you either use CMA for all or you write your own script.


```r
splits <- GenerateLearningsets(y = train.es$LnStatus, method = "CV", fold = 6, 
    strat = TRUE)
```



#### <dim id="7b">[Loop for feature selection and modeling](#7b)

To illustrate how to select a model, I will use the top-50 genes selected by `limma` (within each fold). Note that this number is very arbitrary and other options may make more sense like using a p-value threshold or testing different options with this CV. For this example, I'm using only the top-50 genes as methods like LDA and Logit can not be run on more features than samples. However, other methods like kNN or SVM will do well with more features.

In this example, I will compare 7 different models: kNN for k={1,5,10,15}, LDA, Logit, SVM. Feel free to add other methods to the list!


```r
# Define here the constants that you will not evaluate. For example, I will
# use the top-50 limma genes
ngenes <- 50
nmethod <- 7  #number of methods you plan to compare.

# Define here an output objects to store results
pr.err <- matrix(-1, nfold, nmethod, dimnames = list(paste0("Fold", 1:nfold), 
    c("1NN", "5NN", "10NN", "15NN", "LDA", "Logit", "SVM")))
for (i in 1:nfold) {
    
    # Test Fold for the i-th step
    testdat.fold <- trainDat[, c(fold.pos[i, ], fold.neg[i, ])]
    # I will create a factor of classes for the test set of the i_th fold
    testclass.fold <- train.es$LnStatus[c(fold.pos[i, ], fold.neg[i, ])]
    
    # The rest of the samples are the training set for the i-th step
    traindat.fold <- trainDat[, -c(fold.pos[i, ], fold.neg[i, ])]
    trainclass.fold <- train.es$LnStatus[-c(fold.pos[i, ], fold.neg[i, ])]
    
    # Step 1: feature selection (do you remember limma?)  Note that a different
    # set of genes will be selected for each fold! you can then compare how
    # consistent these sets were.
    limma.dat <- as.data.frame(traindat.fold)
    desMat <- model.matrix(~trainclass.fold, limma.dat)  # design matrix
    trainFit <- lmFit(limma.dat, desMat)
    eBtrainFit <- eBayes(trainFit)
    
    # top-50 limma genes
    top.fold <- topTable(eBtrainFit, coef = which(colnames(coef(trainFit)) != 
        "(Intercept)"), n = ngenes, sort.by = "P")
    
    # Retain the top-50 limma genes from the train and test sets
    traindat.fold <- traindat.fold[rownames(top.fold), ]
    testdat.fold <- testdat.fold[rownames(top.fold), ]
    
    # STEP 2: select a classifier Set a counter for the method tested
    l <- 0
    
    # kNN classifiers
    for (kk in c(1, 5, 10, 15)) {
        # every time you get inside this loop, the l counter gets redefined (i.e.,
        # 1, 2, etc for method 1, method 2, etc)
        l <- l + 1
        
        # knn needs samples in rows
        yhat.knn <- knn(train = t(traindat.fold), test = t(testdat.fold), cl = trainclass.fold, 
            k = kk)
        # Store the prediction error for each kk within this fold
        pr.err[i, l] <- mean(testclass.fold != yhat.knn)
    }  # end of kNN loop
    
    # LDA method. Note that you can change the prior parameter to reflect a
    # different proportion of case and control samples. The default is to use
    # the class proportions from the training set.
    
    m.lda <- lda(x = t(traindat.fold), group = trainclass.fold, prior = c(0.5, 
        0.5))
    yhat.lda <- predict(m.lda, newdata = t(testdat.fold))$class
    pr.err[i, "LDA"] <- mean(testclass.fold != yhat.lda)
    
    # Logit
    glm.dat <- data.frame(t(traindat.fold), group = trainclass.fold)
    m.log <- glm(group ~ ., data = glm.dat, family = binomial)
    pr.log <- predict(m.log, newdata = data.frame(t(testdat.fold)), type = "response")
    pr.cl <- rep(0, length(testclass.fold))
    pr.cl[pr.log > 1/2] <- "pos"
    pr.cl[pr.log <= 1/2] <- "neg"
    pr.cl <- factor(pr.cl)
    pr.err[i, "Logit"] <- mean(pr.cl != testclass.fold)
    
    # SVM
    m.svm <- svm(x = t(traindat.fold), y = trainclass.fold, cost = 1, type = "C-classification", 
        kernel = "linear")
    pr.svm <- predict(m.svm, newdata = t(testdat.fold))
    
    pr.err[i, "SVM"] <- mean(pr.svm != testclass.fold)
}  # end of CV loop
```

```
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```r
## JB: I get 'There were 12 warnings'
```



#### <dim id="8b">[Error Rates](#8b)

Now you can get the average prediction error for all methods. Note that the prediction errors are high! not too much hope for the real test run!


```r
cv.err <- colMeans(pr.err)

# mean - 1 sd (sd of the 6 error rates)
ls <- cv.err - apply(pr.err, 2, sd)

# mean + 1 sd (sd of the 6 error rates)
us <- cv.err + apply(pr.err, 2, sd)

# plot the results
plot(1:nmethod, cv.err, ylim = c(0, 1), xlim = c(1, (nmethod + 0.5)), type = "n", 
    axes = FALSE, xlab = "Classifier", ylab = "Error rate", main = "6-fold CV Error")

for (j in 1:ncol(pr.err)) points(jitter(rep(j, 6), factor = 2), jitter(pr.err[, 
    j]), cex = 0.8, pch = "X", col = "gray")

for (i in 1:nmethod) lines(c(i, i), c(ls[i], us[i]), lwd = 2, col = "gray")
points(1:nmethod, ls, pch = 19, col = "red")
points(1:nmethod, us, pch = 19, col = "green")
points(1:nmethod, cv.err, pch = 19, cex = 1.5, col = "black")
axis(2, ylab = "Error rate")
axis(1, 1:nmethod, colnames(pr.err))

box()
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



### <dim id="9b">[Results of the CV](#9b)

According to these results, 1NN and 10NN may be the better classifier to try in the test data. However, remember that this CV results depend on the first split of the data we did. Thus, we need to repeat this CV.

**Exercise 1**: perform 100 runs of this CV before selecting a model to test! Add at least on model to the list of models, e.g., use genes with a p-val threshold < cutoff.

**Exercise 2**: Use AUC as a criteria to select a model based on the training data! Tip: extract the predicted probabilities from each method and use the roc function in ROCR.


### <dim id="10b">[Testing the selected model](#10b)

Now that we decided on which method we are going to use to classify samples in the test set, we need to train the model using the *FULL* training set and then classify samples of the test set. I will use the 10NN model.


```r
yhat.knn <- knn(train = t(trainDat), test = t(testDat), cl = train.es$LnStatus, 
    k = 10)
# Store the prediction error for each kk within this fold
pr.errTest <- mean(test.es$LnStatus != yhat.knn)
pr.errTest
```


Not good! In real practice, you should not keep trying until we get a good result! However, in this seminar, I encourage you to try different options as an exercise and to see how much the results can change.

### <dim id="11b">[CMA](#11b)
Many steps of the CV defined above can be easily done with CMA. For example, Step 1 in the loop above can also be done using 'CMA' with the function 'GeneSelection', which selects the most informative features (e.g., gene) to build a classifier within each of the splits generated by 'GenerateLearningsets'. Some learning algorithms do better if you only give them "useful" features.


```r
featureScores <- GeneSelection(X = t(trainDat), y = train.es$LnStatus, learningsets = splits, 
    method = "limma")

# Compare list of selected genes using:
toplist(featureScores)
# We can aggregate the results across the 6 splits.
seliter <- numeric()
for (i in 1:nfold) seliter <- c(seliter, toplist(featureScores, iter = i, top = 10, 
    show = FALSE)$index)
sort(table(seliter), dec = TRUE)  # summarize
# Choose the 20 probes which are chosen most commonly in the 6 splits
bestprobes <- as.numeric(names(sort(table(seliter), dec = TRUE)))[1:20]

# examine the annotations. I just selected a few columns from the fData of
# the eSet.
fData(dat)[bestprobes, c("Gene Symbol", "Gene Title", "ENTREZ_GENE_ID", "Representative Public ID")]
```


This looks promising since I get TTC3 and at least a couple of other genes that show up on Table 3 of the paper.

Similarly, you can use CMA to train and test a classifier within each CV fold (learningsets). However, there are things you can not do within CMA or that CMA is not doing right. For example, CMA can not do a full nested cross-validation. Additionally, it is not trivial to train the selected in the full dataset and then test it in the test set. CMA is more designed for CV. Thus, it is good to know how to do this things by hand as well.

Paul solved this problem in the following way: he made a `learningsets` object that has just one "split" defined by the samples in the training set.


```r
m <- matrix(which(dat$Set == "patient type: training set"), 1)

full.learningset <- new("learningsets", learnmatrix = m, method = "my own", 
    ntrain = 96, iter = 1)

fullFeatureScores <- GeneSelection(X = t(exprs(dat)), learningsets = full.learningset, 
    y = dat$LnStatus, method = "t.test")

testclassif <- classification(X = t(exprs(dat)), y = dat$LnStatus, learningsets = full.learningset, 
    genesel = fullFeatureScores, nbgene = 100, classifier = pknnCMA, k = 5)

# Evaluation:
tres <- testclassif[[1]]
ftable(tres)
roc(tres)
```


Note: This optimized classifier did terribly as well.

**Discussion points**  
How good do you think a classifier will have to be to be clinically relevant? What level of specificity or sensitivity do you think is "good enough"? The author's data set has half lymph-node positive and half negative. The incidence of lymph-node-positive breast cancer is about 33% at the time of diagnosis (according to [http://seer.cancer.gov/statfacts/html/breast.html]). How does that affect your opinion about the classifier?

**References**  
[Ref1]: Varma S, Simon R: Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics 2006, 7:91.

[Ref2]: Statnikov A, Aliferis CF, Tsamardinos I, Hardin D, Levy S: A comprehensive evaluation of multicategory classification methods for microarray gene expression cancer diagnosis. Bioinformatics 2005, 21:631-643.
