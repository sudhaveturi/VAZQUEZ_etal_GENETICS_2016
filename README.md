## Integrating multiple Omics for Prediction of BC Survival
The following scripts illustrate how to fit some of the models presented in *Vazquez et al., Genetics, 2016*.

#### (1) Installing BGLR
The code below illustrates how to install and load the necessary package from CRAN using `install.packages()`.
```R
   install.packages(pkg='BGLR')    #1# install BGLR
   library(BGLR); 
 ```   

#### (2) Loading data

**Data**. The code assumes that the user has saved in the file `OMIC_DATA.rda` the objects that contain the phenotypic information, clinical covariates, and omic data. The code assumes that the file `OMIC_DATA.rda` contain the following objects:
   * `XF`: an incidence matrix for clinical covariates.
   * `Xge`: an incidence matrix for gene expression. 
   * `Xmt`: an incidence matrix for methylation values at various sites.
   * `y`: a vector with the response, in this case a 0/1 where 0 denotes alive.

The code below assumes that all the predictors were edited by removing outliers and predictors that did not vary in the sample and transformed and imputed (i.e., no NAs in predictors) if needed.

 
#### (2) Computing similarity matrices
 Some of the models fitted use a similarity matrix of the form G=XX' computed from omics. The following code illustrates how to compute these matrices.
 
 ```R 
  load('OMIC_DATA.rda')

  #Building Gene Expression Similarity Matrix
   Xge<- scale(Xge, scale=FALSE, center=TRUE) #centering and scaling
   Gge<-tcrossprod(Xge)                       #computing crossproductcts
   Gge<-Gge/mean(diag(Gge)                    #scales to an average diagonal value of 1.
```
 
NOTE: for larger data sets it may be more convinient to use the `geG()` function of the [BGData](https://github.com/quantgen/BGData) R-package. This function allows computing G without loading all the data in RAM and offers methods for multi-core computing. 


#### (3)  Fitting a binary regression for (the "fixed effects" of) Clinical Coavariates using BGLR

The following code illustrates how to use BGLR to fit a fixed effects model. The matrix XF is an incidence matrix for clinical covariates. There is no column for intercept in XF because BGLR adds the intercept automatically. The response variable `y` is assumed to be coded with two lables (e.g., 0/1), the argument `response_type` is used to indicate to BGLR that the response is ordinal (the binary case is a special case with only two levels). Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Pérez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf). The code also shows how to retrieve estimates of effects and of success probabilities. In thr example we use 12000 iterations with 2000 iterations discarded for burn-in. 

```R
# Inputs
 XF<- scale(XF, scale=FALSE, center=TRUE) # centering and scaling the incidence matrix for fixed effects.
 ETA.COV<-list( COV=list(X=XF, model='FIXED') )
# Fitting the model
 #(1)  
 fm=BGLR(y=y, ETA=ETA.COV, saveAt='cov_', response_type='ordinal')
# Retrieving estimates
 fm$ETA$COV$b      # posterior means of fixed effects
 fm$ETA$COV$SD.b   # posteriro SD of fixed effects
 head(fm$probs)    # estimated probabilities for the 0/1 outcomes.
 
```

#### (5)  Fitting a binary model for fixed effects and whole genome gene expression (GE) using BGLR

The following code illustrates how to use BGLR to fit a mixed effects model that accomodates both clinical covariates and whole-genome-gene expression.

```R
# Setting up parameters to adjust COV + one omic model
  ETA.COV.GE<-list( list(X=XF, model='FIXED'), list(K=Gge, model='RKHS'))
# Fitting the model
  #(2)
  fm.COV.GE<- BGLR(y=y, ETA=ETA.COV.GE, response_type='ordinal',nIter=)
```

#### (6)  Fitting a binary model for fixed effects covariates and 2 whole genome omics.
The following code shows how to extend the second model (2) to incorporation of two omics. Because the model is more complex that the covariates-only model we use longer MCMC chains.

```R
#Building Methylation Similarity Matrix
#Center the columns of Xmt
Xmt<- scale(Xmt, scale=FALSE, center=TRUE)
Gmt<-tcrossprod(Xmt);tmp<-mean(diag(Gmt)); Gmt<-Gmt/tmp
# Setting up parameters to adjust the COV model + 2 omics
ETA.COV.GE.MT<-list( list(X=XF, model='FIXED'),
           list(K=Gge, model='RKHS'),
           list(K=Gmt, model='RKHS'))

# Fitting models 
# (3) 
fm.COV.GE.MT<- BGLR(y=y, ETA=ETA.COV.GE.MT, 
                 response_type='ordinal', nIter=55000,burnIn=5000)
```

#### (7) Below we demostrate how to test the model (1) in a randomly selected testing set.
```R
#Installing and loading library pROC to compute Area Under the ROC Curve.
install.packages(pkg='pROC')    # install pROC
library(pROC);
 
#  Randomly select a 20% of the data to be the a separated testing set:
n <- length(y)
tst<- runif(n) <0.2
yNA = y
yNA[tst] <-NA

# Fit the model only in the training set
fm.COVtr<- BGLR(y=yNA, ETA=ETA.COV, response_type='ordinal')
# Find probability of survival for the testing set
pred <-fm.COVtr$probs[tst,2]

# Estimate AUC in training and testing sets
AUC_train<-auc(y[-tst],fm.COV2$yHat[-tst])
AUC_test<-auc(y[tst], pred)

#For the first individual, area under the standard normal curve (CDF) of estimated y from full model:
pnorm(fm.COV$yHat[1])
```
