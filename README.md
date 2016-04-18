## Integrating multiple Omics for Prediction of BC Survival
The following scripts illustrate how to fit some of the models presented in [Vazquez et al., Genetics (2016)]().

#### (1) Installing required library
The code below illustrates how to install and load the necessary packages from CRAN using `install.packages()`.
```R
 install.packages(pkg='BGLR')    #1# install BGLR
 install.packages(pkg='pROC')    #2# install pROC
 library(BGLR); 
 library(pROC);
 ```   

#### (2) Loading data.
Data. The code assumes that the user has saved in the file OMIC_DATA.rda the objects that contain the phenotypic, covariates and omic information
•	XF: an incidence matrix with clinical covariates.
•	Xge: an incidence matrix with gene expression. 
•	Xmt: an incidence matrix with methylation values at various sites
•	y: a vector with responses, it may be continuous (time to event type) or binary (0/1), we are assuming a binary response.
It is assumed that the omics data was QC edited removing outliers, constant variables, or columns with too many missing values (e.g. more than 10% missing values). It is also assumed that there are no missing values in the predictors, thus, if there is they should be imputed. 
```R
 load('OMIC_DATA.rda')
 ```   
 
#### (3) Similarity matrices.
 TCGA sata is small enough to be loaded and manipulated in RAM. Here we can computes a similarity matrix of the form G=XX' directly as follows:
 ```R 
 #Building Gene Expression Similarity Matrix
#Center the columns of Xge
 Xge<- scale(Xge, scale=FALSE, center=TRUE)
 Gge<-tcrossprod(Xge);tmp<-mean(diag(Gge)); Gge<-Gge/tmp
 XF<- scale(XF, scale=FALSE, center=TRUE) # centering and scaling the incidence matrix for fixed effects.
```
 
Alternatively, when data sets have higher dimensions and XX' cannot be computed in RAM, the [BGData] R-package  (https://github.com/quantgen/BGData) has an option to compute G in parallel. The function `getG()` from the  computes a similarity matrix. The function has several parameters to manage optional centering and scaling. The function can use multi-core computing. For further details follow the link provided above. With this packages the similarities could be fitted as follows:
```R
 library(BGData)
 Gge<-getG(Xge,scaleCol=T,scaleG=T) # Similarity matrix for gene expression.
 Gmt<-getG(Xmt,scaleCol=T,scaleG=T) # Similarity matrix for methylation. 
```

#### (4)  Fitting a survival model for Fixed effects using BGLR
The following code illustrates how to use BGLR to fit a fixed effects model. The matrix XF is an incidence matrix for effects. There is no column for intercept in XF because BGLR adds the intercept automatically. The response variable `y` is assumed to be coded with two lables (e.g., 0/1), the argument `response_type` is used to indicate to BGLR that the response is ordinal (the binary case is a special case with only two levels). Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Pérez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf) The code also shows how to retrieve estimates of effects and of psuccess probabilities.

```R
 # Inputs
  ETA<-list( COV=list(X=XF, model='FIXED') )
 
 # Fitting the model
  fm=BGLR(y=y, ETA=LP, saveAt='cov_', response_type='ordinal')
 
 # Retrieving estimates
  fm$ETA$COV$b      # posterior means of fixed effects
  fm$ETA$COV$SD.b   # posteriro SD of fixed effects
  head(fm$probs)    # estimated probabilities for the 0/1 outcomes.
```

#### (1) Installing required libraries


The code below illustrates how to install BGLR and BGData from GitHub. BGLR can also be installed from CRAN using `install.packages()`.

```R
 install.packages(pkg='devtools',repos='https://cran.r-project.org/')    #1# install devtools
 library(devtools)                                                       #2# load the library
 install_git('https://github.com/quantgen/BGData/')                      #4# install BGLR from GitHub
```   
