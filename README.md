# Bayesian Network-based Mendelian Randomization (BNMR)
This is an R package to conduct causal estimation between exposure and outcome on GWAS data using BNMR model, which is a two-staged framework to deal with the imperfect genetic instruments. It selects approaciate genetic instruments via random graph forest, an ensemble approach comprises a series Bayesian network structure learning within sampled variants and exposure. It then uses the Bayesian Mendelian randomization with a shrinkage prior to cope with horizontal pleiotropy and obtain a robust estimate. 

![image](https://github.com/sjl-sjtu/bnmr/blob/main/FIG/Fig1.jpg)

## Tutorial
### 1. Installation
You can install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/bnmr")
```
### 2. Beginners' Guide
Here we provide a illustrational step-by-step example with a simulated dataset to demonstrate the usage of `bnmr`.

#### 1) Contructing simulated dataset
Let's start with a small simulated dataset.
```R
n <- 2000
p <- 200
snps <- replicate(p,sample(0:2,n,replace = TRUE))
snps <- apply(snps,2,as.numeric)
snpname <- paste0("g",1:p)
df <- as.data.frame(snps)
colnames(df) <- snpname
truesnp <- paste0("g",sample(1:p,50))
df$x <- as.matrix(df[,truesnp])%*%rnorm(50,0.05,0.05)+rnorm(n,0,1)
df$y <- 0.5*df$x+rnorm(n,0,1)
```

#### 2) Learning stage: RGF for IV selection
The example dataset includes 2,000 samples and 200 SNP loci, as well as two phenotypes, `x` and `y`. We now use a random graph forest (RGF) to select the loci that directly affect the exposed x. This can be done by using the `bn` function of package `bnmr`.
```R
library(bnmr)
rgf_results <- bn(df,snpname,"x",bn_method="hc",repeats=1000,alpha=0.9,nsam=2000,psam=100)
```

The object `rgf_results` is a list containing two objects. The first object, `selectsnp`, is the names of genetic instrumental variables (IVs) selected based on predefined parameters, and the second object, `dfscore`, is the complete RGF result, which is a data frame containing two columns - the first column is the name of the SNP, and the second column is the corresponding adjacency score, arranged in descending order of adjacency scores. 

We provide two criteria for selecting IVs from adjacency scores, the first is to specify the number of IVs directly, which can be achieved by specifying the `selectNum` parameter in the `bn` function; The second is to set a threshold $\alpha$ between 0-1, at which point SNPs with a critical score greater than `alpha*psam/p` (`p` is the total number of SNPs assessed by RGF and `psam` is the number of SNPs selected per sampling) will be selected as IVs, which can be achieved by specifying parameter `alpha` in the `bn` function. Note that the parameter `alpha` will not work when the parameter `selectNum` is specified. The default `alpha` of the program is 0.9. Of course, you can also use the data frame containing full adjacency score of the output to define your own criteria to select IVs.

#### 3) Inference stage: Bayesian estimation
After selecting IVs, the Bayesian Mendelian randomization model can be used to estimate the causal parameters of exposure `x` to outcome `y`, which can be achieved by the function `mr` in the `bnmr` package.
```R
IVs <- rgf_results$selectsnp
mr_results <- mr(df,IVs,"x","y",mr_model="linear",prior="horseshoe",n.iter=5000,n.chain=4)
# show results
print(mr_results$mean,mr_results$se,mr_results$lower,mr_results$upper)
```

We implemented linear (for quantitative outcome, by specifying `mr_model="linear"`) and logistic (for binary outcome, by specifying `mr_model="logit"`) model in Mendelian randomization. Note that for hypothesis testing purposes, a linear model can also be used for dichotomous outcomes. In general, linear models are more efficiency.

The Bayesian estimation is implemented based on `RStan`. The object `mr_results` is a list containing seven elements. The first element `betaList` is a vector containing posterior sampling of causal effect $\beta$. The following four elements are the posterior mean estimate, standard error, the lower and upper bound of 95% credible intervals of $\beta$. The element `Rhat` is an indicator of convergence (Rhat < 1.1 at convergence). `fit_detail` is an S4 class Stan fit object containing details of Bayesian estimation.

The above two steps can be integrated via function `bnmr`.
```R
model <- bnmr(df,snpname,"x","y",bn_method="hc",repeats=1000,alpha=0.9,nsam=2000,psam=100,
              mr_model="linear",prior="horseshoe",n.iter=5000,n.chain=4)
```

### 3. Hyperparameters
TBD

### 4. API
Deatil usage and examples can be found at https://github.com/sjl-sjtu/bnmr/blob/main/bnmr_0.2.1.pdf.

### 5. Adaptation to large-scale biobank-level data
For large-scale dataset (like biobank), we recommend to conduct Bayesian MR analysis using Python Package PyMC with NUTS JAX samplers (NumPyro or BlackJAX) and GPU (https://www.pymc-labs.com/blog-posts/pymc-stan-benchmark/) to achieve faster posterior sampling. A tutorial of PyMC (v5) with JAX and Numba can be found at https://www.pymc.io/projects/examples/en/latest/samplers/fast_sampling_with_jax_and_numba.html. An example in BNMR can be found at https://github.com/sjl-sjtu/bnmr/blob/main/BayesianMR_example_pymc.py.


Supplementary notes, tables, and figures for the paper can be found at https://github.com/sjl-sjtu/bnmr/blob/main/supplementary_notes.pdf.

Latest updation at Nov. 11th, 2023.

Contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)
