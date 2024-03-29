# Bayesian Network-based Mendelian Randomization (BNMR)
This is an R package to conduct causal estimation between exposure and outcome on GWAS data using the Bayesian Network-based Mendelian Randomization (BNMR), which is a two-stage framework to deal with imperfect genetic instruments (horizontal pleiotropy, linkage disequilibrium, epistasis, etc.). It selects variants with a direct effect on the exposure as instrumental variables via random graph forest (RGF), an ensemble approach comprised of a series of Bayesian network structure learning processes within sampled variants and exposure. It then uses the Bayesian Mendelian randomization (BMR) with a shrinkage prior on the nuisance parameters to cope with potential horizontal pleiotropy and obtain a robust estimate. 

Current version: 0.3.1

Latest updation at Feb. 6th, 2024. We have renamed the package to `bnmr2` to avoid potential conflict with another previous package (in Econometrics).

## Tutorial
### 1. Installation
You can install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/bnmr2")
```
Check the installation
```R
library(bnmr2)
```

### 2. Beginners' Guide
Here we provide a illustrational step-by-step example with a simulated dataset to demonstrate the usage of `bnmr2`.

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
The example dataset includes 2,000 samples and 200 SNP loci, as well as two phenotypes, `x` and `y`. We now use a random graph forest (RGF) to select the loci that directly affect the exposed x. This can be done by using the `bn` function of package `bnmr2`.
```R
library(bnmr2)
rgf_results <- bn(df,snpname,"x",bn_method="hc",repeats=1000,alpha=0.5,nsam=2000,psam=100)
```

The object `rgf_results` is a list containing two objects. The first object, `selectsnp`, is the names of genetic instrumental variables (IVs) selected based on predefined parameters, and the second object, `dfscore`, is the complete RGF result, which is a data frame containing two columns - the first column is the name of the SNP, and the second column is the corresponding adjacency score, arranged in descending order of adjacency scores. 

We provide two criteria for selecting IVs from adjacency scores, the first is to specify the number of IVs directly, which can be achieved by specifying the `selectNum` parameter in the `bn` function; The second is to set a threshold $\alpha$ between 0-1, at which point SNPs with a critical score greater than `alpha*psam/p` (`p` is the total number of SNPs assessed by RGF and `psam` is the number of SNPs selected per sampling) will be selected as IVs, which can be achieved by specifying parameter `alpha` in the `bn` function. Note that the parameter `alpha` will not work when the parameter `selectNum` is specified. The default `alpha` of the program is 0.5. Of course, you can also use the data frame containing full adjacency score of the output to define your own criteria to select IVs. Please note that the selection of IVs mentioned here is closely related to the inference stage, so it is necessary to adjust the threshold and select an appropriate number of IVs based on the actual situation.

#### 3) Inference stage: Bayesian estimation
After selecting IVs, the Bayesian Mendelian randomization model with a shrinkage prior on the instruments' horizontal pleiotropic effects can be used to estimate the causal parameters of exposure `x` to outcome `y`, which can be achieved by the function `mr` in the `bnmr2` package.
```R
IVs <- rgf_results$selectsnp
mr_results <- mr(df,IVs,"x","y",mr_model="linear",prior="horseshoe",n.iter=5000,n.chain=4)
# show results
print(c(mr_results$mean,mr_results$se,mr_results$lower,mr_results$upper))
```

We implemented linear (for quantitative outcome, by specifying `mr_model="linear"`) and logistic (for binary outcome, by specifying `mr_model="logit"`) model in Mendelian randomization. Note that for hypothesis testing purposes, a linear model can also be used for binary outcomes like diseases. In general, linear models are more efficiency.

The Bayesian estimation is implemented based on `RStan` (<https://mc-stan.org/users/interfaces/rstan>). The object `mr_results` is a list containing seven elements. The first element `betaList` is a vector containing posterior sampling of causal effect $\beta$. The following four elements `mean`, `se`, `lower`, and `upper` are the posterior mean estimate, standard error (standard deviation of posterior samples of the parameter), the lower and upper bound of 95% credible intervals of $\beta$, respectively. The element `Rhat` is an indicator of convergence (Rhat < 1.1 at convergence). `fit_detail` is an S4 class Stan fit object containing details of Bayesian estimation.

The above two steps can be integrated via function `bnmr2`.
```R
model <- bnmr(df,snpname,"x","y",bn_method="hc",repeats=1000,alpha=0.9,nsam=2000,psam=100,
              mr_model="linear",prior="horseshoe",n.iter=5000,n.chain=4)
```

### 3. Hyperparameters
#### 1) Learning stage:
There are some hyperparameters in RGF that can be specified, including the number of samples per sampling `nsam`, the number of loci selected for each samplings `psam`, the number of sampling in ensemble learning `repeats`, and the algorithm for Bayesian Network (BN) structure learning `bn_method`.
* `nsam`: For small samples (<5000), it is recommended to set nsam to the same as the total sample size and set `sample_replace=TRUE`, which is bootstrapping. For larger samples, nsam recommends setting between 2000-5000, `sample_replace` can be set to `TRUE` or `FALSE`, corresponding to sampling from all samples with and without return, respectively.
* `psam`: We recommend setting the number of variants to be selected between 100 and 150 per sampling. An increase in the number of variants significantly increases the calculation time.
* `repeats`: To ensure adequate sampling, we recommend that both `repeats*nsam/n` and `repeats*psam/p` should be above 100. Here n and p are the total number of samples and the number of loci, respectively.
* `bn_method`: We implement BN learning with algorithms provided by R package `bnlearn` (<https://www.bnlearn.com/>), including contraint-based, score-based, and hybrid approaches. In our simulations, fast.iamb performs best in terms of false discovery rate (FDR) control, while the Hill-Climbing (hc) algorithm is the fastest. Although other methods may be more precise, hc is able to meet the requirements under mundane situations.

#### 2) Inference stage:
The main hyperparameters of the inference phase include the shrinkage prior `prior`, number of IVs, and number of iterations `n.iter`.
* IV number: The number of IVs depends on the results of RGF learning on the one hand, and on the other hand, on the screening criteria based on adjacency scores that we have set artificially (two selection criteria were discussed in Part 2: specifying the number `selectNum` and setting the threshold `alpha` in `bn` function). Simulations manifest that the variance decreases to a stable level as the IV increases, while bias increases when there are too many or too few instruments. We recommend using 20-50 IVs, depending on the actual situation.
* Prior: We implement five types of shrinkage prior (horseshoe, Bayesian lasso, hyperlasso, Bernoulli Spike & Slab, and Uniform Spike & Slab) by specifying parameter `prior` to one of `"horseshoe"`, `"lasso"`, `"hyperlasso"`, `"spikeslabBernoulli"`, and `"spikeslabUniform"`. BNMR estimates are not sensitive to priors in simulations. The Bayesian Lasso prior shows the fastest sampling speed and a deflated standard error, but has a slightly higher bias. Bernoulli spike and slab prior tends to achieve a small standard error. The horseshoe prior, although slightly less efficient, is superior in the performance of convergence due to the lowest Rhat. We recommend to use horseshoe or Bernoulli spike & slab prior if no additional information avaliable.
* Numbers of iterations `n.iter` and chains `n.chain`: We endorse the use of at least 4 chains and at least 2000 iterations per chain in MCMC.

### 4. API
Detail usage and examples can be found at [https://github.com/sjl-sjtu/bnmr/blob/main/bnmr2_0.3.1.pdf](https://github.com/sjl-sjtu/bnmr2/blob/main/bnmr2_0.3.1.pdf).

### 5. Adaptation to large-scale biobank-level data
For Bayesian posterior sampling within large-scale dataset (like biobank), we provide two options:
#### 1) Consolidation of subset posterior sampling
Split the whole dataset into samll subsets and conduct MCMC sampling parallelly and combine the posterior. This can be achieved with function `mr_split` in package `bnmr2`.
```R
mr_results <- mr_split(df,truesnp,"x","y",mr_model="linear",prior="horseshoe",
                       n.iter=5000,n.chain=4,n.split=4)
```
Here we use `n.split=4` to represent that the entire dataset is divided into 4 subsets, sampled separately by MCMC, and then merged the posterior sampling within each subset.

#### 2) Using PyMC with NUTS JAX samplers (Recommended)
We recommend to conduct Bayesian MR analysis within large-scale biobank using Python Package PyMC with NUTS JAX samplers (NumPyro or BlackJAX) and GPU (<https://www.pymc-labs.com/blog-posts/pymc-stan-benchmark/>) to achieve faster posterior sampling. A tutorial of PyMC (v5) with JAX and Numba can be found at <https://www.pymc.io/projects/examples/en/latest/samplers/fast_sampling_with_jax_and_numba.html>. An example in BNMR can be found at https://github.com/sjl-sjtu/bnmr/blob/main/BayesianMR_example_pymc.py.

### * An integrated example using bnmr and PyMC to handle large scale data from UK Biobank
The tutorial for the integrated example is shown at <https://github.com/sjl-sjtu/bnmr2/blob/main/integrated_example.md>.

## Reference
Jianle Sun, et al. Bayesian network-based Mendelian randomization for variant prioritization and phenotypic causal inference, *Human Genetics*, accepted.

Other reference links:
* bnlearn: Bayesian network structure learning (<https://www.bnlearn.com/>)
* stan: statistical modeling and high-performance statistical computation (<https://mc-stan.org/>)
* PyMC: probabilistic programming library for Python (<https://www.pymc.io/welcome.html>)
* Bayesian MR: <https://github.com/carloberzuini/BMR>
* shrinkage priors for Bayesian penalized regression: <https://github.com/sara-vanerp/bayesreg>

Contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)

![image](https://github.com/sjl-sjtu/bnmr/blob/main/FIG/Fig1.jpg)
