# An intergrate exemple of BNMR 

### 1. Data Preparation
Firstly, according to the GWAS summary statistics and the set primary screening criteria, we can prepare a file of sample ID and corresponding phenotypes, and another with list of pre-filtered candidate SNPs. Genomic analysis tools like PLINK (https://www.cog-genomics.org/plink/2.0/) can then be used to extract the data of the pre-screened samples and SNPs. Then, we use the R-package bnmr for ensemble Bayesian network structure learning to select tool variables.

With files converted by PLINK, we need to organize them to prepare for BN learning.
```R
library(tidyverse)
library(data.table)
# prepare data
list_of_files <- list.files(pattern="RBC_chr([1-9]|1[0-9]|2[0-2]).raw",full.names=FALSE)  # files containing dosage (0/1/2) information for each variant
df <- map(.x=set_names(list_of_files), .f=read_delim)
df <- purrr::reduce(df,left_join,by=c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
snps <- grep("^rs",colnames(df),value = TRUE)
df[,snps] <- round(df[,snps])
dftrait <- read_delim("phenotypes.txt")  # files containing the phenotypes you need
df <- df %>% left_join(dftrait,by=c("FID","IID"))
df %>% write_csv("RBC_bind.csv")
```

### 2. Learning Stage
Apply RGF to the datasets and attain the adjacency scores for each loci.
```R
library(bnmr)
ns <- 4000
ps <- 150
r <- 5000

snps <- grep("^rs",colnames(df),value = TRUE)
dfre <- bn(df1,snps,"RBC",bn_method="hc",repeats=r,nsam=ns,psam=ps)
dfs <- dfre$score
dfs %>% write_csv("RBC_score.csv")
```

### 3. Inference Stage
Since this dataset includes more than 200,000 individuals, we used Python library PyMC (https://www.pymc.io/welcome.html) and NUTS JAX samplers to estimate the inference stage. First let me set the environment
```python
import os
import multiprocessing
```
If you are working on CPU and want to multi-core parallel
```python
# if you use cpu
os.environ["XLA_FLAGS"] = "--xla_force_host_platform_device_count={}".format(
    multiprocessing.cpu_count()
)
```
If you want to work on GPU, specify the number of cudas. For example, we can use 2 cudas in computation by specifying `os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"`.
```python
# if you use gpu
os.environ["THEANO_FLAGS"] = "device=cuda"
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"
```

Import neccessary libraries now. Make sure you have PyMC 5.6 and above installed in your environment, as well as at least one of numpyro, blackjax, and nutpie (corresponding to the three types of samplers mentioned below).
```python
import jax
import torch
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
RANDOM_SEED = 0
rng = np.random.default_rng(RANDOM_SEED)
np.set_printoptions(2)
import pymc as pm
import aesara
import numpyro
```

Now load the dataset we want to analysis. Here we need to load two datasets: one containing the individual-level data of exposure, outcome, and genetic instrumental variables (`df` here), and another containing the adjacency scores of each SNP after pre-filtering (`df_score` here). Both datasets have been processed before.
```python
df = pd.read_csv("RBC_bind.csv") # individual-level data of SNPs, exposures and outcomes
df_score = pd.read_csv("RBC_score.csv") # adjacency score of SNP-RBC from random graph forest
df_score = df_score.sort_values(by='score', ascending=False)
IV = df_score['snp'].head(20).tolist()
feature = ['RBC','DBP']
df = df.loc[:, IV + feature].dropna()

exposureName = "RBC"
outcomeName = "DBP"
exposure = df[exposureName].values
outcome = df[outcomeName].values
s = IV
```

Prepare data for PyMC
```python
# Data
N = df.shape[0]
J = len(s)
X = np.array(exposure)
Y = np.array(outcome)
Z = df[s].values.reshape((N, J))
```

Define model and priors (here we show the example of horseshoe prior)
```python
# Define the PyMC3 model
with pm.Model() as shrink_model:
    # Priors
    sigmax = pm.HalfFlat('sigmax') 
    sigmay = pm.HalfFlat('sigmay') 
    sigmaalpha = pm.HalfFlat('sigmaalpha')
    mualpha = pm.Flat('mualpha')
    omegax = pm.Flat('omegax')
    omegay = pm.Flat('omegay')
    deltax = pm.Flat('deltax')
    deltay = pm.Flat('deltay')
    beta = pm.Flat('beta')
    tau = pm.HalfCauchy('tau',beta=1)
    phi = pm.HalfCauchy('phi',beta=1,shape=J)
    gamma = pm.Cauchy('gamma', alpha=0, beta=phi*tau, shape=J)
    alpha = pm.Normal('alpha', mu=mualpha, sigma=sigmaalpha, shape=J)
    u = pm.Normal('u', mu=0, sigma=1, shape=N)
    
    # Likelihoods
    X = pm.Normal('X', mu=omegax + pm.math.dot(Z, alpha) + u * deltax, sigma=sigmax, observed=X)
    Y = pm.Normal('Y', mu=omegay + pm.math.dot(Z, gamma) + X * beta + u * deltay, sigma=sigmay, observed=Y) # for quantitative outcome
    # Y = pm.Bernoulli('Y', p=pm.invlogit(omegay + pm.math.dot(Z, gamma) + X * beta + u * deltay), observed=Y)  # for binary outcome
```

Then sample the posterior distribution using MCMC. Different samplers can be used, like NumPyro JAX NUTS sampler (by specifying `nuts_sampler="numpyro"`, BlackJAX NUTS sampler (by specifying `nuts_sampler="blackjax"`), and Nutpie Rust NUTS sampler (by specifying `nuts_sampler="nutpie"`). We recommend to use NumPyro JAX NUTS sampler here.
```python
with shrink_model:    
    trace = pm.sample(draws=5000, tune=5000, chains=4, cores=4, target_accept=0.9,nuts_sampler="numpyro")
    # Get the samples
    subdata = trace.posterior['beta']
    print(az.summary(subdata))
```
