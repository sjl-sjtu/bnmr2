import os
import multiprocessing

# if you use cpu
os.environ["XLA_FLAGS"] = "--xla_force_host_platform_device_count={}".format(
    multiprocessing.cpu_count()
)

# if you use gpu
# os.environ["THEANO_FLAGS"] = "device=cuda"
# os.environ["CUDA_VISIBLE_DEVICES"] = "0,1,2,3"

import jax
print(jax.default_backend())
print(jax.devices())
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
import aesara.tensor as at 
import aesara
import numpyro

df = pd.read_csv("lym_bind.csv") # data frame of SNPs, exposures and outcomes
df1 = pd.read_csv("lym_score.csv") # adjacency score from RGF (BN)
df1 = df1.sort_values(by='score', ascending=False)
IV = df1['snp'].head(20).tolist()
feature = ['RBC','DBP']

df = df.loc[:, IV + feature].dropna()

exposureName = "RBC"
outcomeName = "DBP"
exposure = df[exposureName].values
outcome = df[outcomeName].values
s = IV

# Data
N = df.shape[0]
J = len(s)
X = np.array(exposure)
Y = np.array(outcome)
Z = df[s].values.reshape((N, J))
# X = at.as_tensor_variable(X)
# Y = at.as_tensor_variable(Y)
# Z = at.as_tensor_variable(Z)


# Define the PyMC3 model (here show the example of horseshoe prior)
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
    Y = pm.Normal('Y', mu=omegay + pm.math.dot(Z, gamma) + X * beta + u * deltay, sigma=sigmay, observed=Y)
    # Y = pm.Bernoulli('Y', p=pm.invlogit(omegay + pm.math.dot(Z, gamma) + X * beta + u * deltay), observed=Y)  # for binary

with shrink_model:    
    trace = pm.sample(draws=5000, tune=5000, chains=4, cores=4, target_accept=0.9,nuts_sampler="blackjax")

    # Get the samples
    subdata = trace.posterior['beta']
    print(az.summary(subdata))

import pickle
with open('RBC_DBP_beta.pkl', 'wb') as f:
    pickle.dump(subdata, f)

