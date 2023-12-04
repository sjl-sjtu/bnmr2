
First let me set the environment
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

Import neccessary libraries now:
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

Now load the dataset we want to analysis. Here we need to load two datasets: one containing the
```python
df = pd.read_csv("lym_bind.csv") # data frame of SNPs, exposures and outcomes
df_score = pd.read_csv("lym_score.csv") # adjacency score of SNP-RBC from random graph fprest
df_score = df_score.sort_values(by='score', ascending=False)
IV = df_score['snp'].head(20).tolist()
feature = ['RBC','DBP']
```
