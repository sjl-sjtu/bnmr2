# Bayesian Network-based Mendelian Randomization (BNMR)
This is an R package to conduct causal estimation between exposure and outcome on GWAS data using BNMR model, which is a two-staged framework to deal with the imperfect genetic instruments. It selects approaciate genetic instruments via random graph forest, which comprises a series Bayesian network structure learning among sampled variants and exposure. It then uses the Bayesian Mendelian randomization with a shrinkage prior to cope with horizontal pleiotropy and obtain a robust estimate. 

![image](https://github.com/sjl-sjtu/bnmr/blob/main/FIG/Fig1.jpg)

You can install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/bnmr")
```

Usage and examples can be found at https://github.com/sjl-sjtu/bnmr/blob/main/bnmr_0.2.1.pdf.

For large-scale dataset (like biobank), we recommend to conduct Bayesian MR analysis using Python Package PyMC with NUTS JAX samplers (NumPyro or BlackJAX) and GPU (https://www.pymc-labs.com/blog-posts/pymc-stan-benchmark/) to achieve faster posterior sampling. A tutorial of PyMC (v5) with JAX and Numba can be found at https://www.pymc.io/projects/examples/en/latest/samplers/fast_sampling_with_jax_and_numba.html. An example in BNMR can be found at https://github.com/sjl-sjtu/bnmr/blob/main/BayesianMR_example_pymc.py.

Supplementary notes, tables, and figures for the paper can be found at https://github.com/sjl-sjtu/bnmr/blob/main/supplementary_notes.pdf.

Latest updation at Nov. 11th, 2023.

Contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)
