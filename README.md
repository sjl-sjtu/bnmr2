# Bayesian Network-based Mendelian Randomization (BNMR)
This is an R package to conduct causal estimation between exposure and outcome on GWAS data using BNMR model, which is a two-staged framework to deal with the imperfect genetic instruments. It selects approaciate genetic instruments via random graph forest, which comprises a series Bayesian network structure learning among sampled variants and exposure. It then uses the Bayesian Mendelian randomization with a shrinkage prior to cope with horizontal pleiotropy and obtain a robust estimate. 

![image](https://github.com/sjl-sjtu/bnmr/blob/main/FIG/overview.jpg)

You can install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/bnmr")
```

Usage and examples can be found at https://github.com/sjl-sjtu/bnmr/blob/main/bnmr_0.2.0.pdf.

Latest updation at Mar. 23rd, 2023.

Contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)
