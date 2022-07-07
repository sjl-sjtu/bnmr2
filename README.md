# bnmr
This is an R package to conduct causal estimation between exposure and outcome on GWAS data using BNMR model, which is a two-staged framework to deal with the imperfect genetic instruments. It uses the adaptive Bayesian Mendelian randomization framework for effect estimations in the presence of horizontal pleiotropy after IV selection by Bayesian network structure learning among variants and exposure. 

You can install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/bnmr")
```
