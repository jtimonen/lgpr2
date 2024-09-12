# lgpr2

## Differences compared to the original lgpr package


- `lgpr2` uses `CmdstanR` as backend, whereas `lgpr` uses `rstan`
- `lgpr2` generates and compiles new Stan code each time a model is specified, whereas `lgpr` uses precompiled Stan code
- `lgpr2` uses approximate GPs (basis function implementation) whereas `lgpr` uses exact ones (covariance matrix implementation)
- `lgpr2` is not on CRAN
- `lgpr2` is a lot faster and scales to bigger data
- `lgpr2` has less features (for example, currently only supports Gaussian observation model)
-  Stan code generated  by `lgpr2` is more readable and can be easily customized
-  `lgpr2` contains the projection predictive method additional model reduction technique


## Requirements
* [CmdStanR](https://mc-stan.org/cmdstanr/)

## Installation

```r
devtools::install_github('jtimonen/lgpr2')
``` 

## Reference

Juho Timonen & Harri Lähdesmäki, **Scalable mixed-domain Gaussian process modeling and model reduction for longitudinal data**, [arXiv](https://arxiv.org/abs/2111.02019) (2024)
