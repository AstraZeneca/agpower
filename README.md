![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)
# agpower
agpower: Recurrent event analysis planning using the Andersen-Gill model in R

Package implements functions useful in prospective planning and monitoring of a clinical trial when a recurrent event endpoint is to be assessed by Lin, Wei, Yang, and Ying (2010) model. The equations developed in Ingel and Jahn-Eimermacher (2014) and their consequences are employed.
In particular functions for sample size planning and assurance calculations are implemented allowing for thorough exploration of the design space (allocation ratio, follow-up time, events targeted, etc.).

The software is implemented in R and tested in version 4.1.0. The software depends on core R, though suggests dplyr, testthat, and tidyr.

After installation, the package may be used as follows:
```{r,  eval = FALSE}
library("agpower")

power.lwyy.test(N = 1000, RR = 0.8, thta = 1, L = 1000, tau = 0.9, alp = 0.05, ar = 0.5)
```

Further details and example usage is provide in the user manual here: [inst/docs/agpower-manual.pdf](inst/docs/agpower-manual.pdf)
