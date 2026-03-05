# Clustering observations with unequal replicates.

The function `repclust()` is used to cluster observations via a finite mixture of matrix-variate normal distributions.

Observations are allowed to have an unequal number of observations, both between observations and between features of hte same observation.

See `?repclust` for more information on the model, and `?generate_data` to see the expect format for data.funciton.

## Installation
With `devtools` (or `remotes`):
```{r}
devtools::install_github("jdberlinski/matclust)
```

## Example usage
```{r}
library(matclust)

simulated_data <- generate_data(1000, 4, 10, 0.2)
res <- repclust(simulated_data$data, 4)
```