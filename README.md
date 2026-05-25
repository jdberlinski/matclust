[![R-CMD-check](https://github.com/jdberlinski/matclust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jdberlinski/matclust/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/jdberlinski/matclust/graph/badge.svg)](https://codecov.io/gh/jdberlinski/matclust)

# Clustering observations with unequal replicates.

The function `repclust()` is used to cluster observations via a finite mixture of matrix-variate normal distributions.

Observations are allowed to have an unequal number of observations, both between observations and between features of the same observation.

See `?repclust` for more information on the model, and `?generate_data` to see the expected format for data.

## Installation
With `pak`:
```{r}
pak::pak("jdberlinski/matclust")
```

## Example usage
```{r}
library(matclust)

simulated_data <- generate_data(1000, 4, 10, 0.2)
res <- repclust(simulated_data$data, 4)
```
