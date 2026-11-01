# Regression: Fourier densities must have probability per unit of the variable,
# and their heights must not change when the plotting grid is refined.
library(DTAM)
psi_normal <- function(u, model) exp(model$mean*u + model$sd^2*u^2/2)
model <- list(mean=.2,sd=.8)
peaks <- numeric(2)
for (i in seq_along(c(201L,401L))) {
  n <- c(201L,401L)[i]
  grid <- seq(-4,4,length.out=n)
  density <- make_pdf(model,grid,psi_normal,max_x=100,nb_x=10000)
  edges <- seq(grid[1],tail(grid,1),length.out=n+1L)
  midpoints <- (head(edges,-1)+tail(edges,-1))/2
  reference <- dnorm(midpoints,model$mean,model$sd)
  stopifnot(length(density)==n,all(is.finite(density)),all(density>=0))
  stopifnot(max(abs(density-reference))<.005)
  mass <- sum(density*diff(edges))
  expected_mass <- diff(pnorm(range(grid),model$mean,model$sd))
  stopifnot(abs(mass-expected_mass)<.005)
  peaks[i] <- max(density)
}
stopifnot(abs(diff(peaks))<.005)
cat('Density normalization and grid-refinement checks passed.\n')
