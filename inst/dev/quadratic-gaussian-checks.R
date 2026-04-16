# Checks and exploratory snippets moved out of package code during cleanup.

mu <- matrix(c(0, 1, .2), 3, 1)
Phi <- .8 * diag(3)
Phi[2, 1] <- .2
Phi[1, 3] <- .1
Phi[3, 1] <- -.1
aux <- matrix(.1 * rnorm(9), 3, 3)
Sigma <- aux %*% t(aux)

# res1 <- vcov_matrix_Quadratic(mu, Phi, Sigma,
#                               BigMats = make_BigMats(dim(Phi)[1]))
# res2 <- make_matrices_cond_mean_variance_Quadratic(
#   list(mu = mu, Phi = Phi, Sigma = Sigma),
#   indic_compute_V = TRUE
# )

model <- list(mu = mu, Phi = Phi, Sigma = Sigma)
n <- dim(Phi)[1]

Vred <- .1 * rnorm(n * (n + 1) / 2)
u <- matrix(c(.1 * rnorm(n), Vred), ncol = 1)

psi_Quad <- psi.GaussianQVAR(u, psi.parameterization = model)

model$n_w <- 9
res <- compute_expect_variance(psi.GaussianQVAR, model, du = 10^(-6))
res2 <- make_matrices_cond_mean_variance_Quadratic(model, indic_compute_V = TRUE)

x <- simul.GVAR(model, nb.sim = 10000)
xx <- (x %x% matrix(1, 1, 3)) * (matrix(1, 1, 3) %x% x)
w <- cbind(x, xx)
Vw <- var(w)

Sigma12 <- matrix(rnorm(4), 2, 2)
aux <- matrix(rnorm(4), 2, 2)
V <- aux + t(aux)
Sigma12 %*% solve(diag(2) - 2 * t(Sigma12) %*% V %*% Sigma12) %*% t(Sigma12)
solve(solve(Sigma12 %*% t(Sigma12)) - 2 * V)
