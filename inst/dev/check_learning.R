# Development checks for the imperfect-information solvers.
# Run after installing DTAM, or source R/learning.R before this file.

if (!exists("solve_learning", mode = "function")) {
  library(DTAM)
}

feedback_model <- list(
  mu = matrix(0),
  Phi = matrix(0.7),
  Sigma12 = matrix(0.1),
  Omega12 = matrix(0.2),
  A = matrix(1),
  B = matrix(0),
  Lambda0 = matrix(0.2),
  Lambda1 = matrix(0.1)
)

feedback_solution <- solve_learning(feedback_model)
stopifnot(
  max(feedback_solution$structural_residuals) < 1e-10,
  max(Mod(eigen(feedback_solution$Phi_ww, only.values = TRUE)$values)) < 1
)

re_model <- c(
  feedback_model,
  list(
    Psi0 = matrix(0.05),
    Psi1 = matrix(0.02),
    Gamma0 = matrix(0.02),
    Gamma1 = matrix(0.01)
  )
)
re_model$mu <- matrix(0.01)
re_model$Lambda0 <- matrix(0.02)
re_model$Lambda1 <- matrix(0.01)

re_solution <- solve_Learning_RE(re_model)
stopifnot(
  all(re_solution$converged),
  max(re_solution$structural_residuals) < 1e-8,
  max(Mod(eigen(re_solution$Phi_ww, only.values = TRUE)$values)) < 1
)

message("All imperfect-information solver checks passed.")
