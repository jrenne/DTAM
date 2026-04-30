\name{psi.RSVAR}
\alias{psi.RSVAR}
\alias{simul.RSVAR}
\title{Regime-switching Gaussian VAR transform and simulation}
\usage{
psi.RSVAR(u, psi.parameterization)

simul.RSVAR(model, nb.sim, y0 = NaN, z0 = NaN)
}
\arguments{
\item{u}{Matrix of Laplace-transform loadings. It must have \code{n + J} rows,
where \code{n} is the dimension of \code{y_t} and \code{J} is the number of regimes.}

\item{psi.parameterization, model}{List containing \code{Phi}, \code{Pi}, and
either \code{Sigma} or \code{Sigma12}. The optional entry \code{mu} is an
\code{n x J} matrix. If omitted, it is set to zero. For \code{simul.RSVAR()}
only, \code{model} may also contain an optional current-regime intercept
\code{M}.}

\item{nb.sim}{Number of simulated dates.}

\item{y0}{Optional initial continuous state.}

\item{z0}{Optional initial regime index.}
}
\value{
\code{psi.RSVAR()} returns a list with affine-transform coefficients \code{a}
and \code{b}. \code{simul.RSVAR()} returns a list with \code{y}, \code{z},
\code{w}, and \code{regime}.
}
\description{
Computes the one-period conditional Laplace transform associated with a
regime-switching Gaussian VAR. The state vector is ordered as
\code{w_t = (y_t^T, z_t^T)^T}, where \code{z_t} is a finite-state Markov chain
represented by a selection vector.
}
\details{
The continuous component follows
\deqn{
y_{t+1} = \mu z_{t+1} + \Phi y_t
  + \Sigma(z_{t+1})^{1/2} eps_{t+1}.
}
The function \code{psi.RSVAR()} implements the conditional log-Laplace transform for the canonical specification used in the book: both the drift and the volatility of \eqn{y_{t+1}} are selected by the next regime \eqn{z_{t+1}}. If an entry \code{M} is supplied in \code{psi.parameterization}, it is deliberately ignored by \code{psi.RSVAR()}.

The simulation helper \code{simul.RSVAR()} simulates the same canonical model when \code{M} is omitted. It also accepts an optional current-regime intercept \code{M}; in that case it simulates the more general dynamics
\deqn{
y_{t+1} = \mu z_{t+1} + M z_t + \Phi y_t
  + \Sigma(z_{t+1})^{1/2} eps_{t+1}.
}
This extension is kept for backward compatibility with other numerical examples.
}
\references{
Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
\emph{Asset Pricing with Discrete-Time Affine Processes}.
}
\examples{
model <- list(
  Phi = matrix(c(0.95, 0.05,
                 0.00, 0.90), 2, 2, byrow = TRUE),
  Pi = matrix(c(0.95, 0.05,
                0.30, 0.70), 2, 2, byrow = TRUE),
  mu = matrix(c(0.001, 0.000,
                0.000, 0.002), 2, 2, byrow = TRUE),
  Sigma = array(c(0.01^2, 0.00,
                  0.00,   0.02^2,
                  0.02^2, 0.00,
                  0.00,   0.04^2), c(2, 2, 2))
)
u <- matrix(c(1, 0, 0, 0), 4, 1)
psi.RSVAR(u, model)
sim <- simul.RSVAR(model, nb.sim = 20, y0 = c(0, 0), z0 = 1)
dim(sim$w)
}
