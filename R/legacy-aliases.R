#' Legacy dot-style aliases
#'
#' These wrappers preserve backward compatibility with earlier DTAM code that
#' used dot-style function names. New code should use the underscore-named
#' replacements documented throughout the package.
#'
#' @param ... Arguments passed to the underscore-named replacement.
#' @name legacy_aliases
NULL

#' @rdname legacy_aliases
#' @export
reverse.MHLT <- function(...) reverse_MHLT(...)

#' @rdname legacy_aliases
#' @export
psi.GaussianVAR <- function(...) psi_GaussianVAR(...)

#' @rdname legacy_aliases
#' @export
psi.GaussianQVAR <- function(...) psi_GaussianQVAR(...)

#' @rdname legacy_aliases
#' @export
psi.VARG <- function(...) psi_VARG(...)

#' @rdname legacy_aliases
#' @export
psi.VARG_Poisson <- function(...) psi_VARG_Poisson(...)

#' @rdname legacy_aliases
#' @export
psi.QPoisson <- function(...) psi_QPoisson(...)

#' @rdname legacy_aliases
#' @export
psi.RSVAR <- function(...) psi_RSVAR(...)

#' @rdname legacy_aliases
#' @export
psi.BansalShaliastovich <- function(...) psi_BansalShaliastovich(...)

#' @rdname legacy_aliases
#' @export
psi.AR_CH <- function(...) psi_AR_CH(...)

#' @rdname legacy_aliases
#' @export
psi.w.TopDown <- function(...) psi_w_TopDown(...)

#' @rdname legacy_aliases
#' @export
psiQ.w.TopDown <- function(...) psiQ_w_TopDown(...)

#' @rdname legacy_aliases
#' @export
psi.y.TopDown <- function(...) psi_y_TopDown(...)

#' @rdname legacy_aliases
#' @export
simul.ARG <- function(...) simul_ARG(...)

#' @rdname legacy_aliases
#' @export
simul.GVAR <- function(...) simul_GaussianVAR(...)

#' @rdname legacy_aliases
#' @export
simul.RSVAR <- function(...) simul_RSVAR(...)

#' @rdname legacy_aliases
#' @export
simul.TopDown <- function(...) simul_TopDown(...)

#' @rdname legacy_aliases
#' @export
simul.compound.poisson <- function(...) simul_compound_poisson(...)

#' @rdname legacy_aliases
#' @export
simul.rating.migration <- function(...) simul_rating_migration(...)

#' @rdname legacy_aliases
#' @export
make.pdf <- function(...) make_pdf(...)

