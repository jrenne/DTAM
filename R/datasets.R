#' List DTAM datasets
#'
#' Returns a compact catalogue of the datasets shipped with DTAM.
#'
#' @return A data frame with object names, underlying data names, broad
#'   frequency, source, and a short description.
#' @export
dtam_datasets <- function() {
  data.frame(
    name = c(
      "ACMTermPremium",
      "Credit_spds",
      "DAT_GSW",
      "YC_LW",
      "YC_LW_FULL",
      "Data_Macro_EA_quarterly",
      "Data_Macro_US_monthly",
      "Data_Macro_US_quarterly",
      "Futures",
      "JSTdataset",
      "Shiller",
      "SPF",
      "YC_Euro",
      "YC_US",
      "YC_US_weekly"
    ),
    data_name = c(
      "ACMTermPremium",
      "Credit_spds",
      "DAT_GSW",
      "DAT_LW",
      "DAT_LW_FULL",
      "Data_Macro_EA_quarterly",
      "Data_Macro_US_monthly",
      "Data_Macro_US_quarterly",
      "Futures",
      "JSTdataset",
      "Shiller",
      "SPF",
      "YC_Euro",
      "YC_US",
      "YC_US_weekly"
    ),
    frequency = c(
      "monthly",
      "annual",
      "monthly",
      "monthly",
      "monthly",
      "quarterly",
      "monthly",
      "quarterly",
      "daily",
      "annual",
      "monthly",
      "quarterly",
      "daily",
      "daily",
      "weekly"
    ),
    source = c(
      "Federal Reserve Bank of New York",
      "FRED and S&P default data",
      "Federal Reserve Board",
      "Liu and Wu yield data",
      "Liu and Wu yield data",
      "Area-Wide Model database",
      "FRED",
      "FRED",
      "Natural gas futures data",
      "Jorda-Schularick-Taylor Macrohistory Database",
      "Robert Shiller data",
      "Philadelphia Fed Survey of Professional Forecasters",
      "ECB",
      "FRED",
      "FRED"
    ),
    description = c(
      "Adrian-Crump-Moench Treasury yields, term premia, and risk-neutral yields.",
      "U.S. corporate credit spreads and default-rate data.",
      "Gurkaynak-Sack-Wright nominal and real zero-coupon Treasury yields.",
      "Selected Liu-Wu U.S. zero-coupon Treasury yields.",
      "Full Liu-Wu U.S. zero-coupon Treasury yield panel.",
      "Euro-area quarterly real GDP, HICP, unemployment, inflation, and GDP growth.",
      "Monthly U.S. macroeconomic and yield data.",
      "Quarterly U.S. macroeconomic and yield data.",
      "Natural gas futures settlement prices.",
      "Long-run macro-financial panel from the Macrohistory Database.",
      "Shiller monthly stock-market, valuation, price, and consumption data.",
      "Survey forecasts for CPI, GDP, and Treasury-bill rates.",
      "Euro-area zero-coupon yields.",
      "Daily U.S. Treasury yields from FRED.",
      "Weekly U.S. rates, policy targets, and Treasury yields from FRED."
    ),
    stringsAsFactors = FALSE
  )
}

#' Return a DTAM dataset
#'
#' Loads one of the datasets shipped with DTAM and returns it as an object.
#'
#' @param name Character scalar. Name of the dataset to return. See
#'   [dtam_datasets()].
#'
#' @return The requested dataset.
#' @export
#'
#' @examples
#' dtam_datasets()
#' YC_US <- dtam_dataset("YC_US")
dtam_dataset <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a character scalar.", call. = FALSE)
  }
  datasets <- dtam_datasets()
  available <- datasets$name
  if (!name %in% available) {
    stop(
      "`name` must be one of: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  env <- new.env(parent = emptyenv())
  data_name <- datasets$data_name[datasets$name == name]
  utils::data(list = data_name, package = "DTAM", envir = env)
  get(name, envir = env, inherits = FALSE)
}

#' Load a DTAM dataset into an environment
#'
#' Convenience wrapper around [data()] that validates the dataset name against
#' [dtam_datasets()].
#'
#' @param name Character scalar. Name of the dataset to load.
#' @param envir Environment in which the dataset should be loaded. Defaults to
#'   the caller environment.
#'
#' @return Invisibly returns the loaded dataset.
#' @export
#'
#' @examples
#' load_dtam_dataset("YC_Euro")
load_dtam_dataset <- function(name, envir = parent.frame()) {
  obj <- dtam_dataset(name)
  assign(name, obj, envir = envir)
  invisible(obj)
}
