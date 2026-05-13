.dtam_dataset_catalogue <- function() {
  data.frame(
    name = c(
      "ACMTermPremium",
      "DAT_GSW",
      "YC_LW",
      "YC_LW_FULL",
      "Data_Macro_EA_quarterly",
      "Data_Macro_US_monthly",
      "Data_Macro_US_quarterly",
      "JSTdataset",
      "Shiller",
      "SPF",
      "YC_Euro",
      "YC_US",
      "YC_US_weekly"
    ),
    data_name = c(
      "ACMTermPremium",
      "DAT_GSW",
      "DAT_LW",
      "DAT_LW_FULL",
      "Data_Macro_EA_quarterly",
      "Data_Macro_US_monthly",
      "Data_Macro_US_quarterly",
      "JSTdataset",
      "Shiller",
      "SPF",
      "YC_Euro",
      "YC_US",
      "YC_US_weekly"
    ),
    frequency = c(
      "monthly",
      "monthly",
      "monthly",
      "monthly",
      "quarterly",
      "monthly",
      "quarterly",
      "annual",
      "monthly",
      "quarterly",
      "daily",
      "daily",
      "weekly"
    ),
    source = c(
      "Federal Reserve Bank of New York",
      "Federal Reserve Board",
      "Liu and Wu yield data",
      "Liu and Wu yield data",
      "Area-Wide Model database",
      "FRED",
      "FRED",
      "Jorda-Schularick-Taylor Macrohistory Database",
      "Robert Shiller data",
      "Philadelphia Fed Survey of Professional Forecasters",
      "ECB",
      "FRED",
      "FRED"
    ),
    description = c(
      "Adrian-Crump-Moench Treasury yields, term premia, and risk-neutral yields.",
      "Gurkaynak-Sack-Wright nominal and real zero-coupon Treasury yields.",
      "Selected Liu-Wu U.S. zero-coupon Treasury yields.",
      "Full Liu-Wu U.S. zero-coupon Treasury yield panel.",
      "Euro-area quarterly real GDP, HICP, unemployment, inflation, and GDP growth.",
      "Monthly U.S. macroeconomic and yield data.",
      "Quarterly U.S. macroeconomic and yield data.",
      "Long-run macro-financial panel from the Macrohistory Database.",
      "Shiller monthly stock-market, valuation, price, and consumption data.",
      "Survey forecasts for CPI, GDP, and Treasury-bill rates.",
      "Euro-area zero-coupon yields.",
      "Daily U.S. Treasury yields from FRED.",
      "Weekly U.S. rates, policy targets, and Treasury yields from FRED."
    ),
    source_url = c(
      "https://www.newyorkfed.org/research/data_indicators/term-premia-tabs",
      "https://www.federalreserve.gov/data/nominal-yield-curve.htm",
      "https://sites.google.com/view/jingcynthiawu/yield-data",
      "https://sites.google.com/view/jingcynthiawu/yield-data",
      "https://eabcn.org/data/area-wide-model",
      "https://fred.stlouisfed.org/",
      "https://fred.stlouisfed.org/",
      "https://www.macrohistory.net/database/",
      "https://www.econ.yale.edu/~shiller/data.htm",
      "https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters",
      "https://data.ecb.europa.eu/",
      "https://fred.stlouisfed.org/",
      "https://fred.stlouisfed.org/"
    ),
    reference = c(
      "Adrian, Crump, and Moench (2013), Pricing the term structure with linear regressions, Journal of Financial Economics 110(1), 110-138.",
      "Gurkaynak, Sack, and Wright (2007), The U.S. Treasury yield curve: 1961 to the present, Journal of Monetary Economics 54(8), 2291-2304; Gurkaynak, Sack, and Wright (2010), The TIPS yield curve and inflation compensation, American Economic Journal: Macroeconomics 2(1), 70-92.",
      "Liu and Wu (2021), Reconstructing the yield curve, Journal of Financial Economics 142(3), 1395-1425.",
      "Liu and Wu (2021), Reconstructing the yield curve, Journal of Financial Economics 142(3), 1395-1425.",
      "Fagan, Henry, and Mestre (2001), An area-wide model (AWM) for the euro area, ECB Working Paper No. 42.",
      "Federal Reserve Bank of St. Louis, FRED, and original source agencies listed in each FRED series.",
      "Federal Reserve Bank of St. Louis, FRED, and original source agencies listed in each FRED series.",
      "Jorda, Schularick, and Taylor (2017), Macrofinancial History and the New Business Cycle Facts, NBER Macroeconomics Annual 31; see also Jorda et al. (2019) for return data and Jorda et al. (2021) for bank balance-sheet ratios.",
      "Shiller (2015), Irrational Exuberance, 3rd ed., Princeton University Press; see also Shiller's online data page.",
      "Federal Reserve Bank of Philadelphia, Survey of Professional Forecasters.",
      "European Central Bank, ECB Data Portal.",
      "Federal Reserve Bank of St. Louis, FRED, and original source agencies listed in each FRED series.",
      "Federal Reserve Bank of St. Louis, FRED, and original source agencies listed in each FRED series."
    ),
    license = c(
      "Federal Reserve Bank of New York public research data; reuse is subject to the New York Fed website terms and attribution requirements.",
      "Federal Reserve Board staff research data; subject to Federal Reserve Board website terms and model disclaimer.",
      "Authors' public research data; cite Liu and Wu (2021) and check the authors' website for current reuse terms.",
      "Authors' public research data; cite Liu and Wu (2021) and check the authors' website for current reuse terms.",
      "AWM database made available for research use by ECB/EABCN; cite Fagan, Henry, and Mestre (2001) and check the AWM page for current terms.",
      "FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions.",
      "FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions.",
      "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0), with required citation and no commercial resale/integration by data providers.",
      "Robert Shiller online data; cite the source page/book and check Yale/Shiller page terms before redistribution.",
      "Federal Reserve Bank of Philadelphia public survey data; reuse subject to Philadelphia Fed website terms and attribution.",
      "Publicly released ESCB/ECB statistics may be reused free of charge if the source is quoted and the statistics/metadata are not modified; third-party data are excluded.",
      "FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions.",
      "FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions."
    ),
    stringsAsFactors = FALSE
  )
}

.dtam_validate_dataset_name <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a character scalar.", call. = FALSE)
  }
  catalogue <- .dtam_dataset_catalogue()
  if (!name %in% catalogue$name) {
    stop(
      "`name` must be one of: ",
      paste(catalogue$name, collapse = ", "),
      call. = FALSE
    )
  }
  catalogue[catalogue$name == name, , drop = FALSE]
}

.dtam_date_range <- function(x) {
  date_candidates <- intersect(c("date", "Date", "DATE", "q", "year", "Year"), names(x))
  if (length(date_candidates) == 0L) {
    return(c(start = NA_character_, end = NA_character_))
  }

  date_col <- x[[date_candidates[1L]]]
  if (inherits(date_col, "Date") || inherits(date_col, "POSIXt")) {
    return(format(range(date_col, na.rm = TRUE)))
  }
  if (is.numeric(date_col) || is.character(date_col)) {
    rng <- range(date_col, na.rm = TRUE)
    return(as.character(rng))
  }
  c(start = NA_character_, end = NA_character_)
}

.dtam_dataset_summary <- function(name) {
  x <- dtam_dataset(name)
  dims <- dim(x)
  if (is.null(dims)) {
    n_obs <- length(x)
    n_vars <- NA_integer_
  } else {
    n_obs <- dims[1L]
    n_vars <- dims[2L]
  }
  range <- if (is.data.frame(x)) .dtam_date_range(x) else c(NA_character_, NA_character_)
  data.frame(
    n_obs = n_obs,
    n_vars = n_vars,
    start = range[1L],
    end = range[2L],
    stringsAsFactors = FALSE
  )
}

#' List DTAM datasets
#'
#' Returns a catalogue of the datasets shipped with DTAM.
#'
#' @param details Logical. If \code{TRUE}, add dimensions and date/sample
#'   coverage by loading each dataset.
#'
#' @return A data frame with object names, underlying data names, broad
#'   frequency, source, description, source URL, reference, and upstream
#'   licence/terms note. With \code{details = TRUE}, the data frame also
#'   includes the number of observations, number of variables, and start/end
#'   sample markers when a date-like column is available.
#' @export
#'
#' @examples
#' dtam_datasets()
#' dtam_datasets(details = TRUE)
dtam_datasets <- function(details = FALSE) {
  catalogue <- .dtam_dataset_catalogue()
  if (!isTRUE(details)) {
    return(catalogue)
  }
  summaries <- do.call(
    rbind,
    lapply(catalogue$name, .dtam_dataset_summary)
  )
  cbind(catalogue, summaries, row.names = NULL)
}

#' Return a DTAM dataset
#'
#' Loads one of the datasets shipped with DTAM and returns it as an object.
#'
#' @param name Character scalar. Name of the dataset to return. See
#'   \code{\link[=dtam_datasets]{dtam_datasets()}}.
#'
#' @return The requested dataset.
#' @export
#'
#' @examples
#' dtam_datasets()
#' YC_US <- dtam_dataset("YC_US")
dtam_dataset <- function(name) {
  dataset <- .dtam_validate_dataset_name(name)
  env <- new.env(parent = emptyenv())
  utils::data(list = dataset$data_name, package = "DTAM", envir = env)
  get(name, envir = env, inherits = FALSE)
}

#' Describe a DTAM dataset
#'
#' Returns metadata, dimensions, sample coverage, and variable-level
#' information for one dataset shipped with DTAM.
#'
#' @param name Character scalar. Name of the dataset to describe. See
#'   \code{\link[=dtam_datasets]{dtam_datasets()}}.
#'
#' @return A list with two elements: \code{metadata}, a one-row data frame with
#'   dataset-level information, and \code{variables}, a data frame with variable
#'   names, classes, and missing-value counts.
#' @export
#'
#' @examples
#' info <- dtam_dataset_info("Data_Macro_EA_quarterly")
#' info$metadata
#' info$variables
dtam_dataset_info <- function(name) {
  .dtam_validate_dataset_name(name)
  metadata <- dtam_datasets(details = TRUE)
  metadata <- metadata[metadata$name == name, , drop = FALSE]

  x <- dtam_dataset(name)
  if (is.data.frame(x)) {
    variables <- data.frame(
      variable = names(x),
      class = vapply(x, function(z) paste(class(z), collapse = "/"), character(1L)),
      n_missing = vapply(x, function(z) sum(is.na(z)), integer(1L)),
      stringsAsFactors = FALSE
    )
  } else {
    variables <- data.frame(
      variable = name,
      class = paste(class(x), collapse = "/"),
      n_missing = sum(is.na(x)),
      stringsAsFactors = FALSE
    )
  }

  list(metadata = metadata, variables = variables)
}

#' Load a DTAM dataset into an environment
#'
#' Convenience wrapper around \code{\link[utils:data]{data()}} that validates
#' the dataset name against \code{\link[=dtam_datasets]{dtam_datasets()}}.
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
