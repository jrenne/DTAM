# Rebuild the three-series euro-area macroeconomic dataset used by DTAM.
# Run this script from the root of the DTAM source repository.

AWM_raw <- read.csv(
  "data-raw/AWMD_Mar2026_extract.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c("q", "YER", "HICP", "URX")
if (!identical(names(AWM_raw), required_columns)) {
  stop(
    "AWMD extract must contain exactly: ",
    paste(required_columns, collapse = ", "),
    call. = FALSE
  )
}

quarter_to_date <- function(q) {
  if (!all(grepl("^[0-9]{4}Q[1-4]$", q))) {
    stop("Invalid AWMD quarter identifier.", call. = FALSE)
  }
  year <- substr(q, 1, 4)
  month <- 1 + 3 * (as.numeric(substr(q, 6, 6)) - 1)
  as.Date(sprintf("%s-%02d-01", year, month))
}

Data_Macro_EA_quarterly <- data.frame(
  date = quarter_to_date(AWM_raw$q),
  q = AWM_raw$q,
  YER = as.numeric(AWM_raw$YER),
  HICP = as.numeric(AWM_raw$HICP),
  URX = as.numeric(AWM_raw$URX),
  stringsAsFactors = FALSE
)

Data_Macro_EA_quarterly$pi <- c(
  NaN,
  diff(log(Data_Macro_EA_quarterly$HICP))
)
Data_Macro_EA_quarterly$dy <- c(
  NaN,
  diff(log(Data_Macro_EA_quarterly$YER))
)

if (any(!complete.cases(Data_Macro_EA_quarterly[, c("YER", "HICP", "URX")]))) {
  stop("The public AWMD extract contains missing observations.", call. = FALSE)
}

save(Data_Macro_EA_quarterly, file = "data/Data_Macro_EA_quarterly.rda")
