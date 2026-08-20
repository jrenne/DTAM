# DTAM

`DTAM` provides R tools for simulation, estimation, filtering, and pricing in
discrete-time affine asset-pricing models. It is the companion package for
*Asset Pricing with Discrete-Time Affine Processes*.

The executable book examples are available on the
[companion Bookdown site](https://jrenne.github.io/MPRR_codes/).

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("jrenne/DTAM")
```

Then load the package:

```r
library(DTAM)
```

## Data included with DTAM

DTAM includes empirical datasets used by the book examples. The package
provides helpers for discovering and inspecting them:

```r
dtam_datasets()
dtam_datasets(details = TRUE)
info <- dtam_dataset_info("JSTdataset")
info$metadata[, c("source_url", "reference", "redistribution_status", "license")]
info$variables
YC_US <- dtam_dataset("YC_US")
```

The MIT licence applies to the DTAM software code. It does **not** replace the
licences, attribution requirements, or reuse restrictions of third-party
datasets. See [`inst/DATA_SOURCES.md`](inst/DATA_SOURCES.md) before
redistributing a dataset or using it commercially. In particular,
`JSTdataset` is subject to CC BY-NC-SA 4.0 rather than the MIT licence.

Restricted futures prices used for selected book figures are not distributed
through DTAM. The companion site identifies the source beside the relevant
figures.

## Rebuilding data from FRED

The data-construction script expects each user to supply their own FRED API
key through the `FRED_API_KEY` environment variable. Never place a key in a
tracked source file.

```r
Sys.setenv(FRED_API_KEY = "your-personal-key")
```

FRED recommends a distinct key for each application and requires users of an
application to use their own key.

## Development status and support

DTAM is under active development and has not yet reached its first stable
release. Please report reproducible software problems through
[GitHub Issues](https://github.com/jrenne/DTAM/issues). Questions about the
book's examples belong in the
[`MPRR_codes` repository](https://github.com/jrenne/MPRR_codes/issues).

Before the first release, the remaining data-permission questions and the
steps in [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) will be resolved.
