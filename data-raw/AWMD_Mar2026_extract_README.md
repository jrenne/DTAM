# AWMD March 2026 extract

`AWMD_Mar2026_extract.csv` is the reproducible input for the DTAM object
`Data_Macro_EA_quarterly`. It contains only the three series used by DTAM and
the companion book examples:

- `YER`: real GDP;
- `HICP`: Harmonised Index of Consumer Prices; and
- `URX`: unemployment rate.

The extract retains the complete available sample, 1970Q1--2025Q4, from the
`AWMD_Updated` (EA21) sheet of the March 2026 workbook. The alternative
Ireland-excluded sheet was not used because it does not contain an
Ireland-excluded HICP series, so mixing the two sheets would produce an
inconsistent three-variable dataset.

## Provenance and transformation

- Source workbook: `AWMD_Mar2026.xlsx`, updated 14 March 2026.
- Source URL: <https://docs.google.com/spreadsheets/d/19yUIA6dMAN9eSuzD5dfUutBS73SFQ4A3/edit>
- EABCN data page: <https://eabcn.org/data/area-wide-model>
- Source-workbook SHA-256: `10c0bd0ff9e90937dcd405dd8cd630f7fac9f59775ea544c5bfdb0f39bbf1829`.
- Extract SHA-256: `652a67cbdf15041999409455373e73dc8ac18a30b672b20a727f7da2bcecc13d`.

The extract was made by selecting `Period`, `YER`, `HICP`, and `URX` from
`AWMD_Updated`, renaming `Period` to `q`, removing the hyphen from quarter
identifiers (for example, `1970-Q1` becomes `1970Q1`), and dropping the final
empty 2026Q1 row. No observations in the three selected series were otherwise
filtered or transformed. `data-raw/make_awmd_dataset.R` constructs dates and
the quarterly log changes `pi` and `dy` used by DTAM.

## Citation and redistribution

Please cite:

İpek, M. S., and B. Kısacıkoğlu (2026), “Estimating Euro Area Output Gap
Dynamics: Evidence from the Updated Area-Wide Model Database,” *European
Economic Review* 181, 105179. <https://doi.org/10.1016/j.euroecorev.2025.105179>

The database maintainers gave the DTAM/book authors written permission on 19
August 2026 to redistribute the data and reported that EABCN had approved the
redistribution. The correspondence is retained privately under
`data-raw/Authorizations/` and is intentionally excluded from the public Git
repository and package source archive. This data extract remains subject to its
own attribution requirements and is not placed under DTAM's MIT software
licence.
