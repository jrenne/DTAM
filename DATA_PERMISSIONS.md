# DTAM data-permissions tracker

This is a working record for the first public release of DTAM. It distinguishes
published terms that appear sufficient for research redistribution from cases
where written clarification should be obtained. It is not legal advice.

Last reviewed: 2026-08-20.

| DTAM data | Published position | Release action |
| --- | --- | --- |
| `ACMTermPremium` | The New York Fed publishes terms permitting use, copying, modification, and distribution subject to attribution and other stated conditions. | No individual permission request currently indicated; retain the required attribution and terms link. |
| `DAT_GSW` | Federal Reserve Board staff research data are publicly downloadable and accompanied by a model disclaimer. | No individual permission request currently indicated; retain the source, papers, and disclaimer. |
| `YC_Euro` | The ECB permits reuse of publicly released ECB/ESCB statistics with source attribution, subject to exclusions for third-party material. | No individual permission request currently indicated; retain ECB attribution and verify that the series contain no third-party component. |
| `SPF` | The Philadelphia Fed makes the Survey of Professional Forecasters files publicly downloadable. | Retain source attribution; confirm the applicable Philadelphia Fed website terms during the final release audit. |
| `JSTdataset` | The Macrohistory Database states CC BY-NC-SA 4.0 terms. Because DTAM accompanies a commercially published book, the application of the NonCommercial restriction should be clarified even though the package itself is freely distributed. | A clarification request is pending. DTAM now retains only the six variables and 16 countries used in the book, without shortening their samples; the complete upstream workbook is ignored and will not be published in the repository. Before the first stable release, obtain written confirmation or replace the bundled extract with user-side download/import instructions. If retained, keep it outside the MIT code licence and preserve attribution, the licence link, modification notices, the NonCommercial restriction, and ShareAlike conditions. |
| `YC_LW`, `YC_LW_FULL` | The authors' page invites downloads and specifies a citation, but does not state an explicit redistribution licence. Cynthia Wu gave written permission on 2026-08-19 to include the requested book extract in correspondence copied to Yan Liu. This correspondence is retained as the permission record for the requested extract. | Retain only the maturities used by the book: the selected curve points in `YC_LW` and the dense one-month-to-ten-year panel in `YC_LW_FULL`. Keep the complete upstream file ignored, preserve the Liu and Wu (2021) citation, and state that the data are outside the MIT software licence. |
| `Shiller` | The Yale page offers downloads but states no explicit redistribution licence; the series also combine material originating from several sources. | Contact the data author or responsible Yale contact before the first release. If permission is not obtained, replace the bundled copy with download/construction instructions. |
| `Data_Macro_EA_quarterly` | On 2026-08-19, the maintainers gave written permission to redistribute the data and stated that EABCN had approved it. They requested citation of İpek and Kısacıkoğlu (2026) for data construction and the EABCN page for the current vintage. | Resolved. Ship only the `YER`, `HICP`, and `URX` extract used by DTAM, keep the full sample, preserve the requested citations, state that the data are outside DTAM's MIT software licence, and retain the authorization correspondence privately in the ignored `data-raw/Authorizations/` folder. |
| `Data_Macro_US_monthly`, `Data_Macro_US_quarterly`, `YC_US`, `YC_US_weekly` | FRED applies series-specific copyright labels. Representative series checked on 2026-08-18 (`DTB3`, `THREEFY1`, `CPIAUCSL`, `BBKMGDP`, `PCE`, `GDPC1`, and `THREEFYTP10`) were labelled “Public Domain: Citation Requested.” | Audit every remaining constituent series against its current FRED label and record the result. The current list is `DTB4WK`, `DTB3`, `DTB6`, `THREEFY1`–`THREEFY10`, `DFF`, `DFEDTAR`, `DFEDTARU`, `DFEDTARL`, `CPIAUCSL`, `BBKMGDP`, `PCE`, `PCEPI`, `THREEFYTP1`, `THREEFYTP2`, `THREEFYTP3`, `THREEFYTP5`, `THREEFYTP7`, `THREEFYTP10`, `GDPPOT`, and `GDPC1`. Contact an original owner only for a series marked as requiring prior approval; otherwise retain the required citation. |

## Decision rule

A dataset should remain bundled in the first stable release only if its
redistribution basis is documented. When the published terms are unclear, use
one of two outcomes:

1. retain the dataset after receiving written permission; or
2. remove the bundled observations and provide a reproducible download/import
   procedure that requires the user to obtain the data from the provider.

Restricted futures prices used in the book are not DTAM data and are handled
separately in the Bookdown repository.
