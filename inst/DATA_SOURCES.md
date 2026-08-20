# DTAM Data Sources and Upstream Terms

DTAM ships empirical datasets used in examples for discrete-time affine
asset-pricing models. The package code is MIT-licensed, but the datasets come
from external providers and remain subject to their own citation requirements
and reuse terms. The notes below summarize the upstream sources known to DTAM;
users should check the current terms of the original provider before
redistributing or using the data commercially.

Inclusion in the package does not place a dataset under the MIT licence and
does not grant permissions that belong to the original provider. Each dataset
is a separately licensed component of the package. In particular, the
`JSTdataset` component is distributed under CC BY-NC-SA 4.0 and remains subject
to that licence's NonCommercial and ShareAlike conditions.

| DTAM object | Source | Reference | Upstream terms |
| --- | --- | --- | --- |
| `ACMTermPremium` | Federal Reserve Bank of New York term-premium data: <https://www.newyorkfed.org/research/data_indicators/term-premia-tabs> | Adrian, Crump, and Moench (2013), "Pricing the term structure with linear regressions," Journal of Financial Economics 110(1), 110-138. | Public research data; reuse is subject to New York Fed website terms and attribution requirements. |
| `DAT_GSW` | Federal Reserve Board nominal and TIPS yield-curve estimates: <https://www.federalreserve.gov/data/nominal-yield-curve.htm> | Gurkaynak, Sack, and Wright (2007), "The U.S. Treasury yield curve: 1961 to the present," Journal of Monetary Economics 54(8), 2291-2304; Gurkaynak, Sack, and Wright (2010), "The TIPS yield curve and inflation compensation," American Economic Journal: Macroeconomics 2(1), 70-92. | Federal Reserve Board staff research data; subject to Federal Reserve Board website terms and model disclaimer. |
| `YC_LW`, `YC_LW_FULL` | Liu-Wu reconstructed Treasury yield curves: <https://sites.google.com/view/jingcynthiawu/yield-data> | Liu and Wu (2021), "Reconstructing the yield curve," Journal of Financial Economics 142(3), 1395-1425. | Cynthia Wu gave written permission on 2026-08-19 to include the requested book extract in correspondence copied to Yan Liu. DTAM retains only the selected curve points and the dense one-month-to-ten-year panel used by the book. Cite Liu and Wu (2021). The data are not covered by DTAM's MIT software licence. |
| `Data_Macro_EA_quarterly` | Updated Area-Wide Model Database: <https://eabcn.org/data/area-wide-model> | İpek and Kısacıkoğlu (2026), "Estimating Euro Area Output Gap Dynamics: Evidence from the Updated Area-Wide Model Database," European Economic Review 181, 105179; Fagan, Henry, and Mestre (2001), "An area-wide model (AWM) for the euro area," ECB Working Paper No. 42. | The database maintainers gave written permission on 2026-08-19 to redistribute the three-series extract and reported EABCN approval. Cite İpek and Kısacıkoğlu (2026) and the EABCN page. The data are not covered by DTAM's MIT software licence. |
| `Data_Macro_US_monthly`, `Data_Macro_US_quarterly`, `YC_US`, `YC_US_weekly` | FRED, with original source agencies as listed in each FRED series: <https://fred.stlouisfed.org/> | Cite FRED and each original source agency where appropriate. | FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions. |
| `JSTdataset` | Jorda-Schularick-Taylor Macrohistory Database: <https://www.macrohistory.net/database/> | Jorda, Schularick, and Taylor (2017), "Macrofinancial History and the New Business Cycle Facts," NBER Macroeconomics Annual 31. Users of return data should also cite Jorda et al. (2019), "The Rate of Return on Everything, 1870-2015," Quarterly Journal of Economics 134(3), 1225-1298. | DTAM contains only the six variables and 16 countries used by the book, with the full available sample. The upstream database states CC BY-NC-SA 4.0 terms; written clarification about bundling this extract with the free companion package is pending. The data are outside DTAM's MIT software licence. |
| `Shiller` | Robert Shiller online stock-market data, augmented in DTAM with FRED series: <https://www.econ.yale.edu/~shiller/data.htm> | Shiller (2015), Irrational Exuberance, 3rd edition, Princeton University Press. | Cite the source page/book and check the Yale/Shiller page terms before redistribution. |
| `SPF` | Federal Reserve Bank of Philadelphia Survey of Professional Forecasters: <https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters> | Federal Reserve Bank of Philadelphia, Survey of Professional Forecasters documentation. | Public survey data; reuse is subject to Philadelphia Fed website terms and attribution. |
| `YC_Euro` | ECB Data Portal euro-area zero-coupon yields: <https://data.ecb.europa.eu/> | European Central Bank, ECB Data Portal. | Publicly released ESCB/ECB statistics may be reused free of charge if the source is quoted and the statistics/metadata are not modified; third-party data are excluded. |

## Restricted book-only inputs

Some figures in the companion book use licensed futures prices obtained through
the authors' institutional data access. Those observations are not part of DTAM.
The Bookdown page identifies the relevant provider at the point of use.
