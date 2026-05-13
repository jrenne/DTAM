# DTAM Data Sources and Upstream Terms

DTAM ships empirical datasets used in examples for discrete-time affine
asset-pricing models. The package code is MIT-licensed, but the datasets come
from external providers and remain subject to their own citation requirements
and reuse terms. The notes below summarize the upstream sources known to DTAM;
users should check the current terms of the original provider before
redistributing or using the data commercially.

| DTAM object | Source | Reference | Upstream terms |
| --- | --- | --- | --- |
| `ACMTermPremium` | Federal Reserve Bank of New York term-premium data: <https://www.newyorkfed.org/research/data_indicators/term-premia-tabs> | Adrian, Crump, and Moench (2013), "Pricing the term structure with linear regressions," Journal of Financial Economics 110(1), 110-138. | Public research data; reuse is subject to New York Fed website terms and attribution requirements. |
| `DAT_GSW` | Federal Reserve Board nominal and TIPS yield-curve estimates: <https://www.federalreserve.gov/data/nominal-yield-curve.htm> | Gurkaynak, Sack, and Wright (2007), "The U.S. Treasury yield curve: 1961 to the present," Journal of Monetary Economics 54(8), 2291-2304; Gurkaynak, Sack, and Wright (2010), "The TIPS yield curve and inflation compensation," American Economic Journal: Macroeconomics 2(1), 70-92. | Federal Reserve Board staff research data; subject to Federal Reserve Board website terms and model disclaimer. |
| `YC_LW`, `YC_LW_FULL` | Liu-Wu reconstructed Treasury yield curves: <https://sites.google.com/view/jingcynthiawu/yield-data> | Liu and Wu (2021), "Reconstructing the yield curve," Journal of Financial Economics 142(3), 1395-1425. | Authors' public research data; cite Liu and Wu (2021) and check the authors' website for current reuse terms. |
| `Data_Macro_EA_quarterly` | Area-Wide Model database: <https://eabcn.org/data/area-wide-model> | Fagan, Henry, and Mestre (2001), "An area-wide model (AWM) for the euro area," ECB Working Paper No. 42. | Made available for research use by ECB/EABCN; check the AWM page for current terms. |
| `Data_Macro_US_monthly`, `Data_Macro_US_quarterly`, `YC_US`, `YC_US_weekly` | FRED, with original source agencies as listed in each FRED series: <https://fred.stlouisfed.org/> | Cite FRED and each original source agency where appropriate. | FRED content is free to access subject to FRED terms; some series are owned by third parties and may carry additional restrictions. |
| `JSTdataset` | Jorda-Schularick-Taylor Macrohistory Database: <https://www.macrohistory.net/database/> | Jorda, Schularick, and Taylor (2017), "Macrofinancial History and the New Business Cycle Facts," NBER Macroeconomics Annual 31. Users of return data should also cite Jorda et al. (2019), "The Rate of Return on Everything, 1870-2015," Quarterly Journal of Economics 134(3), 1225-1298. Users of bank-balance-sheet ratios should also cite Jorda et al. (2021), "Bank capital redux: solvency, liquidity, and crisis," Review of Economic Studies 88(1), 260-286. | Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0), with required citation and restrictions on commercial resale/integration by data providers. |
| `Shiller` | Robert Shiller online stock-market data, augmented in DTAM with FRED series: <https://www.econ.yale.edu/~shiller/data.htm> | Shiller (2015), Irrational Exuberance, 3rd edition, Princeton University Press. | Cite the source page/book and check the Yale/Shiller page terms before redistribution. |
| `SPF` | Federal Reserve Bank of Philadelphia Survey of Professional Forecasters: <https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters> | Federal Reserve Bank of Philadelphia, Survey of Professional Forecasters documentation. | Public survey data; reuse is subject to Philadelphia Fed website terms and attribution. |
| `YC_Euro` | ECB Data Portal euro-area zero-coupon yields: <https://data.ecb.europa.eu/> | European Central Bank, ECB Data Portal. | Publicly released ESCB/ECB statistics may be reused free of charge if the source is quoted and the statistics/metadata are not modified; third-party data are excluded. |

## Data Not Shipped

Earlier local versions of DTAM included a `Credit_spds` object combining
Moody's Seasoned Aaa/Baa corporate bond yields from FRED with annual corporate
default counts and rates from S&P Global Ratings annual default studies. These
data are no longer shipped with DTAM because the Moody's series notes on FRED
identify the data as proprietary and prohibit redistribution without Moody's
prior written consent. The S&P default-study data may also be subject to S&P
terms. The Rbookdown examples that use these data load them, when available,
from a local ignored `private_data/` folder.
