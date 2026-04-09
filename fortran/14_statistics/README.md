# Chapter 14: Statistical Description of Data

Moments, hypothesis tests, contingency tables, and nonparametric correlation measures.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `moment` | 14.1 | Calculate moments of a data set (mean, variance, skewness, kurtosis) |
| `ttest`  | 14.2 | Student's t-test for difference of means |
| `avevar` | 14.2 | Calculate mean and variance of a data set |
| `tutest` | 14.2 | Student's t-test for means, unequal variances (Welch's test) |
| `tptest` | 14.2 | Student's t-test for means, paired data |
| `ftest`  | 14.2 | F-test for difference of variances |
| `chsone` | 14.3 | Chi-square test for difference between data and model |
| `chstwo` | 14.3 | Chi-square test for difference between two data sets |
| `ksone`  | 14.3 | Kolmogorov-Smirnov test of data against a model distribution |
| `kstwo`  | 14.3 | Kolmogorov-Smirnov test between two data sets |
| `probks` | 14.3 | Kolmogorov-Smirnov probability function |
| `cntab1` | 14.4 | Contingency table analysis using chi-square |
| `cntab2` | 14.4 | Contingency table analysis using entropy measure |
| `pearsn` | 14.5 | Pearson's (linear) correlation between two data sets |
| `spear`  | 14.6 | Spearman's rank correlation between two data sets |
| `crank`  | 14.6 | Replace array elements by their rank |
| `kendl1` | 14.6 | Kendall's tau correlation between two data sets |
| `kendl2` | 14.6 | Kendall's tau for contingency table analysis |

## Notes

- `ttest`, `tutest`, and `tptest` cover the three common t-test variants (equal variance, unequal variance, paired).
- `chsone`/`chstwo` test goodness-of-fit for binned data; `ksone`/`kstwo` test for unbinned data.
- `pearsn` measures linear correlation; `spear` and `kendl1`/`kendl2` measure rank-based (nonparametric) correlation.
- `cntab1` and `cntab2` analyze two-way contingency tables using chi-square and entropy, respectively.
