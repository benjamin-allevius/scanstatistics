---
title: 'scanstatistics: space-time anomaly detection using scan statistics'
authors:
- affiliation: 1
  name: Benjamin Allévius
  orcid: 0000-0002-0927-7183
date: "2 May 2018"
bibliography: paper.bib
tags:
- scan statistic
- cluster detection
- anomaly detection
- spatiotemporal
affiliations:
- index: 1
  name: Department of Mathematics, Stockholm University
---

# Summary

The R package `scanstatistics` enables the detection of anomalous space-time 
clusters using the scan statistics methodology. Scan statistics are commonly 
applied in disease surveillance, where they are used to detect disease outbreaks
as they emerge locally. In this setting, cases of a given disease are recorded 
continuously across a country, and are then aggregated spatially to (say) 
district level, and temporally to (say) weekly counts. Scan statistics 
accomplish the detection task by searching the recent records of clusters of 
neighboring districts for patterns that seem anomalous given either past counts 
or the counts outside the cluster currently searched.

The `scanstatistics` package implements several scan statistics, making it a 
partially overlapping complement to existing scan statistic software such as 
[SaTScan](https://www.satscan.org/). For example, the conditional Poisson 
[@Kulldorff2001] and space-time permutation [@Kulldorff2005] scan statistics 
are available in both [SaTScan](https://www.satscan.org/) and `scanstatistics`, 
while only the latter implements scan statistics for zero-inflated data 
[@Allevius2018], count data with overdispersion [@Tango2011], an unconditional
(expectation-based) Poisson scan statistic [@Neill2005], and a Bayesian scan 
statistic [@Neill2006].

The R package `scanstatistics` is available on 
[CRAN](https://cran.r-project.org/package=scanstatistics) and its source code
is available on [GitHub](https://github.com/BenjaK/scanstatistics).


# References
