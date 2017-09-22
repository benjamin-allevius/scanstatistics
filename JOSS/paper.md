---
title: 'scanstatistics: space-time anomaly detection using scan statistics'
authors:
- affiliation: 1
  name: Benjamin Allévius
  orcid: 0000-0002-0927-7183
date: "22 September 2017"
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

The `scanstatistics` package implements several scan statistics, some not 
publically available elsewhere. The key references for these scan statistics 
are:

* @Kulldorff2001
* @Kulldorff2005
* @Neill2005
* @Neill2006
* @Tango2011
* @Allevius2017

The R package `scanstatistics` is available on 
[CRAN](https://cran.r-project.org/package=scanstatistics) and its source code
is available on [GitHub](https://github.com/BenjaK/scanstatistics).


# References
