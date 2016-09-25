
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/BenjaK/scanstatistics.svg?branch=master)](https://travis-ci.org/BenjaK/scanstatistics) <!--
[![Build status](https://ci.appveyor.com/api/projects/status/tbd7gs7n2aaa1yvf/branch/master?svg=true)](https://ci.appveyor.com/project/BenjaK/scanstatistics/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/BenjaK/scanstatistics/badge.svg?branch=master)](https://coveralls.io/github/BenjaK/scanstatistics?branch=master)
-->

scanstatistics
==============

An R package for space-time anomaly detection using scan statistics.

Installing the package
----------------------

To install the development version of this package, type the following:

``` r
devtools::install_github("benjak/scanstatistics")
```

What are scan statistics?
-------------------------

Scan statistics are used to detect anomalous clusters in spatial or space-time data. The gist of the methodology, at least in this package, is this:

1.  Monitor one or more data streams at multiple *locations* over intervals of time.
2.  Form a set of space-time *clusters*, each consisting of (1) a collection of locations, and (2) an interval of time stretching from the present to some number of time periods in the past.
3.  For each cluster, compute a statistic based on both the observed and the expected responses. Report the clusters with the largest statistics.

Main functions
--------------

-   **`scan_poisson`**: computes a scan statistic for data following a Poisson distribution.
-   **`scan_negbin`**: computes a scan statistic for data following a negative binomial distribution.
-   **`scan_zip`**: computes a scan statistic for data following a zero-inflated Poisson distribution.
-   **`knn_zones`**: Creates a set of spatial *zones* (groups of locations) to scan for anomalies. Input is a matrix in which rows are the enumerated locations, and columns the \(k\) nearest neighbors. To create such a matrix, the following two functions are useful:
    -   **`coords_to_knn`**: use `stats::dist` to get the \(k\) nearest neighbors of each location into a format usable by `knn_zones`.
    -   **`dist_to_knn`**: use an already computed distance matrix to get the \(k\) nearest neighbors of each location into a format usable by `knn_zones`.
-   **`score_locations`**: Score each location by how likely it is to have an ongoing anomaly in it. This score is not theoretically motivated.
-   **`top_clusters`**: Get the top \(k\) space-time clusters, either overlapping or non-overlapping in the spatial dimension.

Example: Brain cancer in New Mexico
-----------------------------------

To demonstrate the scan statistics in this package, we will use a dataset of the annual number of brain cancer cases in the counties of New Mexico, for the years 1973--1991. This data was studied by Kulldorff et al. (1998), who detected a probable cluster of cancer cases in the counties Los Alamos and Santa Fe during the years 1986--1989. The data originally comes from the package *rsatscan* (Kleinman 2015), which provides an interface to the program [SatScan](http://www.satscan.org), but has been aggregated and extended for the *scanstatistics* package.

The counties of New Mexico may be plotted on a map as follows:

``` r
library(scanstatistics)
#> Loading required package: data.table
library(ggplot2)

# Load map data
data(NM_map)

# Place the county names at the centroid of the counties. The following function
# calculates the geographical centroids from the polygons in NM_map.
# See: https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
polygon_centroid <- function(x0, y0) {
  x1 <- c(x0[-1], x0[1])
  y1 <- c(y0[-1], y0[1])
  A <- sum(x0 * y1 - x1 * y0) / 2
  Cx <- sum((x0 + x1) * (x0 * y1 - x1 * y0)) / (6 * A)
  Cy <- sum((y0 + y1) * (x0 * y1 - x1 * y0)) / (6 * A)
  data.frame(long = Cx, lat = Cy)
}

# Calculate geographic centroids for each county
centroids <- NM_map[, polygon_centroid(long, lat), by = .(subregion)]

# Plot map with labels at centroids
ggplot() + 
  geom_polygon(data = NM_map,
               mapping = aes(x = long, y = lat, group = group),
               color = "grey", fill = "white") +
  geom_text(data = centroids, 
            mapping = aes(x = long, y = lat, label = subregion)) +
  ggtitle("Counties of New Mexico")
```

![](README_figures/unnamed-chunk-2-1.png)

### Creating spatial zones

The anomalies considered in the *scanstatistics* package have both a temporal and a spatial component. The spatial component, called a zone, consists of one or more locations grouped together according to their similarity across features. In this example, the locations are the counties of New Mexico and the features are the coordinates of the county seats. These are made available in the data table `NM_geo`. Similarity will be measured using the geographical distance between the seats of the counties, taking into account the curvature of the earth. A distance matrix is calculated using the `sp_Dists` function from the *sp* package, which is then passed to `dist_to_knn` (with \(k=15\) neighbors) and on to `knn_zones`:

``` r
library(sp)
library(magrittr)

data(NM_geo)

zones <- NM_geo[, c("long", "lat"), with = FALSE] %>%
  as.matrix %>%
  spDists(x = ., y = ., longlat = TRUE) %>%
  dist_to_knn(k = 15) %>%
  knn_zones
```

### A scan statistic for Poisson data

The Poisson distribution is a natural first option when dealing with (practically) unbounded count data. The *scanstatistics* package provides the function `scan_poisson`, which is a scan statistic for univariate Poisson-distributed data proposed by Neill et al. (2005).

The first argument to `scan_poisson` should be a **data table** with columns 'location', 'duration', 'count' and 'mean'. This table holds data for the period in which we want to detect anomalies. Locations should be encoded as the integers 1, 2, ..., which means that factor variables can be used for this purpose. The duration column counts time backwards, so that a duration of 1 is the most recent time interval, duration 2 is the second most recent, and so on.

We will create such a table by subsetting the `NM_popcas` table, which holds the population and the number of brain cancer cases for each year between \(1973--1991\) and each county of New Mexico. Note that the population numbers are (perhaps poorly) interpolated from the censuses conducted in 1973, 1982, and 1991.

``` r
data(NM_popcas)

tab <- NM_popcas[year >= 1986 & year < 1990, ]
tab[, duration := max(year) - year + 1]
tab[, location := county]
```

We still need to add the column 'mean', which should hold the predicted Poisson expected value parameter \(\mu_{it}\) for each location \(i\) and time interval \(t\). In this example we would like to detect a potential cluster of brain cancer in the counties of New Mexico during the years 1986--1989. Thus, we will use data from the years prior to 1986 to estimate the Poisson parameter for all counties in the years following. A simple generalized linear model (GLM) with a linear time trend and an offset for county population size will suffice to demonstrate the scan statistic. We fit such a model and create the needed column as follows:

``` r
# Fit model and add predicted mean parameters
mod_poisson <- glm(count ~ offset(log(population)) + 1 + I(year - 1985), 
                   data = NM_popcas[year < 1986, ], 
                   family = poisson(link = "log"))

# Add the mean column
tab[, mean := predict(mod_poisson, tab, type = "response")]
```

We can now calculate the Poisson scan statistic. To give us more confidence in our detection results, we will perform 99 Monte Carlo replications, by which data is generated using the parameters from the null hypothesis and a new scan statistic calculated. This is then summarized in a \(p\)-value, calculated as the proportion of times the replicated scan statistics exceeded the observed one. The output of `scan_poisson` is an object of class "scanstatistic", which comes with the print method seen below.

``` r
set.seed(1)
poisson_result <- scan_poisson(tab, zones, n_mcsim = 99)
print(poisson_result)
#> Data distribution:                Poisson
#> Type of scan statistic:           Expectation-based
#> Number of locations considered:   32
#> Maximum duration considered:      4
#> Number of spatial zones:          415
#> Number of Monte Carlo replicates: 99
#> p-value of observed statistic:    0.01
#> Most likely event duration:       4
#> ID of locations in most likely cluster: 15, 26
```

As we can see, the most likely cluster for an anomaly stretches from 1986--1989 and involves the locations numbered 15 and 26, which correspond to the counties

``` r
counties <- as.character(NM_geo$county)
counties[c(15, 26)]
[1] "losalamos" "santafe"  
```

These are the same counties detected by Kulldorff et al. (1998), though their analysis was retrospective rather than prospective as ours was. Ours was also data dredging (adjective) as we used the same study period with hopes of detecting the same cluster.

We can score each county according to how likely it is to be part of a cluster in a heuristic fashion using the function `score_locations`, and visualize the results on a heatmap as follows:

``` r
# Calculate scores and add column with county names
county_scores <- score_locations(poisson_result)
county_scores[, county := counties]

# Create a table for plotting
score_map_df <- merge(NM_map, county_scores, by = "county", all.x = TRUE)
score_map_df[subregion == "cibola", 
             relative_score := county_scores[county == "valencia", relative_score]]

ggplot() + 
  geom_polygon(data = score_map_df,
               mapping = aes(x = long, y = lat, group = group, 
                             fill = relative_score),
               color = "grey") +
  scale_fill_gradient(low = "#e5f5f9", high = "darkgreen",
                      guide = guide_colorbar(title = "Relative\nScore")) +
  geom_text(data = centroids, 
            mapping = aes(x = long, y = lat, label = subregion),
            alpha = 0.5) +
  ggtitle("County scores")
```

![](README_figures/unnamed-chunk-8-1.png)

As we can see, this does not match up entirely with the previous results as Torrance was not part of the most likely cluster.

References
==========

Kleinman, Ken. 2015. *Rsatscan: Tools, Classes, and Methods for Interfacing with SaTScan Stand-Alone Software*. <https://CRAN.R-project.org/package=rsatscan>.

Kulldorff, Martin, William F. Athas, Eric J. Feuer, Barry A. Miller, and Charles R. Key. 1998. “Evaluating Cluster Alarms: A Space-Time Scan Statistic and Brain Cancer in Los Alamos.” *American Journal of Public Health* 88 (9): 1377–80.

Neill, Daniel B, Andrew W Moore, Maheshkumar Sabhnani, and Kenny Daniel. 2005. “Detection of Emerging Space-Time Clusters.” In *Proceedings of the Eleventh ACM SIGKDD International Conference on Knowledge Discovery in Data Mining*, 218–27. ACM. doi:[10.1145/1081870.1081897](https://doi.org/10.1145/1081870.1081897).
