A simple tool for calculating peptide sequence and fragment masses.

## Installation

To install, open R and type:

``` r
install.packages("devtools")
devtools::install_github("jeffsocal/rmstandem")
```

## Get Started

Calculate the mass of a peptide sequence

``` r
library(rmstandem)

sequence_mass("SAMPLE")
#> [1] 647.3069
```

Calculate the mass of a peptide sequence with a mass modification

``` r
sequence_mass("SA[M15.998]PLE")
#> [1] 663.3049
```

Generate a table of expected fragment masses Calculate the mass of a
peptide sequence

``` r
fragments("SAMPLE", charge = 1)
#>     ion        mz z    seq pair pos      type
#> 1   b1+  88.03928 1      S  p01   1         b
#> 2   y1+ 148.06044 1      E  p05   1         y
#> 3   b2+ 159.07638 1     SA  p02   2         b
#> 4   y2+ 261.14454 1     LE  p04   2         y
#> 5   b3+ 290.11688 1    SAM  p03   3         b
#> 6  MH++ 324.15711 2 SAMPLE  p00  NA precursor
#> 7   y3+ 358.19734 1    PLE  p03   3         y
#> 8   b4+ 387.16968 1   SAMP  p04   4         b
#> 9   y4+ 489.23784 1   MPLE  p02   4         y
#> 10  b5+ 500.25378 1  SAMPL  p05   5         b
#> 11  y5+ 560.27494 1  AMPLE  p01   5         y
#> 12  MH+ 647.30694 1 SAMPLE  p00  NA precursor
```
