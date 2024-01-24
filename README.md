A simple tool for calculating peptide sequence and fragment masses.

## Installation

To install, open R and type:

``` r
install.packages("devtools")
devtools::install_github("jeffsocal/mspredictr")
```

## Get Started

Calculate the mass of a peptide sequence

``` r
library(mspredictr)

peptide_mass("SAMPLER")
#> [1] 802.4007
```

Calculate the mass of a peptide sequence with a mass modification

``` r
peptide_mass("SA[M15.998]PLER")
#> [1] 818.3987
```

Generate a table of expected fragment masses

``` r
fragments("SA[M15.998]PLER", charge = 1)
#>     ion        mz z             seq pair pos      type
#> 1   b1+  88.03931 1               S  p01   1         b
#> 2   b2+ 159.07643 1              SA  p02   2         b
#> 3   y1+ 175.11894 1               R  p06   1         y
#> 4   y2+ 304.16153 1              ER  p05   2         y
#> 5   b3+ 306.11491 1     SA[M15.998]  p03   3         b
#> 6   b4+ 403.16767 1    SA[M15.998]P  p04   4         b
#> 7  MH++ 410.20663 2 SA[M15.998]PLER  p00  NA precursor
#> 8   y3+ 417.24559 1             LER  p04   3         y
#> 9   y4+ 514.29836 1            PLER  p03   4         y
#> 10  b5+ 516.25173 1   SA[M15.998]PL  p05   5         b
#> 11  b6+ 645.29432 1  SA[M15.998]PLE  p06   6         b
#> 12  y5+ 661.33684 1   [M15.998]PLER  p02   5         y
#> 13  y6+ 732.37396 1  A[M15.998]PLER  p01   6         y
#> 14  MH+ 819.40599 1 SA[M15.998]PLER  p00  NA precursor
```

Read in spectra

``` r
ms2data <- path_to_example() |> read_spectra()
```

Plot spectrum with peak matching

``` r
ms2data[1,] |> plot_spectrum(peptides = 'HAVSEGTK')
```

![](docs/README-unnamed-chunk-7-1.png)
