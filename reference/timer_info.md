# Source code for the TIMER deconvolution method.

Formats and displays informational messages for timing or logging
purposes. Useful for tracking progress or stages of execution within
scripts.

## Usage

``` r
timer_info(string)
```

## Arguments

- string:

  Character. Message to be displayed.

## Value

None; used for its side effect of printing a message.

## Details

This code is adapted from https://github.com/hanfeisun/TIMER, which
again is an adapted version of the original TIMER source code from
http://cistrome.org/TIMER/download.html.

The method is described in Li et al. Genome Biology 2016;17(1):174.,
PMID 27549193. Display Timer Information Messages

## Author

Bo Li

## Examples

``` r
timer_info("Data processing started.")
#> ℹ Data processing started.
```
