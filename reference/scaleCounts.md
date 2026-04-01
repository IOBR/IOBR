# Scaling raw counts from each sample.

Normalizing the sum of counts from each sample to 1e6.

## Usage

``` r
scaleCounts(counts, sigGenes = NULL, renormGenes = NULL, normFact = NULL)
```

## Details

Function taking a matrix (*genes* x *samples*) of counts as input and
returning the scaled counts for the subset signature genes (or all genes
if it sigGenes is `NULL`), with the scaled counts computed based on all
the 'renormGenes' (or all genes if it is NULL). The renormalization is
made independently for each sample, to have the sum of each columns over
the renormGenes equal to 1e6.

normFact, if not null, is used as the normalization factor instead of
the colSums (used to renormalize the refProfiles.var by the same amount
than the refProfiles).
