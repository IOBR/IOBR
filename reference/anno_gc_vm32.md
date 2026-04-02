# Mouse Genome Annotation (GC/VM32)

Mouse gene annotation dataset for genome assembly GRCm38/mm10,
containing gene features including GC content and gene length
information. Used for mouse transcriptomic data analysis.

## Usage

``` r
data(anno_gc_vm32)
```

## Format

A data frame with columns:

- id:

  Ensembl gene identifier (e.g., "ENSMUSG00000000001").

- eff_length:

  Effective gene length for TPM calculation.

- gc:

  GC content proportion (0-1 range).

- symbol:

  Official gene symbol.

- mgi_id:

  Mouse Genome Informatics identifier.

- gene_type:

  Gene type classification (e.g., "protein_coding", "lncRNA").

- start:

  Genomic start position.

- end:

  Genomic end position.

- transcript_id:

  Transcript identifier (mostly NA in this dataset).

- ont:

  Gene ontology information (mostly NA in this dataset).

## Source

Ensembl database for mouse GRCm38/mm10

## Examples

``` r
data(anno_gc_vm32)
head(anno_gc_vm32)
#>                   id eff_length        gc symbol      mgi_id      gene_type
#> 1 ENSMUSG00000000001       3262 0.4350092  Gnai3   MGI:95773 protein_coding
#> 2 ENSMUSG00000000003        902 0.3481153   Pbsn MGI:1860484 protein_coding
#> 3 ENSMUSG00000000028       3506 0.4962921  Cdc45 MGI:1338073 protein_coding
#> 4 ENSMUSG00000000031       2625 0.5588571    H19   MGI:95891         lncRNA
#> 5 ENSMUSG00000000037       6397 0.4377052  Scml2 MGI:1340042 protein_coding
#> 6 ENSMUSG00000000049       1594 0.5050188   Apoh   MGI:88058 protein_coding
#>       start       end transcript_id  ont
#> 1 108014596 108053462          <NA> <NA>
#> 2  76881507  76897229          <NA> <NA>
#> 3  18599197  18630737          <NA> <NA>
#> 4 142129262 142131886          <NA> <NA>
#> 5 159865521 160041209          <NA> <NA>
#> 6 108234180 108305222          <NA> <NA>
```
