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
