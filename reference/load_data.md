# Load IOBR Datasets

Loads internal datasets from the IOBR package. Supports both sysdata
(internal) and exported data files included in the package.

## Usage

``` r
load_data(name)
```

## Arguments

- name:

  Character string. Name of the dataset to load. Must be a single value.
  Available datasets include: - Expression data: \`"eset_stad"\`,
  \`"imvigor210_eset"\`, \`"melanoma_data"\` - Signatures:
  \`"signature_tme"\`, \`"signature_metabolism"\`,
  \`"signature_collection"\` - Gene sets: \`"hallmark"\`, \`"kegg"\`,
  \`"go_bp"\`, \`"go_cc"\`, \`"go_mf"\` - Cell markers:
  \`"cellmarkers"\`, \`"mcp_genes"\` - Phenotype data: \`"pdata_stad"\`,
  \`"pdata_sig_tme"\`, \`"pdata_acrg"\` - Reference data:
  \`"xCell.data"\`, \`"quantiseq_data"\`, \`"TRef"\`, \`"BRef"\` - Color
  palettes: \`"palette1"\`, \`"palette2"\`, \`"palette3"\`,
  \`"palette4"\`

## Value

Dataset object, typically a \`list\`, \`data.frame\`, or \`matrix\`. The
exact type depends on the requested dataset.

## Examples

``` r
# Load signature collection (stored in sysdata, no download)
sig_tme <- load_data("signature_tme")

# Load color palette (stored in sysdata, no download)
colors <- load_data("palette1")

# Error handling with suggestions for similar names
try(load_data("sign_tme")) # Will suggest "signature_tme"
#> Error in load_data("sign_tme") : 
#>   Dataset "sign_tme" not found in IOBR package.
#> ℹ Available datasets: BRef, PurityDataAffy, SI_geneset, TRef, anno_gc_vm32,
#>   anno_grch38, anno_hug133plus2, anno_illumina, anno_rnaseq, cancer_type_genes,
#>   cellmarkers, common_genes, deg, eset_blca, eset_gse62254, eset_stad,
#>   eset_tme_stad, go_bp, go_cc, go_mf, hallmark, immuneCuratedData,
#>   imvigor210_eset, imvigor210_pdata, imvigor210_sig, ips_gene_set, kegg,
#>   length_ensembl, lm22, mRNA_cell_default, mcp_genes, mcp_probesets,
#>   melanoma_data, msig_immune, msig_sc, mus_human_gene_symbol, null_models,
#>   onco_sig, palette1, palette2, palette3, palette4, panel_for_gene,
#>   panel_for_signature, patterns_to_na, pdata_GSE63557, pdata_acrg,
#>   pdata_sig_tme, pdata_sig_tme_binary, pdata_stad, pdata_tme_binary,
#>   quantiseq_data, reactome, sig_excel, sig_group, sig_stad,
#>   signature_collection, signature_collection_citation, signature_metabolism,
#>   signature_sc, signature_tme, signature_tumor, stad_group, subgroup_data,
#>   tcga_stad_pdata, tcga_stad_sig, tcga_stad_var, xCell.data

if (FALSE) { # \dontrun{
# Load expression data (triggers download from GitHub)
eset <- load_data("eset_stad")
} # }
```
