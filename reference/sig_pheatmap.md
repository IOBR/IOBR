# Generate Heatmap for Signature Data

Creates a heatmap from signature data with grouping variables, offering
flexible options for colors, clustering, and output formats using
ComplexHeatmap.

## Usage

``` r
sig_pheatmap(
  input,
  feas,
  group,
  group2 = NULL,
  group3 = NULL,
  ID = "ID",
  path = NULL,
  cols1 = "random",
  cols2 = "random",
  cols3 = "random",
  seed = 54321,
  show_col = FALSE,
  palette1 = 1,
  palette2 = 2,
  palette3 = 3,
  cluster_cols = TRUE,
  palette_for_heatmape = 6,
  scale.matrix = TRUE,
  cellwidth = 1,
  cellheight = 9,
  show_colnames = FALSE,
  fig.type = "pdf",
  width = 6,
  height = NULL,
  file_name_prefix = 1
)
```

## Arguments

- input:

  Data frame with variables in columns.

- feas:

  Character vector. Feature names (columns) to include in heatmap.

- group:

  Character string. Column name for primary grouping variable.

- group2:

  Character string or \`NULL\`. Optional secondary grouping variable.

- group3:

  Character string or \`NULL\`. Optional tertiary grouping variable.

- ID:

  Character string. Column name for sample identifiers. Default is
  \`"ID"\`.

- path:

  Character string or \`NULL\`. Directory to save output files. Default
  creates \`"Marker-heatmap-average"\`.

- cols1:

  Character vector or \`"random"\` or \`"normal"\`. Colors for primary
  group. Default is \`"random"\`.

- cols2:

  Character vector or \`"random"\` or \`"normal"\`. Colors for secondary
  group. Default is \`"random"\`.

- cols3:

  Character vector or \`"random"\` or \`"normal"\`. Colors for tertiary
  group. Default is \`"random"\`.

- seed:

  Integer. Random seed for color generation. Default is \`54321\`.

- show_col:

  Logical indicating whether to display colors. Default is \`FALSE\`.

- palette1:

  Integer. Palette for primary group. Default is \`1\`.

- palette2:

  Integer. Palette for secondary group. Default is \`2\`.

- palette3:

  Integer. Palette for tertiary group. Default is \`3\`.

- cluster_cols:

  Logical indicating whether to cluster columns. Default is \`TRUE\`.

- palette_for_heatmape:

  Integer. Palette number for heatmap. Default is \`6\`.

- scale.matrix:

  Logical indicating whether to scale the matrix. Default is \`TRUE\`.

- cellwidth:

  Numeric. Width of each cell in points. Default is \`1\`.

- cellheight:

  Numeric. Height of each cell in points. Default is \`9\`.

- show_colnames:

  Logical indicating whether to show column names. Default is \`FALSE\`.

- fig.type:

  Character string. File format for saving. Default is \`"pdf"\`.

- width:

  Numeric. Width of saved figure in inches. Default is \`6\`.

- height:

  Numeric or \`NULL\`. Height of saved figure in inches. Calculated if
  \`NULL\`.

- file_name_prefix:

  Character or numeric. Prefix for saved file name. Default is \`1\`.

## Value

A list containing:

- p_anno:

  Annotation data frame

- p_cols:

  List of cluster colors

- plot:

  ComplexHeatmap object

- eset:

  Transformed expression matrix

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
tcga_stad_sig <- load_data("tcga_stad_sig")
tcga_stad_pdata <- load_data("tcga_stad_pdata")
input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
feas <- grep("MCPcounter", colnames(input), value = TRUE)
sig_pheatmap(
  input = input, feas = feas, group = "subtype",
  scale.matrix = TRUE, path = tempdir()
)
#> ℹ Heatmap palettes: 1 (pheatmap), 2 (peach), 3 (blues), 4 (virids), 5 (reds), 6 (RdBu), 7 (navy_firebrick), 8 (magma)
#> ℹ Random palettes: 1 (palette1), 2 (palette2), 3 (palette3), 4 (palette4)
#> ℹ Using random seed: 54321
#> $subtype
#>         EBV          GS         CIN         MSI 
#> "#8F7700FF"   "#386CB0" "#008B45FF"   "#e31a1c" 
#> 
#> ✔ Heatmap saved to: /tmp/Rtmp8wyU3P/1-pheatmap-subtype.pdf
#> $p_anno
#>              subtype
#> TCGA-3M-AB46    <NA>
#> TCGA-B7-5818     EBV
#> TCGA-B7-A5TI    <NA>
#> TCGA-B7-A5TJ    <NA>
#> TCGA-B7-A5TK    <NA>
#> TCGA-B7-A5TN    <NA>
#> TCGA-BR-4187      GS
#> TCGA-BR-4191     CIN
#> TCGA-BR-4201     MSI
#> TCGA-BR-4253     EBV
#> TCGA-BR-4256     MSI
#> TCGA-BR-4257     MSI
#> TCGA-BR-4267     CIN
#> TCGA-BR-4279      GS
#> TCGA-BR-4280     MSI
#> TCGA-BR-6452     MSI
#> TCGA-BR-6453      GS
#> TCGA-BR-6455     EBV
#> TCGA-BR-6456     CIN
#> TCGA-BR-6457      GS
#> TCGA-BR-6458     CIN
#> TCGA-BR-6563      GS
#> TCGA-BR-6564     CIN
#> TCGA-BR-6565     CIN
#> TCGA-BR-6566     MSI
#> TCGA-BR-6705      GS
#> TCGA-BR-6707     EBV
#> TCGA-BR-6709    <NA>
#> TCGA-BR-6801      GS
#> TCGA-BR-6802     MSI
#> TCGA-BR-6803      GS
#> TCGA-BR-6852     MSI
#> TCGA-BR-7196     EBV
#> TCGA-BR-7197     CIN
#> TCGA-BR-7704    <NA>
#> TCGA-BR-7707     MSI
#> TCGA-BR-7715     CIN
#> TCGA-BR-7716     CIN
#> TCGA-BR-7717     CIN
#> TCGA-BR-7722      GS
#> TCGA-BR-7723     CIN
#> TCGA-BR-7901     CIN
#> TCGA-BR-7957      GS
#> TCGA-BR-7958     EBV
#> TCGA-BR-7959     CIN
#> TCGA-BR-8058      GS
#> TCGA-BR-8059     MSI
#> TCGA-BR-8060    <NA>
#> TCGA-BR-8077     CIN
#> TCGA-BR-8080     CIN
#> TCGA-BR-8081     MSI
#> TCGA-BR-8284     MSI
#> TCGA-BR-8286     CIN
#> TCGA-BR-8289     CIN
#> TCGA-BR-8291     CIN
#> TCGA-BR-8295     CIN
#> TCGA-BR-8296     CIN
#> TCGA-BR-8297     CIN
#> TCGA-BR-8361     MSI
#> TCGA-BR-8364      GS
#> TCGA-BR-8365     CIN
#> TCGA-BR-8366     EBV
#> TCGA-BR-8367      GS
#> TCGA-BR-8368     MSI
#> TCGA-BR-8369     CIN
#> TCGA-BR-8371      GS
#> TCGA-BR-8372     MSI
#> TCGA-BR-8373     CIN
#> TCGA-BR-8380      GS
#> TCGA-BR-8381     EBV
#> TCGA-BR-8382     MSI
#> TCGA-BR-8384      GS
#> TCGA-BR-8483     CIN
#> TCGA-BR-8484     CIN
#> TCGA-BR-8485     CIN
#> TCGA-BR-8486     CIN
#> TCGA-BR-8487     MSI
#> TCGA-BR-8588      GS
#> TCGA-BR-8589     EBV
#> TCGA-BR-8590      GS
#> TCGA-BR-8591     MSI
#> TCGA-BR-8592      GS
#> TCGA-BR-8676     EBV
#> TCGA-BR-8677      GS
#> TCGA-BR-8678     CIN
#> TCGA-BR-8680      GS
#> TCGA-BR-8682     CIN
#> TCGA-BR-8683     CIN
#> TCGA-BR-8686     EBV
#> TCGA-BR-8687     CIN
#> TCGA-BR-8690     CIN
#> TCGA-BR-A44T      GS
#> TCGA-BR-A44U     CIN
#> TCGA-BR-A4CS     CIN
#> TCGA-BR-A4IV      GS
#> TCGA-BR-A4J4     EBV
#> TCGA-BR-A4J5      GS
#> TCGA-BR-A4J6      GS
#> TCGA-BR-A4J7      GS
#> TCGA-BR-A4J8     CIN
#> TCGA-BR-A4J9      GS
#> TCGA-BR-A4PF     CIN
#> TCGA-BR-A4QL     MSI
#> TCGA-CD-5798      GS
#> TCGA-CD-5799     CIN
#> TCGA-CD-5800     CIN
#> TCGA-CD-5801     EBV
#> TCGA-CD-5803      GS
#> TCGA-CD-5804     CIN
#> TCGA-CD-5813      GS
#> TCGA-CD-8524     CIN
#> TCGA-CD-8525     CIN
#> TCGA-CD-8526     CIN
#> TCGA-CD-8527     CIN
#> TCGA-CD-8528     CIN
#> TCGA-CD-8529     CIN
#> TCGA-CD-8530     CIN
#> TCGA-CD-8531      GS
#> TCGA-CD-8532      GS
#> TCGA-CD-8533    <NA>
#> TCGA-CD-8534     CIN
#> TCGA-CD-8535     CIN
#> TCGA-CD-A486     CIN
#> TCGA-CD-A487     CIN
#> TCGA-CD-A489     CIN
#> TCGA-CD-A48A     CIN
#> TCGA-CD-A48C     CIN
#> TCGA-CD-A4MG     MSI
#> TCGA-CD-A4MH     CIN
#> TCGA-CG-4301     CIN
#> TCGA-CG-4305     MSI
#> TCGA-CG-4306     MSI
#> TCGA-CG-4436     CIN
#> TCGA-CG-4437     CIN
#> TCGA-CG-4438     CIN
#> TCGA-CG-4440     CIN
#> TCGA-CG-4441     CIN
#> TCGA-CG-4443     CIN
#> TCGA-CG-4444     CIN
#> TCGA-CG-4460    <NA>
#> TCGA-CG-4465     MSI
#> TCGA-CG-4466     CIN
#> TCGA-CG-4469     CIN
#> TCGA-CG-4475     CIN
#> TCGA-CG-4477     CIN
#> TCGA-CG-5717      GS
#> TCGA-CG-5718     CIN
#> TCGA-CG-5719     CIN
#> TCGA-CG-5720      GS
#> TCGA-CG-5721     MSI
#> TCGA-CG-5722     EBV
#> TCGA-CG-5723     MSI
#> TCGA-CG-5724     CIN
#> TCGA-CG-5725    <NA>
#> TCGA-CG-5726     MSI
#> TCGA-CG-5732     CIN
#> TCGA-CG-5734      GS
#> TCGA-D7-5577     EBV
#> TCGA-D7-5578     CIN
#> TCGA-D7-6519     CIN
#> TCGA-D7-6520     CIN
#> TCGA-D7-6521     CIN
#> TCGA-D7-6522      GS
#> TCGA-D7-6524      GS
#> TCGA-D7-6525     CIN
#> TCGA-D7-6526     CIN
#> TCGA-D7-6527     CIN
#> TCGA-D7-6528     CIN
#> TCGA-D7-6815     CIN
#> TCGA-D7-6818     CIN
#> TCGA-D7-6822     CIN
#> TCGA-D7-8570     EBV
#> TCGA-D7-8572      GS
#> TCGA-D7-8573     EBV
#> TCGA-D7-8574      GS
#> TCGA-D7-8575     CIN
#> TCGA-D7-8576     CIN
#> TCGA-D7-8578     CIN
#> TCGA-D7-8579      GS
#> TCGA-D7-A4YU     CIN
#> TCGA-D7-A4YX     EBV
#> TCGA-D7-A4Z0      GS
#> TCGA-D7-A6EV    <NA>
#> TCGA-D7-A6EX    <NA>
#> TCGA-D7-A6EY    <NA>
#> TCGA-D7-A6EZ    <NA>
#> TCGA-D7-A6F0    <NA>
#> TCGA-D7-A6F2    <NA>
#> TCGA-D7-A747    <NA>
#> TCGA-D7-A748    <NA>
#> TCGA-D7-A74A    <NA>
#> TCGA-EQ-8122     CIN
#> TCGA-F1-6177     MSI
#> TCGA-F1-6874     MSI
#> TCGA-F1-6875     CIN
#> TCGA-F1-A448     MSI
#> TCGA-F1-A72C    <NA>
#> TCGA-FP-7735     CIN
#> TCGA-FP-7829     CIN
#> TCGA-FP-7916     EBV
#> TCGA-FP-7998     EBV
#> TCGA-FP-8099     CIN
#> TCGA-FP-8210      GS
#> TCGA-FP-8211     CIN
#> TCGA-FP-8631     CIN
#> TCGA-FP-A4BF     CIN
#> TCGA-FP-A8CX    <NA>
#> TCGA-FP-A9TM    <NA>
#> TCGA-HF-7132     MSI
#> TCGA-HF-7133     CIN
#> TCGA-HF-7134    <NA>
#> TCGA-HF-A5NB    <NA>
#> TCGA-HJ-7597     MSI
#> TCGA-HU-8238    <NA>
#> TCGA-HU-8244    <NA>
#> TCGA-HU-8249     CIN
#> TCGA-HU-8602     MSI
#> TCGA-HU-8604     CIN
#> TCGA-HU-8608     EBV
#> TCGA-HU-8610     CIN
#> TCGA-HU-A4G2     EBV
#> TCGA-HU-A4G3      GS
#> TCGA-HU-A4G8     MSI
#> TCGA-HU-A4G9     MSI
#> TCGA-HU-A4GC     CIN
#> TCGA-HU-A4GD     CIN
#> TCGA-HU-A4GF     CIN
#> TCGA-HU-A4GH     CIN
#> TCGA-HU-A4GJ    <NA>
#> TCGA-HU-A4GP     CIN
#> TCGA-HU-A4GQ     MSI
#> TCGA-HU-A4GT     MSI
#> TCGA-HU-A4GU     MSI
#> TCGA-HU-A4GX     MSI
#> TCGA-HU-A4GY      GS
#> TCGA-HU-A4H0     EBV
#> TCGA-HU-A4H2     CIN
#> TCGA-HU-A4H3     MSI
#> TCGA-HU-A4H4     MSI
#> TCGA-HU-A4H5     CIN
#> TCGA-HU-A4H6     CIN
#> TCGA-HU-A4H8     MSI
#> TCGA-HU-A4HB    <NA>
#> TCGA-HU-A4HD     CIN
#> TCGA-IN-7806     CIN
#> TCGA-IN-7808     CIN
#> TCGA-IN-8462     CIN
#> TCGA-IN-8663     CIN
#> TCGA-IN-A6RI    <NA>
#> TCGA-IN-A6RJ    <NA>
#> TCGA-IN-A6RL    <NA>
#> TCGA-IN-A6RN    <NA>
#> TCGA-IN-A6RR    <NA>
#> TCGA-IN-A6RS    <NA>
#> TCGA-IN-A7NR    <NA>
#> TCGA-IN-A7NT    <NA>
#> TCGA-IN-A7NU    <NA>
#> TCGA-IN-AB1V    <NA>
#> TCGA-IN-AB1X    <NA>
#> TCGA-IP-7968     CIN
#> TCGA-KB-A6F7    <NA>
#> TCGA-KB-A93G    <NA>
#> TCGA-KB-A93H    <NA>
#> TCGA-KB-A93J    <NA>
#> TCGA-MX-A5UG    <NA>
#> TCGA-MX-A5UJ    <NA>
#> TCGA-MX-A663    <NA>
#> TCGA-MX-A666    <NA>
#> TCGA-R5-A7O7    <NA>
#> TCGA-R5-A7ZE    <NA>
#> TCGA-R5-A7ZF    <NA>
#> TCGA-R5-A7ZI    <NA>
#> TCGA-R5-A7ZR    <NA>
#> TCGA-R5-A805    <NA>
#> TCGA-RD-A7BS    <NA>
#> TCGA-RD-A7BT    <NA>
#> TCGA-RD-A7BW    <NA>
#> TCGA-RD-A7C1    <NA>
#> TCGA-RD-A8MV    <NA>
#> TCGA-RD-A8MW    <NA>
#> TCGA-RD-A8N0    <NA>
#> TCGA-RD-A8N1    <NA>
#> TCGA-RD-A8N4    <NA>
#> TCGA-RD-A8N5    <NA>
#> TCGA-RD-A8N6    <NA>
#> TCGA-RD-A8N9    <NA>
#> TCGA-RD-A8NB    <NA>
#> TCGA-SW-A7EA    <NA>
#> TCGA-SW-A7EB    <NA>
#> TCGA-VQ-A8DT    <NA>
#> TCGA-VQ-A8DU    <NA>
#> TCGA-VQ-A8DV    <NA>
#> TCGA-VQ-A8DZ    <NA>
#> TCGA-VQ-A8E0    <NA>
#> TCGA-VQ-A8E2    <NA>
#> TCGA-VQ-A8E3    <NA>
#> TCGA-VQ-A8E7    <NA>
#> TCGA-VQ-A8P2    <NA>
#> TCGA-VQ-A8P3    <NA>
#> TCGA-VQ-A8P5    <NA>
#> TCGA-VQ-A8P8    <NA>
#> TCGA-VQ-A8PB    <NA>
#> TCGA-VQ-A8PC    <NA>
#> TCGA-VQ-A8PD    <NA>
#> TCGA-VQ-A8PE    <NA>
#> TCGA-VQ-A8PF    <NA>
#> TCGA-VQ-A8PH    <NA>
#> TCGA-VQ-A8PJ    <NA>
#> TCGA-VQ-A8PK    <NA>
#> TCGA-VQ-A8PM    <NA>
#> TCGA-VQ-A8PO    <NA>
#> TCGA-VQ-A8PP    <NA>
#> TCGA-VQ-A8PQ    <NA>
#> TCGA-VQ-A8PU    <NA>
#> TCGA-VQ-A8PX    <NA>
#> TCGA-VQ-A91A    <NA>
#> TCGA-VQ-A91D    <NA>
#> TCGA-VQ-A91E    <NA>
#> TCGA-VQ-A91K    <NA>
#> TCGA-VQ-A91N    <NA>
#> TCGA-VQ-A91Q    <NA>
#> TCGA-VQ-A91S    <NA>
#> TCGA-VQ-A91U    <NA>
#> TCGA-VQ-A91V    <NA>
#> TCGA-VQ-A91X    <NA>
#> TCGA-VQ-A91Y    <NA>
#> TCGA-VQ-A91Z    <NA>
#> TCGA-VQ-A922    <NA>
#> TCGA-VQ-A923    <NA>
#> TCGA-VQ-A924    <NA>
#> TCGA-VQ-A925    <NA>
#> TCGA-VQ-A927    <NA>
#> TCGA-VQ-A928    <NA>
#> TCGA-VQ-A92D    <NA>
#> TCGA-VQ-A94O    <NA>
#> TCGA-VQ-A94P    <NA>
#> TCGA-VQ-A94R    <NA>
#> TCGA-VQ-A94T    <NA>
#> TCGA-VQ-A94U    <NA>
#> TCGA-VQ-AA64    <NA>
#> TCGA-VQ-AA68    <NA>
#> TCGA-VQ-AA69    <NA>
#> TCGA-VQ-AA6A    <NA>
#> TCGA-VQ-AA6D    <NA>
#> TCGA-VQ-AA6F    <NA>
#> TCGA-VQ-AA6G    <NA>
#> TCGA-VQ-AA6J    <NA>
#> TCGA-VQ-AA6K    <NA>
#> TCGA-ZA-A8F6    <NA>
#> TCGA-ZQ-A9CR    <NA>
#> 
#> $p_cols
#> $p_cols$subtype
#>         EBV          GS         CIN         MSI 
#> "#8F7700FF"   "#386CB0" "#008B45FF"   "#e31a1c" 
#> 
#> 
#> $plot

#> 
#> $eset
#>                                    TCGA-3M-AB46 TCGA-B7-5818 TCGA-B7-A5TI
#> T_cells_MCPcounter                   -1.5293778    0.1815390  -0.54458690
#> CD8_T_cells_MCPcounter               -1.2988159    0.3727035  -0.73750686
#> Cytotoxic_lymphocytes_MCPcounter     -1.4976838    0.1669304  -0.54708913
#> NK_cells_MCPcounter                  -1.3918835   -0.5790658  -0.59383342
#> B_lineage_MCPcounter                 -1.3435369   -1.0299282   0.32677556
#> Monocytic_lineage_MCPcounter         -1.5854336    0.6810721   0.03354702
#> Myeloid_dendritic_cells_MCPcounter   -1.4598778    0.2970510  -0.65255430
#> Neutrophils_MCPcounter               -1.1347084   -0.3320885  -1.49655644
#> Endothelial_cells_MCPcounter         -1.3776122   -1.0668018   1.16973320
#> Fibroblasts_MCPcounter               -0.6958131   -0.2625985   0.50803612
#>                                    TCGA-B7-A5TJ TCGA-B7-A5TK TCGA-B7-A5TN
#> T_cells_MCPcounter                   -1.4385270    1.2625325   -0.7143300
#> CD8_T_cells_MCPcounter               -1.4211184    1.9710453   -0.8108444
#> Cytotoxic_lymphocytes_MCPcounter     -1.3869823    2.4308313   -0.8344176
#> NK_cells_MCPcounter                  -1.4522585    2.3950048   -1.0414019
#> B_lineage_MCPcounter                 -1.0941828   -0.4824224   -0.7715137
#> Monocytic_lineage_MCPcounter         -1.2048338    1.3906292    0.1228352
#> Myeloid_dendritic_cells_MCPcounter   -1.0133761    0.2879003   -0.6457672
#> Neutrophils_MCPcounter               -1.1135425   -0.6906884   -0.3962743
#> Endothelial_cells_MCPcounter         -0.6265681    0.1064296    0.2004661
#> Fibroblasts_MCPcounter               -0.5093713    0.8000789    1.4392078
#>                                    TCGA-BR-4187 TCGA-BR-4191 TCGA-BR-4201
#> T_cells_MCPcounter                  -0.07286809    0.9087963    0.2065507
#> CD8_T_cells_MCPcounter               0.62867029    0.4463726   -0.9600922
#> Cytotoxic_lymphocytes_MCPcounter     0.07454271    0.7296823    0.7160903
#> NK_cells_MCPcounter                 -0.27262007    0.9253402    1.7445752
#> B_lineage_MCPcounter                 0.16354694    0.5231614    0.3174053
#> Monocytic_lineage_MCPcounter         1.57507576    2.2492433    1.2689643
#> Myeloid_dendritic_cells_MCPcounter   0.77658420    3.5299029    0.3216347
#> Neutrophils_MCPcounter               1.60569553    1.0031462    1.8783557
#> Endothelial_cells_MCPcounter         1.24735034    1.2575540    0.4631374
#> Fibroblasts_MCPcounter               1.69277731    0.2293574    1.3286410
#>                                    TCGA-BR-4253 TCGA-BR-4256 TCGA-BR-4257
#> T_cells_MCPcounter                   2.14282020  0.801661296  -0.08688099
#> CD8_T_cells_MCPcounter               1.82595997  0.833118093   0.15960153
#> Cytotoxic_lymphocytes_MCPcounter     3.53516365  1.161595076   0.41770432
#> NK_cells_MCPcounter                  1.66763427  0.685363139  -0.26662888
#> B_lineage_MCPcounter                -0.03397195 -0.009119045  -0.98684457
#> Monocytic_lineage_MCPcounter         2.73193815  2.028672660   1.12542959
#> Myeloid_dendritic_cells_MCPcounter   1.05770647  0.253222410   1.31393125
#> Neutrophils_MCPcounter               0.51282480  2.554853720   0.59316390
#> Endothelial_cells_MCPcounter        -0.18147421  1.310472183  -0.40595047
#> Fibroblasts_MCPcounter              -0.78230684  1.440230956   0.58084973
#>                                    TCGA-BR-4267 TCGA-BR-4279 TCGA-BR-4280
#> T_cells_MCPcounter                   -0.8912010 -0.001296296   0.08013704
#> CD8_T_cells_MCPcounter               -0.7605730  0.903684753  -0.19152041
#> Cytotoxic_lymphocytes_MCPcounter     -0.9093449  0.938429071   0.26023941
#> NK_cells_MCPcounter                  -0.9070585  1.186445410  -0.34490345
#> B_lineage_MCPcounter                 -0.1925955 -0.200696949  -0.57689723
#> Monocytic_lineage_MCPcounter          0.3757014  1.088680676   0.16923193
#> Myeloid_dendritic_cells_MCPcounter   -0.6344686  0.940639365  -0.75375316
#> Neutrophils_MCPcounter                1.3675192  1.006503773  -0.03702467
#> Endothelial_cells_MCPcounter         -0.6506859  1.051844979  -1.27659393
#> Fibroblasts_MCPcounter               -0.5092481  1.379247916  -1.01470014
#>                                    TCGA-BR-6452 TCGA-BR-6453 TCGA-BR-6455
#> T_cells_MCPcounter                   0.44666801    1.1026645   0.48530716
#> CD8_T_cells_MCPcounter              -0.03938577    0.5857928   1.26769243
#> Cytotoxic_lymphocytes_MCPcounter     1.26625457   -0.1190170   0.73035145
#> NK_cells_MCPcounter                  1.12461232   -0.6917539   1.62739130
#> B_lineage_MCPcounter                -0.15054849    2.5549877  -0.05177583
#> Monocytic_lineage_MCPcounter         0.82479110    0.2843303   1.20937259
#> Myeloid_dendritic_cells_MCPcounter   1.17999976    1.6666543   1.82996134
#> Neutrophils_MCPcounter              -0.26017806   -1.2556234  -0.17779579
#> Endothelial_cells_MCPcounter        -0.04879316    0.4021866   0.31849559
#> Fibroblasts_MCPcounter               0.48180331    0.1623347   0.23894445
#>                                    TCGA-BR-6456 TCGA-BR-6457 TCGA-BR-6458
#> T_cells_MCPcounter                  0.573167103   0.02289880    0.9134304
#> CD8_T_cells_MCPcounter              0.128034012   0.25811769    0.4939806
#> Cytotoxic_lymphocytes_MCPcounter    0.195213423   0.14786524    0.7920014
#> NK_cells_MCPcounter                -0.327419320  -0.14646287    0.4944324
#> B_lineage_MCPcounter                0.165561126   1.05969214    0.2017638
#> Monocytic_lineage_MCPcounter        0.855268112  -0.07766011    1.3636574
#> Myeloid_dendritic_cells_MCPcounter  0.001841283   0.88425703    0.1695385
#> Neutrophils_MCPcounter              0.808473974  -0.46539653   -0.4774011
#> Endothelial_cells_MCPcounter        1.487256021   1.48103593    0.8293872
#> Fibroblasts_MCPcounter              1.630915954   1.89698874    0.9967863
#>                                    TCGA-BR-6563 TCGA-BR-6564 TCGA-BR-6565
#> T_cells_MCPcounter                    0.7070845  0.207783747    0.6289755
#> CD8_T_cells_MCPcounter                1.3711111  0.087208972    0.4083453
#> Cytotoxic_lymphocytes_MCPcounter      0.7011896 -0.118983160    0.7889203
#> NK_cells_MCPcounter                  -0.1349577 -0.981120478    0.2143991
#> B_lineage_MCPcounter                  1.4702845 -0.006451999    0.4851531
#> Monocytic_lineage_MCPcounter          0.3403385 -0.297761725    0.8343310
#> Myeloid_dendritic_cells_MCPcounter    0.7518076  0.258308181   -0.6658379
#> Neutrophils_MCPcounter               -0.7939110 -0.996702808   -0.7305177
#> Endothelial_cells_MCPcounter          0.4931703  1.360259783    0.1248542
#> Fibroblasts_MCPcounter                1.0542571  1.573761663    0.7386146
#>                                    TCGA-BR-6566 TCGA-BR-6705 TCGA-BR-6707
#> T_cells_MCPcounter                  -0.69882044  -0.06059829    0.5554160
#> CD8_T_cells_MCPcounter              -0.38242361  -0.31142308    1.5747207
#> Cytotoxic_lymphocytes_MCPcounter     0.01985942   0.17901636    1.2348222
#> NK_cells_MCPcounter                  0.84859453  -0.12523843    0.8627753
#> B_lineage_MCPcounter                -0.41060761   0.33626592   -0.1247874
#> Monocytic_lineage_MCPcounter        -0.06122360  -0.10348009    0.3426645
#> Myeloid_dendritic_cells_MCPcounter   0.10206009   0.38769652   -0.8660467
#> Neutrophils_MCPcounter              -1.72521153   0.56790524   -0.9395928
#> Endothelial_cells_MCPcounter        -1.59919755   2.74429823   -0.5382930
#> Fibroblasts_MCPcounter               0.77389964   1.93148857   -0.6460671
#>                                    TCGA-BR-6709 TCGA-BR-6801 TCGA-BR-6802
#> T_cells_MCPcounter                   1.40199651   -1.5351422   0.63986498
#> CD8_T_cells_MCPcounter               1.82097457   -1.2038139   0.85245986
#> Cytotoxic_lymphocytes_MCPcounter     1.33868867   -1.1222615   1.79565328
#> NK_cells_MCPcounter                  0.79064028   -1.4059104   2.00051095
#> B_lineage_MCPcounter                 0.74657668   -1.5628081   0.62583764
#> Monocytic_lineage_MCPcounter         2.15144036   -1.4916590   0.65871192
#> Myeloid_dendritic_cells_MCPcounter   0.87868166   -1.8357242   0.26289623
#> Neutrophils_MCPcounter              -0.08931934   -0.5409907   0.52168669
#> Endothelial_cells_MCPcounter        -0.21855987   -0.7355018   0.03272997
#> Fibroblasts_MCPcounter               0.32301115    0.3304275   0.29283574
#>                                    TCGA-BR-6803 TCGA-BR-6852 TCGA-BR-7196
#> T_cells_MCPcounter                   0.44789955    1.4894734    0.4751209
#> CD8_T_cells_MCPcounter               0.82785847    2.1587729    1.0256156
#> Cytotoxic_lymphocytes_MCPcounter     0.60238534    2.9675873    1.0964372
#> NK_cells_MCPcounter                 -0.18775434    2.6303872    1.6181225
#> B_lineage_MCPcounter                 0.83242258    0.4457565    0.9288396
#> Monocytic_lineage_MCPcounter         0.05303055    1.3508807    0.8396502
#> Myeloid_dendritic_cells_MCPcounter   0.65293755    0.3755038    0.5416161
#> Neutrophils_MCPcounter              -0.84470221   -0.5864361   -0.3919438
#> Endothelial_cells_MCPcounter         1.42930210   -0.1020505    2.0425445
#> Fibroblasts_MCPcounter               1.15285869    0.4647914    1.4687773
#>                                    TCGA-BR-7197 TCGA-BR-7704 TCGA-BR-7707
#> T_cells_MCPcounter                  -1.07794857  1.794657275  -0.60347079
#> CD8_T_cells_MCPcounter              -1.41926305  0.800584983  -0.39740898
#> Cytotoxic_lymphocytes_MCPcounter    -1.58119800  1.385652700   0.41736724
#> NK_cells_MCPcounter                 -1.21661887  0.633561463   0.81203744
#> B_lineage_MCPcounter                -1.09205094  0.613685839  -0.68453711
#> Monocytic_lineage_MCPcounter        -1.39736473  0.405752229  -0.68186841
#> Myeloid_dendritic_cells_MCPcounter  -1.31169929  1.852754468  -1.19067849
#> Neutrophils_MCPcounter              -1.40635570 -0.916340044   0.03022563
#> Endothelial_cells_MCPcounter        -0.34960726 -0.366605250  -0.63977502
#> Fibroblasts_MCPcounter               0.08319044  0.005296758   0.28332018
#>                                    TCGA-BR-7715 TCGA-BR-7716 TCGA-BR-7717
#> T_cells_MCPcounter                   -1.0518739    1.4190761   -1.5367146
#> CD8_T_cells_MCPcounter               -1.1598674    2.6013641    0.1377626
#> Cytotoxic_lymphocytes_MCPcounter     -1.2605295    1.1685892   -1.0708182
#> NK_cells_MCPcounter                  -1.2943659    0.1444851   -1.3176935
#> B_lineage_MCPcounter                 -1.3828910    0.4702353   -0.1228357
#> Monocytic_lineage_MCPcounter         -0.5386820    1.3724100   -1.5490466
#> Myeloid_dendritic_cells_MCPcounter   -1.4807273    0.3864229   -1.4682597
#> Neutrophils_MCPcounter               -1.2885608   -0.1437321   -1.2055586
#> Endothelial_cells_MCPcounter          0.2433338    0.1030076   -0.3025693
#> Fibroblasts_MCPcounter                0.7300876    0.3301410    0.2390101
#>                                    TCGA-BR-7722 TCGA-BR-7723 TCGA-BR-7901
#> T_cells_MCPcounter                  -0.63440065  0.484333901   0.17967437
#> CD8_T_cells_MCPcounter              -0.86049421  1.825937331  -0.03315405
#> Cytotoxic_lymphocytes_MCPcounter    -0.97334817  0.401026879  -0.02389751
#> NK_cells_MCPcounter                 -0.61840445 -0.002982094  -0.50204027
#> B_lineage_MCPcounter                -0.36914836  0.320187033  -0.08344899
#> Monocytic_lineage_MCPcounter        -0.48284810  0.689982055   0.62392797
#> Myeloid_dendritic_cells_MCPcounter  -0.35457510 -0.155287971   1.73279771
#> Neutrophils_MCPcounter              -0.73764760 -0.981949055   0.85679781
#> Endothelial_cells_MCPcounter        -0.25876151  0.297038078   1.10995680
#> Fibroblasts_MCPcounter              -0.06806234  0.087815949   0.93780839
#>                                    TCGA-BR-7957 TCGA-BR-7958 TCGA-BR-7959
#> T_cells_MCPcounter                  -0.46052054   1.74117999   -0.7419238
#> CD8_T_cells_MCPcounter              -0.19882603   2.03011397    0.1473773
#> Cytotoxic_lymphocytes_MCPcounter    -0.05620313   1.99569147   -0.7986916
#> NK_cells_MCPcounter                 -0.45039357   0.91084379   -0.9201722
#> B_lineage_MCPcounter                -0.75677458   0.09273172   -0.5482823
#> Monocytic_lineage_MCPcounter         0.24282620   1.52155339    0.7023745
#> Myeloid_dendritic_cells_MCPcounter   0.34178384   2.35514099   -0.4339359
#> Neutrophils_MCPcounter               1.66127070  -0.70547388    0.3067386
#> Endothelial_cells_MCPcounter         2.23752866   0.11484486    1.2680504
#> Fibroblasts_MCPcounter               1.89216781   0.63277311    1.4046253
#>                                    TCGA-BR-8058 TCGA-BR-8059 TCGA-BR-8060
#> T_cells_MCPcounter                    1.5418653  -1.11677774  0.006729937
#> CD8_T_cells_MCPcounter                0.7786092  -0.97043144 -0.672262898
#> Cytotoxic_lymphocytes_MCPcounter      1.3572579  -0.55381821  0.138412589
#> NK_cells_MCPcounter                   0.2187054  -0.49650847  0.435267475
#> B_lineage_MCPcounter                  1.3287288  -0.57437040  0.539603822
#> Monocytic_lineage_MCPcounter          1.4305193  -0.97848338  1.006992227
#> Myeloid_dendritic_cells_MCPcounter    1.6328197  -1.37996347  0.300636977
#> Neutrophils_MCPcounter                1.3832755  -0.09834598  0.231555795
#> Endothelial_cells_MCPcounter          0.6194573   0.48342144  1.126676260
#> Fibroblasts_MCPcounter                0.6621921   1.43847709  0.722602592
#>                                    TCGA-BR-8077 TCGA-BR-8080 TCGA-BR-8081
#> T_cells_MCPcounter                  0.294680176    0.4056478   1.09508358
#> CD8_T_cells_MCPcounter              2.280187403   -0.1274944   0.48764417
#> Cytotoxic_lymphocytes_MCPcounter    0.181009669    0.6392880   0.84659085
#> NK_cells_MCPcounter                -0.146965312    0.3548939   0.33783596
#> B_lineage_MCPcounter               -0.300079836    1.2989041   1.25579282
#> Monocytic_lineage_MCPcounter        0.089355443    0.7963234   1.34686821
#> Myeloid_dendritic_cells_MCPcounter  0.874064634    0.7214899   0.85466029
#> Neutrophils_MCPcounter             -0.569942571    0.1221878  -0.50482633
#> Endothelial_cells_MCPcounter        0.004115349    1.7256700   0.09695537
#> Fibroblasts_MCPcounter             -0.672068757    1.4790869   0.62868406
#>                                     TCGA-BR-8284 TCGA-BR-8286 TCGA-BR-8289
#> T_cells_MCPcounter                  1.4697302980  -0.22498483   -0.9321558
#> CD8_T_cells_MCPcounter              1.6386523943  -0.10632080   -0.4045234
#> Cytotoxic_lymphocytes_MCPcounter    2.0573877416  -0.02767038   -0.9071876
#> NK_cells_MCPcounter                 0.7062732030   0.01328003   -0.9345912
#> B_lineage_MCPcounter                0.4131924091  -0.58039655   -0.2349604
#> Monocytic_lineage_MCPcounter        1.2696009891   0.20083148    0.1745647
#> Myeloid_dendritic_cells_MCPcounter  1.0747809522   1.49074963   -0.8967442
#> Neutrophils_MCPcounter             -0.0005504517  -0.94764878   -1.6504979
#> Endothelial_cells_MCPcounter        0.9522432050   0.81494527   -0.1679650
#> Fibroblasts_MCPcounter              0.7247819036   0.04518764    0.5992355
#>                                    TCGA-BR-8291 TCGA-BR-8295 TCGA-BR-8296
#> T_cells_MCPcounter                   0.42686797  -1.52241546   1.26062306
#> CD8_T_cells_MCPcounter               0.48341055  -1.01390221   0.19283590
#> Cytotoxic_lymphocytes_MCPcounter     0.69560555  -0.77239612   0.19740314
#> NK_cells_MCPcounter                  0.86507362  -0.22285692   1.32800307
#> B_lineage_MCPcounter                -0.19956781  -0.05770472   2.11594718
#> Monocytic_lineage_MCPcounter         1.36952137  -1.52888713   0.28397497
#> Myeloid_dendritic_cells_MCPcounter   1.86766743  -1.52592344   1.25807702
#> Neutrophils_MCPcounter              -0.04225999  -1.06785330   0.09371892
#> Endothelial_cells_MCPcounter         1.89962916  -0.37317805   0.65960333
#> Fibroblasts_MCPcounter               1.53210353  -1.31062134  -0.17888056
#>                                    TCGA-BR-8297 TCGA-BR-8361 TCGA-BR-8364
#> T_cells_MCPcounter                   -0.1135552    0.3388134   0.59354438
#> CD8_T_cells_MCPcounter               -0.3880230    0.9571575   0.40227894
#> Cytotoxic_lymphocytes_MCPcounter     -0.2884956    1.1055172   0.85341440
#> NK_cells_MCPcounter                  -0.3072789    2.3377186   0.42388592
#> B_lineage_MCPcounter                  0.7035169   -0.4990919   0.59217368
#> Monocytic_lineage_MCPcounter         -0.4569697   -0.3697003   1.26348820
#> Myeloid_dendritic_cells_MCPcounter   -0.2168401   -1.4631623   1.52078036
#> Neutrophils_MCPcounter               -1.1036778   -0.6402848  -0.06347479
#> Endothelial_cells_MCPcounter          0.6519535   -0.5493772   1.41882215
#> Fibroblasts_MCPcounter                1.0976153   -0.7659661   2.27161985
#>                                    TCGA-BR-8365 TCGA-BR-8366 TCGA-BR-8367
#> T_cells_MCPcounter                   0.69176554    1.7577630   -0.3240666
#> CD8_T_cells_MCPcounter               0.42063454    1.7229871   -0.3635681
#> Cytotoxic_lymphocytes_MCPcounter     0.22900528    2.0492441    0.1117893
#> NK_cells_MCPcounter                 -0.47927772    3.4346969    0.1602634
#> B_lineage_MCPcounter                 0.89474860    1.4267474   -0.2047782
#> Monocytic_lineage_MCPcounter         0.40016535    1.8184780    0.6785989
#> Myeloid_dendritic_cells_MCPcounter   0.46304418    0.2906231    0.3898941
#> Neutrophils_MCPcounter              -0.05132177   -0.3276903    0.4148624
#> Endothelial_cells_MCPcounter         1.81879068    0.2413490    1.6380115
#> Fibroblasts_MCPcounter               1.44195342    0.4733461    1.0824236
#>                                    TCGA-BR-8368 TCGA-BR-8369 TCGA-BR-8371
#> T_cells_MCPcounter                   0.22885717   -0.6725619   0.28798458
#> CD8_T_cells_MCPcounter              -0.02968514   -0.6200191  -0.41162437
#> Cytotoxic_lymphocytes_MCPcounter    -0.04499447   -0.9019559  -0.73145641
#> NK_cells_MCPcounter                  0.40395703    0.1152733  -0.78269675
#> B_lineage_MCPcounter                -0.13298017    0.5743837   1.53326500
#> Monocytic_lineage_MCPcounter        -0.91522285   -0.9573186  -1.10298753
#> Myeloid_dendritic_cells_MCPcounter  -0.73110996   -0.7078119  -0.52018810
#> Neutrophils_MCPcounter              -0.81393660   -0.6048896  -1.08872319
#> Endothelial_cells_MCPcounter         0.03154735    0.6600474  -0.05948564
#> Fibroblasts_MCPcounter              -1.09110543    0.4506656   0.76788246
#>                                    TCGA-BR-8372 TCGA-BR-8373 TCGA-BR-8380
#> T_cells_MCPcounter                    0.8017442   0.06711858  -0.23779340
#> CD8_T_cells_MCPcounter                0.7983046   0.16972095  -0.70296721
#> Cytotoxic_lymphocytes_MCPcounter      1.2435244  -0.21272762  -0.48464578
#> NK_cells_MCPcounter                   1.0139598  -0.47956054  -0.08909080
#> B_lineage_MCPcounter                 -0.2418535  -0.51338001   0.33138988
#> Monocytic_lineage_MCPcounter          1.0114457   0.77601774  -0.04263330
#> Myeloid_dendritic_cells_MCPcounter    2.9228154   0.58463992  -0.06907682
#> Neutrophils_MCPcounter               -1.1682452   2.12162390  -0.22952687
#> Endothelial_cells_MCPcounter         -0.8826809   1.06711323   0.83983431
#> Fibroblasts_MCPcounter               -0.5672966   0.97739566   1.85995058
#>                                    TCGA-BR-8381 TCGA-BR-8382 TCGA-BR-8384
#> T_cells_MCPcounter                    1.6029087  -0.52527531    0.9777396
#> CD8_T_cells_MCPcounter                1.0990967  -1.09390884    0.8300289
#> Cytotoxic_lymphocytes_MCPcounter      1.2093536   0.44578130    0.8006242
#> NK_cells_MCPcounter                   0.6719895   1.65163989    0.4689062
#> B_lineage_MCPcounter                  1.1008072  -0.31286933    0.1944507
#> Monocytic_lineage_MCPcounter          1.0891902   0.08613423    0.5190539
#> Myeloid_dendritic_cells_MCPcounter    0.7249968   0.03756457    1.1866681
#> Neutrophils_MCPcounter               -0.9357035   0.45598787   -0.3563972
#> Endothelial_cells_MCPcounter         -0.2291877  -0.50166112    1.8578587
#> Fibroblasts_MCPcounter                0.2015060   0.67394808    1.6619366
#>                                    TCGA-BR-8483 TCGA-BR-8484 TCGA-BR-8485
#> T_cells_MCPcounter                   -1.5038883  1.398985483    1.0162428
#> CD8_T_cells_MCPcounter               -0.8438450  0.361427828    0.8824894
#> Cytotoxic_lymphocytes_MCPcounter     -1.1789675  0.844062035    1.3570818
#> NK_cells_MCPcounter                  -0.7615124 -0.190177201    0.9150149
#> B_lineage_MCPcounter                 -1.0056094  0.714326458    0.3743291
#> Monocytic_lineage_MCPcounter         -1.7989118  1.321498932    0.4303497
#> Myeloid_dendritic_cells_MCPcounter   -1.8192122  1.133180393    0.2639193
#> Neutrophils_MCPcounter               -0.9774519  0.506929206    0.3879810
#> Endothelial_cells_MCPcounter         -1.7069644  0.562806892    1.1343523
#> Fibroblasts_MCPcounter               -0.6471484  0.001741291    0.7004725
#>                                    TCGA-BR-8486 TCGA-BR-8487 TCGA-BR-8588
#> T_cells_MCPcounter                   0.86278947   -0.8586571   0.32254449
#> CD8_T_cells_MCPcounter              -0.99809049   -0.7746100   0.54129531
#> Cytotoxic_lymphocytes_MCPcounter    -0.41721331   -0.3995426   0.49440050
#> NK_cells_MCPcounter                 -0.38903863    0.5587517   0.08144136
#> B_lineage_MCPcounter                -0.02687122   -0.7932142   0.66567629
#> Monocytic_lineage_MCPcounter        -0.03046340   -0.4534154   0.20148469
#> Myeloid_dendritic_cells_MCPcounter   0.07140182   -1.0294591   1.33981775
#> Neutrophils_MCPcounter               0.61494601   -1.8489602   1.05219171
#> Endothelial_cells_MCPcounter         1.36280453   -1.4204676   0.29908371
#> Fibroblasts_MCPcounter               0.29890308   -0.6805069   0.65682371
#>                                    TCGA-BR-8589 TCGA-BR-8590 TCGA-BR-8591
#> T_cells_MCPcounter                    1.4355793    0.7354713   0.18245769
#> CD8_T_cells_MCPcounter                2.1283086    0.2613703  -0.98757925
#> Cytotoxic_lymphocytes_MCPcounter      2.4612022    0.5279011   0.08231699
#> NK_cells_MCPcounter                   2.6638907    0.1814867   0.06513350
#> B_lineage_MCPcounter                 -0.6559567    0.8012439   0.48102119
#> Monocytic_lineage_MCPcounter          1.2216310    1.4200771   0.65176460
#> Myeloid_dendritic_cells_MCPcounter    2.0741179    1.3491494   0.73752023
#> Neutrophils_MCPcounter               -0.5466511    0.3743000  -0.24497388
#> Endothelial_cells_MCPcounter         -0.4273732    0.8547585   0.71286902
#> Fibroblasts_MCPcounter               -1.6224903    1.3959263   0.41943963
#>                                    TCGA-BR-8592 TCGA-BR-8676 TCGA-BR-8677
#> T_cells_MCPcounter                  0.510804607  -0.66766697    1.2119117
#> CD8_T_cells_MCPcounter             -0.206204520   0.58356193    1.4015145
#> Cytotoxic_lymphocytes_MCPcounter    0.134304026  -0.08961074    2.0534560
#> NK_cells_MCPcounter                -0.023883569  -0.10503561    1.1841230
#> B_lineage_MCPcounter                1.348991354  -1.28271978   -0.3194401
#> Monocytic_lineage_MCPcounter        0.170805744  -0.29835457    0.3553119
#> Myeloid_dendritic_cells_MCPcounter  0.244401540  -0.79115067    0.1493601
#> Neutrophils_MCPcounter              0.003119563  -0.99689924   -0.3477128
#> Endothelial_cells_MCPcounter        1.302925103  -1.51900159    1.8681157
#> Fibroblasts_MCPcounter              1.565790495  -2.10838126    1.0646099
#>                                    TCGA-BR-8678 TCGA-BR-8680 TCGA-BR-8682
#> T_cells_MCPcounter                   0.05707897   -1.4314240   0.03241836
#> CD8_T_cells_MCPcounter               2.04584082   -0.9187353  -0.44072668
#> Cytotoxic_lymphocytes_MCPcounter    -0.03297288   -1.0932945  -0.42458570
#> NK_cells_MCPcounter                  0.35112546   -0.1759446  -0.17398051
#> B_lineage_MCPcounter                -0.87452692   -0.7244482   0.30389725
#> Monocytic_lineage_MCPcounter        -0.41093580   -1.4429634   0.99334821
#> Myeloid_dendritic_cells_MCPcounter  -0.91912158   -1.7845362   1.22186466
#> Neutrophils_MCPcounter              -0.37665380   -1.7228362   0.26391909
#> Endothelial_cells_MCPcounter        -0.22332011   -1.6438798   1.05057063
#> Fibroblasts_MCPcounter               0.23522176   -0.2693902   0.84959619
#>                                    TCGA-BR-8683 TCGA-BR-8686 TCGA-BR-8687
#> T_cells_MCPcounter                   0.83037884   1.31344757   -1.0357941
#> CD8_T_cells_MCPcounter               0.47377171   1.22349447   -1.3844592
#> Cytotoxic_lymphocytes_MCPcounter     0.79729804   1.62329942   -1.2463579
#> NK_cells_MCPcounter                 -0.32396894   0.63793601   -0.8939450
#> B_lineage_MCPcounter                 0.04363816   1.21471976   -1.0345544
#> Monocytic_lineage_MCPcounter         0.69318719   1.01362692   -0.6473146
#> Myeloid_dendritic_cells_MCPcounter   0.93086147   0.61228229   -0.4818159
#> Neutrophils_MCPcounter               1.01869924  -0.09473877   -0.4519933
#> Endothelial_cells_MCPcounter         2.09125158   1.48040191   -0.3292561
#> Fibroblasts_MCPcounter               1.14665037   0.37777289    0.2721501
#>                                    TCGA-BR-8690 TCGA-BR-A44T TCGA-BR-A44U
#> T_cells_MCPcounter                   0.38245568    1.8414324   -1.7726041
#> CD8_T_cells_MCPcounter               0.63825118    0.5907348   -1.6163353
#> Cytotoxic_lymphocytes_MCPcounter     1.03589512    0.9515275   -1.3060185
#> NK_cells_MCPcounter                  0.59074743    1.1192489   -1.2370764
#> B_lineage_MCPcounter                -0.23786318    2.9063007   -1.4834255
#> Monocytic_lineage_MCPcounter        -0.05420241   -0.3692508   -1.8525733
#> Myeloid_dendritic_cells_MCPcounter   0.14527860    1.3536345   -1.6816981
#> Neutrophils_MCPcounter              -1.38000860   -1.3669333   -1.9398489
#> Endothelial_cells_MCPcounter        -1.66561556   -0.6805785   -1.1582936
#> Fibroblasts_MCPcounter              -0.11975233    0.2813028   -0.3967287
#>                                    TCGA-BR-A4CS TCGA-BR-A4IV TCGA-BR-A4J4
#> T_cells_MCPcounter                 -1.548306024 -0.634567600   -0.4570279
#> CD8_T_cells_MCPcounter             -1.409236370 -0.934854613   -0.3446681
#> Cytotoxic_lymphocytes_MCPcounter   -1.313066442 -0.258363461   -0.1636457
#> NK_cells_MCPcounter                -1.389377764 -0.084082168    0.6285456
#> B_lineage_MCPcounter               -1.273733023 -0.004010939   -0.3722606
#> Monocytic_lineage_MCPcounter       -0.773851519 -0.421636552    0.4965392
#> Myeloid_dendritic_cells_MCPcounter -0.609718025  0.090606581    0.9087483
#> Neutrophils_MCPcounter              0.003061856 -0.390807435    0.1840728
#> Endothelial_cells_MCPcounter       -0.378974173  1.382406804   -0.5586493
#> Fibroblasts_MCPcounter             -0.451177587  2.439318761    0.1532788
#>                                    TCGA-BR-A4J5 TCGA-BR-A4J6 TCGA-BR-A4J7
#> T_cells_MCPcounter                  -0.39357850  -0.92907712    0.8248828
#> CD8_T_cells_MCPcounter              -0.57509890  -0.82090892    1.6833831
#> Cytotoxic_lymphocytes_MCPcounter    -0.29887751  -0.34743040    1.3981892
#> NK_cells_MCPcounter                 -0.40778824  -1.11748701    1.1889249
#> B_lineage_MCPcounter                 0.04196898  -0.09267531   -0.2161863
#> Monocytic_lineage_MCPcounter         0.23661535  -1.53799507    0.3044178
#> Myeloid_dendritic_cells_MCPcounter   0.48658183  -0.29017873    0.1512023
#> Neutrophils_MCPcounter              -0.52741140  -1.19702205   -0.9395282
#> Endothelial_cells_MCPcounter         0.78414360  -1.05523044    0.4278596
#> Fibroblasts_MCPcounter               1.60856965  -0.03484396    1.2277623
#>                                    TCGA-BR-A4J8 TCGA-BR-A4J9 TCGA-BR-A4PF
#> T_cells_MCPcounter                   -1.1375969   -0.1436717    0.4223669
#> CD8_T_cells_MCPcounter               -0.9751123   -0.7102338    1.1662639
#> Cytotoxic_lymphocytes_MCPcounter     -0.9693338   -0.8648378    0.2827860
#> NK_cells_MCPcounter                  -0.8649120   -0.8912893    0.5574982
#> B_lineage_MCPcounter                 -0.3345676    1.0259095   -0.8711627
#> Monocytic_lineage_MCPcounter         -1.1517813   -1.1060320    0.9997170
#> Myeloid_dendritic_cells_MCPcounter   -1.0524321   -0.7014591   -0.3604579
#> Neutrophils_MCPcounter               -1.1426741   -1.0696085    2.7469645
#> Endothelial_cells_MCPcounter          0.2750863    0.4352139   -1.2657274
#> Fibroblasts_MCPcounter                0.1944077    1.0970571   -0.9303587
#>                                    TCGA-BR-A4QL TCGA-CD-5798 TCGA-CD-5799
#> T_cells_MCPcounter                   -1.2781704   -0.8967933   -1.7089532
#> CD8_T_cells_MCPcounter               -0.1863133   -0.2575776   -1.3174431
#> Cytotoxic_lymphocytes_MCPcounter     -0.6176446   -0.7235021   -1.1807508
#> NK_cells_MCPcounter                  -1.0316031   -0.6309516   -1.4719710
#> B_lineage_MCPcounter                 -1.2149526   -0.9864305   -1.0762495
#> Monocytic_lineage_MCPcounter         -1.8450006    0.1038514   -1.2159068
#> Myeloid_dendritic_cells_MCPcounter   -1.7377493    1.3519377   -1.5884880
#> Neutrophils_MCPcounter               -2.3999355   -0.6029836   -0.7334543
#> Endothelial_cells_MCPcounter         -2.6309560    0.4001000   -0.5175613
#> Fibroblasts_MCPcounter               -1.4163488    1.8363322    0.1524457
#>                                    TCGA-CD-5800 TCGA-CD-5801 TCGA-CD-5803
#> T_cells_MCPcounter                   -0.6671678    0.9087164    2.2560840
#> CD8_T_cells_MCPcounter                0.8253353    1.9304772    1.8936460
#> Cytotoxic_lymphocytes_MCPcounter     -0.8448586    1.5077412    1.7502628
#> NK_cells_MCPcounter                  -0.8372590    0.4399147    0.1677833
#> B_lineage_MCPcounter                  0.1098968   -0.3194557    1.4345963
#> Monocytic_lineage_MCPcounter         -0.7443863    0.1707128    1.3147720
#> Myeloid_dendritic_cells_MCPcounter   -0.1606040   -0.5827033    2.2718017
#> Neutrophils_MCPcounter                0.6284791   -1.4945388   -0.7644462
#> Endothelial_cells_MCPcounter         -0.3796929   -1.8439086    1.9531777
#> Fibroblasts_MCPcounter               -0.1366490   -0.9101390    1.8070855
#>                                    TCGA-CD-5804 TCGA-CD-5813 TCGA-CD-8524
#> T_cells_MCPcounter                   -0.5928931    0.8877552 -0.132284583
#> CD8_T_cells_MCPcounter               -0.2968601    0.9710510  0.405733136
#> Cytotoxic_lymphocytes_MCPcounter     -0.3567097    0.5682984  0.226391310
#> NK_cells_MCPcounter                  -0.8688048   -0.5050646 -0.472147089
#> B_lineage_MCPcounter                 -0.5999154    1.2464694  0.175393584
#> Monocytic_lineage_MCPcounter          0.1969164    1.0949361  0.003249202
#> Myeloid_dendritic_cells_MCPcounter   -0.8648727    0.8872338 -0.827869592
#> Neutrophils_MCPcounter                0.7345711   -0.4640988 -0.213862405
#> Endothelial_cells_MCPcounter          0.3963979    2.0341380  0.734194949
#> Fibroblasts_MCPcounter                0.6129305    1.0999535  0.519480859
#>                                    TCGA-CD-8525 TCGA-CD-8526 TCGA-CD-8527
#> T_cells_MCPcounter                  -0.21889814   0.60322085   -0.9265000
#> CD8_T_cells_MCPcounter               0.83809851   0.37289348    0.2446020
#> Cytotoxic_lymphocytes_MCPcounter     0.65778917   0.03633019   -0.3719317
#> NK_cells_MCPcounter                  0.57170325  -0.17288561   -0.8222250
#> B_lineage_MCPcounter                 0.93648002   0.91865455   -0.8500952
#> Monocytic_lineage_MCPcounter         0.05896432   0.04317719   -0.0612910
#> Myeloid_dendritic_cells_MCPcounter  -0.24393702   1.03597215   -1.2346565
#> Neutrophils_MCPcounter              -0.71467072  -0.01862833   -0.8804843
#> Endothelial_cells_MCPcounter        -0.39444327   1.01963147   -1.1924713
#> Fibroblasts_MCPcounter              -0.26175095  -0.14083031   -0.4565134
#>                                    TCGA-CD-8528 TCGA-CD-8529 TCGA-CD-8530
#> T_cells_MCPcounter                   -1.3255287   1.03416095  -0.32259457
#> CD8_T_cells_MCPcounter               -1.2224285   0.41557369  -0.27114000
#> Cytotoxic_lymphocytes_MCPcounter     -0.5184514   0.79851299  -0.30283660
#> NK_cells_MCPcounter                  -0.6496701   0.39591585   0.18528243
#> B_lineage_MCPcounter                 -1.1131756   1.74797191   0.10454185
#> Monocytic_lineage_MCPcounter          0.4946867   1.77336728  -0.07217013
#> Myeloid_dendritic_cells_MCPcounter   -0.5772879   0.58330249   0.34280214
#> Neutrophils_MCPcounter                0.6581040   0.01352478  -0.13817322
#> Endothelial_cells_MCPcounter         -0.7469464   1.58073542   2.23895369
#> Fibroblasts_MCPcounter                0.3380747   0.93413325   1.16196477
#>                                    TCGA-CD-8531 TCGA-CD-8532 TCGA-CD-8533
#> T_cells_MCPcounter                   -0.9210629   1.58890930 -0.393692029
#> CD8_T_cells_MCPcounter               -1.0935236   0.58597242 -0.923871350
#> Cytotoxic_lymphocytes_MCPcounter     -1.0249925   0.05323882 -1.200807128
#> NK_cells_MCPcounter                  -1.3402637  -0.07390049 -0.228960078
#> B_lineage_MCPcounter                  2.8564871   0.81026247 -0.155022039
#> Monocytic_lineage_MCPcounter          0.4293517   0.75348010 -0.709170065
#> Myeloid_dendritic_cells_MCPcounter   -0.5861025   1.90539537 -0.871340805
#> Neutrophils_MCPcounter               -1.2716515   0.28331288  0.330128122
#> Endothelial_cells_MCPcounter         -0.5207691   0.15074292  0.060860431
#> Fibroblasts_MCPcounter                0.1130669  -0.02577957  0.006148538
#>                                    TCGA-CD-8534 TCGA-CD-8535 TCGA-CD-A486
#> T_cells_MCPcounter                  -0.29519283   -1.1294057  -1.06277404
#> CD8_T_cells_MCPcounter              -0.68195364   -1.7039231  -0.70171100
#> Cytotoxic_lymphocytes_MCPcounter    -0.54155734   -1.0082262  -1.13893592
#> NK_cells_MCPcounter                  0.13658705   -0.8263926  -1.31728785
#> B_lineage_MCPcounter                 0.03114773   -0.7460883  -0.44619428
#> Monocytic_lineage_MCPcounter        -0.72564410   -0.9774655   0.09081195
#> Myeloid_dendritic_cells_MCPcounter  -0.11725097   -1.6605148  -0.81441743
#> Neutrophils_MCPcounter              -0.43731710   -0.7785428   1.11167758
#> Endothelial_cells_MCPcounter         0.07075308   -0.4861317  -0.78153339
#> Fibroblasts_MCPcounter              -0.59484467   -0.1744532   0.02913533
#>                                    TCGA-CD-A487 TCGA-CD-A489 TCGA-CD-A48A
#> T_cells_MCPcounter                   -1.2671832   0.05012195   -1.2060055
#> CD8_T_cells_MCPcounter               -0.9650339   0.28156601   -0.8635372
#> Cytotoxic_lymphocytes_MCPcounter     -1.2293019  -0.17660776   -0.8489806
#> NK_cells_MCPcounter                  -1.1153516  -0.68316658   -0.2876138
#> B_lineage_MCPcounter                 -1.2794083   0.74520910   -0.6703567
#> Monocytic_lineage_MCPcounter         -1.4919744  -0.28152427   -0.5862182
#> Myeloid_dendritic_cells_MCPcounter   -0.6222663  -0.02773554   -0.9323388
#> Neutrophils_MCPcounter                0.9291713  -0.72444995   -1.0677571
#> Endothelial_cells_MCPcounter         -0.4568860   0.10611835   -1.6707145
#> Fibroblasts_MCPcounter                0.6698512   0.94808723   -0.7880924
#>                                    TCGA-CD-A48C TCGA-CD-A4MG TCGA-CD-A4MH
#> T_cells_MCPcounter                   -1.0811846   -1.1473085   -0.6156125
#> CD8_T_cells_MCPcounter                0.2477002   -1.2245552    0.1249126
#> Cytotoxic_lymphocytes_MCPcounter     -0.8912183   -0.5632224   -1.1275098
#> NK_cells_MCPcounter                  -1.0652025   -0.9783657   -1.0450227
#> B_lineage_MCPcounter                 -0.5684946    0.2971801   -0.2092329
#> Monocytic_lineage_MCPcounter         -0.3516482   -0.2892926   -0.7600736
#> Myeloid_dendritic_cells_MCPcounter   -0.7669022    0.3895946   -0.1009900
#> Neutrophils_MCPcounter               -1.4494703   -1.5291448    0.3377263
#> Endothelial_cells_MCPcounter         -1.0473869   -0.8092570   -1.4923966
#> Fibroblasts_MCPcounter                0.1090483    0.1208688   -0.8186821
#>                                    TCGA-CG-4301 TCGA-CG-4305 TCGA-CG-4306
#> T_cells_MCPcounter                   0.02113383    0.9608188    0.1306519
#> CD8_T_cells_MCPcounter              -1.08563605    0.6182981   -0.3786750
#> Cytotoxic_lymphocytes_MCPcounter    -0.56115775    1.4839568    0.3411368
#> NK_cells_MCPcounter                 -0.79139851    1.7222682   -0.3091255
#> B_lineage_MCPcounter                -0.48996267    0.1160292    0.1854086
#> Monocytic_lineage_MCPcounter         1.33143766    1.5015159    0.2786719
#> Myeloid_dendritic_cells_MCPcounter   0.48637705    0.5021430    0.7138691
#> Neutrophils_MCPcounter               0.87739085    0.4933233    1.1867360
#> Endothelial_cells_MCPcounter         0.74389020   -0.3530292   -0.3534321
#> Fibroblasts_MCPcounter               0.25199754    0.4167505    0.3411352
#>                                    TCGA-CG-4436 TCGA-CG-4437 TCGA-CG-4438
#> T_cells_MCPcounter                  -0.70093203   0.01169572   0.78490731
#> CD8_T_cells_MCPcounter               0.02392692   0.29973769   0.09231216
#> Cytotoxic_lymphocytes_MCPcounter    -1.29787825   0.47740622   0.98558820
#> NK_cells_MCPcounter                 -1.28550783  -0.28626374   0.65420673
#> B_lineage_MCPcounter                -1.09804407   0.39077846  -0.72205545
#> Monocytic_lineage_MCPcounter        -0.79372215   0.64579110   1.34568246
#> Myeloid_dendritic_cells_MCPcounter  -1.09185424  -0.69573915  -0.34991923
#> Neutrophils_MCPcounter               1.55433005   0.39879834  -0.22723828
#> Endothelial_cells_MCPcounter        -1.32197759  -0.59101042  -1.03620042
#> Fibroblasts_MCPcounter              -1.57074389   0.40368495  -0.40409934
#>                                    TCGA-CG-4440 TCGA-CG-4441 TCGA-CG-4443
#> T_cells_MCPcounter                  -1.34255832   -0.7282413  -0.95104790
#> CD8_T_cells_MCPcounter              -0.03030937   -0.7082745  -1.66891861
#> Cytotoxic_lymphocytes_MCPcounter    -0.54311442    0.2217449  -1.71219891
#> NK_cells_MCPcounter                 -0.04569009    0.1197627  -0.96332343
#> B_lineage_MCPcounter                -1.21592192   -0.5466225  -1.89109375
#> Monocytic_lineage_MCPcounter         0.58580816    0.6540032  -0.99565183
#> Myeloid_dendritic_cells_MCPcounter  -1.04364964   -0.3857866  -1.82241405
#> Neutrophils_MCPcounter               0.91987893    1.6761055  -0.59555033
#> Endothelial_cells_MCPcounter        -0.78387557    0.1458936   0.03265709
#> Fibroblasts_MCPcounter              -0.46351702    0.6773138  -1.50602318
#>                                    TCGA-CG-4444 TCGA-CG-4460 TCGA-CG-4465
#> T_cells_MCPcounter                   0.02267435 -0.494161516   -0.4754335
#> CD8_T_cells_MCPcounter              -0.73537784  0.676753998   -0.9731606
#> Cytotoxic_lymphocytes_MCPcounter    -0.01738737 -0.210717888    0.1253622
#> NK_cells_MCPcounter                  0.41699542 -0.480214480    0.2549814
#> B_lineage_MCPcounter                 0.53637772 -0.122126342   -0.8179853
#> Monocytic_lineage_MCPcounter         0.66342256 -0.009916442    1.4606410
#> Myeloid_dendritic_cells_MCPcounter   0.54415303  0.417009800    1.2605056
#> Neutrophils_MCPcounter               0.77374503 -0.199880828    1.1082308
#> Endothelial_cells_MCPcounter         0.18088834 -0.208256107   -0.3823893
#> Fibroblasts_MCPcounter              -0.19471277 -0.567918328    0.2373842
#>                                    TCGA-CG-4466 TCGA-CG-4469 TCGA-CG-4475
#> T_cells_MCPcounter                   -1.1266591  -0.59321203  -0.75653746
#> CD8_T_cells_MCPcounter               -1.3438740  -1.39759678  -1.39599588
#> Cytotoxic_lymphocytes_MCPcounter     -0.5967996  -0.98355415  -0.86462739
#> NK_cells_MCPcounter                  -0.7015131   0.05098353  -0.97310808
#> B_lineage_MCPcounter                 -1.2559781  -1.28036508  -0.99861784
#> Monocytic_lineage_MCPcounter         -1.1636697   0.85311288   0.46065572
#> Myeloid_dendritic_cells_MCPcounter   -0.8047296   0.10050193   0.06034939
#> Neutrophils_MCPcounter                0.2730024   1.12923879   2.05399575
#> Endothelial_cells_MCPcounter         -1.1224329  -0.37365389   1.63168832
#> Fibroblasts_MCPcounter               -1.7709144  -1.06433126   1.01737762
#>                                    TCGA-CG-4477  TCGA-CG-5717 TCGA-CG-5718
#> T_cells_MCPcounter                    1.9134107  1.0296742345   -0.3156867
#> CD8_T_cells_MCPcounter                0.2672699  0.1102594356   -0.6164246
#> Cytotoxic_lymphocytes_MCPcounter      1.0495646  0.0008638508   -0.3989879
#> NK_cells_MCPcounter                   1.2112368 -0.1756683589   -0.7961282
#> B_lineage_MCPcounter                  0.5832167  0.7174889391   -0.3572057
#> Monocytic_lineage_MCPcounter          1.1611279 -0.0443336236   -0.1364578
#> Myeloid_dendritic_cells_MCPcounter    1.5865368  2.1015053808   -0.7188617
#> Neutrophils_MCPcounter                1.8407517  1.4252653510    0.9274007
#> Endothelial_cells_MCPcounter          0.2843974  0.5697376864   -0.2396859
#> Fibroblasts_MCPcounter               -0.2908777 -0.6878131654   -0.5903437
#>                                    TCGA-CG-5719 TCGA-CG-5720 TCGA-CG-5721
#> T_cells_MCPcounter                  -0.61978370  -0.17213426   1.48194444
#> CD8_T_cells_MCPcounter              -0.89534232  -0.21837673   1.41571024
#> Cytotoxic_lymphocytes_MCPcounter    -0.55047310  -0.39567656   2.35261925
#> NK_cells_MCPcounter                 -0.51925531   0.33330296   0.97016797
#> B_lineage_MCPcounter                -0.57392795   0.10939953   0.29838501
#> Monocytic_lineage_MCPcounter         0.09390342   0.19840320   1.36323766
#> Myeloid_dendritic_cells_MCPcounter  -0.26284293  -0.24608997   0.67092966
#> Neutrophils_MCPcounter               0.84541094   1.18870801  -0.68229420
#> Endothelial_cells_MCPcounter         1.29238826   0.46564461  -0.77182866
#> Fibroblasts_MCPcounter               1.45328708   0.07872336  -0.06572324
#>                                    TCGA-CG-5722 TCGA-CG-5723 TCGA-CG-5724
#> T_cells_MCPcounter                    1.7456196   0.09157095  -1.04018173
#> CD8_T_cells_MCPcounter                1.4702940   0.95121575  -0.50354295
#> Cytotoxic_lymphocytes_MCPcounter      1.8968437   0.72301826  -0.67128346
#> NK_cells_MCPcounter                   3.1826433   0.33778782  -0.67877756
#> B_lineage_MCPcounter                  1.4280026  -0.12145558  -0.79954503
#> Monocytic_lineage_MCPcounter          0.3883372   0.40600925  -0.21863607
#> Myeloid_dendritic_cells_MCPcounter    2.1143438   1.27349799  -1.11777649
#> Neutrophils_MCPcounter               -0.8547646  -0.05339823   0.02331936
#> Endothelial_cells_MCPcounter          0.4372670  -1.04335328  -1.03976896
#> Fibroblasts_MCPcounter               -0.4386326  -0.31386030   0.30383956
#>                                    TCGA-CG-5725 TCGA-CG-5726 TCGA-CG-5732
#> T_cells_MCPcounter                   -1.1048270   -1.5516022    0.3051658
#> CD8_T_cells_MCPcounter                3.8962295   -1.1896194    0.4752452
#> Cytotoxic_lymphocytes_MCPcounter     -0.9711558   -1.6004627   -0.2357061
#> NK_cells_MCPcounter                  -1.2557731   -1.5029749   -0.2589825
#> B_lineage_MCPcounter                 -0.7532144   -0.9978322    1.3359475
#> Monocytic_lineage_MCPcounter         -1.5708770   -1.3287902   -0.2134644
#> Myeloid_dendritic_cells_MCPcounter   -1.2548625   -0.4792397    1.0812901
#> Neutrophils_MCPcounter               -1.2644746   -1.1658649    1.7970584
#> Endothelial_cells_MCPcounter         -0.5204554   -1.5981959   -0.7843264
#> Fibroblasts_MCPcounter               -0.9681474   -1.4142250   -1.3063817
#>                                    TCGA-CG-5734 TCGA-D7-5577 TCGA-D7-5578
#> T_cells_MCPcounter                   0.20079915    0.6841024  0.036543293
#> CD8_T_cells_MCPcounter              -0.39778328    1.3689568 -0.251088703
#> Cytotoxic_lymphocytes_MCPcounter     0.08686887    1.2978186  0.186716577
#> NK_cells_MCPcounter                  1.09076193    0.1018040 -0.361007372
#> B_lineage_MCPcounter                 0.62663384   -0.7021831  0.772383915
#> Monocytic_lineage_MCPcounter         0.11393224    1.2838623  0.129356933
#> Myeloid_dendritic_cells_MCPcounter   0.62975718   -0.3702751  0.190346480
#> Neutrophils_MCPcounter               0.29731035   -0.2742975  0.005463758
#> Endothelial_cells_MCPcounter        -0.59801587   -1.1228109  0.612431225
#> Fibroblasts_MCPcounter              -1.06452074   -0.4854787  0.328620424
#>                                    TCGA-D7-6519 TCGA-D7-6520 TCGA-D7-6521
#> T_cells_MCPcounter                  -0.96897184   -0.9827702   0.22606274
#> CD8_T_cells_MCPcounter               0.02746441   -1.1365373   0.65112068
#> Cytotoxic_lymphocytes_MCPcounter     0.09709383   -1.1766920  -0.01273718
#> NK_cells_MCPcounter                  0.07004205   -1.1170652  -1.05396592
#> B_lineage_MCPcounter                -1.43941954   -0.7635208   1.17045121
#> Monocytic_lineage_MCPcounter        -0.74268263    0.3331075   0.74522706
#> Myeloid_dendritic_cells_MCPcounter  -0.42415793   -0.7628614  -0.24547310
#> Neutrophils_MCPcounter              -0.56447079   -0.7191119  -0.31576558
#> Endothelial_cells_MCPcounter        -1.11383279   -0.2931479   0.55960593
#> Fibroblasts_MCPcounter              -0.01634743    0.8123838   0.63611712
#>                                    TCGA-D7-6522 TCGA-D7-6524 TCGA-D7-6525
#> T_cells_MCPcounter                   1.54354576   -0.4146482  -0.83545217
#> CD8_T_cells_MCPcounter               1.28202120   -0.9205450  -0.72729028
#> Cytotoxic_lymphocytes_MCPcounter     0.75423368   -0.6490324  -0.23253786
#> NK_cells_MCPcounter                  0.03790814   -0.9326409  -0.06707051
#> B_lineage_MCPcounter                 2.21685780   -0.3942216  -0.79297295
#> Monocytic_lineage_MCPcounter         1.55933263    1.0210014   0.58929969
#> Myeloid_dendritic_cells_MCPcounter   2.89579614    1.6849929  -0.98887510
#> Neutrophils_MCPcounter               0.07767919    2.8554520   0.95813725
#> Endothelial_cells_MCPcounter         1.17372500    2.0107420  -0.63762658
#> Fibroblasts_MCPcounter               1.38871558    1.7164448   0.98971602
#>                                    TCGA-D7-6526 TCGA-D7-6527 TCGA-D7-6528
#> T_cells_MCPcounter                  -1.28141309  -1.06764487  -1.10614683
#> CD8_T_cells_MCPcounter              -1.01495324  -0.97107075  -1.46016466
#> Cytotoxic_lymphocytes_MCPcounter    -1.44777719  -0.99776522  -1.02960615
#> NK_cells_MCPcounter                 -1.23915569  -1.18643777  -0.60675254
#> B_lineage_MCPcounter                -0.65401111  -0.31798069  -1.15593612
#> Monocytic_lineage_MCPcounter        -0.65264032  -0.78039115  -1.56111811
#> Myeloid_dendritic_cells_MCPcounter  -0.07324885  -0.08339207  -2.03015920
#> Neutrophils_MCPcounter               0.76700988   0.28818190   0.05358874
#> Endothelial_cells_MCPcounter        -0.83934020  -0.18607695  -0.96908086
#> Fibroblasts_MCPcounter              -0.47546675  -0.22827114  -0.39742560
#>                                    TCGA-D7-6815 TCGA-D7-6818 TCGA-D7-6822
#> T_cells_MCPcounter                   -0.8138597   -1.1648308  -0.33290192
#> CD8_T_cells_MCPcounter                0.1537651   -1.3550610   2.94813072
#> Cytotoxic_lymphocytes_MCPcounter     -0.7977140   -1.0229644   0.35368666
#> NK_cells_MCPcounter                  -1.0597445   -1.0840240  -0.42156564
#> B_lineage_MCPcounter                 -0.1329639    0.1633415  -1.06137916
#> Monocytic_lineage_MCPcounter         -0.1801805    1.2259506  -0.01041736
#> Myeloid_dendritic_cells_MCPcounter   -0.4491749    0.6178930   0.38405714
#> Neutrophils_MCPcounter                0.5216954    1.7217438   2.15795099
#> Endothelial_cells_MCPcounter         -0.6537292    0.9654743  -0.11491052
#> Fibroblasts_MCPcounter                0.2005152    1.3239512  -0.04978778
#>                                    TCGA-D7-8570 TCGA-D7-8572 TCGA-D7-8573
#> T_cells_MCPcounter                    1.0968553    0.5973789   -0.9901327
#> CD8_T_cells_MCPcounter                0.3302214    1.1836548   -1.0191593
#> Cytotoxic_lymphocytes_MCPcounter      1.4027761    0.8506105   -1.1504053
#> NK_cells_MCPcounter                   1.9759369    0.9679806   -0.7121218
#> B_lineage_MCPcounter                  0.4164165    0.3545102   -0.4020856
#> Monocytic_lineage_MCPcounter          1.8178286    1.1902807   -1.0369814
#> Myeloid_dendritic_cells_MCPcounter    1.6683887    0.2985727   -1.0091941
#> Neutrophils_MCPcounter               -1.2269478    2.3422392   -0.2013310
#> Endothelial_cells_MCPcounter         -0.3280501    0.9644417   -1.2936886
#> Fibroblasts_MCPcounter                0.3475059    1.3648887   -1.1140924
#>                                    TCGA-D7-8574 TCGA-D7-8575 TCGA-D7-8576
#> T_cells_MCPcounter                    2.0998664    0.1805876 -1.038592920
#> CD8_T_cells_MCPcounter                0.9293164   -0.4490820 -0.870138661
#> Cytotoxic_lymphocytes_MCPcounter      0.9278474   -0.3377344 -1.276842166
#> NK_cells_MCPcounter                   0.4218645    0.1797262 -1.039877015
#> B_lineage_MCPcounter                  2.8296070    0.3987395 -0.213097504
#> Monocytic_lineage_MCPcounter          0.8694488    1.2798994  0.840681290
#> Myeloid_dendritic_cells_MCPcounter    0.6401735    1.2462616  0.107553063
#> Neutrophils_MCPcounter               -0.4054409    1.3246589  0.545523995
#> Endothelial_cells_MCPcounter          0.4783036   -0.4782630  0.641743603
#> Fibroblasts_MCPcounter                0.9663961   -0.6441586 -0.009732151
#>                                    TCGA-D7-8578 TCGA-D7-8579 TCGA-D7-A4YU
#> T_cells_MCPcounter                  -0.39927882   -0.1655609    1.8912170
#> CD8_T_cells_MCPcounter               0.34705471   -0.5496914    0.6336385
#> Cytotoxic_lymphocytes_MCPcounter    -0.52663666   -0.4366417    1.6678451
#> NK_cells_MCPcounter                 -0.51722115   -0.2245465    0.7724640
#> B_lineage_MCPcounter                -0.01000906    0.5954404    1.9780373
#> Monocytic_lineage_MCPcounter        -0.05529333    0.5509788    1.5142107
#> Myeloid_dendritic_cells_MCPcounter  -0.70774542    0.8667427    0.6141937
#> Neutrophils_MCPcounter               0.13225540    1.8268006    0.1570001
#> Endothelial_cells_MCPcounter         0.10239778    1.0507243    0.6848253
#> Fibroblasts_MCPcounter               1.67699202    1.0446479    0.2230368
#>                                    TCGA-D7-A4YX TCGA-D7-A4Z0 TCGA-D7-A6EV
#> T_cells_MCPcounter                    0.6003394   0.07021377   -0.9024603
#> CD8_T_cells_MCPcounter                0.6432693  -0.76961619   -0.6775560
#> Cytotoxic_lymphocytes_MCPcounter      1.0954318  -0.40945995   -1.2556734
#> NK_cells_MCPcounter                   0.6332635  -0.52421541   -1.3017176
#> B_lineage_MCPcounter                 -0.8833186   1.10911441   -0.5082253
#> Monocytic_lineage_MCPcounter          0.8699105   0.04547389   -0.8252725
#> Myeloid_dendritic_cells_MCPcounter   -0.8380530   1.14047538   -0.3686719
#> Neutrophils_MCPcounter               -0.2513515   0.59873753   -0.1599892
#> Endothelial_cells_MCPcounter         -0.9435058   0.96384102   -1.0189737
#> Fibroblasts_MCPcounter               -0.8724369   0.20305189   -0.5071298
#>                                    TCGA-D7-A6EX TCGA-D7-A6EY TCGA-D7-A6EZ
#> T_cells_MCPcounter                  -0.39810592    1.1867137   0.62476722
#> CD8_T_cells_MCPcounter              -0.93173658    1.0624265   0.69462706
#> Cytotoxic_lymphocytes_MCPcounter    -0.84699148    0.8134216   1.20965581
#> NK_cells_MCPcounter                 -0.96488208    0.8825964   0.09013565
#> B_lineage_MCPcounter                -0.02263931    2.0719641  -0.85388439
#> Monocytic_lineage_MCPcounter        -1.32834094    0.7745380   1.05619683
#> Myeloid_dendritic_cells_MCPcounter  -0.76904300    1.6020135  -0.58847016
#> Neutrophils_MCPcounter              -1.38127414    1.3304490  -1.15578631
#> Endothelial_cells_MCPcounter        -1.60580340   -0.2243806  -0.24308957
#> Fibroblasts_MCPcounter              -0.16322112    0.1986658   0.12627737
#>                                    TCGA-D7-A6F0 TCGA-D7-A6F2 TCGA-D7-A747
#> T_cells_MCPcounter                   0.75765295   0.18597768  0.604887512
#> CD8_T_cells_MCPcounter               1.35576236   0.26784101  0.042343918
#> Cytotoxic_lymphocytes_MCPcounter     0.32981484  -0.16451923  0.009339174
#> NK_cells_MCPcounter                  0.43169638  -0.55155816  0.101261576
#> B_lineage_MCPcounter                -0.03866028  -0.44435577  1.237585772
#> Monocytic_lineage_MCPcounter         0.70196603   0.05064485 -0.208920231
#> Myeloid_dendritic_cells_MCPcounter   0.52342386   0.08080060  0.131531173
#> Neutrophils_MCPcounter              -0.27873996  -0.24477222 -0.847142825
#> Endothelial_cells_MCPcounter        -0.42715313  -0.52806767 -0.798285579
#> Fibroblasts_MCPcounter              -0.56636272  -0.61016603  1.220197985
#>                                    TCGA-D7-A748 TCGA-D7-A74A TCGA-EQ-8122
#> T_cells_MCPcounter                  -0.56562498   -0.8251468   -1.4101856
#> CD8_T_cells_MCPcounter              -0.67534331   -0.5823989   -0.6427547
#> Cytotoxic_lymphocytes_MCPcounter    -0.63663885   -0.9100246   -1.0733489
#> NK_cells_MCPcounter                 -0.27070114   -1.2383184   -0.8190418
#> B_lineage_MCPcounter                -0.34756686   -0.8436862   -1.1155583
#> Monocytic_lineage_MCPcounter         1.75735069   -1.1590964    0.1066605
#> Myeloid_dendritic_cells_MCPcounter   2.22090502   -0.2625821   -1.1708600
#> Neutrophils_MCPcounter               1.82471855   -1.5856522    0.6697993
#> Endothelial_cells_MCPcounter        -0.04623023   -2.1009817    0.9477580
#> Fibroblasts_MCPcounter               0.88867757   -1.9839327    0.5213346
#>                                    TCGA-F1-6177 TCGA-F1-6874 TCGA-F1-6875
#> T_cells_MCPcounter                   -0.4541380   0.50811697   -0.7288228
#> CD8_T_cells_MCPcounter               -0.9228503  -0.02185808    1.8761707
#> Cytotoxic_lymphocytes_MCPcounter     -0.5780892   1.62238310   -1.3953068
#> NK_cells_MCPcounter                  -0.3383474   2.24886944   -1.3033884
#> B_lineage_MCPcounter                 -1.5262835  -0.95111232   -1.4620664
#> Monocytic_lineage_MCPcounter         -0.6426012   0.97932669   -1.4004970
#> Myeloid_dendritic_cells_MCPcounter   -1.0585504   0.07505051   -2.1571758
#> Neutrophils_MCPcounter               -1.0200340   3.53753965   -0.9133895
#> Endothelial_cells_MCPcounter         -0.6920229   1.59693814   -0.2083479
#> Fibroblasts_MCPcounter               -0.3779553   0.82076712   -0.7642005
#>                                    TCGA-F1-A448 TCGA-F1-A72C TCGA-FP-7735
#> T_cells_MCPcounter                   0.05134008  -0.12965884   -0.4240454
#> CD8_T_cells_MCPcounter              -0.84923872  -0.07193885   -0.5181869
#> Cytotoxic_lymphocytes_MCPcounter     0.28180731   0.21884120   -0.5367173
#> NK_cells_MCPcounter                  0.10569592   0.40340805    0.4212924
#> B_lineage_MCPcounter                 1.07994213   0.39171457    0.2231642
#> Monocytic_lineage_MCPcounter         0.62382523  -0.18170365   -0.6693195
#> Myeloid_dendritic_cells_MCPcounter   1.94124978  -0.20217194   -0.9042991
#> Neutrophils_MCPcounter               1.06046481  -0.60109004   -1.5230696
#> Endothelial_cells_MCPcounter         0.89324915  -0.11391594    0.3981306
#> Fibroblasts_MCPcounter               0.32638406   0.51013860   -0.2468965
#>                                    TCGA-FP-7829 TCGA-FP-7916 TCGA-FP-7998
#> T_cells_MCPcounter                  -1.48487534    1.4351766    2.1228457
#> CD8_T_cells_MCPcounter              -1.82830879    1.0966897    2.6236219
#> Cytotoxic_lymphocytes_MCPcounter    -1.87121704    1.4684714    2.4427354
#> NK_cells_MCPcounter                 -1.43192212    1.2114836    1.6108595
#> B_lineage_MCPcounter                -0.36721494    1.1940060    1.1650121
#> Monocytic_lineage_MCPcounter         0.19213330    1.7497118    2.5097276
#> Myeloid_dendritic_cells_MCPcounter  -0.27857615    0.2287783    0.7039139
#> Neutrophils_MCPcounter               0.04996474   -0.6773431   -0.5695254
#> Endothelial_cells_MCPcounter         0.24294423    0.8179491    1.0450357
#> Fibroblasts_MCPcounter               0.70943617    0.4159725    0.7200308
#>                                    TCGA-FP-8099 TCGA-FP-8210 TCGA-FP-8211
#> T_cells_MCPcounter                   -0.7473233    2.0628388  -0.06111483
#> CD8_T_cells_MCPcounter               -1.4972110    1.1361539  -0.06703657
#> Cytotoxic_lymphocytes_MCPcounter     -0.8863104    0.9715503  -0.39398281
#> NK_cells_MCPcounter                  -0.2525737    1.2224982  -0.88418323
#> B_lineage_MCPcounter                 -0.3385714    2.3741449  -0.46126144
#> Monocytic_lineage_MCPcounter          0.4322177    1.1529236  -0.97136323
#> Myeloid_dendritic_cells_MCPcounter   -0.4392100    2.0786212  -0.93076420
#> Neutrophils_MCPcounter                0.8236051    1.4448777  -0.68729194
#> Endothelial_cells_MCPcounter         -0.6955349    1.8109563  -0.91887339
#> Fibroblasts_MCPcounter                0.3466095    1.5358104  -1.23349771
#>                                    TCGA-FP-8631 TCGA-FP-A4BF TCGA-FP-A8CX
#> T_cells_MCPcounter                   -0.7285998    1.2801217   0.08006747
#> CD8_T_cells_MCPcounter               -1.2491098    0.9880879  -0.52627263
#> Cytotoxic_lymphocytes_MCPcounter     -1.4252275    2.3330325  -0.88017840
#> NK_cells_MCPcounter                  -0.1269108    3.7573305  -0.86652229
#> B_lineage_MCPcounter                 -0.3629054    0.2122430   0.04779480
#> Monocytic_lineage_MCPcounter         -0.4358658    1.1565744  -0.46954227
#> Myeloid_dendritic_cells_MCPcounter   -0.4090448    0.2871034  -0.17542170
#> Neutrophils_MCPcounter                0.4030741   -1.0977204  -0.41130503
#> Endothelial_cells_MCPcounter          0.1152409   -0.6715280  -0.96314633
#> Fibroblasts_MCPcounter                0.4913639    0.7622362  -0.83067490
#>                                    TCGA-FP-A9TM TCGA-HF-7132 TCGA-HF-7133
#> T_cells_MCPcounter                    0.6881702    1.4907649    0.5112719
#> CD8_T_cells_MCPcounter               -0.6739010    0.4741894    0.6010807
#> Cytotoxic_lymphocytes_MCPcounter     -0.4233516    1.1080072    0.3124028
#> NK_cells_MCPcounter                  -0.2228633    1.0176607    0.1508329
#> B_lineage_MCPcounter                  2.6572201    1.4356152    0.3586989
#> Monocytic_lineage_MCPcounter         -0.4741132    1.0090948    1.3959456
#> Myeloid_dendritic_cells_MCPcounter    1.0431611    1.0825442    1.7312792
#> Neutrophils_MCPcounter               -0.1515265    0.7387106    1.1780306
#> Endothelial_cells_MCPcounter         -0.7121093    0.4730183    0.4997663
#> Fibroblasts_MCPcounter               -1.4701923    0.3696295    0.0644028
#>                                    TCGA-HF-7134 TCGA-HF-A5NB TCGA-HJ-7597
#> T_cells_MCPcounter                 -0.748570675   -0.3353602    0.4701723
#> CD8_T_cells_MCPcounter             -0.968370022   -0.3761505    0.7654361
#> Cytotoxic_lymphocytes_MCPcounter   -0.540912615   -0.1917195    1.3447352
#> NK_cells_MCPcounter                -0.002844287   -0.5726316    2.2651543
#> B_lineage_MCPcounter               -0.158621099   -1.4592601    0.5157704
#> Monocytic_lineage_MCPcounter       -0.125509781   -0.7566786   -0.6832066
#> Myeloid_dendritic_cells_MCPcounter -0.656561437   -1.7081019   -1.0758326
#> Neutrophils_MCPcounter              0.807608602   -0.8707341    0.7519758
#> Endothelial_cells_MCPcounter        1.478841817   -1.2541591   -0.7132986
#> Fibroblasts_MCPcounter             -0.554563590   -1.0117658   -0.4576280
#>                                    TCGA-HU-8238 TCGA-HU-8244 TCGA-HU-8249
#> T_cells_MCPcounter                  -0.75395908   -1.3521394   0.41036895
#> CD8_T_cells_MCPcounter              -0.74108437   -1.0896416   0.09353722
#> Cytotoxic_lymphocytes_MCPcounter    -0.98909037   -1.5315892   0.18497355
#> NK_cells_MCPcounter                 -0.04839321   -0.4747373  -0.18967803
#> B_lineage_MCPcounter                 0.09770644   -1.2088777  -0.53534894
#> Monocytic_lineage_MCPcounter        -0.48812599   -1.9305732  -1.08871430
#> Myeloid_dendritic_cells_MCPcounter  -0.09574519   -1.1277485  -1.03159511
#> Neutrophils_MCPcounter              -0.55892388   -1.0967610  -0.62272261
#> Endothelial_cells_MCPcounter        -0.69656967   -1.8511630  -0.41062471
#> Fibroblasts_MCPcounter              -1.00527284   -2.9346374  -0.87478004
#>                                    TCGA-HU-8602 TCGA-HU-8604 TCGA-HU-8608
#> T_cells_MCPcounter                   0.75613632    1.2769055    2.3585496
#> CD8_T_cells_MCPcounter              -0.36792108    0.7135141    1.6360456
#> Cytotoxic_lymphocytes_MCPcounter     0.99057658    2.0299685    2.8166587
#> NK_cells_MCPcounter                  1.19915429    3.4088852    2.2055442
#> B_lineage_MCPcounter                -0.87423589    0.2885302    0.9126950
#> Monocytic_lineage_MCPcounter         1.74569746    0.8771575    2.6187385
#> Myeloid_dendritic_cells_MCPcounter  -0.59545448   -0.7950323    0.1941808
#> Neutrophils_MCPcounter               1.30804853    0.3359939    1.3894505
#> Endothelial_cells_MCPcounter         0.09978377    0.3610968    0.9846486
#> Fibroblasts_MCPcounter               0.18498875    0.4908632   -0.2332137
#>                                    TCGA-HU-8610 TCGA-HU-A4G2 TCGA-HU-A4G3
#> T_cells_MCPcounter                    1.4522425  -0.66180079   -1.6093779
#> CD8_T_cells_MCPcounter               -1.3131652   0.05950909   -1.2074427
#> Cytotoxic_lymphocytes_MCPcounter      0.5676310  -0.08919893   -1.0546113
#> NK_cells_MCPcounter                  -0.1682870   0.44961421   -0.7990202
#> B_lineage_MCPcounter                 -0.3797975  -0.99796461   -0.7236671
#> Monocytic_lineage_MCPcounter          1.3027499  -0.34008440   -1.1041069
#> Myeloid_dendritic_cells_MCPcounter    0.5870276   0.20865585   -0.5915385
#> Neutrophils_MCPcounter                0.3159934  -1.16803577   -0.6551882
#> Endothelial_cells_MCPcounter          1.5837282   1.06328671   -0.5668115
#> Fibroblasts_MCPcounter                0.2705998  -0.90085058   -0.9552537
#>                                    TCGA-HU-A4G8 TCGA-HU-A4G9 TCGA-HU-A4GC
#> T_cells_MCPcounter                   0.30618633    -1.831699   -0.9676661
#> CD8_T_cells_MCPcounter              -0.36181298    -1.692944   -0.3493120
#> Cytotoxic_lymphocytes_MCPcounter     0.33961052    -1.546920   -0.5525097
#> NK_cells_MCPcounter                  0.04469672    -1.476086   -0.4659242
#> B_lineage_MCPcounter                 1.21500382    -1.171196   -0.5477778
#> Monocytic_lineage_MCPcounter         0.04165809    -1.977086   -0.6690931
#> Myeloid_dendritic_cells_MCPcounter   0.67474222    -1.189500   -0.8379038
#> Neutrophils_MCPcounter              -0.73856278    -1.113252   -0.5286090
#> Endothelial_cells_MCPcounter        -0.40293322    -2.028854   -0.5964296
#> Fibroblasts_MCPcounter              -0.65597785    -2.492085   -0.7769645
#>                                    TCGA-HU-A4GD TCGA-HU-A4GF TCGA-HU-A4GH
#> T_cells_MCPcounter                   -0.6343217   -0.1446038  -1.06132317
#> CD8_T_cells_MCPcounter                0.1804857    0.5193891  -0.79603245
#> Cytotoxic_lymphocytes_MCPcounter     -0.4920175   -0.1894998  -0.84682892
#> NK_cells_MCPcounter                   1.6807181   -0.4245306   0.03064024
#> B_lineage_MCPcounter                 -0.5467774   -0.7348713  -0.46724674
#> Monocytic_lineage_MCPcounter         -1.3217796   -0.1246575  -1.10064064
#> Myeloid_dendritic_cells_MCPcounter   -0.6606342   -0.7007439  -0.77768067
#> Neutrophils_MCPcounter               -0.2845816    2.0865081   0.86814531
#> Endothelial_cells_MCPcounter         -1.5216367    0.5498380  -1.58789758
#> Fibroblasts_MCPcounter               -2.1586717   -1.1417784  -2.46692499
#>                                    TCGA-HU-A4GJ TCGA-HU-A4GP TCGA-HU-A4GQ
#> T_cells_MCPcounter                   2.84314874   -1.3936335   -0.6960175
#> CD8_T_cells_MCPcounter               1.27868662   -0.8012554   -0.5893638
#> Cytotoxic_lymphocytes_MCPcounter     1.13185025   -1.3400536   -0.3678990
#> NK_cells_MCPcounter                  0.02027848   -0.7072571   -0.5595463
#> B_lineage_MCPcounter                 3.76137417   -0.1755291   -1.1921747
#> Monocytic_lineage_MCPcounter         0.32433086   -1.2414440   -0.3729079
#> Myeloid_dendritic_cells_MCPcounter   0.84116846   -0.6583102   -0.6478051
#> Neutrophils_MCPcounter              -0.61591184   -1.4926031   -1.0894666
#> Endothelial_cells_MCPcounter         0.06360296    0.1199794   -0.9857721
#> Fibroblasts_MCPcounter              -0.77159023   -1.0368932    0.3878947
#>                                    TCGA-HU-A4GT TCGA-HU-A4GU TCGA-HU-A4GX
#> T_cells_MCPcounter                   -0.9141816   -0.5034956    0.3141285
#> CD8_T_cells_MCPcounter               -0.5173244   -1.1380992    1.1880656
#> Cytotoxic_lymphocytes_MCPcounter     -1.1066362   -0.6649737    0.8747111
#> NK_cells_MCPcounter                  -0.9120861   -1.0417291    0.5496320
#> B_lineage_MCPcounter                 -0.6367667   -1.6254556   -0.0778396
#> Monocytic_lineage_MCPcounter         -1.0214683   -1.1960505   -0.6168687
#> Myeloid_dendritic_cells_MCPcounter   -0.1053461   -1.7876686   -1.5570573
#> Neutrophils_MCPcounter               -0.4907000   -0.6307041   -1.7660847
#> Endothelial_cells_MCPcounter         -1.1689941   -1.3912213   -1.5361169
#> Fibroblasts_MCPcounter               -1.0664019   -0.4139173   -1.3206843
#>                                    TCGA-HU-A4GY TCGA-HU-A4H0 TCGA-HU-A4H2
#> T_cells_MCPcounter                    2.4275517    0.9477997   -0.6572247
#> CD8_T_cells_MCPcounter                0.4800308    0.5426138   -0.5944241
#> Cytotoxic_lymphocytes_MCPcounter      1.0975948    1.0675887   -0.5815858
#> NK_cells_MCPcounter                   0.2865828    0.3535259   -0.6835401
#> B_lineage_MCPcounter                  2.2644614   -0.7097515    0.2293700
#> Monocytic_lineage_MCPcounter          1.7395274    0.9826592   -0.9383022
#> Myeloid_dendritic_cells_MCPcounter    2.0026778   -0.6248765   -0.8560678
#> Neutrophils_MCPcounter               -0.8486565   -1.1796817   -1.3551223
#> Endothelial_cells_MCPcounter          0.6545757   -0.5439526   -1.8846564
#> Fibroblasts_MCPcounter                0.5696756   -1.1567964   -1.1549937
#>                                    TCGA-HU-A4H3 TCGA-HU-A4H4  TCGA-HU-A4H5
#> T_cells_MCPcounter                   -0.4049793   0.82831273 -0.6302370142
#> CD8_T_cells_MCPcounter               -0.3648689   0.40366137 -0.8649489112
#> Cytotoxic_lymphocytes_MCPcounter      0.1486127   1.11412367 -0.2961712904
#> NK_cells_MCPcounter                   0.5357280   0.18756175  0.1127807026
#> B_lineage_MCPcounter                 -0.8952936   0.35107991 -0.1488192879
#> Monocytic_lineage_MCPcounter         -1.1910861   0.36471793 -0.5931728025
#> Myeloid_dendritic_cells_MCPcounter   -1.8221942  -0.22525255 -0.0992595589
#> Neutrophils_MCPcounter               -0.8394884  -0.63834245  0.8264188096
#> Endothelial_cells_MCPcounter         -1.3064783  -0.08476919 -0.0002963556
#> Fibroblasts_MCPcounter               -0.4118089  -0.70922496 -1.0506706189
#>                                    TCGA-HU-A4H6 TCGA-HU-A4H8 TCGA-HU-A4HB
#> T_cells_MCPcounter                  -0.48671541   -1.0568269    1.6840906
#> CD8_T_cells_MCPcounter              -0.69775021   -0.9962819    0.6484671
#> Cytotoxic_lymphocytes_MCPcounter    -0.82040083   -0.5999654    0.4336776
#> NK_cells_MCPcounter                 -0.30245257   -0.5048574    0.5839081
#> B_lineage_MCPcounter                 1.30809748   -0.9142332    2.8182344
#> Monocytic_lineage_MCPcounter        -0.98775701   -1.1604815    0.5236263
#> Myeloid_dendritic_cells_MCPcounter  -0.69122230   -0.6232153    1.2030387
#> Neutrophils_MCPcounter               0.11978191   -1.1899379   -0.2233462
#> Endothelial_cells_MCPcounter        -0.01616177   -1.1125975   -0.3481518
#> Fibroblasts_MCPcounter              -0.88646897   -1.5840700   -0.8765339
#>                                    TCGA-HU-A4HD TCGA-IN-7806 TCGA-IN-7808
#> T_cells_MCPcounter                  -1.11207653  -0.34280623   2.80845014
#> CD8_T_cells_MCPcounter               0.74652023  -1.20527224   1.39357169
#> Cytotoxic_lymphocytes_MCPcounter    -1.29843277  -1.50797533   1.06263226
#> NK_cells_MCPcounter                 -1.19563580  -0.74117300   2.06086071
#> B_lineage_MCPcounter                -0.77459600  -0.20879650   2.98776162
#> Monocytic_lineage_MCPcounter        -0.14982026  -1.06342301  -0.07570153
#> Myeloid_dendritic_cells_MCPcounter  -0.96776989  -0.73248096   1.20038904
#> Neutrophils_MCPcounter               0.05466657  -0.22887894  -1.31349995
#> Endothelial_cells_MCPcounter        -0.43622100   0.05213299  -0.38607694
#> Fibroblasts_MCPcounter              -0.30984120  -0.65923623  -1.23318192
#>                                     TCGA-IN-8462 TCGA-IN-8663 TCGA-IN-A6RI
#> T_cells_MCPcounter                 -0.2348183909   -1.2198737   -1.1024497
#> CD8_T_cells_MCPcounter             -0.9974506249   -1.1851018   -1.1446122
#> Cytotoxic_lymphocytes_MCPcounter   -1.0386407131   -1.1527906   -1.1626423
#> NK_cells_MCPcounter                 0.0008745893   -0.4863508   -1.3255975
#> B_lineage_MCPcounter               -0.4389717293   -1.5737355   -1.0587992
#> Monocytic_lineage_MCPcounter       -1.1318737245   -0.1863068   -1.7319890
#> Myeloid_dendritic_cells_MCPcounter  0.9807729544   -1.3159530   -0.1993816
#> Neutrophils_MCPcounter             -0.3377497501    0.4726686    0.3682701
#> Endothelial_cells_MCPcounter       -0.3675887589   -0.6415687   -1.0363888
#> Fibroblasts_MCPcounter             -0.2559602107   -0.2140173   -1.7414158
#>                                    TCGA-IN-A6RJ TCGA-IN-A6RL TCGA-IN-A6RN
#> T_cells_MCPcounter                   -0.5019799  -0.33702636   -0.9091938
#> CD8_T_cells_MCPcounter               -1.0383447  -0.92237829   -1.2715763
#> Cytotoxic_lymphocytes_MCPcounter     -1.1386257  -0.78112161   -1.3142640
#> NK_cells_MCPcounter                  -0.3273378  -0.92981760   -1.0139760
#> B_lineage_MCPcounter                 -0.7727281   0.08776159   -0.7805431
#> Monocytic_lineage_MCPcounter         -1.9620562  -0.73642195   -2.1251295
#> Myeloid_dendritic_cells_MCPcounter    0.9933210   0.34277291    0.5407614
#> Neutrophils_MCPcounter                0.0677835   0.17137525   -0.4969676
#> Endothelial_cells_MCPcounter         -1.6815205  -0.31451735   -0.9587696
#> Fibroblasts_MCPcounter               -2.1679883  -0.68680894   -0.9502645
#>                                    TCGA-IN-A6RR TCGA-IN-A6RS TCGA-IN-A7NR
#> T_cells_MCPcounter                   -0.9548977  0.146191592   0.38754687
#> CD8_T_cells_MCPcounter               -0.7597316 -0.002770485  -0.09274166
#> Cytotoxic_lymphocytes_MCPcounter     -0.7304528  0.287258641  -0.63951718
#> NK_cells_MCPcounter                  -0.2607966 -1.261195224  -0.69227572
#> B_lineage_MCPcounter                 -0.5417830 -0.657090392   0.54050205
#> Monocytic_lineage_MCPcounter         -1.3862232 -0.618758306  -1.29067585
#> Myeloid_dendritic_cells_MCPcounter   -0.7605720 -0.953006387   0.38040419
#> Neutrophils_MCPcounter               -0.5492989  0.175280260  -0.21735345
#> Endothelial_cells_MCPcounter         -0.4166209  0.040349681  -0.60132091
#> Fibroblasts_MCPcounter               -1.4919253 -0.971510812  -0.74702198
#>                                    TCGA-IN-A7NT TCGA-IN-A7NU TCGA-IN-AB1V
#> T_cells_MCPcounter                 -0.779386482  -0.05363204   0.06899453
#> CD8_T_cells_MCPcounter             -1.691442045  -0.65181254  -0.62307216
#> Cytotoxic_lymphocytes_MCPcounter   -0.013375219  -0.48409294  -0.51956143
#> NK_cells_MCPcounter                -0.330021921  -0.58457497   1.09937683
#> B_lineage_MCPcounter               -1.114499998  -0.61725608   0.24481312
#> Monocytic_lineage_MCPcounter        0.035158077  -0.49209154  -1.85135066
#> Myeloid_dendritic_cells_MCPcounter -0.727762178   0.39445156   0.20906985
#> Neutrophils_MCPcounter              0.007454531   0.18964126   1.37578890
#> Endothelial_cells_MCPcounter       -1.101445350  -0.80638579  -0.86773456
#> Fibroblasts_MCPcounter              0.530893535   0.12932706  -1.51244633
#>                                    TCGA-IN-AB1X TCGA-IP-7968 TCGA-KB-A6F7
#> T_cells_MCPcounter                   2.11943779    0.6147680   0.98408282
#> CD8_T_cells_MCPcounter               1.33487506    0.5409654   1.21842910
#> Cytotoxic_lymphocytes_MCPcounter     1.90413042    1.0947092   1.10066845
#> NK_cells_MCPcounter                  1.68767993    0.3612907   1.18174743
#> B_lineage_MCPcounter                 0.68334825   -0.5114751   0.07718502
#> Monocytic_lineage_MCPcounter         0.46322520    0.1880279   0.56907265
#> Myeloid_dendritic_cells_MCPcounter   0.27438712   -0.4982216   1.10487671
#> Neutrophils_MCPcounter               1.71270631    0.4771820   0.17786197
#> Endothelial_cells_MCPcounter        -0.04458868    1.5890670  -0.60606405
#> Fibroblasts_MCPcounter              -0.76794762    0.7579484  -1.68655799
#>                                    TCGA-KB-A93G TCGA-KB-A93H TCGA-KB-A93J
#> T_cells_MCPcounter                   -0.3109756    -1.016417    0.3710396
#> CD8_T_cells_MCPcounter               -1.3885948     2.399238   -0.2756385
#> Cytotoxic_lymphocytes_MCPcounter     -0.5100531    -1.340131    0.7471930
#> NK_cells_MCPcounter                   0.4589173    -1.490383    2.6420188
#> B_lineage_MCPcounter                 -0.6458285    -1.521649   -0.2763097
#> Monocytic_lineage_MCPcounter          0.6729436    -1.131822    1.2473090
#> Myeloid_dendritic_cells_MCPcounter    0.6137678    -1.987647    0.2338439
#> Neutrophils_MCPcounter                3.7733633    -1.293292    0.6797115
#> Endothelial_cells_MCPcounter          2.3984792    -1.669478   -0.3160977
#> Fibroblasts_MCPcounter                2.6065999    -2.263071   -0.8515199
#>                                    TCGA-MX-A5UG TCGA-MX-A5UJ TCGA-MX-A663
#> T_cells_MCPcounter                    1.8624443   -0.4484992  -0.42798754
#> CD8_T_cells_MCPcounter                0.5608955   -0.5709078  -0.08683090
#> Cytotoxic_lymphocytes_MCPcounter      0.8482950   -0.2682040  -0.25023103
#> NK_cells_MCPcounter                   0.6728872   -0.8874266  -0.01505505
#> B_lineage_MCPcounter                  2.6211689    0.6469748  -0.19472260
#> Monocytic_lineage_MCPcounter          1.2890810   -0.5635708  -0.09013648
#> Myeloid_dendritic_cells_MCPcounter    1.7709396   -0.4960747  -0.51087069
#> Neutrophils_MCPcounter                1.7578771   -0.2243880   0.30909602
#> Endothelial_cells_MCPcounter          2.1316499    0.8578987   1.93643112
#> Fibroblasts_MCPcounter                1.9796845    0.7646741   1.40774562
#>                                    TCGA-MX-A666 TCGA-R5-A7O7 TCGA-R5-A7ZE
#> T_cells_MCPcounter                   0.38018217   -0.2462152   -1.2355186
#> CD8_T_cells_MCPcounter               1.37461700   -1.0183485    0.1157674
#> Cytotoxic_lymphocytes_MCPcounter    -0.73240681   -0.9996277   -1.7361365
#> NK_cells_MCPcounter                 -0.03264818   -0.7068696   -1.1286711
#> B_lineage_MCPcounter                 2.18406750   -0.4349729   -1.5131972
#> Monocytic_lineage_MCPcounter        -0.48127820   -1.6602791   -1.7577920
#> Myeloid_dendritic_cells_MCPcounter  -0.22247974   -1.0578287    0.2098126
#> Neutrophils_MCPcounter               1.44327269   -0.9488383    0.7644849
#> Endothelial_cells_MCPcounter         0.59271403   -0.3018178   -0.5444238
#> Fibroblasts_MCPcounter              -0.63437849   -0.9906582   -1.6913523
#>                                    TCGA-R5-A7ZF TCGA-R5-A7ZI TCGA-R5-A7ZR
#> T_cells_MCPcounter                   -1.0656768   0.35805852   -1.1826910
#> CD8_T_cells_MCPcounter               -1.2658539   1.11497785    0.4223707
#> Cytotoxic_lymphocytes_MCPcounter     -0.9598799   0.73715072   -1.0960189
#> NK_cells_MCPcounter                  -0.4197636   0.21031365   -1.2405352
#> B_lineage_MCPcounter                 -1.1555835   0.14745301   -0.3385025
#> Monocytic_lineage_MCPcounter         -2.3393099  -0.07469863   -0.9347337
#> Myeloid_dendritic_cells_MCPcounter   -1.4928781  -0.24645783   -1.0504828
#> Neutrophils_MCPcounter               -0.9351950  -1.96333194   -0.6950528
#> Endothelial_cells_MCPcounter         -0.8746742  -1.84475015   -1.5163646
#> Fibroblasts_MCPcounter               -1.3549146  -2.02218630   -2.0174471
#>                                    TCGA-R5-A805 TCGA-RD-A7BS TCGA-RD-A7BT
#> T_cells_MCPcounter                  0.007660377   0.58002270   -0.1924783
#> CD8_T_cells_MCPcounter             -0.569005586  -0.41720907   -0.3698278
#> Cytotoxic_lymphocytes_MCPcounter   -0.399398604   0.05516082   -0.7937205
#> NK_cells_MCPcounter                -0.593785077   1.80762600   -1.3991726
#> B_lineage_MCPcounter               -0.548520864   1.28075287   -0.6716172
#> Monocytic_lineage_MCPcounter        0.840349806  -0.53842861   -1.1180342
#> Myeloid_dendritic_cells_MCPcounter -0.672055935   1.15426644   -0.6931298
#> Neutrophils_MCPcounter              2.992710724   0.62722391    0.2934270
#> Endothelial_cells_MCPcounter        0.429383606   0.15906169   -1.1968223
#> Fibroblasts_MCPcounter              0.067068066  -0.16838966   -1.6076389
#>                                    TCGA-RD-A7BW  TCGA-RD-A7C1 TCGA-RD-A8MV
#> T_cells_MCPcounter                    1.3239758  0.3158152779    0.8358551
#> CD8_T_cells_MCPcounter                1.7002711  0.8578124227    0.2584469
#> Cytotoxic_lymphocytes_MCPcounter      1.8338759  0.9941270682    1.0254839
#> NK_cells_MCPcounter                   0.4319312  2.5750850808    1.6469364
#> B_lineage_MCPcounter                 -0.5933689 -0.0535738179    1.2241314
#> Monocytic_lineage_MCPcounter          0.7848779  0.5459316639    0.6300079
#> Myeloid_dendritic_cells_MCPcounter    1.1943860  0.0004231372   -0.1141100
#> Neutrophils_MCPcounter               -0.3161227  0.9349330705    1.6594393
#> Endothelial_cells_MCPcounter          1.0658407 -0.4429114778   -0.3172437
#> Fibroblasts_MCPcounter                1.2775026 -0.8444660003   -0.9941662
#>                                    TCGA-RD-A8MW TCGA-RD-A8N0 TCGA-RD-A8N1
#> T_cells_MCPcounter                    0.2590742   1.58800426   1.89317996
#> CD8_T_cells_MCPcounter                0.8599108   0.89571871   0.77690332
#> Cytotoxic_lymphocytes_MCPcounter     -0.2813855   1.31686191   1.32866368
#> NK_cells_MCPcounter                  -0.5328222   0.97287842   1.88885988
#> B_lineage_MCPcounter                  0.3092608   2.48384853   1.88082274
#> Monocytic_lineage_MCPcounter          0.8722095  -0.03473960   0.76037110
#> Myeloid_dendritic_cells_MCPcounter    0.5356978   0.89516970   1.48030864
#> Neutrophils_MCPcounter               -0.5132824  -0.51773311  -0.03451794
#> Endothelial_cells_MCPcounter         -0.2257081  -0.01012224   0.15326346
#> Fibroblasts_MCPcounter                0.6405262  -0.12334663   0.09294732
#>                                    TCGA-RD-A8N4 TCGA-RD-A8N5 TCGA-RD-A8N6
#> T_cells_MCPcounter                   0.26916175  -1.00570602   -0.7514774
#> CD8_T_cells_MCPcounter               0.20928926  -1.34024562    0.4400325
#> Cytotoxic_lymphocytes_MCPcounter    -0.03920455  -0.72450077   -0.9242042
#> NK_cells_MCPcounter                 -0.39199181  -0.04691583   -0.6104751
#> B_lineage_MCPcounter                 0.44165233  -0.91542049   -1.2308377
#> Monocytic_lineage_MCPcounter        -0.46811214   0.43733576   -0.3139150
#> Myeloid_dendritic_cells_MCPcounter  -0.29895856  -0.56115548    0.2169923
#> Neutrophils_MCPcounter              -0.77498463   2.33449270    0.1993975
#> Endothelial_cells_MCPcounter         0.35768324   0.77902671   -0.2511461
#> Fibroblasts_MCPcounter               0.91649651   1.09300407    0.4271583
#>                                    TCGA-RD-A8N9 TCGA-RD-A8NB TCGA-SW-A7EA
#> T_cells_MCPcounter                   0.28029368   0.75773141   -0.8585477
#> CD8_T_cells_MCPcounter               0.15774657   0.06523145   -0.0310774
#> Cytotoxic_lymphocytes_MCPcounter    -0.03695698   0.89458437    0.4221331
#> NK_cells_MCPcounter                 -0.06095853   2.58120632    0.2452058
#> B_lineage_MCPcounter                 0.53968777   0.93519669   -0.7415811
#> Monocytic_lineage_MCPcounter         0.53054555   0.57094490   -1.0096804
#> Myeloid_dendritic_cells_MCPcounter   1.06253617   0.32914278   -0.4390756
#> Neutrophils_MCPcounter              -0.93630787   1.08918120   -1.1561958
#> Endothelial_cells_MCPcounter         1.64030140   0.74448042   -1.9569691
#> Fibroblasts_MCPcounter               1.31715087   0.70276472    0.2477852
#>                                    TCGA-SW-A7EB TCGA-VQ-A8DT TCGA-VQ-A8DU
#> T_cells_MCPcounter                   0.66449654   -1.1730269  -0.42659718
#> CD8_T_cells_MCPcounter               0.61009765   -0.6825798   0.37211495
#> Cytotoxic_lymphocytes_MCPcounter     0.08030363   -0.9839822  -0.06155575
#> NK_cells_MCPcounter                 -0.28410019   -1.0199636  -0.21494543
#> B_lineage_MCPcounter                 0.94854824   -0.5021647  -0.34335410
#> Monocytic_lineage_MCPcounter        -0.95462700   -1.6825060  -0.62624571
#> Myeloid_dendritic_cells_MCPcounter  -0.24328890   -0.5814271   0.03879639
#> Neutrophils_MCPcounter              -0.01140449   -0.8399314   0.05381766
#> Endothelial_cells_MCPcounter        -0.54056127   -1.5029271  -0.78318650
#> Fibroblasts_MCPcounter              -0.89579964   -1.7028863  -0.11243880
#>                                    TCGA-VQ-A8DV TCGA-VQ-A8DZ TCGA-VQ-A8E0
#> T_cells_MCPcounter                   -1.7576238   -1.1622297  -0.43868205
#> CD8_T_cells_MCPcounter                0.0252193    1.8450835  -1.58675815
#> Cytotoxic_lymphocytes_MCPcounter     -1.2418048   -0.6225566  -1.13891815
#> NK_cells_MCPcounter                  -1.1525129   -0.8612923  -0.57555113
#> B_lineage_MCPcounter                 -1.1064282   -0.7536642  -1.03876603
#> Monocytic_lineage_MCPcounter         -1.6103901   -0.1389235  -0.13573218
#> Myeloid_dendritic_cells_MCPcounter   -1.2866723   -0.5950698   1.06737575
#> Neutrophils_MCPcounter               -0.7049323    2.1498366   1.13214726
#> Endothelial_cells_MCPcounter         -0.8480613    0.3722122   0.03108201
#> Fibroblasts_MCPcounter               -1.8821677   -0.4083560  -1.09185606
#>                                    TCGA-VQ-A8E2 TCGA-VQ-A8E3 TCGA-VQ-A8E7
#> T_cells_MCPcounter                  0.825241649   0.63966806  -1.16811578
#> CD8_T_cells_MCPcounter              1.458099804   0.05947011  -0.44567782
#> Cytotoxic_lymphocytes_MCPcounter    0.534266939   1.78123793  -1.23152001
#> NK_cells_MCPcounter                 0.076775489   4.24287334  -1.18530661
#> B_lineage_MCPcounter               -0.528210121  -0.28104688  -1.11207860
#> Monocytic_lineage_MCPcounter        0.096690774  -0.01416675  -1.47752409
#> Myeloid_dendritic_cells_MCPcounter  0.001160666   0.51307798  -0.81001256
#> Neutrophils_MCPcounter              0.026884587   0.09385457   0.14553339
#> Endothelial_cells_MCPcounter        1.122562105   0.03683041  -0.01911034
#> Fibroblasts_MCPcounter              1.123522839  -0.63262127   0.04038854
#>                                    TCGA-VQ-A8P2 TCGA-VQ-A8P3 TCGA-VQ-A8P5
#> T_cells_MCPcounter                   -1.3306285   -0.6937524  -0.33947737
#> CD8_T_cells_MCPcounter               -1.3215708   -0.3604799  -0.50106650
#> Cytotoxic_lymphocytes_MCPcounter     -0.8922454    0.1839203   0.30425010
#> NK_cells_MCPcounter                  -0.9620793    0.1257733   1.36145526
#> B_lineage_MCPcounter                 -1.2391117   -1.1251757  -0.55595897
#> Monocytic_lineage_MCPcounter         -1.8420496   -1.7235407  -0.06944349
#> Myeloid_dendritic_cells_MCPcounter   -0.6344778   -0.6369254  -0.40372276
#> Neutrophils_MCPcounter               -1.2066465   -0.6194391  -0.11629483
#> Endothelial_cells_MCPcounter         -0.9894670   -0.3234961  -0.21855734
#> Fibroblasts_MCPcounter               -1.1165483    0.1250300   0.49008528
#>                                    TCGA-VQ-A8P8 TCGA-VQ-A8PB TCGA-VQ-A8PC
#> T_cells_MCPcounter                   0.02538795  -0.77973724   0.42777061
#> CD8_T_cells_MCPcounter              -0.67378084  -0.39416111  -0.52804703
#> Cytotoxic_lymphocytes_MCPcounter    -0.85982138   0.32023845  -0.20700396
#> NK_cells_MCPcounter                 -0.78388559  -0.69037588  -0.18754659
#> B_lineage_MCPcounter                 2.04363032  -0.73210182   1.36114591
#> Monocytic_lineage_MCPcounter        -0.74732953  -0.85951699   1.32384273
#> Myeloid_dendritic_cells_MCPcounter   0.54118119  -0.77355358   0.15456893
#> Neutrophils_MCPcounter              -0.76087876  -0.86971261   1.59854139
#> Endothelial_cells_MCPcounter         0.72183500  -0.34133614   0.88686161
#> Fibroblasts_MCPcounter               0.12715029   0.04614483  -0.05378996
#>                                    TCGA-VQ-A8PD TCGA-VQ-A8PE TCGA-VQ-A8PF
#> T_cells_MCPcounter                    1.9043170   0.80748146   1.74947983
#> CD8_T_cells_MCPcounter                0.6136823   1.25005760   1.44705853
#> Cytotoxic_lymphocytes_MCPcounter      1.5645132   0.15000304   2.00438454
#> NK_cells_MCPcounter                   0.6552649   0.65351061   1.13251477
#> B_lineage_MCPcounter                  2.1684814   1.30600006   0.98063353
#> Monocytic_lineage_MCPcounter          0.8253319   0.42085912   1.99142639
#> Myeloid_dendritic_cells_MCPcounter    1.5331811   0.03256724   1.89936010
#> Neutrophils_MCPcounter                1.9083227   1.38456459   0.24376845
#> Endothelial_cells_MCPcounter          1.5174484  -0.10940504   1.48789487
#> Fibroblasts_MCPcounter                0.7289977  -0.34583703   0.03202277
#>                                    TCGA-VQ-A8PH TCGA-VQ-A8PJ TCGA-VQ-A8PK
#> T_cells_MCPcounter                    0.4085762  -1.13950079  -0.69044345
#> CD8_T_cells_MCPcounter                1.4301148   1.26936578   0.79807456
#> Cytotoxic_lymphocytes_MCPcounter     -0.3020192  -1.16222774  -0.06303456
#> NK_cells_MCPcounter                  -1.0784503  -0.66150760   1.45907532
#> B_lineage_MCPcounter                  0.2569651  -0.74666765  -0.72404679
#> Monocytic_lineage_MCPcounter         -1.2042254  -0.70366686  -0.37723783
#> Myeloid_dendritic_cells_MCPcounter   -0.3370582   0.01434138  -0.88373451
#> Neutrophils_MCPcounter               -0.7562184   1.00232665   0.62016213
#> Endothelial_cells_MCPcounter         -1.1767931   0.16078296  -0.43113759
#> Fibroblasts_MCPcounter               -1.2788268  -0.86397407   0.58445504
#>                                    TCGA-VQ-A8PM TCGA-VQ-A8PO TCGA-VQ-A8PP
#> T_cells_MCPcounter                   0.39084241    0.9849134    0.4354941
#> CD8_T_cells_MCPcounter              -0.25477244    0.6347043    0.2862798
#> Cytotoxic_lymphocytes_MCPcounter     0.06177224    0.8586062    1.1700466
#> NK_cells_MCPcounter                 -0.18003425    0.3160514   -0.1124643
#> B_lineage_MCPcounter                 1.54053497    1.4043230    0.7848588
#> Monocytic_lineage_MCPcounter        -0.02482044    0.1307744    0.1444482
#> Myeloid_dendritic_cells_MCPcounter   0.23080349    0.1260239   -0.5993433
#> Neutrophils_MCPcounter               0.17233656   -0.9816933   -0.5405659
#> Endothelial_cells_MCPcounter         0.90054250   -0.4969906   -0.2906012
#> Fibroblasts_MCPcounter              -0.08763784    0.1675277    0.9149664
#>                                    TCGA-VQ-A8PQ TCGA-VQ-A8PU TCGA-VQ-A8PX
#> T_cells_MCPcounter                    3.9194964   -1.8135576   0.31695493
#> CD8_T_cells_MCPcounter                1.6844221   -1.6905255  -0.18317533
#> Cytotoxic_lymphocytes_MCPcounter      1.3375325   -1.5333928  -0.40451776
#> NK_cells_MCPcounter                   1.4740865   -0.7708626  -1.06058797
#> B_lineage_MCPcounter                  3.0859800   -1.2599498  -0.04774783
#> Monocytic_lineage_MCPcounter          1.2176885   -1.4297399  -1.04663374
#> Myeloid_dendritic_cells_MCPcounter    2.0506523   -1.4749110   0.36061338
#> Neutrophils_MCPcounter                0.9620854    0.4670867   1.20778232
#> Endothelial_cells_MCPcounter          1.2340100   -0.7986368  -0.67206932
#> Fibroblasts_MCPcounter                0.3753419   -1.1030014  -1.99380731
#>                                    TCGA-VQ-A91A TCGA-VQ-A91D TCGA-VQ-A91E
#> T_cells_MCPcounter                   -0.8854020  -0.57953968    0.5710609
#> CD8_T_cells_MCPcounter               -1.3471680  -0.98355577    0.2584934
#> Cytotoxic_lymphocytes_MCPcounter     -0.4493499  -0.06781445    0.9525780
#> NK_cells_MCPcounter                  -0.2109057   0.71931988    0.4430208
#> B_lineage_MCPcounter                 -0.5853957  -1.14097547    0.9824180
#> Monocytic_lineage_MCPcounter          1.1158066   0.01879422   -0.6900950
#> Myeloid_dendritic_cells_MCPcounter   -0.2915582  -1.36550131    0.6364613
#> Neutrophils_MCPcounter                1.2427508  -0.44243017    0.3361117
#> Endothelial_cells_MCPcounter          1.6181221  -0.87947532   -0.9718754
#> Fibroblasts_MCPcounter                1.5886857   0.32965115   -1.2939760
#>                                    TCGA-VQ-A91K TCGA-VQ-A91N TCGA-VQ-A91Q
#> T_cells_MCPcounter                   0.05799075  -0.38964998   -0.6322600
#> CD8_T_cells_MCPcounter              -0.19433126  -0.62095746   -1.4797608
#> Cytotoxic_lymphocytes_MCPcounter    -0.18327051  -0.04744502   -0.3971273
#> NK_cells_MCPcounter                 -0.83814177   0.22139538   -0.3288821
#> B_lineage_MCPcounter                -0.91047418   0.20118542   -0.5090544
#> Monocytic_lineage_MCPcounter         0.77800108  -1.07120830   -0.6079447
#> Myeloid_dendritic_cells_MCPcounter  -0.72303626   0.62230895    0.4633617
#> Neutrophils_MCPcounter              -0.28373957   0.51524216    1.2278852
#> Endothelial_cells_MCPcounter        -0.57283885  -0.29916188    1.6444886
#> Fibroblasts_MCPcounter               0.35649092  -0.89021250    0.9093250
#>                                    TCGA-VQ-A91S TCGA-VQ-A91U TCGA-VQ-A91V
#> T_cells_MCPcounter                   0.52702223   0.56382862   -1.1408926
#> CD8_T_cells_MCPcounter               1.82835982   0.36536388   -1.0030820
#> Cytotoxic_lymphocytes_MCPcounter     1.16152324  -0.07835340   -1.1332786
#> NK_cells_MCPcounter                  0.67534790  -0.23760076   -0.4829409
#> B_lineage_MCPcounter                -0.93014816   0.16908074   -0.7621869
#> Monocytic_lineage_MCPcounter        -0.17205727  -0.18469439   -1.4007454
#> Myeloid_dendritic_cells_MCPcounter  -0.07529549  -0.55647831   -0.7424123
#> Neutrophils_MCPcounter               0.90732638  -0.69587838    0.2981380
#> Endothelial_cells_MCPcounter         0.31219546   0.09092133   -0.6965436
#> Fibroblasts_MCPcounter              -0.74985136  -0.78083794   -1.4905419
#>                                    TCGA-VQ-A91X TCGA-VQ-A91Y TCGA-VQ-A91Z
#> T_cells_MCPcounter                   -1.3899727  1.013839396   -1.0185866
#> CD8_T_cells_MCPcounter               -1.5737090 -0.001782453    1.1085102
#> Cytotoxic_lymphocytes_MCPcounter     -1.7775097  0.765928251   -1.5640511
#> NK_cells_MCPcounter                  -1.4712716  0.791524182   -1.4108329
#> B_lineage_MCPcounter                 -0.5601535  0.359047888   -1.2578040
#> Monocytic_lineage_MCPcounter         -2.4000914  1.370450648   -1.9391015
#> Myeloid_dendritic_cells_MCPcounter   -1.8494242  0.217610929   -1.4025527
#> Neutrophils_MCPcounter               -0.1474874 -0.631330208    0.2028961
#> Endothelial_cells_MCPcounter         -2.7068557  0.621042804   -1.6034610
#> Fibroblasts_MCPcounter               -2.3930615  1.306091760   -1.3704546
#>                                    TCGA-VQ-A922 TCGA-VQ-A923 TCGA-VQ-A924
#> T_cells_MCPcounter                  -0.64728122    1.9883658   0.28183745
#> CD8_T_cells_MCPcounter              -0.20459875    1.5574287   0.10221624
#> Cytotoxic_lymphocytes_MCPcounter    -1.06330085    1.8832406   0.78970509
#> NK_cells_MCPcounter                 -0.35099259    0.9286879   0.61468854
#> B_lineage_MCPcounter                 0.05070682    1.3150796  -0.66021810
#> Monocytic_lineage_MCPcounter        -0.35151169    0.9870137   0.39291368
#> Myeloid_dendritic_cells_MCPcounter  -0.77811166    1.5423870  -0.09712303
#> Neutrophils_MCPcounter               0.25114027   -0.5319076  -0.96334667
#> Endothelial_cells_MCPcounter         1.67919846    0.9027571  -0.78684477
#> Fibroblasts_MCPcounter               1.15440855   -0.4424755   0.32428157
#>                                    TCGA-VQ-A925 TCGA-VQ-A927 TCGA-VQ-A928
#> T_cells_MCPcounter                   -0.1183874  -0.23837078   -0.7071312
#> CD8_T_cells_MCPcounter               -0.2745005  -0.81989850   -1.3109649
#> Cytotoxic_lymphocytes_MCPcounter     -0.3626886  -0.58952857   -0.9011235
#> NK_cells_MCPcounter                  -0.7539923   0.46317842   -0.6463232
#> B_lineage_MCPcounter                 -0.4478396  -0.12595077   -0.6217727
#> Monocytic_lineage_MCPcounter         -0.9706010  -0.96135403   -0.5516056
#> Myeloid_dendritic_cells_MCPcounter    0.6241392  -0.62654080   -1.0030163
#> Neutrophils_MCPcounter                0.1204876   0.48223956    0.4417220
#> Endothelial_cells_MCPcounter          0.2331063   0.80712125    0.7196718
#> Fibroblasts_MCPcounter               -0.4855902  -0.02327551    0.8304147
#>                                    TCGA-VQ-A92D TCGA-VQ-A94O TCGA-VQ-A94P
#> T_cells_MCPcounter                 -1.444886224   0.75750930  0.099571389
#> CD8_T_cells_MCPcounter             -1.220539837  -0.37419509  0.153581903
#> Cytotoxic_lymphocytes_MCPcounter   -1.188391417   0.11182320  0.306657426
#> NK_cells_MCPcounter                -0.842002850  -0.78652401  0.120646960
#> B_lineage_MCPcounter               -0.887764998  -0.57606798  0.001626009
#> Monocytic_lineage_MCPcounter       -1.386261344  -0.04247667 -0.114362084
#> Myeloid_dendritic_cells_MCPcounter -1.666959752  -0.96205842 -0.008797279
#> Neutrophils_MCPcounter             -1.186982264   0.52976588 -0.454117986
#> Endothelial_cells_MCPcounter       -0.272234223   0.34454483  1.931238342
#> Fibroblasts_MCPcounter              0.002508071  -0.57584604  1.759783493
#>                                    TCGA-VQ-A94R TCGA-VQ-A94T TCGA-VQ-A94U
#> T_cells_MCPcounter                   0.19772519   -1.7651996  -0.73509199
#> CD8_T_cells_MCPcounter               0.02580576   -1.5808421  -1.11505480
#> Cytotoxic_lymphocytes_MCPcounter    -0.27230367   -0.5192142  -1.19363235
#> NK_cells_MCPcounter                 -0.03351967    0.2726331  -0.08342231
#> B_lineage_MCPcounter                 0.16893489   -0.9571955  -0.82893554
#> Monocytic_lineage_MCPcounter         0.23453819   -0.8307437  -1.06498224
#> Myeloid_dendritic_cells_MCPcounter  -0.34844409   -0.4876545   0.06769843
#> Neutrophils_MCPcounter              -0.43444321    0.5172288  -0.19679327
#> Endothelial_cells_MCPcounter         0.54034324    0.3197640   1.41386142
#> Fibroblasts_MCPcounter               0.56371734   -0.4003001   0.62100835
#>                                    TCGA-VQ-AA64 TCGA-VQ-AA68 TCGA-VQ-AA69
#> T_cells_MCPcounter                  -0.43344878   0.30100552   0.28159045
#> CD8_T_cells_MCPcounter              -0.20651712  -0.02714418   0.05947647
#> Cytotoxic_lymphocytes_MCPcounter    -0.64449230  -0.05531259   0.62223271
#> NK_cells_MCPcounter                 -0.59398485  -0.45050517  -0.48460428
#> B_lineage_MCPcounter                 0.04253756   0.81020505  -1.39442660
#> Monocytic_lineage_MCPcounter        -0.50091195  -0.82761638  -0.23651703
#> Myeloid_dendritic_cells_MCPcounter  -0.80956657   0.19785496  -1.33320214
#> Neutrophils_MCPcounter               0.23188797   0.10625567  -0.11764652
#> Endothelial_cells_MCPcounter         1.09663984  -0.79599310  -0.87024446
#> Fibroblasts_MCPcounter               1.57517523  -0.78428084  -1.47691352
#>                                    TCGA-VQ-AA6A TCGA-VQ-AA6D TCGA-VQ-AA6F
#> T_cells_MCPcounter                   -0.8237120   -0.7872447   0.22228054
#> CD8_T_cells_MCPcounter               -1.3883873    1.1486359  -0.40376641
#> Cytotoxic_lymphocytes_MCPcounter     -1.0051090   -1.3118081  -0.09080977
#> NK_cells_MCPcounter                  -0.4756096   -1.2519015   0.36817597
#> B_lineage_MCPcounter                 -0.6722927   -0.9753169   1.60771986
#> Monocytic_lineage_MCPcounter         -0.8530319   -0.6201404  -0.48196602
#> Myeloid_dendritic_cells_MCPcounter   -0.5403599   -1.3396724  -0.01231494
#> Neutrophils_MCPcounter               -0.3750900    1.1227178   0.81077608
#> Endothelial_cells_MCPcounter         -0.7685050   -1.0604408   0.28221025
#> Fibroblasts_MCPcounter               -0.7311119   -1.0843382  -0.82701498
#>                                    TCGA-VQ-AA6G TCGA-VQ-AA6J TCGA-VQ-AA6K
#> T_cells_MCPcounter                   -1.0657632   1.07235070   0.02612623
#> CD8_T_cells_MCPcounter               -1.4257562   0.86307904   0.31360601
#> Cytotoxic_lymphocytes_MCPcounter     -1.1960533   0.08052617   0.64121313
#> NK_cells_MCPcounter                  -0.7594583  -0.37469887  -0.09462378
#> B_lineage_MCPcounter                 -0.5193979   1.54078874  -0.44436004
#> Monocytic_lineage_MCPcounter         -0.5271932   1.40973253  -0.63003169
#> Myeloid_dendritic_cells_MCPcounter   -0.5207764   1.25100190  -0.64118712
#> Neutrophils_MCPcounter                0.1350932   0.06711491  -0.02073155
#> Endothelial_cells_MCPcounter          0.3000447   0.18324288   0.02472273
#> Fibroblasts_MCPcounter               -0.8569854   0.01658955   0.16593439
#>                                    TCGA-ZA-A8F6 TCGA-ZQ-A9CR
#> T_cells_MCPcounter                   -0.3292301  -0.29934761
#> CD8_T_cells_MCPcounter               -0.8844827  -0.13962526
#> Cytotoxic_lymphocytes_MCPcounter     -0.7320959   0.47456564
#> NK_cells_MCPcounter                  -0.1752463  -0.32894159
#> B_lineage_MCPcounter                 -0.5891289  -0.42281958
#> Monocytic_lineage_MCPcounter          0.8619831  -0.47593928
#> Myeloid_dendritic_cells_MCPcounter    0.1994664  -0.23693035
#> Neutrophils_MCPcounter               -0.1168008   0.03737702
#> Endothelial_cells_MCPcounter          0.4382820   1.14822654
#> Fibroblasts_MCPcounter                0.6898235   1.01548718
#> attr(,"scaled:center")
#>                 T_cells_MCPcounter             CD8_T_cells_MCPcounter 
#>                          2.1272195                          2.3339357 
#>   Cytotoxic_lymphocytes_MCPcounter                NK_cells_MCPcounter 
#>                          1.4836188                          0.5761658 
#>               B_lineage_MCPcounter       Monocytic_lineage_MCPcounter 
#>                          3.0253606                          3.1534554 
#> Myeloid_dendritic_cells_MCPcounter             Neutrophils_MCPcounter 
#>                          1.7693935                          2.0361769 
#>       Endothelial_cells_MCPcounter             Fibroblasts_MCPcounter 
#>                          2.9390371                          6.7645463 
#> attr(,"scaled:scale")
#>                 T_cells_MCPcounter             CD8_T_cells_MCPcounter 
#>                          0.7315805                          1.2765545 
#>   Cytotoxic_lymphocytes_MCPcounter                NK_cells_MCPcounter 
#>                          0.6648504                          0.3270494 
#>               B_lineage_MCPcounter       Monocytic_lineage_MCPcounter 
#>                          1.3620574                          0.8511569 
#> Myeloid_dendritic_cells_MCPcounter             Neutrophils_MCPcounter 
#>                          0.6446256                          0.4470047 
#>       Endothelial_cells_MCPcounter             Fibroblasts_MCPcounter 
#>                          0.6509202                          1.2078533 
#> 
# }
```
