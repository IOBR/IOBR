# Create Nested Output Folders

Creates one to three nested folders (if not existing) under the current
working directory and returns their names and absolute paths.

## Usage

``` r
creat_folder(f1, f2 = NULL, f3 = NULL, return = NULL)
```

## Arguments

- f1:

  Character. First-level folder name.

- f2:

  Character or \`NULL\`. Second-level folder name. Default is \`NULL\`.

- f3:

  Character or \`NULL\`. Third-level folder name. Default is \`NULL\`.

- return:

  Deprecated (not used). Kept for backward compatibility.

## Value

List with elements:

- folder_name:

  Relative path to the created folder

- abspath:

  Absolute path ending with '/'

## Examples

``` r
creat_folder(file.path(tempdir(), "1-result"))
#> $folder_name
#> [1] "/tmp/RtmpYX5KuO/1-result"
#> 
#> $abspath
#> [1] "/tmp/RtmpYX5KuO/1-result/"
#> 
creat_folder(file.path(tempdir(), "1-result"), "figures", "correlation")
#> $folder_name
#> [1] "/tmp/RtmpYX5KuO/1-result/figures/correlation"
#> 
#> $abspath
#> [1] "/tmp/RtmpYX5KuO/1-result/figures/correlation/"
#> 
```
