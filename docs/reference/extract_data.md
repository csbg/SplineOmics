# Extract a rectangular block from a dataframe and make it to a matrix

This function extracts a rectangular block from a dataframe using
user-specified top/bottom row indices and left/right column identifiers
(numeric or Excel-style letters). It ensures the block contains only
numeric values or NAs, and returns a cleaned matrix.

## Usage

``` r
extract_data(
  data,
  bottom_row,
  right_col,
  top_row = 1,
  left_col = 1,
  feature_name_columns = NULL,
  use_row_index = FALSE
)
```

## Arguments

- data:

  `data.frame`: A dataframe containing the full input, including
  annotation columns and the numeric block to extract.

- bottom_row:

  `integer(1)`: Specifies the last (bottom) row of the numeric data
  block. Must be \>= `top_row`.

- right_col:

  `integer(1)` \| `character(1)`: Same format as `left_col`. Specifies
  the right-most column of the numeric block. Must be \>= `left_col`
  after conversion.

- top_row:

  `integer(1)`: Specifies the first (top) row of the numeric data block.
  Row indexing is 1-based.

- left_col:

  `integer(1)` \| `character(1)`: Column specifier for the left-most
  column of the data block. Can be either:

  - An integer index (e.g., 2 for the second column), or

  - A character string using Excel-style letters (e.g., "A", "AB").

  Column names (e.g., "age") are **not** allowed here. Only letters or
  numeric indices are accepted.

- feature_name_columns:

  [`character()`](https://rdrr.io/r/base/character.html) \| `NULL`:
  Optional character vector specifying columns in `data` to be used as
  row (feature) names in the output. If `NA`, generic feature names are
  used. These row names are used everywhere to label the features, such
  as to label the plots in the
  [`cluster_hits()`](https://csbg.github.io/SplineOmics/reference/cluster_hits.md)
  function report.

- use_row_index:

  `logical(1)`: If `TRUE`, prepend the row index to each combined
  feature name to ensure uniqueness. Defaults to `FALSE`.

## Value

A numeric matrix with cleaned values and appropriate column names.

## Examples

``` r
# Tiny demo table with two header rows, feature columns, and numeric block
df <- data.frame(
    feat_id = c(NA, NA, "g1", "g2", "g3"),
    feat_sym = c(NA, NA, "TP53", "EGFR", "BAX"),
    A = c("cond", "t0", 1, 2, 3),
    B = c("cond", "t1", 4, 5, 6),
    C = c("ctrl", "t0", 7, 8, 9),
    D = c("ctrl", "t1", 10, 11, 12),
    check.names = FALSE
)

# Example 1: extract numeric block using Excel letters, build headers from
# the two rows above (they get collapsed like "cond_t0", "ctrl_t1", ...)
m1 <- extract_data(
    data = df,
    top_row = 3,
    bottom_row = 5,
    left_col = "A",
    right_col = "D",
    feature_name_columns = c("feat_id", "feat_sym"),
    use_row_index = FALSE
)
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
m1
#>         X1 X2 cond_t0 cond_t1
#> g1_TP53 NA NA       1       4
#> g2_EGFR NA NA       2       5
#> g3_BAX  NA NA       3       6

# Example 2: same extraction but with numeric column indices and row index
# prepended to ensure uniqueness of feature names
m2 <- extract_data(
    data = df,
    top_row = 3,
    bottom_row = 5,
    left_col = 3,
    right_col = 6,
    feature_name_columns = c("feat_id", "feat_sym"),
    use_row_index = TRUE
)
m2
#>           cond_t0 cond_t1 ctrl_t0 ctrl_t1
#> 1_g1_TP53       1       4       7      10
#> 2_g2_EGFR       2       5       8      11
#> 3_g3_BAX        3       6       9      12
```
