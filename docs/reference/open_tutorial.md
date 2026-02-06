# open_tutorial()

This function opens the \`tutorial.Rmd\` file in RStudio for interactive
use. Users can then run each code chunk step by step.

## Usage

``` r
open_tutorial()
```

## Value

If successful, opens the \`tutorial.Rmd\` file in RStudio for the user
to interact with. If \`rstudioapi\` is not installed or available, or
the tutorial file is not found, an error is thrown with a corresponding
message.
