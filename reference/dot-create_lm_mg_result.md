# Create a marker gene result object for linear modelling approach

This function creates a structured output object named 'glm_mg_result'
for storing the marker gene results. The object contains two data
frames: top results and full results.

## Usage

``` r
.create_lm_mg_result(top_result, full_result)
```

## Arguments

- top_result:

  A data frame containing top results.

- full_result:

  A data frame containing full results.

## Value

An S3 object of class 'glm_mg_result' which includes both results data
frames.
