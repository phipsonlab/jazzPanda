# Convert SingleCellExperiment/SpatialExperiment/SpatialFeatureExperiment objects to list object for jazzPanda.

This function takes an object of class SingleCellExperiment,
SpatialExperiment or SpatialFeatureExperimentreturns and returns a list
object that is expected for the `get_vector` functions.

## Usage

``` r
.convert_data(x, sample_names, test_genes)
```

## Arguments

- x:

  a SingleCellExperiment or SpatialExperiment or
  SpatialFeatureExperiment object

- sample_names:

  a vector of strings giving the sample names

- test_genes:

  A vector of strings giving the name of the genes you want to create
  gene vector.

## Value

outputs a list object with the following components

- trans_lst :

  A list of named dataframes. Each dataframe refers to one sample and
  shows the transcript detection coordinates for each gene. The name
  matches the input sample_names

- cm_lst :

  A list of named dataframes containing the count matrix for each
  sample. The name matches the input sample_names
