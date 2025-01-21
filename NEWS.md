# jazzPanda 0.2.1
*Updated on January 20, 2025*

Changes to the package include:

- **New Accessors for the Correlation Approach Results**:
    - Add the function `get_cor()`
    - Add the function `get_perm_p()`
    - Add the function `get_perm_adjp()`
    
- **New Accessors for Results from the Linear Modeling Approach**:
    - Add the function `get_top_mg()` 
    - Add the function `get_full_mg()`



# jazzPanda 0.2.0
*Updated on November 16, 2024*

To enhance interoperability with Bioconductor, the jazzPanda package has 
been updated to depend on the following data classes: 
SingleCellExperiment, SpatialExperiment, or SpatialFeatureExperiment.

Changes to the package include:

- **Removed Function**:
  - The `convert_data()` function has been removed from the export list.

- **function `get_vectors()`**:
    - Input `trans_lst` is renamed to `x`.
    - Input `cm_lst` is removed.
    - `x` must be a `SingleCellExperiment`, `SpatialExperiment`, 
        or `SpatialFeatureExperiment` object.
    - Renamed `all_genes` to `test_genes` to store genes 
        for building spatial vectors.

- **function `compute_permp()`**:
    - Input `data` is renamed to `x`.
    - `n.cores` is renamed to `n_cores`.
    - `x` must be a `SingleCellExperiment`, `SpatialExperiment`,
        or `SpatialFeatureExperiment` object.

- **function `create_genesets()`**:
    - Input `data_lst` is renamed to `x`.
    - `x` must be a `SingleCellExperiment`, `SpatialExperiment`, 
        or `SpatialFeatureExperiment` object.


# jazzPanda 0.1.1
*Updated on November 12, 2024*

* Add the `convert_data()` function to expand interoperability with 
Bioconductor. Support input objects from 
SingleCellExperiment/SpatialExperimentSpatialFeatureExperiment class 
* Modify one input type for function `get_vectors()`. 
Renamed input parameter `data_lst` as `trans_lst` specifically for storing
transcript coordinates only. 


# jazzPanda 0.0.1
* This package is renamed as jazzPanda 

# spaceMarker 0.0.2

* Remove the `get_data()` function from the package. 
* The `get_data()` function has been renamed as `get_xenium_data()` and 
    can be used as a help function from vignette. 

# spaceMarker 0.0.2

* Add a parallel option for the `get_vectors()` function. 

# spaceMarker 0.0.1

* First version of the spaceMarker package contains functions to detect marker
genes for image-based spatial transcriptomics data. 
