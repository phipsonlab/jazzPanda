# jazzPanda 0.3.0
*Updated on June 24, 2025*

Changes to the package include:

- **function `get_vectors()`**:  
  - Parameters `w_x` and `w_y` have been removed. Spatial windows are now 
    automatically determined **per sample** using a 5% buffer around the 
    coordinate range. 
  - Improved handling of genes that are only detected in a subset of samples. 
    If a gene is not present in a sample, a zero-filled vector is used for that
    sample to maintain consistent vector lengths across samples. 

- **function `create_genesets()`**:  
  - Parameters `w_x` and `w_y` have been removed. Spatial windows are now 
    computed **per sample** using a 5% buffer based on input coordinates.  
  - Improved handling of genes that are only detected in a subset of samples. 
    If a gene is not present in a sample, a zero-filled vector is used for that
    sample to maintain consistent vector lengths across samples.
    
- **function `compute_permp()`**:  
  - Parameters `w_x` and `w_y` are no longer required. Spatial boundaries 
    are now automatically derived **per sample** using a 5% buffer.  


# jazzPanda 0.2.2
*Updated on February 23, 2025*

Changes to the package include:

- **function `get_vectors()`**:
    - Input `x` can also be a named list now. This approach enables users to 
    directly input transcript detection coordinates using a dataframe. 
    This helps prevent integer overflow issues that may arise when 
    creating the SpatialExperiment object for a large dataset
    
    - If invalid cluster/sample names are detected, the function `make.names()` 
    will be employed to generate valid names. A message will also be displayed
    to indicate this change.
    
- **function `create_genesets()`**:
    - Input `x` can also be a named list now as in `get_vectors()`


    
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
