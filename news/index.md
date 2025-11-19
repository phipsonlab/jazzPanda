# Changelog

## jazzPanda 1.1.2

*Updated on Oct 14, 2025*

- Minor bug fix for Bioconductor devel.

## jazzPanda 1.0.1

*Updated on Sep 1, 2025*

Changes to the package include: - **function
[`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)**:  
- Buffer used for auto-detecting spatial windows has been changed from
5% to a fixed value of `1e-6`.  
- Added a new parameter to allow exporting the **calculated bounding
box** for each sample.

- **function
  [`compute_permp()`](https://phipsonlab.github.io/jazzPanda/reference/compute_permp.md)**:
  - Buffer for automatic spatial boundaries is now set to `1e-6`.  
  - Fixed spatial bin indexing issues to ensure correct alignment when
    calcualting spatial vectors with count matrix.

## jazzPanda 1.0.1

- Backported updates from version 1.1.1 on Bioconductor devel

## jazzPanda 1.1.1

*Updated on June 24, 2025*

Changes to the package include:

- **function
  [`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)**:
  - Parameters `w_x` and `w_y` have been removed. Spatial windows are
    now automatically determined **per sample** using a 5% buffer around
    the coordinate range.
  - Improved handling of genes that are only detected in a subset of
    samples. If a gene is not present in a sample, a zero-filled vector
    is used for that sample to maintain consistent vector lengths across
    samples.
- **function
  [`create_genesets()`](https://phipsonlab.github.io/jazzPanda/reference/create_genesets.md)**:
  - Parameters `w_x` and `w_y` have been removed. Spatial windows are
    now computed **per sample** using a 5% buffer based on input
    coordinates.  
  - Improved handling of genes that are only detected in a subset of
    samples. If a gene is not present in a sample, a zero-filled vector
    is used for that sample to maintain consistent vector lengths across
    samples.
- **function
  [`compute_permp()`](https://phipsonlab.github.io/jazzPanda/reference/compute_permp.md)**:
  - Parameters `w_x` and `w_y` are no longer required. Spatial
    boundaries are now automatically derived **per sample** using a 5%
    buffer.

## jazzPanda 0.2.2

*Updated on February 23, 2025*

Changes to the package include:

- **function
  [`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)**:
  - Input `x` can also be a named list now. This approach enables users
    to directly input transcript detection coordinates using a
    dataframe. This helps prevent integer overflow issues that may arise
    when creating the SpatialExperiment object for a large dataset

  - If invalid cluster/sample names are detected, the function
    [`make.names()`](https://rdrr.io/r/base/make.names.html) will be
    employed to generate valid names. A message will also be displayed
    to indicate this change.
- **function
  [`create_genesets()`](https://phipsonlab.github.io/jazzPanda/reference/create_genesets.md)**:
  - Input `x` can also be a named list now as in
    [`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)

## jazzPanda 0.2.1

*Updated on January 20, 2025*

Changes to the package include:

- **New Accessors for the Correlation Approach Results**:
  - Add the function
    [`get_cor()`](https://phipsonlab.github.io/jazzPanda/reference/get_cor.md)
  - Add the function
    [`get_perm_p()`](https://phipsonlab.github.io/jazzPanda/reference/get_perm_p.md)
  - Add the function
    [`get_perm_adjp()`](https://phipsonlab.github.io/jazzPanda/reference/get_perm_adjp.md)
- **New Accessors for Results from the Linear Modeling Approach**:
  - Add the function
    [`get_top_mg()`](https://phipsonlab.github.io/jazzPanda/reference/get_top_mg.md)
  - Add the function
    [`get_full_mg()`](https://phipsonlab.github.io/jazzPanda/reference/get_full_mg.md)

## jazzPanda 0.2.0

*Updated on November 16, 2024*

To enhance interoperability with Bioconductor, the jazzPanda package has
been updated to depend on the following data classes:
SingleCellExperiment, SpatialExperiment, or SpatialFeatureExperiment.

Changes to the package include:

- **Removed Function**:
  - The `convert_data()` function has been removed from the export list.
- **function
  [`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)**:
  - Input `trans_lst` is renamed to `x`.
  - Input `cm_lst` is removed.
  - `x` must be a `SingleCellExperiment`, `SpatialExperiment`, or
    `SpatialFeatureExperiment` object.
  - Renamed `all_genes` to `test_genes` to store genes for building
    spatial vectors.
- **function
  [`compute_permp()`](https://phipsonlab.github.io/jazzPanda/reference/compute_permp.md)**:
  - Input `data` is renamed to `x`.
  - `n.cores` is renamed to `n_cores`.
  - `x` must be a `SingleCellExperiment`, `SpatialExperiment`, or
    `SpatialFeatureExperiment` object.
- **function
  [`create_genesets()`](https://phipsonlab.github.io/jazzPanda/reference/create_genesets.md)**:
  - Input `data_lst` is renamed to `x`.
  - `x` must be a `SingleCellExperiment`, `SpatialExperiment`, or
    `SpatialFeatureExperiment` object.

## jazzPanda 0.1.1

*Updated on November 12, 2024*

- Add the `convert_data()` function to expand interoperability with
  Bioconductor. Support input objects from
  SingleCellExperiment/SpatialExperimentSpatialFeatureExperiment class
- Modify one input type for function
  [`get_vectors()`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md).
  Renamed input parameter `data_lst` as `trans_lst` specifically for
  storing transcript coordinates only.

## jazzPanda 0.0.1

- This package is renamed as jazzPanda
