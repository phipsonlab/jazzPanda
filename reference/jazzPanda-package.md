# jazzPanda: A hybrid approach to find spatially relevant marker genes in image-based spatial transcriptomics data

jazzPanda pacakge provides hybrid approaches to prioritize marker genes
that uses the spatial coordinates of gene detections and cells making up
clusters. We propose a binning approach
[`get_vectors`](https://phipsonlab.github.io/jazzPanda/reference/get_vectors.md)
that summarises the number of genes and cells within a cluster as
spatial vectors. We have developed two approaches to detect and
prioritize marker genes. The first approach
[`compute_permp`](https://phipsonlab.github.io/jazzPanda/reference/compute_permp.md)
is based on correlation coefficients between genes and cluster spatial
vectors, where significance of the marker genes are assessed through
permutation. The second approach
[`lasso_markers`](https://phipsonlab.github.io/jazzPanda/reference/lasso_markers.md)
is based on lasso regularisation and linear modeling of our defined
spatial vectors. This second approach is more flexible and can account
for multiple samples and background noise.

## Author

Melody Jin <jin.m@wehi.edu.au>
