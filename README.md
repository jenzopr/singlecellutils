# Utility functions for data from single-cell experiments

The singlecellutils package provides auxillary functions to work with data from single cell RNAseq data.

### Features

- Easy and versatile sample/feature filtering with `filter_samples` and `filter_features`
- Implements TF-IDF weighting, e.g. for data from single-cell ATAC-seq
- Provides easy access to new dimension reduction methods, like *randomized SVD* and *UMAP*
- Provides access to the powerful hierarchical density-based clustering *HDBSCAN*
- Most methods are suitable for magrittr *pipelines*

### Installation

Install the latest development version from **GitHub**:

```r
remotes::install_github("jenzopr/singlecellutils")
```

### Contribution

Please feel free to contribute to this project - either by opening an [issue](https://github.com/jenzopr/singlecellutils/issues) or a [pull request](https://github.com/jenzopr/singlecellutils/pulls). 

### License

The project is licensed under the MIT license.
