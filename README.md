# *jaccard* R package

This R package enables statistical testing of similarity between binary data using the Jaccard/Tanimoto similarity coefficient -- the ratio of intersection to union. Biochemical fingerprints, genomic intervals, and ecological communities are some examples of binary data in life sciences. For examples, competition between two different operational taxonomic units (OTUs) are often evaluated by a Jaccard/Tanimoto coefficient between their absence/presence vectors across multiple bioregions.

We provide 4 methods of computing statistical significance of such similarity coefficients for binary data: the exact solution, the asymptotic approximation, the bootstrap method, and the measure concentration algorithm. We recommand using either the bootstrap method or the measure concentration algorithm, since the exact solution can be slow and the asymptotic approximation could be inaccurate depending on the data size.

# Basic Usage

To install this package:

```R
install.packages("devtools")
library("devtools")
install_github("ncchung/jaccard")
```

To load this package:

```R
library("jaccard")
```

To compute the exact p-value of similarity between two binary vectors, x and y:

```R
jaccard.test(x,y,method="exact")$pvalue
```

When a length of a binary vector is moderately long, the bootstrap and the measure concentration algorithm (mca) are much faster while maintaining high accuracy:

```R
jaccard.test(x,y,method="bootstrap",B=1000)
jaccard.test(x,y,method="mca",accuracy=1e-05)
```

Help documents can be loaded in R, such as:

```R
? jaccard.test.bootstrap
```


# License

[GNU General Public License 2](https://opensource.org/licenses/GPL-2.0)
