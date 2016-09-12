Query DGIdb using the rDGIdb package
=======

Queries DGIdb through the R/Bioconductor rDGIdb package and saves the results to text files. 

## Installation

rDGIdb is available on [Bioconductor](http://bioconductor.org/packages/rDGIdb/).

To install the rDGIdb package, type this into your R console:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("rDGIdb")
```

## Usage

In the terminal, use:

```bash
Rscript query_dgidb.r path/to/input path/to/output min-db-support
```

## Licence

This script is release under the MIT licence.

## Contact

thomas.thurnherr (at) bsse.ethz.ch
