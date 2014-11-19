# `"bigmemoryExt"

`bigmemoryExt` is an R package designed to provide extended features to 'big.matrix' objects from the [bigmemory package](http://cran.r-project.org/web/packages/bigmemory/index.html)

## Motivation
As wonderful as the bigmemory package is, there currently is only limited functionality for the analysis of these objects.  There are other packages extending the capabilities such as bigalgebra and biganalytics but as with any project there is more to done.  This is my humble attempt at creating some extensions to these big.matrix objects that I currently could not find anywhere else.  Ideally this will be a useful package in itself or merged with on of the other bigmemory extension packages (TBD).  As of the moment, these are functions I find useful to apply to these objects.

I have made a conscious effort to heavily document my code so hopefully those who wish to either contribute here or develop there own extenstions will find this useful.

## Installation
Not currently submitted to CRAN

Currently, this project depends on my personal copy of the bigmemory package that is also hosted on [github](https://github.com/cdeterman/bigmemory).  This is clearly not the best scenario but I needed some methods that did not previously exist.  I have added these methods to the big.matrix class (specifically the sharedbigmatrix class) that make some functions herein possible (e.g. in-place transpose).  As such, you must also install my version of bigmemory until such additions are either merged with the authoritative package or improved upon and implemented separately.

Development version on [github](https://github.com/cdeterman/bigmemoryExt)
```r
# development version
library(devtools)

# install my version of bigmemory
install_github('cdeterman/bigmemory')

# install bigmemoryExt
install_github('cdeterman/bigmemoryExt')
```