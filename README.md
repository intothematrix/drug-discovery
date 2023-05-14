# drug-discovery
Drug discovery bio project

## setup
Assuming we are using MacOS...

1. Install [Homebrew](https://brew.sh) (MacOS package/software manager):
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

2. Install R and RStudio
[R](https://www.r-project.org) is a programming language and environment for statistical, data analysis and other research. [RStudio](https://posit.co/downloads/) is the most popular development environment for R.

```Bash
brew install --cask r
brew install --cask rstudio
```

3. Get familiar with R programming and RStudio
   * [R for Data Science (a free online book)](https://r4ds.had.co.nz/index.html)
   * [An Introduction to R (another free online book)](https://intro2r.com)
   * [Data Analysis and Visualization Using R (video tutorials)](http://varianceexplained.org/RData/)
   * [swirl: an interactive tutorial to learn R in R environment](https://swirlstats.com/students.html)
   * [Data Visualization with R](https://rkabacoff.github.io/datavis/)

Other Learning Materials:
   * [Awesome R: a curated list](https://github.com/qinwf/awesome-R)
   * [R Data Science Tutorials: another curated list](https://github.com/ujjwalkarn/DataScienceR)

4. Set up DepMap

In RStudio, create a new project.

Then, run these commands:
```R
install.packages("BiocManager")
BiocManager::install("depmap", force = TRUE)
```

Try out one example following instructions documented [here](https://bioconductor.org/packages/release/data/experiment/vignettes/depmap/inst/doc/using_depmap.html).

Other References:
   * [depmap](https://bioconductor.org/packages/devel/data/experiment/html/depmap.html)
   * [depmap data](https://bioconductor.org/packages/devel/data/experiment/vignettes/depmap/inst/doc/depmap.html)
   * [using the depmap data](https://bioconductor.org/packages/release/data/experiment/vignettes/depmap/inst/doc/using_depmap.html)


## project (TBD)

