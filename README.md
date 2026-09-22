# ST-COS Tutorial

This repository contains the source and rendered pages for the
[ST-COS tutorial](https://yiyannnnn.github.io/STCOS_tutorial/).

## Build the site

Install the current ST-COS package, then render all pages from the repository
root:

```r
devtools::install_github("Yiyannnnn/ST-COS")
rmarkdown::render_site(encoding = "UTF-8")
```

The navigation and shared output settings are defined in `_site.yml`. Edit the
`.Rmd` source files rather than the generated `.html` files.
