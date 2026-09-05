
# Install packages


```{bash}
Rscript --vanilla pkgs.R
```


# Fix/deal with metadata

```{bash}
Rscript --vanilla 04-FixMetadata.R
```

# Annotations

```{r}
#in R
rmarkdown::render("04-annotations-host-new.Rmd", "04-annotations-host-new.html")
rmarkdown::render("04-annotations-virus-new.Rmd", "04-annotations-virus-new.html")
```
# EdgeR analyses

```{r}
#in R (edit the file to select the desired sampSet and then run:)

rmarkdown::render("05-salmon-to-edgeR.Rmd", glue::glue("salmon-to-edge-R-{sampSet}.html"))

```

# deseq2 analyses

```{r}
#in R (edit the file to select the desired sampSet and then run:)

rmarkdown::render("06-salmon-to-deseq2.Rmd", glue::glue("salmon-to-deseq2-{sampSet}.html"))

```
