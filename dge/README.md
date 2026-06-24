#install packages
Rscript --vanilla pkgs.R

# run
Rscript --vanilla 04-FixMetadata.R
##in R
rmarkdown::render("04-Annotations-Host.R", "04-Annotations-Host.html")

# Salmon-to-EdgeR

##in R:

rmarkdown::render("05-salmon-to-edgeR.Rmd", "salmon-to-edge-R")



