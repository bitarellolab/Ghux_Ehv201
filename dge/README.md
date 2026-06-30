#install packages
Rscript --vanilla pkgs.R

# run
Rscript --vanilla 04-FixMetadata.R
##in R
rmarkdown::render("04-annotations-host-new.Rmd", "04-annotations-host-new.html")
rmarkdown::render("04-annotations-virus-new.Rmd", "04-annotations-virus-new.html")

# Salmon-to-EdgeR

##in R (edit the file to select sampSet and then run:)

rmarkdown::render("05-salmon-to-edgeR.Rmd", "salmon-to-edge-R-sampSet.html")



