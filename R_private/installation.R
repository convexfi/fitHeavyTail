##
## User installation
##
install.packages(file.choose(), repos = NULL, type="source")
library(covHeavyTail)
help(package="covHeavyTail")
?momentsStudentt




##
## Developer commands (http://r-pkgs.had.co.nz/)
##
#devtools::create("covHeavyTail")
devtools::load_all()  #or Ctrl-Shift-L
#devtools::use_package("Gmedian")
devtools::document()
devtools::install()
devtools::build()  # to generate the installation file
#devtools::install_github("yourusername/myfirstpackage")  # for GitHub
#devtools::release()  #for CRAN
#devtools::use_readme_rmd()  # to create the README file
#rmarkdown::render("README.Rmd")  # to create .Rd from .Rmd

#code checking
lintr::lint_package()
devtools::check()

covr::package_coverage()  #coverage of tests
goodpractice::gp()
