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
library(devtools)
#devtools::create("covHeavyTail")
devtools::load_all()  #or Ctrl-Shift-L
#devtools::use_package("Gmedian")
#devtools::use_package("mvtnorm")
devtools::document()
devtools::install()
devtools::build()  # to generate the installation file
#devtools::use_readme_rmd()  # to create the README file
#devtools::use_data_raw()  # to set up the raw-data folder

# code checking
lintr::lint_package()
devtools::check()

# code tests
#devtools::use_testthat()
devtools::test()
covr::package_coverage()  #coverage of tests


# overall checks:
#goodpractice::gp()

# CRAN check
rcmdcheck::rcmdcheck()
#R CMD check --as-cran  # this is before submission to CRAN



# to upload to CRAN or GitHub
#devtools::install_github("yourusername/myfirstpackage")  # for GitHub
#devtools::release()  #for CRAN
