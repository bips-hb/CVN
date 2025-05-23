## TEST and BUILD

install.packages(c("devtools", "roxygen2", "usethis", "testthat"))
install.packages("rcmdcheck")
library(devtools)
library(roxygen2)
library(rcmdcheck)

## Run roxygen to generate/update NAMESPACE and Rd files
roxygen2::roxygenise()
devtools::spell_check()

## Optional: build vignettes
devtools::build_vignettes()

## Run Tests (if using testthat)
# devtools::test()
# usethis::use_testthat()
# usethis::use_test("core_function")

## Run comprehensive tests
# devtools::check()
devtools::check(vignettes = FALSE, cran = TRUE)
devtools::check(clean = TRUE, cran = TRUE)
rcmdcheck::rcmdcheck()

## Final Cleanup
# Make sure you have:
#   A clean .Rbuildignore file
#   No temp files or system folders in your package directory
#   Valid DESCRIPTION and NAMESPACE files
#   README.md and inst/CITATION if needed


# Build package and submit to
# https://cran.r-project.org/submit.html
devtools::build()
devtools::build_manual()
devtools::build(ignore = TRUE)  # Dateien, die ignoriert werden
