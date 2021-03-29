

install:
	R CMD INSTALL .

test_all:
	R --slave -e "library(testthat); test_dir('./tests/testthat/')"

test_lpSolve:
	R --slave -e "library(testthat); test_file('./tests/testthat/test-CARNIVAL_pipeline.R')"

vignettes/CARNIVAL.html:
	R --slave -e "library(rmarkdown); rmarkdown::render('vignettes/CARNIVAL.Rmd', output_format = 'html_vignette')"

vignettes/CARNIVAL.pdf:
	R --slave -e "library(rmarkdown); rmarkdown::render('vignettes/CARNIVAL.Rmd', output_format = 'pdf_document')"

clean:
	rm -rf vignettes/CARNIVAL.html vignettes/CARNIVAL.pdf vignettes/example*/
