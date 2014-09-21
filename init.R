library(devtools)
load_all()

run_check <- function(){
  check(args = '-no-manual')
}

make_readme <- function() {
  library(rmarkdown)
  render('README.Rmd', 'md_document')
}