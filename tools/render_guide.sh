#!/usr/bin/env bash
R --slave --vanilla -e "rmarkdown::render('GUIDE.Rmd', output_file = 'GUIDE.pdf', quiet = TRUE)"
