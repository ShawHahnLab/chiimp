#!/usr/bin/env bash

set -e

VERSION=$1

chiimp_check='x<-devtools::check();quit(save="no",status=length(c(x$errors,x$warnings,x$notes)))'

sed -i "s/Version: .*$/Version: $VERSION/" DESCRIPTION
R --slave --vanilla -e "$chiimp_check"
R --slave --vanilla -e "rmarkdown::render('GUIDE.Rmd', output_file = 'GUIDE.pdf', quiet = TRUE)"
