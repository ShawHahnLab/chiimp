#!/usr/bin/env bash

set -e

VERSION=$1

chiimp_check='x<-devtools::check();quit(save="no",status=length(c(x$errors,x$warnings,x$notes)))'

# Update version in download link in README
VER_MSG="The most recent released version is"
TAG_URL="https\\://github.com/ShawHahnLab/chiimp/releases/tag"
SED_README="s:$VER_MSG \\[[0-9.]+\\]\\($TAG_URL/[0-9.]+\\)\\.:$VER_MSG [$VERSION]($TAG_URL/$VERSION).:"
sed -i -r "$SED_README" README.md

# Update version in DESCRIPTION
sed -i "s/Version: .*$/Version: $VERSION/" DESCRIPTION

R --slave --vanilla -e "$chiimp_check"
R --slave --vanilla -e "rmarkdown::render('GUIDE.Rmd', output_file = 'GUIDE.pdf', quiet = TRUE)"
