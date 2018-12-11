#!/usr/bin/env bash

# NOTE: assumes working directory is the project directory

set -e

VERSION=$1

chiimp_check='x<-devtools::check();quit(save="no",status=length(c(x$errors,x$warnings)))'

# Run lint script
./.utils/lint.R

# Update version in download link in README
VER_MSG="The most recent released version is"
TAG_URL="https\\://github.com/ShawHahnLab/chiimp/releases/tag"
SED_README="s:$VER_MSG \\[[0-9.]+\\]\\($TAG_URL/[0-9.]+\\)\\.:$VER_MSG [$VERSION]($TAG_URL/$VERSION).:"
sed -i -r "$SED_README" README.md

# Update version in DESCRIPTION and NEWS.md
sed -i "s/Version: .*$/Version: $VERSION/" DESCRIPTION
sed -i "s/# chiimp dev/# chiimp $VERSION/" NEWS.md

R --slave --vanilla -e "$chiimp_check"
R --slave --vanilla -e "rmarkdown::render('GUIDE.Rmd', output_file = 'GUIDE.pdf', quiet = TRUE)"

# Create bundled ZIP and TGZ versions without hidden top level files (such as
# the git and travis stuff) and with the GUIDE.pdf.
pushd ..
zip -r chiimp-v${VERISON}.zip chiimp/*
tar czvf chiimp-v${VERSION}.tgz chiimp/*
popd

echo "REMINDER BEFORE TAGGING RELEASE $VERSION:"
echo
echo " * Run full test on Mac OS, Windows, and Linux"
echo " * Update NEWS.md with all updates under a heading matching this version"
echo " * Make sure GUIDE.Rmd is up-to-date and rendered GUIDE.pdf is correct"
echo
