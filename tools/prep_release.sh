#!/usr/bin/env bash

# NOTE: assumes working directory is the project directory

set -e

VERSION=$1
SEP="==="

chiimp_check='x<-devtools::check();quit(save="no",status=length(c(x$errors,x$warnings)))'

echo "$SEP Running spell check"
./tools/check_spelling.R

# Run lint script
echo "$SEP Running lint check"
./tools/check_lint.R

if [[ $VERSION != "" ]]; then
	# Update version in download link in README
	VER_MSG="The most recent released version is"
	TAG_URL="https\\://github.com/ShawHahnLab/chiimp/releases/tag"
	SED_README="s:$VER_MSG \\[[0-9.]+\\]\\($TAG_URL/[0-9.]+\\)\\.:$VER_MSG [$VERSION]($TAG_URL/$VERSION).:"
	sed -i -r "$SED_README" README.md

	# Update version in DESCRIPTION and NEWS.md
	sed -i "s/Version: .*$/Version: $VERSION/" DESCRIPTION
	sed -i "s/# chiimp dev/# chiimp $VERSION/" NEWS.md
fi

echo "$SEP Running devtools::check()"
R --slave --vanilla -e "$chiimp_check"

echo "$SEP Rendering user guide"
R --slave --vanilla -e "rmarkdown::render('GUIDE.Rmd', output_file = 'GUIDE.pdf', quiet = TRUE)"

# Create bundled ZIP and TGZ versions without hidden top level files (such as
# the git and travis stuff) and with the GUIDE.pdf.
if [[ $VERSION != "" ]]; then
	echo "$SEP Creating release archives"
	pushd ..
	zip -r chiimp-v${VERSION}.zip chiimp/*
	tar czvf chiimp-v${VERSION}.tgz chiimp/*
	popd
fi

echo
echo "REMINDER BEFORE TAGGING RELEASE $VERSION:"
echo
echo " * Run full test on Mac OS, Windows, and Linux"
echo " * Update NEWS.md with all updates under a heading matching this version"
echo " * Check README.md for link to this version"
echo " * Make sure GUIDE.Rmd is up-to-date and rendered GUIDE.pdf is correct"
echo
echo "ALSO:"
echo " * Draft release from tag on github including archive files with bundled"
echo "   GUIDE.pdf"
echo " * Merge release-### into master, dev, and gh-pages"
echo
