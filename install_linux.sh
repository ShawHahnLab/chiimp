#!/usr/bin/env bash

# Install CHIIMP on Linux.

# A base R install is assumed to already be present, but all dependencies
# should be installed automatically here.

rexe=$(which R)
pkgdir=$(readlink -f $(dirname $BASH_SOURCE))
pkgdir_r=$pkgdir:\=\\%\\.

devtools_setup="install.packages('devtools',repos='https://cloud.r-project.org')"
bioclite_setup="source('https://bioconductor.org/biocLite.R');biocLite();biocLite('msa')"
deps_setup="devtools::install_deps('$pkgdir_r',dependencies=TRUE)"
chiimp_test="quit(save='no',status=sum(as.data.frame(devtools::test('$pkgdir_r'))\$failed))"
chiimp_setup="devtools::install('$pkgdir_r')"
chiimp_get_path="cat(system.file('bin','chiimp',package='chiimp'))"

"$rexe" --version
echo
echo "### Installing devtools"
echo
"$rexe" --slave -e "$devtools_setup"
echo
echo "### Installing Bioconductor and MSA"
echo
"$rexe" --slave -e "$bioclite_setup"
echo
echo "### Installing dependencies"
echo
"$rexe" --slave -e "$deps_setup"
echo
echo "### Testing CHIIMP"
echo
if "$rexe" --slave -e "$chiimp_test"
then
	echo
	echo
	echo  "   Warning: Tests indicated failures."
	echo
	echo
fi
echo
echo "### Installing CHIIMP"
echo
"$rexe" --slave -e "$chiimp_setup"
