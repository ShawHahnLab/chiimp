#!/usr/bin/env bash

# Install microsat on Linux.

# A base R install is assumed to already be present, but all dependencies
# should be installed automatically here.

rexe=$(which R)
pkgdir=$(readlink -f $(dirname $BASH_SOURCE))
pkgdir_r=$pkgdir:\=\\%\\.

devtools_setup="install.packages('devtools',repos='https://cloud.r-project.org')"
bioclite_setup="source('https://bioconductor.org/biocLite.R');biocLite();biocLite('msa')"
deps_setup="devtools::install_deps('$pkgdir_r',dependencies=TRUE)"
microsat_test="quit(save='no',status=sum(as.data.frame(devtools::test('$pkgdir_r'))\$failed))"
microsat_setup="devtools::install('$pkgdir_r')"
microsat_get_path="cat(system.file('bin','microsat',package='microsat'))"

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
echo "### Testing microsat"
echo
if "$rexe" --slave -e "$microsat_test"
then
	echo
	echo
	echo  "   Warning: Tests indicated failures."
	echo
	echo
fi
echo
echo "### Installing microsat"
echo
"$rexe" --slave -e "$microsat_setup"
