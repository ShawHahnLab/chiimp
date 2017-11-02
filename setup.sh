#!/usr/bin/env bash

# Set up an install on Linux from scratch.
#
# This uses Miniconda to provide R and encapsulate packages.
# OS-provided packages: texlive, pandoc

PKG="microsat"
CONDA="$HOME/miniconda3"
# Using a specific version until Anaconda gets this issue fixed
# https://stackoverflow.com/questions/46450912/
MINICONDA="Miniconda3-4.3.21-Linux-x86_64.sh"

if [ ! -d "$CONDA" ]
then
	d=$(mktemp -d)
	pushd $d
	wget https://repo.continuum.io/miniconda/"$MINICONDA"
	bash "$MINICONDA" -b -p "$CONDA"
	popd
fi

# Contrast this more general approach with the more exact package matching
# with:
#     conda create --name microsat --file spec-file.txt
$CONDA/bin/conda env list | grep -q "$PKG" ||
	$CONDA/bin/conda env create --file environment.yml

source $CONDA/bin/activate "$PKG"

# Install some fonts needed for R CMD Rd2pdf (via R CMD check)
# https://stackoverflow.com/questions/10819959
LATEX="$CONDA/envs/$PKG/lib/R/share/texmf/tex/latex"
if [ ! -e $LATEX/inconsolata.sty ]
then
	pushd $LATEX
	wget http://mirrors.ctan.org/install/fonts/inconsolata.tds.zip
	unzip inconsolata.tds.zip tex/latex/inconsolata/inconsolata.sty
	mv tex/latex/inconsolata/inconsolata.sty .
	popd
fi
if [ ! -e $LATEX/upquote.sty ]
then
	pushd $LATEX
	wget http://mirrors.ctan.org/macros/latex/contrib/upquote.zip
	unzip upquote.zip upquote/upquote.sty
	mv upquote/upquote.sty .
	popd
fi


# devtools trips up when trying to extract archives unless we explicitly point
# it to the correct tar
# https://github.com/hadley/devtools/issues/379
export TAR=$(which tar)
# Likewise, tools::texi2pdf() by default sets the texi2dvi option to a
# non-existant file, even though it's actually right there on the path.
export R_TEXI2DVICMD=$(which texi2dvi)

# Special install procedure for msa
R -q -e '"msa" %in% installed.packages() || {source("https://bioconductor.org/biocLite.R"); biocLite("msa")}'
# Install microsat itself
R -q -e '"microsat" %in% installed.packages() || devtools::install(".")'

# Add the microsat executable to the conda environment path
[ -e $CONDA/envs/$PKG/bin/microsat ] ||
	ln -s $CONDA/envs/$PKG/lib/R/library/microsat/bin/microsat $CONDA/envs/$PKG/bin
