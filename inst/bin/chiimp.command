#!/usr/bin/env bash

# Mac OS wrapper to the CHIIMP command-line script.
#
# This script will wait for a keypress before exiting since it presumably
# opened its own terminal window.

# This directory
dir=$(dirname $BASH_SOURCE)
cfg_dir=$(dirname "$1")

which pandoc > /dev/null || export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc/

if [[ ! "$dir" =~ ^/ ]]; then
	rel=$(pwd)
fi
cd "$cfg_dir"
Rscript "$rel/$dir/chiimp" $*
read -p "Press any key to continue... " -n1 -s
echo
