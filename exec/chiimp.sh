#!/usr/bin/env bash

# Linux wrapper to the CHIIMP command-line script.
#
# This script will wait for a keypress before exiting since it presumably
# opened its own terminal window.

# This directory
dir=$(readlink -f $(dirname $BASH_SOURCE))
cfg_dir=$(dirname "$1")

if [[ $# -eq 0 ]]; then
	echo "To run CHIIMP, drag and drop a configuration file onto this icon."
	echo
	echo "For more information see the user guide bundled with the program or here:"
	echo "https://shawhahnlab.github.io/chiimp/GUIDE.pdf"
	echo
else
	export RSTUDIO_PANDOC=$(Rscript "$dir/find_pandoc.R")
	cd "$cfg_dir"
	Rscript "$dir/chiimp" $*
fi
read -p "Press any key to continue... " -n1 -s
echo
