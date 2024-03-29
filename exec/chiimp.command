#!/usr/bin/env bash

# Mac OS wrapper to the CHIIMP command-line script.
#
# This script will wait for a keypress before exiting since it presumably
# opened its own terminal window and Terminal may or may not be configured to
# auto-close windows when commands finish.  If CHIIMP_AUTOCLOSE is set to yes,
# this is skipped.

# This directory
dir=$(dirname $BASH_SOURCE)
cfg_dir=$(dirname "$1")

if [[ $# -eq 0 ]]; then
	echo "To run CHIIMP, drag and drop a configuration file onto this icon."
	echo
	echo "For more information see the user guide bundled with the program or here:"
	echo "https://shawhahnlab.github.io/chiimp/GUIDE.pdf"
	echo
else
	if [[ ! "$dir" =~ ^/ ]]; then
		rel=$(pwd)
	fi
	export RSTUDIO_PANDOC=$(Rscript "$rel/$dir/find_pandoc.R")
	cd "$cfg_dir"
	Rscript "$rel/$dir/chiimp" "$@"
fi
if [[ "$CHIIMP_AUTOCLOSE" != "yes" ]]; then
	read -p "Press any key to continue... " -n1 -s
	echo
fi
