#!/usr/bin/env bash

# Linux wrapper to the CHIIMP command-line script.
#
# This script will wait for a keypress before exiting since it presumably
# opened its own terminal window.

# This directory
dir=$(readlink -f $(dirname $BASH_SOURCE))
cfg_dir=$(dirname "$1")

cd "$cfg_dir"
Rscript "$dir/chiimp" $*
read -p "Press any key to continue... " -n1 -s
echo
