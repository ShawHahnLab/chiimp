#!/usr/bin/env bash

# Linux/Mac wrapper to the microsat command-line script.
#
# This script will wait for a keypress before exiting since it presumably
# opened its own terminal window.

# This directory
dir=$(readlink -f $(dirname $BASH_SOURCE))

Rscript "$dir/microsat" $*
read -p "Press any key to continue... " -n1 -s
echo