#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "chiimp" executable
# R script and the test data.
#
# Special case: all empty input files

if [[ "$2" == "strict" ]]; then
	set -e
fi

cd "$(dirname $BASH_SOURCE)"
dir=$(pwd -P)
inst="../inst"
[ -d "$inst" ] || inst=".."
scratch=$1
if [[ "$scratch" == "" ]]; then
	scratch=$(mktemp -d)
else
	mkdir -p "$scratch"
fi
mkdir -p "$scratch/str-dataset"
for samp in {1..3}; do
	for locus in A B 1 2; do
		touch "$scratch/str-dataset/Replicate1-Sample${samp}-${locus}.fasta"
	done
done
cp "$dir/$inst/example_locus_attrs.csv" "$scratch/locus_attrs.csv"
cp "$dir/$inst/example_config.csv" "$scratch/config.csv"
cd "$scratch"
"$dir/chiimp" "config.csv"
echo "Demo files and output stored in $(pwd -P)"
