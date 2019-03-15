#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "chiimp" executable
# R script and the test data.
#
# Special case: all empty input files

cd "$(dirname $BASH_SOURCE)"
dir=$(pwd -P)
scratch=$(mktemp -d)
mkdir -p "$scratch/str-dataset"
for samp in {1..3}; do
	for locus in A B 1 2; do
		touch "$scratch/str-dataset/Replicate1-Sample${samp}-${locus}.fasta"
	done
done
cp "$dir/../example_locus_attrs.csv" "$scratch/locus_attrs.csv"
cp "$dir/../example_config.yml" "$scratch/config.yml"
cd "$scratch"
"$dir/chiimp" "config.yml"
echo "Demo files and output stored in $(pwd -P)"
