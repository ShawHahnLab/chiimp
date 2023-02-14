#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "chiimp" executable
# R script and the test data.

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
R --vanilla -q -e "devtools::load_all('..', quiet=T); test_data\$write_seqs(test_data\$seqs, '$scratch/str-dataset', 'Replicate1-Sample%s-%s.fasta')" > /dev/null
cp "$dir/$inst/example_locus_attrs.csv" "$scratch/locus_attrs.csv"
cp "$dir/$inst/example_config.csv" "$scratch/config.csv"
cd "$scratch"
"$dir/chiimp" "config.csv"
echo "Demo files and output stored in $(pwd -P)"
