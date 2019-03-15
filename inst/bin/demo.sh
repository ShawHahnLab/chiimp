#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "chiimp" executable
# R script and the test data.

cd "$(dirname $BASH_SOURCE)"
dir=$(pwd -P)
scratch=$(mktemp -d)
R --vanilla -q -e "devtools::load_all('..', quiet=T); test_data\$write_seqs(test_data\$seqs, '$scratch/str-dataset', 'Replicate1-Sample%s-%s.fasta')" > /dev/null
cp "$dir/../example_locus_attrs.csv" "$scratch/locus_attrs.csv"
cp "$dir/../example_config.yml" "$scratch/config.yml"
cd "$scratch"
"$dir/chiimp" "config.yml"
echo "Demo files and output stored in $(pwd -P)"
