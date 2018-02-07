#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "microsat" executable
# R script and the test data.

dir=$(readlink -f $(dirname $BASH_SOURCE))
scratch=$(mktemp -d)
R --vanilla -q -e "devtools::load_all('.', quiet=T); test_data\$write_seqs(test_data\$seqs, '$scratch/str-dataset', 'Replicate1-Sample%s-%s.fasta')"
cp "$dir/../example_locus_attrs.csv" "$scratch/locus_attrs.csv"
cp "$dir/../example_config.yml" "$scratch/config.yml"
cd "$scratch"
"$dir/microsat" "config.yml"
