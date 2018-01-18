#!/usr/bin/env bash

# Wrapper script to set up and execute a demo, using the "microsat" executable
# R script and the test data.

dir=$(readlink -f $(dirname $BASH_SOURCE))
scratch=$(mktemp -d)
R --vanilla -q -e "devtools::load_all('.', quiet=T); write_seqs(seqs, '$scratch/str-dataset', 'Replicate1-Sample%s-%s.fasta')"
cp "$dir/../example_locus_attrs.tsv" "$scratch/locus_attrs.tsv"
cp "$dir/../example_config.yml" "$scratch/config.yml"
cd "$scratch"
"$dir/microsat" "config.yml"
