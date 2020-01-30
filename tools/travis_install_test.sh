#!/usr/bin/env bash

# This is a helper script for Travis testing.
# Assumes working directory is the project directory.

# make nonzero exit codes trigger failure
set -e
# print commands as they are executed
set -x

# (These scripts pause for input hence the pipe from the yes command.)
if   [[ $TRAVIS_OS_NAME == "linux"   ]]; then
	mkdir -p ~/Desktop
	yes | ./install_linux.sh
	# Desktop icon should be a regular file
	test -f ~/Desktop/CHIIMP.desktop
	# TODO can we test the icon on the command line?
	# Maybe xdg-open, but this bug has been open since forever:
	# https://bugs.launchpad.net/ubuntu/+source/glib2.0/+bug/378783
elif [[ $TRAVIS_OS_NAME == "osx"     ]]; then
	mkdir -p ~/Desktop
	yes | ./install_mac.command
	# Desktop icon should be a symbolic link to an existent directory
	test -d ~/Desktop/CHIIMP
	test -h ~/Desktop/CHIIMP
	# Desktop icon should run chiimp.command when a config file is dragged
	# onto it
	# First, set up inputs in a temp location.
	# TODO This is largely a repeat of what's in demo.sh
	cd "$(dirname $BASH_SOURCE)"
	dir=$(pwd -P)
	inst="../inst"
	scratch=$(mktemp -d)
	R --vanilla -q -e "devtools::load_all('..', quiet=T); test_data\$write_seqs(test_data\$seqs, '$scratch/str-dataset', 'Replicate1-Sample%s-%s.fasta')" > /dev/null
	cp "$dir/$inst/example_locus_attrs.csv" "$scratch/locus_attrs.csv"
	cp "$dir/$inst/example_config.yml" "$scratch/config.yml"
	cd "$scratch"
	# Simulate drag-and-drop onto the icon, and specify that it should exit
	# automatically when finished so Travis can continue
	CHIIMP_AUTOCLOSE=yes open -W -a ~/Desktop/CHIIMP config.yml
	# TODO why do we need to do this?  open -W should wait for the command to complete.
	timeout=60
	while [[ $timeout -gt 0 ]]; do
		test -e str-results/summary.csv && break || true
		echo "waiting for str-results/summary.csv... (${timeout}s)"
		timeout=$((timeout - 5))
		sleep 5
	done
	# In case we timed out run one final check (and fail)
	# Check that the results are there
	test -e str-results/summary.csv
elif [[ $TRAVIS_OS_NAME == "windows" ]]; then
	# (If and when enabled in the Travis config.)
	yes | ./install_windows.cmd
else
	exit 1
fi
