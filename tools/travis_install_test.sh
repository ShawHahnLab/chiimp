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
elif [[ $TRAVIS_OS_NAME == "osx"     ]]; then
	yes | ./install_mac.command
	# Desktop icon should be a symbolic link to an existent file
	test -f ~/Desktop/CHIIMP
	test -h ~/Desktop/CHIIMP
elif [[ $TRAVIS_OS_NAME == "windows" ]]; then
	# (If and when enabled in the Travis config.)
	yes | ./install_windows.cmd
else
	exit 1
fi
