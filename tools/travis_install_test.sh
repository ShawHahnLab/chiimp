#!/usr/bin/env bash

# This is a helper script for Travis testing.
# Assumes working directory is the project directory.

# (These scripts pause for input hence the pipe from the yes command.)
if   [[ $TRAVIS_OS_NAME -eq "linux"   ]]; then
	yes | ./install_linux.sh
elif [[ $TRAVIS_OS_NAME -eq "osx"     ]]; then
	yes | ./install_mac.command
elif [[ $TRAVIS_OS_NAME -eq "windows" ]]; then
	# (If and when enabled in the Travis config.)
	yes | ./install_windows.cmd
else
	exit 1
fi
