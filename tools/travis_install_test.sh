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
	# TODO we may be able to test the drag-and-drop behavior with a command like this:
	# open -W -a ~/Desktop/CHIIMP .../config.yml
	# Currently the .app exits and leaves the Terminal window open, though,
	# so -W doesn't work as expected.  Changing the AppleScript .app export
	# options might resolve this.
elif [[ $TRAVIS_OS_NAME == "windows" ]]; then
	# (If and when enabled in the Travis config.)
	yes | ./install_windows.cmd
else
	exit 1
fi
