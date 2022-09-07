#!/usr/bin/env bash
# This is a helper script for CircleCI testing.
# Assumes working directory is the project directory.

set -e # make nonzero exit codes trigger failure
set -x # print commands as they are executed

# (These scripts pause for input hence the pipe from the yes command.)
mkdir -p ~/Desktop
yes | ./install_linux.sh
# Desktop icon should be a regular file
test -f ~/Desktop/CHIIMP.desktop
# TODO can we test the icon on the command line?
# Maybe xdg-open, but this bug has been open since forever:
# https://bugs.launchpad.net/ubuntu/+source/glib2.0/+bug/378783
