#!/usr/bin/env bash

# Install CHIIMP on Mac OS.

# A base R install is assumed to already be present, but all dependencies
# should be installed automatically here.

rscript=$(which Rscript)
pkgdir=$(dirname $BASH_SOURCE)
"$rscript" --vanilla "$pkgdir/tools/install.R"
read -n 1 -s -p "Press any key to continue . . ."; echo # see "pause" in cmd.exe
