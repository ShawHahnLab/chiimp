#!/usr/bin/env bash

# Install CHIIMP on Linux.

# A base R install is assumed to already be present, but all dependencies
# should be installed automatically here.

rscript=$(which Rscript)
pkgdir=$(readlink -f $(dirname $BASH_SOURCE))
"$rscript" --vanilla "$pkgdir/inst/installer/install.R"
read -n 1 -s -p "Press any key to continue . . ."; echo # see "pause" in cmd.exe
