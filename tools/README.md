# Package Development Tools

These are files related to CHIIMP package development and initial installation.
These aren't relevant within an installation or for regular users.

 * `chiimp.AppleScript`: Source code for the Desktop icon wrapper application
   on Mac OS.  See comments in the script for usage.
 * `conda_environment.yml`: Conda environment definition file.
 * `conda_install.sh`: Conda install wrapper.
 * `install.R`: Small wrapper script to call installation code from within the
   package.
 * `travis_install_test.sh` Script to manage testing the various installers
   under Travis.
 * `check_lint.R`: Checks for R code quality
 * `check_spelling.R` Checks for spelling in documentation.  Add extra words to
   `inst/WORDLIST`.  Note that it's clever about code encapsulated in backticks
   and won't flag it.
 * `prep_release.sh` New version release preparation script

## Manual Testing Steps

On each OS test that:

 1. From-scratch installation (no existing Desktop icon or user R library)
    works as expected
 1. Re-installation works as expected (package and Desktop icon are replaced)
 1. Desktop icon shows usage message when double-clicked
 1. Desktop icon runs the demo when the config is dragged onto it
 1. The report is generated as expected
