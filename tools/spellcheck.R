#!/usr/bin/env Rscript

# Spell-check the documentation files.  Note they'll have to be updated e.g.
# with devtools::document() first.

ignore <- read.table("tools/wordlist.txt",
                     header = FALSE,
                     stringsAsFactors = FALSE)[, 1]
results <- devtools::spell_check(ignore = ignore)
if (length(results) > 0) {
  results
}
