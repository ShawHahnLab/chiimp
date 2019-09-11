#!/usr/bin/env Rscript

# Spell-check the documentation files.  Note they'll have to be updated e.g.
# with devtools::document() first.

results <- devtools::spell_check()
if (length(results) > 0) {
  results
}
