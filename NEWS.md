# chiimp 0.2.1

 * Minor improvements to release process ([#14]).
 * Fixed install script for Mac OS ([#13]).
 * Fixed file-saving on Windows ([#12]).
 * Fixed installation on Windows for usernames with spaces ([#11]).
 * Added automatic categorization of genotyping results for samples from known
 individuals ([#8]).
   * Added function to pair samples with known correct genotypes,
   `match_known_genotypes`.
   * Added function to categorize results of genotyping for known individuals,
   `categorize_genotype_results`.
   * Enabled categorization features in `summarize_dataset` when Name column is
   supplied in results summary data frame.

[#14]: https://github.com/ShawHahnLab/chiimp/issues/14
[#13]: https://github.com/ShawHahnLab/chiimp/issues/13
[#12]: https://github.com/ShawHahnLab/chiimp/issues/12
[#11]: https://github.com/ShawHahnLab/chiimp/issues/11
[#8]: https://github.com/ShawHahnLab/chiimp/issues/8

# chiimp 0.2.0

 * Restructured code to avoid analyzing multiplexed samples more than once ([#3]).
 * Reorganized output into per-sequence-file (full) and per-sample (filtered)
   sections ([#5]).
   * Added a new saving function, `save_seqfile_data`, to save a directory tree
     of per-sequence-file output files starting from the first shared directory
     in the input file paths.
   * Moved functionality from `analyze_sample` into a new `analyze_seqs`
     function to be used as a separate step (enabling shared processing between
     multiplexed samples in a single data file for [#3]).
   * Split data list from `analyze_dataset` output into two separate lists
     called files and samples.
 * Added sequence name matching to `analyze_dataset`.  The summary data frame
   now has Allele1Name and Allele2Name columns, and the sample data frames a
   SeqName column, matching any sequence recognized as a called allele from any
   sample in the current analysis (or a previous analysis if the allele table
   is provided).
 * Improved `histogram` function to recognize more categories of unique
   sequences including sequences called as alleles elsewhere, and return the
   counts-by-length data by category in a list.
   * Removed `histogram2` to centralize functionality in `histogram`.
 * Fixed bugs causing failure of report generation for completely blank
   dataset analysis results ([#7]).
 * Removed `summarize_sample_by_length` function.
 * Clarified behavior of `summarize_sample` functions to allow any combination
   of TRUE/FALSE values in the Ambiguous/Stutter/Artifact entries.  Previously
   only the first (highest-count) case would be flagged.
 * Added features to track and filter ambiguous sequences ([#4]).
   * Added column named Ambiguous to `analyze_sample` output to flag sequences
     with non-ACTG characters.
   * Added entry named Ambiguous to `summarize_sample` output to track
     filtering of sequences with non-ACTG characters.

[#7]: https://github.com/ShawHahnLab/chiimp/issues/7
[#5]: https://github.com/ShawHahnLab/chiimp/issues/5
[#4]: https://github.com/ShawHahnLab/chiimp/issues/4
[#3]: https://github.com/ShawHahnLab/chiimp/issues/3

# chiimp 0.1.0

 * Initial release