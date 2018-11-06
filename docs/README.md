# Worked Examples

These examples assume familiarity with the basics in the
[user guide](../GUIDE.pdf) and a working CHIIMP installation.  Windows and Mac
OS are given as examples but all steps should also apply under Linux.

We have a zip file with all the required inputs to CHIIMP (using one dataset
from the publication as an example) available for download here:

[https://upenn.box.com/chiimp-demo](https://upenn.box.com/chiimp-demo)

## Sample Details via Spreadsheet

The simplest and most flexible method to define the samples to be analyzed is
to create a spreadsheet that explicitly defines the sample, replicate, and
locus identifiers for each filename.  Here we'll download existing input files
prepared in this way and run the analysis.

**Download the [chiimp-demo zip file] and extract its contents.  Open
`samples.csv` in a spreadsheet editor like Microsoft Excel.**

The columns, from left to right:

 * **Sample** defines what biological sample each file is associated with.
 * There may be multiple sequencer samples for each biological sample, so the
   next column, **Replicate**, allows an integer identifier to distinguish
   these.  Here we don't have multiple replicates so the column is blank, but
   still required.
 * **GMCode** is not a sample attribute used by CHIIMP, so the software will
   carry it along with sample attributes but will otherwise ignore it.  (This
   allows custom sample attributes in existing spreadsheets to be kept  in
   place.)  In this case it's one of our sample identifiers.
 * The **Locus** column has identifiers corresponding to those in the
   `locus_attrs.csv`.  This defines which microsatellite loci should have
   amplicons present in a given file.  (If amplicons for multiple loci were
   pooled in each file, there can be multiple rows per file, with only the
   Locus varying.)
 * The **Filename** column gives the location of each file.  If this is
   relative (like folder/file.fastq.gz) it will be interpreted relative to the
   config file (more on this below).  If absolute (like C:\folder\file.fastq.gz
   or /home/user/file.fastq.gz) it will be taken as-is.

The other three spreadsheets (`known_alleles.csv`, `known_genotypes.csv`, and
`locus_attrs.csv`) are laid out similarly and are described more in the
[user guide].  The `data/` folder contains the sequence data (adapter-trimmed
forward reads from a MiSeq run).

The configuration file, `config.yml`, is a simple text file that will tell the
software where to find these spreadsheet and data files and how to do the
analysis.  The built-in text editor in RStudio (installed during the CHIIMP
setup) can edit this format and will provide some features like color-coding
and automatic indentation for sub-sections.  Other text editors are fine too,
though Windows doesn't come with a particularly useful text editor, so it may
be helpful to install one yourself such as [Notepad++].  After a default
install you can simply right-click on a CHIIMP config file and click "Edit with
Notepad++". On Mac OS TextEdit will work if kept in plain text mode.  For now
you don't need to worry about the contents of the configuration file.

**Drag and drop `config.yml` onto the CHIIMP icon on the Desktop.**

A terminal window should appear showing the progress of the analysis.  It will
remain open after completion so you can refer to the output text if needed.
The text should look something like this:

    Loading dataset: samples.csv...                              2018-11-06 16:21:38
    Loading locus attrs: locus_attrs.csv...                      2018-11-06 16:21:38
    Analyzing samples...                                         2018-11-06 16:21:38
    Summarizing results...                                       2018-11-06 16:24:19
    Saving output files...                                       2018-11-06 16:24:25
    Creating report...                                           2018-11-06 16:25:47
    Done.                                                        2018-11-06 16:25:58
    Press any key to continue . . .

If the text shows "Done" near the bottom the software ran successfully.  A new
folder should appear alongside the example spreadsheets and data folder called
`results`.

(**If no results folder appears or there are errors, please [report an issue]
and we can help investigate.**)

**Inside the results folder, double-click `report.html`.**

The report is an HTML document summarizing the results of an analysis.  It
should open automatically in your default web browser.  The most basic summary
is the table of alleles across samples and loci at the top.  This is the core
of CHIIMP's output, since it identifies a genotype for each biological sample.
Additional sections below that describe additional aspects as detailed in the
[user guide].

**Open the `summary.csv` spreadsheet.**

The summary spreadsheet contains one row per analyzed sample per expected
locus, and many columns of output data. The first few columns are just the same
sample attributes that were given as input.  The remaining columns describe
what alleles were identified and any notable features flagged during analysis.
In this case a table of names for alleles was configured (`known_alleles.csv`)
so a short name is given in both the report and this spreadsheet for recognized
alleles.  In this spreadsheet the alleles are *not* duplicated for homozygous
cases, as the columns specifically refer to the unique sequences identified in
the data for the first and second most prominent candidate sequence.  (In the
report, candidate homozygous cases have the same allele shown twice.)

The other output files and folders provide *much* more detail, down to the
level of each unique sequence identified.  Except for the alignment and
histogram images, the data is saved in FASTA and CSV formats so it's easy to
work with from this point on. `results.rds` was defined in this configuration
file to save the full set of output data in an R-compatible binary format, so
any later analysis in R can easily load the results in one step.  Alternatively
the individual files can be loaded in Geneious, Excel, etc for post-processing.

## Sample Details via Filename Matching

*To add here: instructions for using CHIIMP without a samples.csv spreadsheet.*

[chiimp-demo zip file]: https://upenn.box.com/chiimp-demo
[user guide]: ../GUIDE.pdf
[report an issue]: https://github.com/ShawHahnLab/chiimp/issues/new
[Notepad++]: https://notepad-plus-plus.org/
