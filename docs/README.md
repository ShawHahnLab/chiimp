# Worked Examples

These examples assume familiarity with the basics in the
[user guide](../GUIDE.pdf) and a working CHIIMP installation.

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
   `locus_attrs.csv`.  This defines which loci should have amplicons present in
   a given file.  (If amplicons for multiple loci were pooled in each file,
   there can be multiple rows per file, with only the Locus varying.)
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

**Drag and drop the file onto the CHIIMP icon on the Desktop.**

[chiimp-demo zip file]: https://upenn.box.com/chiimp-demo
[user guide]: ../GUIDE.pdf
[Notepad++]: https://notepad-plus-plus.org/
