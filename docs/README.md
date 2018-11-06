# Worked Examples

These examples assume familiarity with the basics in the
[user guide](../GUIDE.pdf) and a working CHIIMP installation.

We have a zip file with all the required inputs to CHIIMP (using one dataset
from the publication as an example) available for download here:

https://upenn.box.com/chiimp-demo

## Sample Details via Spreadsheet

The simplest and most flexible method to define the samples to be analyzed is
to create a spreadsheet that explicitly defines the sample, replicate, and
locus identifiers for each filename.  Download the [chiimp-demo zip file] and
extract its contents.  Inside there are four files ending in .csv.  Open
`samples.csv` in a spreadsheet editor like Microsoft Excel.

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
[user guide].

[chiimp-demo zip file]: https://upenn.box.com/chiimp-demo
[user guide]: ../GUIDE.pdf
