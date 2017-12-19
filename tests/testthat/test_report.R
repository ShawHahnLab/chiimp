
# test plot_alignment -----------------------------------------------------

test_that("plot_alignment plots alignments for a locus", {
  with(results_summary_data, {
    alignments <- align_alleles(results$summary)
    fp_png <- tempfile()
    png(fp_png)
    plot_alignment(alignments[["A"]])
    dev.off()
    data_png <- readBin(fp_png, n=1e6, raw())
    hash <- openssl::md5(data_png)
    known <- unlist(strsplit("62:57:60:b9:da:13:a8:24:b2:0a:ce:ab:49:1b:62:0a",
                           ":"))
    expect_true(all(hash == known))
  })
})
