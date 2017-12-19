
# test plot_alignment -----------------------------------------------------

test_that("plot_alignment plots alignments for a locus", {
  with(results_summary_data, {
    alignments <- align_alleles(results$summary)
    fp_png <- tempfile()
    png(fp_png, width = 480, height = 480, antialias = "none")
    plot_alignment(alignments[["A"]])
    dev.off()
    data_png <- readBin(fp_png, n=1e6, raw())
    hash <- openssl::md5(data_png)
    known <- unlist(strsplit("88:89:3e:cc:32:69:0b:9a:34:26:45:fe:36:c1:dc:9c",
                           ":"))
    expect_true(all(hash == known))
  })
})
