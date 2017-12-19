
# test plot_alignment -----------------------------------------------------

test_that("plot_alignment plots alignments for a locus", {
  with(results_summary_data, {
    alignments <- align_alleles(results$summary)
    fp_png <- tempfile()
    png(fp_png,
        width = 480, height = 480, res = NA, type = "cairo", antialias = "none")
    plot_alignment(alignments[["A"]])
    dev.off()
    data_png <- readBin(fp_png, n=1e6, raw())
    hash <- openssl::md5(data_png)
    #known_txt <- "22:5b:4a:e9:dd:c7:65:fe:78:30:37:ec:83:d5:b5:ac"
    known_txt <- "88:89:3e:cc:32:69:0b:9a:34:26:45:fe:36:c1:dc:9c"
    known <- unlist(strsplit(known_txt, ":"))
    cat(paste0("\n\n\n", hash, "\n\n\n"), file="/dev/stderr")
    expect_true(all(hash == known))
  })
})
