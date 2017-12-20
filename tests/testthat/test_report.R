
# test plot_alignment -----------------------------------------------------

test_that("plot_alignment plots alignments for a locus", {
  with(results_summary_data, {
    alignments <- align_alleles(results$summary)
    fp_png <- tempfile()
    bmp(fp_png,
        width = 100, height = 100, units = "px", pointsize = 2, bg="white",
        res = NA, type = "cairo", antialias = "none")
    plot_alignment(alignments[["A"]])
    dev.off()
    data_png <- readBin(fp_png, n=1e6, raw())
    hash <- openssl::md5(data_png)
    #known_txt <- "22:5b:4a:e9:dd:c7:65:fe:78:30:37:ec:83:d5:b5:ac" # png travis
    #known_txt <- "88:89:3e:cc:32:69:0b:9a:34:26:45:fe:36:c1:dc:9c" # png
    #known_txt <- "20:d6:97:61:f8:fd:26:dc:f4:d3:35:16:03:4b:59:b2" # bmp travis
    known_txt <- "26:6e:41:90:5d:51:00:4c:f9:9c:2d:b9:c2:cd:64:e0" # bmp
    known <- unlist(strsplit(known_txt, ":"))

    cat(paste0("\n\n\n", hash, "\n\n\n"), file = "/dev/stderr")

    stderr <- file('binary.log', 'w')
    for (i in seq(1, length(data_png), 16)) {
      cat(data_png[i:min(length(data_png),i+15)], "\n", file = binarylog)
    }
    flush(binarylog)
    close(binarylog)

    expect_true(all(hash == known))
  })
})
