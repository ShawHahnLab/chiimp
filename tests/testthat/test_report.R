
# test plot_alignment -----------------------------------------------------

test_that("plot_alignment plots alignments for a locus", {
  with(results_summary_data, {
    alignments <- align_alleles(results$summary)
    fp_img <- "alignment.svg.log"
    svg(fp_img)
    plot_alignment(alignments[["A"]])
    dev.off()
    data_img <- readBin(fp_img, n=1e6, raw())
    hash <- openssl::md5(data_img)
    known_txt <- "25:95:62:2f:b7:ee:d6:5d:41:bf:17:c1:24:16:d5:e5"
    known <- unlist(strsplit(known_txt, ":"))
    cat(paste0("\n\n\n", hash, "\n\n\n"), file = "/dev/stderr")
    expect_true(all(hash == known))
  })
})
