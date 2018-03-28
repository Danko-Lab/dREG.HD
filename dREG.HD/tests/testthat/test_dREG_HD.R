context("dREG_HD")

## Needed to load snow options.
library(snow)

## Test fixtures.
bundled_file <- function(file_name) {
    system.file("extdata", file_name, package = "dREG.HD")
}
bed_path <- bundled_file("k562.chr21.predictions.bed")
bigwig_plus <- bundled_file("K562.chr21.plus.bw")
bigwig_minus <- bundled_file("K562.chr21.minus.bw")
load(bundled_file("dREG_HD.model.rdata"))

## Tests.
test_that("dREG_HD runs null value", {
    ret <- dREG_HD(bed_path, bigwig_plus, bigwig_minus, model, ncores = 12,
                   use_rgtsvm = FALSE)
    expect_null(ret)
})
