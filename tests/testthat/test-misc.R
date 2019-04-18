
context("Testing miscellaneous operations")

# library(jackalope)
# library(testthat)


test_that("jackalope:::using_openmp() works", {
    expect_is(jackalope:::using_openmp(), "logical")
    expect_length(jackalope:::using_openmp(), 1L)

})

