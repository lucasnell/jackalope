


# library(jackalope)
# library(testthat)

context("Very basic tests of sequencing output")

# For more in-depth tests, see the diagnostics folder


dir <- tempdir()

ref <- create_genome(5, 100)
tr <- ape::rcoal(4)
vars <- create_variants(ref, method = "phy", method_info = tr,
                        mevo_obj = create_mevo(ref, sub = list(model = "JC69",
                                                               lambda = 0.1)))


# ================================================================================`
# ================================================================================`

# >> Illumina -----

# ================================================================================`
# ================================================================================`

# ------*
# __Reference -----
# ------*

# single ----

test_that("no weirdness with Illumina single-end reads on ref. genome", {

    illumina(ref, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = FALSE,
             frag_mean = 400, frag_sd = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})


# paired ----

test_that("no weirdness with Illumina paired-end reads on ref. genome", {

    illumina(ref, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = TRUE,
             frag_mean = 400, frag_sd = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))
    expect_true(sprintf("%s_R2.fq", "test") %in% list.files(dir))

    fasta1 <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))
    fasta2 <- readLines(sprintf("%s/%s_R2.fq", dir, "test"))

    expect_length(fasta1, 200L)
    expect_length(fasta2, 200L)
    expect_true(all(grepl("^@", fasta1[seq(1, 200, 4)])))
    expect_true(all(grepl("^@", fasta2[seq(1, 200, 4)])))
    expect_identical(fasta1[seq(3, 200, 4)], rep("+", 50))
    expect_identical(fasta2[seq(3, 200, 4)], rep("+", 50))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))
    file.remove(sprintf("%s/%s_R2.fq", dir, "test"))

})


illumina(vars, out_prefix = sprintf("%s/%s", dir, "test"),
         n_reads = 100, read_length = 100, paired = FALSE,
         frag_mean = 400, frag_sd = 100)



# ------*
#  __Variants -----
# ------*



# single ----

test_that("no weirdness with Illumina single-end reads on variants", {

    illumina(vars, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = FALSE,
             frag_mean = 400, frag_sd = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})


# paired ----

test_that("no weirdness with Illumina paired-end reads on variants", {

    illumina(vars, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = TRUE,
             frag_mean = 400, frag_sd = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))
    expect_true(sprintf("%s_R2.fq", "test") %in% list.files(dir))

    fasta1 <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))
    fasta2 <- readLines(sprintf("%s/%s_R2.fq", dir, "test"))

    expect_length(fasta1, 200L)
    expect_length(fasta2, 200L)
    expect_true(all(grepl("^@", fasta1[seq(1, 200, 4)])))
    expect_true(all(grepl("^@", fasta2[seq(1, 200, 4)])))
    expect_identical(fasta1[seq(3, 200, 4)], rep("+", 50))
    expect_identical(fasta2[seq(3, 200, 4)], rep("+", 50))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))
    file.remove(sprintf("%s/%s_R2.fq", dir, "test"))

})




# ================================================================================`
# ================================================================================`

# >> PacBio -----

# ================================================================================`
# ================================================================================`



test_that("no weirdness with PacBio reads on ref. genome", {

    pacbio(ref, out_prefix = sprintf("%s/%s", dir, "test"),
           n_reads = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})



test_that("no weirdness with PacBio reads on variants", {

    pacbio(vars, out_prefix = sprintf("%s/%s", dir, "test"),
           n_reads = 100)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})

