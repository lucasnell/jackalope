

# library(jackalope)
# library(testthat)

context("Testing FASTA file input/output")

dir <- tempdir(check = TRUE)



chroms <- jackalope:::rando_chroms(10, 100)

ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))





# ----------*
# Errors ----
# ----------*

test_that("Read/writing FASTA files produces errors when nonsense is input", {

    fa_fn <- sprintf("%s/%s.fa", dir, "test")

    expect_error(write_fasta("ref", fa_fn),
                 regexp = "argument `obj` must be a \"ref_genome\" or \"haplotypes\"")

    expect_error(write_fasta(ref, 3),
                 regexp = "argument `out_prefix` must be a single string")

    expect_error(write_fasta(ref, fa_fn, compress = "3"),
                 regexp = paste("argument `compress` must be a single logical",
                                "or integer from 1 to 9"))

    expect_error(write_fasta(ref, fa_fn, text_width = "y"),
                 regexp = "argument `text_width` must be a single integer >= 1")


    expect_error(read_fasta(c()),
                 regexp = "argument `fasta_files` must be a character vector")

    expect_error(read_fasta(fa_fn, NA),
                 regexp = paste("argument `fai_files` must be NULL or a",
                                "character vector of the same length as the",
                                "`fasta_files` argument"))

    expect_error(read_fasta(rep(fa_fn, 2), "NA"),
                 regexp = paste("argument `fai_files` must be NULL or a",
                                "character vector of the same length as the",
                                "`fasta_files` argument"))

    expect_error(read_fasta(fa_fn, cut_names = "yeah"),
                 regexp = "argument `cut_names` must be a single logical")

})




# ================================================================================`
# ================================================================================`

# >>> Non-indexed ----

# ================================================================================`
# ================================================================================`


# ----------*
# Single ----
# ----------*


test_that("Read/writing single non-indexed FASTA files works with uncompressed output", {

    fa_fn <- sprintf("%s/%s", dir, "test")

    write_fasta(ref, fa_fn, compress = FALSE, overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa", dir, "test")
    new_ref <- read_fasta(fa_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing single non-indexed FASTA files works with gzipped output", {

    fa_fn <- sprintf("%s/%s", dir, "test")

    write_fasta(ref, fa_fn, compress = TRUE, comp_method = "gzip", overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa.gz", dir, "test")
    new_ref <- read_fasta(fa_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing single non-indexed FASTA files works with bgzipped output", {

    fa_fn <- sprintf("%s/%s", dir, "test")

    write_fasta(ref, fa_fn, compress = TRUE, comp_method = "bgzip", overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa.gz", dir, "test")
    new_ref <- read_fasta(fa_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})




# ----------*
# Multiple ----
# ----------*


ref1 <- ref_genome$new(jackalope:::make_ref_genome(chroms[1:5]))
ref2 <- ref_genome$new(jackalope:::make_ref_genome(chroms[6:10]))


test_that("Read/writing multiple non-indexed FASTA files works with uncompressed output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = FALSE, overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = FALSE, overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing multiple non-indexed FASTA files works with gzipped output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = TRUE, comp_method = "gzip", overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = TRUE, comp_method = "gzip", overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa.gz", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})

test_that("Read/writing multiple non-indexed FASTA files works with bgzipped output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = TRUE, comp_method = "bgzip", overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = TRUE, comp_method = "bgzip", overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa.gz", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})








# ================================================================================`
# ================================================================================`

# >>> Indexed ----

# ================================================================================`
# ================================================================================`




# ----------*
# Single ----
# ----------*


# Making my own fasta index file:
fa_index_df <-
    data.frame(name = ref$chrom_names(),
               nbases = ref$sizes(),
               byte_index = cumsum(c(nchar(ref$chrom_names()[1]) + 2,
                                     nchar(ref$chrom_names()[-1]) + 2 +
                                         utils::head(ref$sizes(), -1) +
                                         ceiling(utils::head(ref$sizes(), -1) /
                                                     formals(write_fasta)$text_width))),
               bases_perline = formals(write_fasta)[["text_width"]],
               bytes_perline = formals(write_fasta)[["text_width"]] + 1)


utils::write.table(fa_index_df, sprintf("%s/%s.fa.fai", dir, "test"), quote = FALSE,
                   sep = "\t", row.names = FALSE, col.names = FALSE)


test_that("Read/writing single indexed FASTA files works with uncompressed output", {

    fa_fn <- sprintf("%s/%s", dir, "test")
    fai_fn <- sprintf("%s/%s.fa.fai", dir, "test")

    write_fasta(ref, fa_fn, compress = FALSE, overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa", dir, "test")
    new_ref <- read_fasta(fa_fn, fai_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing single indexed FASTA files works with gzipped output", {

    fa_fn <- sprintf("%s/%s", dir, "test")
    fai_fn <- sprintf("%s/%s.fa.fai", dir, "test")

    write_fasta(ref, fa_fn, compress = TRUE, comp_method = "gzip", overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa.gz", dir, "test")
    new_ref <- read_fasta(fa_fn, fai_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing single indexed FASTA files works with bgzipped output", {

    fa_fn <- sprintf("%s/%s", dir, "test")
    fai_fn <- sprintf("%s/%s.fa.fai", dir, "test")

    write_fasta(ref, fa_fn, compress = TRUE, comp_method = "bgzip", overwrite = TRUE)

    fa_fn <- sprintf("%s/%s.fa.gz", dir, "test")
    new_ref <- read_fasta(fa_fn, fai_fn)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})



# ----------*
# Multiple ----
# ----------*


# Making my own fasta index file:
fa_index_df1 <- fa_index_df[1:5,]
fa_index_df2 <- fa_index_df[6:10,]
fa_index_df2$byte_index <- fa_index_df2$byte_index + ((nchar(ref$chrom_names()[1]) + 2) -
                                                          fa_index_df2$byte_index[1])


ref1 <- ref_genome$new(jackalope:::make_ref_genome(chroms[1:5]))
ref2 <- ref_genome$new(jackalope:::make_ref_genome(chroms[6:10]))

utils::write.table(fa_index_df1, sprintf("%s/%s1.fa.fai", dir, "test"),
                   quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
utils::write.table(fa_index_df2, sprintf("%s/%s2.fa.fai", dir, "test"),
                   quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


test_that("Read/writing multiple indexed FASTA files works with uncompressed output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)
    fai_fns <- sprintf("%s/%s%i.fa.fai", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = FALSE, overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = FALSE, overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns, fai_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})


test_that("Read/writing multiple indexed FASTA files works with gzipped output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)
    fai_fns <- sprintf("%s/%s%i.fa.fai", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = TRUE, comp_method = "gzip", overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = TRUE, comp_method = "gzip", overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa.gz", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns, fai_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})



test_that("Read/writing multiple indexed FASTA files works with bgzipped output", {

    fa_fns <- sprintf("%s/%s%i", dir, "test", 1:2)
    fai_fns <- sprintf("%s/%s%i.fa.fai", dir, "test", 1:2)

    write_fasta(ref1, fa_fns[1], compress = TRUE, comp_method = "bgzip", overwrite = TRUE)
    write_fasta(ref2, fa_fns[2], compress = TRUE, comp_method = "bgzip", overwrite = TRUE)

    fa_fns <- sprintf("%s/%s%i.fa.gz", dir, "test", 1:2)
    new_ref <- read_fasta(fa_fns, fai_fns)

    expect_identical(ref$n_chroms(), new_ref$n_chroms())

    for (i in 1:ref$n_chroms()) {
        expect_identical(ref$chrom(i), ref$chrom(i))
    }

})






# ___ Writing haplotypes -----

haps <- create_haplotypes(ref, haps_theta(0.1, 2), sub_JC69(0.001))

test_that("Writing FASTA files with haplotypes", {

    fa_fn <- sprintf("%s/%s", dir, "test")

    write_fasta(haps, fa_fn, compress = FALSE, overwrite = TRUE)

    fa_fn <- sprintf("%s/%s__%s.fa", dir, "test", haps$hap_names())
    new_ref1 <- read_fasta(fa_fn[1])
    new_ref2 <- read_fasta(fa_fn[2])

    expect_identical(haps$n_chroms(), new_ref1$n_chroms())
    expect_identical(haps$n_chroms(), new_ref2$n_chroms())

    expect_identical(sapply(1:haps$n_chroms(), function(i) haps$chrom(1, i)),
                     sapply(1:haps$n_chroms(), function(i) new_ref1$chrom(i)))
    expect_identical(sapply(1:haps$n_chroms(), function(i) haps$chrom(2, i)),
                     sapply(1:haps$n_chroms(), function(i) new_ref2$chrom(i)))

})

