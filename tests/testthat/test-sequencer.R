


# library(jackalope)
# library(testthat)

context("Very basic tests of sequencing output")


dir <- tempdir(check = TRUE)

ref <- create_genome(5, 100)
tr <- ape::rcoal(4)
haps <- create_haplotypes(ref, haps_phylo(tr), sub = sub_JC69(0.1))


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
             overwrite = TRUE)

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
             overwrite = TRUE)

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



# Profile with no msimatches
profile_df <- expand.grid(nucleo = c("T", "C", "A", "G"),
                          pos = 0:99,
                          qual = c(255L, 1000L), stringsAsFactors = FALSE)
profile_df <- profile_df[order(profile_df$nucleo, profile_df$pos, profile_df$qual),]
write.table(profile_df, file = sprintf("%s/%s", dir, "test_prof.txt"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)



test_that("proper pairs created with Illumina paired-end reads on ref. genome", {

    # 1 chromosome of length 200
    chrom <- paste(c(rep('C', 25), rep('N', 150), rep('T', 25)), collapse = "")

    poss_pairs <- c(paste(c(rep('C', 25), rep('N', 75)), collapse = ""),
                    paste(c(rep('A', 25), rep('N', 75)), collapse = ""))

    # Make ref_genome object from a pointer to a RefGenome object based on `chrom`
    rg <- ref_genome$new(jackalope:::make_ref_genome(chrom))

    illumina(rg, out_prefix = paste0(dir, "/test"),
             n_reads = 10e3, read_length = 100,
             # Paired-end reads:
             paired = TRUE, matepair = FALSE,
             # Fragments will always be of length 200:
             frag_len_min = 200, frag_len_max = 200,
             # No sequencing errors:
             ins_prob1 = 0, del_prob1 = 0,
             ins_prob2 = 0, del_prob2 = 0,
             profile1 = paste0(dir, "/test_prof.txt"),
             profile2 = paste0(dir, "/test_prof.txt"),
             overwrite = TRUE)

    fq1 <- readLines(paste0(dir, "/test_R1.fq"))
    fq2 <- readLines(paste0(dir, "/test_R2.fq"))
    reads1 <- fq1[seq(2, length(fq1), 4)]
    reads2 <- fq2[seq(2, length(fq2), 4)]

    # Should both be true:
    expect_identical(sort(unique(reads1)), sort(poss_pairs))
    expect_identical(sort(unique(reads2)), sort(poss_pairs))

})



test_that("proper pairs created with Illumina mate-pair reads on ref. genome", {

    # 1 chromosome of length 200
    chrom <- paste(c(rep('C', 25), rep('N', 150), rep('T', 25)), collapse = "")

    poss_pairs <- c(paste(c(rep('N', 75), rep('T', 25)), collapse = ""),
                    paste(c(rep('N', 75), rep('G', 25)), collapse = ""))

    # Make ref_genome object from a pointer to a RefGenome object based on `chrom`
    rg <- ref_genome$new(jackalope:::make_ref_genome(chrom))

    illumina(rg, out_prefix = paste0(dir, "/test"),
             n_reads = 10e3, read_length = 100,
             # Mate-pair reads:
             paired = TRUE, matepair = TRUE,
             # Fragments will always be of length 200:
             frag_len_min = 200, frag_len_max = 200,
             # No sequencing errors:
             ins_prob1 = 0, del_prob1 = 0,
             ins_prob2 = 0, del_prob2 = 0,
             profile1 = paste0(dir, "/test_prof.txt"),
             profile2 = paste0(dir, "/test_prof.txt"),
             overwrite = TRUE)

    fq1 <- readLines(paste0(dir, "/test_R1.fq"))
    fq2 <- readLines(paste0(dir, "/test_R2.fq"))
    reads1 <- fq1[seq(2, length(fq1), 4)]
    reads2 <- fq2[seq(2, length(fq2), 4)]

    # Should both be true:
    expect_identical(sort(unique(reads1)), sort(poss_pairs))
    expect_identical(sort(unique(reads2)), sort(poss_pairs))

})


# ------*
#  __Variants -----
# ------*



# single ----

test_that("no weirdness with Illumina single-end reads on haplotypes", {

    illumina(haps, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = FALSE,
             overwrite = TRUE)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})


test_that("no weirdness with Illumina single-end reads on haplotypes w/ sep. files", {

    illumina(haps, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = FALSE,
             overwrite = TRUE, sep_files = TRUE)

    fns <- sprintf("%s_%s_R1.fq", "test", haps$hap_names())

    expect_true(all(fns %in% list.files(dir)))

    fns <- paste0(dir, "/", fns)

    fastas <- lapply(fns, readLines)

    expect_identical(sum(sapply(fastas, length)), 400L)
    expect_true(all(sapply(fastas,
                           function(fa) all(grepl("^@", fa[seq(1, length(fa), 4)])))))
    expect_identical(do.call(c, lapply(fastas, function(fa) fa[seq(3, length(fa), 4)])),
                     rep("+", 100))

    file.remove(fns)

})


# paired ----

test_that("no weirdness with Illumina paired-end reads on haplotypes", {

    illumina(haps, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = TRUE,
             overwrite = TRUE)

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




test_that("no weirdness with Illumina paired-end reads on haplotypes w/ sep. files", {

    illumina(haps, out_prefix = sprintf("%s/%s", dir, "test"),
             n_reads = 100, read_length = 100, paired = TRUE,
             overwrite = TRUE, sep_files = TRUE)

    fns <- lapply(haps$hap_names(),
                  function(x) sprintf("%s_%s_R%i.fq", "test", x, 1:2))

    expect_true(all(unlist(fns) %in% list.files(dir)))

    fns <- lapply(fns, function(x) paste0(dir, "/", x))

    fastas <- lapply(fns, function(x) lapply(x, readLines))

    expect_true(all(sapply(fastas, function(x) length(x[[1]]) == length(x[[2]]))))

    fastas <- unlist(fastas, recursive = FALSE)

    expect_length(unlist(fastas), 400L)
    expect_true(all(sapply(fastas,
                           function(fa) all(grepl("^@", fa[seq(1, length(fa), 4)])))))
    expect_identical(do.call(c, lapply(fastas, function(fa) fa[seq(3, length(fa), 4)])),
                     rep("+", 100))

    file.remove(unlist(fns))

})


# ================================================================================`
# ================================================================================`

# >> PacBio -----

# ================================================================================`
# ================================================================================`



test_that("no weirdness with PacBio reads on ref. genome", {

    pacbio(ref, out_prefix = sprintf("%s/%s", dir, "test"),
           n_reads = 100, overwrite = TRUE)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})



test_that("no weirdness with PacBio reads on haplotypes", {

    pacbio(haps, out_prefix = sprintf("%s/%s", dir, "test"),
           n_reads = 100, overwrite = TRUE)

    expect_true(sprintf("%s_R1.fq", "test") %in% list.files(dir))

    fasta <- readLines(sprintf("%s/%s_R1.fq", dir, "test"))

    expect_length(fasta, 400L)
    expect_true(all(grepl("^@", fasta[seq(1, 400, 4)])))
    expect_identical(fasta[seq(3, 400, 4)], rep("+", 100))

    file.remove(sprintf("%s/%s_R1.fq", dir, "test"))

})

