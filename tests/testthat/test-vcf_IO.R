


# library(jackalope)
# library(testthat)

# SETUP -----

context("Testing VCF file input/output")

options(stringsAsFactors = FALSE)

dir <- tempdir(check = TRUE)

chroms <- rep("TCAGTCAGTC", 2)
ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))
haps <- haplotypes$new(jackalope:::make_hap_set(ref$ptr(), 4), ref$ptr())

# First chromosome combines deletions and substitutions
{
    chrom <- 1

    haps$add_sub(1, chrom, 6, "T")
    haps$add_sub(2, chrom, 6, "A")
    haps$add_sub(3, chrom, 6, "A")
    haps$add_sub(3, chrom, 7, "T")
    haps$add_sub(4, chrom, 6, "G")
    haps$add_sub(4, chrom, 8, "T")

    haps$add_del(1, chrom, 7, 1)
    haps$add_del(2, chrom, 7, 2)
    haps$add_del(4, chrom, 9, 1)

    haps$add_del(1, chrom, 1, 3)
    haps$add_del(2, chrom, 2, 3)
    haps$add_del(3, chrom, 1, 3)
    haps$add_del(4, chrom, 3, 2)

}
# Second chromosome combines insertions, deletions, and substitutions
{
    chrom <- 2

    haps$add_del(1, chrom, 9, 1)

    haps$add_ins(1, chrom, 8, "A")

    haps$add_sub(1, chrom, 6, "A")
    haps$add_sub(2, chrom, 6, "A")
    haps$add_sub(3, chrom, 6, "T")
    haps$add_del(4, chrom, 6, 1)

    haps$add_ins(1, chrom, 5, "TT")
    haps$add_ins(2, chrom, 5, "TT")
    haps$add_ins(3, chrom, 5, "T")
    haps$add_ins(4, chrom, 5, "C")

    haps$add_sub(4, chrom, 3, "T")

    haps$add_ins(2, chrom, 2, "AG")
    haps$add_del(3, chrom, 2, 2)
    haps$add_ins(4, chrom, 2, "AG")

    haps$add_del(1, chrom, 1, 1)

}


# From these mutations, I know what the ref and alt strings should be, as well as the
# genotypes for each haplotype:
vcf_info <-
    rbind(data.frame(chrom = ref$chrom_names()[1],
                     pos = c(1, 6),
                     ref = c("TCAG", "CAGT"),
                     alt = c("G,T,TC", "TGT,AT,ATGT,GAT"),
                     gt1 = c(1, 1),
                     gt2 = c(2, 2),
                     gt3 = c(1, 3),
                     gt4 = c(3, 4)),
          data.frame(chrom = ref$chrom_names()[2],
                     pos = c(1, 5, 8),
                     ref = c("TCA", "TC", "GT"),
                     alt = c("CA,TCAGA,T,TCAGT", "TTTA,TTT", "GA"),
                     gt1 = c(1, 1, 1),
                     gt2 = c(2, 1, 0),
                     gt3 = c(3, 2, 0),
                     gt4 = c(4, 0, 0)))



# ===============================================================`
# ===============================================================`

#               WRITING -----

# ===============================================================`
# ===============================================================`


# ------------------------*
# Haploid version -----
# (i.e., each haplotype is a separate sample)
# ------------------------*

write_vcf(haps, paste0(dir, "/test"), overwrite = TRUE)
vcf <- readLines(paste0(dir, "/test.vcf"))


test_that("VCF writing produces error when attempting to overwrite files", {
    expect_error(write_vcf(haps, paste0(dir, "/test")), "test.vcf already exists")
})





test_that("VCF file header is accurate for haploid samples", {

    header <-
        c("##fileformat=VCFv4.3",
          sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
          "##source=jackalope",
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[1], ref$sizes()[1]),
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[2], ref$sizes()[2]),
          "##phasing=full",
          paste("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of",
                "Samples With Data\">"),
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GenotypeQuality\">",
          paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                 paste(haps$hap_names(), collapse = "\t")))

    expect_identical(vcf[1:length(header)], header)

})

# Remove header lines
vcf <- vcf[!grepl("^#", vcf)]

test_that("VCF file data lines are accurate for haploid samples", {

    # String to insert VCF info into:
    test_str <- paste0(c("%s\t%i\t.\t%s\t%s\t441453\tPASS\tNS=4\tGT:GQ",
                         rep("%i:441453", 4)), collapse = "\t")

    data_lines <- lapply(1:2,
                         function(i) {
                             sprintf(
                                 rep(test_str, sum(vcf_info$chrom == ref$chrom_names()[i])),
                                 ref$chrom_names()[i],
                                 vcf_info$pos[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$ref[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$alt[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt1[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt2[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt3[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt4[vcf_info$chrom == ref$chrom_names()[i]])
                         })

    expect_identical(data_lines[[1]], vcf[grepl(paste0("^", ref$chrom_names()[1]), vcf)])
    expect_identical(data_lines[[2]], vcf[grepl(paste0("^", ref$chrom_names()[2]), vcf)])

})





# Haploid compressed -----


write_vcf(haps, paste0(dir, "/test"), compress = TRUE, overwrite = TRUE)
vcf_gz <- gzfile(paste0(dir, "/test.vcf.gz"), "rt")
vcf <- readLines(vcf_gz)
close(vcf_gz)

test_that("VCF file header is accurate for haploid samples", {

    header <-
        c("##fileformat=VCFv4.3",
          sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
          "##source=jackalope",
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[1], ref$sizes()[1]),
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[2], ref$sizes()[2]),
          "##phasing=full",
          paste("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of",
                "Samples With Data\">"),
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GenotypeQuality\">",
          paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                 paste(haps$hap_names(), collapse = "\t")))

    expect_identical(vcf[1:length(header)], header)

})


# Remove header lines
vcf <- vcf[!grepl("^#", vcf)]

test_that("VCF file data lines are accurate for haploid samples", {

    # String to insert VCF info into:
    test_str <- paste0(c("%s\t%i\t.\t%s\t%s\t441453\tPASS\tNS=4\tGT:GQ",
                         rep("%i:441453", 4)), collapse = "\t")

    data_lines <- lapply(1:2,
                         function(i) {
                             sprintf(
                                 rep(test_str, sum(vcf_info$chrom == ref$chrom_names()[i])),
                                 ref$chrom_names()[i],
                                 vcf_info$pos[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$ref[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$alt[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt1[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt2[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt3[vcf_info$chrom == ref$chrom_names()[i]],
                                 vcf_info$gt4[vcf_info$chrom == ref$chrom_names()[i]])
                         })

    expect_identical(data_lines[[1]], vcf[grepl(paste0("^", ref$chrom_names()[1]), vcf)])
    expect_identical(data_lines[[2]], vcf[grepl(paste0("^", ref$chrom_names()[2]), vcf)])

})



# ------------------------*
# Diploid version -----
# (i.e., 2 haplotypes represent one sample)
# ------------------------*



sample_mat <- matrix(1:haps$n_haps(), ncol = 2, byrow = TRUE)

write_vcf(haps, paste0(dir, "/test"), sample_matrix = sample_mat, overwrite = TRUE)
vcf <- readLines(paste0(dir, "/test.vcf"))



test_that("VCF writing produces error with nonsense sample matrix", {

    expect_error(write_vcf(haps, paste0(dir, "/test"), sample_matrix = sample_mat * -1,
                           overwrite = TRUE),
                 "there are values < 1.")

    sample_mat2 <- sample_mat
    sample_mat2[nrow(sample_mat2),ncol(sample_mat2)] <- haps$n_haps() + 1

    expect_error(write_vcf(haps, paste0(dir, "/test"), sample_matrix = sample_mat2,
                           overwrite = TRUE),
                 "there are values > the number of haplotypes")

    sample_mat2 <- sample_mat
    sample_mat2[nrow(sample_mat2),ncol(sample_mat2)] <- 1

    expect_error(write_vcf(haps, paste0(dir, "/test"), sample_matrix = sample_mat2,
                           overwrite = TRUE),
                 "contained duplicates")

})



test_that("VCF file header is accurate for diploid samples", {

    sample_names <- apply(sample_mat, 1, function(x) paste(haps$hap_names()[x],
                                                           collapse = "__"))

    header <-
        c("##fileformat=VCFv4.3",
          sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
          "##source=jackalope",
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[1], ref$sizes()[1]),
          sprintf("##contig=<ID=%s,length=%i>", ref$chrom_names()[2], ref$sizes()[2]),
          "##phasing=full",
          paste("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of",
                "Samples With Data\">"),
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GenotypeQuality\">",
          paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                 paste(sample_names, collapse = "\t")))

    expect_identical(vcf[1:length(header)], header)

})



test_that("VCF file data lines are accurate for diploid samples", {

    # String to insert VCF info into:
    test_str <- paste0(c("%s\t%i\t.\t%s\t%s\t441453\tPASS\tNS=%i\tGT:GQ",
                         rep("%s:441453", nrow(sample_mat))), collapse = "\t")

    data_lines <- lapply(1:2,
                         function(i) {
                             vcf_info_ <- vcf_info[vcf_info$chrom == ref$chrom_names()[i],]
                             n_lines <- nrow(vcf_info_)
                             gt_mat <- vcf_info_[, paste0("gt", 1:4)]
                             gt_mat <- as.matrix(gt_mat)
                             gt_strings <- lapply(
                                 1:nrow(sample_mat),
                                 function(j) {
                                     sapply(1:n_lines,
                                           function(k) paste(gt_mat[k,][sample_mat[j,]],
                                                             collapse = "|"))
                                 })
                             arg_list <- c(list(rep(test_str, n_lines),
                                                rep(ref$chrom_names()[i], n_lines),
                                                vcf_info_$pos,
                                                vcf_info_$ref,
                                                vcf_info_$alt,
                                                rep(nrow(sample_mat), n_lines)),
                                           gt_strings)

                             do.call(sprintf, arg_list)
                         })

    expect_identical(data_lines[[1]], vcf[grepl(paste0("^", ref$chrom_names()[1]), vcf)])
    expect_identical(data_lines[[2]], vcf[grepl(paste0("^", ref$chrom_names()[2]), vcf)])

})







# ===============================================================`
# ===============================================================`

#           READING -----

# ===============================================================`
# ===============================================================`







test_that("reading haploid haplotype info from VCF produces proper output", {

    write_vcf(haps, out_prefix = sprintf("%s/%s", dir, "test"), overwrite = TRUE)

    haps2 <- create_haplotypes(ref, haps_info = haps_vcf(sprintf("%s/%s.vcf", dir, "test")))


    expect_identical(haps$n_haps(), haps2$n_haps())

    for (i in 1:haps$n_haps()) {
        expect_identical(sapply(1:ref$n_chroms(), function(j) haps$chrom(i, j)),
                         sapply(1:ref$n_chroms(), function(j) haps2$chrom(i, j)))
    }


})


test_that("reading diploid haplotype info from VCF produces proper output", {

    sample_mat <- matrix(1:4, 2, 2, byrow = TRUE)

    write_vcf(haps, out_prefix = sprintf("%s/%s", dir, "test"),
              sample_matrix = sample_mat, overwrite = TRUE)

    haps2 <- create_haplotypes(ref, haps_info = haps_vcf(sprintf("%s/%s.vcf", dir, "test")))

    expect_identical(haps$n_haps(), haps2$n_haps())

    for (i in 1:haps$n_haps()) {
        expect_identical(sapply(1:ref$n_chroms(), function(j) haps$chrom(i, j)),
                         sapply(1:ref$n_chroms(), function(j) haps2$chrom(i, j)))
    }

})



test_that("reading haplotype info from VCF produces proper output when chromosomes mixed", {

    # --------------*
    # Haploid version:
    # --------------*

    write_vcf(haps, out_prefix = sprintf("%s/%s", dir, "test"), overwrite = TRUE)

    vcf_fn <- sprintf("%s/%s.vcf", dir, "test")

    # Mix up the first two chromosomes in the VCF file:
    vcf <- readLines(vcf_fn)
    c0 <- vcf[grepl(paste0("^", haps$chrom_names()[1]), vcf)]
    c1 <- vcf[grepl(paste0("^", haps$chrom_names()[2]), vcf)]
    vcf <- c(vcf[!grepl(paste0("^", haps$chrom_names()[1]), vcf) &
                     !grepl(paste0("^", haps$chrom_names()[2]), vcf)],
             c1, c0)
    writeLines(paste(vcf, collapse = "\n"), vcf_fn)

    # Now create haplotypes and check output
    haps2 <- create_haplotypes(ref, haps_info = haps_vcf(vcf_fn))

    expect_identical(haps$n_haps(), haps2$n_haps())

    for (i in 1:haps$n_haps()) {
        expect_identical(sapply(1:ref$n_chroms(), function(j) haps$chrom(i, j)),
                         sapply(1:ref$n_chroms(), function(j) haps2$chrom(i, j)))
    }


    # --------------*
    # Diploid version:
    # --------------*

    sample_mat <- matrix(1:4, 2, 2, byrow = TRUE)

    write_vcf(haps, out_prefix = sprintf("%s/%s", dir, "test"),
              sample_matrix = sample_mat, overwrite = TRUE)

    # Mix up the first two chromosomes in the VCF file:
    vcf <- readLines(vcf_fn)
    c0 <- vcf[grepl(paste0("^", haps$chrom_names()[1]), vcf)]
    c1 <- vcf[grepl(paste0("^", haps$chrom_names()[2]), vcf)]
    vcf <- c(vcf[!grepl(paste0("^", haps$chrom_names()[1]), vcf) &
                     !grepl(paste0("^", haps$chrom_names()[2]), vcf)],
             c1, c0)
    writeLines(paste(vcf, collapse = "\n"), vcf_fn)

    haps2 <- create_haplotypes(ref, haps_info = haps_vcf(vcf_fn))

    expect_identical(haps$n_haps(), haps2$n_haps())

    for (i in 1:haps$n_haps()) {
        expect_identical(sapply(1:ref$n_chroms(), function(j) haps$chrom(i, j)),
                         sapply(1:ref$n_chroms(), function(j) haps2$chrom(i, j)))
    }


})





test_that("out-of-order VCF file returns error", {

    write_vcf(haps, out_prefix = sprintf("%s/%s", dir, "test"), overwrite = TRUE)

    vcf_fn <- sprintf("%s/%s.vcf", dir, "test")

    # Reverse lines for first two chromosomes in the VCF file:
    vcf <- readLines(vcf_fn)
    c0 <- vcf[grepl(paste0("^", haps$chrom_names()[1]), vcf)]
    c1 <- vcf[grepl(paste0("^", haps$chrom_names()[2]), vcf)]
    c0 <- rev(c0)
    c1 <- rev(c1)
    vcf <- c(vcf[!grepl(paste0("^", haps$chrom_names()[1]), vcf) &
                     !grepl(paste0("^", haps$chrom_names()[2]), vcf)],
             c1, c0)
    writeLines(paste(vcf, collapse = "\n"), vcf_fn)

    # Now attempt to create haplotypes and check output
    expect_error(create_haplotypes(ref, haps_info = haps_vcf(vcf_fn)),
                 regexp = paste("Positions are sorted numerically, in",
                                "increasing order, within each reference",
                                "sequence CHROM"))

})


