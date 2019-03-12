


library(gemino)


rg <- create_genome(10, 100e3, 100, pi_tcag = c(0, 0, 0, 1))

var <- create_variants(rg, "theta", c(theta = 0.01, n_vars = 3),
                       make_mevo(rg, sub = list(model = "JC69", lambda = 0.4)))

t0 <- Sys.time()
illumina(var, out_prefix = "~/Desktop/test",
         n_reads = 10e3, read_length = 100, paired = FALSE,
         frag_mean = 400, frag_sd = 100, barcodes = c("TTTT", "CCCC", "AAAA"))
print(Sys.time() - t0)

z <- readLines("~/Desktop/test_R1.fq")
mean(grepl("^TTTT", z)) * 4
mean(grepl("^CCCC", z)) * 4
mean(grepl("^AAAA", z)) * 4





