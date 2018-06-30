
library(dplyr)
library(readr)
library(rvest)


# Species that had pi_s instead of theta_s for estimate of population mutation rate
pis_s_spp <- c("Agrobacterium tumefaciens", "Pseudomonas aeruginosa", "Mus musculus")

evo_rates <- read_html("http://www.g3journal.org/content/6/8/2583.figures-only") %>%
    html_node("table") %>%
    html_table() %>%
    setNames(c("species", "label", "Ge", "Gc_Gnc",
               "indels", "subs", "pop_mutation_rate", "Ne")) %>%
    slice(c(-1,-9)) %>%
    mutate_at(vars(Ge:Ne),
              function(x) {
                  gsub("[a-zA-Z]|[^[:^punct:].-]", "", x, perl = TRUE) %>%
                      as.numeric()
              }) %>%
    mutate(domain = c(rep("Bacteria", 7), rep("Eukarya", 8))) %>%
    mutate(pi_s = ifelse(species %in% pis_s_spp, pop_mutation_rate, NA),
           theta_s = ifelse(species %in% pis_s_spp, NA, pop_mutation_rate)) %>%
    select(domain, species, everything(), -label, -pop_mutation_rate)



write_csv(evo_rates, "data-raw/evo_rates.csv")
devtools::use_data(evo_rates, overwrite = TRUE)

