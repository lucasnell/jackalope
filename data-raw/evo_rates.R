
library(dplyr)
library(readr)
library(rvest)

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
    select(domain, species, everything(), -label)

write_csv(evo_rates, "data-raw/evo_rates.csv")
devtools::use_data(evo_rates, overwrite = TRUE)

