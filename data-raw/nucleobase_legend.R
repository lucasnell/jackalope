

library(dplyr)
library(readr)
library(tidyr)
library(rvest)


page <- read_html(paste0("https://en.wikipedia.org/wiki/List_of_restriction_enzyme",
                         "_cutting_sites:_A#Whole_list_navigation"))


nucleobase_legend <- page %>%
    html_nodes("table") %>%
    .[[1]] %>%
    html_table() %>%
    setNames(c("code", "nucleotides")) %>%
    slice(-1:-5) %>%
    mutate(nucleotides = strsplit(nucleotides, " or |, ")) %>%
    unnest()

write_csv(nucleobase_legend, "data-raw/nucleobase_legend.csv")
devtools::use_data(nucleobase_legend, overwrite = TRUE)

