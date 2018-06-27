
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(rvest)



# Altered version of rvest::html_table that doesn't remove first row:
html_table2 <- function (x, header = NA, trim = TRUE, fill = FALSE, dec = ".") {
    stopifnot(html_name(x) == "table")
    rows <- html_nodes(x, "tr")
    n <- length(rows)
    cells <- lapply(rows, "html_nodes", xpath = ".//td|.//th")
    ncols <- lapply(cells, html_attr, "colspan", default = "1")
    ncols <- lapply(ncols, as.integer)
    nrows <- lapply(cells, html_attr, "rowspan", default = "1")
    nrows <- lapply(nrows, as.integer)
    p <- unique(vapply(ncols, sum, integer(1)))
    maxp <- max(p)
    if (length(p) > 1 & maxp * n != sum(unlist(nrows)) & maxp *
        n != sum(unlist(ncols))) {
        if (!fill) {
            stop("Table has inconsistent number of columns. ",
                 "Do you want fill = TRUE?", call. = FALSE)
        }
    }
    values <- lapply(cells, html_text, trim = trim)
    out <- matrix(NA_character_, nrow = n, ncol = maxp)
    for (i in seq_len(n)) {
        row <- values[[i]]
        ncol <- ncols[[i]]
        col <- 1
        for (j in seq_len(length(ncol))) {
            out[i, col:(col + ncol[j] - 1)] <- row[[j]]
            col <- col + ncol[j]
        }
    }
    # Not sure why, but this part removes the first row:
    # for (i in seq_len(maxp)) {
    #     for (j in seq_len(n)) {
    #         rowspan <- nrows[[j]][i]
    #         colspan <- ncols[[j]][i]
    #         if (!is.na(rowspan) & (rowspan > 1)) {
    #             if (!is.na(colspan) & (colspan > 1)) {
    #                 nrows[[j]] <- c(utils::head(nrows[[j]], i),
    #                                 rep(rowspan, colspan - 1), utils::tail(nrows[[j]],
    #                                                       length(rowspan) - (i + 1)))
    #                 rowspan <- nrows[[j]][i]
    #             }
    #             for (k in seq_len(rowspan - 1)) {
    #                 l <- utils::head(out[j + k, ], i - 1)
    #                 r <- utils::tail(out[j + k, ], maxp - i + 1)
    #                 out[j + k, ] <- utils::head(c(l, out[j, i],
    #                                               r), maxp)
    #             }
    #         }
    #     }
    # }
    if (is.na(header)) {
        header <- all(html_name(cells[[1]]) == "th")
    }
    if (header) {
        col_names <- out[1, , drop = FALSE]
        out <- out[-1, , drop = FALSE]
    }
    else {
        col_names <- paste0("X", seq_len(ncol(out)))
    }
    df <- lapply(seq_len(maxp), function(i) {
        utils::type.convert(out[, i], as.is = TRUE, dec = dec)
    })
    names(df) <- col_names
    class(df) <- "data.frame"
    attr(df, "row.names") <- .set_row_names(length(df[[1]]))
    if (length(unique(col_names)) < length(col_names)) {
        warning("At least two columns have the same name")
    }
    df
}



binding_sites <- paste0("https://www.neb.com/tools-and-resources/selection-charts/",
                        "alphabetized-list-of-recognition-specificities") %>%
    read_html() %>%
    html_node("table") %>%
    html_table2() %>%
    setNames(c("sequence", "enzyme")) %>%
    tbl_df() %>%
    mutate(enzyme = strsplit(gsub("®", "", enzyme), " ")) %>%
    unnest() %>%
    mutate(sequence = ifelse(grepl("BstZ17I", x = enzyme),
                             sprintf("%.3s/%s", sequence, substr(sequence,4,6)),
                             sequence))




process_isos <- function(isos) {
    strsplit(isos, ", ") %>%
        map(function(x) {
            xx <- keep(x, ~ !grepl("\\^$", x = .x))
            return(xx)
        })
}


other_isos <- "https://www.neb.com/tools-and-resources/selection-charts/isoschizomers" %>%
    read_html() %>%
    html_node("table") %>%
    html_table2() %>%
    tbl_df() %>%
    select(sequence = `Sequence`, enzyme = `Enzyme`, isos = `Other Isoschizomers`) %>%
    filter(enzyme != "", sequence != "", isos != "") %>%
    mutate(isos = process_isos(isos),
           enzyme = gsub(" x$", "", enzyme)) %>%
    filter(map_lgl(isos, ~ length(.x) > 0))



binding_sites <- full_join(binding_sites, other_isos,
                           by = c("enzyme", "sequence")) %>%
    mutate(isos = map(isos, ~ {if (length(.x) == 0) character(0) else .x})) %>%
    group_by(enzyme) %>%
    summarize(sequence = sequence[nchar(sequence) == max(nchar(sequence))][1],
              isos = isos[sapply(isos, length) == max(sapply(isos, length))]) %>%
    ungroup() %>%
    group_by(sequence) %>%
    summarize(enzyme = list(unique(c(enzyme, isos, recursive = TRUE)))) %>%
    unnest() %>%
    select(enzyme, sequence) %>%
    as.data.frame()


write_csv(binding_sites, "data-raw/binding_sites.csv")
devtools::use_data(binding_sites, overwrite = TRUE)

