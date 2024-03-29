library(DT)
library(tidyverse)

read.abricate <- function(x){
    t <- data.frame(read_delim(x, delim = "\t", col_names = TRUE, col_types = "ffnnfffffnnffff"))
    colnames(t) <- c("isolate","sequence", "start", "end",
                     "strand", "gene", "coverage", "coverage_map", "gaps", "percent_coverage",
                     "percent_identity", "db", "db_accession", "product", "resistance")
    t$resistance <- gsub(";", "; ", t$resistance)
    t$product <- gsub(":", ": ", t$product)
    return(t)
}

# parse args from snakemake rule
# all messages in the log
abricate_res_dir <- snakemake@params[[1]]
abricate_summary <- snakemake@output[[1]]
log <- file(snakemake@log[[1]], open="wt")
sink(log)

# craft the table
abricate.results.files <- list.files(path = abricate_res_dir, pattern = "*_*.out", full.names = TRUE)
#print(abricate.results.files)
abricate.results <- data.frame(do.call("rbind", lapply(abricate.results.files, read.abricate))) %>% mutate(resistance=tolower(resistance))

abricate.results %>%
    DT::datatable(filter = list(position="top",
                  clear=TRUE,
                  plain=FALSE),
    extensions = c("Buttons",
                   "ColReorder",
                   "FixedHeader",
                   "Scroller"),
    options = list(
        autoWidth = TRUE,
        dom = "Blfrtip",
        colReorder = TRUE,
        fixedHeader = FALSE,
        pageLength = 50,
        lengthMenu = c(5, 10, 25, 50, 100, 200, 500, 1000),
        #scrollY = 450,
        #scrollCollapse = TRUE,
        buttons = list(
            # list(extend = "colvis", columns = 1:ncol(.)),
            c('copy', 'csv', 'excel'))
    ),
    rownames = FALSE) %>%
	htmlwidgets::saveWidget(abricate_summary, selfcontained=TRUE)
