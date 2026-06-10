setwd("~/Documents/R/NAO")  # adapt path
#inspect the first lines
raw <- readLines("~/downloads/nao_index-2.tim", n = 20)
cat(raw, sep = "\n")

fn <- "nao_index-2.tim"

# Read everything, then locate the header row
all_lines <- readLines(paste0("~/downloads/",fn))

# Find the first line that looks like "YEAR MONTH INDEX"
hdr_row <- grep("YEAR *MONTH *INDEX", all_lines, ignore.case = TRUE)

if (length(hdr_row) == 0) {
  stop("No line matching YEAR MONTH INDEX found; adjust pattern.")
}

# Read starting from that header row; this becomes the header
nao <- read.table(
  paste0("~/downloads/",fn),
  header = FALSE,
  skip = hdr_row,
  stringsAsFactors = FALSE
)

# Set correct column names
names(nao) <- c("YEAR", "MONTH", "INDEX")

library(dplyr)

nao <- nao |>
  mutate(
    YEAR  = as.integer(YEAR),
    MONTH = as.integer(MONTH),
    date  = as.Date(sprintf("%d-%02d-01", YEAR, MONTH))
  ) |>
  arrange(date)

head(nao)
tail(nao)
all_lines[1:8]
