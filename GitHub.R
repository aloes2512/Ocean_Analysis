# Find the exact web address of your GitHub repository
url <- system("git config --get remote.origin.url", intern = TRUE)

# Convert git format to web URL if needed and open in your default browser
url <- gsub("git@github.com:", "https://github.com/", url)
url <- gsub("\\.git$", "", url)
browseURL(url)
