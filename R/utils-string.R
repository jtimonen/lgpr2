# Remove certain suffixes from strings
replicate_without_suffix <- function(strings, sfx) {
  regex <- paste0("_", sfx, "$")
  # Identify strings that end with suffix
  has_suffix <- grepl(regex, strings)

  # Remove suffix and replicate the string if it has the suffix, else
  # return original string
  ifelse(has_suffix, gsub(regex, "", strings), strings)
}

# Colorize string
colorize_string <- function(x, col) {
  if (interactive()) {
    x <- paste0(col, x, "\u001b[0m")
  }
  x
}

# Number string
number_string <- function(x) {
  col <- "\u001b[34;1m" # bold blue
  colorize_string(x, col)
}

# Highlight string
hl_string <- function(x) {
  col <- "\u001b[33m" # orange
  colorize_string(x, col)
}

# Highlight class
hl_class <- function(x) {
  col <- "\u001b[31m" # red
  colorize_string(x, col)
}

# Main class name
class_name <- function(x) {
  class(x)[1]
}

# Main class name
class_name_hl <- function(x) {
  hl_class(class_name(x))
}
