# Autoformat a 'Stan' code string
autoformat_stancode <- function(code) {
  tryCatch(
    {
      file <- cmdstanr::write_stan_file(code)
      model <- cmdstanr::cmdstan_model(file, compile = FALSE)
      res <- processx::run(
        file.path(cmdstanr::cmdstan_path(), "bin", "stanc"),
        args = c(model$stan_file(), "--auto-format")
      )
      return(res$stdout)
    },
    error = function(e) {
      stop(e$stderr)
    }
  )
}

# Read stan code from inst/functions dir
read_stan_function <- function(name) {
  checkmate::assert_character(name, min.chars = 1)
  fname <- file.path("functions", paste0(name, ".stan"))
  filepath <- system.file(fname, package = "lgpr2")
  code <- read_file_lines(filepath)
  paste0("  // --------------- ", name, " ---------------\n", code, "\n")
}

# Read lines from file
read_file_lines <- function(file) {
  a <- readLines(file)
  paste(a, collapse = "\n")
}
