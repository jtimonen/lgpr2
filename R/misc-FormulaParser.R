# Define the FormulaParser class
FormulaParser <- R6::R6Class("FormulaParser",
  private = list(
    terms = NULL,
    y_name = NULL
  ),
  public = list(
    initialize = function(formula) {
      str <- as.character(formula)
      private$y_name <- str[2]
      private$terms <- create_rhs(str[3])
    },
    num_terms = function() {
      length(private$terms)
    },
    get_y_name = function() {
      private$y_name
    },
    parse_term = function(idx) {
      checkmate::assert_integerish(idx, lower = 0, upper = self$num_terms())
      t <- private$terms[idx]
      pt <- parse_formulaterm(t)
      if (pt$type == "gp") {
        check_no_hierarchy(pt, "gp")
        out <- GPTerm$new(x_name = pt$covariates[1], z_name = pt$covariates[2])
      } else if (pt$type == "offset") {
        out <- create_offsetterm(pt$covariates, pt$hierarchy)
      } else {
        stop("allowed terms are gp() or offset()!")
      }
      return(out)
    },
    parse_terms = function() {
      J <- self$num_terms()
      terms <- list()
      for (j in seq_len(J)) {
        terms[[j]] <- self$parse_term(j)
      }
      gn <- function(x) x$stanname_base()
      tn <- sapply(terms, gn)
      names(terms) <- tn
      terms
    }
  )
)

extract_function_and_argument <- function(input_string) {
  # Regular expression to match function name and arguments
  regex <- "([^\\(]+)\\((.*)\\)$"

  # Apply regular expression to input string
  matches <- regmatches(input_string, regexec(regex, input_string))

  # Extract matches, if any
  if (length(matches[[1]]) > 2) {
    func_name <- matches[[1]][2]
    args <- matches[[1]][3]
    return(list(func_name = func_name, args = args))
  } else {
    return(NULL)
  }
}

split_expression <- function(input_string) {
  # Initialize variables to keep track of parentheses depth and positions for splitting
  paren_depth <- 0
  start_pos <- 1
  parts <- list()

  # Iterate through the string to find top-level plus signs
  for (i in 1:nchar(input_string)) {
    char <- substr(input_string, i, i)
    if (char == "(") {
      paren_depth <- paren_depth + 1
    } else if (char == ")") {
      paren_depth <- paren_depth - 1
    } else if (char == "+" && paren_depth == 0) {
      # Split at top-level plus, excluding the plus sign itself
      parts <- c(parts, substr(input_string, start_pos, i - 1))
      start_pos <- i + 1
    }
  }

  # Add the last part of the string if there's any remaining
  if (start_pos <= nchar(input_string)) {
    parts <- c(parts, substr(input_string, start_pos, nchar(input_string)))
  }

  return(parts)
}


# Parse formula terms
create_rhs <- function(rhs) {
  checkmate::assert_character(rhs)
  parts <- split_expression(rhs)
  trimws(unlist(parts))
}

# Parse formula terms
parse_formulaterm <- function(t) {
  checkmate::assert_character(t)
  a <- extract_function_and_argument(t)
  type <- trimws(a$func_name)
  rem <- strsplit(a$args, split = "[|]")[[1]]
  covs <- trimws(strsplit(rem[1], ",", fixed = TRUE)[[1]])
  if (length(rem) > 2) {
    hier <- paste(rem[2:length(rem)], collapse = "|")
  } else {
    hier <- rem[2]
  }
  hier <- parse_hier_or_formula(hier)
  list(
    type = type,
    covariates = covs,
    hierarchy = hier
  )
}

# helper
parse_hier_or_formula <- function(hier) {
  # split at comma but not inside parentheses
  pattern <- ",\\s*(?![^()]*\\))"
  hier <- trimws(strsplit(hier, pattern, perl = TRUE)[[1]])
  if (length(hier) < 2) {
    if (is.na(hier)) {
      hier <- NULL
    }
  }
  hier
}

# Make sure that error is thrown if hierarchy argument is given
# but will not be used
check_no_hierarchy <- function(parsed_term, type) {
  msg <- paste0(type, " terms cannot have hierarchy!")
  if (!is.null(parsed_term$hierarchy)) {
    stop(msg)
  }
  invisible(NULL)
}
