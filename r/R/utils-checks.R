# Argument checks shared by every exported function. Each check returns the
# value in the type the engine expects, or stops with a message that names
# the argument in backticks and says what would be accepted (call. = FALSE,
# since the caller's name adds nothing for the user). See "Missing values
# and invalid input" in ?net_frame for the rules these implement.

# A single TRUE or FALSE.
check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

# For arguments whose `NULL` default means "decided elsewhere".
check_flag_or_null <- function(value, name) {
  if (!is.null(value) && (!is.logical(value) || length(value) != 1L || is.na(value))) {
    stop(sprintf("`%s` must be TRUE, FALSE, or NULL.", name), call. = FALSE)
  }
  value
}

# match.arg() for arguments that are not formals with a vector of choices.
check_choice <- function(value, choices, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop(sprintf("`%s` must be one of %s.", name,
                 paste0("\"", choices, "\"", collapse = ", ")), call. = FALSE)
  }
  value
}

# A single string, or NULL when `null_ok`.
check_string <- function(value, name, null_ok = FALSE) {
  if (is.null(value) && null_ok) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be a single string%s.", name, if (null_ok) " or NULL" else ""),
         call. = FALSE)
  }
  value
}

# A single file path.
check_path <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("`%s` must be a single file path.", name), call. = FALSE)
  }
  value
}

# A single whole number of at least `min`, returned as an integer.
check_count <- function(value, name, min = 0L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < min ||
      value != floor(value) || value > .Machine$integer.max) {
    what <- if (min <= 0L) {
      "a single non-negative whole number"
    } else {
      sprintf("a single whole number of at least %d", as.integer(min))
    }
    stop(sprintf("`%s` must be %s.", name, what), call. = FALSE)
  }
  as.integer(value)
}

# A single number within [lower, upper] (or the open interval when
# `open`), returned as a double. `finite = FALSE` admits Inf at an
# unbounded end.
check_number <- function(value, name, lower = -Inf, upper = Inf, open = FALSE,
                         finite = TRUE) {
  ok <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    (!finite || is.finite(value)) &&
    (if (open) value > lower && value < upper else value >= lower && value <= upper)
  if (!ok) {
    stop(sprintf("`%s` must be %s.", name, number_text(lower, upper, open, finite)),
         call. = FALSE)
  }
  as.double(value)
}

number_text <- function(lower, upper, open, finite) {
  if (lower == 0 && upper == Inf) {
    return(if (open) "a single positive number" else "a single non-negative number")
  }
  if (lower == -Inf && upper == Inf) {
    return(if (finite) "a single finite number" else "a single number")
  }
  if (is.finite(lower) && is.finite(upper)) {
    return(sprintf("a single number %s %s and %s", if (open) "strictly between" else "between",
                   format(lower), format(upper)))
  }
  if (is.finite(lower)) {
    return(sprintf("a single number %s %s", if (open) "greater than" else "of at least",
                   format(lower)))
  }
  sprintf("a single number %s %s", if (open) "less than" else "of at most", format(upper))
}

# A single probability.
check_probability <- function(p, name) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
    stop(sprintf("`%s` must be a single probability between 0 and 1.", name), call. = FALSE)
  }
  as.double(p)
}

# A vector of non-negative whole numbers (degrees, block sizes).
check_counts <- function(d, name) {
  if (!is.numeric(d) || anyNA(d) || any(d < 0) || any(d != floor(d)) ||
      any(d > .Machine$integer.max)) {
    stop(sprintf("`%s` must be a vector of non-negative whole numbers.", name),
         call. = FALSE)
  }
  as.integer(d)
}

# The convergence controls of iterative measures.
check_iteration <- function(tol, max_iter) {
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !isTRUE(max_iter >= 1) ||
      max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive whole number.", call. = FALSE)
  }
  invisible(NULL)
}

# A seed for the engine's random streams: a whole number between 0 and
# 2^53, as a double. `NULL` draws one from R's generator (two draws give 52
# random bits, so set.seed() controls the result). Every function that
# takes `seed` uses this, so a seed means the same thing everywhere.
check_seed <- function(seed) {
  if (is.null(seed)) {
    return(floor(stats::runif(1L) * 2^26) * 2^26 + floor(stats::runif(1L) * 2^26))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
      seed > 2^53 || seed != floor(seed)) {
    stop("`seed` must be a whole number between 0 and 2^53.", call. = FALSE)
  }
  as.double(seed)
}

# Resolution parameter of modularity-based methods.
check_resolution <- function(resolution) {
  if (!is.numeric(resolution) || length(resolution) != 1L || is.na(resolution) ||
      !is.finite(resolution) || resolution < 0) {
    stop("`resolution` must be a single finite, non-negative number.", call. = FALSE)
  }
  as.double(resolution)
}

# A snapping or matching distance.
check_tolerance <- function(tolerance) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance < 0) {
    stop("`tolerance` must be a single non-negative distance.", call. = FALSE)
  }
  invisible(tolerance)
}

# ---- Tables and columns ------------------------------------------------------

# A column of `df` given by name or position, as its name. `role` is the
# argument being resolved and `what` the table, for the message.
resolve_column <- function(df, col, role, what) {
  if (is.numeric(col) && length(col) == 1L && !is.na(col) && col >= 1 &&
      col <= ncol(df) && col == floor(col)) {
    return(names(df)[[col]])
  }
  if (is.character(col) && length(col) == 1L && !is.na(col) && col %in% names(df)) {
    return(col)
  }
  shown <- if (ncol(df)) {
    sprintf(" It has %s: %s.", plural(ncol(df), "column"),
            paste0("`", utils::head(names(df), 10L), "`", collapse = ", "))
  } else {
    " It has no columns."
  }
  stop(sprintf(paste0("Column for `%s` not found in `%s`; give a column name or ",
                      "position.%s"), role, what, shown), call. = FALSE)
}

# Checks an identifier column and returns it in the form used for
# matching: factors as character, 64-bit integers from bit64 as character
# (their storage is not meaningful as doubles), other classed vectors such
# as dates unchanged, and plain vectors without attributes.
check_ids <- function(v, col, table, what) {
  if (!is.atomic(v) || is.matrix(v)) {
    stop(sprintf("%s (column `%s` of `%s`) must be an atomic vector, not %s.",
                 what, col, table, describe_class(v)), call. = FALSE)
  }
  v <- if (inherits(v, "integer64")) {
    as.character(v)
  } else if (is.factor(v)) {
    as.vector(v)
  } else if (is.object(v)) {
    v
  } else {
    as.vector(v)
  }
  bad <- which(is.na(v))
  if (length(bad)) {
    stop(sprintf(paste0("%s must not contain missing values (NA); %s of `%s` %s one ",
                        "in column `%s`. Remove those rows or fill in the identifiers."),
                 what, rows_text(bad), table, if (length(bad) == 1L) "has" else "have",
                 col), call. = FALSE)
  }
  if (is.character(v)) {
    bad <- which(!nzchar(v))
    if (length(bad)) {
      stop(sprintf(paste0("%s must not be empty strings; %s of `%s` %s one in column ",
                          "`%s`. Remove those rows or fill in the identifiers."),
                   what, rows_text(bad), table, if (length(bad) == 1L) "has" else "have",
                   col), call. = FALSE)
    }
  }
  v
}

is_classed <- function(v) is.object(v) && !is.factor(v)

# ---- Text for messages -------------------------------------------------------

plural <- function(k, what) sprintf("%d %s%s", k, what, if (k == 1L) "" else "s")

# "row 3", "rows 3 and 7", or "rows 3, 7, 9, 12, 15, and 35 more". `total`
# is the number of rows affected when `idx` lists only the first of them.
rows_text <- function(idx, total = length(idx)) {
  k <- length(idx)
  shown <- format(utils::head(idx, 5L), scientific = FALSE, trim = TRUE)
  if (k == 1L && total == 1L) return(paste("row", shown))
  if (k == 2L && total == 2L) return(paste("rows", shown[1], "and", shown[2]))
  if (k <= 5L && total == k) {
    return(paste0("rows ", paste(shown[-k], collapse = ", "), ", and ", shown[k]))
  }
  sprintf("rows %s, and %d more", paste(shown, collapse = ", "), total - length(shown))
}

describe_class <- function(v) {
  if (is.null(v)) return("NULL")
  sprintf("an object of class `%s`", class(v)[[1]])
}

# Stops unless a suggested package is installed. `what` names the feature
# that needs it ("net_to_sf()", "Writing a GeoPackage").
require_pkg <- function(pkg, what = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- if (is.null(what)) {
      sprintf("The %s package is required", pkg)
    } else {
      sprintf("%s needs the %s package", what, pkg)
    }
    stop(sprintf("%s; install it with install.packages(\"%s\").", msg, pkg), call. = FALSE)
  }
  invisible(TRUE)
}
