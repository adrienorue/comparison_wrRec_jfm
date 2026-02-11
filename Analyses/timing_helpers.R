format_elapsed <- function(seconds) {
  minutes <- floor(seconds / 60)
  secs <- seconds - minutes * 60
  sprintf("%d:%05.2f", minutes, secs)
}

time_expr <- function(expr) {
  start <- Sys.time()
  value <- eval.parent(substitute(expr))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  list(value = value, elapsed = elapsed)
}

print_timing <- function(label, seconds) {
  cat(sprintf("[Duration] %s: %.2f s (%s)\n", label, seconds, format_elapsed(seconds)))
}
