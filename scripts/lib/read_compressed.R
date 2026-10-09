read_compressed_table <- function(path, ...) {
  if (!grepl("[.]gz$", path, ignore.case = TRUE)) {
    return(data.table::fread(path, ...))
  }

  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(data.table::fread(path, ...))
  }

  gzip <- Sys.which("gzip")
  if (!nzchar(gzip)) {
    stop("Reading .gz files requires either the R.utils package or the gzip command.")
  }

  command <- paste(shQuote(gzip), "-dc", shQuote(path))
  data.table::fread(cmd = command, ...)
}
