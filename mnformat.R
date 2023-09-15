#' Write matrix to a file in binary format.
mnbin.write <- function(data) {
  con <- file("sample.mn2", open="wb")
  writeBin(ncol(data), con, size=4, endian='big')
  for (i in 1:ncol(data)) {
    writeBin(as.integer(round(1)), con, size=4, endian='big')
  }
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      writeBin(as.integer(round(data[i, j])), con, size=4, endian='big')
    }
  }
  close(con)
}
