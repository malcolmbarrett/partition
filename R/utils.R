#  variables used in various NSE calls
utils::globalVariables(
  c(
    "params"
  ))

corr <- function(x, y = NULL, spearman = FALSE) {
  if (is.null(y)) {
      dim_names <- names(x)
      if (spearman) {
        x <- apply_rank(as.matrix(x))
      }
      
      correlation <- corr_c_mat(as.matrix(x))
      
      attr(correlation, "dimnames") <- list(dim_names, dim_names)
  } else {
      if (spearman) {
        x <- rank_c(x)
        y <- rank_c(y)
      }
    
      correlation <- corr_c_2vec(x, y)
    }
  
  correlation
}

quietly <- function(.f) {
  function(...) capture_output(.f(...))
}

capture_output <- function(code) {
  warnings <- character()
  wHandler <- function(w) {
    warnings <<- c(warnings, w$message)
    invokeRestart("muffleWarning")
  }

  messages <- character()
  mHandler <- function(m) {
    messages <<- c(messages, m$message)
    invokeRestart("muffleMessage")
  }

  temp <- file()
  sink(temp)
  on.exit({
    sink()
    close(temp)
  })

  result <- withCallingHandlers(
    code,
    warning = wHandler,
    message = mHandler
  )

  output <- paste0(readLines(temp, warn = FALSE), collapse = "\n")

  list(
    result = result,
    output = output,
    warnings = warnings,
    messages = messages
  )
}

quiet_kmeans <- function(data, clusters) {
  rslts <- quietly(kmeans_c)(as.matrix(data), clusters)$result
  rslts$cluster <- as.numeric(rslts$cluster) + 1
  rslts
}

`%nin%` <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}
