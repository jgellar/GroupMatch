################################################################################
# Utility functions used throughout the package
################################################################################

### toZ: turn various objects into a logical vector indicating treatment
### accepts a variety of inputs and keeps names/rownames if present

setGeneric("toZ", function(x) {
  if (is.data.frame(x) | is.matrix(x)) {
    if (dim(x)[2] > 1) {
      stop("Treatment indicators must be vectors or single column
      matrices/data.frames")
    }

    nms <- rownames(x)
    x <- x[,1]
    names(x) <- nms
  }

  nms <- names(x)
  tmp <- standardGeneric("toZ")
  names(tmp) <- nms
  return(tmp)
})

setMethod("toZ", "logical", function(x) {
  u <- unique(x)
  if (!(TRUE %in% u)) {
    stop("There must be at least one treatment unit.")
  }
  if (!(FALSE %in% u)) {
    stop("There must be at least one control unit.")
  }
  x
})

setMethod("toZ", "numeric", function(x) {
  u <- unique(x)
  if (any(!(u[!is.na(u)] %in% 0:1))) {
    stop("Numeric treatment indicators can only take on values 1 (treatment) and 0 (control).")
  }
  if (!(1 %in% u)) {
    stop("There must be at least one treatment unit.")
  }
  if (!(0 %in% u)) {
    stop("There must be at least one control unit.")
  }
  as.logical(x)
})

setMethod("toZ", "character", function(x) {
  stop(paste("Character treatment indicators no longer supported.\n",
             "Convert into a numeric or logical vector."))
})

setMethod("toZ", "factor", function(x) {
  stop(paste("Factor treatment indicators no longer supported.\n",
             "Convert into a numeric or logical vector."))
})

#' (Internal) Remove the call before digesting a distance so things
#' like omitting caliper and calling caliper=NULL give the same digest
#' #
#' @param dist Distance object to hash. Must be one of
#'   \code{InfinitySparseMatrix}, \code{BlockedInfinitySparseMatrix},
#'   \code{DenseMatrix}, \code{matrix}, or \code{distmatch.dlist}.
#' @return Hash on the distance object with a null \code{call}
#' @import digest
dist_digest <- function(dist) {
  if (class(dist)[1] %in% c("InfinitySparseMatrix", "BlockedInfinitySparseMatrix", "optmatch.dlist", "DenseMatrix", "matrix")) {
    csave <- attr(dist, "call")
    attr(dist, "call") <- NULL
    out <- digest::digest(dist)
    attr(dist, "call") <- csave
    return(out)
  }
  stop("Must pass distance object")
}

#' (Internal) If the x argument does not exist for
#' match_on, fullmatch, or pairmatch, use this function
#' to print a helpful message.
#'
#' @param x_str x as a string of code, usually deparse(substitute(x))
#' @param data_str data arg as string of code
#' @param ... will look for 'z = <stuff>' in the extra args of caller
#' @return string a helpful error message
#' @author Josh Buckner
missing_x_msg <- function(x_str, data_str, ...) {
  if(data_str == "NULL")
    data_str <- "<data argument>"

  # use a way of finding z that doesn't require z
  # to exist as an actual R obj
  z <- ""
  extra_args_str <- deparse(substitute(list(...)))
  z_regex <- "z = (\\S+)[,\\)]"
  z_search <- regexpr(z_regex, extra_args_str, perl=TRUE)
  if(z_search != -1) {
    z_match <- regmatches(extra_args_str, z_search)[1]
    z <- sub(z_regex, "\\1", z_match, perl=TRUE)
  }

  msg_tail <- if(z != "")
                paste("or ", z, "~", x_str, sep="")
              else
                ""

  paste("Can't find",
        paste(x_str, ".", sep=""),
        "If it lives within the data frame provided",
        "as the data argument, try",
        paste(data_str, "$", x_str, sep=""),
        msg_tail)
}

#' Convert GroupMatch match matrix to matching weights attached to original dataframe
#' @param mat GroupMatch match matrix
#' @param dat Dataset GroupMatch was performed on 
#' @return Dataset with matching weights
matrix_to_weights <- function(mat, dat){
  # convert matrix to dataframe and get weights for each treated unit
  match_df <- data.frame(mat) %>% mutate(wt = 1/rowSums(!is.na(.)))
  options(warn = -1) # suppress warnings temporarily
  match_df <- melt(match_df, id.vars = ncol(match_df))
  options(warn = getOption("warn"))
  # get total weight for each control unit
  match_df <- match_df %>% group_by(value) %>% 
    dplyr::summarize(match_wt = sum(wt))
  # add ID column to dat
  dat$uniqueID <- as.character(1:nrow(dat))
  # join weights to the dataset
  dat <- dat %>% left_join(match_df, by = c("uniqueID" = "value")) %>%
    mutate(match_wt = ifelse(treat == 1, 1, match_wt),
           match_wt = ifelse(is.na(match_wt) == T, 0, match_wt))
  
  return(dat)
}


#' Convert GroupMatch match matrix to optmatch setID format
#' @param match GroupMatch output
#' @return Optmatch setID vector
matrix_to_set <- function(match){
  if(length(unique(as.vector(match$cells))) == nrow(match$cells)*ncol(match$cells)){
    unit_names <- names(attr(match, "group")[[1]])
    output <- rep(0, length(unit_names))
    names(output) <- unit_names
    set_name <- seq(1, 2, length.out = nrow(match$cells))
    for(i in 1:nrow(match$cells)){
      output[which(names(output) %in% match$cells[i,] | names(output) == rownames(match$cells)[i])] <- set_name[i]
    }
  } else {
    output <- "Function only works for matching without replacement!"
  }
  return(output)
}
