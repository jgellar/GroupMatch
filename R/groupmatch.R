#' (Internal) Sets up option to try recovery in \code{groupmatch}.
#'
#' @return NULL
setTryRecovery <- function() {
  options("groupmatch_try_recovery" = TRUE)
}


#' Optimal full matching with control groups
#'
#' This is an adaption of \code{\link[optmatch]{fullmatch}} to allow for
#' restrictions when control observations are "grouped". The motivating use
#' case is when there are multiple observations of control data for each
#' control subject. In this case, the grouping variable is the subject. We
#' may want to place restrictions, for example that only one observation of
#' a subject can be matched, or in the case of one:many matching, a given
#' control subject can only be matched to a given treated subject once.
#'
#'
#' @param x Any valid input to \code{match_on}. \code{groupmatch} will use
#' \code{x} and any optional arguments to generate a distance before performing
#' the matching.
#'
#' If \code{x} is a numeric vector, there must also be passed a vector \code{z}
#' indicating grouping. Both vectors must be named.
#'
#' Alternatively, a precomputed distance may be entered. A matrix of
#' non-negative discrepancies, each indicating the permissibility and
#' desirability of matching the unit corresponding to its row (a 'treatment') to
#' the unit corresponding to its column (a 'control'); or, better, a distance
#' specification as produced by \code{\link{match_on}}.
#'
#' @param group Grouping variable for control group.
#'
#' @param allow_duplicates Whether we allow duplicates (more to come....)
#'
#' @param min.controls The minimum ratio of controls to treatments that is to
#' be permitted within a matched set: should be non-negative and finite.  If
#' \code{min.controls} is not a whole number, the reciprocal of a whole number,
#' or zero, then it is rounded \emph{down} to the nearest whole number or
#' reciprocal of a whole number.
#'
#' When matching within subclasses (such as those created by
#' \code{\link{exactMatch}}), \code{min.controls} may be a named numeric vector
#' separately specifying the minimum permissible ratio of controls to treatments
#' for each subclass.  The names of this vector should include names of all
#' subproblems \code{distance}.
#'
#' @param max.controls The maximum ratio of controls to treatments that is
#' to be permitted within a matched set: should be positive and numeric.
#' If \code{max.controls} is not a whole number, the reciprocal of a
#' whole number, or \code{Inf}, then it is rounded \emph{up} to the
#' nearest whole number or reciprocal of a whole number.
#'
#' When matching within subclasses (such as those created by
#' \code{\link{exactMatch}}), \code{max.controls} may be a named numeric vector
#' separately specifying the maximum permissible ratio of controls to treatments
#' in each subclass.
#'
#' @param omit.fraction Optionally, specify what fraction of controls or treated
#' subjects are to be rejected.  If \code{omit.fraction} is a positive fraction
#' less than one, then \code{groupmatch} leaves up to that fraction of the control
#' reservoir unmatched.  If \code{omit.fraction} is a negative number greater
#' than -1, then \code{groupmatch} leaves up to |\code{omit.fraction}| of the
#' treated group unmatched.  Positive values are only accepted if
#' \code{max.controls} >= 1; negative values, only if \code{min.controls} <= 1.
#' If neither \code{omit.fraction} or \code{mean.controls} are specified, then
#' only those treated and control subjects without permissible matches among the
#' control and treated subjects, respectively, are omitted.
#'
#' When matching within subclasses (such as those created by
#' \code{\link{exactMatch}}), \code{omit.fraction} specifies the fraction of
#' controls to be rejected in each subproblem, a parameter that can be made to
#' differ by subclass by setting \code{omit.fraction} equal to a named numeric
#' vector of fractions.
#'
#' At most one of \code{mean.controls} and \code{omit.fraction} can be non-\code{NULL}.
#'
#' @param mean.controls Optionally, specify the average number of controls per
#' treatment to be matched. Must be no less than than \code{min.controls} and no
#' greater than the either \code{max.controls} or the ratio of total number of
#' controls versus total number of treated. Some controls will likely not be
#' matched to ensure meeting this value. If neither \code{omit.fraction} or
#' \code{mean.controls} are specified, then only those treated and control
#' subjects without permissible matches among the control and treated subjects,
#' respectively, are omitted.
#'
#' When matching within subclasses (such as those created by
#' \code{\link{exactMatch}}), \code{mean.controls} specifies the average number of
#' controls per treatment per subproblem, a parameter that can be made to
#' differ by subclass by setting \code{mean.controls} equal to a named numeric
#' vector.
#'
#' At most one of \code{mean.controls} and \code{omit.fraction} can be non-\code{NULL}.
#'
#' @param tol Because of internal rounding, \code{groupmatch} may
#' solve a slightly different matching problem than the one
#' specified, in which the match generated by
#' \code{groupmatch} may not coincide with an optimal solution of
#' the specified problem.  \code{tol} times the number of subjects
#' to be matched specifies the extent to
#' which \code{groupmatch}'s output is permitted to differ from an
#' optimal solution to the original problem, as measured by the
#' sum of discrepancies for all treatments and controls placed
#' into the same matched sets.
#'
#' @param data Optional \code{data.frame} or \code{vector} to use to get order
#' of the final matching factor. If a \code{data.frame}, the \code{rownames}
#' are used. If a vector, the \code{names} are first tried, otherwise the contents
#' is considered to be a character vector of names. Useful to pass if you want to
#' combine a match (using, e.g., \code{cbind}) with the data that were used to
#' generate it (for example, in a propensity score matching).
#'
#' @param ... Additional arguments, including \code{within}, which may be passed to \code{match_on}.
#'
#' @return A \code{\link{optmatch}} object (\code{factor}) indicating matched groups.
#'
#' @references
#' Our paper?
#'
#' @example inst/examples/groupmatch.R
#' @keywords nonparametric optimize
#' @export
groupmatch <- function(x, group = NULL,
                       allow_duplicates = FALSE,
                       return_style = "vector",
                       min.controls = 0,
                       max.controls = Inf,
                       replace_value = FALSE,
                       omit.fraction = NULL,
                       mean.controls = NULL,
                       tol = .001,
                       data = NULL,
                       ...) {

  # if x does not exist then print helpful error msg
  x_str <- deparse(substitute(x))
  data_str <- deparse(substitute(data))
  tryCatch(x, error = function(e) {
    stop(missing_x_msg(x_str, data_str, ...))})

  cl <- match.call()
  if (is.null(data)) {
    if (is(x, "InfinitySparseMatrix") |
        is(x, "matrix") |
        is(x, "optmatch.dlist") )
      warning("Without 'data' argument the order of the match is not guaranteed
    to be the same as your original data.")
  }
  UseMethod("groupmatch")
}

#' @export
groupmatch.default <- function(x, group = NULL,
                               allow_duplicates = FALSE,
                               return_style = "vector",
                               min.controls = 0,
                               max.controls = Inf,
                               omit.fraction = NULL,
                               mean.controls = NULL,
                               replace_value = FALSE,
                               tol = .001,
                               data = NULL,
                               within = NULL,
                               ...) {

  if (!inherits(x, gsub("match_on.","",methods("match_on")))) {
    stop("Invalid input, must be a potential argument to match_on")
  }

  mfd <- if (!is.null(data)) {
    model.frame(data, na.action=na.pass)
  } else {
    if (inherits(x, "function")) {
      stop("A data argument must be given when passing a function")
    }
    model.frame(x, na.action=na.pass)
  }
  if (!class(mfd) == "data.frame") {
    stop("Please pass data argument")
  }
  m <- match_on(x, within=within, data=mfd, ...)
  out <- groupmatch(m, group, allow_duplicates,
                    return_style = return_style,
                    min.controls=min.controls,
                    max.controls=max.controls,
                    omit.fraction=omit.fraction,
                    mean.controls=mean.controls,
                    replace_value,
                    tol=tol,
                    data=mfd,
                    ...)
  if (!exists("cl")) cl <- match.call()
  attr(out, "call") <- cl
  out
}

#' @export
groupmatch.numeric <- function(x, group = NULL, allow_duplicates = FALSE,
                               return_style = "vector",
                               min.controls = 0,
                               max.controls = Inf,
                               omit.fraction = NULL,
                               mean.controls = NULL,
                               replace_value = FALSE,
                               tol = .001,
                               data = NULL,
                               z,
                               within = NULL,
                               ...) {

  m <- match_on(x, within=within, z=z, ...)
  out <- groupmatch(m, group, allow_duplicates,
                    return_style=return_style,
                    min.controls=min.controls,
                    max.controls=max.controls,
                    omit.fraction=omit.fraction,
                    mean.controls=mean.controls,
                    replace_value,
                    tol=tol,
                    data=data,
                    ...)
  if (!exists("cl")) cl <- match.call()
  attr(out, "call") <- cl
  out
}

#' @export
groupmatch.matrix <- function(x, group = NULL, allow_duplicates = FALSE,
                              return_style = "vector",
                              min.controls = 0,
                              max.controls = Inf,
                              omit.fraction = NULL,
                              mean.controls = NULL,
                              replace_value = FALSE,
                              tol = .001,
                              data = NULL,
                              within = NULL,
                              ...) {
  ###LNV: add a check here that max.controls is not greater than the number of groups divided by number of Ts?
  ###      Also possibly a check that min.controls >=1 (i.e. no replacement)?

  ### check of the class of group argument -- vector of group IDs (current set-up); quoted variable name, tilde variable name

  # Check return_style
  if (!return_style %in% c("vector", "matrix")) {
    stop("return style must be 'vector' or 'matrix'")
  }

  if (allow_duplicates & return_style == "vector") {
    warning("Matching with replacement requires matrix-style return object; setting return_style to 'matrix'")
    return_style <- "matrix"
  }





  ### Checking Input ###

  # Convert group to a vector of group IDs
  group <- if (is.character(group) & length(group)==1) {
    if (is.null(data)) stop("data required when group is a character")
    data[[group]]
  } else if (class(group) == "formula") {
    if (is.null(data)) stop("data required when group is a formula")
    model.frame(group, data = data, na.action = na.pass)[[1]]
  } else if (is.vector(group)) {
    group
  } else {
    stop("Unrecognized type for group")
  }

  # Convert to a factor
  group <- as.factor(group)

  # Check/assign group names
  if (!is.null(data)) {
    # Name the groups with the T/C name (rowname of data)
    names(group) <- rownames(data)
  } else if (is.null(names(group))) {
    stop("group must be named if data is not supplied")
  } else if (!(names(group) %in% c(rownames(x), colnames(x)))) {
    stop("group names do not match dimnames of distance matrix")
  }

  # Checks for group argument
  if (length(group) != sum(dim(x))) {
    stop("length of group does not match dimensions of difference matrix")
  }


  # this will throw an error if not valid
  validDistanceSpecification(x)

  # note: we might want to move these checks to validDistSpec
  dnms <- dimnames(x)
  if (is.null(dnms) | is.null(dnms[[1]]) | is.null(dnms[[2]])) {
    stop("argument \'x\' must have dimnames")
  }

  if (any(duplicated(unlist(dnms)))){
    stop("dimnames of argument \'x\' contain duplicates")
  }

  if (!is.null(within)) warning("Ignoring non-null 'within' argument.  When using 'groupmatch' with\n pre-formed distances, please combine them using '+'.")

  nmtrt <- dnms[[1]]
  nmctl <- dnms[[2]]

  # note: this next _should_ be unnecessary, the objects should do this
  # but better safe than sorry
  if (!isTRUE(all.equal(dim(x), c(length(nmtrt), length(nmctl))))) {
    stop("argument \'x\' dimensions do not match row and column names")
  }

  if (!is.numeric(min.controls)) {
    stop("argument \'min.controls\' must be numeric")
  }
  if (!is.numeric(max.controls)) {
    stop("argument \'max.controls\' must be numeric")
  }
  if (!is.null(omit.fraction)) {
    # A vector of all NA's is logical, not numeric, so the first condition is needed.
    if (all(is.na(omit.fraction))) {
      omit.fraction <- NULL
    } else if (any(abs(omit.fraction) > 1, na.rm = TRUE) | !is.numeric(omit.fraction)) {
      stop("omit.fraction must be NULL or numeric between -1 and 1")
    }
  }
  if (!is.null(mean.controls)) {
    if (all(is.na(mean.controls))) {
      mean.controls <- NULL
    } else if (any(mean.controls <= 0, na.rm = TRUE) | !is.numeric(mean.controls)) {
      stop("mean.controls must be NULL or numeric greater than 0")
    }
  }

  if (!is.null(omit.fraction) & !is.null(mean.controls)) {
    stop("omit.fraction and mean.controls cannot both be specified")
  }


  # Issue #56: Checking for sane input in data
  if (!is.null(data)) {
    if (!is.vector(data)) {
      dnames <- rownames(data)
    } else {
      dnames <- names(data)
    }
    if (any(!unlist(dimnames(x)) %in% dnames)) {
      stop("Some elements of the distance matrix are not found in the data argument.")
    }
  }

  # problems is guaranteed to be a list of DistanceSpecifictions
  # it may only have 1 entry
  problems <- findSubproblems(x)

  # the number of problems should match the argument lengths for
  # min, max, and omit

  np <- length(problems)
  if (length(min.controls) > 1 & np != length(min.controls)) {
    stop(paste("Length of \'min.controls\' arg must be same ",
               "as number of subproblems [", np, "]", sep = ""))
  }
  if (length(max.controls) > 1 & np != length(max.controls)) {
    stop(paste("Length of \'max.controls\' arg must be same ",
               "as number of subproblems [", np, "]", sep = ""))
  }
  if (!is.null(omit.fraction) & length(omit.fraction) > 1 & np !=
      length(omit.fraction)) {
    stop(paste("Length of \'omit.fraction\' arg must be same ",
               "as number of subproblems [", np, "]", sep = ""))
  }
  if (!is.null(mean.controls) & length(mean.controls) > 1 & np !=
      length(mean.controls)) {
    stop(paste("Length of \'mean.controls\' arg must be same ",
               "as number of subproblems [", np, "]", sep = ""))
  }

  # reset the arguments to be the right length if they are not
  if (length(min.controls) == 1) {
    min.controls <- rep(min.controls, np)
  }
  if (length(max.controls) == 1) {
    max.controls <- rep(max.controls, np)
  }

  if (is.null(omit.fraction)) {
    omit.fraction <- NA
  }
  if (length(omit.fraction) == 1) {
    omit.fraction <- rep(omit.fraction, np)
  }
  denom <- if (allow_duplicates) ncol(x) else length(unique(group[colnames(x)]))
  if (min.controls == max.controls & is.na(omit.fraction) & (min.controls*nrow(x) < denom)) {
    #  Feasible, fixed-ratio matching. Take a shortcut by setting the omf.
    omit.fraction <- 1 - min.controls * nrow(x)/denom
  }

  if (is.null(mean.controls)) {
    mean.controls <- NA
  }
  if (length(mean.controls) == 1) {
    mean.controls <- rep(mean.controls, np)
  }

  if (!is.list(group)) {
    group <- list(group)
  }
  if (length(group) == 1) {
    group <- rep(group, np)
  }

  if (any(mean.controls < min.controls, na.rm=TRUE)) {
    stop("mean.controls cannot be smaller than min.controls")
  }

  if (any(mean.controls > max.controls, na.rm=TRUE)) {
    stop("mean.controls cannot be larger than max.controls")
  }

  if (any(!is.na(mean.controls))) {
    if (any(mean.controls > lapply(problems, function(p) {x <- subdim(p)[[1]] ;  x[2]/x[1]}), na.rm=TRUE)) {
      stop("mean.controls cannot be larger than the ratio of number of controls to treatments")
    }
  }

  if (any(omit.fraction > 0 & max.controls <= .5, na.rm=TRUE)) {
    stop("positive \'omit.fraction\' with \'max.controls\' <= 1/2 not permitted")
  }

  if (any(omit.fraction < 0 & min.controls >= 2, na.rm=TRUE)) {
    stop("negative \'omit.fraction\' with \'min.controls\' >= 2 not permitted")
  }


  user.input.mean.controls <- FALSE

  # ?????Do we need to fix this check?????
  if (any(!is.na(mean.controls) & is.na(omit.fraction))) {
    user.input.mean.controls <- TRUE
    omit.fraction <- 1 - mapply(function(x,y) {z <- subdim(y)[[1]] ; x*z[1]/z[2]}, mean.controls, problems)
  }

  total.n <- sum(dim(x))

  TOL <- tol * total.n

  # a helper to handle a single matching problem. all args required.
  # input error checking happens in the public groupmatch function.
  .groupmatch <- function(d, mnctl, mxctl, omf, g, ad, r) {

    # if the subproblem is completely empty, short circuit
    if (length(d) == 0 || all(is.infinite(d))) {
      x <- dim(d)
      cells.a <- rep(NA, x[1])
      cells.b <- rep(NA, x[2])
      names(cells.a) <- rownames(d)
      names(cells.b) <- colnames(d)
      tmp <- list(cells = c(cells.a, cells.b), maxerr = -1)
      return(tmp)
    }

    ncol <- dim(d)[2]
    nrow <- dim(d)[1]
    ncg  <- length(unique(g[colnames(d)]))

    tol.frac <- (nrow + ncol - 2)/(total.n - 2 * np)

    # if omf is specified (i.e. not NA), see if is non-negative
    # if omf is not specified, check to see if mxctl is > .5
    if (switch(1 + is.na(omf), omf >= 0,  mxctl > .5)) {
      maxc <- min(mxctl, ncg)
      minc <- max(mnctl, 1/nrow)
      omf.calc <- omf

    } else {
      maxc <- min(1/mnctl, ncol)
      minc <- max(1/mxctl, 1/nrow)
      omf.calc <- -1 * omf
      d <- t(d)
    }

    temp <- gSubDivStrat(rownames = rownames(d),
                         colnames = colnames(d),
                         distspec = d,
                         group = g,
                         allow_duplicates = ad,
                         return_style = return_style,
                         max.cpt = maxc,
                         min.cpt = minc,
                         replace_value = r,
                         tolerance = TOL * tol.frac,
                         omit.fraction = if(!is.na(omf)) { omf.calc }) # passes NULL for NA

    return(temp)
  }

  # a second helper function, that will attempt graceful recovery in situations where the match
  # is infeasible with the given max.controls
  .groupmatch.with.recovery <- function(d.r, mnctl.r, mxctl.r, omf.r,
                                        g.r, ad.r, r.r) {
    denom <- if (ad.r) ncol(d.r) else length(unique(g.r[colnames(d.r)]))

    # if (mnctl.r == mxctl.r & is.na(omf.r) & (mnctl.r*nrow(d.r) < denom)) {
    #   #  Feasible, fixed-ratio matching. Take a shortcut by setting the omf.
    #   omf.r <- 1 - mnctl.r * nrow(d.r)/denom
    # }

    # if the subproblem isn't clearly infeasible, try to get a match
    if (mxctl.r * dim(d.r)[1] >= prod(denom, 1-omf.r, na.rm=TRUE)) {
      tmp <- .groupmatch(d.r, mnctl.r, mxctl.r, omf.r, g.r, ad.r, r.r)
      if (!all(is.na(tmp[1]$cells)) && !all(tmp[1]$cells == 'NA')) {
        # subproblem is feasible with given constraints, no need to recover
        new.omit.fraction <<- c(new.omit.fraction, omf.r)
        return(tmp)
      }
    }
    # if max.control is in [1, Inf), and we're infeasible
    if(is.finite(mxctl.r) & mxctl.r >= 1) {
      # Re-solve with no max.control
      # Since max.control is capped at number of control groups, ensure we
      # omit any excess flow first.
      ncg  <- length(unique(g.r[colnames(d.r)]))
      if(denom > nrow(d.r)*ncg){
        omf.r.excess <- (denom - nrow(d.r)*ncg)/denom
        if(!is.na(omf.r) && omf.r > omf.r.excess){
          omf.r.excess <- omf.r
        }else{
          new.omit.fraction <<- c(new.omit.fraction, omf.r.excess)
        }
      }else{
        omf.r.excess <- omf.r
      }
      tmp2 <- list(.groupmatch(d.r, mnctl.r, Inf, omf.r.excess, g.r, ad.r, r.r))
      #need to remove code that depends on the makeOptmatch command, no longer
      #compatible with our "allow duplicates" setting.
      # tmp2.optmatch <- makeOptmatch(d.r, tmp2, match.call(), data)
      # trial.ss <- stratumStructure(tmp2.optmatch)
      #treats <- as.numeric(unlist(lapply(strsplit(names(trial.ss), ":"),"[",1)))
      #ctrls <- as.numeric(unlist(lapply(strsplit(names(trial.ss), ":"),"[",2)))
      ctrls <- apply(tmp2[[1]]$cells, 1, function(x) length(na.omit(x)))
      #num.controls <- sum((pmin(ctrls, mxctl.r)*trial.ss)[treats > 0])
      num.controls <- sum(pmin(ctrls, mxctl.r))
      #if(num.controls == 0) {
      if (all(is.na(tmp2[[1]]$cells)) || all(tmp2[[1]]$cells == 'NA') ||
          length(tmp2[[1]]$cells) == 0){
        # infeasible anyways
        if (!exists("tmp")) {
          tmp <- .groupmatch(d.r, mnctl.r, mxctl.r, omf.r, g.r, ad.r, r.r)
        }
        new.omit.fraction <<- c(new.omit.fraction, omf.r)
        return(tmp)
      }
      new.omf.r <- 1 - num.controls/denom

      # feasible with the new omit fraction
      new.omit.fraction <<- c(new.omit.fraction, new.omf.r)
      return(.groupmatch(d.r, mnctl.r, mxctl.r, new.omf.r, g.r, ad.r, r.r))
    } else {
      # subproblem is infeasible, but we can't try to fix because no max.controls
      if (!exists("tmp")) {
        tmp <- .groupmatch(d.r, mnctl.r, mxctl.r, omf.r, g.r, ad.r, r.r)
      }

      new.omit.fraction <<- c(new.omit.fraction, omf.r)
      return(tmp)
    }
  }

  # In case we need to try and recover from infeasible, save the new.omit.fraction's used for output to user
  new.omit.fraction <- numeric(0)

  if (is.null(options()$groupmatch_try_recovery)) {
    warning("The flag groupmatch_try_recovery is unset, setting to TRUE")
    setTryRecovery()
  }

  if (options()$groupmatch_try_recovery) {
    solutions <- mapply(.groupmatch.with.recovery, problems,
                        min.controls, max.controls, omit.fraction,
                        group, allow_duplicates, replace_value,
                        SIMPLIFY = FALSE)
  } else {
    solutions <- mapply(.groupmatch, problems,
                        min.controls, max.controls, omit.fraction,
                        group, allow_duplicates, replace_value,
                        SIMPLIFY = FALSE)
  }

  # g for matrix output
  #mout <- makeOptmatch(x, solutions, match.call(), data)
  mout <- solutions[[1]]

  names(min.controls) <- names(problems)
  names(max.controls) <- names(problems)
  attr(mout, "min.controls") <- min.controls
  attr(mout, "max.controls") <- max.controls
  attr(mout, "group") <- group
  attr(mout, "allow_duplicates") <- allow_duplicates

  # length(new.omit.fraction) will be strictly positive if we ever entered .groupmatch.with.recovery
  if(length(new.omit.fraction) > 0) {
    out.omit.fraction <- new.omit.fraction
  } else {
    out.omit.fraction <- omit.fraction
  }
  out.mean.controls <- mapply(function(x,y) (1 - x)*y[2]/y[1], out.omit.fraction, subdim(x))

  names(out.mean.controls) <- names(problems)
  names(out.omit.fraction) <- names(problems)

  if(user.input.mean.controls) {
    attr(mout, "mean.controls") <- out.mean.controls
  } else {
    attr(mout, "omit.fraction") <- out.omit.fraction
  }

  if(length(new.omit.fraction) > 0 & !identical(new.omit.fraction, omit.fraction) & !all(is.na(new.omit.fraction))) {
    if(!any(is.na(new.omit.fraction)) & all(new.omit.fraction == 1)) {
      # If we never got a feasible subproblem
      warning("The problem appears infeasible with the given constraints.")
    } else {
      warning("The problem is infeasible with the given constraints; some units were omitted to allow a match.")
    }
  }

  # save hash of distance
  attr(mout, "hashed.distance") <- dist_digest(x)

  if (!exists("cl")) cl <- match.call()
  attr(mout, "call") <- cl
  return(mout)
}


#' @export
groupmatch.optmatch.dlist <- groupmatch.matrix
#' @export
groupmatch.InfinitySparseMatrix <- groupmatch.matrix
#' @export
groupmatch.BlockedInfinitySparseMatrix <- groupmatch.matrix
