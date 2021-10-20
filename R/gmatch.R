gmatch <- function(distance, group, allow_duplicates = FALSE,
                   max.row.units, max.col.units, min.col.units = 1, f = 1,
                   stability.increment=1L, num_inst_use = 10, replace_value = FALSE)
  
  #NOTE: f is 1-omit.fraction, i.e. the proportion of Cs to keep.  Defaults to 1, i.e. full matching.
{
  if(!inherits(distance, "data.frame") & !all(colnames("data.frame") %in% c("treated", "control", "distance"))) {
    stop("Distance argument is not a canonical matching problem (an adjacency list of the graph): A data.frame with columns `treated`, `control`, `distance`.")
  }
  
  if (is.null(distance$group)) {
    stop("group not assigned before gmatch")
  }
  if (any(is.na(distance$group))) {
    stop("Missing at least one group name")
  }
  
  # NB: ORDER OF ARGUMENTS SWITCHED FROM PREV VERSION!
  mxc <- round(max.col.units) #  (formerly kt)
  mnc <- round(min.col.units) #  (formerly ktl)
  mxr <- round(max.row.units)
  
  if (mnc > 1) {
    mxr <- 1
  }
  
  # Check that matching problem is well-specified
  if (mxc < mnc) {
    stop("min.col.units may not exceed max.col.units")
  }
  
  if (any(c(mxc, mnc, mxr) < 1)) {
    stop("max and min constraints must be 1 or greater")
  }
  
  if (!is.numeric(f) | f > 1 | f < 0) {
    stop("f must be a number in [0,1]")
  }
  
  # Ok, let's see how a big a problem we are working on
  # this returns a "canonical" matching: a data.frame with
  # three columns: treated, control, distance. Where the first two
  # are factors and the last is numeric.
  treated.units <- unique(distance$treated)
  control.units <- unique(distance$control)
  group.units   <- unique(distance$group)
  # group.units   <- factor(group.units, levels = group.units)
  nt <- length(treated.units)
  nc <- length(control.units)
  ncg <- length(group.units)
  narcs <- dim(distance)[1]
  problem.size <- narcs + nt + nc + ncg
  
  if (problem.size > getMaxProblemSize()) {
    stop(paste('Maximum matching problem may have only',
               signif(getMaxProblemSize(), 2), "- (nrows + ncols + 2) finite entries;",
               problem.size - getMaxProblemSize(), 'too many.',
               "Set 'options(\"optmatch_max_problem_size\" = Inf)' to disable this check."),
         call. = FALSE)
  }
  
  
  if (any(as.integer(distance$distance) != distance$distance |
          distance$distance < 0)) {
    stop("distance should be nonnegative integer")
  }
  
  # SHOULD PROBABLY DISABLE NEXT TWO WARNINGS
  if (mxc != max.col.units | mnc!=min.col.units |
      (mxr == round(max.row.units) & mxr != max.row.units)) {
    warning("fmatch coerced one or more constraints to integer")
  }
  
  if (mnc > 1 & round(max.row.units) > 1) {
    warning("since min.col.units > 1, fmatch coerced max.row.units to 1")
  }
  
  
  #################################################
  # fmatch code, leaving it in here for reference #
  #################################################
  # # set up the problem for the Fortran algorithm
  # # each node has a integer ID number
  # # startn indicates where each arc starts (using ID num)
  # # endn indicates where each arc ends (using ID num)
  # # nodes 1:nt are the treated units
  # # nodes (nt + 1):nc are the control units
  # # we use the levels of the treated and control factors to generate the ID numbers
  # # the capacity of these arcs is 1
  #
  # dists <- as.vector(distance$distance) + stability.increment
  # startn <- as.numeric(distance$treated)
  # endn <- nt + as.numeric(distance$control)
  # ucap <- rep(1, narcs)
  #
  # # Add entries for "end" and "sink" nodes
  # # "end" is node nt+nc+1; "sink" is node nt+nc+2
  # dists <- c(dists, rep(0, nc + nt), rep(0, nc))
  # startn <- c(startn, 1:(nt + nc), nt + 1:nc)
  # endn <- c(endn, rep(nt + nc + 1, nc + nt), rep(nt + nc + 2, nc))
  # ucap <- c(ucap, rep(mxc - mnc, nt), rep(mxr - 1, nc), rep(1, nc))
  #
  # # supply
  # b <- c(rep(mxc, nt), rep(0, nc), -(mxc * nt - round(f * nc)), -round(f * nc))
  #
  
  
  
  ######################################
  # Set up the problem for callrelax() #
  ######################################
  
  if (!allow_duplicates) {
    # Scenario 1
    # f is defined relative to ncg
    target <- round(f * ncg)
    
    nodes <-
      do.call(what = "rbind",
              list(data.frame(type = "T",
                              name = treated.units,
                              b    = mxc,
                              id   = as.numeric(treated.units)),
                   data.frame(type = "C",
                              name = control.units,
                              b    = 0,
                              id   = as.numeric(control.units) + nt),
                   data.frame(type = "G",
                              name = group.units,
                              b    = 0,
                              id   = as.numeric(group.units) + nt + nc),
                   data.frame(type = c("S", "O"),
                              name = c("Sink", "Overflow"),
                              b    = c(-1*target, target - nt*mxc),
                              id   = 1:2 + nt + nc + ncg)))
    udistance <- distance[!duplicated(distance$control), ]
    
    edges <-
      do.call(what = "rbind",
              list(data.frame(type   = "TC",
                              startn = as.numeric(distance$treated),
                              endn   = as.numeric(distance$control) + nt,
                              ucap   = 1,
                              dist   = distance$distance),
                   data.frame(type   = "CG",
                              startn = as.numeric(udistance$control) + nt,
                              endn   = as.numeric(udistance$group) + nt + nc,
                              ucap   = 1,
                              dist   = 0),
                   data.frame(type   = "TO",
                              startn = nodes$id[nodes$type == "T"],
                              endn   = nodes$id[nodes$type == "O"],
                              ucap   = mxc - mnc,
                              dist   = 0),
                   data.frame(type   = "GS",
                              startn = nodes$id[nodes$type == "G"],
                              endn   = nodes$id[nodes$type == "S"],
                              ucap   = 1,
                              dist   = 0)))
  } else {
    # Scenario 2
    # f defined relative to nc
    target <- round(f * nc)
    
    distance$T.id <- as.numeric(distance$treated)
    distance$G.id <- as.numeric(distance$group)
    distance$TG   <- factor(paste(distance$T.id, distance$G.id, sep = "."))
    TG.units <- unique(distance$TG)
    ntcg <- length(TG.units)
    
    nodes <-
      do.call(what = "rbind",
              list(data.frame(type = "T",
                              name = treated.units,
                              b    = mxc,
                              id   = as.numeric(treated.units)),
                   data.frame(type = "G",
                              name = TG.units,
                              b    = 0,
                              id   = as.numeric(TG.units) + nt),
                   data.frame(type = "C",
                              name = control.units,
                              b    = 0,
                              id   = as.numeric(control.units) + nt + ntcg),
                   data.frame(type = c("S", "O"),
                              name = c("Sink", "Overflow"),
                              b    = c(-1*target, target - nt*mxc),
                              id   = 1:2 + nt + ntcg + nc)))
    
    udistance <- distance[!duplicated(distance$TG), ]
    if(replace_value == FALSE){
      ucap_value <- 1
    } else {
      ucap_value <- 100
    }
    
    edges <-
      do.call(what = "rbind",
              list(data.frame(type   = "GC",  # One for each row in distance
                              startn = as.numeric(distance$TG) + nt,
                              endn   = as.numeric(distance$control) + nt + ntcg,
                              ucap   = 1,
                              dist   = distance$distance),
                   
                   data.frame(type   = "TG",  # One for each TG
                              startn = as.numeric(udistance$treated),
                              endn   = as.numeric(udistance$TG) + nt,
                              ucap   = 1,
                              dist   = 0),
                   data.frame(type   = "TO",
                              startn = nodes$id[nodes$type == "T"],
                              endn   = nodes$id[nodes$type == "O"],
                              ucap   = mxc - mnc,
                              dist   = 0),
                   data.frame(type   = "CS",
                              startn = nodes$id[nodes$type == "C"],
                              endn   = nodes$id[nodes$type == "S"],
                              ucap   = ucap_value, # 1 is without replacement, increase for with replacement
                              dist   = 0)))
  }
  
  net <- list(startn = as.integer(edges$startn),
              endn   = as.integer(edges$endn),
              ucap   = as.integer(edges$ucap),
              cost   = as.integer(edges$dist),
              b      = as.integer(nodes$b))
  fnet <- rcbalance::callrelax(net)
  
  # # If the user specifies using the old version of the relax algorithm. The `if` will be
  # # FALSE if use_fallback_optmatch_solver is anything but TRUE, including NULL.
  # # We have to duplicate the .Fortran code to make R CMD Check not complain about "registration" problems
  # if (identical(options()$use_fallback_optmatch_solver, TRUE)) {
  #   fop <- .Fortran("relaxalgold",
  #                   as.integer(nc + nt + 2),
  #                   as.integer(length(startn)),
  #                   as.integer(startn),
  #                   as.integer(endn),
  #                   as.integer(dists),
  #                   as.integer(ucap),
  #                   as.integer(b),
  #                   x1=integer(length(startn)),
  #                   crash1=as.integer(0),
  #                   large1=as.integer(.Machine$integer.max/4),
  #                   feasible1=integer(1),
  #                   NAOK = FALSE,
  #                   DUP = TRUE,
  #                   PACKAGE = "optmatch")
  # } else {
  #   fop <- .Fortran("relaxalg",
  #                   as.integer(nc + nt + 2),
  #                   as.integer(length(startn)),
  #                   as.integer(startn),
  #                   as.integer(endn),
  #                   as.integer(dists),
  #                   as.integer(ucap),
  #                   as.integer(b),
  #                   x1=integer(length(startn)),
  #                   crash1=as.integer(0),
  #                   large1=as.integer(.Machine$integer.max/4),
  #                   feasible1=integer(1),
  #                   NAOK = FALSE,
  #                   DUP = TRUE,
  #                   PACKAGE = "optmatch")
  # }
  #
  
  if (fnet$crash==0 & fnet$feasible==1) {
    edges$flow <- fnet$x
  }
  
  feas <- fnet$feasible & ((mnc*nt <= target & mxc*nt >= target) |
                             (target <= nt & target*mxr >= nt))
  x <- feas * fnet$x - (1 - feas)
  ans <- numeric(narcs)
  ans <- x[1:narcs]
  return(cbind(distance, solution = ans))
}
