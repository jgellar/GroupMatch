### Tests for the code ###
# Load in packages
library(tidyverse)

# set working directory
#setwd("~/Desktop/GroupMatch/R")

# source code
# source('zzzDistanceSpecification.R')
# source('zzz.R')
# source('utilities.R')
# source('rollingMatch.R')
# source('Optmatch.R')
# source('mdist.R')
# source('gSubDivStrat.R')
# source('groupmatch.R')
# source('gmatch.R')
# source('feasible.R')

# Units 1 & 7 and 2 & 6 should be matched, and all have wt 1
dat <- data.frame(treat = c(1, 1, 0, 0, 0, 0, 0),
                  personID = c(1, 2, 3, 3, 4, 4, 4),
                  X1 = c(1, 2, 0, 0, 0, 2, 1),
                  X2 = c(1, 2, 0, 0, 0, 2, 1))
psDist <- match_on(treat ~ X1 + X2, data = dat)
match <- groupmatch(psDist, min.controls = 1, max.controls = 2,
                    allow_duplicates = as.logical(1), group = ~personID, data = dat, 
                    replace_value = TRUE)
match$cells
matrix_to_weights(match$cells, dat)
matrix_to_set(match)

# Change the ordering of the treated and control units to make sure ordering is preserved
dat <- data.frame(treat = c(1, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 1),
                  personID = c(1, 4,
                               3, 3, 3, 3, 3,
                               4, 4, 4, 4,4,2),
                  X1 = c(2, 2.5,
                         -1, 2, 1, 1.5, 1,
                         0, 1, 1.6, -1, 2, 2.5))
psDist <- match_on(treat ~ X1, data = dat)
match <- groupmatch(psDist, min.controls = 2, max.controls = 2,
                    allow_duplicates = as.logical(1), group = ~personID, 
                    data = dat, replace_value = FALSE)
match$cells
matrix_to_weights(match$cells, dat)
matrix_to_set(match)

# Match each treated unit to 2 controls
dat <- data.frame(treat = c(1, 1,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0),
                  personID = c(1, 2,
                               3, 3, 3, 3, 3,
                               4, 4, 4, 4, 4, 4),
                  X1 = c(1.6748082, 1.0067069,
                         -0.9902554, 0.8125199, -0.8709979, -1.4816648, 0.7157801,
                         0.4030891, 0.7158147, 0.7277648, -0.7860486, 3.4868057, 0.1183520),
                  X2 = c(0.9452371, 1.8184389,
                         2.0007598, -1.7643497, -0.8343905, -0.1206871, 0.4981403,
                         -0.5475328, -1.4133694, 0.3800757, 0.5239266, -0.5836660, -0.5212903))
psDist <- match_on(treat ~ X1 + X2, data = dat)
match <- groupmatch(psDist, min.controls = 2, max.controls = 2,
                    allow_duplicates = as.logical(1), group = ~personID, data = dat,
                    replace_value = TRUE)
match$cells
matrix_to_weights(match$cells, dat)

# try with variable number of controls
matchOut <- groupmatch(psDist, min.controls = 1, max.controls = 3,
                    allow_duplicates = as.logical(1), group = ~personID, data = dat,
                    replace_value = TRUE)
matchOut$cells

# change allow duplicates to false
match <- groupmatch(psDist, min.controls = 1, max.controls = 1,
                    allow_duplicates = as.logical(0), group = ~personID, data = dat,
                    replace_value = TRUE)
match$cells
matrix_to_weights(match$cells, dat)
