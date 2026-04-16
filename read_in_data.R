############################# read in data #####################################

### load packages --------------------------------------------------------------

library(showtext)
library(data.table)
library(readxl)
library(lubridate)
library(OptimalDesign)
library(fitdistrplus)
library(matrixcalc)
library(reshape2)
library(viridis)
library(writexl)
#library(asreml) # commented out because asreml is not available on CRAN and requires a license: More info in the readme file of this repository
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rstan)
library(MCMCpack)
library(MASS)
library(coda)
library(tidyr)
library(ggmcmc)
library(GGally)
library(SpATS)
library(fitdistrplus)
library(patchwork)
library(purrr)
library(stringr)
library(htmltools)
library(knitr)
library(grid)
library(base64enc)


find_project_root <- function() {
  candidates <- character()

  add_candidate <- function(path) {
    if (!is.character(path) || length(path) == 0 || is.na(path) || !nzchar(path)) {
      return(invisible(NULL))
    }

    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    candidates <<- c(candidates, normalized, file.path(normalized, "project"))
  }

  add_candidate(getwd())
  add_candidate("/home/rstudio/project")

  file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_args) > 0) {
    add_candidate(dirname(sub("^--file=", "", file_args[[1]])))
  }

  frame_files <- Filter(
    Negate(is.null),
    lapply(sys.frames(), function(frame) {
      if (exists("ofile", envir = frame, inherits = FALSE)) {
        get("ofile", envir = frame, inherits = FALSE)
      } else {
        NULL
      }
    })
  )
  for (file_path in frame_files) {
    add_candidate(dirname(file_path))
  }

  expanded_candidates <- unique(unlist(lapply(candidates, function(path) {
    if (!nzchar(path)) {
      return(character())
    }

    parent_1 <- dirname(path)
    parent_2 <- dirname(parent_1)
    c(path, parent_1, parent_2)
  }), use.names = FALSE))

  for (path in expanded_candidates) {
    has_data <- file.exists(file.path(path, "0 Data", "yield_winter_medium.csv"))
    has_functions <- file.exists(file.path(path, "1 Functions", "My Residuals Function.R"))
    has_project <- file.exists(file.path(path, "read_in_data.R")) ||
      length(list.files(path, pattern = "\\.Rproj$", full.names = TRUE)) > 0

    if (has_data && has_functions && has_project) {
      return(path)
    }
  }

  stop(
    paste(
      "Could not locate the project root containing '0 Data' and '1 Functions'.",
      "Checked candidates:",
      paste(expanded_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

project_root <- find_project_root()
setwd(project_root)

source_project_file <- function(...) {
  source(file.path(project_root, ...))
}


### load functions -------------------------------------------------------------

## Frequentist asreml-package pipeline

source_project_file("1 Functions", "My Residuals Function.R")
source_project_file("1 Functions", "Grid of Optimal Design Function.R")

## Bayesian rstan-package pipeline

source_project_file("1 Functions", "Fit general multiple Cycles Function.R")
source_project_file("1 Functions", "Grid of Optimal Bayes Design Function.R")
source_project_file("1 Functions", "Render Bayesian Cycle Report Function.R")


### read in data ---------------------------------------------------------------

yield_winter_medium <- fread(file.path(project_root, "0 Data", "yield_winter_medium.csv"))

## set variable classes
yield_winter_medium[ , ID := factor(ID)]
yield_winter_medium[ , Release_year := factor(Release_year)]
yield_winter_medium[ , year := factor(year)]
yield_winter_medium[ , Zone := factor(Zone)]
yield_winter_medium[ , Location := factor(Location)]
yield_winter_medium[ , Genotype := factor(Genotype)]
yield_winter_medium[ , Group := factor(Group)]
yield_winter_medium[ , Rep := factor(Rep)]
yield_winter_medium[ , environment := factor(environment)]
yield_winter_medium[ , Loc_Count := factor(Loc_Count)]
yield_winter_medium[ , yield := as.numeric(yield)]

## As dataframe
yield_winter_medium <- na.omit(yield_winter_medium)
yield_winter_medium <- as.data.frame(yield_winter_medium)
