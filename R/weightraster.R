#' Succinct title 8 words max
#'
#' Description paragraph: Automates delta log normal boosted regression trees abundance prediction.
#' Loops through all permutations of parameters provided (learning
#' rate, tree complexity, bag fraction), chooses the best, then simplifies it.
#' Generates line, dot and bar plots, and outputs these and the predictions
#' and a report of all variables used, statistics for tests, variable
#' interactions, predictors used and dropped, etc. If selected, generates
#' predicted abundance maps, and Unrepresentativeness surfaces.
#' See www.GitHub.com/SimonDedman/gbm.auto for issues, feedback, and development
#' suggestions. See SimonDedman.com for links to walkthrough paper, and papers
#' and thesis published using this package.
#' 
#' @export
#' @param locations No terminal slash.
#' @param weightings No terminal slash.
#' @param saveloc No terminal slash.
# Mo receiver array doublescale####
# scale to each receiver array (loop of 3) then scale those together
# 1. run dBBMM.build.R (change this filename to dBBMMhomeRange.R ?) for each array
# 2. run scaleraster on each array
# new script: # 3. Pull All_Rasters_Scaled.asc from each array/scaled subfolder & weight each (of 3) scaleraster outputs by number of cells containing any receiver (a number the user calculates). Save to a single folder.
weightraster <- function(locations = c("array1", "array2"), # assumes they have subfolders called Scaled with files called All_Rasters_Scaled.asc
                         weightings = c(1, 2),
                         saveloc = "location") # respective weightings per array
  # 4. run scaleraster on the 3 outputs (All_Rasters_Scaled.asc (which you weighted) * 3)
  # 5. run plot

    tmp <- 123

# for (f in list.files("/home/simon/Dropbox/Galway/Analysis/R/dBBMMhomeRange/R", full.names = TRUE)) {
#   print(f)
#   parse(f)
# }
