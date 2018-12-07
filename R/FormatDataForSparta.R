# Function to format data for sparta. THIS NEEDS DOCUMENTING

FormatDataForSparta <- function(All.Obs) {
  require(sparta)
  
  All.Obs <- All.Obs
  # visit: name combined to represent "site_date_year"
  All.Obs$Site <- paste(All.Obs$Site, "_", rownames(All.Obs), "_", All.Obs$year, sep = "")
  colnames(All.Obs)[1] <- "visit"
  # change 1s and 0s to TRUE and FALSE
  All.Obs[, 3:ncol(All.Obs)] <- All.Obs[, 3:ncol(All.Obs)] == 1
  # extract this now before removing the year column
  occDetData <- All.Obs[, c("visit", "year")]
  # remove year column
  All.Obs <- All.Obs[, c(1, 3:ncol(All.Obs))]
  # data element 2: occDetData - dataframe giving the site, list length and TP (time period)
  occDetData <- cbind(occDetData, apply(X = All.Obs[, 3:ncol(All.Obs)], MARGIN = 1, FUN = sum))
  
  colnames(occDetData) <- c("visit", "TP", "L")
  occDetData$site <- sub("_.*", "", occDetData$visit)
  
  # currently there are visits with no observations, don't have these in real datasets
  # remove these
  occDetData <- occDetData[!occDetData$L == 0, ]
  
  # filter off these visits with 0 LL from the All.Obs
  All.Obs <- All.Obs[All.Obs$visit %in% occDetData$visit, ]
  
  # check that these fit into the sparta function
  output <- list(occDetData = occDetData, spp_vis = All.Obs)
  output
}
