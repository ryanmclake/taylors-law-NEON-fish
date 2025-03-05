# clear your environment and memory
rm(list=ls())
gc()


# download and load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(neonUtilities,tidyverse,nlme,MuMIn,FSA,reshape2,trend, minpack.lm, lme4)

# download and unpack neon fish data
fish <- try(neonUtilities::loadByProduct(
  dpID='DP1.20107.001',
  check.size= F, 
  startdate = "2014-01",
  enddate = "2024-10",
  release = "LATEST",
  package='expanded',
  include.provisional = T,
  token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJtY2NsdXJlckBiYXR0ZWxsZWVjb2xvZ3kub3JnIiwic2NvcGUiOiJyYXRlOnVubGltaXRlZCByZWFkOnJlbGVhc2VzIHJlYWQ6cmVsZWFzZXMtbGF0ZXN0IiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4OTIyMzc5NzIsImlhdCI6MTczNDU1Nzk3MiwiZW1haWwiOiJtY2NsdXJlckBiYXR0ZWxsZWVjb2xvZ3kub3JnIn0.a2CvXWQ60JBuCooRUn3_uBbtouzJdGk5LIsCBQEqZJ65wpiGueO_mpG9_AXxkw937W69MWCnQBcNL7UCNHpk6w"))

if(is.null(attr(fish, "class"))){
  invisible(list2env(fish, envir=.GlobalEnv))
}else{
  cat("There are no data in this date range to check, no script outputs generated.")
  knitr::knit_exit()
}

