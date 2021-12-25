# rob_lemons_run_script
# 2021-12-22 Simon Dedman & Mo VZB

library(dBBMM_HomeRange) # change as required
DET <- read.csv("../data/TRACKS1.csv")

DET %<>%
  mutate(Datetime = as.POSIXct(Datetime,
                               format = "%m/%d/%y %H:%M",
                               tz = dat.TZ
  ),
  Shark = make.names(Shark) # prefixes "X" to numerical-named sharks to avoid issues later
  ) %>%
  rename(
    Lat = N,
    Lon = W,
    T.Ph = Tidal.Phase
  ) %>%
  select(Datetime, Shark, T.Ph, Lat, Lon) %>% # dropped 2 variables (Date, Time)
  arrange(Shark, Datetime) # df size: 1308 x 5

write.csv(x = DET, file = "TracksCleaned.csv", row.names = FALSE) # Could load this directly here

dBBMM_HomeRange(
  data = DET,
)

scaleraster(
  
)

dBBMM_plot(
  
)