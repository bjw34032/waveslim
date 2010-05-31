exchange <-
  ts(as.matrix(read.table("exchange.txt", col.names=c("DEM.USD", "JPY.USD"))),
               start=1970, freq=12)
