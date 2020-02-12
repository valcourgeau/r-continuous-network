# loading the packages
source("package_to_load.R")

# We do not use the load data but solar since load is too correlated and not very random

# Functions and procedures to clean the data
data_path <- "~/GitHub/r-continuous-network/data/re-europe/"

n_df_load <- 1000
n_nodes <- 70
df_load <- read.csv(paste(data_path, "Nodal_TS/wind_signal_COSMO.csv", sep=""), nrows = n_df_load+10)[,2:(n_nodes+1)]
df_load <- df_load[-c(1:10),]

clean_wind_data <- CleanData(df_load)
core_wind <- clean_wind_data$remainders

# EXPLORATION PLOT
par(mfrow=c(1,1))
clrs <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(n_nodes))
plot(core_wind[,1], type="l", ylim=c(-0.4,0.53), ylab="Hourly load in MWh",
     main="Load across nodes", col = clrs$color[1])
for(i in 2:n_nodes){
  lines(core_wind[,i], type="l", col=clrs$color[i])
}


