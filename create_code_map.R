# Get a list of all vegan functions and their dependencies
source("https://raw.githubusercontent.com/MangoTheCat/remotes/master/install-github.R")$value("mangothecat/functionMap")
library("functionMap")

fun_map <- map_r_package(path = "/Users/lschreib/Downloads/vegan")

# Get global overview 
print(fun_map)

# This shows the indivvidual functions and their associated calls
fun_map$edge_df

# Export as tsv file so that we can work through these dependencies for our functions of interest
install.packages("readr")
library("readr")

write_tsv(as.data.frame(fun_map$edge_df), file = "vegan_2.6-10_funs-deps.tsv")