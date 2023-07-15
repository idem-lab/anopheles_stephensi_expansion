to_unzip <- list.files("data/internal_migration_flow", full.names = TRUE)

lapply(as.list(to_unzip), unzip, exdir = "data/internal_migration_flow")


to_unzip <- list.files("data/internal_migration_flow", full.names = TRUE,pattern = "_AdminUnit_Centroids.zip|_AdminUnit_Centroidss.zip")
#Eritrea data spelling error 


unzip_centroid <- function(x){
  path.x <- substr(x,1,52)
  unzip(x,exdir = path.x)
}

lapply(as.list(to_unzip), unzip_centroid)


