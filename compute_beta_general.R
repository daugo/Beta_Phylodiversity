rm(list=ls())
#####Packages######
library(raster)
library(picante)
library(multicore)

#Set Working directory
setwd("~/Dropbox/IAvH/betaphylodiversity/")

#Load community file and clean NA values
community<-read.table("communities.txt",h=T)
community<-na.omit(community)

#Load phylogeny
phylogeny<-read.nexus("consenso_1k_sumtrees.nex")

#Create Raster containing pixels numbers 
#For the Raster to have an appropiate extent and resolution, raster is created from an available species distribution map in .asc
r<-raster("Aburria_aburri.asc")
pix_ras<-r
names(pix_ras)<-"Grid"
values(pix_ras)<-NA
#2- Asignar al raster los valores de PD que corresponden a cada pixel
pix_ras[as.numeric(rownames(community))]<-as.numeric(rownames(community))

######BetaPD################
### funcion con indice 1
beta1<-function(x){
  pixel<-x[length(x)/2+0.5]
  if(is.na(pixel)==F && sum(!is.na(x)) > 1) {
    beta_phylosor<-as.matrix( phylosor( subset(community, rownames(community) %in% x[!is.na(x)]), phylogeny) )
    row <- as.character(pixel)
    value<-vector()
    for (i in 1:length(x)) {
        if( i!= x[length(x)/2+0.5]) {
          column<-as.character(x[i])
          if(is.na(column)==F){
            value[i]<- beta_phylosor[row, column]
          }
        }
    }
    beta<-mean(value,na.rm=T)
    return(beta)
  }
  else return(NA)
}


######Get ext meshes#######
get_exts_meshes <- function(pix_ras=pix_ras, mesh1=mesh1,  mesh2=mesh2) {
  
  pix_col_mesh1 <- floor(seq(0,ncol(pix_ras),length.out=ceiling(mesh1[1])+1))
  pix_row_mesh1 <- floor(seq(0,nrow(pix_ras),length.out=ceiling(mesh1[2])+1))
  
  pix_col_mesh2 <- floor(seq(0,ncol(pix_ras),length.out=ceiling(mesh2[1])+1))
  pix_row_mesh2 <- floor(seq(0,nrow(pix_ras),length.out=ceiling(mesh2[2])+1))
  
  rep_col <- setdiff(intersect(pix_col_mesh1, pix_col_mesh2),c(0,ncol(pix_ras)))
  rep_row <- setdiff(intersect(pix_row_mesh1, pix_row_mesh2),c(0,nrow(pix_ras)))
  
  if (length(rep_col) != 0) {
    positions2change <- which(pix_col_mesh2 %in% rep_col)
    for(i in positions2change) {
      pix_col_mesh2[i] <- pix_col_mesh2[i] + 5
    }
  }
  
  pixels2split_l <- list(col_mesh1=pix_col_mesh1, row_mesh1=pix_row_mesh1, col_mesh2=pix_col_mesh2, row_mesh2=pix_row_mesh2)
  return(pixels2split_l)
}

####Crop Rasters###########
crop_rasters <- function(pix_ras=pix_ras, col_ext, row_ext) {
  c <- 1
  extents_jobs <- vector()
  c_i <- 1
  for (i in col_ext ) {
    c_j <- 1
    for ( j in row_ext ) {
      if ( !is.na(col_ext[c_i+1]) && !is.na(row_ext[c_j+1]) ) {
        extents_jobs <- rbind(extents_jobs, c(i, col_ext[c_i+1], j, row_ext[c_j+1] ) )
      }
      c_j = c_j +1
      c = c + 1
    }
    c_i = c_i + 1
  }
  pix_ras_clusters <- apply(extents_jobs,1, crop,x=pix_ras)
  return(pix_ras_clusters)
}


####Partition Jobs#########
partition_jobs <- function(pix_ras=pix_ras, mesh1=mesh1,  mesh2=mesh2) {
  
  pixels2split_list <- get_exts_meshes(pix_ras, mesh1,  mesh2)
  
  col_ext_m1 <- xmin(pix_ras) + res(pix_ras)[2]* pixels2split_list$col_mesh1
  col_ext_m2 <- xmin(pix_ras) + res(pix_ras)[2]* pixels2split_list$col_mesh2
  row_ext_m1 <- ymin(pix_ras) + res(pix_ras)[1]* pixels2split_list$row_mesh1
  row_ext_m2 <- ymin(pix_ras) + res(pix_ras)[1]* pixels2split_list$row_mesh2
  
  pix_ras_clusters_m1 <- crop_rasters(pix_ras, col_ext_m1, row_ext_m1)
  pix_ras_clusters_m2 <- crop_rasters(pix_ras, col_ext_m2, row_ext_m2)
  pix_ras_clusters <- list(mesh1_cl=pix_ras_clusters_m1, mesh2_cl=pix_ras_clusters_m2)
  return(pix_ras_clusters)
}


####Neighbours_complete######
neighbours_complete <- function(x) {
  n <- length(which( (x >= 0) == T) )
  n<- n + length( which(is.na(x) == T ) )
  if (n == length(x) ) { return(T) }
  else { return(F) }
}


###Pixel2Keep#####
pixel2keep <- function(x,y,k,m) {
  if ( k == 1 || m == 1 ) {
    if (k == 1 && m ==1 && x == y) {
      return (x)
    }
    else if (k == 1 && m != 1 ) {
      return (x)
    }
    else if (k != 1 && m == 1 ) {
      return (y)
    }
  }
  else {
    warning( "both windows overlay at a corner: try resizing one of the meshes" )
  }
}

#####It overlays? #####
it_overlays <- function( k , m) {
  suppressWarnings(try(
    if (class(intersect( k, m) ) == "RasterLayer") { 
        return(T)
      }
    , silent=T))
  return(F)
}

merge_partitions <- function (k , m) {
  
}

check_borders <- function(x,y) {
  if (x == y) {
    return(x)
    } 
  else if (is.na(x)) {
    return(y)
    } 
  else {
    return(x)
    }
}

#### La funcion focal recorre el raster y para cada pixel genera una ventana alrededor de las dimensiones especificadas en w.  En este caso una ventana de 3*3 pixeles y para la ventana seleccionada aplica la funcion beta

#copute betaphylodiversity by windows using neighbours geometry define in w
#ras_beta_phylosor<-raster::focal(pix_ras,w=matrix(1,nrow=3,ncol=3),fun=beta1)
#piplot(ras_beta_phylosor)

#get total number of cores on machine
total_cores <- multicore:::detectCores()
#establish number of cores to use, we define to use 4/5 of the total number of cores
cores2use <- ceiling(total_cores*4/5)
#mesh geometry for data partitioning
mesh1 <- c(cores2use/2, 2)
mesh2 <- c(cores2use/3, 3)

pix_ras_clusters <- partition_jobs(pix_ras, mesh1, mesh2)

ras_beta_phylosor_clusters1 <- mclapply(pix_ras_clusters$mesh1_cl, raster::focal, w=matrix(1,nrow=3,ncol=3), fun=beta1)
#ras_neigh_clusters1 <- mclapply(pix_ras_clusters$mesh1_cl, raster::focal, w=matrix(1,nrow=3,ncol=3), fun=neighbours_complete)
ras_beta_phylosor_clusters2 <- mclapply(pix_ras_clusters$mesh2_cl, raster::focal, w=matrix(1,nrow=3,ncol=3), fun=beta1)
#ras_neigh_clusters2 <- mclapply(pix_ras_clusters$mesh2_cl, raster::focal, w=matrix(1,nrow=3,ncol=3), fun=neighbours_complete)

#beta_clu1 <- list(beta=ras_beta_phylosor_clusters1, neigh=ras_neigh_clusters1)
#beta_clu2 <- list(beta=ras_beta_phylosor_clusters2, neigh=ras_neigh_clusters2)
all_rasters <- c(ras_beta_phylosor_clusters1, ras_beta_phylosor_clusters2)
rasters.mosaicargs <- all_rasters
rasters.mosaicargs$fun <- mean
mos <- do.call(mosaic, rasters.mosaicargs)

#beta_good <- list()
#c_i <- 1
#for ( i in beta_clu1$beta ) {
#  c_j <- 1
#  for ( j in beta_clu2$beta ) {
#    if ( it_overlays(i, j) ) {
#      beta_good <- c(beta_good, overlay( i, j, fun=pixel2keep, k= beta_clu1$neigh[[c_i]], m= beta_clu2$neigh[[c_j]] ) )
#    }
#    c_i = c_i + i
#  }
#  c_j = c_j + j
#}


