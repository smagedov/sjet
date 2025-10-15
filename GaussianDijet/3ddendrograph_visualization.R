library(rgl)
library(dendextend)

df <- list()

i = 0
for(i in 0:1){
  file_prefix <- "/home/smagedov/local/src/sjet-master/GaussianDijet/dijetHistFile/dijetHistFile_213_"
  file_prefix <- paste0(file_prefix, as.character(i))
  file_prefix <- paste0(file_prefix, ".txt")
  maxdist <- as.double(readLines(file_prefix, n = 1))
  temp <- read.table(file_prefix, header=TRUE, sep=" ", skip=1)
  df <- append(df, list(temp))
}

cluster_list <- list()
particle_list <- list()

i = 0
for(i in 1:length(df)){
  nclusts <- nrow(df[[i]])
  particle_df <- df[[1]][FALSE,]
  cluster_df <- df[[1]][FALSE,]
  j = 0
  for(j in 1:nclusts) {
    if(df[[i]][j,2]==-1 && df[[i]][j,3]==-1) {
      particle_df <- rbind(particle_df, df[[i]][j,])
      } else {
        cluster_df <- rbind(cluster_df, df[[i]][j,])
      }
  }
  particle_list <- append(particle_list, list(particle_df))
  cluster_list <- append(cluster_list, list(cluster_df))
}

# Recursive function to get 3D segments

get_segments_3d <- function(l, m, ptot, x_offset = 0) {
  px = cluster_list[[l]][m, 6]
  py = cluster_list[[l]][m, 7]
  pz = cluster_list[[l]][m, 4]
  
  parent1 <- as.integer(cluster_list[[l]][m, 2])
  parent2 <- as.integer(cluster_list[[l]][m, 3])
  node = 0
  if(parent1 <= (ptot-1)){
    for(n in 1:nrow(particle_list[[l]])){
      if(particle_list[[l]][n,1] == parent1){
        node = n
      }
    }
    x1 = particle_list[[l]][node, 6]
    y1 = particle_list[[l]][node, 7]
    z1 = particle_list[[l]][node, 4]
  } else {    
    for(n in 1:nrow(cluster_list[[l]])){
      if(cluster_list[[l]][n,1] == parent1){
        node = n
      }
    }
    x1 = cluster_list[[l]][node, 6]
    y1 = cluster_list[[l]][node, 7]
    z1 = cluster_list[[l]][node, 4]
  }
  
  node = 0
  if(parent2 <= (ptot-1)){
    for(n in 1:nrow(particle_list[[l]])){
      if(particle_list[[l]][n,1] == parent2){
        node = n
      }
    }
    x2 = particle_list[[l]][node, 6]
    y2 = particle_list[[l]][node, 7]
    z2 = particle_list[[l]][node, 4]
  } else {
    for(n in 1:nrow(cluster_list[[l]])){
      if(cluster_list[[l]][n,1] == parent2){
        node = n
      }
    }
    x2 = cluster_list[[l]][node, 6]
    y2 = cluster_list[[l]][node, 7]
    z2 = cluster_list[[l]][node, 4]
  }

  left_step1 <- list(rbind(c(x1, y1, z1), c(x1, y1, pz)))
  left_step2 <- list(rbind(c(x1, y1, pz), c(px, py, pz)))
    
  # Right connector similarly
  right_step1 <- list(rbind(c(x2, y2, z2), c(x2, y2, pz)))
  right_step2 <- list(rbind(c(x2, y2, pz), c(px, py, pz)))
    
  # Add segments to list
  # Save segments: child -> parent
  segments_list <<- append(segments_list, left_step1)
  segments_list <<- append(segments_list, left_step2)
  segments_list <<- append(segments_list, right_step1)
  segments_list <<- append(segments_list, right_step2)
}

full_seg_list <- list()

# Build the dendrogram and gather segments
particle_tot = 0
for(i in 1:length(particle_list)){
  particle_tot = particle_tot + nrow(particle_list[[i]])
}

for(i in 1:length(cluster_list)){
  segments_list <- list()
  for(j in 1:nrow(cluster_list[[i]])){
    get_segments_3d(i, j, particle_tot)
  }
  final_step <- rbind(c(cluster_list[[i]][1,6],cluster_list[[i]][1,7],cluster_list[[i]][1,4]),c(cluster_list[[i]][1,6],cluster_list[[i]][1,7],maxdist))
  segments_list <- append(segments_list, list(final_step))
  full_seg_list <- append(full_seg_list, list(segments_list))
}

open3d()
bg3d("white")

# Plot all segments (edges)
for (seg in full_seg_list[[1]]) {
  lines3d(seg, col = "blue", lwd = 2)
}
for (seg in full_seg_list[[2]]) {
   lines3d(seg, col = "magenta", lwd = 2)
 }

# Plot leaf points
points3d(particle_list[[1]]$eta, particle_list[[1]]$phi, particle_list[[1]]$dist, col = "red", size = 8)
points3d(cluster_list[[1]]$eta, cluster_list[[1]]$phi, cluster_list[[1]]$dist, col = "red", size = 8)

points3d(particle_list[[2]]$eta, particle_list[[2]]$phi, particle_list[[2]]$dist, col = "pink", size = 8)
points3d(cluster_list[[2]]$eta, cluster_list[[2]]$phi, cluster_list[[2]]$dist, col = "pink", size = 8)

axes3d(edges = c("x--", "y--", "z--"), col = "black", lwd = 2)
title3d(main = "3D Dendrogram of Jet Clustering", xlab = "Eta", ylab = "Phi", zlab = "Distance")

