library(rgl)
library(dendextend)

df <- list()

i = 0
for(i in 0:1){
  file_prefix <- "/home/smagedov/local/src/sjet-master/GaussianDijet/dijetHistFile/dijetHistFile_0_"
  file_prefix <- paste0(file_prefix, as.character(i))
  file_prefix <- paste0(file_prefix, ".txt")
  temp <- read.table(file_prefix, header=TRUE, sep=" ")
  df <- append(df, list(temp))
}

cluster_list <- list()
particle_list <- list()
particle_tot_list <- list()

i = 0
for(i in 1:length(df)){
  nclusts <- nrow(df[[i]])
  particle_tot <- 0
  j = 0
  for(j in 1:nclusts) {
    if(df[[i]][j,1]!=-1 && df[[i]][j,2]!=-1) {
      break
      } else {
        particle_tot = j
      }
  }
  particle_tot_list <- append(particle_tot_list, particle_tot)
  
  particles <- df[[i]][1:particle_tot, ]
  clusters <- df[[i]][(particle_tot + 1):nclusts, ]
  
  particle_list <- append(particle_list, list(particles))
  cluster_list <- append(cluster_list, list(clusters))
}

# Recursive function to get 3D segments

get_segments_3d <- function(l, m, x_offset = 0) {
  ptot <- as.integer(particle_tot_list[l])
  px = cluster_list[[l]][m, 5]
  py = cluster_list[[l]][m, 6]
  pz = cluster_list[[l]][m, 3]
  
  parent1 <- as.integer(cluster_list[[l]][m, 1]) + 1
  parent2 <- as.integer(cluster_list[[l]][m, 2]) + 1
  if(parent1 <= (ptot+1)){
    x1 = particle_list[[l]][parent1, 5]
    y1 = particle_list[[l]][parent1, 6]
    z1 = particle_list[[l]][parent1, 3]
  } else {
    x1 = cluster_list[[l]][parent1-ptot, 5]
    y1 = cluster_list[[l]][parent1-ptot, 6]
    z1 = cluster_list[[l]][parent1-ptot, 3]
  }
  
  if(parent2 <= (ptot+1)){
    x2 = particle_list[[l]][parent2, 5]
    y2 = particle_list[[l]][parent2, 6]
    z2 = particle_list[[l]][parent2, 3]
  } else {
    x2 = cluster_list[[l]][parent2-ptot, 5]
    y2 = cluster_list[[l]][parent2-ptot, 6]
    z2 = cluster_list[[l]][parent2-ptot, 3]
  }
  
  left_step1 <- list(rbind(c(x1, y1, z1), c(x1, y1, pz/2)))
  left_step2 <- list(rbind(c(x1, y1, pz/2), c(px, py, pz/2)))
    
  # Right connector similarly
  right_step1 <- list(rbind(c(x2, y2, z2), c(x2, y2, pz/2)))
  right_step2 <- list(rbind(c(x2, y2, pz/2), c(px, py, pz/2)))
  
  final_step <- list(rbind(c(px, py, pz/2), c(px, py, pz)))
    
  # Add segments to list
  # Save segments: child -> parent
  segments_list <<- append(segments_list, left_step1)
  segments_list <<- append(segments_list, left_step2)
  segments_list <<- append(segments_list, right_step1)
  segments_list <<- append(segments_list, right_step2)
  segments_list <<- append(segments_list, final_step)
}

full_seg_list <- list()

# Build the dendrogram and gather segments
for(i in 1:length(cluster_list)){
  segments_list <- list()
  for(j in 1:(nrow(cluster_list[[i]])-1)){
    get_segments_3d(i,j)
  }
  full_seg_list <- append(full_seg_list, list(segments_list))
}

open3d()
bg3d("white")

# Plot all segments (edges)
for (seg in full_seg_list[[1]]) {
  lines3d(seg, col = "blue", lwd = 2)
}
for (seg in full_seg_list[[2]]) {
  lines3d(seg, col = "yellow", lwd = 2)
}

# Plot leaf points
points3d(particle_list[[1]]$eta, particle_list[[1]]$phi, particle_list[[1]]$dist, col = "red", size = 8)
points3d(cluster_list[[1]]$eta, cluster_list[[1]]$phi, cluster_list[[1]]$dist, col = "red", size = 8)

# Plot leaf points
points3d(particle_list[[2]]$eta, particle_list[[2]]$phi, particle_list[[2]]$dist, col = "green", size = 8)
points3d(cluster_list[[2]]$eta, cluster_list[[2]]$phi, cluster_list[[2]]$dist, col = "green", size = 8)

axes3d(edges = c("x--", "y--", "z--"), col = "black", lwd = 2)
title3d(xlab = "Eta", ylab = "Phi", zlab = "Distance")

