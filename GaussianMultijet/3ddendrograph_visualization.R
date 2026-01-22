library(rgl)
library(dendextend)

df <- list()

i = 0
for(i in 0:3){
  file_prefix <- "/home/smagedov/local/src/sjet-master/GaussianMultijet/multijetHistFile/multijetHistFile_0_"
  file_prefix <- paste0(file_prefix, as.character(i))
  file_prefix <- paste0(file_prefix, ".txt")
  maxdist <- as.double(readLines(file_prefix, n = 1))
  temp <- read.table(file_prefix, header=TRUE, sep=" ", skip=1)
  df <- append(df, list(temp))
}

lines <- readLines("/home/smagedov/local/src/sjet-master/GaussianMultijet/multijetHistFile/dijetStudyExampleOutput.txt")

# Find the line numbers of the event markers
start <- grep("^# Event 0", lines)
end   <- grep("^# Event 1", lines)

# Extract the lines between them (exclusive)
data_block <- lines[(start + 1):(end - 1)]

# Convert to a data frame (assuming space-separated numeric data)
clustering_data <- read.table(text = data_block)

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

ref_jet <- data.frame()
probe_one <- data.frame()
probe_two <- data.frame()
probe_three <- data.frame()

i = 0
j = 0
for(i in 1:length(particle_list)){
  for(j in 1:nrow(particle_list[[i]])){
    node = particle_list[[i]][j, 1]+1
    if (clustering_data[node,4]==0){
      ref_jet <- rbind(ref_jet, particle_list[[i]][j,])
    } else if (clustering_data[node,4]==1) {
      probe_one <- rbind(probe_one, particle_list[[i]][j,])
    } else if (clustering_data[node,4]==2) {
      probe_two <- rbind(probe_two, particle_list[[i]][j,])
    } else {
      probe_three <- rbind(probe_three, particle_list[[i]][j,])
    }
  }
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

tot_pt <- list()
for(i in 1:length(cluster_list)){
  tot_pt <<- append(tot_pt, cluster_list[[i]][1,5])
}


ref_jet_pt <- ref_jet$pt/tot_pt[[1]]
probe_one_pt <- probe_one$pt/tot_pt[[2]]
probe_two_pt <- probe_two$pt/tot_pt[[3]]
probe_three_pt <- probe_three$pt/tot_pt[[4]]
clust1_pt <- cluster_list[[1]]$pt/tot_pt[[1]]
clust2_pt <- cluster_list[[2]]$pt/tot_pt[[2]]
clust3_pt <- cluster_list[[3]]$pt/tot_pt[[3]]
clust4_pt <- cluster_list[[4]]$pt/tot_pt[[4]]

open3d()
bg3d("black")

# Plot all segments (edges)
for (seg in full_seg_list[[1]]) {
  lines3d(seg, col = "blue", lwd = 2)
}
for (seg in full_seg_list[[2]]) {
   lines3d(seg, col = "red", lwd = 2)
 }
for (seg in full_seg_list[[3]]) {
  lines3d(seg, col = "yellow", lwd = 2)
}
for (seg in full_seg_list[[4]]) {
  lines3d(seg, col = "magenta", lwd = 2)
}
# Plot leaf points
size_scale = 0.15

spheres3d(ref_jet$eta, ref_jet$phi, ref_jet$dist, col = "red", radius = (ref_jet_pt)^(1/3)*size_scale)
spheres3d(probe_one$eta, probe_one$phi, probe_one$dist, col = "blue", radius = (probe_one_pt)^(1/3)*size_scale)
spheres3d(probe_two$eta, probe_two$phi, probe_two$dist, col = "yellow", radius = (probe_two_pt)^(1/3)*size_scale)
spheres3d(probe_three$eta, probe_three$phi, probe_three$dist, col = "magenta", radius = (probe_three_pt)^(1/3)*size_scale)
spheres3d(cluster_list[[1]]$eta, cluster_list[[1]]$phi, cluster_list[[1]]$dist, col = "blue", radius = (clust1_pt)^(1/3)*size_scale)
spheres3d(cluster_list[[2]]$eta, cluster_list[[2]]$phi, cluster_list[[2]]$dist, col = "red", radius = (clust2_pt)^(1/3)*size_scale)
spheres3d(cluster_list[[3]]$eta, cluster_list[[3]]$phi, cluster_list[[3]]$dist, col = "yellow", radius = (clust3_pt)^(1/3)*size_scale)
spheres3d(cluster_list[[4]]$eta, cluster_list[[4]]$phi, cluster_list[[4]]$dist, col = "magenta", radius = (clust4_pt)^(1/3)*size_scale)

labels = c("Reference Jet", "Probe One", "Probe Two", "Probe Three")
colors = c("red", "blue", "yellow", "magenta")

legend_x <- 3
legend_y <- 0
legend_z <- seq(0.3, 0.4, length.out = length(labels))

points3d(rep(legend_x, 2), rep(legend_y, 2), legend_z, col = colors, size = 6)
text3d(rep(legend_x + 0.15, 2), rep(legend_y, 2), legend_z, texts = labels, adj = 0)

aspect3d(1, 1, 1.1)

axes3d(edges = c("x--", "y--", "z--"), col = "white", lwd = 2)
title3d(main = "3D Dendrogram of Jet Clustering", xlab = "Eta", ylab = "Phi", zlab = "Distance")

