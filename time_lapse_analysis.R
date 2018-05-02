### Code written by Kristoffer T. Baek, 2018 ###

# Load packages -----------------------------------------------------------

library(plyr)
library(ape) # Tree plotting
library(dplyr)
library(magrittr)
library(png)
library(tidyr)

# Source functions --------------------------------------------------------

source("Time_lapse_tracking_functions.R")

# Description of input files ----------------------------------------------

## The analysis uses the two output files from TrackMate as input. The spot file contains information about each 'spot', such as spatial location in the frame, and frame number. In the case of cell tracking, a spot is the instance of one cell in one frame. The link file tells us if two spots in consecutive frames are connected (i.e. if either it's the same cell in the two frames, or if a cell in one frame became two cells in the next frame). Each row in the link file represents a link defined by a source spot (the first frame) and a target spot (the next frame). 

# Import files  --------------------------------------

WT_links <- read.table("160225_WT/TrackMate_output/Links_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

WT_spots <- read.table("160225_WT/TrackMate_output/Spots_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

WT_OX_links <- read.table("160225_WT_OXA/TrackMate_output/Links_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

WT_OX_spots <- read.table("160225_WT_OXA/TrackMate_output/Spots_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

X_A_links <- read.table("160225_X/TrackMate_output/Links_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

X_A_spots <- read.table("160225_X/TrackMate_output/Spots_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

X_A_OX_links <- read.table("160225_X_OXA/TrackMate_output/Links_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

X_A_OX_spots <- read.table("160225_X_OXA/TrackMate_output/Spots_in_tracks_statistics.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Start of analysis -------------------------------------------------------

WT <- analyze_tracks(links = WT_links, spots = WT_spots) 

WT_OX <- analyze_tracks(links = WT_OX_links, spots = WT_OX_spots)

X_A <- analyze_tracks(links = X_A_links, spots = X_A_spots)

X_OX_A <- analyze_tracks(links = X_A_OX_links, spots = X_A_OX_spots)

# Create Figure 2c --------------------------------------------------------

png("growth_curve_all.png", width = 7, height = 20, units = "cm", res = 600)
par(mar = c(2,2,2,2), mfcol = c(4,1), oma = c(4,4,2,2), cex = 0.9)

plot_growth(WT, "WT", 600)
plot_growth(X_A, substitute(paste(italic("clpX"))), 600)
plot_growth(WT_OX, "WT oxacillin", 600)
plot_growth(X_OX_A, substitute(paste(italic("clpX"), " oxacillin")), 600)

mtext("Time (min)", side = 1, outer = TRUE, line = 1)
mtext("Number of cells", side = 2, outer = TRUE, line = 1)

dev.off()

# Cut data to include only the first 8 hours ----------------------------------------------------------------
WT %<>%
  filter(POSITION_T < 28800)

WT_OX %<>%
  filter(POSITION_T < 28800)

X_A %<>%
  filter(POSITION_T < 28800)

X_OX_A %<>%
  filter(POSITION_T < 28800)


# Calculate WT pedigree -----------------------------------------------------

wt <- calc_pedigree(WT, 12.91)

# There were a few errors in the tracking, and manual corrections were therefore necessary. At the time of analysis it was not super easy to do that in the tracking plugin itself so they are made in the following:

wt %<>%
  mutate(d1 = replace(d1, cell == 41, 39)) %>%
  mutate(d2 = replace(d2, cell == 41, 36)) %>%
  mutate(d1 = replace(d1, cell == 58, NA)) %>%
  mutate(d2 = replace(d2, cell == 58, NA)) %>%
  mutate(lifetime = replace(lifetime, cell == 58, 180)) %>%
  mutate(lifetime = replace(lifetime, cell == 41, 164)) %>%
  mutate(lifetime = replace(lifetime, cell == 120, 136)) %>%
  mutate(lifetime = replace(lifetime, cell == 190, 272)) %>%
  mutate(lifetime = replace(lifetime, cell == 155, 64)) %>%
  mutate(d1 = replace(d1, cell == 190, NA)) %>%
  mutate(d2 = replace(d2, cell == 190, NA)) %>%
  mutate(d1 = replace(d1, cell == 155, 122)) %>%
  mutate(d2 = replace(d2, cell == 155, 128)) %>%
  mutate(dead = replace(dead, cell == 120, FALSE))

# Calculate WT oxa pedigree -----------------------------------------------

# First I need to correct for some minor bug in calc_pedigree function that I can't seem to fix easily.

WT_OX %<>%
  mutate(generation = replace(generation, POSITION_T == 0 & TRACK_ID == 1, 211)) %>%
  mutate(generation = replace(generation, POSITION_T == 0 & TRACK_ID == 0, 212)) 

wt_ox <- calc_pedigree(WT_OX, 12.91)

# Corrections:

wt_ox %<>%
  mutate(d1 = replace(d1, cell == 93, 188)) %>%
  mutate(lifetime = replace(lifetime, cell == 188, 84)) 

# Calculate clpX pedigree -----------------------------------------------

x_a <- calc_pedigree(X_A, 12.91)

# Corrections:

x_a %<>%
  mutate(d2 = replace(d2, cell == 34, 36)) %>%
  mutate(lifetime = replace(lifetime, cell == 35, 12)) %>%
  mutate(sibling = replace(sibling, cell == 35, 36))

df <- x_a %>%
  filter(cell == 35) %>%
  mutate(cell = replace(cell, cell == 35, 36)) %>%
  mutate(sibling = replace(sibling, cell == 36, 35))

x_a <- rbind_list(x_a, df)

# Calculate clpX oxa pedigree -------------------------------------------

x_ox_a <- calc_pedigree(X_OX_A, 12.91)

# Convert pedigree to Newick format-------------------------------------------------------

wt_nwk_0 <- get_newick(wt, Track = 0)

wt_ox_nwk_0 <- get_newick(wt_ox, Track = 0)

wt_nwk_1<- get_newick(wt, Track = 1)

wt_ox_nwk_1 <- get_newick(wt_ox, Track = 1)

x_a_nwk_0 <- get_newick(x_a, Track = 0)

x_a_nwk_1 <- get_newick(x_a, Track = 1)

x_ox_a_nwk <- get_newick(x_ox_a, Track = 0)

# Create Supplemental Figure 2a -------------------------------------------

png("lineages_WT_1.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,0,2,0), cex = 0.7)

plot_lineage(wt, wt_nwk_0, "WT cell 1", y = 500 - root_edge_length(wt, 0), x = 80)

dev.off()

png("lineages_WT_2.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,0,2,0), cex = 0.7)

plot_lineage(wt, wt_nwk_1, "WT cell 2", y = 500 - root_edge_length(wt, 1), x = 80)

dev.off()

png("lineages_WT_3.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,0,2,0), cex = 0.7)

plot_lineage(wt_ox,wt_ox_nwk_0, "WT + OXA cell 1", y = 500 - root_edge_length(wt_ox, 0), x = 80)

dev.off()

png("lineages_WT_4.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,0,2,0), cex = 0.7)

plot_lineage(wt_ox, wt_ox_nwk_1, "WT + OXA cell 2", y = 500 - root_edge_length(wt_ox, 1), x = 80)

dev.off()


png("lineages_X_1.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,2,2,2),  cex = 0.7)

plot_lineage(x_a, x_a_nwk_0, "clpX movie 1, cell 1", x = 80, y = 500 - root_edge_length(x_a, 0))

dev.off()

png("lineages_X_2.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,2,2,2), cex = 0.7)

plot_lineage(x_a, x_a_nwk_1, "clpX movie 1, cell 2", x = 80, y = 500 - root_edge_length(x_a, 1))

dev.off()

png("lineages_X_3.png", width = 10, height = 10, units = "cm", res = 600)
par(mar = c(2,2,2,2), cex = 0.7)


plot_lineage(x_ox_a, x_ox_a_nwk, "clpX + OXA movie 3", x = 80, y = 500 - root_edge_length(x_ox_a, 0))

dev.off()

# Create Supplemental Figure 2b -------------------------------------------

png("all_gen_times_WT.png", width = 6, height = 6, units = "cm", res = 600)
par(mar = c(4,4,2,2), cex = 0.9)

plot_gen_times_3(wt, "WT", 12.91, 480)

dev.off()


png("all_gen_times_X.png", width = 6, height = 6, units = "cm", res = 600)
par(mar = c(4,4,2,2), cex = 0.9)

plot_gen_times_3(x_a, substitute(paste(italic("clpX"))), 12.91, 480)

dev.off()


png("all_gen_times_WT_OXA.png", width = 6, height = 6, units = "cm", res = 600)
par(mar = c(4,4,2,2), cex = 0.9)

plot_gen_times_3(wt_ox, "WT oxacillin", 12.91, 480)

dev.off()


png("all_gen_times_X_OXA.png", width = 6, height = 6, units = "cm", res = 600)
par(mar = c(4,4,2,2), cex = 0.9)

plot_gen_times_3(x_ox_a, substitute(paste(italic("clpX"), " oxacillin")), 12.91, 480)

dev.off()

# Summarize data for Supplemental Table 1 ---------------------------------------------
 
summary(wt)
summary(wt_ox)
summary(x_a)
summary(x_ox_a)
















