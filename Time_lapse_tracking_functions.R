### Code written by Kristoffer T. Baek, 2018 ###

### Functions used for analyzing tracking output from TrackMate (Fiji).

## Calculate for each spot whether it is a mother (links to two target spots in next frame), a daughter (linked to a mother in the previous frame), or neither. The function identifies all spots belonging to a cell, and assigns a unique cell ID number (called generation ID in the code) to each spot. It also parses additional information from the input spot file, such as spatial location in the frame.  

analyze_tracks <- function(links, spots) {
  
  # Remove irrelevant information
  spots %<>%
    select(ID, TRACK_ID, POSITION_T, POSITION_X, POSITION_Y, FRAME) 
  
  # Define mothers as source spots having two target spots
  mother <- links$SPOT_SOURCE_ID[duplicated(links$SPOT_SOURCE_ID)]
  
  # Define founders as the spots present in the first frame
  founder_cells <- spots %>% 
    filter(FRAME == 0) %>%
    select(ID)
  
  # Define the end of a cell generation 
  links %<>%
    select(SPOT_SOURCE_ID, SPOT_TARGET_ID) %>%
    set_colnames(c("source", "target")) %>%
    mutate(generation_end = ifelse(source %in% mother, TRUE, FALSE)) 
  
  # Define daugthers
  daughter <- links %>%
    filter(generation_end) %>%
    select(target) 
  
  daughter <- c(daughter$target, founder_cells$ID)
  
  # Define start of cell generations
  links %<>%
    mutate(generation_start = ifelse(source %in% daughter, TRUE, FALSE))
  
  if(links %>% filter(generation_end, generation_start) %>% nrow() > 0) {
    message("one frame generations")
  }
  
  # Assign a unique cell ID to each daughter spot
  cell_ID <- links %>%
    filter(!generation_end | !generation_start) %>%
    filter(generation_start) %>%
    mutate(generation = seq_len(length.out = nrow(.))) %>%
    select(source, generation)
  
  links %<>%
    join(cell_ID, by = c("source")) %>% tbl_df ## Non-daughter cells will for now have NAs as gen_ID
  
  i <- nrow(links)
  
  # Apply the 'grow' function to the data until no new generation-IDs are added
  
  repeat{
    
    links <- grow(links)
    
    n <- links %>% 
      filter(is.na(generation)) %>%
      nrow()
    
    if(n == i) break
    
    i <- n
    
  }
  
  # Print the number of unassigned spots to the console. Unassigned spots are usually noise (not cells) occurring for a few consecutive frames outside the microcolony
  message(paste(n, "out of", nrow(links), "spots could not be assigned to a cell generation", sep = " "))
  
  
  # Arrange data for output and include information from spotfile
  links %<>%  
    arrange(generation, source) %>% 
    rename(ID = source) %>% 
    join(spots, by = "ID") %>% tbl_df
  
  return(links)
  
  
}


## Grow each cell generation by one frame. This function takes the spots that already have an assigned cell-ID, and assigns the same ID to their targets. If the ID of a target is already known, no change is made. The function is used inside the 'analyze_tracks' function in a loop. 

grow <- function(data) {
  
  cell_ID <- data %>%
    filter(!is.na(generation)) %>% 
    select(target, generation) %>%
    rename(source = target,
           y = generation) 
  
  data %<>%
    join(cell_ID, by = c("source")) %>%
    rename(x = generation) %>%
    mutate(generation = ifelse(is.na(x),
                               ifelse(!is.na(y), y, x),
                               x)
    ) %>%
    select(-x, -y)
  
  return(data)
}


## Plot the number of cells over time

plot_growth <- function(data, title, maxx) {
  
  x <- data %>% 
    group_by(POSITION_T) %>%
    summarise(n = n())
  
    plot(x$POSITION_T/60, x$n,
         pch = 21,
         ylim = c(0,170),
         xlim = c(0,maxx),
         xlab = "Time (min)",
         ylab = "Number of cells",
         main = title)
 
  box(which = "plot")
  
}


## Calculate the pedigree relationship between all cells. The function also outputs information such as generation time, whether the cell dies etc. 

calc_pedigree <- function(data, dim) {
  
  # Define timepoint for end of experiment
  EOE <- max(data$POSITION_T) 
  
  x <- data %>%
    # Calculate if a cell is located near the edge of the frame
    mutate(edge = ifelse(POSITION_X > dim * 0.9 |
                           POSITION_X < dim * 0.1 |
                           POSITION_Y > dim * 0.9 |
                           POSITION_Y < dim * 0.1, 
                         TRUE, FALSE)) %>%
    group_by(generation) %>%
    summarise(start_T = min(POSITION_T), #time of birth of cell
              end_T = max(POSITION_T), #time of end of cell generation
              duration = (end_T - start_T)/60, #generation time in minutes
              generation_end = unique(generation_end[POSITION_T == end_T]), #does the generation end by dividing?
              generation_start = unique(generation_start[POSITION_T == start_T]), #does the generation start as a daughter cell (including founder cells)?
              dead = ifelse(end_T != max(data$POSITION_T)  & !generation_end, TRUE, FALSE),
              outside = unique(edge[POSITION_T == end_T]),
              first_cell = unique(ID[POSITION_T == start_T]), #ID of first spot
              last_cell = unique(ID[POSITION_T == end_T]), #ID of last spot
              daughter_1 = min(target[ID == last_cell]),
              daughter_2 = max(target[ID == last_cell]),
              end_of_exp = ifelse(end_T == EOE, TRUE, FALSE), #cell present at the EOE?
              origin = ifelse(start_T == 0, TRUE, FALSE),  #cell present at the beginning?
              track = unique(TRACK_ID)
    ) 
  
  # Create data frame with daughter spot IDs for each cell.
  z <- x %>%
    select(generation, first_cell, last_cell, daughter_1, daughter_2) %>%
    gather(d, ID, daughter_1, daughter_2) %>% arrange(generation)
  
  # Create data frame with first spot of each cell generation and rename columns for later use
  y <- z %>%
    select(generation, first_cell) %>%
    mutate(ID = first_cell,
           daughter_gen = generation) %>%
    select(-generation, -first_cell) %>%
    distinct()
  
  # Combine the two above data frames to establish link between cell generation IDs and their daughter generation IDs
  w <- z %>%
    join(y, by = "ID") %>%
    select(generation, daughter_gen) 
  
  
  # Add this information to the full data frame 
  p <- x %>% select(generation, duration, dead, end_of_exp, outside, origin, track, start_T, end_T) %>%
    join(w, by = "generation") %>%
    rename(cell = generation,
           daughter_cell = daughter_gen,
           lifetime = duration) %>%
    mutate(n = rep(c("d1","d2"),n()/2)) %>%
    spread(n,daughter_cell)
  
  # Create data frame of sibling pairs
  a1 <- p %>%
    select(d1, d2) %>%
    filter(!is.na(d1)) %>%
    rename(cell = d1,
           sibling = d2)
  
  # Create data frame of sibling pairs but with reverse column naming
  a2 <- p %>%
    select(d1, d2) %>%
    filter(!is.na(d1)) %>%
    rename(cell = d2,
           sibling = d1) 
  
  # Combine these two data frames to get a data frame giving the sibling cell for every cell
  a <- rbind_list(a1, a2) %>%
    distinct
  
  # Create data frame giving the mother cell for every cell
  m <- p %>% 
    gather(x, daughter, d1, d2) %>%
    select(cell, daughter) %>%
    rename(mother = cell,
           cell = daughter) %>%
    select(cell, mother) %>%
    filter(!is.na(cell)) %>%
    arrange(cell)
  
  
  # Add mother and sibling cell information to the full data frame
  p %<>%
    join(a, by = "cell") %>%
    join(m, by = "cell")
  
  return(p)
  
}


## The following functions converts the calculated pedigree to Newick format, which can be used to plot simple good looking lineage trees using the ape package. 

get_text <- function(p, i, frametime)  {
  
  A <- p %>%
    filter(cell == i) %>%
    select(cell) %>%
    as.numeric()
  
  B <- p %>%
    filter(cell == i) %>%
    select(lifetime) %>%
    as.numeric()
  
  B <- B + frametime
  
  text <- paste(A, ":", B, sep= "")
  
  return(text)
  
}

get_sister <- function(p, i) {
  
  C <- p %>%
    filter(cell == i) %>%
    select(sibling) %>%
    as.numeric()
  
  return(C)
  
}

get_mother <- function(p, i) {
  
  D <- p %>%
    filter(cell == i) %>%
    select(mother) %>%
    as.numeric()
  
  return(D)
  
}

get_daughter <- function(p, i) {
  
  d1 <- p %>%
    filter(cell == i) %>%
    select(d1) %>%
    as.numeric()
  
  d2 <- p %>%
    filter(cell == i) %>%
    select(d2) %>%
    as.numeric()
  
  return(c(d1, d2))
  
}

siblings <- function(p, i, frametime) {
  
  paste(get_text(p, i, frametime), ",", get_text(p, get_sister(p, i), frametime), sep = "")
  
}

branch_down <- function(p, i, frametime) {
  
  paste("(", siblings(p, get_daughter(p, i)[1], frametime), ")", get_text(p, i, frametime), ";", sep = "")
  
}

df2newick_2 <- function(p, vector, string, frametime) {
  
  c <- vector[1]
  
  d <- get_daughter(p, c)
  
  if(!any(is.na(d))) {
    
    d1 <- get_text(p, d[1], frametime)
    d2 <- get_text(p, d[2], frametime)
    
    sib <- paste("(", d1, ",", d2, ")", sep = "")
    
    newick <- sub(paste("(.*\\D)(", c ,":.*)", sep = ""), paste("\\1", sib, "\\2", sep = ""), x = string) #replace text after (
    
    
  }else{
    
    newick <- string
    
  }
  
  
  
  vector <- vector[-1]
  
  if(any(!is.na(d)))  vector <- c(vector, d)
  
  
  
  
  return(list(vector, newick))
  
}

get_newick <- function(p, Track = 0, frametime = 4) {
  
  p %<>%
    filter(track == Track)
  
  ori <- p %>%
    filter(origin) %>%
    select(cell) %>%
    as.numeric()
  
  newick <- df2newick_2(p, get_daughter(p, ori), branch_down(p, ori, frametime), frametime)
  
  
  
  repeat{
    
    newick <- df2newick_2(p, newick[[1]], newick[[2]], frametime)
    
    if(length(newick[[1]]) == 0) break
    
  }
  
  return(newick)
  
}


## Plot lineage trees

root_edge_length <- function(p, Track) {
  
  p %<>% filter(origin, track == Track) %>%
    select(lifetime) %>% 
    as.numeric()
  
  return(p)
  
}

plot_lineage <- function(p, nwk, title, x = 50, y = 600) {
  
  dead <- p %>% filter(dead) %>% select(cell) 
  dead <- dead$cell
  tip_labels <- as.numeric(collapse.singles(read.tree(text = nwk[[2]]))$tip.label)
  
  plot.phylo(collapse.singles(read.tree(text = nwk[[2]])), 
             show.tip.label= FALSE, 
             root.edge = TRUE, 
             direction = "downwards",
             main = "",
             x.lim = c(0, x),
             y.lim = c(0, y)
  )
  
  tiplabels(pch = ifelse(tip_labels %in% dead, 19, NA), 
            col = "red",
            bg = "red"
            
  )
  
}


## Plot generation times as function of elapsed time

plot_gen_times_3 <- function(data, title, dim, maxx) {
  
  s <- data %>% filter(!dead, !is.na(d2), start_T != 0)
  
  plot(s$start_T/60, s$lifetime, 
       ylim = c(0,400), 
       pch = 21,
       col= rgb(red = 0, green = 0, blue = 0.0, alpha = 1),
       bg = rgb(red = 0, green = 0, blue = 0.0, alpha = 0.4),
       xlim = c(0,maxx),
       xlab = "Time (min)",
       ylab = "Time between divisions (min)",
       main = title
  )
  
  box(which = "plot")
  
}


## Compile summary of generation times for divisions taking place in 1 h intervals
summary <- function(data) {
  
  data_divisions <- data %>% filter(!dead, !is.na(d2), start_T != 0)
  
  data_divisions %<>% mutate(start_T_min = start_T/60)
  
  data_divisions %<>% mutate(interval = ifelse(start_T_min > 420, 8, 8),
                             interval = ifelse(start_T_min > 360 & start_T_min < 421, 7, interval),
                             interval = ifelse(start_T_min > 300 & start_T_min < 361, 6, interval),
                             interval = ifelse(start_T_min > 240 & start_T_min < 301, 5, interval),
                             interval = ifelse(start_T_min > 180 & start_T_min < 241, 4, interval),
                             interval = ifelse(start_T_min > 120 & start_T_min < 181, 3, interval),
                             interval = ifelse(start_T_min > 60 & start_T_min < 121, 2, interval),
                             interval = ifelse(start_T_min < 60, 1, interval))
  
  
  summary <- data_divisions %>% group_by(interval) %>%
    summarize(number_of_divisions = n(),
              mean_generation_time_min = round(mean(lifetime),0),
              SD = round(sd(lifetime),1)) %>%
    rename(hour_interval = interval) %>%
    rbind(list("All_cells", 
               nrow(data_divisions),
               round(mean(data_divisions$lifetime), 0), 
               round(sd(data_divisions$lifetime),1)) )
  
  return(summary)
  
}