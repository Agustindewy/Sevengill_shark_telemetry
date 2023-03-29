
# De Wysiecki et al. Habitat use patterns of broadnose sevengill sharks (Notorynchus cepedianus) in Caleta Valdés

# Analyses and figures
# Each figure is self-sufficient, so expect high repetition of lines in script
library(raster)
library(pheatmap)
library(viridis)
library(circular)
library(ggplot2)
library(rgdal)
library(RSP)
library(actel)
library(gridExtra)
library(geosphere)
library(suncalc)
library(ggsn)
library(tidyverse)
library(ggforce)
library(lubridate)
library(mgcv)
library(MuMIn)
library(lme4)
library(tidymv)
library(investr) 

setwd('SET YOUR WORKING DIRECTORY')

# Study area coastline (spatial polygons) both at high (18-09-2013) and low (15-10-2018) tides
coastline <- 'DOWNLOAD, READ AND SET PATH' # Available at Dryad repository
coast <- readOGR(dsn = coastline, layer = 'Caleta_marea_alta')

# Raw detection data
dat_raw <- read.csv('Detections.csv') # Available at Dryad repository
dat_raw$Long_date <- as.POSIXct(dat_raw$Long_date, format = '%d/%m/%Y %H:%M', tz = 'UTC')
dat_raw$Date <- date(dat_raw$Long_date)
dat_raw$Year <- year(dat_raw$Long_date)
dat_raw$Month <- month(dat_raw$Long_date)
dat_raw$Day <- day(dat_raw$Long_date)
dat_raw$Hour <- hour(dat_raw$Long_date)
dat_raw <- dat_raw %>% filter(Full_tag != 'A69-9001-5966') %>% # tag detached from shark (constantly detected throughout the study period)
  filter(Date >= as.Date('2019-11-01')) # filter tagging days at the beginning

# Daily detection data
dat_day <- dat_raw %>% group_by(Date, Year, Month, Tag) %>%
  summarize(Raw.det = length(Long_date), Roam = length(unique(Receiver)), Sex = unique(Sex), TL = unique(TL_cm)) %>%
  mutate(Binary.det = ifelse(Raw.det > 1, 1, 0)) %>% filter(Binary.det > 0)

# Hourly detection data
dat_hour <- dat_raw %>% group_by(Date, Year, Month, Hour, Tag) %>%
  summarize(Raw.det = length(Long_date), Roam = length(unique(Receiver)), Sex = unique(Sex), TL = unique(TL_cm))
dat_hour <- dat_hour %>% right_join(dat_day[, c('Date', 'Tag')], by = c('Date', 'Tag')) %>% mutate(Hour = as.integer(Hour))


#--------------------------------- Figure 1 ----------------------------------

# Low-high tide contrast
coast1 <- readOGR(dsn = coastline, layer = 'Caleta_marea_baja') # Available at Dryad repository
coast2 <- readOGR(dsn = coastline, layer = 'Caleta_marea_alta') # Available at Dryad repository
df_exclude <- data.frame(xmin = -63.66, xmax = -63.72, ymin = -42.45, ymax = -42.60)

g1 <- ggplot()+
  geom_polygon(data = coast1, aes(x = long, y = lat, group = group), col = 'white', size = 0.01) +
  coord_equal(xlim = c(-63.6875, -63.595), ylim = c(-42.51, -42.2575), expand = 0) + 
  ggtitle('Low tide') + xlab('Longitude (W)') + ylab('Latitude (S)') +
  geom_rect(data = df_exclude, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey20') +
  scale_y_continuous(name = NULL, breaks = c(-42.49, -42.42, -42.35, -42.28), labels = c('42.49°', '42.42°', '42.35°', '42.28°')) +
  scale_x_continuous(name = NULL, breaks = c(-63.66, -63.62), labels = c('63.66°', '63.62°')) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5),
        plot.title = element_text(size = 12), 
        axis.title = element_text(size = 10)) 

g2 <- ggplot()+
  geom_polygon(data = coast2, aes(x = long, y = lat, group = group), col = 'white', size = 0.01) +
  coord_equal(xlim = c(-63.6875, -63.595), ylim = c(-42.51, -42.2575), expand = 0) + 
  ggtitle('High tide') + xlab('Longitude (W)') + ylab('Latitude (S)') +
  geom_rect(data = df_exclude, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey20') +
  scale_y_continuous(name = NULL, breaks = c(-42.49, -42.42, -42.35, -42.28), labels = c('42.49°', '42.42°', '42.35°', '42.28°')) +
  scale_x_continuous(name = NULL, breaks = c(-63.66, -63.62), labels = c('63.66°', '63.62°')) +
  north(location = 'topright', symbol = 3, x.min = -63.6875, x.max = -63.598, y.min = -42.51, y.max = -42.255) +
  scalebar(x.min = -63.6875, x.max = -63.595, y.min = -42.51, y.max = -42.2575, transform = T, 
                 dist = 1, st.size = 3, height = 0.01, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.6, anchor = c(x = -63.66, y = -42.5)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5),
        plot.title = element_text(size = 12), 
        axis.title = element_text(size = 10)) 

ggg <- grid.arrange(g1, g2, nrow = 1)
ggsave('Figure 1.pdf', ggg, width = 13, height = 15, units = 'cm')


#--------------------------------- Table 1 ----------------------------------

# First and last detection dates and detection duration for each shark
first_last <- data.frame()
for(i in sort(unique(dat_day$Tag))){
  sub <- subset(dat_day, Tag == i)
  first_last <- rbind(first_last, data.frame(ID = i, first = min(sub$Date), last = max(sub$Date)))
} # tagging date coincides with first detection date in all cases
first_last$duration <- (first_last$last - first_last$first) + 1 # detection duration

# Number of days each shark was detected more than once in the study period
num_days <- table(dat_day$Tag, dat_day$Binary.det)

# Maximum detection window (consecutive days) until last detection
study_period <- seq(as.Date(min(dat_day$Date)), as.Date(max(dat_day$Date)), 1)

all_tags <- data.frame()
for(i in sort(unique(dat_day$Tag))){
  sub <- data.frame(Tag = rep(i, length(study_period)), Date = study_period)
  all_tags <- rbind(all_tags, sub)
}
all_tags$Binary.det <- 0

max_merge <- rbind(dat_day[, c('Tag', 'Date', 'Binary.det')], all_tags)
max_aggr <- aggregate(Binary.det ~ Tag + Date, data = max_merge, FUN = sum)
max_detection <- max_aggr[which(max_aggr$Binary.det < 1), ]

max_detection_df <- data.frame()
for(i in sort(unique(max_detection$Tag))){
  sub <- subset(max_detection, Tag == i)
  max_det <- vector()
  for(j in 1:dim(sub)[1]){
    max_det <- c(max_det, (sub[j, ]$Date - sub[j - 1, ]$Date - 1))
  }
  max_det <- max(max_det)
  max_detection_df <- rbind(max_detection_df, data.frame(Tag = i, Value = max_det))
}

# Maximum non-detection window (consecutive days) until last detection
max_detection_df <- data.frame()
for(i in sort(unique(dat_day$Tag))){
  sub <- subset(dat_day, Tag == i)
  max_det <- vector()
  for(j in 1:dim(sub)[1]){
    max_det <- c(max_det, (sub[j, ]$Date - sub[j - 1, ]$Date - 1))
  }
  max_det <- max(max_det)
  max_detection_df <- rbind(max_detection_df, data.frame(Tag = i, Value = max_det))
}

# Number of total detections for each shark
df <- data.frame(table(dat_raw$Tag))
names(df) <- c('ID', 'N_detections')

# Residency index for each shark based on full study period (minimum residency)
round(rowSums(num_days) / as.numeric((max(dat_day$Date) - min(dat_day$Date) + 1)), 3)

# Residency index for each shark based on first and last detection (maximum residency)
rmax <- data.frame()
for(i in sort(unique(dat_day$Tag))){
  sub <- subset(dat_day, Tag == i)
  num_days0 <- table(sub$Tag, sub$Date)
  value <- round(rowSums(num_days0) / as.numeric((max(sub$Date) - min(sub$Date) + 1)), 3)
  rmax <- rbind(rmax, data.frame(rmax = value))
}

# Mean daily Roaming index
roam_df <- data.frame()
for(i in sort(unique(dat_day$Tag))){
  sub <- subset(dat_day, Tag == i)
  roam_mean <- mean(sub$Roam / 5)
  roam_sd <- sd(sub$Roam / 5)
  roam_df <- rbind(roam_df, data.frame(Tag = i, Roaming.mean = round(roam_mean, 3), Roaming.sd = round(roam_sd, 3)))
}

# Minimum distance travelled and area used by each shark 
# Refining the Shortest Paths (RSP) and Actel packages --> read guide in https://github.com/YuriNiella/RSP
setwd('~/Actel')

# Create water raster of study area (NA for land and open sea, 1 for study area)
Ext <- extent(-63.69, -63.595, -42.42, -42.26)
efective_coast <- crop(coast, Ext)
writeOGR(efective_coast, '.', 'Study_area', driver = 'ESRI Shapefile')
coast_raster <- loadShape(shape = 'Study_area.shp', size = 0.0002, buffer = 0.07, type = 'land', coord.x = 'Longitude', coord.y = 'Latitude')

# Shortest in-water paths
trans <- transitionLayer(coast_raster, directions = 16)
run1 <- explore(tz = 'America/Buenos_Aires', report = F)
rsp_run <- runRSP(input = run1, t.layer = trans, coord.x = 'Longitude', coord.y = 'Latitude', distance = 250,
                  time.step = 10, min.time = 10, max.time = 24)

# RSP minimum distance travelled by each shark
rsp_distances <- getDistances(rsp_run)
distances_df <- data.frame()
for(i in unique(rsp_distances$Animal.tracked)){
  sub <- subset(rsp_distances, Animal.tracked == i & Loc.type == 'RSP')
  distance_value <- sum(sub$Dist.travel)
  distances_df <- rbind(distances_df, cbind(shark_id = i, distance_travelled = round(distance_value/1000, 1)))
}

# RSP area used by each shark
area_df <- data.frame()
for(i in rsp_run$bio$Transmitter){
  dbbmm_run <- dynBBMM(input = rsp_run, base.raster = coast_raster, UTM = 20, tags = i,
                       start.time = '2019-10-29 17:36:00', stop.time = '2021-03-12 11:28:00')
  rsp_areas <- getAreas(input = dbbmm_run, type = 'group', breaks = c(0.5, 0.95))
  area_df <- rbind(area_df, cbind(shark_id = substr(i, 0, 13), area_used_0.5 = round(rsp_areas$areas$area.5 / 1000 ^ 2, 1),
                                  area_used_0.95 = round(rsp_areas$areas$area.95 / 1000 ^ 2, 1)))
}
area_df


#--------------------------------- Figure 2 ----------------------------------

# Monthly distribution of sharks in the array per receiver 

dat <- dat_raw %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != 5973) %>% # Tag 5973 was detected in 5 or less days in total so it is not considered
  group_by(Date, Receiver, Month, Tag) %>% summarize(dummie = 1)

dat_day <- dat_day %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != 5973) # Tag 5973 was detected in 5 or less days in total so it is not considered
  
dat <- dat %>% right_join(dat_day, by = c('Date', 'Month', 'Tag')) %>% select(Date, Month, Receiver, Tag, dummie)

Sitas <- dat %>% group_by(Date, Month, Receiver) %>% summarize(Sita = length(unique(Tag))) # Sharks In The Array

dat <- dat %>% group_by(Date, Month, Tag) %>% summarize(dummie = 1) %>% select(Date, Month, Tag, dummie)

my.fun <- function(X) {data.frame(Tag = X[1], Date = seq(ymd(X[2]), ymd(X[3]), by = 'day'))} # function to create day sequence until last detection
dat_dates <- dat %>% group_by(Tag) %>% summarize(minDate = '2019-11-01', maxDate = max(Date)) 
seq_list <- apply(X = data.frame(tag = dat_dates$Tag, first = dat_dates$minDate, last = dat_dates$maxDate), 1, FUN = my.fun) 
seq_df <- bind_rows(seq_list)
seq_df$Month <- month(seq_df$Date)

final_df <- seq_df %>% mutate(Tag = as.integer(Tag)) %>% left_join(dat, by = c('Tag', 'Date', 'Month')) 
final_df[is.na(final_df)] <- 0 # NAs replacement by zeros
final_df <- final_df %>% group_by(Date) %>% summarize(nTags = length(dummie)) %>% right_join(Sitas, by = 'Date') 
final_df$Proportion <- round(final_df$Sita / final_df$nTags * 100, 1)
final_df <- final_df %>% group_by(Month, Receiver) %>% summarize(Mean = mean(Proportion), SD = sd(Proportion)) %>%
  mutate(Month = as.factor(Month), Receiver = as.factor(Receiver))
final_df[is.na(final_df)] <- 0 # NAs replacement by zeros

# Bubble plot
ggplot(data = final_df, aes(x = Month, y = Receiver)) + 
  geom_point(aes(size = Mean, color = Mean)) + 
  scale_colour_viridis_c() +
  scale_y_discrete(name = 'Receiver', breaks = c(6:10), labels = c('6', '7', '8', '9', '10')) +
  scale_x_discrete(name = 'Month', breaks = c(1:12), 
                   labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure 2.pdf', width = 12, height = 8.6, units = 'cm')


#--------------------------------- Figure 3 ----------------------------------

# Read raw detection data again to include all days
dat_raw <- read.csv('Detections.csv') # Available at Dryad repository
dat_raw$Date <- as.Date(dat_raw$Date, format = '%d/%m/%Y')
dat_raw <- dat_raw %>% filter(Full_tag != 'A69-9001-5966') # tag detached from shark (constantly detected throughout the study period)

# Percentage of tagged sharks detected
agg_monthly <- aggregate(Tag ~ Year + Month, data = dat_raw, FUN = unique)
agg_monthly$Sum <- lapply(agg_monthly$Tag, FUN = length)  
agg_monthly$Sum <- unlist(agg_monthly$Sum)
agg_monthly$Percentage <- agg_monthly$Sum * 100 / length(unique(dat_raw$Tag))
agg_monthly <- agg_monthly[with(agg_monthly, order(Year, Month)), ]
agg_monthly$Date <- seq.Date(as.Date('2019-10-01'), as.Date('2021-03-01'), by = 'month')

ggplot(agg_monthly, aes(x = Date, y = Percentage)) + 
  geom_line(color = '#FDE725FF') + geom_point(color = '#440154FF') + ylim(0, 100) + 
  xlab('Date') + ylab('% of sharks detected') +
  scale_x_date(name = NULL, breaks = seq(as.Date('2019-10-01'), as.Date('2021-03-01'), 'months'), limits = as.Date(c('2019-10-01', '2021-03-12')), 
               labels = c('', 'Nov 2019', '', '', 'Feb 2020', '', '', 'May 2020', '', '', 'Aug 2020', '', '', 'Nov 2020', '', '', 'Feb 2021', '')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA)) 
ggsave('Figure 3a.pdf', width = 12, height = 4, units = 'cm')

# Detection plot at all receivers
dat_raw$Dummie <- 1
agg_daily <- aggregate(Dummie ~ Date + Tag + Sex, data = dat_raw, FUN = sum)
agg_daily$Tag <- as.factor(agg_daily$Tag)
agg_daily$Sex <- as.factor(agg_daily$Sex)
agg_daily <- agg_daily[with(agg_daily, order(Sex)), ]
agg_daily$Tag <- factor(agg_daily$Tag, levels = rev(unique(agg_daily$Tag[order(agg_daily$Sex)])))

last_detec <- data.frame()
for(i in sort(unique(agg_daily$Tag))){
  sub <- subset(agg_daily, Tag == i)
  last_detec <- rbind(last_detec, data.frame(ID = i, last = max(sub$Date)))
} # tagging date coincides with first detection date in all cases

ggplot() +
  geom_point(data = agg_daily, aes(x = Date, y = Tag, color = Sex), size = 0.5) + 
  geom_point(data = last_detec, aes(x = last, y = ID), size = 1.25, col = '#21908CFF', shape = 8) +
  scale_colour_viridis_d() + xlab('Date') + ylab('Tag ID') +
  scale_x_date(name = NULL, breaks = seq(as.Date('2019-10-01'), as.Date('2021-03-01'), 'months'), limits = as.Date(c('2019-10-01', '2021-03-12')), 
               labels = c('', 'Nov 2019', '', '', 'Feb 2020', '', '', 'May 2020', '', '', 'Aug 2020', '', '', 'Nov 2020', '', '', 'Feb 2021', '')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA),
        legend.position = 'none') 
ggsave('Figure 3b.pdf', width = 12, height = 8, units = 'cm')


#--------------------------------- Figure 4 ----------------------------------

dat <- dat_day %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != 5973) # # Tag 5973 was detected in the first 5 or less days in total so it is not considered

# Seasonal residency and roaming indices for each shark (minimum seasonal residency)
dat$Season <- as.character('') # Set season column based on peak of abundance season (October-February period)
dat$Season[which(dat$Month %in% c(11, 12, 1, 2))] <- 'High season 19/20'
dat$Season[which(dat$Month %in% c(3:9))] <- 'Low season 20'
dat$Season[which(dat$Month %in% c(10, 11, 12) & dat$Year %in% 2020)] <- 'High season 20/21'
dat$Season[which(dat$Month %in% c(1, 2) & dat$Year %in% 2021)] <- 'High season 20/21'

# Minimum residency
residency_df <- data.frame()
sub <- subset(dat, Season == 'High season 19/20')
num_days <- table(sub$Tag, sub$Date) 
residency0 <- round(rowSums(num_days) / 120, 3)
residency_df0 <- data.frame(Tag = names(residency0), Season = 'High season 19/20', Residency.min = residency0)
mean_roam <- aggregate(Roam / 5 ~ Tag, data = sub, FUN = max) # Mean Roaming index of each shark detected in the season
residency_df0$Roaming <- round(mean_roam$`Roam/5`, 3) # Mean daily Roaming index only considering detection days
residency_df <- rbind(residency_df, residency_df0)

sub <- subset(dat, Season == 'Low season 20')
num_days <- table(sub$Tag, sub$Date) 
residency0 <- round(rowSums(num_days) / 213, 3)
residency_df0 <- data.frame(Tag = names(residency0), Season = 'Low season 20', Residency.min = residency0)
mean_roam <- aggregate(Roam / 5 ~ Tag, data = sub, FUN = max) # Mean Roaming index of each shark detected in the season
residency_df0$Roaming <- round(mean_roam$`Roam/5`, 3) # Mean daily Roaming index only considering detection days
residency_df <- rbind(residency_df, residency_df0)

sub <- subset(dat, Season == 'High season 20/21')
num_days <- table(sub$Tag, sub$Date) 
residency0 <- round(rowSums(num_days) / 150, 3)
residency_df0 <- data.frame(Tag = names(residency0), Season = 'High season 20/21', Residency.min = residency0)
mean_roam <- aggregate(Roam / 5 ~ Tag, data = sub, FUN = max) # Mean Roaming index of each shark detected in the season
residency_df0$Roaming <- round(mean_roam$`Roam/5`, 3) # Mean daily Roaming index only considering detection days
residency_df <- rbind(residency_df, residency_df0)

residency_df <- transform(residency_df, Season = factor(Season, levels = c('High season 19/20', 'High season 20/21', 'Low season 20')))

# Plot
ggplot(data = residency_df) +
  geom_vline(xintercept = c(0.1, 0.5), linetype = 'dashed', color = 'grey20', size = 0.4) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey20', size = 0.4) +
  geom_point(aes(x = Residency.min, y = Roaming, color = Season, fill = Season, size = Residency.min), shape = 21) +
  geom_smooth(aes(x = Residency.min, y = Roaming, color = Season, fill = Season), se = F, span = 5, size = 1) +
  scale_colour_viridis_d() + scale_fill_viridis_d() +
  xlab('Residency index') + ylab('Roaming index') + xlim(0, 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3)) +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2), labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'))
ggsave('Figure 4.pdf', width = 15, height = 9.6, units = 'cm')


#--------------------------------- Data preparation for figures 5 & 6 ------------------------------

# Create daily detection data

dat <- dat_day %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != 5973) # Tag 5973 was detected in 5 or less days in total so it is not considered

# Detection data
my.fun <- function(X) {data.frame(Tag = X[1], Date = seq(ymd(X[2]), ymd(X[3]), by = 'day'))} # function to create day sequence until last detection
dat_dates <- dat %>% group_by(Tag, Sex, TL) %>% summarize(minDate = '2019-11-01', maxDate = max(Date)) 
seq_list <- apply(X = data.frame(tag = dat_dates$Tag, first = dat_dates$minDate, last = dat_dates$maxDate), 1, FUN = my.fun) 
seq_df <- bind_rows(seq_list)
seq_df$Month <- month(seq_df$Date)
seq_df$Julian <- yday(seq_df$Date) # Julian day
seq_df$Year <- year(seq_df$Date) # Julian day

Sitas <- dat %>% group_by(Date) %>% summarize(Sita = length(unique(Tag))) # Sharks In The Array
Fitas <- dat %>% group_by(Date) %>% filter(Sex == 'female') %>% summarize(Fita = length(unique(Tag))) # Females In The Array
Mitas <- dat %>% group_by(Date) %>% filter(Sex == 'male') %>% summarize(Mita = length(unique(Tag))) # Males In The Array

model_input <- seq_df %>% mutate(Tag = as.integer(Tag)) %>%
  left_join(as.data.frame(dat)[, -3], by = c('Tag', 'Date')) %>%
  select(Tag, Date, Julian, Year.x, Month, Raw.det, Binary.det, Roam) %>%
  left_join(as.data.frame(dat_dates[, c('Tag', 'Sex', 'TL')]), by = 'Tag') %>%
  mutate(Tag = as.factor(Tag), Sex = as.factor(Sex), TL = as.numeric(TL)) %>%
  left_join(as.data.frame(Sitas), by = 'Date') %>%
  left_join(as.data.frame(Fitas), by = 'Date') %>%
  left_join(as.data.frame(Mitas), by = 'Date') %>%
  mutate(Season = cut(Date, breaks = as.Date(c('2019-11-01', '2020-03-01', '2020-10-01', '2021-03-12')),
                      labels = c('high season', 'low season', 'high season'),
                      include.lowest = T, right = F))

colnames(model_input)[4] <- 'Year'
model_input[is.na(model_input)] <- 0 # NAs replacement by zeros
model_input$fail.Roam <- 5 - model_input$Roam # number of non-visited receivers 

# Proportion of sharks, females and males present relative to the number of sharks still active each day, i.e. have not being last detected yet
Prop_sharks <- model_input %>% group_by(Date) %>% summarize(Active_sharks = length(unique(Tag)))
Prop_females <- model_input %>% group_by(Date) %>% filter(Sex == 'female') %>% summarize(Active_females = length(unique(Tag)))
Prop_males <- model_input %>% group_by(Date) %>% filter(Sex == 'male') %>% summarize(Active_males = length(unique(Tag)))
model_input <- model_input %>% left_join(as.data.frame(Prop_sharks), by = 'Date')
model_input <- model_input %>% left_join(as.data.frame(Prop_females), by = 'Date')
model_input <- model_input %>% left_join(as.data.frame(Prop_males), by = 'Date')
model_input$Prop_sharks <- model_input$Sita / model_input$Active_sharks  
model_input$Prop_females <- model_input$Fita / model_input$Active_females  
model_input$Prop_males <- model_input$Mita / model_input$Active_males 

# Predictors 
study_period <- seq(as.Date(min(model_input$Date)), as.Date(max(model_input$Date)), 1)

# Sea Surface Temperature (°C) - downloaded 10-2019 to 03-2021 daily, ~4 km
sst <- read.csv('Predictors/sst_2019-2021.csv', skip = 2, header = F) # Available at Dryad repository
names(sst) <- c('Date', 'Latitude', 'Longitude', 'SST')
sst$Date <- as.Date(substr(sst$Date, 0, 10), format = '%Y-%m-%d')
grid <- SpatialPoints(expand.grid(unique(sst$Longitude), unique(sst$Latitude), KEEP.OUT.ATTRS = F), proj4string = CRS(proj4string(coast))) # Determine some nearby water area to consider for daily sst description
longs <- c(-63.56250, -63.52083, -63.47917, -63.43750, -63.39583, -63.35417, -63.31250, -63.27083, -63.22917, -63.18750, -63.14583) # chosen longitudes
lats <- c(-42.18750, -42.22918, -42.27083, -42.31250, -42.35418, -42.39583, -42.43750, -42.47918, -42.52083) # chosen latitudes
sst_subset <- subset(sst, Longitude %in% longs) # Subset desired area
sst_subset <- subset(sst_subset, Latitude %in% lats)
sst_final <- aggregate(SST ~  Date, data = sst_subset, FUN = mean, na.action = na.omit) # Get a mean sst daily value for this area
smooth_model <- smooth.spline(sst_final, spar = 0.75) # Create a simple model of SST to get estimates of SST in 'NA' days due to cloud cover
sst_prediction <- predict(smooth_model, as.numeric(study_period))
predictor_df <- data.frame(Date = as.Date(study_period, format = '%Y-%m-%d'), SST = round(sst_prediction$y, 2))

# Tide amplitude
tide_amplitude <- read.table('Predictors/Comodoro_2019-10-01_2021-03-31.txt') # Available at Dryad repository
tide_amplitude <- data.frame(Date = tide_amplitude$V1, Tide = tide_amplitude$V5)
tide_amplitude$Date <- as.Date(tide_amplitude$Date)
tide_high <- aggregate(Tide ~ Date, data = tide_amplitude, FUN = max)
tide_low <- aggregate(Tide ~ Date, data = tide_amplitude, FUN = min)
tide_amplitude <- data.frame(Date = tide_high$Date, Amplitude = tide_high$Tide - tide_low$Tide)
tide_amplitude <- tide_amplitude %>% filter(Date >= as.Date('2019-11-01') & Date <= as.Date('2021-02-28'))
predictor_df$TA <- tide_amplitude$Amplitude

# Paste final set of predictors to model input
model_input <- model_input %>% left_join(predictor_df, by = c('Date' = 'Date'))

source('HighstatLib.R') # for VIF multicollinearity analysis

# Prdictor correlacion
vars <- c('Julian', 'TL', 'SST', 'TA', 'Prop_females', 'Prop_males')
corvif(model_input[, vars])
pairs(model_input[, vars], lower.panel = panel.cor, cex.labels = 1, labels = vars)

write.csv(model_input, 'Daily_data.csv', row.names = F)


#--------------------------------- Figures 5, S2a and S2b ----------------------------

# Residency index - detection probability

model_input <- read.csv('Daily_data.csv')
model_input <- model_input %>% mutate(Sex = factor(Sex), Tag = factor(Tag), Year = factor(Year), Date = as.Date(Date))

# Do the tags as random effects improve the fit?
Mod1 <- gam(Binary.det ~ Sex + s(TL)+ s(Month, bs = 'cc', k = 4) + s(Month, bs = 'cc', by = Sex, k = 4) + s(TA) + 
              s(SST, bs = 'cc', k = 4), select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod1)
Mod2 <- gam(Binary.det ~ Sex + s(TL) + s(Month, bs = 'cc', k = 4) + s(Month, bs = 'cc', by = Sex, k = 4) + s(TA) +  
              s(SST, bs = 'cc', k = 4) + s(Tag, bs = 're'), select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod2)
Mod3 <- gam(Binary.det ~ Sex + s(TL) + s(Month, bs = 'cc', k = 4) + s(Month, bs = 'cc', by = Sex, k = 4) + s(TA) + 
              s(SST, bs = 'cc', k = 4) + s(Tag, bs = 're', by = Sex), select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod3)

# Hypothesis of sexual segregation: does the Month:Sex interaction improve the fit?
Mod4 <- gam(Binary.det ~ Sex + s(TL) + s(Month, bs = 'cc', k = 4) + s(TA) + s(SST, bs = 'cc', k = 4) + s(Tag, bs = 're'), 
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod4)

AIC(Mod1, Mod2, Mod3, Mod4)
# Tags as random effects improve fit, but not the interaction Tag by Sex
# Interactions improve the fit
# Mod2 is the best model
# The term s(TL) has no effect

# Simplest best model
final_model <- gam(Binary.det ~ Sex + s(Month, bs = 'cc', k = 4) + s(Month, bs = 'cc', by = Sex, k = 4) + s(TA) + 
                     s(SST, bs = 'cc', k = 4) + s(Tag, bs = 're'), select = T, family = binomial, data = model_input, method = 'REML') 
summary(final_model)
AIC(final_model)

# Detection probability
model_pred <- data.frame(Month = rep(seq(1, 12, length.out = 500), 2),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = mean(model_input$SST),
                         TA = mean(model_input$TA),
                         Tag = '5970')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = Month, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + ylim(c(0, 1)) +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  scale_x_continuous(name = 'Month', breaks = seq(1, 12, 1), 
                     labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure 5.pdf', width = 15, height = 9.6, units = 'cm')

# Figure S1a - surface temperature
model_pred <- data.frame(Month = mean(model_input$Month),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = rep(seq(min(model_input$SST), max(model_input$SST), length.out = 500), 2),
                         TA = mean(model_input$TA),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = SST, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') + 
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + xlab('Surface temperature (°C)') +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure S2a.pdf', width = 15, height = 9.6, units = 'cm')

# Figure S1b - tide amplitude
model_pred <- data.frame(Month = mean(model_input$Month),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = mean(model_input$SST),
                         TA = rep(seq(min(model_input$TA), max(model_input$TA), length.out = 500), 2),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = TA, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + xlab('Tide amplitude (m)') +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure S2b.pdf', width = 15, height = 9.6, units = 'cm')


#--------------------------------- Figures 6 and S3 ----------------------------

# Roaming index - Roaming patterns

model_input <- read.csv('Daily_data.csv')
model_input <- model_input %>% mutate(Sex = factor(Sex), Tag = factor(Tag), Season = factor(Season), Date = as.Date(Date)) %>% filter(Binary.det > 0)
model_input <- model_input[complete.cases(model_input), ]

# Do the tags as random effects improve the fit?
Mod1 <- gam(cbind(Roam, fail.Roam) ~ Sex + Season + s(TL) + s(Prop_females, by = Sex) + s(Prop_males, by = Sex) + s(TA) + s(SST), 
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod1)
Mod2 <- gam(cbind(Roam, fail.Roam) ~ Sex + Season + s(TL) + s(Prop_females, by = Sex) + s(Prop_males, by = Sex) + s(TA) + s(SST) + s(Tag, bs = 're'), 
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod2)
Mod3 <- gam(cbind(Roam, fail.Roam) ~ Sex + Season + s(TL) + s(Prop_females, by = Sex) + s(Prop_males, by = Sex) + s(TA) + s(SST) + s(Tag, bs = 're', by = Sex), 
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod3)

# Hypothesis of sexual segregation: does the Prop_females:Sex and Prop_males:Sex interactions improve the fit? 
Mod4 <- gam(cbind(Roam, fail.Roam) ~ Sex + Season + s(TL) + s(Prop_females) + s(Prop_males) + s(TA) + s(SST) + s(Tag, bs = 're'), 
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod4)

AIC(Mod1, Mod2, Mod3, Mod4)
# Tags as random effects improve fit, but not the interaction Tag by Sex
# Interactions improve the fit
# Mod2 is the best model
# The terms s(TA) and Season have no effect

# Simplest best model
final_model <- gam(cbind(Roam, fail.Roam) ~ Sex + s(Prop_females, k = -1) + s(Prop_females, by = Sex, k = -1) + 
                     s(Prop_males, k = -1) + s(Prop_males, by = Sex, k = -1) + s(SST) + s(Tag, bs = 're'), 
                   select = T, family = binomial, data = model_input, method = 'REML')
summary(final_model)

# Roaming patterns - prop females
model_pred <- data.frame(Prop_females = rep(seq(0, 1, length.out = 500), 2),
                         Prop_males = mean(model_input$Prop_males),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = mean(model_input$SST),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = Prop_females, y = fit, group = Sex, fill = Sex)) + ylab('Roaming index') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + ylim(c(0, 1)) +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  scale_x_continuous(name = 'Proportion of females') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3)) 
ggsave('Figure 6a.pdf', width = 15, height = 9.6, units = 'cm')

# Roaming patterns - prop males
model_pred <- data.frame(Prop_males = rep(seq(0, 1, length.out = 500), 2),
                         Prop_females = mean(model_input$Prop_females),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = mean(model_input$SST),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = Prop_males, y = fit, group = Sex, fill = Sex)) + ylab('Roaming index') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + ylim(c(0, 1)) +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  scale_x_continuous(name = 'Proportion of males') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3)) 
ggsave('Figure 6b.pdf', width = 15, height = 9.6, units = 'cm')

# Figure S3 - surface temperature
model_pred <- data.frame(Prop_females = mean(model_input$Prop_females),
                         Prop_males = mean(model_input$Prop_males),
                         Sex = c(rep('female', 500), rep('male', 500)), 
                         SST = rep(seq(min(model_input$SST), max(model_input$SST), length.out = 500), 2),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = SST, y = fit, group = Sex, fill = Sex)) + ylab('Roaming index') + 
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + xlab('Surface temperature (°C)') +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3)) 
ggsave('Figure S3.pdf', width = 15, height = 9.6, units = 'cm')


#--------------------------------- Figure 7 ----------------------------------

# Space use by sex

# Refining the Shortest Paths (RSP) and Actel packages --> read guide in https://github.com/YuriNiella/RSP
setwd('~/Actel')

# Sex differences in space use considering abundance season
dat_raw <- read.csv('Detections.csv') # Available at Dryad repository
dat <- subset(dat_raw, Timestamp <= as.Date('2021-02-28'))
dat <- subset(dat, Timestamp >= as.Date('2019-11-01'))
dat <- subset(dat, Signal != '5973')

# Create water raster of study area (NA for land and open sea, 1 for study area)
Ext <- extent(-63.69, -63.595, -42.42, -42.26)
efective_coast <- crop(coast, Ext)
writeOGR(efective_coast, '.', 'Study_area', driver = 'ESRI Shapefile')
coast_raster <- loadShape(shape = 'Study_area.shp', size = 0.0002, buffer = 0.07, type = 'land', coord.x = 'Longitude', coord.y = 'Latitude')
df_sea <- matrix(c(-63.6705, -63.6073, -63.6023, -63.6102, -63.6087, -63.4477, -63.5232, -63.6705,
                   -42.2373, -42.3178, -42.3583, -42.4095, -42.4204, -42.4098, -42.2203, -42.2373), 8, 2) # polygon to get rid of the open sea parts
df_sea <- SpatialPolygons(list(Polygons(list(Polygon(df_sea)), ID = 1)), proj4string = crs(coast))

# Shortest in-water paths
trans <- transitionLayer(coast_raster, directions = 16)
run1 <- explore(tz = 'America/Buenos_Aires', report = F)
rsp_run <- runRSP(input = run1, t.layer = trans, coord.x = 'Longitude', coord.y = 'Latitude', distance = 250,
                  time.step = 10, min.time = 10, max.time = 24)

# By abundance seasons

# Space use areas: Dynamic Brownian Bridge Movement Models
dbbmm_run_1 <- dynBBMM(input = rsp_run, base.raster = coast_raster, UTM = 20,
                       start.time = '2019-11-01 00:00:00', stop.time = '2020-02-29 23:59:00')
dbbmm_run_2 <- dynBBMM(input = rsp_run, base.raster = coast_raster, UTM = 20,
                       start.time = '2020-03-01 00:00:00', stop.time = '2020-09-30 23:59:00')
dbbmm_run_3 <- dynBBMM(input = rsp_run, base.raster = coast_raster, UTM = 20,
                       start.time = '2020-10-01 00:00:00', stop.time = '2021-02-28 23:59:00')

# In-water areas in RSP tracks
areas_run_1 <- getAreas(dbbmm_run_1, type = 'group', breaks = 0.5)
areas_run_2 <- getAreas(dbbmm_run_2, type = 'group', breaks = 0.5)
areas_run_3 <- getAreas(dbbmm_run_3, type = 'group', breaks = 0.5)

# Overlap between sex
overlap_run_1 <- getOverlaps(areas_run_1)
overlap_run_2 <- getOverlaps(areas_run_2)
overlap_run_3 <- getOverlaps(areas_run_3)

source('plotOverlaps_adw.R')
p1 <- plotOverlaps_adw(overlaps = overlap_run_1, areas = areas_run_1, base.raster = coast_raster, group = c('hembra', 'macho'), 
                       level = 0.5, land.col = 'white', title = 'High season 19/20')
p2 <- plotOverlaps_adw(overlaps = overlap_run_2, areas = areas_run_2, base.raster = coast_raster, group = c('hembra', 'macho'), 
                       level = 0.5, land.col = 'white', title = 'Low season 20')
p3 <- plotOverlaps_adw(overlaps = overlap_run_3, areas = areas_run_3, base.raster = coast_raster, group = c('hembra', 'macho'), 
                       level = 0.5, land.col = 'white', title = 'High season 20/21')

ggsave('Figure 7.pdf', width = 15, height = 7.5, units = 'cm', plot = grid.arrange(p1, p2, p3, ncol = 3))

# Percentage of female areas use overlapped with male area use, and viceversa

#1
pixel_male <- table(values(areas_run_1$rasters$macho$`0.5`) > 0)[[2]]
pixel_female <- table(values(areas_run_1$rasters$hembra$`0.5`) > 0)[[2]]
porc_female <- overlap_run_1$areas$`0.5`$percentage[2] * 100
porc_male <- pixel_female * overlap_run_1$areas$`0.5`$percentage[2] / pixel_male * 100
porc_male
porc_female

#2
pixel_male <- table(values(areas_run_2$rasters$macho$`0.5`) > 0)[[2]]
pixel_female <- table(values(areas_run_2$rasters$hembra$`0.5`) > 0)[[2]]
porc_female <- overlap_run_2$areas$`0.5`$percentage[2] * 100
porc_male <- pixel_female * overlap_run_2$areas$`0.5`$percentage[2] / pixel_male * 100
porc_male
porc_female

#3
pixel_male <- table(values(areas_run_3$rasters$macho$`0.5`) > 0)[[2]]
pixel_female <- table(values(areas_run_3$rasters$hembra$`0.5`) > 0)[[2]]
porc_male <- overlap_run_3$areas$`0.5`$percentage[2] * 100
porc_female <- pixel_male * overlap_run_3$areas$`0.5`$percentage[2] / pixel_female * 100
porc_male
porc_female


#--------------------------------- Figure 8 ----------------------------------

dat <- dat_raw %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != '5973') # # Tag 5973 was detected in 5 or less days in total so it is not considered

# Is there variation in the hourly detection patterns among individuals?

# Transform data to a detection matrix
dat$Tag <- as.character(dat$Tag)

detec_list <- list()
for(i in sort(unique(dat$Tag))) {
  sub <- subset(dat, Tag == i)
  detec_list[[i]] <- as.numeric(sub$Hour)
}

# Rao's test of homogeneity: are differences among hourly patterns of sharks?
rao.test(circular(detec_list[[1]], units = 'hours'), circular(detec_list[[2]], units = 'hours'),
         circular(detec_list[[3]], units = 'hours'), circular(detec_list[[4]], units = 'hours'),
         circular(detec_list[[5]], units = 'hours'), circular(detec_list[[6]], units = 'hours'),
         circular(detec_list[[7]], units = 'hours'), circular(detec_list[[8]], units = 'hours'),
         circular(detec_list[[9]], units = 'hours'), circular(detec_list[[10]], units = 'hours'),
         circular(detec_list[[11]], units = 'hours'), circular(detec_list[[12]], units = 'hours'),
         circular(detec_list[[13]], units = 'hours'), circular(detec_list[[14]], units = 'hours'),
         circular(detec_list[[15]], units = 'hours'), circular(detec_list[[16]], units = 'hours'),
         circular(detec_list[[17]], units = 'hours'), circular(detec_list[[18]], units = 'hours'),
         alpha = 0.05)

# Rao's test of uniformity: does any shark show a uniform hourly pattern?
for(i in sort(unique(dat$Tag))) {
  each_shark <- rao.spacing.test(circular(detec_list[[i]], units = 'hours'), alpha = 0.1)
  print(each_shark)
  print(i)
}

# Heat map of proportion of detections per hour for each shark
hourly_detec <- table(dat$Hour, dat$Tag)
hourly_detec <- hourly_detec[rev(rownames(hourly_detec)), ]
hourly_prop <- t(round(prop.table(hourly_detec, margin = 2), 4))
dd <- dist(hourly_prop, method = 'euclidean', diag = F, upper = F, p = 2)
clusters <- hclust(dd, method = 'average', members = NULL)
plot(clusters)

pdf('Figure 8a.pdf', width = 5.5, height = 4)
pheatmap(hourly_prop, cluster_rows = F, cluster_cols = T, border_color = NA, color = viridis(124), 
         breaks = seq(0, 0.124, 0.001), clustering_method = 'ward.D', angle_col = 45, cutree_cols = 2, treeheight_col = 0)
dev.off()

# Plots proportions
nocturnals <- 5970
no_pattern_female <- c(5962, 5969, 5965, 5963, 5968, 5979, 5980, 5975, 5981, 5978)
no_pattern_male <- c(5964, 5967, 5971, 5972, 5974, 5976, 5977)
prop_df = as.data.frame(hourly_prop)
names(prop_df) <- c('Hour', 'Tag', 'Proportion')
prop_df$Type <- ifelse(prop_df$Tag %in% nocturnals, 'Odd pattern', ifelse(prop_df$Tag %in% no_pattern_female, 'Female', 'Male'))
prop_df$Proportion.sd <- prop_df$Proportion
prop_df <- prop_df %>% group_by(Type, Hour) %>% summarize(Proportion = mean(Proportion), Proportion.sd = sd(Proportion.sd)) %>%
  mutate(Hour = as.numeric(Hour))

ggplot(prop_df, aes(x = rev(Hour), y = Proportion, group = Type, col = Type, fill = Type)) +
  geom_errorbar(aes(ymin = Proportion - Proportion.sd, ymax = Proportion + Proportion.sd), position = position_dodge(0.5), width = 0) +
  geom_point() +
  geom_line() +
  scale_fill_viridis_d(option = 'D') + scale_colour_viridis_d(option = 'D') + 
  scale_x_continuous(name = 'Time of day', breaks = c(1, 7, 13, 19, 25), labels = c('00:00', '06:00', '12:00', '18:00', '23:59')) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3)) 
ggsave('Figure 8b.pdf', width = 15, height = 9.6, units = 'cm')


#--------------------------------- Data preparation for figure 9 ------------------------------

# Intra-day tidal effect on detection probability

dat <- dat_hour %>% filter(Date <= as.Date('2021-02-28')) %>% # avoid underrepresented month at the end of study period
  filter(Tag != 5973) %>%# # Tag 5973 was detected in 5 or less days in total so it is not considered
  mutate(Season = cut(Date, breaks = as.Date(c('2019-11-01', '2020-03-01', '2020-10-01', '2021-03-12')),
                      labels = c('high season', 'low season', 'high season'),
                      include.lowest = T, right = F))

# Prepare data
dat.x <- dat %>% group_by(Date, Tag) %>% summarize(Season = unique(Season), Year = unique(Year), Month = unique(Month), Sex = unique(Sex), TL = unique(TL))
hours <- rep(0:23, each = dim(dat.x)[1])
dat.x <- dat.x[rep(seq_len(nrow(dat.x)), 24), ]
dat.x$Hour <- hours
dat <- dat %>% right_join(dat.x, by = c('Date', 'Tag', 'Hour'))
dat <- dat[, c('Date', 'Year.y', 'Month.y', 'Hour', 'Tag', 'Raw.det', 'Roam', 'Sex.y', 'TL.y', 'Season.y')]
names(dat) <- c(names(dat_hour), 'Season')
dat[is.na(dat)] <- 0

# Create tide predictor for study period
# Tidal data from nearest station (validation period 08-2015 to 08-2016)
tides_old <- read.table('Predictors/MADRYN_prediccion_2015-2016.txt') # astronomic prediction for the validation period available at the Dryad repository
names(tides_old) <- c('Date', 'Hour', 'Tide')
tides_old$Date <- as.Date(tides_old$Date, format = '%d/%m/%Y')
tides_old$Hour <- substr(tides_old$Hour, 1, 2)
tides_old <- tides_old %>% mutate(Hour = as.integer(Hour)) %>% filter(Date <= as.Date('2016-03-28'))
tides_old$Long_date <- paste(tides_old$Date, tides_old$Hour)
tides_old$Long_date <- as.POSIXct(tides_old$Long_date, format = '%Y-%m-%d %H')

# Tidal data from the array as validation data
tides_valid <- read.csv('Predictors/Data_loggers.csv') # in situ tidal data from data loggers available at the Dryad repository
names(tides_valid) <- c('Date', 'Hour', 'Tide')
tides_valid$Date <- as.Date(tides_valid$Date, format = '%d/%m/%Y')
tides_valid$Hour <- as.POSIXct(tides_valid$Hour, format = '%H:%M:%S')
tides_valid$Hour <- substr(tides_valid$Hour, 12, 13)
tides_valid <- tides_valid %>% mutate(Hour = as.integer(Hour))
tides_valid <- tides_valid %>% group_by(Date, Hour) %>% summarize(Tide = mean(Tide))
tides_valid$Long_date <- paste(tides_valid$Date, tides_valid$Hour)
tides_valid$Long_date <- as.POSIXct(tides_valid$Long_date, format = '%Y-%m-%d %H')

# Find best tide correction
first_corrA <- tides_valid %>% filter(Date %in% c(as.Date('2015-08-15'), as.Date('2015-08-16')))
first_corrB <- tides_old %>% filter(Date %in% c(as.Date('2015-08-15'), as.Date('2015-08-16')))
plot(first_corrB$Long_date + 3.5*60*60, first_corrB$Tide)
points(first_corrA$Long_date, first_corrA$Tide - 12, pch = 20, col = 'blue')

second_corrA <- tides_valid %>% filter(Date %in% c(as.Date('2015-11-14'), as.Date('2015-11-15')))
second_corrB <- tides_old %>% filter(Date %in% c(as.Date('2015-11-14'), as.Date('2015-11-15')))
plot(second_corrB$Long_date + 3.5*60*60, second_corrB$Tide)
points(second_corrA$Long_date, second_corrA$Tide - 12, pch = 20, col = 'blue')

third_corrA <- tides_valid %>% filter(Date %in% c(as.Date('2016-03-12'), as.Date('2016-03-13')))
third_corrB <- tides_old %>% filter(Date %in% c(as.Date('2016-03-12'), as.Date('2016-03-13')))
plot(third_corrB$Long_date + 3.5*60*60, third_corrB$Tide)
points(third_corrA$Long_date, third_corrA$Tide - 12, pch = 20, col = 'blue')

# Roughly +3.5 hours to the tide prediction of the nearest port (Puerto Madryn) matches tidal data in all three validation cases

# Tidal data from nearest station (study period)
tides <- read.table('Predictors/MADRYN_prediccion_2019-2021.txt') # astronomic prediction for the study period available at the Dryad repository
names(tides) <- c('Date', 'Hour', 'Tide')
tides$Date <- as.Date(tides$Date, format = '%d/%m/%Y')
tides$Hour <- substr(tides$Hour, 1, 2)
tides <- tides %>% mutate(Hour = as.integer(Hour))
tides$Long_date <- paste(tides$Date, tides$Hour)
tides$Long_date <- as.POSIXct(tides$Long_date, format = '%Y-%m-%d %H')
tides$Long_date <- tides$Long_date + 3.5*60*60 # correct data by +3.5 hours
tides$Date <- date(tides$Long_date)
tides$Hour <- hour(tides$Long_date)
tides <- tides %>% filter(Date < as.Date('2021-04-01'))
tides <- tides[which(!is.na(tides$Date)), ]

# To understand the global effect of tide create three categorical variables 
# Tide height: is it a high, intermediate or low tidal moment of the day? 
tides$Tide_height <- NA
q25 <- quantile(tides$Tide, probs = 0.25)
q75 <- quantile(tides$Tide, probs = 0.75)
tides <- tides %>%
  mutate(Tide_height = case_when(Tide < q25 ~ 'low',
                              Tide > q75 ~ 'high',
                              Tide >= q25 & Tide <= q75 ~ 'intermediate'))

# Tide direction: is it an inflow or outflow tidal moment of the day?
tides$Tide_direction <- NA
for(i in 1:dim(tides)[1]) {
  if(tides$Tide[i] > tides$Tide[i + 1]) {
    tides$Tide_direction[i + 1] <- 'outflow'
  }
  if(tides$Tide[i] < tides$Tide[i + 1]) {
    tides$Tide_direction[i + 1] <- 'inflow'
  }
  if(tides$Tide[i] == tides$Tide[i + 1]) {
    tides$Tide_direction[i] <- 'inflow'
    tides$Tide_direction[i + 1] <- 'inflow'
  }
}
tides$Tide_direction[1] <- 'outflow'

# Tide strength: is it a weak or strong tidal moment of the day? 
tides$Tide_strength <- 'strong'
for(i in 1:dim(tides)[1]) {
    if(tides$Tide_direction[i] != tides$Tide_direction[i + 1]) {
      tide_dif <- tides$Tide[i + 1] - tides$Tide[i]
      if(tide_dif > 0 & tide_dif < 0.5) {
        tides$Tide_strength[i] <- 'weak'
        tides$Tide_strength[i + 1] <- 'weak'
      }
      if(tide_dif < 0 & tide_dif > -0.5) {
        tides$Tide_strength[i] <- 'weak'
        tides$Tide_strength[i + 1] <- 'weak'
      }
      tide_dif <- tides$Tide[i - 1] - tides$Tide[i]
      if(tide_dif > 0 & tide_dif < 0.5) {
        tides$Tide_strength[i] <- 'weak'
        tides$Tide_strength[i - 1] <- 'weak'
      }
      if(tide_dif < 0 & tide_dif > -0.5) {
        tides$Tide_strength[i] <- 'weak'
        tides$Tide_strength[i - 1] <- 'weak'
      }
    }
}  

# Sunlight position
daytime <- getSunlightTimes(as.Date(unique(tides$Date)), lat = -42.398, lon = -63.613, tz = "America/Buenos_Aires",
                        keep = c('sunriseEnd', 'sunsetStart', 'dawn', 'dusk', 'night', 'nightEnd')) # sunlight times
daytime <- daytime[, -2]
daytime <- daytime[, -2]

tides <- tides %>% left_join(daytime, by = c('Date' = 'date'))

tides$Daytime <- NA
for(i in 1:dim(tides)[1]) {
  if(tides$Long_date[i] > tides$nightEnd[i] & tides$Long_date[i] < tides$sunriseEnd[i] + 3600) {
    tides$Daytime[i] <- "dawn"
  }
  if(tides$Long_date[i] > tides$sunriseEnd[i] + 3600 & tides$Long_date[i] < tides$sunsetStart[i] - 3600) {
    tides$Daytime[i] <- "daylight"
  }
  if(tides$Long_date[i] > tides$sunsetStart[i] - 3600 & tides$Long_date[i] < tides$night[i]) {
    tides$Daytime[i] <- "dusk"
  }
  if(tides$Long_date[i] < tides$nightEnd[i]) {
    tides$Daytime[i] <- "night"
  }
  if(tides$Long_date[i] > tides$night[i]) {
    tides$Daytime[i] <- "night"
  }
}
tides <- tides[, c(1:7, 14)]

# Paste tide predictor to model input
model_input <- dat %>% left_join(tides[, -4], by = c('Date' = 'Date', 'Hour' = 'Hour'))
model_input$Binary.det <- ifelse(model_input$Raw.det > 0, 1, 0)
model_input$fail.Roam <- 5 - model_input$Roam # number of non-visited receivers 

source('HighstatLib.R') # for VIF multicollinearity analysis

# Prdictor correlacion
vars <- c('Hour', 'Tide')
corvif(model_input[, vars])
pairs(model_input[, vars], lower.panel = panel.cor, cex.labels = 1, labels = vars)

write.csv(model_input, 'Intraday_data.csv', row.names = F)


#--------------------------------- Figures 9 and S4 ------------------------------

# Residency index - detection probability

model_input <- read.csv('Intraday_data.csv')
model_input <- model_input %>% mutate(Sex = factor(Sex), Tag = factor(Tag), Date = as.Date(Date), Season = factor(Season),
                                      Tide_height = factor(Tide_height), Tide_direction = factor(Tide_direction),
                                      Tide_strength = factor(Tide_strength))

# Simpler models
Mod1 <- gam(Binary.det ~ Sex + Season + s(TL) + s(Hour, bs = 'cc') + Tide_height + Tide_direction + Tide_strength,
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod1)
Mod2 <- gam(Binary.det ~ Sex + s(TL) + s(Hour, bs = 'cc') + Tide_height + Tide_direction,
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod2)
# no variables artificially inflate deviance, correlation expected to be low
# vars Tide_strength and Season have no strong effect on detection probability

# Does Sex interactions improve the fit?
Mod3 <- gam(Binary.det ~ Sex + s(TL) + s(Hour, bs = 'cc') + s(Hour, bs = 'cc', by = Sex) + Tide_height + Tide_direction,
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod3)
Mod4 <- gam(Binary.det ~ Sex + s(TL) + s(Hour, bs = 'cc') + s(Hour, bs = 'cc', by = Sex) + Tide_height * Sex + Tide_direction * Sex,
            select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod4)
  # All interactions improve the fit

# Random effects
Mod5 <- gam(Binary.det ~ Sex + s(TL) + s(Hour, bs = 'cc') + s(Hour, bs = 'cc', by = Sex) + Tide_height * Sex + Tide_direction * Sex + 
              s(Tag, bs = 're'), select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod5)
Mod6 <- gam(Binary.det ~ Sex + s(TL) + s(Hour, bs = 'cc') + s(Hour, bs = 'cc', by = Sex) + Tide_height * Sex + Tide_direction * Sex +  
              s(Tag, bs = 're', by = Sex), select = T, family = binomial, data = model_input, method = 'REML')
summary(Mod6)
# Tags as random effects improve fit, but not the interaction Tag by Sex
# The term s(TL) has no effect when Tag is present
# Mod5 is the best model

# Simplest best model
final_model <- gam(Binary.det ~ Sex + s(Hour, bs = 'cc') + s(Hour, bs = 'cc', by = Sex) + 
                     Tide_height * Sex + Tide_direction * Sex + s(Tag, bs = 're'), 
                   select = T, family = binomial, data = model_input, method = 'REML')
summary(final_model)
AIC(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, final_model)

# Detection probability ~ Tide_height
model_pred <- data.frame(Sex = c(rep('female', 1500), rep('male', 1500)), 
                         Hour = mean(model_input$Hour),
                         Tide_height = rep(c(rep('low', 500), rep('intermediate', 500), rep('high', 500)), 2),
                         Tide_direction = 'outflow',
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      
model_pred <- transform(model_pred, Tide_height = factor(Tide_height, levels = c('low', 'intermediate', 'high'))) # ordering of factors

ggplot(data = model_pred, aes(x = Tide_height, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') + xlab('Tide height') +
  geom_pointrange(aes(ymin = ICinf, ymax = ICsup, color = Sex),
    position = position_dodge(0.3)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure 9a.pdf', width = 15, height = 9.6, units = 'cm')

# Detection probability ~ Tide_direction
model_pred <- data.frame(Sex = c(rep('female', 1000), rep('male', 1000)), 
                         Hour = mean(model_input$Hour),
                         Tide_height = 'intermediate',
                         Tide_direction = rep(c(rep('inflow', 500), rep('outflow', 500)), 2),
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = Tide_direction, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') + xlab('Tide direction') +
  geom_pointrange(aes(ymin = ICinf, ymax = ICsup, color = Sex),
                  position = position_dodge(0.3)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure 9b.pdf', width = 15, height = 9.6, units = 'cm')

# Figure S4 ~ Hour of day
model_pred <- data.frame(Sex = c(rep('female', 1500), rep('male', 1500)), 
                         Tide = mean(model_input$Tide),
                         Hour = rep(seq(min(model_input$Hour), max(model_input$Hour), length.out = 500), 2),
                         Tide_height = 'intermediate',
                         Tide_direction = 'inflow',
                         Tag = '5962')
prediction <- predict.gam(final_model, newdata = model_pred, type = 'response', se.fit = T, exclude = 's(Tag)')
model_pred$fit <- as.vector(prediction$fit)
link.pred <- predict.gam(final_model, newdata = model_pred, type = 'link', se.fit = T, exclude = 's(Tag)')     
model_pred$ICinf <- as.vector(1 / (1 + exp(-((link.pred$fit - 1.96 * link.pred$se.fit)))))       
model_pred$ICsup <- as.vector(1 / (1 + exp(-((link.pred$fit + 1.96 * link.pred$se.fit)))))      

ggplot(data = model_pred, aes(x = Hour, y = fit, group = Sex, fill = Sex)) + ylab('Dectection probability') + 
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, alpha = Sex)) + xlab('Time of day') +
  geom_line(aes(color = Sex)) + scale_alpha_discrete(range = c(0.15, 0.5)) + ylim(c(0, 1)) +
  scale_fill_viridis_d(option = 'D') + scale_color_viridis_d(option = 'D') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = 'grey90', size = 0.3),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.3))
ggsave('Figure S4.pdf', width = 15, height = 9.6, units = 'cm')


#--------------------------------- Figure S1 ------------------------------

# Coastline (spatial polygons) - taken from GSHHG coastline database of 15 June 2017
# Dowloaded from https://www.soest.hawaii.edu/pwessel/gshhg/
coastline0 <- 'DOWNLOAD, READ AND SET PATH' # Available at Dryad repository
coast0 <- readOGR(dsn = coastline, layer = 'GSHHS_f_L1_SouthAmerica')

# Recaptures in the region to date
recaptures <- data.frame(lon = c(-63.61, -67.21, -64.91), 
                         lat = c(-42.39, -45.57, -42.62))

# Plot
ggplot() +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_point(data = recaptures, aes(x = lon, y = lat), size  =  0.0001, color = '#FFEA46FF', fill = '#FFEA46FF') +
  scale_y_continuous(name = NULL, breaks = c(-47, -45, -43, -41), labels = c('47º', '45º', '43º', '41º')) + 
  scale_x_continuous(name = NULL, breaks = c(-67, -65, -63), labels = c('67º', '65º', '63º')) +
  ggsn::scalebar(x.min = -62, x.max = -67.8, y.min = -47.2, y.max = -40.5, transform = T,
                 dist = 50, st.size = 2.5, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -62.4, y = -46.94)) +
  coord_equal(xlim = c(-67.8, -62), ylim = c(-47.2, -40.5), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure S1.pdf', width = 9.6, height = 15, units = 'cm')


#-------------------------------------------- END ----------------------------------
