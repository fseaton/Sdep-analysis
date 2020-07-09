# Code for plotting out data for exploratory purposes
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(gganimate)
library(mice)
library(janitor)
library(patchwork)

# atmospheric data plots ####
# Sdep maps
str(Sdep_avg)
Sdep_avg_long <- melt(Sdep_avg, id.vars = c("x","y"))
Sdep_avg_long$year <- substring(Sdep_avg_long$variable, 9,12)

ggplot(Sdep_avg_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Sdep grid average kgS / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Sdep grid average maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

# animated version
Sdep_avg_long %>%
  mutate(year = as.integer(year)) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Sdep (kgS/ha)") +
  coord_fixed() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title = 'Year: {frame_time}', x = '', y = '') +
  transition_time(year) +
  ease_aes('linear')

str(Sdep_for)
Sdep_for_long <- melt(Sdep_for, id.vars = c("x","y"))
Sdep_for_long$year <- substring(Sdep_for_long$variable, 8,11)

ggplot(Sdep_for_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Sdep forest kgS / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Sdep forest maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

str(Sdep_moo)
Sdep_moo_long <- melt(Sdep_moo, id.vars = c("x","y"))
Sdep_moo_long$year <- substring(Sdep_moo_long$variable, 6,9)

ggplot(Sdep_moo_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Sdep moorland kgS / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Sdep moorland maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

# Ndep maps
str(Ndep_avg)
Ndep_avg_long <- melt(Ndep_avg, id.vars = c("x","y"))
Ndep_avg_long$year <- substring(Ndep_avg_long$variable, 9,12)

ggplot(Ndep_avg_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Ndep grid average kgN / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Ndep grid average maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

str(Ndep_for)
Ndep_for_long <- melt(Ndep_for, id.vars = c("x","y"))
Ndep_for_long$year <- substring(Ndep_for_long$variable, 8,11)

ggplot(Ndep_for_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Ndep forest kgN / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Ndep forest maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

str(Ndep_moo)
Ndep_moo_long <- melt(Ndep_moo, id.vars = c("x","y"))
Ndep_moo_long$year <- substring(Ndep_moo_long$variable, 6,9)

ggplot(Ndep_moo_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "") +
  coord_fixed() +
  facet_wrap(~year, nrow = 5) +
  labs(x = "",y = "", title = "Ndep moorland kgN / ha") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("Ndep moorland maps.png", path = "Outputs/Graphs/",
       width = 28, height = 28, units = "cm")

# Plot Sdep and Ndep against each other
colnames(Sdep_avg_long) <- c("x","y","variable","Sdep","year")
colnames(Ndep_avg_long) <- c("x","y","variable","Ndep","year") 
AtDep_avg <- merge(Sdep_avg_long, Ndep_avg_long, by = c("x","y","variable","year"))
str(AtDep_avg)
AtDep_avg$year <- as.integer(AtDep_avg$year)
AtDep_avg <- na.omit(AtDep_avg)
AtDep_avg$square <- paste(AtDep_avg$x, AtDep_avg$y, sep = "_")


ggplot(AtDep_avg, aes(x = Sdep, y = Ndep)) + 
  geom_point(aes(group = square)) + 
  labs(title = 'Year: {frame_time}', x = 'Sdep (kgS/ha)', y = 'Ndep (kgN/ha)') +
  transition_time(year) +
  ease_aes('linear')

AtDep_avg %>%
  pivot_longer(Sdep:Ndep, names_to = "Element",
               values_to = "measure") %>%
  mutate(Element = ifelse(Element == "Sdep", "Sdep (kgS/ha)", "Ndep (kgN/ha)")) %>%
  ggplot(aes(x = x, y = y, fill = measure)) +
  geom_tile() +
  scale_fill_viridis_c(name = "", na.value = "white") +
  coord_fixed() +
  facet_wrap(~Element) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  labs(title = 'Year: {frame_time}', x = '', y = '') +
  transition_time(year) +
  ease_aes('linear')

# forest
colnames(Sdep_for_long) <- c("x","y","variable","Sdep","year")
colnames(Ndep_for_long) <- c("x","y","variable","Ndep","year") 
AtDep_for <- merge(Sdep_for_long, Ndep_for_long, by = c("x","y","variable","year"))
str(AtDep_for)
AtDep_for$year <- as.integer(AtDep_for$year)
AtDep_for <- na.omit(AtDep_for)
AtDep_for$square <- paste(AtDep_for$x, AtDep_for$y, sep = "_")

ggplot(AtDep_for, aes(x = Sdep, y = Ndep)) + 
  geom_point(aes(group = square)) + 
  labs(title = 'Year: {frame_time}', x = 'Sdep (kgS/ha)', y = 'Ndep (kgN/ha)') +
  transition_time(year) +
  ease_aes('linear')

AtDep_for %>%
  pivot_longer(Sdep:Ndep, names_to = "Element",
               values_to = "measure") %>%
  mutate(Element = ifelse(Element == "Sdep", "Sdep (kgS/ha)", "Ndep (kgN/ha)")) %>%
  ggplot(aes(x = x, y = y, fill = measure)) +
  geom_tile() +
  scale_fill_viridis_c(name = "", na.value = "white") +
  coord_fixed() +
  facet_wrap(~Element) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  labs(title = 'Year: {frame_time}', x = '', y = '') +
  transition_time(year) +
  ease_aes('linear')


# moorland
colnames(Sdep_moo_long) <- c("x","y","variable","Sdep","year")
colnames(Ndep_moo_long) <- c("x","y","variable","Ndep","year") 
AtDep_moo <- merge(Sdep_moo_long, Ndep_moo_long, by = c("x","y","variable","year"))
str(AtDep_moo)
AtDep_moo$year <- as.integer(AtDep_moo$year)
AtDep_moo <- na.omit(AtDep_moo)
AtDep_moo$square <- paste(AtDep_moo$x, AtDep_moo$y, sep = "_")

ggplot(AtDep_moo, aes(x = Sdep, y = Ndep)) + 
  geom_point(alpha = 0.2, aes(group = square)) + 
  labs(title = 'Year: {frame_time}', x = 'Sdep (kgS/ha)', y = 'Ndep (kgN/ha)') +
  transition_time(year) +
  ease_aes('linear')

AtDep_moo %>%
  pivot_longer(Sdep:Ndep, names_to = "Element",
               values_to = "measure") %>%
  mutate(Element = ifelse(Element == "Sdep", "Sdep (kgS/ha)", "Ndep (kgN/ha)")) %>%
  ggplot(aes(x = x, y = y, fill = measure)) +
  geom_tile() +
  scale_fill_viridis_c(name = "", na.value = "white") +
  coord_fixed() +
  facet_wrap(~Element) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  labs(title = 'Year: {frame_time}', x = '', y = '') +
  transition_time(year) +
  ease_aes('linear')

# Soil pH data ####
# data manipulation
str(CS78_PH)
str(CS98_PH)
str(CS07_PH)
str(CS16_PH)
str(UK19_PH)

CS78_PH$REP_ID <- paste(CS78_PH$SQUARE_NUM,CS78_PH$REP_NUM, sep = "X")
CS98_PH$REP_ID <- paste(CS98_PH$SQUARE_NUM,CS98_PH$REP_NUM, sep = "X")
CS07_PH$REP_ID <- paste(CS07_PH$SQUARE_NUM,CS07_PH$REP_NUM, sep = "X")
CS16_PH$REP_ID <- paste(CS16_PH$SQUARE_NUM,CS16_PH$REP_NUM, sep = "X")
UK19_PH$REP_ID <- paste(UK19_PH$SQUARE_NUM,UK19_PH$REP_NUM, sep = "X")

PH <- full_join(select(CS78_PH, REP_ID, PH1978),
                select(CS98_PH, REP_ID, PH2000 = PHF2000)) %>%
  full_join(select(CS07_PH, REP_ID, PH2007 = PH2007_IN_WATER)) %>%
  full_join(select(CS16_PH, REP_ID, PH2016 = PH_DIW)) %>%
  full_join(select(UK19_PH, REP_ID, PH2019 = PH_DIW))
str(PH)
summary(PH)

mice::md.pattern(PH)

# histograms
PH_long <- pivot_longer(PH, starts_with("PH"),
                        values_to = "pH",
                        values_drop_na = TRUE)
ggplot(PH_long, aes(x = pH)) + 
  geom_histogram() +
  facet_wrap(~name, scales = "free_y")

PH_long$year <- as.integer(substring(PH_long$name, 3,6))
ggplot(PH_long, aes(x = year, y = pH, group = REP_ID)) +
  geom_line(alpha = 0.2, col = "dodgerblue2")+
  geom_jitter(alpha = 0.2, width = 1, height = 0, shape = 16) 

# calculate differences between survey years 
PH <- PH %>%
  mutate(diff7898 = PH2000 - PH1978,
         diff7807 = PH2007 - PH1978,
         diff7816 = PH2016 - PH1978,
         diff7819 = PH2019 - PH1978,
         diff9807 = PH2007 - PH2000,
         diff9816 = PH2016 - PH2000,
         diff9819 = PH2019 - PH2000,
         diff0716 = PH2016 - PH2007,
         diff0719 = PH2019 - PH2007) %>%
  mutate(diff0718 = ifelse(!is.na(diff0719), diff0719,
                           ifelse(!is.na(diff0716), diff0716, NA)),
         diff7818 = ifelse(!is.na(diff7819), diff7819,
                           ifelse(!is.na(diff7816), diff7816, NA)),
         diff9818 = ifelse(!is.na(diff9819), diff9819,
                           ifelse(!is.na(diff9816), diff9816, NA)))
summary(PH)

PH_diff_long <- PH %>%
  select(REP_ID, starts_with("diff")) %>%
  pivot_longer(starts_with("diff"),
               values_to = "pH",
               values_drop_na = TRUE) %>%
  mutate(name = as.factor(name)) %>%
  mutate(name = forcats::fct_inorder(name))
ggplot(PH_diff_long, aes(x = pH)) + 
  geom_histogram() +
  facet_wrap(~name, scales = "free_y") +
  geom_vline(xintercept = 0)

# select only most recent change and convert into wide format for plotting
PH_Diff_wide <- select(PH, REP_ID, diff0718) %>%
  na.omit() %>%
  left_join(select(CS07_PLOTS, REP_ID, POINT_X, POINT_Y))
summary(PH_Diff_wide)

ggplot(PH_Diff_wide, aes(x = POINT_X, y = POINT_Y, colour = diff0718)) +
  geom_jitter(width = 5000, height = 5000) +
  coord_fixed() +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_Diff_wide$diff0718)),
                         name = "pH change", na.value = "white") +
  theme_dark()

# ** pH maps ####
library(sf)
library(leaflet)

# convert to sf object
CS_PH_loc <- PH_Diff_wide %>%
  select(POINT_X, POINT_Y) %>%
  as.matrix() %>%
  st_multipoint(dim="XY") %>%
  st_sfc(crs = 27700) %>% 
  st_transform(crs = 4326) %>%
  st_cast("POINT") 

CS_PH_loc <- st_sf(cbind(select(PH_Diff_wide, REP_ID, pH_change = diff0718),CS_PH_loc))


# Create variable for colouring of points. first cut the continuous variable
# into bins - these bins are now factors
CS_PH_loc$pH_lev <- cut(CS_PH_loc$pH_change, 
                        c(-3,-1.5,-1,-0.5,0,0.5,1,1.5,3))

pHCol <- colorFactor(palette = 'RdBu', CS_PH_loc$pH_lev)

# add random jitter to points so not overplotting 
CS_PH_loc_jitter <- st_jitter(CS_PH_loc, factor = 0.005)

# read in UK boundary shapefile
UK_boundary <- st_read("../../../GBR_adm/GBR_adm0.shp")

# plot interactively
leaflet() %>%
  addPolygons(data = UK_boundary, stroke = FALSE, 
              color = "black") %>%
  addCircleMarkers(data = CS_PH_loc_jitter, radius = 5,
                   label = CS_PH_loc$REP_ID, 
                   color = ~pHCol(CS_PH_loc$pH_lev),
                   fillOpacity = 1, stroke = FALSE) %>%
  addLegend('topright', pal = pHCol, values = CS_PH_loc$pH_lev,
            title = 'pH change',
            opacity = 1)

# plot histograms of difference between survey years wrapping together 16 and 19
PH_diff_long %>% filter(name %in%
                          c("diff7807","diff9807","diff7898",
                            "diff7818","diff9818","diff0718")) %>%
  ggplot(aes(x = pH)) + 
  geom_histogram() +
  facet_wrap(~name) +
  geom_vline(xintercept = 0)
ggsave("pH change histograms facet by survey comparison.png",
       path = "Outputs/Graphs/",
       width = 20, height = 12, units = "cm")

# remove 18 variables for consistency later in script
PH_diff_long <- filter(PH_diff_long,
                       names %in% c("diff7898",
                                    "diff7807",
                                    "diff7816",
                                    "diff7819",
                                    "diff9807",
                                    "diff9816",
                                    "diff9819",
                                    "diff0716",
                                    "diff0719"))

# ** breakdown by AVC data ####
# AVC data manipulation
hab07 <- select(CS07_IBD, REP_ID = REP_ID07, AVC07) %>%
  unique()
hab98 <- select(CS98_IBD, REP_ID = REP_ID98, AVC98) %>%
  unique()
hab78 <- select(CS78_IBD, REP_ID = REP_ID78, AVC78) %>%
  unique()

# create combined AVC variable, if 07 has AVC use that otherwise use 98 then 78.
# There are only 3 sites with no AVC data and I can't see how to get theirs as
# they don't appear in 2016/19.
hab <- full_join(hab07, hab98) %>% full_join(hab78) %>%
  mutate_if(is.factor, as.character) %>%
  mutate(AVC = ifelse(!is.na(AVC07), AVC07,
                     ifelse(!is.na(AVC98), AVC98,
                            ifelse(!is.na(AVC78), AVC78, NA)))) %>%
  mutate(AVC_desc = recode(AVC,
                           `1` = "Crops/Weeds",
                           `2` = "Tall herb/grass",
                           `3` = "Fertile grassland",
                           `4` = "Infertile grassland",
                           `5` = "Lowland wooded",
                           `6` = "Upland wooded",
                           `7` = "Moorland grass/mosaic",
                           `8` = "Heath/bog"))

# calculate total change in pH over survey years
PH$change_dir <- rowSums(select(PH, diff7898, diff9807, diff0716, diff0719), na.rm = TRUE)
summary(PH$change_dir)
filter(PH,change_dir == 0) %>% select(starts_with("diff")) %>%
  summary()

PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7807), PH$diff7807, PH$change_dir) 
PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7816), PH$diff7816, PH$change_dir) 

# Combine pH and AVC data and convert to long format
PH_long_hab <- left_join(PH, select(hab, REP_ID, AVC = AVC_desc)) %>%
  droplevels() %>%
  pivot_longer(starts_with("PH"),
               values_to = "pH",
               values_drop_na = TRUE) %>%
  mutate(year = as.integer(substring(name, 3,6))) %>%
  mutate(year = ifelse(year == 2000, 1998, year))

# plots of pH change over time
PH_long_hab %>% 
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = year, y = pH, group = REP_ID)) +
  geom_line(alpha = 0.5, aes(colour = change_dir) )+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.8,
              aes(colour = change_dir)) +
  facet_wrap(~AVC, nrow = 2) + 
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  theme_dark()
ggsave("pH change over time facetted by AVC.png", path = "Outputs/Graphs/",
       width = 28, height = 12, units = "cm")

PH_long_hab %>% 
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = year, y = pH)) +
  geom_line(alpha = 0.5, aes(colour = change_dir, group = REP_ID))+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.5,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA) +
  facet_wrap(~AVC, nrow = 2) + 
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("pH change over time boxplots facetted by AVC.png", path = "Outputs/Graphs/",
       width = 28, height = 15, units = "cm")

# combine ph difference and AVC data
PH_diff_long <- left_join(PH_diff_long, 
                          select(hab, REP_ID, AVC = AVC_desc)) %>%
  droplevels()
table(PH_diff_long$AVC)
get_dupes(PH_diff_long, REP_ID, name)

PH_diff_long %>%
  filter(!is.na(AVC)) %>% 
  filter(name %in% c("diff7807","diff7898","diff9807","diff7818","diff9818","diff0718")) %>%
  ggplot(aes(x = pH)) + 
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(AVC ~ name, scales = "free_y")
ggsave("pH difference histograms facetted by AVC and year.png",
       path = "Outputs/Graphs/", width = 28, height = 24, units = "cm")

# ** Soil pH in CaCl2 ####
# Only have pH in CaCl2 data for 2007 onwards
str(CS78_PH)
str(CS98_PH)
str(CS07_PH)
str(CS16_PH)
str(UK19_PH)

# data manipulation
PHC <- full_join(select(CS07_PH, REP_ID, PHC2007 = PH2007_IN_CACL2),
                select(CS16_PH, REP_ID, PHC2016 = PH_CACL2)) %>%
  full_join(select(UK19_PH, REP_ID, PHC2019 = PH_CACL)) %>%
  mutate(pH_change = ifelse(!is.na(PHC2019), PHC2019 - PHC2007,
                            ifelse(!is.na(PHC2016), PHC2016 - PHC2007, NA))) %>%
  left_join(unique(select(hab, REP_ID, AVC = AVC_desc)))
str(PHC)
summary(PHC)
md.pattern(PHC)

# histograms of pH CaCl2 change
PHC %>% filter(!is.na(AVC)) %>%
  ggplot(aes(x = pH_change)) + geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC, nrow = 2)
ggsave("pH CaCl2 change 07 to 1619 facet by AVC.png", path = "Outputs/Graphs/",
       width = 28, height = 12, units = "cm")

p1 <-PHC %>% 
  ggplot(aes(x = pH_change)) + geom_histogram() +
  geom_vline(xintercept = 0)+
  labs(x = "pH change", title = bquote("pH (CaCl"[2]*")")) +
  scale_x_continuous(limits = c(-3,3))+
  scale_y_continuous(limits = c(0,110), expand = c(0,0))
p2 <- PH_diff_long %>% filter(name %in% c("diff0716","diff0719")) %>%
  ggplot(aes(x = pH)) + geom_histogram() + 
  geom_vline(xintercept = 0) +
  labs(x = "", title = "pH (DIW)")+
  scale_x_continuous(limits = c(-3,3))+
  scale_y_continuous(limits = c(0,110), expand = c(0,0))
p2/p1
ggsave("pH change 07 to 1619 DIW and CaCl2.png", path = "Outputs/Graphs/",
       width = 15, height = 18, units = "cm")

# data manipulation into long format
PHC_long <- PHC %>%
  pivot_longer(starts_with("PHC"), names_to = "year",
               names_prefix = "PHC", values_to = "pH_CaCl2",
               values_drop_na = TRUE)
str(PHC_long)

# boxplot/line/scatter plot of pH CaCl2 change over time
PHC_long %>% mutate(year = as.numeric(year)) %>%
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = year, y = pH_CaCl2)) +
  geom_jitter(shape = 16, size = 0.5, alpha = 0.5,
              width = 1, height = 0) +
  geom_boxplot(aes(group = year), fill = NA) +
  geom_line(aes(group = REP_ID, colour = pH_change), alpha = 0.5) +
  facet_wrap(~AVC, nrow = 2) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*abs(max(PHC_long$pH_change)))
ggsave("pH CaCl2 over time boxplots facet by AVC.png",
       path = "Outputs/Graphs/",
       width =28, height = 15, units = "cm")


# combine pH in CaCl2 and DIW and plot against each other
phc_wide_diff <- PH %>%
  mutate(pH_diw_change = ifelse(!is.na(diff0719),diff0719,
                                ifelse(!is.na(diff0716), diff0716, NA))) %>%
  select(REP_ID, PH2007, PH2016, PH2019, pH_diw_change) %>%
  full_join(PHC)

ggplot(phc_wide_diff, aes(x = pH_diw_change, y = pH_change)) + 
  geom_abline(intercept = 0,slope = 1, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_point() +
  # geom_smooth(method = "lm") +
  labs(x = "pH (DIW) change", y = bquote("pH (CaCl"[2]*") change"))
ggsave("pH change over time DIW vs CaCl2 scatterplot.png", 
       path = "Outputs/Graphs/",
       width = 15, height = 15, units = "cm")

# there is one sample with a NA for AVC so removing
phc_wide_diff %>% filter(!is.na(AVC)) %>%
  ggplot(aes(x = pH_diw_change, y = pH_change)) + 
  geom_abline(intercept = 0,slope = 1, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  facet_wrap(~AVC, nrow = 2) +
  geom_point() +
  # geom_smooth(method = "lm") +
  labs(x = "pH (DIW) change", y = bquote("pH (CaCl"[2]*") change"))
ggsave("pH change over time DIW vs CaCl2 scatterplots facet by AVC.png", 
       path = "Outputs/Graphs/",
       width = 28, height = 15, units = "cm")

# Plant Ellenberg scores ####
# Data manipulation
str(CS07_IBD)
str(CS98_IBD)
str(CS78_IBD)
str(GM16_IBD)

# get GMEP data to have CS REP_ID 
GMEP_CS_match <- CS16_PH %>%
  select(SQUARE_NUM, GMEP_NUM, PLOT_TYPE, REP_NUM) %>%
  filter(!is.na(GMEP_NUM)) %>%
  mutate(CS_REP_ID = paste0(SQUARE_NUM, PLOT_TYPE, REP_NUM),
         GMEP_REP_ID = paste0(GMEP_NUM, PLOT_TYPE, REP_NUM)) %>%
  select(CS_REP_ID, GMEP_REP_ID)

# Get GMEP data into similar format to CS data
GM16_IBD <- GM16_IBD %>%
  right_join(GMEP_CS_match, by = c("REP_ID" = "GMEP_REP_ID")) %>%
  select(REP_ID16 = CS_REP_ID, R16 = PH, N16 = FERT, L16 = LIGHT, F16 = WET)

# Combine IBD files for the different years
IBD_comb <- full_join(CS07_IBD, CS98_IBD, by = c("REP_ID07" = "REP_ID98")) %>%
  full_join(CS78_IBD, by = c("REP_ID07" = "REP_ID78")) %>%
  full_join(GM16_IBD, by = c("REP_ID07" = "REP_ID16"))

# Use AVC data from 2007 if it is there, otherwise 98 or 78
IBD_comb$AVC <- ifelse(!is.na(IBD_comb$AVC07), IBD_comb$AVC07,
                           ifelse(!is.na(IBD_comb$AVC98), IBD_comb$AVC98,
                                  IBD_comb$AVC78))
summary(IBD_comb$AVC)
# get plot type from REP_ID
IBD_comb$PLOT_TYPE <- gsub("[^a-zA-Z]", "", IBD_comb$REP_ID07)
summary(as.factor(IBD_comb$PLOT_TYPE))

# Calculate difference in Ell R over the years
ELL <- IBD_comb %>%
  mutate(diff7898 = R98 - R78,
         diff7807 = R07 - R78,
         diff9807 = R07 - R98,
         diff7816 = R16 - R78,
         diff9816 = R16 - R98,
         diff0716 = R16 - R07,
         AVC_desc = recode(AVC,
                           `1` = "Crops/Weeds",
                           `2` = "Tall herb/grass",
                           `3` = "Fertile grassland",
                           `4` = "Infertile grassland",
                           `5` = "Lowland wooded",
                           `6` = "Upland wooded",
                           `7` = "Moorland grass/mosaic",
                           `8` = "Heath/bog")) 
# Calculate overall Ell R change
ELL$Rchange <- rowSums(select(ELL, diff7898, diff9807, diff0716), na.rm = TRUE)
summary(ELL$Rchange)
filter(ELL, Rchange == 0) %>% select(starts_with("diff")) %>%
  summary()

ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff7807), ELL$diff7807, ELL$Rchange) 
ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff9807), ELL$diff9807, ELL$Rchange) 
ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff7898), ELL$diff7898, ELL$Rchange) 
summary(ELL$Rchange)

# Convert Ell R change into long format
ELL_diff_long <- ELL %>%
  select(REP_ID = REP_ID07, AVC = AVC_desc, starts_with("diff")) %>%
  droplevels %>%
  pivot_longer(starts_with("diff"), values_to = "Ell_R") %>%
  filter(!is.na(Ell_R))

ELL_diff_long %>%
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = Ell_R)) + 
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(name ~ AVC, scales = "free_y")
ggsave("Ellenberg R change histograms all plots facetted AVC year.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units ="cm")

# Convert Ellenberg R scores file to long format
ELL_R_LONG <- ELL %>%  select(REP_ID = REP_ID07,
                              PLOT_TYPE = PLOT_TYPE.x,
                              AVC = AVC_desc, R07, R98, R78, R16, Rchange) %>%
  filter(PLOT_TYPE == "X") %>%
  droplevels() %>%
  select(-PLOT_TYPE) %>%
  pivot_longer(cols = c(R07,R98,R78,R16), names_to = "year",
               names_prefix = "R") %>%
  mutate(year = ifelse(year == "07", 2007,
                       ifelse(year == "98", 1998,
                              ifelse(year == "78", 1978, 
                                     ifelse(year == "16", 2016, NA)))))
str(ELL_R_LONG)  
summary(ELL_R_LONG$AVC)

# Ellenberg R score change over time boxplot/scatter/line graph
ELL_R_LONG %>%
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = year, y = value)) +
  geom_line(alpha = 0.5, aes(group = REP_ID, colour = Rchange))+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.5,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA, width = 3) +
  facet_wrap(~AVC, nrow = 2) + 
  labs(y = "Ellenberg R") +
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "Ell R change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("Ellenberg R change over time X plots boxplots facetted by AVC.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units = "cm")



# Combined pH and ellenberg R ####
# change graphs - combine difference stats
ph_ell_comb <- ELL_diff_long %>%
  select(-AVC) %>%
  filter(grepl("X", REP_ID)) %>%
  full_join(filter(PH_diff_long, name %in% unique(ELL_diff_long$name)))
str(ph_ell_comb)
md.pattern(ph_ell_comb)

# scatter plot of pH change against Ellenberg R change facetted by survey year
# comparison and AVC
ph_ell_comb %>%
  filter(!is.na(AVC)) %>%
  ggplot(aes(x = pH, y = Ell_R)) +
  geom_point() +
  facet_grid(name~AVC) +
  geom_smooth(method = "lm")
ggsave("Ellenberg R vs pH change by AVC and survey.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units = "cm")


# pH in CaCl2 and Ellenberg R 
phc_ell_comb <- ELL_diff_long %>%
  filter(grepl("X", REP_ID)) %>%
  filter(name %in% c("diff0716","diff0719")) %>%
  full_join(na.omit(select(PHC, REP_ID, pH_change)))
str(phc_ell_comb)
md.pattern(phc_ell_comb)

# scatter plot of pH CaCl2 vs Ellenberg R change, no facets
phc_ell_comb %>%
  # filter(!is.na(AVC)) %>%
  ggplot(aes(x = pH_change, y = Ell_R)) +
  geom_point() +
  # facet_wrap(~AVC) +
  geom_smooth(method = "lm")

# combine pH and Ellenberg R scores in wide format
str(PH)
str(ELL)

# differences
ph_ell_wide_diff <- ELL %>%
  filter(REP_ID07 %in% PH$REP_ID) %>%
  select(REP_ID = REP_ID07,  diff7807, diff7898, diff9807) %>%
  full_join(select(PH, REP_ID,diff7807, diff7898, diff9807),
                         suffix = c("_ellR","_ph"), by = "REP_ID")
md.pattern(ph_ell_wide_diff)

psych::pairs.panels(select_if(ph_ell_wide_diff, is.numeric))

# actual pH and Ellenberg R values
ph_ell_wide <- ELL %>%
  filter(REP_ID07 %in% PH$REP_ID) %>%
  select(REP_ID = REP_ID07,  R78, R98, R07, R16) %>%
  full_join(select(PH, REP_ID, PH1978, PH2000, PH2007, PH2016, PH2019),
            suffix = c("_ellR","_phw"), by = "REP_ID")
md.pattern(ph_ell_wide)

psych::pairs.panels(select_if(ph_ell_wide, is.numeric), 
                    ellipses = FALSE, rug = FALSE,
                    method = "spearman")

# strongest correlation tends to be current year
ph_ell_long <- ph_ell_wide %>%
  pivot_longer(ends_with(c("78","98","00","07","16","19")),
               values_drop_na = TRUE) %>%
  mutate(year = sapply(strsplit(name, "[A-Z]{1,2}"),"[",2),
         variable = sapply(strsplit(name, "[0-9]{1,4}"),"[",1)) %>%
  mutate(year = as.numeric(recode(year, "98" = "1998",
                                  "78" = "1978",
                                  "07" = "2007",
                                  "16" = "2018",
                                  "19" = "2018",
                                  "2000" = "1998",
                                  "2016" = "2018",
                                  "2019" = "2018")),
         variable = recode(variable, 
                           "PH" = "Soil_pH",
                           "R" = "Ell_R")) %>%
  select(-name) %>%
  unique() %>%
  pivot_wider(names_from = variable,
              values_from = value)

ggplot(ph_ell_long, aes(x = Soil_pH, y = Ell_R)) + 
  geom_point(size = 1) +
  facet_wrap(~year) +
  labs(x = "Soil pH", y = "Ellenberg R")
ggsave("Ellenberg R vs Soil pH by year.png", path = "Outputs/Graphs",
       width = 15, height = 15, units = "cm")

# checking if sigmoidal curve seems appropriate
x <- seq(3.5,9,0.1)
c1 <- 4.5
c2 <- 1.5
c3 <- 5
c4 <- 2
y <- c1/(1 + exp(-c2*(x - c3))) + c4
c1 <- 4.5
c2 <- 2
c3 <- 4.5
c4 <- 2
y2 <- c1/(1 + exp(-c2*(x - c3))) + c4
dat <- data.frame(x,y,y2)

# asymmetrical sigmoidal curve
c1 <- 4
c2 <- 1.5
c3 <- 3.5
c4 <- 2.5
c5 <- 6
y3 <- c1/((1 + exp(-c2*(x - c3)))^c5) + c4
dat <- data.frame(x,y,y2,y3)


ggplot(ph_ell_long, aes(x = Soil_pH, y = Ell_R)) +
  geom_point(size = 1) +
  facet_wrap(~year) + 
  geom_smooth() +
  geom_line(data = dat, aes(x = x, y = y), 
            colour = "purple", size = 1.5) +
  geom_line(data = dat, aes(x = x, y = y2), 
            colour = "red", size = 1.5) +
  geom_line(data = dat, aes(x = x, y = y3), 
            colour = "orange", size = 1.5) +
  labs(x = "Soil pH", y = "Ellenberg R")

# seems reasonable - interestingly in 1978 the Ellenberg R value for any given
# pH is higher. So in model need to make sure that c3 varies by year, and not
# sure about the other parameters.

# include CaCl2 but only recent years for graphical simplicity
phr_ell_wide <- ELL %>%
  filter(REP_ID07 %in% PH$REP_ID) %>%
  select(REP_ID = REP_ID07,  R07, R16) %>%
  full_join(select(PH, REP_ID, PH2007, PH2016),
            suffix = c("_ellR","_phw"), by = "REP_ID") %>%
  full_join(select(PHC, REP_ID, PHC2007, PHC2016))

psych::pairs.panels(select_if(phr_ell_wide, is.numeric), 
                    ellipses = FALSE, rug = FALSE,
                    method = "spearman")


# Soil moisture ####
UK19_WET <- UK19_WET %>%
  mutate(REP_ID = paste0(`Square number...1`,`X plot...2`)) %>%
  select(-starts_with("Square")) %>%
  select(-starts_with("X "))

MOISTURE <- CS_tier4 %>%
  mutate(REP_ID = ifelse(!is.na(REP_ID07), REP_ID07,
                         ifelse(!is.na(REP_ID98), REP_ID98,
                                ifelse(!is.na(REP_ID78), REP_ID78, NA)))) %>%
  select(REP_ID, MOISTURE_CONTENT_07, MOISTURE_CONTENT_98) %>%
  full_join(select(UK19_WET, REP_ID, MOISTURE_CONTENT_19 = `g water/wet weight of soil`)) %>%
  # mutate(VWC_19 = 100*VWC_19) %>%
  pivot_longer(starts_with("MOISTURE"), names_to = "variable", 
               values_to = "Moisture") %>%
  mutate(Year = ifelse(variable == "MOISTURE_CONTENT_07", 2007,
                       ifelse(variable == "MOISTURE_CONTENT_98", 1998,
                              ifelse(variable == "MOISTURE_CONTENT_19", 2019, NA))))

ggplot(MOISTURE, aes(x = Moisture)) + 
  geom_histogram() +
  facet_wrap(~Year, nrow = 3, scales = "free_y") +
  labs(x = "Soil moisture (%)") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
ggsave("Soil moisture g per wet soil histograms.png", path = "Outputs/Graphs/",
       width = 12, height = 15, units = "cm")

MOISTURE_PH <- PH_long %>%
  mutate(year = ifelse(year == 2000, 1998, year)) %>%
  full_join(rename(MOISTURE, year = Year)) %>%
  filter(year %in% c(1998, 2007, 2019))
ggplot(MOISTURE_PH, aes(x = Moisture, y = pH,
                        colour = as.factor(year), fill = as.factor(year))) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(formula = 'y ~ log(x) + log(x)', method = "lm") +
  scale_colour_brewer(palette = "Dark2", name = "Year") +
  scale_fill_brewer(palette = "Dark2", name = "Year") +
  labs(x = "Soil moisture (g/g wet soil)") +
  scale_y_continuous(limits = c(3.3,9.2))
ggsave("Soil moisture vs pH log line.png", path = "Outputs/Graphs",
       width =14, height = 12, units = "cm")  


ggplot(MOISTURE_PH, aes(x = Moisture, y = pH)) +
  geom_point(alpha = 0.5, size = 0.8, colour = "dodgerblue2") +
    geom_smooth(colour = "black", formula = 'y ~ log(x)', method = "lm") +
  facet_wrap(~year) +
  labs(x = "Soil moisture (g/g wet soil)") +
  scale_y_continuous(limits = c(3.3,9.2))
ggsave("Soil moisture vs pH facet log line.png", path = "Outputs/Graphs",
       width =18, height = 12, units = "cm")  


library(brms)
MOISTURE_PH$Year_cat <- as.factor(MOISTURE_PH$year)
get_prior(bf(pH ~ a*exp(b*Moisture) + c*exp(d*Moisture) + e, 
             a + b + c + d + e ~1, nl = TRUE),
          data = MOISTURE_PH)
