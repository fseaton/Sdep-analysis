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

# Cumulative deposition per square ####
AtDepavg_cumdep_sq <- AtDep_avg %>%
  mutate(year_cat = ifelse(year %in% 1970:1978, "70_78",
                           ifelse(year %in% 1979:1998, "79_98",
                                  ifelse(year %in% 1999:2007, "99_07",
                                         ifelse(year %in% 2008:2018, "08_18",NA))))) %>%
  group_by(x,y,square, year_cat) %>%
  summarise(Sdep = sum(Sdep),
            Ndep = sum(Ndep)) %>%
  ungroup() %>%
  pivot_wider(names_from = year_cat, 
              values_from = c(Sdep, Ndep)) 

psych::pairs.panels(select(AtDepavg_cumdep_sq,-x,-y,-square))

# merge with CS plots
dep_x <- AtDepavg_cumdep_sq$x
dep_y <- AtDepavg_cumdep_sq$y

CS_m <- CS07_PLOTS %>% select(plot_x = POINT_X, 
                              plot_y = POINT_Y, REP_ID)
for(i in 1:nrow(CS_m)) {
  CS_m[i,"x"] <- dep_x[which.min(abs(dep_x - CS_m$plot_x[i]))]
  CS_m[i,"y"] <- dep_y[which.min(abs(dep_y - CS_m$plot_y[i]))]
}

CS_Atdep <- left_join(CS_m, AtDepavg_cumdep_sq)
summary(CS_Atdep)
psych::multi.hist(select_if(CS_Atdep, is.numeric))

AtDepfor_cumdep_sq <- AtDep_for %>%
  mutate(year_cat = ifelse(year %in% 1970:1978, "70_78",
                           ifelse(year %in% 1979:1998, "79_98",
                                  ifelse(year %in% 1999:2007, "99_07",
                                         ifelse(year %in% 2008:2018, "08_18",NA))))) %>%
  group_by(x,y,square, year_cat) %>%
  summarise(Sdep = sum(Sdep),
            Ndep = sum(Ndep)) %>%
  ungroup() %>%
  pivot_wider(names_from = year_cat, 
              values_from = c(Sdep, Ndep)) 

# difference since 1970
AtDepavg_diff_sq <- AtDep_avg %>%
  mutate(year_cat = ifelse(year %in% 1970:1978, "70_78",
                           ifelse(year %in% 1979:1998, "79_98",
                                  ifelse(year %in% 1999:2007, "99_07",
                                         ifelse(year %in% 2008:2018, "08_18",NA))))) %>%
  group_by(x,y,square, year_cat) %>%
  summarise(Sdep = sum(Sdep),
            Ndep = sum(Ndep)) %>%
  ungroup() %>%
  pivot_wider(names_from = year_cat, 
              values_from = c(Sdep, Ndep)) 

psych::pairs.panels(select(AtDepavg_cumdep_sq,-x,-y,-square))

# merge with CS plots
dep_x <- AtDepavg_cumdep_sq$x
dep_y <- AtDepavg_cumdep_sq$y
  
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
PH$change_dir <- rowSums(select(PH, diff7898, diff9807, diff0719), na.rm = TRUE)
summary(PH$change_dir)
filter(PH,change_dir == 0) %>% select(starts_with("diff")) %>%
  summary()

PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7807), PH$diff7807, PH$change_dir) 
PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7819), PH$diff7819, PH$change_dir) 

# Combine pH and AVC data and convert to long format
PH_long_hab <- left_join(PH, select(BH_IMP, REP_ID, Management)) %>%
  droplevels() %>%
  select(-starts_with("diff")) %>%
  pivot_longer(starts_with("PH"),
               names_to = c("Variable","year"),
               names_sep = "_",
               values_to = "pH",
               values_drop_na = TRUE) %>%
  filter(!is.na(Management )) %>%
  mutate(year = as.numeric(year))

# plots of pH change over time
PH_long_hab %>% 
  ggplot(aes(x = year, y = pH, group = REP_ID)) +
  geom_line(alpha = 0.5, aes(colour = change_dir) )+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.8,
              aes(colour = change_dir)) +
  facet_wrap(~Management, nrow = 2) + 
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  theme_dark()
ggsave("pH change over time facetted by Management.png", path = "Outputs/Graphs/",
       width = 12, height = 12, units = "cm")

PH_long_hab %>% 
  ggplot(aes(x = year, y = pH)) +
  geom_line(alpha = 0.5, aes(colour = change_dir, group = REP_ID))+
  geom_jitter(size = 0.2, width = 1, height = 0,
              shape = 16, alpha = 0.1,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA) +
  facet_wrap(~Management, nrow = 2) + 
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("pH change over time boxplots facetted by management.png", path = "Outputs/Graphs/",
       width = 12, height = 15, units = "cm")

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

PH %>%
  select(REP_ID, PHC_2007, PH_2007, PH_2019, PHC_2019) %>%
  pivot_longer(starts_with("PH"), 
               names_to = c("Variable","Year"),
               names_sep = "_") %>%
  na.omit() %>%
  pivot_wider(names_from = "Variable",
              values_from = "value") %>%
  ggplot(aes(x = PH, y = PHC, colour = Year)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  facet_wrap(~Year) +
  labs(x = "pH (DIW)", y = bquote("pH CaCl"[2]))

# Plant Ellenberg scores ####
str(CS19_SP)
table(CS19_SP$PLOT_TYPE)
unique(CS19_SP[CS19_SP$PLOT_TYPE=="XX","REP_ID"])
table(CS19_SP[CS19_SP$PLOT_TYPE=="X","PLOTYEAR"])

str(SPECIES_LIB_TRAITS)
filter(SPECIES_LIB_CODES, COLUMN_NAME == "GROWTH_FORM")

CS18_ELL <- filter(CS19_SP, PLOT_TYPE %in% c("X","XX")) %>%
  mutate(REP_ID = paste0(SQUARE,PLOT_TYPE,PLOT_NUMBER)) %>%
  mutate(REP_ID = gsub("XX","X",REP_ID)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"),
                   GROWTH_FORM)) %>%
  filter(GROWTH_FORM %in% c("f","fe","g","m","s","ss","w")) %>% # filter to vascular plants
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE,
                   .names = "{col}18")) %>%
  rename_with(~gsub("EBERG","",.x))
summary(CS18_ELL)
test <- full_join(CS18_ELL, GM16_IBD, by = c("REP_ID" = "REP_ID16"))
plot(N18 ~ N16, test);abline(0,1)

CS98_ELL <- CS98_SP %>%
  select(REP_ID, BRC_NUMBER, TOTAL_COVER) %>%
  unique() %>%
  filter(TOTAL_COVER > 0) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"),
                   GROWTH_FORM)) %>%
  # filter(GROWTH_FORM %in% c("f","fe","g","m","s","ss","w")) %>% # filter to vascular plants
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(REP_ID) %>%
  summarise(across(starts_with("EBER"), function(x) sum(x, na.rm=TRUE)/length(na.omit(EBERGN)),
                   .names = "{col}98_new")) %>%
  rename_with(~gsub("EBERG","",.x))
test <- full_join(CS98_ELL, CS98_IBD, by = c("REP_ID" = "REP_ID98"))
#par(mfrow=c(2,2))
plot(R98_new ~ R98, test);abline(0,1)
plot(N98_new ~ N98, test);abline(0,1)
plot(W98_new ~ F98, test);abline(0,1)
plot(L98_new ~ L98, test);abline(0,1)
par(mfrow=c(1,1))
summary(CS18_ELL)

X_Ell_comp <- full_join(X_Ell_inner, X_Ell_whole) %>%
  left_join(hab) %>%
  mutate(R_diff = SM_R - WH_R,
         N_diff = SM_N - WH_N,
         W_diff = SM_W - WH_W,
         L_diff = SM_L - WH_L)

p1 <- ggplot(X_Ell_comp, aes(x = WH_R, y = SM_R)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg R 400m"^2~"plot"),
       y = bquote("Ellenberg R 4m"^2~"plot"))

ggplot(X_Ell_comp, aes(x = WH_N, y = SM_N)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc)

ggplot(X_Ell_comp, aes(x = R_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

ggplot(X_Ell_comp, aes(x = WH_R, y = SM_R)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(Year~AVC_desc) +
  theme_bw()

ggplot(X_Ell_comp, aes(x = N_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

X_Ell_comp %>%
  select(Year, REP_ID, ends_with("diff"), AVC_desc) %>%
  pivot_longer(ends_with("diff"), names_to = "Ellenberg") %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  facet_grid(Ellenberg~AVC_desc)

# weighted Ellenberg comparison
X_wEll_comp <- full_join(X_wEll_inner, X_wEll_whole) %>%
  left_join(hab) %>%
  mutate(R_diff = SM_R - WH_R,
         N_diff = SM_N - WH_N,
         W_diff = SM_W - WH_W,
         L_diff = SM_L - WH_L)

p2 <- ggplot(filter(X_wEll_comp, Year != 1978), aes(x = WH_R, y = SM_R)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg R 400m"^2~"plot"),
       y = bquote("Ellenberg R 4m"^2~"plot"))

ggplot(filter(X_wEll_comp, Year != 1978), aes(x = R_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

X_wEll_comp %>%
  select(Year, REP_ID, ends_with("diff"), AVC_desc) %>%
  pivot_longer(ends_with("diff"), names_to = "Ellenberg") %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  facet_grid(Ellenberg~AVC_desc)


t.test(X_wEll_comp$SM_R,X_wEll_comp$WH_R)
# t = -0.47853, df = 18816, p-value = 0.6323
t.test(X_Ell_comp$SM_R,X_Ell_comp$WH_R)
# t = -3.3001, df = 19015, p-value = 0.0009682

x <- na.omit(unique(X_wEll_comp$AVC_desc))
for(i in 1:length(x)){
  
  dat <- filter(X_wEll_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_R,dat$WH_R)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.85153"
# [1] "Crops/Weeds p = 0.04266"
# [1] "Fertile grassland p = 0.92217"
# [1] "Heath/bog p = 2e-05"
# [1] "Moorland grass/mosaic p = 0.12709"
# [1] "Upland wooded p = 0.74556"
# [1] "Infertile grassland p = 0.92811"
# [1] "Lowland wooded p = 0.39651"

x <- na.omit(unique(X_Ell_comp$AVC_desc))
for(i in 1:length(x)){
  dat <- filter(X_Ell_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_R,dat$WH_R)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.91431"
# [1] "Crops/Weeds p = 0.06944"
# [1] "Fertile grassland p = 0.05292"
# [1] "Heath/bog p = 0"
# [1] "Moorland grass/mosaic p = 0"
# [1] "Upland wooded p = 0.00329"
# [1] "Infertile grassland p = 0.06137"
# [1] "Lowland wooded p = 0.96831"

p1 + ggtitle("Unweighted") + p2 + ggtitle("Cover weighted")
ggsave("Ellenberg R plot size comparison.png", path = "Outputs/Graphs/",
       width = 24, height = 12, units = "cm")

p1 <- ggplot(filter(X_Ell_comp, Year != 1978), aes(x = WH_N, y = SM_N)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg N 400m"^2~"plot"),
       y = bquote("Ellenberg N 4m"^2~"plot")) +
  ggtitle("Unweighted")
p2 <- ggplot(filter(X_wEll_comp, Year != 1978), aes(x = WH_N, y = SM_N)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg N 400m"^2~"plot"),
       y = bquote("Ellenberg N 4m"^2~"plot")) +
  ggtitle("Cover weighted")
p1 + p2
ggsave("Ellenberg N plot size comparison.png", path = "Outputs/Graphs/",
       width = 24, height = 12, units = "cm")

t.test(X_wEll_comp$SM_N,X_wEll_comp$WH_N)
# t = -0.12149, df = 18823, p-value = 0.9033
t.test(X_Ell_comp$SM_N,X_Ell_comp$WH_N)
# t = -1.621, df = 19043, p-value = 0.105

# correlations
x <- na.omit(unique(X_wEll_comp$AVC_desc))
for(i in 1:length(x)){
  
  dat <- filter(X_wEll_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_N,dat$WH_N)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.8639"
# [1] "Crops/Weeds p = 0.04501"
# [1] "Fertile grassland p = 0.83626"
# [1] "Heath/bog p = 0.03957"
# [1] "Moorland grass/mosaic p = 0.07602"
# [1] "Upland wooded p = 0.63365"
# [1] "Infertile grassland p = 0.38431"
# [1] "Lowland wooded p = 0.31157"

x <- na.omit(unique(X_Ell_comp$AVC_desc))
for(i in 1:length(x)){
  dat <- filter(X_Ell_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_N,dat$WH_N)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.91795"
# [1] "Crops/Weeds p = 0.05538"
# [1] "Fertile grassland p = 0.38703"
# [1] "Heath/bog p = 0"
# [1] "Moorland grass/mosaic p = 0"
# [1] "Upland wooded p = 0.01399"
# [1] "Infertile grassland p = 0.41287"
# [1] "Lowland wooded p = 0.62069"


# Data manipulation
str(CS07_IBD)
str(CS98_IBD)
str(CS78_IBD)
str(GM16_IBD)
str(CS18_ELL)

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
  full_join(GM16_IBD, by = c("REP_ID07" = "REP_ID16")) %>%
  full_join(CS18_ELL, by = c("REP_ID07" = "REP_ID"))

# Use AVC data from 2007 if it is there, otherwise 98 or 78
IBD_comb$AVC <- ifelse(!is.na(IBD_comb$AVC07), IBD_comb$AVC07,
                           ifelse(!is.na(IBD_comb$AVC98), IBD_comb$AVC98,
                                  IBD_comb$AVC78))
summary(IBD_comb$AVC)
# get plot type from REP_ID
IBD_comb$PLOT_TYPE <- gsub("[^a-zA-Z]", "", IBD_comb$REP_ID07)
summary(as.factor(IBD_comb$PLOT_TYPE))

# Calculate difference in Ell R over the years
ELL <- X_Ell %>%
  select(Year, REP_ID, contains("_R_")) %>%
  pivot_longer(contains("_R_"), names_to = "Ellenberg") %>%
  mutate(Year = as.character(Year)) %>%
  pivot_wider(names_from = Year,
              names_prefix = "R") %>%
  mutate(diff7890 = R1990 - R1978,
         diff9098 = R1998 - R1990,
         diff9807 = R2007 - R1998,
         diff0719 = R2019 - R2007
  ) 
  

# Calculate overall Ell R change
ELL$Rchange <- rowSums(select(ELL, diff7890, diff9098, 
                              diff9807, diff0719), na.rm = TRUE)
summary(ELL$Rchange)
filter(ELL, Rchange == 0) %>% select(starts_with("diff")) %>%
  summary()

# Convert Ell R change into long format
ELL_diff_long <- ELL %>%
  select(REP_ID, starts_with("diff")) %>%
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
ELL %>%
  select(-starts_with("diff")) %>%
  pivot_longer(R1978:R2019,
               names_to = "year",
               names_prefix = "R",
               names_transform = list(year = as.integer)) %>%
  left_join(BH_IMP) %>%
  filter(!is.na(Management)) %>%
  mutate(Management = recode(Management,
                             "High" = "High intensity",
                             "Low" = "Low intensity"),
         Ellenberg = recode(Ellenberg,
                            "SM_R_UW" = "Small unweighted",
                            "SM_R_W" = "Small weighted",
                            "WH_R_UW" = "Full unweighted",
                            "WH_R_W"= "Full weighted")) %>%
  ggplot(aes(x = year, y = value)) +
  geom_line(alpha = 0.5, aes(group = REP_ID, colour = Rchange))+
  geom_jitter(size = 0.2, width = 1, height = 0, 
              shape = 16, alpha = 0.1,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA, width = 3) +
  facet_grid(Ellenberg~Management) + 
  labs(y = "Ellenberg R") +
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(ELL$Rchange)),
                         name = "Ell R change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("Ellenberg R change over time X plots boxplots facetted by Management.png",
       path = "Outputs/Graphs/", width = 15, height = 20, units = "cm")

str(X_Ell)
str(BH_IMP)

plot_dat <- left_join(X_Ell, BH_IMP) %>%
  filter(!is.na(Management)) %>%
  select(REP_ID, Year, contains("_R_"), Management) %>%
  pivot_longer(contains("_R_"), 
               names_to = c("PlotSize","Score","Weighting"),
               names_sep = "_") %>%
  mutate(PlotSize = recode(PlotSize,
                           "SM" = "Small",
                           "WH" = "Full"),
         Weighting = recode(Weighting,
                            "UW" = "Unweighted",
                            "W" = "Weighted"),
         Management = recode(Management,
                             "High" = "High intensity management",
                             "Low" = "Low intensity management")) %>%
  pivot_wider(names_from = PlotSize,
              values_from = value)
ggplot(plot_dat, aes(x = Full, y = Small)) + 
  geom_point(alpha = 0.25, colour = "dodgerblue3") + 
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() + 
  facet_grid(Management~Weighting) +
  labs(x = "Full size plot",  y = "Smaller size plot")
ggsave("Ellenberg R plot size comparison by weighting and management.png",
       path = "Outputs/Graphs/", width = 12, height = 12, units ="cm")

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

# other monotonic curves instead
# monomolecular/Mitscherlich law/von Bertalanffy law.
a <- 7
b <- 0
c <- 0.3
ym1 <- a - (a-b)*exp(-c*x)
plot(x,ym1)

# von bertalanffy
a <- 7
b <- 0.5
c <- -2
yb <- a*(1-exp(-b*(x-c)))
plot(x,yb)

# Michaelis Menten
a <- 9
b <- 4
ym <- a*x/(b + x)


dat <- data.frame(x,ym,y2,y3,yb)


ggplot(ph_ell_long, aes(x = Soil_pH, y = Ell_R)) +
  geom_point(size = 1) +
  facet_wrap(~year) + 
  geom_smooth() +
  geom_line(data = dat, aes(x = x, y = yb), 
            colour = "purple", size = 1.5) +
  geom_line(data = dat, aes(x = x, y = y2), 
            colour = "red", size = 1.5) +
  geom_line(data = dat, aes(x = x, y = y3), 
            colour = "orange", size = 1.5) +
  labs(x = "Soil pH", y = "Ellenberg R")

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
MOISTURE <- CS_tier4 %>%
  mutate(REP_ID = ifelse(!is.na(REP_ID07), REP_ID07,
                         ifelse(!is.na(REP_ID98), REP_ID98,
                                ifelse(!is.na(REP_ID78), REP_ID78, NA)))) %>%
  select(REP_ID, MOISTURE_CONTENT_07, MOISTURE_CONTENT_98) %>%
  full_join(select(UK19_WET, REP_ID, MOISTURE_CONTENT_19 = `g_water/wet_weight_of_soil`)) %>%
  # mutate(VWC_19 = 100*VWC_19) %>%
  pivot_longer(starts_with("MOISTURE"), names_to = "variable", 
               values_to = "Moisture") %>%
  mutate(Year = ifelse(variable == "MOISTURE_CONTENT_07", 2007,
                       ifelse(variable == "MOISTURE_CONTENT_98", 1998,
                              ifelse(variable == "MOISTURE_CONTENT_19", 2019, NA)))) %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y07)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)

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


# CHess soil moisture
soil_moist_long <- soil_moist %>%
  # `colnames<-`(c(paste0(LETTERS[1:14],colnames(soil_moist)))) %>%
  pivot_longer(AJan:LDec, names_to = "Month", values_to = "Moisture") %>%
  mutate(month_num = as.numeric(as.factor(Month)))

ggplot(soil_moist_long, aes(x = Measting, y = Nnorthing, fill = Moisture)) +
  geom_tile() + 
  coord_fixed() + 
  facet_wrap(~month_num) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave("JULES soil moisture maps 2007.png",
       path = "Outputs/Graphs/", width = 30, height = 30, units = "cm")

# Rainfall data ####
ggplot(cs_survey_rainfall, aes(x = mean_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year) +
  labs(x = "Average monthly rainfall for 4 months pre-survey")

p1 <- ggplot(cs_survey_rainfall, aes(x = mean_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Average monthly rainfall for 4 months pre-survey")
p2 <- ggplot(cs_survey_rainfall, aes(x = sum_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Total rainfall for 4 months pre-survey")
p1 + p2

ggplot(cs_survey_rainfall, aes(x = mean_rainfall, y = sum_rainfall)) + 
  geom_point() +
  geom_abline(slope=4, intercept = 0) +
  facet_wrap(~Year)
# basically the same

p1
ggsave("Average monthly rainfall for 4 months pre-survey.png",
       path = "Outputs/Graphs/",width = 12, height = 18, units = "cm")


# Rainfall and soil moisture ####
rain_moist <- cs_survey_rainfall %>% ungroup() %>%
  mutate(Year = as.numeric(Year)) %>%
  full_join(MOISTURE)
summary(rain_moist)
mice::md.pattern(rain_moist)
rain_moist %>% filter(is.na(mean_rainfall)) %>% .$Year %>% table()

ggplot(rain_moist, aes(x = mean_rainfall, y = Moisture, 
                       colour = as.factor(Year))) + 
  geom_point()
ggplot(rain_moist, aes(x = mean_rainfall, y = Moisture)) + 
  geom_point() +
  facet_wrap(~Year)

rain_moist_hab <- left_join(rain_moist, BH_comb)
summary(rain_moist_hab)
ggplot(rain_moist_hab, aes(x = mean_rainfall, y = Moisture)) + 
  geom_point() +
  facet_wrap(~BH_DESC)

rain_moist_man <- left_join(rain_moist, BH_IMP)
summary(rain_moist_man)
ggplot(na.omit(rain_moist_man), 
       aes(x = mean_rainfall, y = Moisture)) + 
  geom_point(aes(colour = Management, fill = Management)) +
  geom_smooth(method="lm", colour = "black") +
  geom_smooth(method="lm", aes(colour = Management, fill = Management)) +
  facet_wrap(~Year)

# compare rainfall to sample and EO soil moisture
rain_moist_hab <- left_join(rain_moist_hab, cs_loc_moist07_long)

ggplot(rain_moist_hab, aes(x = mean_moisture, y = Moisture)) + 
  geom_point(aes(colour = month))

rain_moist_hab07 <- rain_moist_hab %>%
  filter(Year == 2007) %>%
  full_join(cs_loc_moist07_sample_month)
ggplot(rain_moist_hab07, aes(x = eo_moisture, y = Moisture)) + 
  geom_point(aes(colour = month))
ggplot(rain_moist_hab07, aes(x = eo_moisture, y = Moisture)) + 
  geom_point() +
  labs(x = "JULES soil moisture", y = "Field soil moisture")
ggsave("Field moisture compared to JULES soil moisture.png",
       path = "Outputs/Graphs/", width = 20, height = 12, units = "cm")
ggplot(rain_moist_hab07, aes(x = eo_moisture, y = mean_moisture)) + 
  geom_point()
ggplot(rain_moist_hab07, aes(x = mean_rainfall, y = Moisture)) + 
  geom_point(aes(colour = eo_moisture))

rain_moist_hab07 <- left_join(rename(rain_moist_hab07, REPEAT_PLOT_ID = REP_ID), 
                              filter(plot_locations, YEAR == "y07"))
ggplot(rain_moist_hab07, aes(x = eo_moisture, y = Moisture)) + 
  geom_point(aes(colour = N_10_FIG_1M))
p1_eo <- ggplot(filter(rain_moist_hab07, !is.na(eo_moisture)), 
                aes(x = E_10_FIG_1M, y = N_10_FIG_1M)) + 
  geom_point(aes(colour = eo_moisture)) +
  coord_fixed() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  scale_color_continuous(name = "") +
  ggtitle("JULES soil moisture")
p2_meas <- ggplot(filter(rain_moist_hab07, !is.na(Moisture)),
                  aes(x = E_10_FIG_1M, y = N_10_FIG_1M)) + 
  geom_point(aes(colour = Moisture)) +
  coord_fixed() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  scale_color_continuous(name = "") +
  ggtitle("Field soil moisture")
p3_rain <- ggplot(filter(rain_moist_hab07, !is.na(mean_rainfall)),
                  aes(x = E_10_FIG_1M, y = N_10_FIG_1M)) + 
  geom_point(aes(colour = mean_rainfall)) +
  coord_fixed() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  scale_color_continuous(name = "") +
  ggtitle("Mean rainfall")

p1_eo + p2_meas + p3_rain
ggsave("Map of Soil moisture and rainfall over CS 2007.png",
       path = "Outputs/Graphs", width = 30, height = 15, units = "cm")

p1_eo <- ggplot(rain_moist_hab07, aes(x = eo_moisture)) +
  geom_histogram(bins = 40)
p2_meas <- ggplot(rain_moist_hab07, aes(x = Moisture)) +
  geom_histogram(bins = 40)
p3_rain <- ggplot(rain_moist_hab07, aes(x = mean_rainfall)) +
  geom_histogram(bins = 40)

p1_eo + p2_meas + p3_rain

psych::pairs.panels(select(rain_moist_hab07, eo_moisture, Moisture,
                           mean_rainfall), rug=FALSE)

summary(lm(Moisture ~ eo_moisture + mean_rainfall, data = rain_moist_hab07))

summary(lm(Moisture ~ mean_moisture + mean_rainfall, data = rain_moist_hab07))

ggplot(PH_moisture,
       aes(x = mean_moisture, y = PH_2007)) + 
  geom_point()

ggplot(PH_moisture,
       aes(x = Moisture, y = PH_2007)) + 
  geom_point()

summary(lm(Moisture ~ mean_moisture + mean_rainfall, data = rain_moist_hab07))


PH_moisture <- left_join(PH_moisture,
                         select(CS07_CN, REP_ID = REP_ID07,
                                C_PERCENT, N_PERCENT))
summary(lm(PH_2007 ~ log(C_PERCENT) + eo_moisture, 
           data = PH_moisture))
summary(lm(PHC_2007 ~ Moisture + log(C_PERCENT), 
           data = PH_moisture))

ggplot(PH_moisture,
       aes(x = C_PERCENT, y = PH_2007)) + 
  geom_point() +
  scale_x_log10()

# pH and atmospheric deposition ####
ph_atdep <- PH %>% select(REP_ID, diff7807) %>%
  mutate(Year = 2007) %>%
  left_join(CS_plot_atdep)

ggplot(ph_atdep, aes(x = Sdep, y = diff7807))+ 
  geom_point()

ph_atdep <- PH %>% 
  select(REP_ID, diff7819) %>%
  na.omit() %>%
  mutate(Year = 2018) %>%
  left_join(CS_plot_atdep)

ggplot(ph_atdep, aes(x = Sdep, y = diff7819))+ 
  geom_point()

ph_atdep <- PH %>% 
  mutate(diffc0719 = PHC_2019 - PHC_2007) %>%
  select(REP_ID, diff0719, diffc0719) %>%
  filter(!is.na(diff0719)) %>%
  mutate(Year = 2018) %>%
  left_join(CS_plot_atdep)

ggplot(ph_atdep, aes(x = Sdep, y = diff0719))+ 
  geom_point()
ggplot(ph_atdep, aes(x = Sdep, y = diffc0719))+ 
  geom_point()


ph_atdep <- PH %>% select(REP_ID, PH_2007, PHC_2007) %>%
  mutate(Year = 2007) %>%
  left_join(CS_plot_atdep)

ggplot(ph_atdep, aes(x = Sdep, y = PH_2007))+ 
  geom_point()

ggplot(ph_atdep, aes(x = Sdep, y = PHC_2007))+ 
  geom_point()

# soil N against N deposition
summary(CS07_CN)
CN_atdep <- CS07_CN %>% 
  mutate(CN_ratio = C_PERCENT/N_PERCENT,
         Year = 2007,
         REP_ID = REP_ID07) %>%
  select(REP_ID, Year, CN_ratio) %>%
  full_join(mutate(CS98_CN,
                   CN_ratio = C_PERCENT/N_PERCENT,
                   Year = 1998,
                   REP_ID = REP_ID98)) %>%
  left_join(CS_plot_atdep)


ggplot(CN_atdep, 
       aes(x = Ndep, y = CN_ratio)) +
  geom_point() +
  labs(x = "Cumulative N deposition", y = "C:N") +
  facet_wrap(~Year + Habitat, nrow = 3)


# Soil nitrogen ### 

# Investigating which total N based metric is most related to mineralisable
# N/nitrate. Options are N% or N:C, there are multiple metrics of mineralisable
# N but nitrate is known to be related to NPP.
Nitrogen <- CS07_MINN %>%
  mutate(REP_ID = paste0(SQUARE_NUM, PLOT_TYPE, REP_NUM)) %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y07)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID) %>%
  left_join(CS07_CN)

ggplot(Nitrogen, aes(x = N_PERCENT, y = NE_NMINTOT_SOIL)) + 
  geom_point() +
  labs(x = "Total N (%)", y = "Total Mineralisable N (mg N / g dry soil)")
p3 <- ggplot(Nitrogen, aes(x = N_PERCENT, y = NE_NMINTOT_SOM)) + 
  geom_point() +
  labs(x = "Total N (%)", y = "Total Mineralisable N (mg N / g LOI)") +
  scale_y_log10()

cor(Nitrogen$N_PERCENT, Nitrogen$NE_NMINTOT_SOM, method = "spearman",
    use = "complete.obs")
# [1] -0.5087721

ggplot(Nitrogen, aes(x = N_PERCENT, y = NE_NH4N_SOM)) + 
  geom_point() +
  labs(x = "Total N (%)", y = "Total NH4 (mg N / g LOI)") +
  scale_y_log10()
p2<- ggplot(Nitrogen, aes(x = N_PERCENT, y = NE_NO3N_SOM)) + 
  geom_point() +
  labs(x = "Total N (%)", y = "Total NO3 (mg N / g LOI)") +
  scale_y_log10()

cor(Nitrogen$N_PERCENT, Nitrogen$NE_NO3N_SOM, method = "spearman",
    use = "complete.obs")
# [1] -0.5469838

Nitrogen$NC_ratio <- Nitrogen$N_PERCENT/Nitrogen$C_PERCENT

ggplot(Nitrogen, aes(x = NC_ratio, y = NE_NMINTOT_SOIL)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Total Mineralisable N (mg N / g dry soil)") +
  scale_y_log10()

p4 <- ggplot(Nitrogen, aes(x = NC_ratio, y = NE_NMINTOT_SOM)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Total Mineralisable N (mg N / g LOI)") +
  scale_y_log10()

cor(Nitrogen$NC_ratio, Nitrogen$NE_NMINTOT_SOM, method = "spearman",
    use = "complete.obs")
# 0.5645039

ggplot(Nitrogen, aes(x = NC_ratio, y = NE_NH4N_SOM)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Total NH4 (mg N / g LOI)") +
  scale_y_log10()

cor(Nitrogen$NC_ratio, Nitrogen$NE_NH4N_SOM, method = "spearman",
    use = "complete.obs")
# [1] 0.1507093

p1 <- ggplot(Nitrogen, aes(x = NC_ratio, y = NE_NO3N_SOM)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Total NO3 (mg N / g LOI)") +
  scale_y_log10()

cor(Nitrogen$NC_ratio, Nitrogen$NE_NO3N_SOM, method = "spearman",
    use = "complete.obs")
# 0.608635

p2+p1+p3+p4
ggsave("Total N and mineralisable N by LOI.png",
       path = "Outputs/Graphs/", width = 25, height = 20, units = "cm")

Nitrogen <- left_join(Nitrogen, select(CS07_PH, -BATCH_NUM))

ggplot(Nitrogen, aes(x = NC_ratio, y = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total N:C", y = "pH") +
  scale_y_log10()

ggplot(Nitrogen, aes(x = N_PERCENT, y = C_PERCENT)) + 
  geom_point() +
  labs(x = "Total N", y = "Total C") +
  scale_y_log10()


ggplot(Nitrogen, aes(x = NC_ratio, y = NE_NMINTOT_SOM, 
                     colour = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Total Mineralisable N (mg N / g LOI)") +
  scale_y_log10()


coplot(NC_ratio ~ log(NE_NMINTOT_SOM) | PH2007_IN_WATER, Nitrogen)
coplot(log(NE_NMINTOT_SOM) ~ NC_ratio | PH2007_IN_WATER, Nitrogen)


p1 <- ggplot(Nitrogen, aes(x = NC_ratio, y = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total N:C", y = "pH") +
  geom_smooth(method = "lm", fill = "#3366FF")
p2 <- ggplot(Nitrogen, aes(x = NE_NMINTOT_SOM, y = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total Mineralisable N (mg N / g LOI)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", fill = "#3366FF")
p3 <- ggplot(Nitrogen, aes(x = NE_NMINTOT_SOIL, y = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total Mineralisable N (mg N / g soil)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", fill = "#3366FF")
p4 <- ggplot(Nitrogen, aes(x = NE_NO3N_SOM, y = PH2007_IN_WATER)) + 
  geom_point() +
  labs(x = "Total NO3-N (mg N / g LOI)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", fill = "#3366FF")
p1 + p2 + p3 + p4
ggsave("pH and N measurements.png", path = "Outputs/Graphs",
       width = 25, height = 20, units = "cm")

Nitrogen2 <- Nitrogen %>% 
  left_join(filter(BH_IMP, Year == 2007)) %>%
  filter(!is.na(Management))

p1 <- ggplot(Nitrogen2, aes(x = NC_ratio, y = PH2007_IN_WATER, colour = Management)) + 
  geom_point() +
  labs(x = "Total N:C", y = "pH") +
  geom_smooth(method = "lm", aes(fill = Management))
p2 <- ggplot(Nitrogen2, aes(x = NE_NMINTOT_SOM, y = PH2007_IN_WATER, colour = Management)) + 
  geom_point() +
  labs(x = "Total Mineralisable N (mg N / g LOI)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", aes(fill = Management))
p3 <- ggplot(Nitrogen2, aes(x = NE_NMINTOT_SOIL, y = PH2007_IN_WATER, colour = Management)) + 
  geom_point() +
  labs(x = "Total Mineralisable N (mg N / g soil)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", aes(fill = Management))
p4 <- ggplot(Nitrogen2, aes(x = NE_NO3N_SOM, y = PH2007_IN_WATER, colour = Management)) + 
  geom_point() +
  labs(x = "Total NO3-N (mg N / g LOI)", y = "pH") +
  scale_x_log10() +
  geom_smooth(method = "lm", aes(fill = Management))
p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave("pH and N measurements by Management type.png", path = "Outputs/Graphs",
       width = 25, height = 20, units = "cm")

Nitrogen3 <- Nitrogen2 %>%
  left_join(filter(X_Ell, Year == 2007) %>%
              select(REP_ID, contains("_R_")))
p1 <- ggplot(Nitrogen3, aes(x = NC_ratio, y = WH_R_W, colour = Management)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Ellenberg R (full X weighted)") +
  geom_smooth(method = "lm", aes(fill = Management))
p2 <- ggplot(Nitrogen3, aes(x = NC_ratio, y = WH_R_UW, colour = Management)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Ellenberg R (full X unweighted)") +
  geom_smooth(method = "lm", aes(fill = Management))
p3 <- ggplot(Nitrogen3, aes(x = NC_ratio, y = SM_R_W, colour = Management)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Ellenberg R (small X weighted)") +
  geom_smooth(method = "lm", aes(fill = Management))
p4 <- ggplot(Nitrogen3, aes(x = NC_ratio, y = SM_R_UW, colour = Management)) + 
  geom_point() +
  labs(x = "Total N:C", y = "Ellenberg R (small X unweighted)") +
  geom_smooth(method = "lm", aes(fill = Management))
p1+p2+p3+p4+ plot_layout(guides="collect")
ggsave("Ellenberg R and NC by Management type.png", path = "Outputs/Graphs",
       width = 25, height = 20, units = "cm")

# Moisture and rainfall and pH ####

# Look at whether soil moisture actually predicted by change in rainfall to
# evaluate whether change in rainfall is a useful metric
rain_moist <- cs_survey_rainfall %>% ungroup() %>%
  mutate(Year = as.numeric(Year)) %>%
  full_join(MOISTURE)

test <- rain_moist %>% select(-sum_rainfall) %>% 
  pivot_wider(names_from = Year, values_from = c(mean_rainfall, Moisture)) %>% 
  mutate(Moisture_diff9807 = Moisture_2007 - Moisture_1998,
         Moisture_diff0719 = Moisture_2019 - Moisture_2007,
         rain_diff9807 = mean_rainfall_2007 - mean_rainfall_1998,
         rain_diff0719 = mean_rainfall_2019 - mean_rainfall_2007) %>%
  select(REP_ID, contains("diff")) %>%
  pivot_longer(contains("diff"), names_to = c("Variable","Time_period"),
               names_sep = "_diff") %>%
  pivot_wider(names_from = Variable, values_from = value) %>%
  na.omit()

# Also combine with LOI as moisture highly dependent on SOM
test <- LOI %>%
  mutate(LOI_2019 = ifelse(!is.na(LOI_2019), LOI_2019, LOI_2016)) %>%
  mutate(LOI_diff9807 = LOI_2007 - LOI_1998,
         LOI_diff0719 = LOI_2019 - LOI_2007) %>%
  select(REP_ID, contains("diff")) %>%
  pivot_longer(contains("diff"), names_to = c("Variable","Time_period"),
               names_sep = "_diff") %>%
  pivot_wider(names_from = Variable, values_from = value) %>%
  right_join(test)

ggplot(test, aes(x = rain, y = Moisture)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(test, aes(x = LOI, y = Moisture)) +
  geom_point() +
  geom_smooth(method = "lm")

summary(lm(Moisture ~ rain+LOI, test))
# Moisture difference predicted by rain + LOI difference (19%)

# however we want to predict change in pH
test2 <- PH_diff_long %>%
  mutate(Time_period = gsub("diff","",as.character(name))) %>%
  right_join(test)

ggplot(test2, aes(x = rain, y = pH)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(test2, aes(x = LOI, y = pH)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(test2, aes(x = Moisture, y = pH)) +
  geom_point() +
  geom_smooth(method = "lm")

summary(lm(pH ~ Moisture*rain*LOI, test2))
summary(lm(pH ~ Moisture, test2))
summary(lm(pH ~ rain, test2))
summary(lm(pH ~ LOI, test2))
# pretty low variance explained, change in rainfall explained a lot more than
# change in soil moisture or LOI but still low

# Ran the above for different time windows of rainfall calculation and looked at
# correlation. This supported initial assumption that 4 months was best based on
# evidence provided by Don about the upwater monitoring network

# 2 month calculation
cor(test2$pH, test2$rain, use = "complete.obs")
# 0.05777657

# 3 month calculation
cor(test2$pH, test2$rain, use = "complete.obs")
# [1] 0.08958804

# 4 month calculation
cor(test2$pH, test2$rain, use = "complete.obs")
#[1] 0.1373099

# 5 month calculation
cor(test2$pH, test2$rain, use = "complete.obs")
# [1] 0.1050855


# Species functional graphs ####
str(Sp_Ell)

Sp_Ell %>%
  na.omit() %>%
  pivot_longer(BUTTLARVFOOD:Low_grass, names_to = "Function",
               values_to = "Count") %>%
  mutate(Count = ifelse(Count == 1, "Yes","No")) %>%
  group_by(EBERGR, Function) %>%
  count(Count) %>%
  ungroup() %>% group_by(Function, EBERGR) %>%
  mutate(prop = n/sum(n, na.rm = TRUE)) %>%
  filter(Count == "Yes") %>%
  select(-n,-Count) %>%
  pivot_wider(names_from = Function, values_from = prop)

p1 <- Sp_Ell %>%
  unique() %>%
  na.omit() %>%
  pivot_longer(BUTTLARVFOOD:Low_grass, names_to = "Function",
               values_to = "Count") %>%
  mutate(Count = ifelse(Count == 1, "Yes","No"),
         Function = recode(Function, 
                           "BUTTLARVFOOD" = "Butterfly larvae food",
                           "KgSughacovyr" = "Nectar producing",
                           "Low_grass" = "Lowland grass indicators")) %>%
  ggplot(aes(x = EBERGR, fill = Count)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("grey","black")) +
  facet_wrap(~Function, ncol = 1) +
  labs(x = "Ellenberg R", y = "Number of plant species") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  NULL

p2 <- Ell_F %>%
  select(Year, REP_ID, WH_R, starts_with("F_")) %>%
  pivot_longer(starts_with("F_"), names_to = "Function",
               values_to = "value") %>%
  mutate(Function = recode(Function, 
                           "F_Butt" = "Butterfly larvae food",
                           "F_Nectar" = "Nectar producing",
                           "F_Lgrass" = "Lowland grass indicators")) %>%
  ggplot(aes(x = WH_R, y = value))+
  geom_point(alpha = 0.1, colour = "#0072B2") +
  facet_wrap(~Function, ncol = 1) +
  labs(x = "Ellenberg R", y = "Proportion of plant species") +
  theme_minimal() +
  # geom_smooth() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid = element_blank()) +
  NULL

p3 <- Ell_F %>%
  select(Year, REP_ID, WH_R, starts_with("Fr")) %>%
  pivot_longer(starts_with("Fr"), names_to = "Function",
               values_to = "value") %>%
  mutate(Function = recode(Function, 
                           "Fr_Butt" = "Butterfly larvae food",
                           "Fr_Nectar" = "Nectar producing",
                           "Fr_Lgrass" = "Lowland grass indicators")) %>%
  ggplot(aes(x = WH_R, y = value))+
  geom_point(alpha = 0.1, colour = "#0072B2") +
  facet_wrap(~Function, ncol = 1) +
  labs(x = "Ellenberg R", y = "Group richness") +
  theme_minimal() +
  # geom_smooth() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid = element_blank()) +
  NULL

p1 + p2 + p3
ggsave("Ellenberg R and plant functions.png",
       path = "Outputs/Graphs/",
       width = 25, height = 15, units=  "cm")

ggplot(X_Ell_nect, aes(x = WH_R, y = Nectar)) + 
  geom_point(alpha = 0.1) +
  scale_y_log10()

library(brms)

mod_pr <- c(prior(normal(0,1), class = "b"),
            prior(student_t(3, 0, 1), class = "Intercept"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 1), class = "sigma"),
            prior(normal(0,1), class = "ar"))

str(Ell_F)
Ell_F <- Ell_F %>% ungroup() %>%
  mutate(SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         YR = as.factor(Year),
         YRnm = as.integer(YR))

buttmod <- brm(F_Butt ~ WH_R + (YR|YR*SQUARE) + 
                 ar(time = YRnm, gr = REP_ID),
               data = Ell_F, prior = mod_pr, cores = 4,
               iter = 4000)
