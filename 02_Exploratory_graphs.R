# Code for plotting out data for exploratory purposes
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(gganimate)

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
         diff0719 = PH2019 - PH2007)
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


# breakdown by habitat data
# habitat data manipulation
hab07 <- select(CS07_IBD, REP_ID = REP_ID07, BH_CODE = BH07) %>%
  left_join(select(BHPH_NAMES, BH07 = BROAD_HABITAT, BH_CODE)) %>%
  select(-BH_CODE) %>%
  unique()
hab98 <- select(CS98_IBD, REP_ID = REP_ID98, BH_CODE = BH98) %>%
  left_join(select(BHPH_NAMES, BH98 = BROAD_HABITAT, BH_CODE))%>%
  select(-BH_CODE)%>%
  unique()
hab78 <- select(CS78_IBD, REP_ID = REP_ID78, BH_CODE = BH78) %>%
  left_join(select(BHPH_NAMES, BH78 = BROAD_HABITAT, BH_CODE))%>%
  select(-BH_CODE)%>%
  unique()

# create combined broad habitat variable, if 07 has broad habitat use that
# otherwise use 98 then 78. There are only 3 sites with no habitat data and I
# can't see how to get theirs as they don't appear in 2016/19.
hab <- full_join(hab07, hab98) %>% full_join(hab78) %>%
  mutate_if(is.factor, as.character) %>%
  mutate(BH = ifelse(!is.na(BH07), BH07,
                     ifelse(!is.na(BH98), BH98,
                            ifelse(!is.na(BH78), BH78, NA))))

# calculate total change in pH over survey years
PH$change_dir <- rowSums(select(PH, diff7898, diff9807, diff0716, diff0719), na.rm = TRUE)
summary(PH$change_dir)
filter(PH,change_dir == 0) %>% select(starts_with("diff")) %>%
  summary()

PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7807), PH$diff7807, PH$change_dir) 
PH$change_dir <- ifelse(PH$change_dir == 0 & !is.na(PH$diff7816), PH$diff7816, PH$change_dir) 

# Combine pH and habitat data and convert to long format
PH_long_hab <- left_join(PH, select(hab, REP_ID, BH07 = BH)) %>%
  droplevels() %>%
  mutate(BH07 = recode_factor(BH07,
                              "Boundary and Linear Features" = "Other",
                              "Littoral Sediment" = "Coastal",
                              "Sea" = "Coastal",
                              "Supra-littoral Rock" = "Coastal",
                              "Supra-littoral Sediment" = "Coastal",
                              "Urban" = "Other",
                              "Standing Open Waters and Canals" = "Other",
                              "Montane" = "Other",
                              "Inland Rock" = "Other",
                              "Calcareous Grassland" = "Other",
                              "Arable and Horticulture" = "Arable",
                              "Broadleaved Mixed and Yew Woodland" = "Broadleaf",
                              "Coniferous Woodland" = "Conifer",
                              "Dwarf Shrub Heath" = "Heath",
                              "Fen, Marsh, Swamp" = "Marsh",
                              "Acid Grassland" = "Acid Gr",
                              "Improved Grassland" = "Improved Gr",
                              "Neutral Grassland" = "Neutral Gr")) %>%
  mutate(BH07 = replace_na(BH07, "Other")) %>%
  pivot_longer(starts_with("PH"),
                      values_to = "pH",
                      values_drop_na = TRUE) %>%
  mutate(year = as.integer(substring(name, 3,6))) 

# plots of pH change over time
ggplot(PH_long_hab, aes(x = year, y = pH, group = REP_ID)) +
  geom_line(alpha = 0.5, aes(colour = change_dir) )+
  geom_jitter(alpha = 0.2, width = 1, height = 0, shape = 16, colour = change_dir) +
  facet_wrap(~BH07) + 
  scale_colour_distiller(palette = "RdBu", direction = -1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change") +
  theme_dark()


# combine ph difference and habitat data
PH_diff_long <- left_join(PH_diff_long, hab07) %>%
  droplevels()
summary(PH_diff_long$BH07)
PH_diff_long <- PH_diff_long %>%
  mutate(BH07 = recode_factor(BH07,
                              "Boundary and Linear Features" = "Other",
                              "Littoral Sediment" = "Other",
                              "Sea" = "Other",
                              "Supra-littoral Rock" = "Other",
                              "Supra-littoral Sediment" = "Other",
                              "Urban" = "Other",
                              "Standing Open Waters and Canals" = "Other",
                              "Bracken" = "Other",
                              "Inland Rock" = "Other",
                              "Calcareous Grassland" = "Other",
                              "Arable and Horticulture" = "Arable",
                              "Broadleaved Mixed and Yew Woodland" = "Broadleaf",
                              "Coniferous Woodland" = "Conifer",
                              "Dwarf Shrub Heath" = "Heath",
                              "Fen, Marsh, Swamp" = "Marsh",
                              "Acid Grassland" = "Acid Gr",
                              "Improved Grassland" = "Improved Gr",
                              "Neutral Grassland" = "Neutral Gr")) %>%
  mutate(BH07 = replace_na(BH07, "Other"))

ggplot(PH_diff_long, aes(x = pH)) + 
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(BH07~name, scales = "free_y")
