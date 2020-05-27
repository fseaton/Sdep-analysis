# Code for plotting out data for exploratory purposes
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(gganimate)
library(mice)

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
  mutate(year = as.integer(substring(name, 3,6))) %>%
  mutate(year = ifelse(year == 2000, 1998, year))

# plots of pH change over time
ggplot(PH_long_hab, aes(x = year, y = pH, group = REP_ID)) +
  geom_line(alpha = 0.5, aes(colour = change_dir) )+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.8,
              aes(colour = change_dir)) +
  facet_wrap(~BH07) + 
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  theme_dark()
ggsave("pH change over time facetted by habitat.png", path = "Outputs/Graphs/",
       width = 28, height = 15, units = "cm")

ggplot(PH_long_hab, aes(x = year, y = pH)) +
  geom_line(alpha = 0.5, aes(colour = change_dir, group = REP_ID))+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.5,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA) +
  facet_wrap(~BH07) + 
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "pH change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("pH change over time boxplots facetted by habitat.png", path = "Outputs/Graphs/",
       width = 28, height = 20, units = "cm")

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
ggsave("pH difference histograms facetted by habitat and year.png",
       path = "Outputs/Graphs/", width = 28, height = 24, units = "cm")



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

# Use habitat data from 2007 if it is there, otherwise 98 or 78
IBD_comb$BH_CODE <- ifelse(!is.na(IBD_comb$BH07), IBD_comb$BH07,
                           ifelse(!is.na(IBD_comb$BH98), IBD_comb$BH98,
                                  IBD_comb$BH78))
summary(IBD_comb$BH_CODE)
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
         diff0716 = R16 - R07) 
# Calculate overall Ell R change
ELL$Rchange <- rowSums(select(ELL, diff7898, diff9807, diff0716), na.rm = TRUE)
summary(ELL$Rchange)
filter(ELL, Rchange == 0) %>% select(starts_with("diff")) %>%
  summary()

ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff7807), ELL$diff7807, ELL$Rchange) 
ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff9807), ELL$diff9807, ELL$Rchange) 
ELL$Rchange <- ifelse(ELL$Rchange == 0 & !is.na(ELL$diff7898), ELL$diff7898, ELL$Rchange) 
summary(ELL$Rchange)

# Convert Ell R change into long format with better named habitats
ELL_diff_long <- ELL %>%
  select(REP_ID = REP_ID07, BH_CODE, starts_with("diff")) %>%
  left_join(select(BHPH_NAMES, BH_rec = BROAD_HABITAT, BH_CODE))%>%
  mutate(BH_rec = recode_factor(BH_rec,
                                "Boundary and Linear Features" = "Other",
                                "Littoral Sediment" = "Other",
                                "Sea" = "Other",
                                "Supra-littoral Rock" = "Other",
                                "Supra-littoral Sediment" = "Other",
                                "Littoral Rock" = "Other",
                                "Urban" = "Other",
                                "Standing Open Waters and Canals" = "Other",
                                "Rivers and Streams" = "Other",
                                "Montane" = "Other",
                                "Bracken" = "Other",
                                "No Allocation" = "Other",
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
  mutate(BH_rec = replace_na(BH_rec, "Other")) %>%
  droplevels %>%
  pivot_longer(starts_with("diff"), values_to = "Ell_R") %>%
  filter(!is.na(Ell_R))

ggplot(ELL_diff_long, aes(x = Ell_R)) + 
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(name ~ BH_rec, scales = "free_y")
ggsave("Ellenberg R change histograms facetted habitat year.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units ="cm")

# Convert Ellenberg R scores file to long format
ELL_R_LONG <- ELL %>%  select(REP_ID = REP_ID07,
                              PLOT_TYPE = PLOT_TYPE.x,
                              BH_CODE, R07, R98, R78, R16, Rchange) %>%
  left_join(select(BHPH_NAMES, BH_rec = BROAD_HABITAT, BH_CODE))%>%
  select(-BH_CODE) %>%
  filter(PLOT_TYPE == "X") %>%
  droplevels() %>%
  select(-PLOT_TYPE) %>%
  pivot_longer(cols = c(R07,R98,R78,R16), names_to = "year",
               names_prefix = "R") %>%
  mutate(year = ifelse(year == "07", 2007,
                       ifelse(year == "98", 1998,
                              ifelse(year == "78", 1978, 
                                     ifelse(year == "16", 2016, NA)))),
         BH_rec = recode_factor(BH_rec,
                                "Boundary and Linear Features" = "Other",
                                "Littoral Sediment" = "Coastal",
                                "Sea" = "Coastal",
                                "Supra-littoral Rock" = "Coastal",
                                "Supra-littoral Sediment" = "Coastal",
                                "Littoral Rock" = "Coastal",
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
  mutate(BH_rec = replace_na(BH_rec, "Other"))
str(ELL_R_LONG)  
summary(ELL_R_LONG$BH_rec)

ggplot(ELL_R_LONG, aes(x = year, y = value)) +
  geom_line(alpha = 0.5, aes(group = REP_ID, colour = Rchange))+
  geom_jitter(size = 0.2, width = 1, height = 0, shape = 16, alpha = 0.5,
              colour = "grey50") +
  geom_boxplot(fill= NA, aes(group = year), outlier.shape = NA, width = 3) +
  facet_wrap(~BH_rec) + 
  labs(y = "Ellenberg R") +
  # geom_smooth(formula = y ~ poly(x,3)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-1,1)*max(abs(PH_long_hab$change_dir)),
                         name = "Ell R change", na.value = "white") +
  # theme_dark() +
  NULL
ggsave("Ellenberg R change over time boxplots facetted by habitat.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units = "cm")



# Combined pH and ellenberg R ####
# change graphs - combine difference stats
ph_ell_comb <- ELL_diff_long %>%
  select(-BH_CODE) %>%
  filter(grepl("X", REP_ID)) %>%
  full_join(filter(PH_diff_long, name %in% unique(ELL_diff_long$name)))
str(ph_ell_comb)
mice::md.pattern(ph_ell_comb)

ggplot(ph_ell_comb, aes(x = pH, y = Ell_R)) +
  geom_point() +
  facet_grid(name~BH_rec) +
  geom_smooth(method = "lm")
ggsave("Ellenberg R vs pH change by habitat and survey.png",
       path = "Outputs/Graphs/", width = 28, height = 20, units = "cm")


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
            suffix = c("_ellR","_ph"), by = "REP_ID")
md.pattern(ph_ell_wide)

psych::pairs.panels(select_if(ph_ell_wide, is.numeric), 
                    ellipses = FALSE, rug = FALSE,
                    method = "spearman")
