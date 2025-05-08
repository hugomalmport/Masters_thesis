# Load necessary libraries
library(googlesheets4)
library(emmeans)
library(tidyverse)
library(agricolae)
library(ggh4x)

# Clean the environment
rm(list = ls())
graphics.off()

# Read all sheets
kelp_data_4 <- read_sheet("https://docs.google.com/spreadsheets/d/1FHVlBrssyievdxxegrkq6g7MS1HEL3AAirMl7akIjMo", 
                          sheet = "digitata")
kelp_data_5 <- read_sheet("https://docs.google.com/spreadsheets/d/1FHVlBrssyievdxxegrkq6g7MS1HEL3AAirMl7akIjMo", 
                          sheet = "saccharina & digitata")
kelp_data_6 <- read_sheet("https://docs.google.com/spreadsheets/d/1FHVlBrssyievdxxegrkq6g7MS1HEL3AAirMl7akIjMo", 
                          sheet = "hyperborea")
kelp_data_7 <- read_sheet("https://docs.google.com/spreadsheets/d/1FHVlBrssyievdxxegrkq6g7MS1HEL3AAirMl7akIjMo", 
                          sheet = "juvenile_hyperborea")

# Set species, location, and pigment color
species_colours <- c("L. digitata" = "seagreen3",
                     "S. latissima" = "skyblue1",
                     "L. hyperborea" = "orangered2")

# Function to apply all transformations to each dataset
# Add column treatment and experiment, and then factorize those and species
process_data <- function(df, experiment_name) {
  df <- df %>%
    mutate(
      # Add experiment column
      experiment = experiment_name
    ) %>%
    # Add treatment column
    mutate(treatment = paste0(salinity, "PSU & ", temperature, "°C")) %>%
    # Make sure both treatment and species are factors
    mutate(
      treatment = factor(treatment, 
                         levels = c("10PSU & 10°C", "10PSU & 17°C",
                                    "10PSU & 20°C", "20PSU & 10°C",
                                    "20PSU & 17°C", "20PSU & 20°C",
                                    "30PSU & 10°C", "30PSU & 20°C")),
      species = factor(species, 
                       levels = c("L. digitata", "S. latissima", "L. hyperborea")),
      location = factor(location,
                        levels = c("Ängklåvbukten", "Lökholmen", "Klövskär", "Klåvningarna"))
    )
  
  return(df)
}

# Apply the function to each dataset
kelp_data_4 <- process_data(kelp_data_4, "Experiment 4")
kelp_data_5 <- process_data(kelp_data_5, "Experiment 5")
kelp_data_6 <- process_data(kelp_data_6, "Experiment 6")
kelp_data_7 <- process_data(kelp_data_7, "Experiment 7")

# Mortality:
# Bind the data into one dataframe and filter for treatment, experiment and mortality
kelp_data_mortality <- bind_rows(kelp_data_4, kelp_data_5, kelp_data_6, kelp_data_7) %>%
  select(Mortality, experiment, treatment, species)

# Set baseline species and treatment
# I'll use 10PSU & 10°C as I want it in numerical order and L. digitata as baseline species
kelp_data_mortality$treatment <- relevel(kelp_data_mortality$treatment, ref = "10PSU & 10°C")  
kelp_data_mortality$species <- relevel(kelp_data_mortality$species, ref = "L. digitata")

# Count the number of mortalities per treatment, experiment, and species
mortality_counts <- kelp_data_mortality %>%
  filter(Mortality == "yes") %>%
  group_by(experiment, treatment, species) %>%
  summarise(count = n(), .groups = "drop")

# Plot
mortality <- ggplot(mortality_counts, aes(x = treatment, y = count, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = species_colours) +
  facet_grid(~experiment, scales = "free_x", space = "free_x") +  # Align treatments under each experiment
  scale_y_continuous(breaks = seq(0, max(mortality_counts$count, na.rm = TRUE), by = 2)) +  # Adjust y-axis
  labs(title = "Mortality",
       x = "Treatment",
       y = "Count",
       fill = "Species") +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "top",
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(mortality)

## Maybe adjust the requirement for mortality, depending on own requirements
## It's a matter of own definitions also!
## But it looks nice!

# Data processing and set up for statistical analysis
# Calculate SGR_total and replace NAs with mean SGR_total of treatment and species
mean_sgr <- function(df, experiment_id) {
  # Ensuring all weight columns exist in the dataset, adding NA if missing
  weight_cols <- paste0("weight", 1:5, "_g")
  missing_cols <- setdiff(weight_cols, colnames(df))
  df[missing_cols] <- NA
  
  df <- df %>%
    group_by(experiment, treatment, species) %>%
    mutate(
      # Adjusting SGR_total based on available final weight
      SGR_total = case_when(
        experiment == "Experiment 4" & !is.na(weight4_g) ~ (log(weight4_g) - log(weight1_g)) / 12 * 100,  # Experiment 4, 12 days
        experiment == "Experiment 5" & !is.na(weight5_g) ~ (log(weight5_g) - log(weight1_g)) / 14 * 100,  # Experiment 5, 14 days
        experiment %in% c("Experiment 6", "Experiment 7") & !is.na(weight5_g) ~ (log(weight5_g) - log(weight1_g)) / 16 * 100,  # Experiments 6 & 7, 16 days
        
        TRUE ~ NA_real_  # If no weight is available, return NA
      ),
      # Replace NA SGR_total values with the mean for that treatment and species only if SGR_total is NA
      SGR_total = ifelse(is.na(SGR_total), 
                         mean(SGR_total, na.rm = TRUE),
                         SGR_total)
    ) %>%
    ungroup()
  
  return(df)
}

# Now the same with weights 
## This needs to be done after so the SGR_total isn't calculated based on this
mean_weight <- function(df) {
  df <- df %>%
    group_by(experiment, treatment, species) %>%
    mutate(
      weight2_g = ifelse(is.na(weight2_g), mean(weight2_g, na.rm = TRUE), weight2_g),
      weight3_g = ifelse(is.na(weight3_g), mean(weight3_g, na.rm = TRUE), weight3_g),
      weight4_g = ifelse(is.na(weight4_g), mean(weight4_g, na.rm = TRUE), weight4_g),
      # Check if weight5_g exists and if NA values are present before calculating the mean
      weight5_g = if("weight5_g" %in% colnames(df)) {
        ifelse(is.na(weight5_g), mean(weight5_g, na.rm = TRUE), weight5_g)
      } else {
        weight4_g  # If weight5_g doesn't exist, use weight4_g as the last weight
      }
    ) %>%
    ungroup()
  
  return(df)
}

# Same function but with dry weight ratio instead!
mean_dw <- function(df) {
  df <- df %>%
    group_by(experiment, treatment, species) %>%
    mutate(
      DW_ratio = ifelse(is.na(DW_ratio), mean(DW_ratio, na.rm = TRUE), DW_ratio)
    ) %>%
    ungroup()
  
  return(df)
}

# And lastly the pigments
mean_pigments <- function(df) {
  # Ensure that necessary columns exist in the dataset, adding NA if missing
  pigment_cols <- c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g", "Total_carotenoid1_mg_g")
  missing_cols <- setdiff(pigment_cols, colnames(df))
  df[missing_cols] <- NA
  
  df <- df %>%
    group_by(experiment, treatment, species) %>%
    mutate(
      # Chlorophyll_a_mg_g: Replace NA with group mean
      Chlorophyll_a_mg_g = ifelse(is.na(Chlorophyll_a_mg_g),
                                  mean(Chlorophyll_a_mg_g, na.rm = TRUE),
                                  Chlorophyll_a_mg_g),
      
      # Chlorophyll_c_mg_g: Replace NA with group mean
      Chlorophyll_c_mg_g = ifelse(is.na(Chlorophyll_c_mg_g),
                                  mean(Chlorophyll_c_mg_g, na.rm = TRUE),
                                  Chlorophyll_c_mg_g),
      
      # Handling Total_carotenoid_mg_g: Only replace NA values with group mean
      Total_carotenoid_mg_g = ifelse(is.na(Total_carotenoid_mg_g),
                                     mean(Total_carotenoid_mg_g, na.rm = TRUE),
                                     Total_carotenoid_mg_g),
      
      # Handling Total_carotenoid1_mg_g (only for Experiment 4): Only replace NA values with group mean
      Total_carotenoid1_mg_g = ifelse(is.na(Total_carotenoid1_mg_g) & experiment == "Experiment 4",
                                      mean(Total_carotenoid1_mg_g, na.rm = TRUE),
                                      Total_carotenoid1_mg_g)
    ) %>%
    ungroup()
  
  return(df)
}

# Apply all four functions with the days argument for SGR calculation to all dataframes
# Experiment 4
kelp_data_4 <- kelp_data_4 %>%
  mean_sgr() %>%
  mean_weight() %>%
  mean_dw() %>%
  mean_pigments()

## Also there's this one sample (S1T1-5) that is deviating a lot in pigment values
### Mainly really high chlorophyll values even though it was dissolving
## So this is only for my data set
### Let's replace that one with mean pigment values
kelp_data_4 <- kelp_data_4 %>%
  mutate(
    # Replace S1T1-5 with treatment mean values for pigments
    Chlorophyll_a_mg_g = ifelse(label == "S1T1-5", mean(Chlorophyll_a_mg_g, na.rm = TRUE), Chlorophyll_a_mg_g),
    Chlorophyll_c_mg_g = ifelse(label == "S1T1-5", mean(Chlorophyll_c_mg_g, na.rm = TRUE), Chlorophyll_c_mg_g),
    Total_carotenoid_mg_g = ifelse(label == "S1T1-5", mean(Total_carotenoid1_mg_g, na.rm = TRUE), Total_carotenoid1_mg_g)
  )
  
# Experiment 5                   
kelp_data_5 <- kelp_data_5 %>%
  mean_sgr() %>%
  mean_weight() %>%
  mean_dw() %>%
  mean_pigments()

# Experiment 6
kelp_data_6 <- kelp_data_6 %>%
  mean_sgr() %>%
  mean_weight() %>%
  mean_dw() %>%
  mean_pigments()

# Experiment 7
kelp_data_7 <- kelp_data_7 %>%
  mean_sgr() %>%
  mean_weight() %>%
  mean_dw() %>%
  mean_pigments()

# Factorize all numeric factors used in the experiments
kelp_data_4$salinity <- factor(kelp_data_4$salinity, levels = c(10, 20, 30))
kelp_data_4$locationnumber <- factor(kelp_data_4$locationnumber, levels = c(1, 2, 3))
kelp_data_4$temperature <- factor(kelp_data_4$temperature, levels = c(10, 20))
kelp_data_5$salinity <- factor(kelp_data_5$salinity, levels = c(10, 20))
kelp_data_5$temperature <- factor(kelp_data_5$temperature, levels = c(10, 20))
kelp_data_5$locationnumber <- factor(kelp_data_5$locationnumber, levels = c(1, 2))
kelp_data_6$salinity <- factor(kelp_data_6$salinity, levels = c(10, 20))
kelp_data_6$temperature <- factor(kelp_data_6$temperature, levels = c(10, 20))
kelp_data_6$locationnumber <- factor(kelp_data_6$locationnumber, levels = c(4))
kelp_data_7$salinity <- factor(kelp_data_7$salinity, levels = c(10, 20))
kelp_data_7$temperature <- factor(kelp_data_7$temperature, levels = c(10, 17))
kelp_data_7$locationnumber <- factor(kelp_data_7$locationnumber, levels = c(4))

# Combine all datasets into one
kelp_data_combined <- bind_rows(kelp_data_4, kelp_data_5, kelp_data_6, kelp_data_7)

# Set a baseline of treatment and species
## Same one here :)
kelp_data_combined$treatment <- relevel(kelp_data_combined$treatment, ref = "10PSU & 10°C")  
kelp_data_combined$species <- relevel(kelp_data_combined$species, ref = "L. digitata")
kelp_data_combined$location <- relevel(kelp_data_combined$location, ref = "Ängklåvbukten")

## Let's figure out what the plots should look like based on our factors 
### And which are significant - cause no need to include things that don't matter anyway
# Statistical analysis: two-, three-, and four-factor ANOVAs
## Starting with Specific growth rate

# Experiment 4
sgr_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
sgr_model_4 <- aov(SGR_total ~ salinity * temperature * locationnumber, data = sgr_data_4)

summary(sgr_model_4)

# Perform SNK tests using existing ANOVA models, only use the significant factors!!
## These letters will then be added to the plots to differentiate the groups
# And since the interaction is significant, we need to create an interaction term
sgr_data_4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  mutate(sal_temp = interaction(salinity, temperature))

sgr_model_4 <- aov(SGR_total ~ sal_temp, data = sgr_data_4)

sgr_snk_4 <- SNK.test(sgr_model_4, "sal_temp")
sgr_snk_4$groups

## For experiment 5 we do separate ANOVAs for the two species

# Experiment 5 - L. digitata
sgr_dig_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "L. digitata")
sgr_dig_5 <- aov(SGR_total ~ salinity * temperature * locationnumber, data = sgr_dig_5)

summary(sgr_dig_5)

# Experiment 5 - S. latissima
sgr_sac_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "S. latissima")
sgr_sac_5 <- aov(SGR_total ~ salinity * temperature * locationnumber, data = sgr_sac_5)

summary(sgr_sac_5)

# Perform SNK tests
# L. digitata
sgr_dig_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "L. digitata") %>%
  mutate(sal_temp = interaction(salinity, temperature))

sgr_snk_dig_5 <- aov(SGR_total ~ sal_temp, data = sgr_dig_5)
sgr_snk_dig_5 <- SNK.test(sgr_snk_dig_5, "sal_temp")
sgr_snk_dig_5$groups

# S. latissima
sgr_sac_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "S. latissima") %>%
  mutate(sal_temp = interaction(salinity, temperature))

sgr_snk_sac_5 <- aov(SGR_total ~ sal_temp, data = sgr_sac_5)
sgr_snk_sac_5 <- SNK.test(sgr_snk_sac_5, "sal_temp")
sgr_snk_sac_5$groups

# Experiment 6
sgr_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
sgr_model_6 <- aov(SGR_total ~ salinity * temperature, data = sgr_data_6)

summary(sgr_model_6)

# Perform SNK test
sgr_data_6 <- kelp_data_combined %>%
  filter(experiment == "Experiment 6") %>%
  mutate(sal_temp = interaction(salinity, temperature))

sgr_model_6_int <- aov(SGR_total ~ sal_temp, data = sgr_data_6)
sgr_snk_6 <- SNK.test(sgr_model_6_int, "sal_temp")
sgr_snk_6$groups

# Experiment 7
sgr_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
sgr_model_7 <- aov(SGR_total ~ salinity * temperature, data = sgr_data_7)

summary(sgr_model_7)

# Perform SNK test
sgr_snk_7 <- SNK.test(sgr_model_7, c("salinity", "temperature"))
sgr_snk_7$groups

## Dry to wet weight ratio

# Experiment 4
dw_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
dw_model_4 <- aov(DW_ratio ~ salinity * temperature * locationnumber, data = dw_data_4)

summary(dw_model_4)

# No significance so no need to do SNK test here!

# Experiment 5 - L. digitata
dw_dig_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "L. digitata")
dw_dig_5 <- aov(DW_ratio ~ salinity * temperature * locationnumber, data = dw_dig_5)

summary(dw_dig_5)

# Experiment 5 - S. latissima
dw_sac_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "S. latissima")
dw_sac_5 <- aov(DW_ratio ~ salinity * temperature * locationnumber, data = dw_sac_5)

summary(dw_sac_5)

# Perform SNK tests
# L. digitata - temperature:location was almost significant, so let's do that here
dw_dig_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "L. digitata") %>%
  mutate(temp_loc = interaction(temperature, locationnumber))

dw_snk_dig_5 <- aov(DW_ratio ~ temp_loc, data = dw_dig_5)
dw_snk_dig_5 <- SNK.test(dw_snk_dig_5, "temp_loc")
dw_snk_dig_5$groups

# S. latissima
dw_sac_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "S. latissima") %>%
  mutate(sal_loc = interaction(salinity, locationnumber))

dw_snk_sac_5 <- aov(DW_ratio ~ sal_loc, data = dw_sac_5)
dw_snk_sac_5 <- SNK.test(dw_snk_sac_5, "sal_loc")
dw_snk_sac_5$groups

## Alright so no difference between groups even though the interaction is significant, interesting

# Experiment 6
dw_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
dw_model_6 <- aov(DW_ratio ~ salinity * temperature, data = dw_data_6)

summary(dw_model_6)

# Perform SNK tests
dw_model_6 <- kelp_data_combined %>%
  filter(experiment == "Experiment 6") %>%
  mutate(sal_temp = interaction(salinity, temperature))

dw_snk_6 <- aov(DW_ratio ~ sal_temp, data = dw_model_6)
dw_snk_6 <- SNK.test(dw_snk_6, "sal_temp")
dw_snk_6$groups

# Experiment 7
dw_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
dw_model_7 <- aov(DW_ratio ~ salinity * temperature, data = dw_data_7)

summary(dw_model_7)

# Perform SNK tests
dw_snk_7 <- SNK.test(dw_model_7, c("temperature"))
dw_snk_7$groups

## Chlorophyll and carotenoid statistics!!!
## Let's do very pigment for itself

# Experiment 4
# Chlorophyll a
chla_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
chla_model_4 <- aov(Chlorophyll_a_mg_g ~ salinity * temperature * locationnumber, data = chla_data_4)

summary(chla_model_4)

# Perform SNK tests
chla_snk_4 <- SNK.test(chla_model_4, c("temperature"))
chla_snk_4$groups

# Chlorophyll c
chlc_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
chlc_model_4 <- aov(Chlorophyll_c_mg_g ~ salinity * temperature * locationnumber, data = chlc_data_4)

summary(chlc_model_4)

# Nothing signficant, no need for SNK test

# Carotenoids
car_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
car_model_4 <- aov(Total_carotenoid1_mg_g ~ salinity * temperature * locationnumber, data = car_data_4)

summary(car_model_4)

# Perform SNK test
car_data_4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  mutate(sal_temp_loc = interaction(salinity, temperature, locationnumber))

car_model_4 <- aov(Total_carotenoid1_mg_g ~ sal_temp_loc, data = car_data_4)
snk_car_4 <- SNK.test(car_model_4, "sal_temp_loc")
snk_car_4$groups

# Experiment 5
# Chlorophyll a - L. digitata
chla_dig_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "L. digitata")
chla_dig_5 <- aov(Chlorophyll_a_mg_g ~ salinity * temperature * locationnumber, data = chla_dig_5)

summary(chla_dig_5)

# Chlorophyll a - S. latissima
chla_sac_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "S. latissima")
chla_sac_5 <- aov(Chlorophyll_a_mg_g ~ salinity * temperature * locationnumber, data = chla_sac_5)

summary(chla_sac_5)

# Perform SNK test
# L. digitata
chla_dig_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "L. digitata") %>%
  mutate(sal_temp = interaction(salinity, temperature))

chla_snk_dig_5 <- aov(Chlorophyll_a_mg_g ~ sal_temp, data = chla_dig_5)
chla_snk_dig_5 <- SNK.test(chla_snk_dig_5, "sal_temp")
chla_snk_dig_5$groups

# S. latissima
chla_sac_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "S. latissima") %>%
  mutate(temp_loc = interaction(temperature, locationnumber))

chla_snk_sac_5 <- aov(Chlorophyll_a_mg_g ~ temp_loc, data = chla_sac_5)
chla_snk_sac_5 <- SNK.test(chla_snk_sac_5, "temp_loc")
chla_snk_sac_5$groups

# Chlorophyll c - L. digitata
chlc_dig_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "L. digitata")
chlc_dig_5 <- aov(Chlorophyll_c_mg_g ~ salinity * temperature * locationnumber, data = chlc_dig_5)

summary(chlc_dig_5)

# Experiment 5
# Chlorophyll c - S. latissima
chlc_sac_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "S. latissima")
chlc_sac_5 <- aov(Chlorophyll_c_mg_g ~ salinity * temperature * locationnumber, data = chlc_sac_5)

summary(chlc_sac_5)

# Perform SNK test
# None for L. digitata
# Salinity and temperature as main factors for S. latissima
chlc_snk_sac_5 <- SNK.test(chlc_sac_5, c("salinity", "temperature"))
chlc_snk_sac_5$groups

# Carotenoids - L. digitata
car_dig_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "L. digitata")
car_dig_5 <- aov(Total_carotenoid_mg_g ~ salinity * temperature * locationnumber, data = car_dig_5)

summary(car_dig_5)

# Carotenoids - S. latissima
car_sac_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5", species == "S. latissima")
car_sac_5 <- aov(Total_carotenoid_mg_g ~ salinity * temperature * locationnumber, data = car_sac_5)

summary(car_sac_5)

# Perform SNK test
# L. digitata
car_dig_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "L. digitata") %>%
  mutate(sal_temp = interaction(salinity, temperature))

car_snk_dig_5 <- aov(Total_carotenoid_mg_g ~ sal_temp, data = car_dig_5)
car_snk_dig_5 <- SNK.test(car_snk_dig_5, "sal_temp")
car_snk_dig_5$groups

# S. latissima
car_sac_5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5", species == "S. latissima") %>%
  mutate(temp_loc = interaction(temperature, locationnumber))

car_snk_sac_5 <- aov(Total_carotenoid_mg_g ~ temp_loc, data = car_sac_5)
car_snk_sac_5 <- SNK.test(car_snk_sac_5, "temp_loc")
car_snk_sac_5$groups

# Experiment 6
# Chlorophyll a
chla_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
chla_model_6 <- aov(Chlorophyll_a_mg_g ~ salinity * temperature, data = chla_data_6)

summary(chla_model_6)

# Perform SNK test
chla_snk_6 <- SNK.test(chla_model_6, c("salinity", "temperature"))
chla_snk_6$groups

# Chlorophyll c
chlc_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
chlc_model_6 <- aov(Chlorophyll_c_mg_g ~ salinity * temperature, data = chlc_data_6)

summary(chlc_model_6)

# Perform SNK test
chlc_snk_6 <- SNK.test(chlc_model_6, c("salinity"))
chlc_snk_6$groups

# Carotenoids
car_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
car_model_6 <- aov(Total_carotenoid_mg_g ~ salinity * temperature, data = car_data_6)

summary(car_model_6)

# Experiment 7
# Chlorophyll a
chla_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
chla_model_7 <- aov(Chlorophyll_a_mg_g ~ salinity * temperature, data = chla_data_7)

summary(chla_model_7)

# Perform SNK test
chla_data_7 <- kelp_data_combined %>%
  filter(experiment == "Experiment 7") %>%
  mutate(sal_temp = interaction(salinity, temperature))

chla_snk_7 <- aov(Chlorophyll_a_mg_g ~ sal_temp, data = chla_data_7)
chla_snk_7 <- SNK.test(chla_snk_7, "sal_temp")
chla_snk_7$groups

# Chlorophyll c
chlc_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
chlc_model_7 <- aov(Chlorophyll_c_mg_g ~ salinity * temperature, data = chlc_data_7)

summary(chlc_model_7)

# Perform SNK test
chlc_snk_7 <- SNK.test(chlc_model_7, c("salinity"))
chlc_snk_7$groups

# Carotenoids
car_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
car_model_7 <- aov(Total_carotenoid_mg_g ~ salinity * temperature, data = car_data_7)

summary(car_model_7)

# Perform SNK test
car_snk_7 <- SNK.test(car_model_7, c("salinity"))
car_snk_7$groups

# PLOTS
## Line plot of how the weight changes over time!
# Global factor order definition
factor_order <- c(
  # Salinity & Temperature
  "10PSU & 10°C", "10PSU & 17°C", "10PSU & 20°C", "20PSU & 10°C", "20PSU & 17°C", "20PSU & 20°C", "30PSU & 10°C", "30PSU & 20°C",
  # Temperature & Location
  "10°C & Ängklåvbukten", "10°C & Lökholmen", "20°C & Ängklåvbukten", "20°C & Lökholmen", "30°C & Ängklåvbukten", "30°C & Lökholmen",
  # Salinity & Location
  "10PSU & Ängklåvbukten", "10PSU & Lökholmen", "20PSU & Ängklåvbukten", "20PSU & Lökholmen", "30PSU & Ängklåvbukten", "30PSU & Lökholmen"
)

# Summarize mean weight and standard error
weight_summary_mean <- kelp_data_combined %>%
  pivot_longer(cols = c(weight1_g, weight2_g, weight3_g, weight4_g, weight5_g), 
               names_to = "day", 
               values_to = "weight") %>%
  mutate(days = case_when(
    day == "weight1_g" ~ 0,
    day == "weight2_g" ~ 4,
    day == "weight3_g" & experiment == "Experiment 5" ~ 7,
    day == "weight3_g" & experiment %in% c("Experiment 4", "Experiment 6", "Experiment 7") ~ 8,
    day == "weight4_g" & experiment == "Experiment 5" ~ 10,
    day == "weight4_g" & experiment %in% c("Experiment 4", "Experiment 6", "Experiment 7") ~ 12,
    day == "weight5_g" & experiment == "Experiment 5" ~ 14,
    day == "weight5_g" & experiment %in% c("Experiment 6", "Experiment 7") ~ 16,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(days)) %>%
  group_by(days, treatment, experiment) %>%
  summarise(
    mean_weight = mean(weight, na.rm = TRUE),
    se_weight = sd(weight, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Create the plot with straight trendlines
weight_plot <- ggplot(weight_summary_mean, aes(x = days, y = mean_weight, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = mean_weight - se_weight, ymax = mean_weight + se_weight), width = 0.4, linewidth = 0.7) +  
  geom_point(size = 2, alpha = 0.8) +  
  geom_line(linewidth = 1) +
  facet_wrap(~experiment, scales = "free_y") +  
  labs(
    title = "Weight change over time",
    x = "Days",
    y = "Mean Weight (g)",
    color = "Treatment"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  ) +
  scale_color_manual(values = c(
    "30PSU & 20°C" = "darkmagenta",
    "30PSU & 10°C" = "darkseagreen",
    "20PSU & 20°C" = "brown1",
    "20PSU & 10°C" = "cornflowerblue",
    "10PSU & 20°C" = "burlywood1",
    "10PSU & 10°C" = "darkslategray1",
    "10PSU & 17°C" = "burlywood1",
    "20PSU & 17°C" = "brown1"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

# Display the plot
print(weight_plot)

## SGR plot!
# Calculate mean SGR_total and standard error for each treatment
sgr_summary <- kelp_data_combined %>%
  group_by(temperature, salinity, treatment, species, experiment) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),  # Calculate mean SGR_total
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Add total days the experiments ran for and create experiment labels
sgr_summary <- sgr_summary %>%
  mutate(experiment_label = case_when(
    experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
    experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
    experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
    experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
    experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
  ))


# Create the plot with only SGR_total
sgr_barplot <- ggplot(sgr_summary, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
  scale_fill_manual(values = species_colours) +
  geom_errorbar(
    aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  labs(
    title = "Specific Growth rate",
    x = "Treatment",
    y = expression("SGR (% day"^-1*")"),
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(sgr_barplot)

# DW-ratio!!
# Add total days and create a customized x-label
dw_summary <- kelp_data_combined %>%
  mutate(
    experiment_label = case_when(
      experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
      experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
      experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
      experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
      experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
    ),
    x_label = case_when(
      experiment == "Experiment 5" & species == "L. digitata" ~ paste0(temperature, "°C & ", location),
      experiment == "Experiment 5" & species == "S. latissima" ~ paste0(salinity, "PSU & ", location),
      TRUE ~ paste0(salinity, "PSU & ", temperature, "°C")
    )
  ) %>%
  group_by(experiment_label, x_label, species) %>%
  summarise(
    mean_DW_ratio = mean(DW_ratio, na.rm = TRUE),
    se_DW_ratio = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Apply the custom label
dw_summary$x_label <- factor(dw_summary$x_label, levels = factor_order)

# Create the plot
dw_barplot <- ggplot(dw_summary, aes(x = x_label, y = mean_DW_ratio, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_DW_ratio - se_DW_ratio, ymax = mean_DW_ratio + se_DW_ratio), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  labs(
    title = "Dry to wet weight-ratio",
    x = "Treatment",
    y = "Mean D:W ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(dw_barplot)

# Chlorophyll a time!
# Same here, apply number of days and add custom label
chla_summary <- kelp_data_combined %>%
  mutate(
    location = factor(location, levels = c("Ängklåvbukten", "Lökholmen", "Klövskär", "Klåvningarna")),
    experiment_label = case_when(
      experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
      experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
      experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
      experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
      experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
    ),
    x_label = case_when(
      experiment == "Experiment 5" & species == "L. digitata" ~ paste0(salinity, "PSU & ", temperature, "°C"),
      experiment == "Experiment 5" & species == "S. latissima" ~ paste0(temperature, "°C & ", location),
      TRUE ~ paste0(salinity, "PSU & ", temperature, "°C")
    ),
    x_label = factor(x_label, levels = unique(x_label))  # Convert x_label to a factor with explicit order
  ) %>%
  group_by(experiment_label, x_label, species) %>%
  summarise(
    mean_chla = mean(`Chlorophyll_a_mg_g`, na.rm = TRUE),
    se_chla = sd(`Chlorophyll_a_mg_g`, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Apply the custom label
chla_summary$x_label <- factor(chla_summary$x_label, levels = factor_order)

# Create the plot
chla_barplot <- ggplot(chla_summary, aes(x = x_label, y = mean_chla, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_chla - se_chla, ymax = mean_chla + se_chla), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  labs(
    title = "Chlorophyll a",
    x = "Treatment",
    y = "Mean Chlorophyll a (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(chla_barplot)

# Chlorophyll c time!
# Calculate mean Chlorophyll c and standard error
chlc_summary <- kelp_data_combined %>%
  group_by(treatment, species, experiment) %>%
  summarise(
    mean_chlc = mean(`Chlorophyll_c_mg_g`, na.rm = TRUE),  
    se_chlc = sd(`Chlorophyll_c_mg_g`, na.rm = TRUE) / sqrt(n()),  
    .groups = "drop"
  )

# Add number of days, no need for the custom label here
chlc_summary <- chlc_summary %>%
  mutate(experiment_label = case_when(
    experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
    experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
    experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
    experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
    experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
  ))

# Create the plot
chlc_barplot <- ggplot(chlc_summary, aes(x = treatment, y = mean_chlc, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_chlc - se_chlc, ymax = mean_chlc + se_chlc), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  labs(
    title = "Chlorophyll c1 + c2",
    x = "Treatment",
    y = "Mean Chlorophyll c (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(chlc_barplot)

# Carotenoids!
# First alternative:
# Add number of days and custom label :)
car_summary <- kelp_data_combined %>%
  mutate(
    location = factor(location, levels = c("Ängklåvbukten", "Lökholmen", "Klövskär", "Klåvningarna")),
    experiment_label = case_when(
      experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
      experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
      experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
      experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
      experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
    ),
    x_label = case_when(
      experiment == "Experiment 5" & species == "L. digitata" ~ paste0(salinity, "PSU & ", temperature, "°C"),
      experiment == "Experiment 5" & species == "S. latissima" ~ paste0(temperature, "°C & ", location),
      TRUE ~ paste0(salinity, "PSU & ", temperature, "°C")
    )
  ) %>%
  group_by(experiment_label, x_label, species) %>%
  summarise(
    mean_car = mean(`Total_carotenoid_mg_g`, na.rm = TRUE),
    se_car = sd(`Total_carotenoid_mg_g`, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Apply the custom label!
car_summary$x_label <- factor(car_summary$x_label, levels = factor_order)

# Create the plot
car_barplot1 <- ggplot(car_summary, aes(x = x_label, y = mean_car, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_car - se_car, ymax = mean_car + se_car), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  labs(
    title = "Total Carotenoids",
    x = "Treatment",
    y = "Mean carotenoid content (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(car_barplot1)

# Second alternative:
# Add number of days and custom label
car_summary <- kelp_data_combined %>%
  mutate(
    experiment_label = case_when(
      experiment == "Experiment 4" ~ "Experiment 4\n(12 days)",
      experiment == "Experiment 5" & species == "L. digitata" ~ "Experiment 5\n(14 days)",
      experiment == "Experiment 5" & species == "S. latissima" ~ "Experiment 5\n(14 days)\u00A0",
      experiment == "Experiment 6" ~ "Experiment 6\n(16 days)",
      experiment == "Experiment 7" ~ "Experiment 7\n(16 days)"
    ),
    location_label = ifelse(experiment == "Experiment 4", as.character(location), ""),
    x_label = case_when(
      experiment == "Experiment 5" & species == "S. latissima" ~ paste0(temperature, "°C & ", location),
      experiment == "Experiment 4" ~ paste0(salinity, "PSU & ", temperature, "°C"),
      TRUE ~ paste0(salinity, "PSU & ", temperature, "°C")
    )
  ) %>%
  mutate(
    location_label = factor(location_label, 
                            levels = c("Ängklåvbukten", "Lökholmen", "Klövskär", "Klåvningarna", ""))
  ) %>%
  group_by(experiment_label, location_label, x_label, species) %>%
  summarise(
    mean_car = mean(`Total_carotenoid_mg_g`, na.rm = TRUE),
    se_car = sd(`Total_carotenoid_mg_g`, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Apply the x-axis order
car_summary$x_label <- factor(car_summary$x_label, levels = factor_order)

# Create the barplot
car_barplot2 <- ggplot(car_summary, aes(x = x_label, y = mean_car, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_car - se_car, ymax = mean_car + se_car), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_nested(
    . ~ experiment_label + location_label, 
    scales = "free_x",
    nest_line = TRUE,
    remove_labels = "y",
    drop = TRUE 
  ) +
  labs(
    title = "Total Carotenoids",
    x = "Treatment",
    y = "Mean carotenoid content (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display the plot
print(car_barplot2)


# Third alternative:
# Here we'll create two separate plots, one for Experiment 4 and then one for the rest
# Find the common y-axis limit for all experiments, so that we have the same for both plots
y_max <- max(car_summary$mean_car + car_summary$se_car, na.rm = TRUE)

# Plot for Experiment 4
car_summary_exp4 <- car_summary %>%
  filter(experiment_label == "Experiment 4\n(12 days)")

# Apply the custom x-labels
car_summary_exp4$x_label <- factor(car_summary_exp4$x_label, levels = factor_order)

# Create the plot for Experiment 4
car_barplot3_1 <- ggplot(car_summary_exp4, aes(x = x_label, y = mean_car, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_car - se_car, ymax = mean_car + se_car),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_nested(. ~ location_label, scales = "free_x", nest_line = TRUE) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(
    title = "Experiment 4: Total Carotenoids",
    x = "Treatment",
    y = "Mean carotenoid content (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Plot for Experiments 5 – 7
car_summary_exp5to7 <- car_summary %>%
  filter(experiment_label != "Experiment 4\n(12 days)")

# Apply the custom x-labels
car_summary_exp5to7$x_label <- factor(car_summary_exp5to7$x_label, levels = factor_order)

# Create the plot
car_barplot3_2 <- ggplot(car_summary_exp5to7, aes(x = x_label, y = mean_car, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_car - se_car, ymax = mean_car + se_car),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment_label, scales = "free_x") +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(
    title = "Experiments 5 – 7: Total Carotenoids",
    x = "Treatment",
    y = "Mean carotenoid content (mg/g)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text.x = element_text(face = "bold")
  )

# Display plots separately
print(car_barplot3_1)
print(car_barplot3_2)


# Salinity fluctuations
# Identify all salinity measurement columns
salinity_cols <- grep("^sal_", colnames(kelp_data_combined), value = TRUE)
salinity_cols <- sort(salinity_cols, decreasing = FALSE) # Ensure correct order

# Pivot to longer format
kelp_salinity <- kelp_data_combined %>%
  pivot_longer(
    cols = all_of(salinity_cols),  
    names_to = "timepoint",
    values_to = "salinity_value"
  ) %>%
  mutate(
    timepoint = gsub("sal_", "", timepoint),  
    timepoint = as.numeric(timepoint),  
    salinity_value = as.numeric(salinity_value)
  ) %>%
  filter(!is.na(timepoint) & !is.na(salinity_value))

# Create the plot for salinity fluctuations!!
sal_plot <- ggplot(kelp_salinity, aes(x = timepoint, y = salinity_value, 
                                      group = individual, color = treatment)) +  
  geom_line(alpha = 1) +
  facet_wrap(~ experiment) +
  scale_color_manual(
    values = c(
      "30PSU & 10°C" = "deepskyblue", "30PSU & 20°C" = "firebrick",
      "20PSU & 10°C" = "deepskyblue", "20PSU & 20°C" = "firebrick",
      "10PSU & 10°C" = "deepskyblue", "10PSU & 20°C" = "firebrick",
      "10PSU & 17°C" = "red", "20PSU & 17°C" = "red"
    )
  ) + 
  labs(
    title = "Salinity fluctuations",
    x = "Time (hours)", y = "Salinity (psu)", color = "Treatment") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text.x = element_text(face = "bold"))

# Display the plot
print(sal_plot)

## Saving the plots!!
# Define output directory
output_dir <- "C:/Users/hugom/Desktop/master thesis/Programming/Plots"

# Common width, height and dpi settings:
w <- 10      # width in inches
h <- 7      # height in inches
res <- 500  # dpi

# Save the plots in the new order
ggsave(file.path(output_dir, "mortality.png"), plot = mortality, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "weight_plot.png"), plot = weight_plot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sgr_barplot.png"), plot = sgr_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "dw_barplot.png"), plot = dw_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "chla_barplot.png"), plot = chla_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "chlc_barplot.png"), plot = chlc_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "car_barplot1.png"), plot = car_barplot1, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "car_barplot2.png"), plot = car_barplot2, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "car_barplot3_1.png"), plot = car_barplot3_1, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "car_barplot3_2.png"), plot = car_barplot3_2, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sal_plot.png"), plot = sal_plot, width = w, height = h, dpi = res)


stop()---------------- :)


## Rests - leave in case we need to look into experiments one by one

# Run ANOVA with the combined treatment
sgr_data_4$treatment <- with(sgr_data_4, interaction(salinity, temperature, locationnumber, sep = ":"))
sgr_model_4_treat <- aov(SGR_total ~ treatment, data = sgr_data_4)

# Tukey HSD with agricolae
tukey_test <- HSD.test(sgr_model_4_treat, "treatment", group = TRUE)

# Get group letters dataframe
groups_df <- tukey_test$groups
groups_df$treatment <- rownames(groups_df)
groups_df <- groups_df[, c("treatment", "groups")]

# Perform the SNK test for salinity * temperature * locationnumber
snk_sgr4 <- SNK.test(sgr_model_4, c("salinity", "temperature", "locationnumber"))

# Extract the SNK results for plotting
snk_sgr4 <- snk_sgr4$groups

# Calculate mean SGR_total and standard error for each treatment
sgr_exp4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = location,
    group_id = paste(salinity, temperature, locationnumber, sep = ":")  # Group ID to join SNK
  )

# Format SNK results for joining
sgr_snk4 <- snk_sgr4 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id")

# Merge SNK results into plot data
sgr_exp4 <- left_join(sgr_exp4, sgr_snk4, by = "group_id")

# Define positions for vertical lines (between salinities)
sal_break4 <- c(2.5, 4.5)  # Adjust based on the number of salinity levels

# Create the plot
sgr_plot4 <- ggplot(sgr_exp4, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  # Add vertical lines between different salinities
  geom_vline(xintercept = sal_break4, linetype = "dashed", color = "gray50") +
  facet_grid(. ~ location_label) +
  # Add the SNK letters above the bars
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) + 
  labs(
    title = "SGR - Experiment 4",
    x = "Treatment",
    y = expression("Mean SGR (% day"^-1*")"),
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")  # Add spacing between temperature facets
  )

# Display the plot
print(sgr_plot4)

## Let's move onto experiment 5!
# Experiment 5
# Perform SNK tests for each species in Experiment 5
# L. digitata
# Perform SNK tests using existing ANOVA models
snk_dig_5 <- SNK.test(sgr_dig_5, c("salinity", "temperature", "locationnumber"))
snk_sac_5 <- SNK.test(sgr_sac_5, c("salinity", "temperature", "locationnumber"))

# Extract and modify SNK groups
snk_dig_df <- snk_dig_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. digitata")

snk_sac_df <- snk_sac_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(
    species = "S. latissima",
    groups = chartr("abc", "def", groups)  # Shift letters to avoid overlap
  )

# Combine SNK groupings
snk_combined_5 <- bind_rows(snk_dig_df, snk_sac_df)

# Calculate mean SGR_total and standard error for each treatment
sgr_exp5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = location,
    group_id = paste(salinity, temperature, locationnumber, sep = ":")
  )

# Merge SNK results into plot data
sgr_exp5 <- left_join(sgr_exp5, snk_combined_5, by = c("group_id", "species"))

# Define position for vertical line between salinity levels
sal_break <- c(2.5)

# Create the plot
sgr_plot5 <- ggplot(sgr_exp5, aes(x = treatment, y = mean_SGR, fill = location_label)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = location_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  facet_grid(. ~ species) +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "SGR - Experiment 5",
    x = "Treatment",
    y = expression("Mean SGR (% day"^-1*")"),
    fill = "Location"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(sgr_plot5)

# Experiment 6
# Perform SNK test for Experiment 6
snk_result_6 <- SNK.test(sgr_model_6, c("salinity", "temperature"))

# Extract the SNK results for plotting
snk_sgr6 <- snk_result_6$groups
snk_sgr6

# Prepare data for plotting
sgr_exp6 <- kelp_data_combined %>%
  filter(experiment == "Experiment 6") %>%
  group_by(temperature, salinity, treatment, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    group_id = paste(salinity, temperature, sep = ":")  # Group ID for SNK matching
  )

# Format SNK results for joining
snk_df6 <- snk_sgr6 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id")

# Merge SNK results into plot data
sgr_exp6 <- left_join(sgr_exp6, snk_df6, by = "group_id")

# Create the plot for Experiment 6
sgr_plot6 <- ggplot(sgr_exp6, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  # Add SNK letters above the bars
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) + 
  labs(
    title = "SGR - Experiment 6",
    x = "Treatment",
    y = expression("Mean SGR (% day"^-1*")"),
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")  # Add spacing between temperature facets
  )

# Display the plot
print(sgr_plot6)

# Experiment 7
# Perform SNK test for Experiment 7
snk_result_7 <- SNK.test(sgr_model_7, c("salinity", "temperature"))

# Extract the SNK results for plotting
snk_sgr7 <- snk_result_7$groups
snk_sgr7

# Prepare data for plotting
sgr_exp7 <- kelp_data_combined %>%
  filter(experiment == "Experiment 7") %>%
  group_by(temperature, salinity, treatment, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    group_id = paste(salinity, temperature, sep = ":")  # Group ID for SNK matching
  )

# Format SNK results for joining
snk_df7 <- snk_sgr7 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id")

# Merge SNK results into plot data
sgr_exp7 <- left_join(sgr_exp7, snk_df7, by = "group_id")

# Create the plot for Experiment 7
sgr_plot7 <- ggplot(sgr_exp7, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  # Add SNK letters above the bars
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) + 
  labs(
    title = "SGR - Experiment 7",
    x = "Treatment",
    y = expression("Mean SGR (% day"^-1*")"),
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")  # Add spacing between temperature facets
  )

# Display the plot
print(sgr_plot7)

## Dry:wet weight ratio plots!!
## Let's check at experiment 4 and 5 separately - with location included!

# Experiment 4
# Perform the SNK test for Experiment 4 (Dry-to-wet weight ratio)
snk_dw_4 <- SNK.test(dw_model_4, c("salinity", "temperature", "locationnumber"))

# Extract the SNK results for plotting
snk_dw_4_groups <- snk_dw_4$groups

# Calculate mean and standard error for each treatment
dw_exp4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = paste0(location),
    group_id = paste(salinity, temperature, locationnumber, sep = ":")  # Group ID for joining SNK
  )

# Format SNK results for joining
dw_snk4 <- snk_dw_4_groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id")

# Merge SNK results into plot data
dw_exp4 <- left_join(dw_exp4, dw_snk4, by = "group_id")

# Create the plot with SNK letters
dw_plot4 <- ggplot(dw_exp4, aes(x = treatment, y = mean_DW, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break4, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  # Add the SNK letters above the bars
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) + 
  labs(
    title = "Dry to Wet Weight Ratio - Experiment 4",
    x = "Treatment",
    y = "Mean D:W Ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")  # Add spacing between temperature facets
  )

# Display the plot
print(dw_plot4)

# Experiment 5
# Perform SNK tests using existing ANOVA models
snk_dw_dig_5 <- SNK.test(dw_dig_5, c("salinity", "temperature", "locationnumber"))
snk_dw_sac_5 <- SNK.test(dw_sac_5, c("salinity", "temperature", "locationnumber"))

# Extract and modify SNK groups
snk_dw_dig_df <- snk_dw_dig_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. digitata")

snk_dw_sac_df <- snk_dw_sac_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(
    species = "S. latissima",
    groups = chartr("abc", "def", groups)  # Shift letters to avoid overlap
  )

# Combine SNK groupings
snk_dw_combined_5 <- bind_rows(snk_dw_dig_df, snk_dw_sac_df)

# Calculate mean DW_ratio and standard error for each treatment
dw_exp5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = location,
    group_id = paste(salinity, temperature, locationnumber, sep = ":")
  )

# Merge SNK results into plot data
dw_exp5 <- left_join(dw_exp5, snk_dw_combined_5, by = c("group_id", "species"))

# Create the plot
dw_plot5 <- ggplot(dw_exp5, aes(x = treatment, y = mean_DW, fill = location_label)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = location_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  facet_grid(. ~ species) +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Dry to Wet Weight Ratio - Experiment 5",
    x = "Treatment",
    y = "Mean D:W Ratio",
    fill = "Location"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(dw_plot5)

# Experiment 6
# Perform SNK test using existing ANOVA model for Experiment 6
snk_dw_6 <- SNK.test(dw_model_6, c("salinity", "temperature"))

# Extract and modify SNK groups
snk_dw_6_df <- snk_dw_6$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. hyperborea")  # Or adjust based on your species for this experiment

# Calculate mean DW_ratio and standard error for each treatment
dw_exp6 <- kelp_data_combined %>%
  filter(experiment == "Experiment 6") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"), 
    location_label = location,
    group_id = paste(salinity, temperature, sep = ":")
  )

# Merge SNK results into plot data
dw_exp6 <- left_join(dw_exp6, snk_dw_6_df, by = c("group_id", "species"))

# Create the plot
dw_plot6 <- ggplot(dw_exp6, aes(x = treatment, y = mean_DW, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  facet_grid(. ~ species) +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Dry to Wet Weight Ratio - Experiment 6",
    x = "Treatment",
    y = "Mean D:W Ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(dw_plot6)

# Experiment 7
# Perform SNK test using existing ANOVA model for Experiment 7
snk_dw_7 <- SNK.test(dw_model_7, c("salinity", "temperature"))

# Extract and modify SNK groups
snk_dw_7_df <- snk_dw_7$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. hyperborea")  # Or adjust based on your species for this experiment

# Calculate mean DW_ratio and standard error for each treatment
dw_exp7 <- kelp_data_combined %>%
  filter(experiment == "Experiment 7") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"), 
    location_label = location,
    group_id = paste(salinity, temperature, sep = ":")
  )

# Merge SNK results into plot data
dw_exp7 <- left_join(dw_exp7, snk_dw_7_df, by = c("group_id", "species"))

# Create the plot
dw_plot7 <- ggplot(dw_exp7, aes(x = treatment, y = mean_DW, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  facet_grid(. ~ species) +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Dry to Wet Weight Ratio - Experiment 7",
    x = "Treatment",
    y = "Mean D:W Ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(dw_plot7)


## Cool!! Moving on to pigments now :)
# Experiment 4
# SNK for Chlorophyll a
snk_chla4 <- SNK.test(chla_model_4, trt = c("salinity", "temperature", "locationnumber"))
groups_chla4 <- snk_chla4$groups
groups_chla4$treatment <- rownames(groups_chla4)
groups_chla4$pigment <- "Chlorophyll_a_mg_g"

# SNK for Chlorophyll c
snk_chlc4 <- SNK.test(chlc_model_4, trt = c("salinity", "temperature", "locationnumber"))
groups_chlc4 <- snk_chlc4$groups
groups_chlc4$treatment <- rownames(groups_chlc4)
groups_chlc4$pigment <- "Chlorophyll_c_mg_g"

# SNK for Carotenoids
snk_car4 <- SNK.test(car_model_4, trt = c("salinity", "temperature", "locationnumber"))
groups_car4 <- snk_car4$groups
groups_car4$treatment <- rownames(groups_car4)
groups_car4$pigment <- "Total_carotenoid_mg_g"

# Combine all SNK results into one dataframe
groups_all4 <- bind_rows(groups_chla4, groups_chlc4, groups_car4)

# Format SNK groups to match the treatment and locationnumber structure
groups_all4 <- groups_all4 %>%
  mutate(
    raw_treatment = treatment,
    treatment = paste0(gsub(":", "PSU & ", sub(":\\d+$", "", raw_treatment)), "°C"),
    locationnumber = sub(".*:(\\d+)$", "\\1", raw_treatment),
    locationnumber = as.character(locationnumber)
  ) %>%
  select(-raw_treatment)

# Pivot the data to long format for pigment columns
pigm_long4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  pivot_longer(
    cols = c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g"),
    names_to = "pigment",
    values_to = "pigment_concentration"
  )

# Calculate mean pigment concentration and standard error
pigm_exp4 <- pigm_long4 %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species, pigment) %>%
  summarise(
    mean_pig = mean(pigment_concentration, na.rm = TRUE),
    se_pig = sd(pigment_concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = location,
    treatment_label = paste(salinity, "PSU &", temperature, "°C")
  )

# Manually join the SNK groupings with the pigment data
pigm_exp4 <- pigm_exp4 %>%
  left_join(groups_all4, by = c("treatment", "locationnumber", "pigment"))

# Adjust SNK group letters to avoid overlap
pigm_exp4 <- pigm_exp4 %>%
  mutate(groups = case_when(
    pigment == "Chlorophyll_a_mg_g" ~ "a",  # Chlorophyll a remains all "a"
    pigment == "Chlorophyll_c_mg_g" ~ "b",  # Chlorophyll c is all "b"
    pigment == "Total_carotenoid_mg_g" ~ case_when(
      groups == "a" ~ "c",  # 'a' becomes 'c' for Total carotenoids
      groups == "b" ~ "d",  # 'b' becomes 'd' for Total carotenoids
      TRUE ~ groups  # leave other groups unchanged
    )
  ))

# Create the plot
pigment_labels <- c(
  "Chlorophyll_a_mg_g" = "Chlorophyll a",
  "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
  "Total_carotenoid_mg_g" = "Total Carotenoids"
)

pigm_plot4 <- ggplot(pigm_exp4, aes(x = treatment_label, y = mean_pig, fill = location_label)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = location_colours) +
  geom_vline(xintercept = sal_break4, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = groups, y = 0), 
            position = position_dodge(0.7), vjust = -0.1, size = 4,
            color = "black", show.legend = FALSE) +
  facet_wrap(~ pigment, labeller = labeller(pigment = pigment_labels)) +
  labs(
    title = "Pigment Concentrations - Experiment 4",
    x = "Treatment",
    y = "Mean Concentration (mg/g)",
    fill = "Location"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )

# Display the plot
print(pigm_plot4)

## Alright so the SNK is a bit confusing here cause we have so many different SNK tests in Experiment 5
## So we define a function here to help with applying different letters for different species and pigments
# Define a function to apply new groupings
apply_new_groups <- function(groups, species, pigment) {
  # Mapping for L. digitata
  if (species == "L. digitata") {
    if (pigment == "Chlorophyll_a_mg_g") {
      groups <- gsub("a", "a", groups)
      groups <- gsub("b", "b", groups)
      groups <- gsub("ab", "ab", groups)  # for cases where ab exists
    } else if (pigment == "Chlorophyll_c_mg_g") {
      groups <- gsub("a", "c", groups)
    } else if (pigment == "Total_carotenoid_mg_g") {
      groups <- gsub("a", "d", groups)
      groups <- gsub("b", "e", groups)
    }
  }
  
  # Mapping for S. latissima
  if (species == "S. latissima") {
    if (pigment == "Chlorophyll_a_mg_g") {
      groups <- gsub("a", "g", groups)
      groups <- gsub("b", "h", groups)
      groups <- gsub("c", "i", groups)
      groups <- gsub("ab", "gh", groups)
      groups <- gsub("abc", "ghi", groups)
      groups <- gsub("bc", "ij", groups)
    } else if (pigment == "Chlorophyll_c_mg_g") {
      groups <- gsub("a", "j", groups)
      groups <- gsub("b", "k", groups)
    } else if (pigment == "Total_carotenoid_mg_g") {
      groups <- gsub("a", "l", groups)
    }
  }
  return(groups)
}

# Experiment 5
# Perform SNK tests using existing ANOVA models
snk_chla_dig_5 <- SNK.test(chla_dig_5, c("salinity", "temperature", "locationnumber"))
snk_chlc_dig_5 <- SNK.test(chlc_dig_5, c("salinity", "temperature", "locationnumber"))
snk_car_dig_5 <- SNK.test(car_dig_5, c("salinity", "temperature", "locationnumber"))

snk_chla_sac_5 <- SNK.test(chla_sac_5, c("salinity", "temperature", "locationnumber"))
snk_chlc_sac_5 <- SNK.test(chlc_sac_5, c("salinity", "temperature", "locationnumber"))
snk_car_sac_5 <- SNK.test(car_sac_5, c("salinity", "temperature", "locationnumber"))

# Extract and modify SNK groups for each species and pigment
snk_chla_dig_df <- snk_chla_dig_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. digitata", pigment = "Chlorophyll_a_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

snk_chlc_dig_df <- snk_chlc_dig_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. digitata", pigment = "Chlorophyll_c_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

snk_car_dig_df <- snk_car_dig_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "L. digitata", pigment = "Total_carotenoid_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

snk_chla_sac_df <- snk_chla_sac_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "S. latissima", pigment = "Chlorophyll_a_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

snk_chlc_sac_df <- snk_chlc_sac_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "S. latissima", pigment = "Chlorophyll_c_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

snk_car_sac_df <- snk_car_sac_5$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(species = "S. latissima", pigment = "Total_carotenoid_mg_g") %>%
  mutate(groups_new = mapply(apply_new_groups, groups = .$groups, species = .$species, pigment = .$pigment))

# Combine SNK groupings for all pigments and species
snk_all <- bind_rows(snk_chla_dig_df, snk_chlc_dig_df, snk_car_dig_df,
                     snk_chla_sac_df, snk_chlc_sac_df, snk_car_sac_df)

# Reshape data for Experiment 5 pigment concentrations
pigm_long5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5") %>%
  pivot_longer(cols = c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g"),
               names_to = "pigment", values_to = "pigment_concentration") %>%
  group_by(temperature, salinity, treatment, location, locationnumber, species, pigment) %>%
  summarise(
    mean_pig = mean(pigment_concentration, na.rm = TRUE),
    se_pig = sd(pigment_concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),
    location_label = location,
    group_id = paste(salinity, temperature, locationnumber, sep = ":")
  )

# Define custom labels for pigment names
pigment_labels <- c(
  "Chlorophyll_a_mg_g" = "Chlorophyll a",
  "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
  "Total_carotenoid_mg_g" = "Total Carotenoids"
)

# Merge SNK results with the reshaped pigment data
pigm_long5 <- left_join(pigm_long5, snk_all, by = c("group_id", "species", "pigment"))

# Create the plot using ggh4x for facet background colors
pigm_plot5 <- ggplot(pigm_long5, aes(x = treatment, y = mean_pig, fill = location_label)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = location_colours) +  # Location colors
  
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  
  # Correct custom labels for pigments
  facet_wrap(~ species + pigment, labeller = labeller(pigment = pigment_labels), nrow = 2) +  # Keep the two-row layout
  
  # Apply custom colors using text labels
  geom_text(aes(label = groups_new, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Pigment Concentrations - Experiment 5",
    x = "Treatment",
    y = "Mean Pigment Concentration",
    fill = "Location"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines"),
    
    # Bold facet labels
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold", color = "black"),
    
    # White background for all facet headers
    strip.background = element_rect(
      fill = "white",  # Set white background for facet headers
      color = "black"
    )
  )

# Display the plot
print(pigm_plot5)

# Experiment 6
# Run SNK tests for each pigment using the ANOVA models
snk_chla_6 <- SNK.test(chla_model_6, trt = c("salinity", "temperature"))
snk_chlc_6 <- SNK.test(chlc_model_6, trt = c("salinity", "temperature"))
snk_car_6  <- SNK.test(car_model_6, trt = c("salinity", "temperature"))

# Extract SNK group letters and add pigment labels
snk_chla_df_6 <- snk_chla_6$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Chlorophyll_a_mg_g")

snk_chlc_df_6 <- snk_chlc_6$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Chlorophyll_c_mg_g")

snk_car_df_6 <- snk_car_6$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Total_carotenoid_mg_g")

# Combine SNK results for all pigments
snk_all_6 <- bind_rows(snk_chla_df_6, snk_chlc_df_6, snk_car_df_6)

# Manually remap SNK group letters by pigment (no overlaps)
snk_all_6 <- snk_all_6 %>%
  mutate(groups = case_when(
    pigment == "Chlorophyll_c_mg_g" ~ gsub("a", "d", gsub("b", "e", groups)),
    pigment == "Total_carotenoid_mg_g" ~ gsub("[a-z]", "f", groups),
    TRUE ~ groups  # Keep original for Chlorophyll a
  ))

# Summarize pigment data for Experiment 6 and prepare for plotting
pigm_long6 <- kelp_data_combined %>%
  filter(experiment == "Experiment 6") %>%
  pivot_longer(cols = c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g"),
               names_to = "pigment", values_to = "pigment_concentration") %>%
  group_by(temperature, salinity, treatment, species, pigment) %>%
  summarise(
    mean_pig = mean(pigment_concentration, na.rm = TRUE),
    se_pig = sd(pigment_concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    group_id = paste(salinity, temperature, sep = ":"),
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C")
  )

# Merge SNK group letters into the summarized pigment data
pigm_long6 <- left_join(pigm_long6, snk_all_6, by = c("group_id", "pigment"))

# Define human-readable pigment labels
pigment_labels <- c(
  "Chlorophyll_a_mg_g" = "Chlorophyll a",
  "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
  "Total_carotenoid_mg_g" = "Total Carotenoids"
)

# Create pigment concentration plot for Experiment 6
pigm_plot6 <- ggplot(pigm_long6, aes(x = treatment, y = mean_pig, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Pigment Concentrations - Experiment 6",
    x = "Treatment",
    y = "Mean Pigment Concentration",
    fill = "Species"
  ) +
  theme_bw() +
  facet_wrap(~ pigment, labeller = labeller(pigment = pigment_labels)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )

# Display the plot
print(pigm_plot6)

# Experiment 7
# Run SNK tests for each pigment using the ANOVA models
snk_chla_7 <- SNK.test(chla_model_7, trt = c("salinity", "temperature"))
snk_chlc_7 <- SNK.test(chlc_model_7, trt = c("salinity", "temperature"))
snk_car_7  <- SNK.test(car_model_7, trt = c("salinity", "temperature"))

# Extract SNK group letters and add pigment labels
snk_chla_df_7 <- snk_chla_7$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Chlorophyll_a_mg_g")

snk_chlc_df_7 <- snk_chlc_7$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Chlorophyll_c_mg_g")

snk_car_df_7 <- snk_car_7$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_id") %>%
  mutate(pigment = "Total_carotenoid_mg_g")

# Combine SNK results for all pigments
snk_all_7 <- bind_rows(snk_chla_df_7, snk_chlc_df_7, snk_car_df_7)

# Manually remap SNK group letters by pigment (adjustments as per your instructions)
snk_all_7 <- snk_all_7 %>%
  mutate(groups = case_when(
    pigment == "Chlorophyll_c_mg_g" ~ gsub("a", "c", gsub("b", "d", groups)),
    pigment == "Total_carotenoid_mg_g" ~ gsub("a", "e", gsub("b", "f", groups)),
    TRUE ~ groups  # Keep original for Chlorophyll a
  ))

# Summarize pigment data for Experiment 7 and prepare for plotting
pigm_long7 <- kelp_data_combined %>%
  filter(experiment == "Experiment 7") %>%
  pivot_longer(cols = c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g"),
               names_to = "pigment", values_to = "pigment_concentration") %>%
  group_by(temperature, salinity, treatment, species, pigment) %>%
  summarise(
    mean_pig = mean(pigment_concentration, na.rm = TRUE),
    se_pig = sd(pigment_concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    group_id = paste(salinity, temperature, sep = ":"),
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C")
  )

# Merge SNK group letters into the summarized pigment data
pigm_long7 <- left_join(pigm_long7, snk_all_7, by = c("group_id", "pigment"))

# Define human-readable pigment labels
pigment_labels <- c(
  "Chlorophyll_a_mg_g" = "Chlorophyll a",
  "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
  "Total_carotenoid_mg_g" = "Total Carotenoids"
)

# Create pigment concentration plot for Experiment 7
pigm_plot7 <- ggplot(pigm_long7, aes(x = treatment, y = mean_pig, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig),
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = groups, y = 0), position = position_dodge(0.7), vjust = -0.1, size = 4) +
  labs(
    title = "Pigment Concentrations - Experiment 7",
    x = "Treatment",
    y = "Mean Pigment Concentration",
    fill = "Species"
  ) +
  theme_bw() +
  facet_wrap(~ pigment, labeller = labeller(pigment = pigment_labels)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.spacing.x = unit(2, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )

# Display the plot
print(pigm_plot7)