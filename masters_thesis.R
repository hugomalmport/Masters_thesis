# Load necessary libraries
library(googlesheets4)
library(emmeans)
library(tidyverse)

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
# Set species color
species_colours <- c("L. digitata" = "darkgreen",
                     "S. latissima" = "lightgreen",
                     "L. hyperborea" = "#33C33C")

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
                       levels = c("L. digitata", "S. latissima", "L. hyperborea"))
    )
  
  return(df)
}

# Apply the function to each dataset
kelp_data_4 <- process_data(kelp_data_4, "Experiment 4")
kelp_data_5 <- process_data(kelp_data_5, "Experiment 5")
kelp_data_6 <- process_data(kelp_data_6, "Experiment 6")
kelp_data_7 <- process_data(kelp_data_7, "Experiment 7")

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

# Plot with experiments as bottom axis and treatments above
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
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(mortality)

## Maybe adjust the requirement for mortality,
## But so far it looks nice!

# Data processing and set up for statistical analysis
# Calculate SGR_total and replace NAs with mean SGR_total of treatment and species
## This kind of replaces mortalities also
### But it's to not have any NA values while still keeping the same variance 
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
### Mainly really high chlorophyll values even though it was dissolving??
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

# Make sure 30 PSU is included as a factor and is represented
# Since there's no 30 PSU in the other experiments it could be automatically
# And also good to do the same with the locations!
kelp_data_4$salinity <- factor(kelp_data_4$salinity, levels = c(10, 20, 30))
kelp_data_4$locationnumber <- factor(kelp_data_4$locationnumber, levels = c(1, 2, 3))
# And now we need to factorize the other as well :)
kelp_data_5$salinity <- factor(kelp_data_5$salinity, levels = c(10, 20))
kelp_data_5$locationnumber <- factor(kelp_data_5$locationnumber, levels = c(1, 2))
# For experiment 6 and 7 we only had 1 location! But we still need to factorize it!
kelp_data_6$salinity <- factor(kelp_data_6$salinity, levels = c(10, 20))
kelp_data_6$locationnumber <- factor(kelp_data_6$locationnumber, levels = c(4))
kelp_data_7$salinity <- factor(kelp_data_7$salinity, levels = c(10, 20))
kelp_data_7$locationnumber <- factor(kelp_data_7$locationnumber, levels = c(4))

# Combine all datasets into one
kelp_data_combined <- bind_rows(kelp_data_4, kelp_data_5, kelp_data_6, kelp_data_7)

# Set a baseline of treatment and species
## Same one here :)
kelp_data_combined$treatment <- relevel(kelp_data_combined$treatment, ref = "10PSU & 10°C")  
kelp_data_combined$species <- relevel(kelp_data_combined$species, ref = "L. digitata")

## Let's figure out what the plots should look like based on our factors 
### And which are significant - cause no need to include things that don't matter anyway
# Statistical analysis: three-factor, four-factor, and two factor ANOVAs

# Experiment 4
sgr_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
sgr_model_4 <- aov(SGR_total ~ temperature * salinity * locationnumber, data = sgr_data_4)

summary(sgr_model_4)

## So here the interaction between temperature and salinity is very significant!
### SGR_total is heavily affected by the treatment they get!!
## So since there's a significance of the interaction - let's do a post-hoc to check where

emmeans(sgr_model_4, pairwise ~ salinity * temperature)

## The strongest differences are between the treatments that differ in both temperature and salinity
### But there's also a negative effect of higher salinity in combination with high temperature
### So they do worse in high salinity when the temperature is high
## Like they do worse in 30PSU & 20°C than in 20PSU & 20°C
### They deteriorate quicker in higher salinity when the temperature is high

# Experiment 5
sgr_data_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5")
sgr_model_5 <- aov(SGR_total ~ temperature * salinity * species * locationnumber, data = sgr_data_5)

summary(sgr_model_5)

## Firstly, the interaction between temperature and species is significant
### So species react differently strongly to temperature
## Then the temperature:salinity is significant
### So there's an effect of both salinity and temperature together
## Let's do post-hoc tests to check the interactions

emmeans(sgr_model_5, pairwise ~ salinity * temperature)

## Except for 10PSU & 10°C compared to 20PSU & 10°C, all the treatments seem to differ greatly
## A lot of significance - the smallest being between 10PSU & 10°C compared to 10PSU & 20°C
### Although the difference between 10PSU & 10°C and 10PSU & 20°C isn't too significant

## Temperature and species
emmeans(sgr_model_5, pairwise ~ temperature * species)

## The species reactions differ greatly, except for when both of them are treated to 20°C
### Which both have a very negative reaction to!

# Experiment 6
sgr_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
sgr_model_6 <- aov(SGR_total ~ temperature * salinity, data = sgr_data_6)

summary(sgr_model_6)

## So here the interaction between temperature and salinity is also significant
## Let's do post-hoc!

emmeans(sgr_model_6, pairwise ~ salinity * temperature)

## A lot of significance
## One interesting thing is 10PSU & 10°C seems to be doing better than the control?
### But that's not really significant - but close to, so if you increased the sample size it might be?

## Let's check Experiment 7, cause there we used juveniles instead!

# Experiment 7
sgr_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
sgr_model_7 <- aov(SGR_total ~ temperature * salinity, data = sgr_data_7)

summary(sgr_model_7)

## So this is different
## Here the interaction is not significant anymore!!
### Only temperature and salinity as main effects are significant
## And salinity isn't as strong as temperature 
## It would be cool to see if adults and juveniles are different
## I can't really analyse this with an ANOVA though, since both the temperatures and lifestages are different
### But it could say something about the usage of disks - maybe isn't the best when looking at salinities
## But either way, when it comes to SGR - temperature and salinity is always a significant factor!!
### Meaning they heavily affected the SGR in my experiments

emmeans(sgr_model_7, pairwise ~ salinity * temperature)

## The interaction isn't significant but there is a huge significant between 20PSU & 10°C compared to 10PSU & 17°C
### And also a significant difference between 20PSU & 10°C compared to 20PSU & 17°C

## Dry to wet weight ratio - statistical analyses!!

# Experiment 4
dw_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
dw_model_4 <- aov(DW_ratio ~ temperature * salinity * locationnumber, data = dw_data_4)

summary(dw_model_4)

## Nothing is significant?
### Which is weird - I don't know what to think about that
## Let's do the others as well

# Experiment 5
dw_data_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5")
dw_model_5 <- aov(DW_ratio ~ temperature * salinity * locationnumber * species, data = dw_data_5)

summary(dw_model_5)

## Hmm, here temperature, species and the interaction between salinity:location is significant
## Interesting...
### So the interaction could mean that there's a difference response to salinity between the two locations
### And then the main factors could be that in the higher temperature they accumulate more water - cause they do worse?
## And S. latissima probably accumulates more water than L. digitata
## But the interaction is interesting
## Let's look closer into that...

# Run pairwise comparisons for the significant interaction salinity:location
emmeans(dw_model_5, pairwise ~ salinity * locationnumber)

## So the interaction of salinity and location is almost significant within location 1 for salinities
## But the interaction becomes the most significant when looking at salinity 20 between the locations

## I also wanna look into how the treatments differ
### Why isn't the effect of salinity and temperature significant

emmeans(dw_model_5, pairwise ~ salinity * temperature)

## So the most significant difference - between 20PSU & 10°C (the control), and 10PSU & 20°C
### This makes sense since that treatment is the most different from the control
## Then there's another significant difference between the control and 20°C
### Which is also seen in the significance of temperature as a main factor

## Post-hoc test on temperature:location

emmeans(dw_model_5, pairwise ~ temperature * locationnumber)

## The differences here is between the degrees on Ängklåvbukten
## And the different degrees between the locations as well - 10°C & Ängklåvbukten compared to 20°C & Lökholmen

## Let's check the other two experiments as well - with only two factors

# Experiment 6
dw_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
dw_model_6 <- aov(DW_ratio ~ temperature * salinity, data = dw_data_6)

summary(dw_model_6)

## So here the DW_ratio depends on temperature and salinity as main factors
### But the interaction is close to being significant
## Let's do a post-hoc anyways!
emmeans(dw_model_6, pairwise ~ salinity * temperature)

## So the reason the interaction isn't significant is because there was no difference between 10PSU & 10°C and 20PSU & 20°C
## But the main factors are very significant

# Experiment 7
dw_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
dw_model_7 <- aov(DW_ratio ~ temperature * salinity, data = dw_data_7)

summary(dw_model_7)

## And here only temperature affected the dry wet weight ratio, hmmm

## Chlorophyll and carotenoid statistics!!!
## Let's do very pigment for itself - not grouping them

# Experiment 4
# Chlorophyll a
chla_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
chla_model_4 <- aov(Chlorophyll_a_mg_g ~ temperature * salinity * locationnumber, data = chla_data_4)

summary(chla_model_4)

## Temperature is significant by itself, salinity is almost

# Chlorophyll c
chlc_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
chlc_model_4 <- aov(Chlorophyll_c_mg_g ~ temperature * salinity * locationnumber, data = chlc_data_4)

summary(chlc_model_4)

## No significance

# Carotenoids
car_data_4 <- kelp_data_combined %>% filter(experiment == "Experiment 4")
car_model_4 <- aov(Total_carotenoid_mg_g ~ temperature * salinity * locationnumber, data = car_data_4)

summary(car_model_4)

## No significance

# Experiment 5
# Chlorophyll a
chla_data_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5")
chla_model_5 <- aov(Chlorophyll_a_mg_g ~ temperature * salinity * locationnumber * species, data = chla_data_5)

summary(chla_model_5)

## Firstly, the interaction temperature:location:species is a little significant
## Secondly, temperature:species is a little more significant
## And then salinity:location is almost significant
## Species by itself is very significant - S. latissima has a lot more chlorophyll a
## location is almost significant
## Salinity is significant
## And finally temperature is very significant
### Post-hoc, checking the interaction of temperature and species between locations!

## Visualize the data with emmip:
emmip_chla_5 <- emmip(chla_model_5, temperature ~ salinity | locationnumber | species, CIs = TRUE) +
  ggtitle("Chlorophyll a - Experiment 5")
emmip_chla_5

## S. latissima shows greater differences than L. digitata
## Let's check the numbers

emm_chla_5 <- emmeans(chla_model_5, pairwise ~ locationnumber * temperature * species)
pairs(emm_chla_5)

## A lot of significance 
### But one notable thing is that the comparisons between L. digitata is never significant

## Salinity:location is also significant
emmeans(chla_model_5, pairwise ~ salinity * locationnumber)

## There is a difference between salinities in location 1
### But not when looking at 10PSU between the locations
### And location 2 doesn't seem to have a difference between salinities
## When comparing 20PSU between the locations it is almost significant

# Chlorophyll c
chlc_data_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5")
chlc_model_5 <- aov(Chlorophyll_c_mg_g ~ temperature * salinity * locationnumber * species, data = chlc_data_5)

summary(chlc_model_5)

## The four-way interaction is almost significant, and so is salinity:location:species
## temperature:species is significant
## Species is very significant by itself - meaning there's a big difference in pigment concentration here as well
### S. latissima has more chlorophyll c as well
## Salinity is a little significant
## Temperature is very significant
## Post-hoc!

## Post-hoc on the four factors - but since species is very significant we don't need to look more into that
emm_chlc_5 <- emmeans(chlc_model_5, ~ salinity * temperature * locationnumber | species)

## Since we have many factors we use tukey
pairs(emm_chlc_5, adjust = "tukey")

## So there are no significant comparisons for L. digitata
### So the difference between S. latissima and L. digitata could be the deciding factor
## But there are some significant differences for S. latissima

## Visualizing the differences!!
emmip_chlc_5 <- emmip(chlc_model_5, temperature ~ salinity | locationnumber | species, CIs = TRUE) +
  ggtitle("Chlorophyll c - Experiment 5")
emmip_chlc_5

## We can see that L. digitata doesn't vary a lot between salinities and temperature
## They have a pretty linear response - especially in location 1
### But in location 2 they have more pigment in 20PSU and 10°C than in the other, although CI are overlapping
### This most likely means that they aren't different at all
## S. latissima has a very different response though
### They do so much better in 10°C but they have more pigments in 20PSU for location 1
## And for location 2 the response is kind of linear - although they have more pigments in 10°C
### But it seems to be a little more in 10PSU there for some reason...

# Carotenoids
car_data_5 <- kelp_data_combined %>% filter(experiment == "Experiment 5")
car_model_5 <- aov(Total_carotenoid_mg_g ~ temperature * salinity * locationnumber * species, data = car_data_5)

summary(car_model_5)

## Four-way interaction is significant
## All the three-way interactions including species is either significant or very close to
### So species likely plays a big role here as well
## Temperature:location is significant
## Species as a main factor is significant
## And salinity as well - but not temperature this time?
### So carotenoids might not be affected by temperature the same way chlorophyll is?

## Post-hoc - but here species isn't as significant so we could start looking at the emmip

emmip_car_5 <- emmip(car_model_5, temperature ~ salinity | locationnumber | species, CIs = TRUE) +
  ggtitle("Carotenoid content - Experiment 5")
emmip_car_5

## So there are two instances where the CIs don't overlap - and that's for S. latissima
### Something that's weird is that S. latissima has more carotenoids in 20°C than in 10°C for location 1
## Let's do the same emmeans - L. digitata doesn't seem to differ a lot

emm_car_5 <- emmeans(car_model_5, pairwise ~ salinity * temperature * locationnumber | species)
pairs(emm_car_5)

## Most of the significant comparisons come from S. latissima, but L. digitata has one comparison
### It's between salinities within location and temperature - which is almost significant
## Otherwise S. latissima has the most significant ones
### But they're not strongly significant, and some are only close to being
## But most of them are between locations

# Experiment 6
# Chlorophyll a
chla_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
chla_model_6 <- aov(Chlorophyll_a_mg_g ~ temperature * salinity, data = chla_data_6)

summary(chla_model_6)

## Salinity is very significant, while temperature is a little significant
## Not the interaction though

# Chlorophyll c
chlc_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
chlc_model_6 <- aov(Chlorophyll_c_mg_g ~ temperature * salinity, data = chlc_data_6)

summary(chlc_model_6)

## Here only salinity is significant

# Carotenoids
car_data_6 <- kelp_data_combined %>% filter(experiment == "Experiment 6")
car_model_6 <- aov(Total_carotenoid_mg_g ~ temperature * salinity, data = car_data_6)

summary(car_model_6)

## Here nothing is significant?

# Experiment 7
# Chlorophyll a
chla_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
chla_model_7 <- aov(Chlorophyll_a_mg_g ~ temperature * salinity, data = chla_data_7)

summary(chla_model_7)

## Salinity is significant
## The interaction is almost significant
## Let's check that with post-hoc test, tukey

emm_chla_7 <- emmeans(chla_model_7, pairwise ~ salinity * temperature)
pairs(emm_chla_7)

emmip_chla_7 <- emmip(chla_model_7, temperature ~ salinity, CIs = TRUE) +
  ggtitle ("Chlorophyll a - Experiment 7")
emmip_chla_7

## So the differences appear when comparing salinities - 10PSU & 10°C with 20PSU & 10°C and 10PSU & 17°C with 20PSU & 17°C
### We can see that they have more chlorophyll a in higher salinities (higher means - negative estimates)
## And there is also a significance when comparing the interactions - 10PSU & 10°C with 20PSU & 17°C and 20PSU & 10°C with 10PSU & 17°C
## This is probably why the interaction is close to being significant. But it only close since the temperature comparisons aren't significant I guess?
### But we can see that they have more chlorophyll a in high salinity and high temperature than in low salinity low temperature
## So salinity is most likely the deciding factor when looking at chlorophyll a

# Chlorophyll c
chlc_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
chlc_model_7 <- aov(Chlorophyll_c_mg_g ~ temperature * salinity, data = chlc_data_7)

summary(chlc_model_7)

## Salinity is significant

# Carotenoids
car_data_7 <- kelp_data_combined %>% filter(experiment == "Experiment 7")
car_model_7 <- aov(Total_carotenoid_mg_g ~ temperature * salinity, data = car_data_7)

summary(car_model_7)

## Salinity is significant
### We should visualize these results also

## SGR plots!

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
    experiment == "Experiment 5" ~ "Experiment 5\n(14 days)",
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
  facet_grid(. ~ experiment_label, scales = "free_x") +  # Facet only by experiment
  labs(
    title = "SGR",
    x = "Treatment",
    y = expression("SGR (% day"^-1*")"),
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Display the plot
print(sgr_barplot)

## Location is not included here, but I want to look into that!
### So let's do two plots for experiment 4 and 5, with locations included
## Let's start with experiment 4!
# Calculate mean SGR_total and standard error for each treatment
sgr_exp4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  group_by(temperature, salinity, treatment, location, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = paste0(location)
  )

# Ensure the x-axis is treated as a factor (to maintain order)
sgr_exp4$sal_temperature <- factor(sgr_exp4$sal_temp, levels = unique(sgr_exp4$sal_temp))

# Define positions for vertical lines (between salinities)
sal_break41 <- c(2.5, 4.5)  # Adjust based on the number of salinity levels

sgr_plot4 <- ggplot(sgr_exp4, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  # Add vertical lines between different salinities
  geom_vline(xintercept = sal_break41, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
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

## I kinda like this look, but changing temperature and location also works

# Let's move onto experiment 5!
sgr_exp5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5") %>%
  group_by(temperature, salinity, treatment, location, species) %>%
  summarise(
    mean_SGR = mean(SGR_total, na.rm = TRUE),
    se_SGR = sd(SGR_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = paste0(location)
  )
# Ensure the x-axis is treated as a factor (to maintain order)
sgr_exp5$sal_temperature <- factor(sgr_exp5$sal_temp, levels = unique(sgr_exp5$sal_temp))

# Define positions for vertical lines (between salinities)
sal_break51 <- c(2.5)

sgr_plot5 <- ggplot(sgr_exp5, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break51, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  labs(
    title = "SGR - Experiment 5",
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
print(sgr_plot5)

## This looks so much nicer - cause this one is easier to read :)

## Line plot of how the weight changes over time!!

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
    legend.position = "right"
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

## Dry:wet weight ratio plots!!

# Calculate mean W:D-ratio and standard error for each treatment, species, and experiment
dw_summary <- kelp_data_combined %>%
  group_by(treatment, species, experiment) %>%
  summarise(
    mean_DW_ratio = mean(`DW_ratio`, na.rm = TRUE),  # Calculate mean wet:dry ratio
    se_DW_ratio = sd(`DW_ratio`, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the plot for wet:dry weight ratio (D:W-ratio)
dw_ratio_barplot <- ggplot(dw_summary, aes(x = treatment, y = mean_DW_ratio, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) + 
  geom_errorbar(
    aes(ymin = mean_DW_ratio - se_DW_ratio, ymax = mean_DW_ratio + se_DW_ratio), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  scale_fill_manual(values = species_colours) +
  facet_grid(. ~ experiment, scales = "free_x") +
  labs(
    title = "Dry to wet weight-ratio",
    x = "Treatment",
    y = "Mean D:W ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Display the plot
print(dw_ratio_barplot)

## This looks weird
### It doesn't look like the treatments are different - even though we might expect the lower salinity to accumulate more water
### But S. latissima is different from digitata in Experiment 5 - clearly
### And 30PSU & 10°C seems to be pretty different from the other treatments in Experiment 4 though

## Let's check at experiment 4 and 5 separately!

# Experiment 4
dw_exp4 <- kelp_data_combined %>%
  filter(experiment == "Experiment 4") %>%
  group_by(temperature, salinity, treatment, location, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = paste0(location)
  )
# Ensure the x-axis is treated as a factor (to maintain order)
dw_exp4$sal_temperature <- factor(dw_exp4$sal_temp, levels = unique(dw_exp4$sal_temp))

# Define positions for vertical lines (between salinities)
sal_break42 <- c(2.5, 4.5)

dw_plot4 <- ggplot(dw_exp4, aes(x = treatment, y = mean_DW, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break42, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  labs(
    title = "Dry to wet weight-ratio - Experiment 4",
    x = "Treatment",
    y = "Mean D:W ratio",
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
dw_exp5 <- kelp_data_combined %>%
  filter(experiment == "Experiment 5") %>%
  group_by(temperature, salinity, treatment, location, species) %>%
  summarise(
    mean_DW = mean(DW_ratio, na.rm = TRUE),
    se_DW = sd(DW_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = paste0(location)
  )
# Ensure the x-axis is treated as a factor (to maintain order)
dw_exp5$sal_temperature <- factor(dw_exp5$sal_temp, levels = unique(dw_exp5$sal_temp))

# Define positions for vertical lines (between salinities)
sal_break52 <- c(2.5)

dw_plot5 <- ggplot(dw_exp5, aes(x = treatment, y = mean_DW, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_DW - se_DW, ymax = mean_DW + se_DW), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_fill_manual(values = species_colours) +
  geom_vline(xintercept = sal_break52, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  labs(
    title = "Dry to wet weight-ratio - Experiment 5",
    x = "Treatment",
    y = "Mean D:W ratio",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")  # Add spacing between temperature facets
  )

# Display the plot
print(dw_plot5)

## Cool!! Moving on to pigments now :)
# Reshape data: Convert pigments into long format
pigm_long <- kelp_data_combined %>%
  pivot_longer(cols = c("Chlorophyll_a_mg_g", "Chlorophyll_c_mg_g", "Total_carotenoid_mg_g"),
               names_to = "pigment",
               values_to = "value")

# Summarize data: Calculate mean and standard error for each pigment type
pigm_summary <- pigm_long %>%
  group_by(treatment, species, experiment, pigment) %>%
  summarise(
    mean_pig = mean(value, na.rm = TRUE),
    se_pig = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define custom colors for pigments
pigment_colors <- c("Chlorophyll_a_mg_g" = "aquamarine4",
                    "Chlorophyll_c_mg_g" = "chartreuse",
                    "Total_carotenoid_mg_g" = "yellow")

# Function to create one plot per experiment
plot_experiment <- function(exp_num) {
  ggplot(filter(pigm_summary, experiment == exp_num), aes(x = treatment, y = mean_pig, fill = pigment, color = species)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6, linewidth = 1) +
    geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig),
                  width = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values = pigment_colors, 
                      labels = c("Chlorophyll_a_mg_g" = "Chlorophyll a",
                                 "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
                                 "Total_carotenoid_mg_g" = "Total Carotenoids")) +
    scale_color_manual(values = species_colours) +
    labs(
      title = paste("Pigment Concentrations - ", exp_num),
      x = "Treatment",
      y = "Mean Concentration (mg/g)",
      fill = "Pigment",
      color = "Species"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

# Generate plots for each experiment
experiments <- unique(pigm_summary$experiment)
pig_plots <- lapply(experiments, plot_experiment)

# Display all plots
pig_plots

## Also, we should do a location based plot for pigments as well!
### That'd be cool :)

# Experiment 4
# Filter for Experiment 4 & summarize data
pigm_summary_exp4 <- pigm_long %>%
  filter(experiment == "Experiment 4") %>%
  group_by(temperature, salinity, treatment, location, species, pigment) %>%
  summarise(
    mean_pig = mean(value, na.rm = TRUE),
    se_pig = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = as.character(location)
  )

# Ensure correct x-axis order
pigm_summary_exp4$sal_temperature <- factor(pigm_summary_exp4$sal_temp, levels = unique(pigm_summary_exp4$sal_temp))

# Define vertical line positions
sal_break43 <- c(2.5, 4.5)

# Create the plot
pigm_exp4 <- ggplot(pigm_summary_exp4, aes(x = treatment, y = mean_pig, fill = pigment, color = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_color_manual(values = species_colours) + 
  scale_fill_manual(values = pigment_colors, 
                    labels = c("Chlorophyll_a_mg_g" = "Chlorophyll a",
                               "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
                               "Total_carotenoid_mg_g" = "Total Carotenoids")) +
  geom_vline(xintercept = sal_break43, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  labs(
    title = "Pigment Concentrations - Experiment 4",
    x = "Treatment",
    y = "Mean Concentration (mg/g)",
    fill = "Pigment",
    color = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(pigm_exp4)

# Experiment 5
# Filter for Experiment 5 & summarize data
pigm_summary_exp5 <- pigm_long %>%
  filter(experiment == "Experiment 5") %>%
  group_by(temperature, salinity, location, treatment, species, pigment) %>%
  summarise(
    mean_pig = mean(value, na.rm = TRUE),
    se_pig = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    sal_temp = paste0(salinity, "PSU - ", temperature, "°C"),  # Add "PSU" to salinity
    location_label = as.character(location)
  )

# Ensure correct x-axis order
pigm_summary_exp5$sal_temperature <- factor(pigm_summary_exp5$sal_temp, levels = unique(pigm_summary_exp5$sal_temp))

# Define vertical line positions
sal_break53 <- c(2.5)

# Create the plot
pigm_exp5 <- ggplot(pigm_summary_exp5, aes(x = treatment, y = mean_pig, fill = pigment, color = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pig - se_pig, ymax = mean_pig + se_pig), 
                width = 0.2, position = position_dodge(0.7)) +
  scale_color_manual(values = species_colours) + 
  scale_fill_manual(values = pigment_colors, 
                    labels = c("Chlorophyll_a_mg_g" = "Chlorophyll a",
                               "Chlorophyll_c_mg_g" = "Chlorophyll c1 + c2",
                               "Total_carotenoid_mg_g" = "Total Carotenoids")) +
  geom_vline(xintercept = sal_break53, linetype = "dashed", color = "gray50") +
  
  facet_grid(. ~ location_label) +
  labs(
    title = "Pigment Concentrations - Experiment 5",
    x = "Treatment",
    y = "Mean Concentration (mg/g)",
    fill = "Pigment",
    color = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    panel.spacing.x = unit(2, "lines")
  )

# Display the plot
print(pigm_exp5)

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
    timepoint = gsub("sal_", "", timepoint),  # Remove 'sal_' prefix
    timepoint = as.numeric(timepoint),  # Attempt to convert to numeric
    salinity_value = as.numeric(salinity_value)  # Ensure values are numeric
  ) %>%
  filter(!is.na(timepoint) & !is.na(salinity_value))

# Plot the salinity fluctuations!!
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
  theme_bw()

print(sal_plot)
## We can clearly see that temperature does not have a big effect on the salinity increases
### since there is clear overlap of the two temperatures (blue and red)
### And it looks pretty much the same for all the experiments (except 4 ig :) )

## Saving the plots!!
# Define output directory
output_dir <- "C:/Users/hugom/Desktop/master thesis/Programming/Plots"

# Common width, height and dpi settings:
w <- 10      # width in inches
h <- 7      # height in inches
res <- 500  # dpi

# Save the plots
ggsave(file.path(output_dir, "emmip_chla_5.png"), plot = emmip_chla_5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "emmip_chlc_5.png"), plot = emmip_chlc_5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "emmip_car_5.png"), plot = emmip_car_5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "emmip_chla_7.png"), plot = emmip_chla_7, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sgr_barplot.png"), plot = sgr_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sgr_plot4.png"), plot = sgr_plot4, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sgr_plot5.png"), plot = sgr_plot5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "weight_plot.png"), plot = weight_plot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "dw_ratio_barplot.png"), plot = dw_ratio_barplot, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "dw_plot4.png"), plot = dw_plot4, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "dw_plot5.png"), plot = dw_plot5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "pigm_exp4.png"), plot = pigm_exp4, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "pigm_exp5.png"), plot = pigm_exp5, width = w, height = h, dpi = res)
ggsave(file.path(output_dir, "sal_plot.png"), plot = sal_plot, width = w, height = h, dpi = res)

## Since there are many plots under "pig_plots"
# Loop through the list of plots and save each one
for (i in 1:length(pig_plots)) {
  ggsave(file.path(output_dir, paste0("pig_plots_exp", experiments[i], ".png")), 
         plot = pig_plots[[i]], width = w, height = h, dpi = res)
}

