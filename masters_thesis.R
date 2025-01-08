# Load necessary libraries
library(googlesheets4)
library(emmeans)
library(car)
library(dplyr)
library(lmerTest)
library(lme4)
library(ggplot2)
library(tidyr)
library(tidyverse)

# Clean the environment
rm(list = ls())

# Load data from Google Sheets
kelp_data <- read_sheet("https://docs.google.com/spreadsheets/d/1FHVlBrssyievdxxegrkq6g7MS1HEL3AAirMl7akIjMo",
                        sheet = "kelp")

# Calculate Specific Growth Rate (SGR) for each growth period
kelp_data <- kelp_data %>%
  mutate(
    SGR_1_2 = (log(weight2_g / weight1_g) / 4) * 100,  # Day 4
    SGR_2_3 = (log(weight3_g / weight2_g) / 3) * 100,  # Day 7
    SGR_3_4 = (log(weight4_g / weight3_g) / 3) * 100,  # Day 10
    SGR_4_5 = (log(weight5_g / weight4_g) / 4) * 100   # Day 14
  )

# Line plot: salinity changes over time with the correct starting value and timepoints
salinity_plot <- kelp_data %>%
  # Pivot to long format, keeping salinity values as separate rows
  pivot_longer(
    cols = c("salinity", "sal_45", "sal_46", "sal_86", "sal_87", "sal_154", "sal_155", "sal_227", "sal_228", "sal_274", "sal_275"), 
    names_to = "timepoint", 
    values_to = "salinity_value"
  ) %>%
  # Recode timepoints to hours (e.g., sal_45 becomes 45 hours, sal_46 becomes 46 hours, etc.)
  mutate(
    timepoint = case_when(
      timepoint == "salinity" ~ 0,   # Initial salinity is the starting point (0 hours)
      timepoint == "sal_45" ~ 45,
      timepoint == "sal_46" ~ 46,
      timepoint == "sal_86" ~ 86,
      timepoint == "sal_87" ~ 87,
      timepoint == "sal_154" ~ 154,
      timepoint == "sal_155" ~ 155,
      timepoint == "sal_227" ~ 227,
      timepoint == "sal_228" ~ 228,
      timepoint == "sal_274" ~ 274,
      timepoint == "sal_275" ~ 275,
      TRUE ~ NA_real_
    ),
    timepoint = as.numeric(timepoint)
  ) %>%
  # Remove rows with NA values in the salinity_value column
  filter(!is.na(salinity_value))

# Plotting salinity changes over time
salinity_plot_final <- ggplot(salinity_plot, aes(x = timepoint, y = salinity_value, group = individual, color = factor(temperature))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Salinity change over time",
    x = "Hours",
    y = "Salinity (PSU)",
    color = "Temperature (\u00B0C)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(salinity_plot_final)

## Statistical analysis of the salinity fluctuations
# Subset the data to focus on temperature and salinity values
kelp_salinity_data <- kelp_data %>%
  pivot_longer(
    cols = c("salinity", "sal_45", "sal_46", "sal_86", "sal_87", "sal_154", "sal_155", "sal_227", "sal_228", "sal_274"), 
    names_to = "timepoint", 
    values_to = "salinity_value"
  ) %>%
  mutate(
    timepoint = case_when(
      timepoint == "salinity" ~ 0,
      timepoint == "sal_45" ~ 45,
      timepoint == "sal_46" ~ 46,
      timepoint == "sal_86" ~ 86,
      timepoint == "sal_87" ~ 87,
      timepoint == "sal_154" ~ 154,
      timepoint == "sal_155" ~ 155,
      timepoint == "sal_227" ~ 227,
      timepoint == "sal_228" ~ 228,
      timepoint == "sal_274" ~ 274,
      TRUE ~ NA_real_
    ),
    timepoint = as.numeric(timepoint)
  ) %>%
  filter(!is.na(salinity_value))

## Is it the temperature that affects the salinity fluctuations??
# ANOVA: Test if temperature treatment affects salinity fluctuations
salinity_aov <- aov(salinity_value ~ temperature + Error(individual/temperature), data = kelp_salinity_data)
summary(salinity_aov)

## It is significant! 
### So based on the information and factors we have, temperature significantly affects the salinity fluctuations

# Pairwise comparison
emmeans_result <- emmeans(salinity_aov, pairwise ~ temperature)
summary(emmeans_result)

## SGR Plot

# Prepare the data by adding a new treatment column
kelp_data_treatment <- kelp_data %>%
  mutate(
    treatment = case_when(
      temperature == 20 & salinity == 20 ~ "20 PSU, 20\u00B0C",
      temperature == 10 & salinity == 20 ~ "20 PSU, 10\u00B0C",
      temperature == 20 & salinity == 10 ~ "10 PSU, 20\u00B0C",
      temperature == 10 & salinity == 10 ~ "10 PSU, 10\u00B0C"
    )
  )

# Calculate mean SGR and standard error for each treatment
sgr_summary <- kelp_data_treatment %>%
  pivot_longer(
    cols = c("SGR_1_2", "SGR_2_3", "SGR_3_4", "SGR_4_5"),  # All growth periods
    names_to = "growth_period",
    values_to = "SGR_value"
  ) %>%
  group_by(treatment, growth_period, species) %>%
  summarise(
    mean_SGR = mean(SGR_value, na.rm = TRUE),  # Calculate mean SGR
    se_SGR = sd(SGR_value, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Plot the results as a barplot with error bars
sgr_barplot <- ggplot(sgr_summary, aes(x = treatment, y = mean_SGR, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  facet_wrap(~growth_period) +
  labs(
    title = "Mean SGR between treatments",
    x = "Treatment",
    y = "Mean Specific Growth Rate (%)",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(sgr_barplot)

# ANOVA for Specific Growth Rate (SGR)
# Test the effect of species, location, temperature, salinity, and individual on SGR

# We'll focus on one growth period as an example, such as SGR_1_2 (Day 4)
sgr_aov <- aov(SGR_4_5 ~ species * locationnumber * temperature * salinity * individual, data = kelp_data)
summary(sgr_aov)

## Temperature seems to be the most affecting factor
### But also the interaction between salinity and temperature as well as temperature and species are significant

# Prepare the data for the line plot
sgr_time_plot_data <- kelp_data %>%
  pivot_longer(
    cols = c("SGR_1_2", "SGR_2_3", "SGR_3_4", "SGR_4_5"),
    names_to = "growth_period",
    values_to = "SGR_value"
  ) %>%
  mutate(
    growth_period = factor(growth_period, 
                           levels = c("SGR_1_2", "SGR_2_3", "SGR_3_4", "SGR_4_5"),
                           labels = c("Day 4", "Day 7", "Day 10", "Day 14"))
  )

# Refine the SGR time plot
sgr_time_plot <- ggplot(sgr_time_plot_data, aes(x = growth_period, y = SGR_value, group = individual, color = species)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +  # Thinner lines with transparency
  scale_color_brewer(palette = "Set1") +  # Professional color palette
  labs(
    title = "SGR changes over time",
    x = "Time",
    y = "SGR (%)",
    color = "Species"
  ) +
  theme_minimal(base_size = 14) +  # Slightly larger font for readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),  # Lighter gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position = "top",  # Move legend to top for clarity
    plot.title = element_text(face = "bold", hjust = 0.5)  # Center bold title
  )

# Display the updated plot
print(sgr_time_plot)

## Saccharina and Laminaria's reactions seem to be different - why is that??
### Let's look into it?

## SGR changes over time depending on species

# Plot the results as a barplot with error bars
sgr_barplot <- ggplot(sgr_summary, aes(x = treatment, y = mean_SGR, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR), 
    width = 0.2, 
    position = position_dodge(0.7)  # Error bars
  ) +
  facet_grid(species ~ growth_period) +  # Facet by species (rows) and growth period (columns)
  labs(
    title = "Mean SGR between treatments",
    x = "Treatment",
    y = "Mean specific growth rate (%)",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(sgr_barplot)

## But is SGR dependent on species as well?
### I mean, it was in an interaction with temperature earlier!
#### Soooooo?


## We need to look into how it looks without the ones that have dissolved....
### So let's remove those who have negative SGR? Or is that right? Let's try it
### Or maybe we divided them into positive, negative and neutral SGR? That'd be cool right?

# Categorize SGR into Negative, Positive, and Neutral
kelp_data_sgr <- kelp_data %>%
  pivot_longer(
    cols = c("SGR_1_2", "SGR_2_3", "SGR_3_4", "SGR_4_5"),
    names_to = "growth_period",
    values_to = "SGR_value"
  ) %>%
  mutate(
    category = case_when(
      SGR_value > 0 ~ "Positive",
      SGR_value < 0 ~ "Negative",
      SGR_value == 0 ~ "Neutral"
    ),
    growth_period = factor(
      growth_period,
      levels = c("SGR_1_2", "SGR_2_3", "SGR_3_4", "SGR_4_5"),
      labels = c("Day 4", "Day 7", "Day 10", "Day 14")
    )
  )

# Count the number of samples in each category per growth period
sgr_kelp_data <- kelp_data_sgr %>%
  group_by(growth_period, category) %>%
  summarise(count = n(), .groups = "drop")

# Create a bar plot showing the count of each category per time point
sgr_category_plot <- ggplot(sgr_kelp_data, aes(x = growth_period, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Distribution of SGR over time",
    x = "Time",
    y = "Count",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(sgr_category_plot)


## Barplot of positive and negative SGR for each species and treatment
# Pivot the data to create SGR_value column
kelp_data_treatment_long <- kelp_data_treatment %>%
  pivot_longer(
    cols = starts_with("SGR_"),  # Adjust to match the column names for SGR
    names_to = "growth_period",
    values_to = "SGR_value"
  )

# Categorize SGR values and calculate summary statistics
sgr_summary <- kelp_data_treatment_long %>%
  mutate(
    SGR_category = case_when(
      SGR_value < 0 ~ "Negative",
      SGR_value == 0 ~ "Neutral",
      SGR_value > 0 ~ "Positive"
    )
  ) %>%
  group_by(treatment, species, growth_period, SGR_category) %>%
  summarise(
    mean_SGR = mean(SGR_value, na.rm = TRUE),
    se_SGR = sd(SGR_value, na.rm = TRUE) / sqrt(n()),
    count = n(),  # Count of samples in each category
    .groups = "drop"
  )

# Filter for Positive and Negative SGR categories
sgr_positive <- filter(sgr_summary, SGR_category == "Positive")
sgr_negative <- filter(sgr_summary, SGR_category == "Negative")

# Define a function to generate the plots
plot_sgr_category <- function(data, title_suffix) {
  ggplot(data, aes(x = treatment, y = mean_SGR, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_SGR - se_SGR, ymax = mean_SGR + se_SGR),
      width = 0.2,
      position = position_dodge(0.7)
    ) +
    geom_text(
      aes(label = count, y = mean_SGR + se_SGR + 1),  # Position text above error bars
      position = position_dodge(0.7),
      size = 3
    ) +
    facet_wrap(~growth_period) +
    labs(
      title = paste("Mean SGR -", title_suffix),
      x = "Treatment",
      y = "Mean SGR (%)",
      fill = "Species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Generate and display plots for Positive and Negative SGR
sgr_plot_positive <- plot_sgr_category(sgr_positive, "Positive")
sgr_plot_negative <- plot_sgr_category(sgr_negative, "Negative")

# Display the plots
print(sgr_plot_positive)
print(sgr_plot_negative)

# Check for missing or problematic rows in sgr_positive
problematic_positive <- sgr_positive %>%
  filter(
    is.na(mean_SGR) | is.na(se_SGR) | is.na(mean_SGR + se_SGR + 1)
  )

# Check for missing or problematic rows in sgr_negative
problematic_negative <- sgr_negative %>%
  filter(
    is.na(mean_SGR) | is.na(se_SGR) | is.na(mean_SGR + se_SGR + 1)
  )

# Print the problematic rows
print("Problematic rows in Positive SGR:")
print(problematic_positive)

print("Problematic rows in Negative SGR:")
print(problematic_negative)
