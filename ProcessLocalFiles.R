# Load required libraries
library(tidyverse)
library(nlme)        # For mixed models with autocorrelation (AR1) structures
library(broom.mixed) # For tidying model outputs
library(RColorBrewer)

# ---------------------------------------------------------
# 1. Data Ingestion and Preparation
# ---------------------------------------------------------

# Locate all CSV files in the Data subfolder
file_list <- list.files(path = "Data", pattern = "*.csv", full.names = TRUE)


# Target speed in km/h
target_speed_kmh <- 80 
# 80 km/h converted to m/s
target_speed_ms <- 80 * (1000 / 3600) 

read_sim_data <- function(file_path) {
   df <- read_csv(file_path, show_col_types = FALSE)
   
   filename <- basename(file_path)
   part_id <- str_extract(filename, "^[0-9]+")
   
   df %>%
      mutate(
         # Force to numeric and convert m/s to km/h
         Speed_kmh = as.numeric(Speed) * 3.6,
         
         ParticipantID = part_id,
         FileID = filename,
         IsStraight = if_else(RoadDirection == "Straight", 1, 0),
         RoadAttribute = as.factor(RoadAttribute),
         TimeIndex = row_number(),
         
         # Calculate deviation from 80 km/h
         Speed_Deviation_kmh = Speed_kmh - target_speed_kmh
      )
}

# Bind all files into one master dataframe
df_all <- map_dfr(file_list, read_sim_data)

# OPTIONAL BUT RECOMMENDED: Downsample data to improve AR(1) processing time
# Keeping every 10th row effectively reduces 50Hz to 5Hz
df_analysis <- df_all %>%
   filter(TimeIndex %% 10 == 0) %>%
   group_by(FileID) %>%
   mutate(
      # Re-index time after downsampling
      TimeIndex = row_number(),
      
      # Relevel RoadSlope inside the mutate function
      RoadSlope = relevel(as.factor(RoadSlope), ref = "Flat")
   ) %>%
   ungroup()

# ---------------------------------------------------------
# 2. Extracting Individual Effect Sizes
# ---------------------------------------------------------

# We map over each participant to fit individual Generalized Least Squares (GLS) 
# models containing an AR(1) structure nested within each file.
individual_effects <- df_analysis %>%
   group_by(ParticipantID) %>%
   nest() %>%
   mutate(
      model = map(data, ~ gls(
         Speed_Deviation_kmh ~ RoadAttribute * PropOfTotalHeight + RoadSlope + IsStraight,
         data = .x,
         # Autocorrelation handled per file for the participant
         correlation = corAR1(form = ~ TimeIndex | FileID),
         na.action = na.omit
      )),
      coefs = map(model, tidy) # Extract coefficients cleanly
   ) %>%
   unnest(coefs) %>%
   select(ParticipantID, term, estimate, std.error, p.value)

# View individual effect sizes
print(individual_effects)

effects_summary <- individual_effects %>%
   group_by(term) %>%
   summarise(
      Min = min(estimate),
      Median = median(estimate),
      Mean = mean(estimate),
      Max = max(estimate),
      SD = sd(estimate),
      .groups = "drop"
   ) %>%
   arrange(desc(Mean))

print(effects_summary)
library(ggplot2)

ggplot(individual_effects, aes(x = estimate, y = reorder(term, estimate, FUN = median))) +
   geom_boxplot(fill = "gray90", outlier.shape = NA, width = 0.6) +
   geom_jitter(aes(color = ParticipantID), width = 0, height = 0.15, size = 2, alpha = 0.9) +
   geom_vline(xintercept = 0, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
   # Force a high-contrast palette supporting up to 12 categories
   scale_color_brewer(palette = "Paired") + 
   theme_minimal() +
   labs(
      title = "Distribution of Individual Effect Sizes",
      x = "Estimated Effect on Speed Deviation (km/h)",
      y = NULL
   ) +
   theme(
      legend.position = "right", # Turn legend on to verify 10 distinct IDs
      axis.text.y = element_text(face = "bold")
   )

# ---------------------------------------------------------
# 3. Global Mixed-Effects Model
# ---------------------------------------------------------

# For population-level inference across all participants, we use a Linear Mixed Model (LMM).
# We include random intercepts and PropHeigh slope for participants.

# 1. Aggressive Downsampling (1Hz)
df_global <- df_all %>%
   # Keep every 50th row (reduces 50Hz to 1Hz)
   filter(TimeIndex %% 50 == 0) %>% 
   group_by(FileID) %>%
   mutate(
      # Sequential index for AR(1) at 1Hz
      TimeIndex_1Hz = row_number(),
      # Ensure Flat is the reference
      RoadSlope = relevel(as.factor(RoadSlope), ref = "Flat"),
      # Ensure Normal is the reference
      RoadAttribute = relevel(as.factor(RoadAttribute), ref = "Normal")
   ) %>%
   ungroup()

# 2. Fit the Updated Global LMM
global_model <- lme(
   fixed = Speed_Deviation_kmh ~ RoadAttribute * PropOfTotalHeight + RoadSlope + IsStraight,
   
   # Random slope + intercept for Participant; Random intercept only for File
   random = list(ParticipantID = ~ PropOfTotalHeight, FileID = ~ 1), 
   
   correlation = corAR1(form = ~ TimeIndex_1Hz | ParticipantID / FileID),
   data = df_global,
   na.action = na.omit,
   control = lmeControl(opt = "optim", maxIter = 500, msMaxIter = 500) 
)

summary(global_model)
