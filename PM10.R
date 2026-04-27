# PM10 mirror analysis (same structure as your PM2.5 workflow), with speed tweaks
# Save as a NEW script and run top-to-bottom.

# =========================
# 0) SETTINGS
# =========================
fast_mode <- FALSE          # TRUE = faster dev run; FALSE = final run
run_imputation <- TRUE     # set FALSE to skip mice while testing

pm10_dir <- "/Users/aakashparekh/Desktop/Epi internship/PM10 Speciation values"
scd_file <- "/Users/aakashparekh/Desktop/Epi internship/R_cognitive_decline_by_state_2015_to_2023.csv"
pm10_cache <- "/Users/aakashparekh/Desktop/Epi internship/pm10_all_states_years.csv"

# fast-mode knobs
m_imp <- if (fast_mode) 3 else 10
maxit_imp <- if (fast_mode) 5 else 10

# =========================
# 1) LIBRARIES
# =========================
library(tidyverse)
library(lubridate)
library(janitor)
library(zoo)
library(mice)
library(lmerTest)
library(gtsummary)
library(broom.mixed)

# =========================
# 2) BUILD OR LOAD PM10 STATE-YEAR DATA
# =========================
if (file.exists(pm10_cache)) {
  message("Loading cached PM10 state-year file...")
  pm10_all_states_years <- read_csv(pm10_cache, show_col_types = FALSE)
} else {
  message("Cache not found. Building PM10 state-year file from hourly CSVs...")
  pm10_files <- list.files(
    pm10_dir,
    pattern = "^hourly_PM10SPEC_[0-9]{4}\\.csv$",
    full.names = TRUE
  )
  if (length(pm10_files) == 0) stop("No PM10 files found.")
  
  # Faster: summarize each file immediately (memory-safe)
  pm10_state_year_list <- lapply(pm10_files, function(f) {
    message("Reading: ", basename(f))
    df <- read_csv(f, show_col_types = FALSE) %>% clean_names()
    
    out <- df %>%
      transmute(
        StateName = state_name,
        DateLocalRaw = as.character(date_local),
        SampleMeasurement = as.numeric(sample_measurement)
      ) %>%
      mutate(
        DateLocal = suppressWarnings(as.Date(DateLocalRaw, format = "%m/%d/%y")),
        DateLocal = coalesce(DateLocal, suppressWarnings(as.Date(DateLocalRaw, format = "%m/%d/%Y"))),
        DateLocal = coalesce(DateLocal, suppressWarnings(as.Date(DateLocalRaw, format = "%Y-%m-%d"))),
        Year = as.integer(format(DateLocal, "%Y"))
      ) %>%
      filter(Year >= 2010, Year <= 2023) %>%
      group_by(StateName, Year) %>%
      summarise(Mean_PM10 = mean(SampleMeasurement, na.rm = TRUE), .groups = "drop")
    
    rm(df); gc()
    out
  })
  
  pm10_all_states_years <- bind_rows(pm10_state_year_list) %>%
    group_by(StateName, Year) %>%
    summarise(Mean_PM10 = mean(Mean_PM10, na.rm = TRUE), .groups = "drop") %>%
    arrange(StateName, Year) %>%
    group_by(StateName) %>%
    mutate(
      PM10_lag1 = lag(Mean_PM10, 1),
      PM10_lag2 = lag(Mean_PM10, 2),
      PM10_lag3 = lag(Mean_PM10, 3)
    ) %>%
    ungroup()
  
  write_csv(pm10_all_states_years, pm10_cache)
  message("Saved cache: ", pm10_cache)
}

# =========================
# 3) SCD PREP (same as PM2.5 workflow)
# =========================
Total_Cog_Decline <- read_csv(scd_file, show_col_types = FALSE)

Total_Cog_Decline_long <- Total_Cog_Decline %>%
  clean_names() %>%
  rename(StateName = state_name) %>%
  pivot_longer(
    cols = starts_with("x"),
    names_to = "Year",
    names_prefix = "x",
    values_to = "MeanCognitiveDecline"
  ) %>%
  mutate(
    Year = as.integer(Year),
    StateName = str_squish(str_to_title(StateName))
  )

pm10_all_states_years <- pm10_all_states_years %>%
  mutate(StateName = str_squish(str_to_title(StateName)))

state_data_counts <- Total_Cog_Decline_long %>%
  group_by(StateName) %>%
  summarise(years_with_data = sum(!is.na(MeanCognitiveDecline)), .groups = "drop")

states_with_enough_data <- state_data_counts %>%
  filter(years_with_data > 2) %>%
  pull(StateName)

Total_Cog_Decline_long_filtered <- Total_Cog_Decline_long %>%
  filter(StateName %in% states_with_enough_data)

Total_Cog_Decline_long_interp <- Total_Cog_Decline_long %>%
  group_by(StateName) %>%
  arrange(Year) %>%
  mutate(InterpolatedLine = na.approx(MeanCognitiveDecline, Year, na.rm = FALSE)) %>%
  ungroup()

# =========================
# 4) MERGE PM10 (same structure as PM2.5 workflow)
# =========================
merged_long_data_pm10 <- Total_Cog_Decline_long_interp %>%
  left_join(pm10_all_states_years, by = c("StateName", "Year"))

filtered_merged_long_data_pm10 <- Total_Cog_Decline_long_filtered %>%
  left_join(pm10_all_states_years, by = c("StateName", "Year"))

# =========================
# 5) RAW MODELS (same style as your PM2.5 chunks)
# =========================
mod1_full_PM10 <- merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag1) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

mod2_full_PM10 <- merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag2) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

mod3_full_PM10 <- merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag3) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

summary(mod1_full_PM10)
summary(mod2_full_PM10)
summary(mod3_full_PM10)

mod1_filt_PM10 <- filtered_merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag1) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

mod2_filt_PM10 <- filtered_merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag2) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

mod3_filt_PM10 <- filtered_merged_long_data_pm10 %>%
  rename(PM10_lag = PM10_lag3) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + PM10_lag + (1 | StateName), data = .)

summary(mod1_filt_PM10)
summary(mod2_filt_PM10)
summary(mod3_filt_PM10)

# =========================
# 6) IMPUTATION PIPELINE (same structure; fast-mode uses fewer iterations)
# =========================
if (run_imputation) {
  wide_merged_data_pm10 <- merged_long_data_pm10 %>%
    select(-InterpolatedLine) %>%
    pivot_wider(
      values_from = c("MeanCognitiveDecline", "Mean_PM10", "PM10_lag1", "PM10_lag2", "PM10_lag3"),
      names_from = "Year"
    )
  
  ini <- mice(wide_merged_data_pm10, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  meth[] <- ""
  meth[grep("^MeanCognitiveDecline_", names(wide_merged_data_pm10))] <- "pmm"
  
  imputed_data_pm10 <- mice(
    wide_merged_data_pm10,
    m = m_imp,
    maxit = maxit_imp,
    method = meth,
    predictorMatrix = pred,
    seed = 8675309,
    printFlag = TRUE
  )
  
  imputed_complete_pm10 <- complete(imputed_data_pm10, action = "long", include = FALSE)
  
  cog_long_pm10 <- imputed_complete_pm10 %>%
    pivot_longer(
      cols = starts_with("MeanCognitiveDecline_"),
      names_to = "Year",
      names_prefix = "MeanCognitiveDecline_",
      values_to = "MeanCognitiveDecline"
    )
  
  pm10_long <- imputed_complete_pm10 %>%
    pivot_longer(
      cols = starts_with("Mean_PM10_"),
      names_to = "Year",
      names_prefix = "Mean_PM10_",
      values_to = "Mean_PM10"
    )
  
  pm10_lag1_long <- imputed_complete_pm10 %>%
    pivot_longer(cols = starts_with("PM10_lag1_"), names_to = "Year", names_prefix = "PM10_lag1_", values_to = "PM10_lag1")
  
  pm10_lag2_long <- imputed_complete_pm10 %>%
    pivot_longer(cols = starts_with("PM10_lag2_"), names_to = "Year", names_prefix = "PM10_lag2_", values_to = "PM10_lag2")
  
  pm10_lag3_long <- imputed_complete_pm10 %>%
    pivot_longer(cols = starts_with("PM10_lag3_"), names_to = "Year", names_prefix = "PM10_lag3_", values_to = "PM10_lag3")
  
  imputed_long_pm10 <- cog_long_pm10 %>%
    left_join(pm10_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(pm10_lag1_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(pm10_lag2_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(pm10_lag3_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    select(.imp, .id, StateName, Year, MeanCognitiveDecline, Mean_PM10, PM10_lag1, PM10_lag2, PM10_lag3) %>%
    mutate(Year = as.integer(Year))
  
  # Fit per imputation and pool (no as.mids needed)
  imp_list <- split(imputed_long_pm10, imputed_long_pm10$.imp)
  
  fit_lag1_PM10 <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + PM10_lag1 + (1 | StateName), data = d)
  })
  fit_lag2_PM10 <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + PM10_lag2 + (1 | StateName), data = d)
  })
  fit_lag3_PM10 <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + PM10_lag3 + (1 | StateName), data = d)
  })
  
  print(summary(pool(as.mira(fit_lag1_PM10))))
  print(summary(pool(as.mira(fit_lag2_PM10))))
  print(summary(pool(as.mira(fit_lag3_PM10))))
  
}
